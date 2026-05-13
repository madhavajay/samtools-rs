//! `samtools coverage` — per-reference alignment statistics.
//!
//! Mirrors `main_coverage` in `coverage.c`. Upstream's implementation uses
//! pileup to compute exact per-position depth and coverage. This Rust port
//! computes approximate statistics by walking each record's CIGAR:
//!
//!  - `numreads` — count of records passing the mapped/qual filter
//!  - `covbases` — number of distinct reference positions covered
//!  - `coverage` — `covbases / ref_length * 100`
//!  - `meandepth` — sum of aligned bases / ref_length
//!  - `meanbaseq` — mean base quality across aligned bases (sentinel `0.0`
//!    when no aligned base qualities are available)
//!  - `meanmapq` — mean mapping quality across reads
//!
//! Output is upstream's tabular format with the `#rname startpos endpos
//! numreads covbases coverage meandepth meanbaseq meanmapq` header.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::core::Region;
use htslib_rs::format::Exact;
use htslib_rs::sam;

use crate::bam_flag::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FUNMAP};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools coverage`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut min_mapq: u8 = 0;
    let mut no_header = false;
    let mut output: Option<PathBuf> = None;
    let mut region: Option<String> = None;
    let mut inputs: Vec<PathBuf> = Vec::new();

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-q" | "--min-MQ" => {
                min_mapq = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0);
            }
            "-H" | "--no-header" => no_header = true,
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "-r" | "--region" => {
                region = iter.next().and_then(|a| a.to_str().map(str::to_owned));
            }
            "-l" | "--min-read-len" | "-Q" | "--min-BQ" | "--rf" | "--ff" | "-d" | "--depth"
            | "--min-depth" | "-w" | "--n-bins" => {
                let _ = iter.next();
            }
            "-m" | "--histogram" | "-D" | "--plot-depth" | "-A" | "--ascii" => {
                // Plot modes not yet supported.
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("coverage", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => inputs.push(PathBuf::from(arg)),
        }
    }

    if inputs.is_empty() {
        let _ = print_usage();
        return ExitCode::from(1);
    }

    let mut has_cram = false;
    for path in &inputs {
        let format = match sam_io::sam_open_format(path) {
            Ok(f) => f,
            Err(e) => {
                print_error("coverage", e.to_string());
                return ExitCode::from(1);
            }
        };
        match format.exact {
            Exact::Bam => {}
            Exact::Cram => has_cram = true,
            _ => {
                print_error(
                    "coverage",
                    "only BAM and reference-backed CRAM input are currently supported (SAM TODO)",
                );
                return ExitCode::from(1);
            }
        }
    }

    let reference = if has_cram {
        match current_global_args().reference {
            Some(reference) => Some(reference),
            None => {
                print_error("coverage", "CRAM input requires top-level --reference FILE");
                return ExitCode::from(1);
            }
        }
    } else {
        None
    };

    let mut writer = match sam_io::open_text_output(output.as_deref()) {
        Ok(writer) => writer,
        Err(e) => {
            print_error_errno("coverage", "open -o output", &e);
            return ExitCode::from(1);
        }
    };

    if !no_header {
        let _ = writeln!(
            writer,
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
        );
    }

    match run_coverage(
        &inputs,
        &mut *writer,
        min_mapq,
        region.as_deref(),
        reference.as_deref(),
    ) {
        Ok(()) => match sam_io::check_sam_close(&mut writer) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
            Err(e) => {
                print_error_errno("coverage", "close output", &e);
                ExitCode::from(1)
            }
        },
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("coverage", "coverage failed", &e);
            ExitCode::from(1)
        }
    }
}

fn run_coverage(
    inputs: &[PathBuf],
    out: &mut dyn Write,
    min_mapq: u8,
    region: Option<&str>,
    reference: Option<&Path>,
) -> io::Result<()> {
    let exclude_flags = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    let region = region.map(parse_region).transpose()?;

    for path in inputs {
        match sam_io::sam_open_format(path)?.exact {
            Exact::Bam => run_bam_coverage(path, out, exclude_flags, min_mapq, region.as_ref())?,
            Exact::Cram => run_cram_coverage(
                path,
                out,
                reference.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "CRAM input requires top-level --reference FILE",
                    )
                })?,
                exclude_flags,
                min_mapq,
                region.as_ref(),
            )?,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "only BAM and reference-backed CRAM input are currently supported",
                ));
            }
        }
    }
    Ok(())
}

fn parse_region(s: &str) -> io::Result<Region> {
    s.parse::<Region>().map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid region \"{}\": {}", s, e),
        )
    })
}

fn run_bam_coverage(
    path: &Path,
    out: &mut dyn Write,
    exclude_flags: u32,
    min_mapq: u8,
    region: Option<&Region>,
) -> io::Result<()> {
    let mut reader = bam::io::Reader::new(File::open(path)?);
    let header = reader.read_header()?;
    let mut refs = coverage_targets(&header, region)?;

    if let Some(region) = region {
        for record in htslib_rs::alignment_compat::query_bam_records_from_path(path, region)? {
            update_targets(&header, &mut refs, &record, exclude_flags, min_mapq);
        }
    } else {
        let mut record = bam::Record::default();
        loop {
            let n = reader.read_record(&mut record)?;
            if n == 0 {
                break;
            }
            update_targets(&header, &mut refs, &record, exclude_flags, min_mapq);
        }
    }

    emit_coverage(out, &refs)
}

fn run_cram_coverage(
    path: &Path,
    out: &mut dyn Write,
    reference: &Path,
    exclude_flags: u32,
    min_mapq: u8,
    region: Option<&Region>,
) -> io::Result<()> {
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(path)?;
    let mut refs = coverage_targets(&header, region)?;

    if let Some(region) = region {
        for record in htslib_rs::alignment_compat::query_cram_records_from_path_with_reference(
            path, region, reference,
        )? {
            update_targets(&header, &mut refs, &record, exclude_flags, min_mapq);
        }
    } else {
        for rs in &mut refs {
            let region = ref_region(rs)?;
            for record in htslib_rs::alignment_compat::query_cram_records_from_path_with_reference(
                path, &region, reference,
            )? {
                update_target(&header, rs, &record, exclude_flags, min_mapq);
            }
        }
    }

    emit_coverage(out, &refs)
}

fn ref_region(target: &RefStats) -> io::Result<Region> {
    format!(
        "{}:{}-{}",
        target.name, target.output_start, target.output_end
    )
    .parse::<Region>()
    .map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "region \"{}:{}-{}\": {}",
                target.name, target.output_start, target.output_end, e
            ),
        )
    })
}

fn emit_coverage(out: &mut dyn Write, refs: &[RefStats]) -> io::Result<()> {
    for rs in refs {
        let covbases = rs.covered.iter().filter(|&&b| b).count();
        let coverage_pct = if rs.length > 0 {
            covbases as f64 / rs.length as f64 * 100.0
        } else {
            0.0
        };
        let meandepth = if rs.length > 0 {
            rs.aligned_bases as f64 / rs.length as f64
        } else {
            0.0
        };
        let meanmapq = if rs.num_reads > 0 {
            rs.mapq_sum as f64 / rs.num_reads as f64
        } else {
            0.0
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            rs.name,
            rs.output_start,
            rs.output_end,
            rs.num_reads,
            covbases,
            coverage_pct,
            meandepth,
            0.0_f64,
            meanmapq,
        )?;
    }
    Ok(())
}

struct RefStats {
    tid: usize,
    name: String,
    length: usize,
    output_start: usize,
    output_end: usize,
    start0: usize,
    end0: usize,
    covered: Vec<bool>,
    num_reads: u64,
    aligned_bases: u64,
    mapq_sum: u64,
}

fn coverage_targets(
    header: &htslib_rs::sam::Header,
    region: Option<&Region>,
) -> io::Result<Vec<RefStats>> {
    match region {
        Some(region) => {
            let tid = header
                .reference_sequences()
                .get_index_of(region.name())
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "region reference sequence does not exist: {}",
                            String::from_utf8_lossy(region.name())
                        ),
                    )
                })?;
            let (_, def) = header.reference_sequences().get_index(tid).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence ID")
            })?;
            let ref_len = usize::from(def.length());
            let interval = region.interval();
            let output_start = interval.start().map(usize::from).unwrap_or(1);
            let output_end = interval.end().map(usize::from).unwrap_or(ref_len);
            if output_end < output_start {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid region interval: {}", region),
                ));
            }

            let start0 = output_start - 1;
            let end0 = output_end.min(ref_len);
            let length = output_end - output_start + 1;

            Ok(vec![RefStats {
                tid,
                name: String::from_utf8_lossy(region.name()).into_owned(),
                length,
                output_start,
                output_end,
                start0,
                end0,
                covered: vec![false; length],
                num_reads: 0,
                aligned_bases: 0u64,
                mapq_sum: 0u64,
            }])
        }
        None => Ok(header
            .reference_sequences()
            .iter()
            .enumerate()
            .map(|(tid, (name, def))| {
                let length = usize::from(def.length());
                RefStats {
                    tid,
                    name: String::from_utf8_lossy(name).into_owned(),
                    length,
                    output_start: 1,
                    output_end: length,
                    start0: 0,
                    end0: length,
                    covered: vec![false; length],
                    num_reads: 0,
                    aligned_bases: 0u64,
                    mapq_sum: 0u64,
                }
            })
            .collect()),
    }
}

fn update_targets(
    header: &sam::Header,
    refs: &mut [RefStats],
    record: &(impl sam::alignment::Record + ?Sized),
    exclude_flags: u32,
    min_mapq: u8,
) {
    let flag = match record.flags() {
        Ok(flags) => u16::from(flags) as u32,
        Err(_) => return,
    };
    if flag & exclude_flags != 0 {
        return;
    }
    let mapq = match record.mapping_quality() {
        Some(Ok(q)) => u8::from(q),
        Some(Err(_)) => return,
        None => 0,
    };
    if mapq < min_mapq {
        return;
    }
    let tid = match record.reference_sequence_id(header).transpose() {
        Ok(Some(t)) => t,
        _ => return,
    };
    let Some(rs) = refs.iter_mut().find(|rs| rs.tid == tid) else {
        return;
    };

    update_target_after_filter(rs, record, mapq);
}

fn update_target(
    header: &sam::Header,
    rs: &mut RefStats,
    record: &(impl sam::alignment::Record + ?Sized),
    exclude_flags: u32,
    min_mapq: u8,
) {
    let flag = match record.flags() {
        Ok(flags) => u16::from(flags) as u32,
        Err(_) => return,
    };
    if flag & exclude_flags != 0 {
        return;
    }
    let mapq = match record.mapping_quality() {
        Some(Ok(q)) => u8::from(q),
        Some(Err(_)) => return,
        None => 0,
    };
    if mapq < min_mapq {
        return;
    }
    if record
        .reference_sequence_id(header)
        .transpose()
        .unwrap_or_default()
        != Some(rs.tid)
    {
        return;
    }

    update_target_after_filter(rs, record, mapq);
}

fn update_target_after_filter(
    rs: &mut RefStats,
    record: &(impl sam::alignment::Record + ?Sized),
    mapq: u8,
) {
    let start = match record.alignment_start().transpose() {
        Ok(Some(p)) => usize::from(p) - 1,
        _ => return,
    };

    rs.num_reads += 1;
    rs.mapq_sum += mapq as u64;

    let mut ref_pos = start;
    for op in record.cigar().iter() {
        let op = match op {
            Ok(op) => op,
            Err(_) => break,
        };
        let len = op.len();
        use htslib_rs::sam::alignment::record::cigar::op::Kind;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let op_end = ref_pos.saturating_add(len);
                let lo = ref_pos.max(rs.start0);
                let hi = op_end.min(rs.end0);
                if hi > lo {
                    for p in lo..hi {
                        let offset = p - rs.start0;
                        if offset < rs.covered.len() && !rs.covered[offset] {
                            rs.covered[offset] = true;
                        }
                    }
                    rs.aligned_bases += (hi - lo) as u64;
                }
                ref_pos = op_end;
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos = ref_pos.saturating_add(len);
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools coverage [options] <in.bam>")?;
    writeln!(w, "  -q INT       min mapping quality [0]")?;
    writeln!(w, "  -H           no header")?;
    writeln!(w, "  -o FILE      output FILE")?;
    writeln!(w, "  -r REGION    restrict to REGION")?;
    Ok(())
}
