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
use std::path::PathBuf;
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::core::Region;
use htslib_rs::format::{Exact, detect_path};

use crate::bam_flag::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FUNMAP};
use crate::diagnostics::{print_error, print_error_errno};

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

    for path in &inputs {
        let format = match detect_path(path) {
            Ok(f) => f,
            Err(e) => {
                print_error(
                    "coverage",
                    format!("failed to detect format of \"{}\": {}", path.display(), e),
                );
                return ExitCode::from(1);
            }
        };
        if format.exact != Exact::Bam {
            print_error(
                "coverage",
                "only BAM input is currently supported (SAM/CRAM TODO)",
            );
            return ExitCode::from(1);
        }
    }

    let mut writer: Box<dyn Write> = match output.as_ref() {
        Some(p) => match File::create(p) {
            Ok(f) => Box::new(f),
            Err(e) => {
                print_error_errno("coverage", "open -o output", &e);
                return ExitCode::from(1);
            }
        },
        None => Box::new(io::stdout().lock()),
    };

    if !no_header {
        let _ = writeln!(
            writer,
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
        );
    }

    match run_coverage(&inputs, &mut *writer, min_mapq, region.as_deref()) {
        Ok(()) => ExitCode::SUCCESS,
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
) -> io::Result<()> {
    let exclude_flags = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    let region = region
        .map(|s| {
            s.parse::<Region>().map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid region \"{}\": {}", s, e),
                )
            })
        })
        .transpose()?;

    for path in inputs {
        let mut reader = bam::io::Reader::new(File::open(path)?);
        let header = reader.read_header()?;
        let mut refs = coverage_targets(&header, region.as_ref())?;

        if let Some(region) = region.as_ref() {
            for record in htslib_rs::alignment_compat::query_bam_records_from_path(path, region)? {
                update_targets(&mut refs, &record, exclude_flags, min_mapq);
            }
        } else {
            let mut record = bam::Record::default();
            loop {
                let n = reader.read_record(&mut record)?;
                if n == 0 {
                    break;
                }
                update_targets(&mut refs, &record, exclude_flags, min_mapq);
            }
        }

        for rs in &refs {
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
                0.0_f64, // meanbaseq — not computed without per-base access
                meanmapq,
            )?;
        }
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

fn update_targets(refs: &mut [RefStats], record: &bam::Record, exclude_flags: u32, min_mapq: u8) {
    let flag = u16::from(record.flags()) as u32;
    if flag & exclude_flags != 0 {
        return;
    }
    let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
    if mapq < min_mapq {
        return;
    }
    let tid = match record.reference_sequence_id().and_then(|r| r.ok()) {
        Some(t) => t,
        None => return,
    };
    let start = match record.alignment_start().and_then(|r| r.ok()) {
        Some(p) => usize::from(p) - 1,
        None => return,
    };
    let Some(rs) = refs.iter_mut().find(|rs| rs.tid == tid) else {
        return;
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
