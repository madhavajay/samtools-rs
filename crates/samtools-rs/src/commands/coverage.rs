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
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::core::Region;
use htslib_rs::format::Exact;
use htslib_rs::sam;

use crate::bam_flag::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FUNMAP, str_to_flag};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools coverage`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut min_mapq: u8 = 0;
    let mut min_baseq: u8 = 0;
    let mut min_depth: u32 = 1;
    let mut max_depth: Option<u32> = None;
    let mut min_read_len: usize = 0;
    let mut no_header = false;
    let mut output: Option<PathBuf> = None;
    let mut region: Option<String> = None;
    let mut inputs: Vec<PathBuf> = Vec::new();
    let mut histogram = false;
    let mut n_bins: usize = 80;
    let mut exclude_flags = default_exclude_flags();
    let mut include_any_flags = 0;

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
            "-Q" | "--min-BQ" => {
                min_baseq = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0);
            }
            "-d" | "--depth" => {
                max_depth = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok());
            }
            "--min-depth" => {
                min_depth = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(1);
            }
            "-l" | "--min-read-len" => {
                min_read_len = iter
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
            "-b" | "--bam-list" => {
                let Some(path) = iter.next().map(PathBuf::from) else {
                    print_error("coverage", "option -b requires an argument");
                    return ExitCode::from(1);
                };
                match read_input_list(&path) {
                    Ok(listed_inputs) => inputs.extend(listed_inputs),
                    Err(e) => {
                        print_error_errno("coverage", "read -b input list", &e);
                        return ExitCode::from(1);
                    }
                }
            }
            "--rf" | "--incl-flags" => match parse_flag_value(iter.next(), s) {
                Ok(flags) => include_any_flags = flags,
                Err(()) => return ExitCode::from(1),
            },
            "--ff" | "--excl-flags" => match parse_flag_value(iter.next(), s) {
                Ok(flags) => exclude_flags = flags,
                Err(()) => return ExitCode::from(1),
            },
            "-w" | "--n-bins" => {
                if let Some(v) = iter.next().and_then(|a| a.to_str())
                    && let Ok(parsed) = v.parse::<usize>()
                    && parsed > 0
                {
                    n_bins = parsed;
                }
            }
            "-m" | "--histogram" | "-A" | "--ascii" => {
                histogram = true;
            }
            "-D" | "--plot-depth" => {
                // Currently routed through the same ASCII histogram as -m;
                // the upstream depth plot variant is not yet implemented.
                histogram = true;
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
            Exact::Sam | Exact::Bam => {}
            Exact::Cram => has_cram = true,
            _ => {
                print_error(
                    "coverage",
                    "only SAM, BAM, and reference-backed CRAM input are currently supported",
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

    if !no_header && !histogram {
        let _ = writeln!(
            writer,
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
        );
    }

    match run_coverage(
        &inputs,
        &mut *writer,
        CoverageConfig {
            min_mapq,
            min_baseq,
            min_depth,
            max_depth,
            min_read_len,
            exclude_flags,
            include_any_flags,
            histogram,
            n_bins,
        },
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

fn default_exclude_flags() -> u32 {
    BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP
}

fn parse_flag_value(value: Option<&OsString>, option: &str) -> Result<u32, ()> {
    let Some(raw) = value.and_then(|a| a.to_str()) else {
        print_error("coverage", format!("option {option} requires an argument"));
        return Err(());
    };
    let Some(flags) = str_to_flag(raw) else {
        print_error("coverage", format!("could not parse {option} {raw}"));
        return Err(());
    };
    Ok(flags as u32)
}

fn read_input_list(path: &Path) -> io::Result<Vec<PathBuf>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut inputs = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end_matches('\r').trim();
        if !trimmed.is_empty() {
            inputs.push(PathBuf::from(trimmed));
        }
    }
    Ok(inputs)
}

#[derive(Clone, Copy)]
struct CoverageConfig {
    min_mapq: u8,
    min_baseq: u8,
    min_depth: u32,
    max_depth: Option<u32>,
    min_read_len: usize,
    exclude_flags: u32,
    include_any_flags: u32,
    histogram: bool,
    n_bins: usize,
}

fn run_coverage(
    inputs: &[PathBuf],
    out: &mut dyn Write,
    config: CoverageConfig,
    region: Option<&str>,
    reference: Option<&Path>,
) -> io::Result<()> {
    let region = region.map(parse_region).transpose()?;

    let mut merged_refs: Option<Vec<RefStats>> = None;
    for path in inputs {
        let refs = match sam_io::sam_open_format(path)?.exact {
            Exact::Sam => collect_sam_coverage(path, config, region.as_ref())?,
            Exact::Bam => collect_bam_coverage(path, config, region.as_ref())?,
            Exact::Cram => collect_cram_coverage(
                path,
                reference.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "CRAM input requires top-level --reference FILE",
                    )
                })?,
                config,
                region.as_ref(),
            )?,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "only SAM, BAM, and reference-backed CRAM input are currently supported",
                ));
            }
        };

        if let Some(merged_refs) = &mut merged_refs {
            merge_coverage_refs(merged_refs, refs)?;
        } else {
            merged_refs = Some(refs);
        }
    }
    if let Some(refs) = merged_refs {
        if config.histogram {
            emit_histogram(
                out,
                &refs,
                config.min_depth,
                config.max_depth,
                config.n_bins,
            )?;
        } else {
            emit_coverage(out, &refs, config.min_depth, config.max_depth)?;
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

fn collect_bam_coverage(
    path: &Path,
    config: CoverageConfig,
    region: Option<&Region>,
) -> io::Result<Vec<RefStats>> {
    let mut reader = bam::io::Reader::new(File::open(path)?);
    let header = reader.read_header()?;
    let mut refs = coverage_targets(&header, region)?;

    if let Some(region) = region {
        for record in htslib_rs::alignment_compat::query_bam_records_from_path(path, region)? {
            update_targets(&header, &mut refs, &record, config);
        }
    } else {
        let mut record = bam::Record::default();
        loop {
            let n = reader.read_record(&mut record)?;
            if n == 0 {
                break;
            }
            update_targets(&header, &mut refs, &record, config);
        }
    }

    Ok(refs)
}

fn collect_sam_coverage(
    path: &Path,
    config: CoverageConfig,
    region: Option<&Region>,
) -> io::Result<Vec<RefStats>> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(path)?));
    let header = reader.read_header()?;
    let mut refs = coverage_targets(&header, region)?;

    for result in reader.records() {
        let record = result?;
        update_targets(&header, &mut refs, &record, config);
    }

    Ok(refs)
}

fn collect_cram_coverage(
    path: &Path,
    reference: &Path,
    config: CoverageConfig,
    region: Option<&Region>,
) -> io::Result<Vec<RefStats>> {
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(path)?;
    let mut refs = coverage_targets(&header, region)?;

    if let Some(region) = region {
        for record in htslib_rs::alignment_compat::query_cram_records_from_path_with_reference(
            path, region, reference,
        )? {
            update_targets(&header, &mut refs, &record, config);
        }
    } else {
        for rs in &mut refs {
            let region = ref_region(rs)?;
            for record in htslib_rs::alignment_compat::query_cram_records_from_path_with_reference(
                path, &region, reference,
            )? {
                update_target(&header, rs, &record, config);
            }
        }
    }

    Ok(refs)
}

fn merge_coverage_refs(merged: &mut [RefStats], next: Vec<RefStats>) -> io::Result<()> {
    if merged.len() != next.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "coverage inputs have incompatible reference dictionaries",
        ));
    }

    for (left, right) in merged.iter_mut().zip(next) {
        if left.name != right.name
            || left.length != right.length
            || left.output_start != right.output_start
            || left.output_end != right.output_end
            || left.start0 != right.start0
            || left.end0 != right.end0
            || left.depths.len() != right.depths.len()
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "coverage inputs have incompatible reference dictionaries",
            ));
        }

        for (left_depth, right_depth) in left.depths.iter_mut().zip(right.depths) {
            *left_depth = left_depth.saturating_add(right_depth);
        }
        left.num_reads = left.num_reads.saturating_add(right.num_reads);
        left.aligned_bases = left.aligned_bases.saturating_add(right.aligned_bases);
        left.baseq_sum = left.baseq_sum.saturating_add(right.baseq_sum);
        left.baseq_count = left.baseq_count.saturating_add(right.baseq_count);
        left.mapq_sum = left.mapq_sum.saturating_add(right.mapq_sum);
    }

    Ok(())
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

fn emit_coverage(
    out: &mut dyn Write,
    refs: &[RefStats],
    min_depth: u32,
    max_depth: Option<u32>,
) -> io::Result<()> {
    for rs in refs {
        let covbases = rs
            .depths
            .iter()
            .map(|&d| cap_depth(d, max_depth))
            .filter(|&d| d >= min_depth)
            .count();
        let coverage_pct = if rs.length > 0 {
            covbases as f64 / rs.length as f64 * 100.0
        } else {
            0.0
        };
        let meandepth = if rs.length > 0 {
            rs.depths
                .iter()
                .map(|&d| cap_depth(d, max_depth) as u64)
                .sum::<u64>() as f64
                / rs.length as f64
        } else {
            0.0
        };
        let meanmapq = if rs.num_reads > 0 {
            rs.mapq_sum as f64 / rs.num_reads as f64
        } else {
            0.0
        };
        let meanbaseq = if rs.baseq_count > 0 {
            rs.baseq_sum as f64 / rs.baseq_count as f64
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
            meanbaseq,
            meanmapq,
        )?;
    }
    Ok(())
}

/// Emits an ASCII histogram per reference/region in the spirit of upstream
/// `samtools coverage -m`. Each row spans `n_bins` columns; the y-axis runs
/// over 10 rows of "percent covered" with `.` and `:` as level marks and `|`
/// as the side border. Byte-for-byte parity with the C tool is **not**
/// pursued here — the UTF-8 box-drawing path, summary sidebar text, and
/// x-axis labels are deliberately omitted; the goal is a serviceable plot
/// that uses the per-position depths we already compute.
fn emit_histogram(
    out: &mut dyn Write,
    refs: &[RefStats],
    min_depth: u32,
    max_depth: Option<u32>,
    n_bins: usize,
) -> io::Result<()> {
    const N_ROWS: usize = 10;

    for rs in refs {
        let region_len = rs.length.max(1);
        let bin_width = ((region_len as f64) / n_bins as f64).max(1.0);

        // Count positions per bin that meet the minimum-depth threshold.
        // hist[i] is in units of "covered positions in bin i".
        let mut hist = vec![0u64; n_bins];
        for (i, &depth) in rs.depths.iter().enumerate() {
            let depth = cap_depth(depth, max_depth);
            if depth >= min_depth {
                let bin = ((i as f64) / bin_width).floor() as usize;
                let bin = bin.min(n_bins.saturating_sub(1));
                hist[bin] += 1;
            }
        }
        // hist_data[i] is the percentage of positions covered within bin i.
        let mut hist_data = vec![0.0f64; n_bins];
        let mut max_val = 0.0f64;
        for i in 0..n_bins {
            hist_data[i] = 100.0 * hist[i] as f64 / bin_width;
            if hist_data[i] > max_val {
                max_val = hist_data[i];
            }
        }

        writeln!(
            out,
            "{} ({}bp)",
            rs.name,
            readable_bp(rs.length as u64).trim_end()
        )?;

        let row_size = (max_val / N_ROWS as f64).max(f64::MIN_POSITIVE);
        for row in (0..N_ROWS).rev() {
            let current = row as f64 * row_size;
            write!(out, ">{:7.2}% |", current)?;
            for &value in hist_data.iter().take(n_bins) {
                let diff = ((value - current) / row_size).round() as i64;
                let glyph = if diff <= 0 {
                    ' '
                } else if diff == 1 {
                    '.'
                } else {
                    ':'
                };
                write!(out, "{}", glyph)?;
            }
            writeln!(out, "|")?;
        }
    }
    Ok(())
}

fn cap_depth(depth: u32, max_depth: Option<u32>) -> u32 {
    max_depth.map_or(depth, |max| depth.min(max))
}

fn readable_bp(bp: u64) -> String {
    if bp >= 1_000_000_000 {
        format!("{:.1}G", bp as f64 / 1_000_000_000.0)
    } else if bp >= 1_000_000 {
        format!("{:.1}M", bp as f64 / 1_000_000.0)
    } else if bp >= 1_000 {
        format!("{:.1}K", bp as f64 / 1_000.0)
    } else {
        format!("{}", bp)
    }
}

struct RefStats {
    tid: usize,
    name: String,
    length: usize,
    output_start: usize,
    output_end: usize,
    start0: usize,
    end0: usize,
    depths: Vec<u32>,
    num_reads: u64,
    aligned_bases: u64,
    baseq_sum: u64,
    baseq_count: u64,
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
                depths: vec![0u32; length],
                num_reads: 0,
                aligned_bases: 0u64,
                baseq_sum: 0u64,
                baseq_count: 0u64,
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
                    depths: vec![0u32; length],
                    num_reads: 0,
                    aligned_bases: 0u64,
                    baseq_sum: 0u64,
                    baseq_count: 0u64,
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
    config: CoverageConfig,
) {
    let flag = match record.flags() {
        Ok(flags) => u16::from(flags) as u32,
        Err(_) => return,
    };
    if !flag_passes(flag, config) {
        return;
    }
    let mapq = match record.mapping_quality() {
        Some(Ok(q)) => u8::from(q),
        Some(Err(_)) => return,
        None => 0,
    };
    if mapq < config.min_mapq {
        return;
    }
    if config.min_read_len != 0
        && read_length_used(record.cigar().iter()).unwrap_or_default() < config.min_read_len
    {
        return;
    }
    let tid = match record.reference_sequence_id(header).transpose() {
        Ok(Some(t)) => t,
        _ => return,
    };
    let Some(rs) = refs.iter_mut().find(|rs| rs.tid == tid) else {
        return;
    };

    update_target_after_filter(rs, record, mapq, config.min_baseq);
}

fn update_target(
    header: &sam::Header,
    rs: &mut RefStats,
    record: &(impl sam::alignment::Record + ?Sized),
    config: CoverageConfig,
) {
    let flag = match record.flags() {
        Ok(flags) => u16::from(flags) as u32,
        Err(_) => return,
    };
    if !flag_passes(flag, config) {
        return;
    }
    let mapq = match record.mapping_quality() {
        Some(Ok(q)) => u8::from(q),
        Some(Err(_)) => return,
        None => 0,
    };
    if mapq < config.min_mapq {
        return;
    }
    if config.min_read_len != 0
        && read_length_used(record.cigar().iter()).unwrap_or_default() < config.min_read_len
    {
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

    update_target_after_filter(rs, record, mapq, config.min_baseq);
}

fn flag_passes(flag: u32, config: CoverageConfig) -> bool {
    if flag & config.exclude_flags != 0 {
        return false;
    }
    if config.include_any_flags != 0 && flag & config.include_any_flags == 0 {
        return false;
    }
    true
}

fn read_length_used(
    cigar: impl Iterator<Item = io::Result<htslib_rs::sam::alignment::record::cigar::Op>>,
) -> io::Result<usize> {
    use htslib_rs::sam::alignment::record::cigar::op::Kind;

    let mut len = 0usize;
    for op in cigar {
        let op = op?;
        match op.kind() {
            Kind::Match | Kind::Insertion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                len = len.saturating_add(op.len());
            }
            Kind::Deletion | Kind::Skip | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
    Ok(len)
}

fn update_target_after_filter(
    rs: &mut RefStats,
    record: &(impl sam::alignment::Record + ?Sized),
    mapq: u8,
    min_baseq: u8,
) {
    let start = match record.alignment_start().transpose() {
        Ok(Some(p)) => usize::from(p) - 1,
        _ => return,
    };

    rs.num_reads += 1;
    rs.mapq_sum += mapq as u64;

    let quality_scores = record.quality_scores();
    let qualities = if quality_scores.is_empty() {
        None
    } else {
        quality_scores.iter().collect::<io::Result<Vec<_>>>().ok()
    };

    let mut ref_pos = start;
    let mut query_pos = 0usize;
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
                        if offset >= rs.depths.len() {
                            continue;
                        }
                        let passes_baseq = if let Some(qualities) = qualities.as_ref() {
                            let query_offset = query_pos + (p - ref_pos);
                            match qualities.get(query_offset) {
                                Some(&quality) if quality >= min_baseq => {
                                    rs.baseq_sum += quality as u64;
                                    rs.baseq_count += 1;
                                    true
                                }
                                Some(_) => false,
                                None => min_baseq == 0,
                            }
                        } else {
                            min_baseq == 0
                        };
                        if passes_baseq {
                            rs.depths[offset] = rs.depths[offset].saturating_add(1);
                            rs.aligned_bases += 1;
                        }
                    }
                }
                ref_pos = op_end;
                query_pos = query_pos.saturating_add(len);
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos = ref_pos.saturating_add(len);
            }
            Kind::Insertion | Kind::SoftClip => {
                query_pos = query_pos.saturating_add(len);
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools coverage [options] <in.bam>")?;
    writeln!(w, "  -q INT       min mapping quality [0]")?;
    writeln!(w, "  -Q INT       min base quality [0]")?;
    writeln!(w, "  -d INT       max per-position depth cap")?;
    writeln!(w, "  -H           no header")?;
    writeln!(w, "  -o FILE      output FILE")?;
    writeln!(w, "  -r REGION    restrict to REGION")?;
    Ok(())
}
