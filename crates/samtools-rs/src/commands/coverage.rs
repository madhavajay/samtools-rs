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
use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
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
    let mut output_mode = CoverageOutputMode::Tabular;
    let mut full_utf = true;
    let mut n_bins: Option<isize> = None;
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
                n_bins = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse::<isize>().ok());
                if matches!(output_mode, CoverageOutputMode::Tabular) {
                    output_mode = CoverageOutputMode::Histogram;
                }
            }
            "-m" | "--histogram" => {
                if matches!(output_mode, CoverageOutputMode::Tabular) {
                    output_mode = CoverageOutputMode::Histogram;
                }
            }
            "-A" | "--ascii" => {
                if matches!(output_mode, CoverageOutputMode::Tabular) {
                    output_mode = CoverageOutputMode::Histogram;
                }
                full_utf = false;
            }
            "-D" | "--plot-depth" => {
                output_mode = CoverageOutputMode::DepthPlot;
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
        if path.as_os_str() != "-" && !path.exists() {
            print_hts_open_missing(path);
            print_error(
                "coverage",
                format!(
                    "Could not open \"{}\": No such file or directory",
                    path.display()
                ),
            );
            return ExitCode::from(1);
        }
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
                    "only SAM, BAM, and CRAM input are currently supported",
                );
                return ExitCode::from(1);
            }
        }
    }
    let reference = has_cram.then(|| current_global_args().reference).flatten();

    let mut writer = match sam_io::open_text_output(output.as_deref()) {
        Ok(writer) => writer,
        Err(e) => {
            print_error_errno("coverage", "open -o output", &e);
            return ExitCode::from(1);
        }
    };

    if !no_header && matches!(output_mode, CoverageOutputMode::Tabular) {
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
            output_mode,
            full_utf,
            n_bins: histogram_bins(n_bins),
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

fn histogram_bins(explicit: Option<isize>) -> usize {
    if let Some(n) = explicit
        && n > 0
    {
        return n as usize;
    }

    let columns = std::env::var("COLUMNS")
        .ok()
        .and_then(|s| s.parse::<isize>().ok())
        .unwrap_or(0);
    if columns > 60 {
        (columns - 40) as usize
    } else {
        40
    }
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
    output_mode: CoverageOutputMode,
    full_utf: bool,
    n_bins: usize,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum CoverageOutputMode {
    Tabular,
    Histogram,
    DepthPlot,
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
            Exact::Cram => collect_cram_coverage(path, reference, config, region.as_ref())?,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "only SAM, BAM, and CRAM input are currently supported",
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
        if matches!(
            config.output_mode,
            CoverageOutputMode::Histogram | CoverageOutputMode::DepthPlot
        ) {
            emit_histogram(
                out,
                &refs,
                config.min_depth,
                config.max_depth,
                config.n_bins,
                config.full_utf,
                config.output_mode == CoverageOutputMode::DepthPlot,
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
    reference: Option<&Path>,
    config: CoverageConfig,
    region: Option<&Region>,
) -> io::Result<Vec<RefStats>> {
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(path)?;
    let mut refs = coverage_targets(&header, region)?;

    if let Some(region) = region {
        for record in query_cram_coverage_records(path, region, reference)? {
            update_targets(&header, &mut refs, &record, config);
        }
    } else {
        for rs in &mut refs {
            let region = ref_region(rs)?;
            for record in query_cram_coverage_records(path, &region, reference)? {
                update_target(&header, rs, &record, config);
            }
        }
    }

    Ok(refs)
}

fn query_cram_coverage_records(
    path: &Path,
    region: &Region,
    reference: Option<&Path>,
) -> io::Result<Vec<sam::alignment::RecordBuf>> {
    if let Some(reference) = reference {
        htslib_rs::alignment_compat::query_cram_records_from_path_with_reference(
            path, region, reference,
        )
    } else {
        htslib_rs::alignment_compat::query_cram_records_from_path_synthesizing_reference(
            path, region,
        )
    }
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
            || left.reference_length != right.reference_length
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
        for (l, r) in left.baseq_pos.iter_mut().zip(right.baseq_pos) {
            *l = l.saturating_add(r);
        }
        for (l, r) in left.baseqn_pos.iter_mut().zip(right.baseqn_pos) {
            *l = l.saturating_add(r);
        }
        left.num_reads = left.num_reads.saturating_add(right.num_reads);
        left.total_reads = left.total_reads.saturating_add(right.total_reads);
        left.aligned_bases = left.aligned_bases.saturating_add(right.aligned_bases);
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

/// C `printf("%.*g")` with `prec` significant digits (default-style, no
/// `#`): pick `%f` when the decimal exponent X satisfies `-4 <= X < prec`,
/// else `%e`; strip trailing zeros and a trailing point. Matches glibc's
/// round-half-to-even (as Rust's float formatting also does).
fn c_printf_g(value: f64, prec: usize) -> String {
    let p = prec.max(1);
    if value == 0.0 {
        return "0".to_string();
    }
    if !value.is_finite() {
        return format!("{value}");
    }

    // Round to `p` significant digits via %e to get the true exponent.
    let sci = format!("{:.*e}", p - 1, value);
    let (mantissa, exp) = sci.split_once('e').unwrap();
    let x: i32 = exp.parse().unwrap();

    let trim = |s: &str| -> String {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s.to_string()
        }
    };

    if x >= -4 && x < p as i32 {
        let decimals = (p as i32 - 1 - x).max(0) as usize;
        trim(&format!("{value:.decimals$}"))
    } else {
        let m = trim(mantissa);
        format!("{}e{}{:02}", m, if x < 0 { "-" } else { "+" }, x.abs())
    }
}

fn emit_coverage(
    out: &mut dyn Write,
    refs: &[RefStats],
    min_depth: u32,
    max_depth: Option<u32>,
) -> io::Result<()> {
    // coverage.c prints references as their pileup columns are reached,
    // flushing any reference with no selected reads at the very end (still
    // in tid order). Reproduce: references with reads first, then the
    // empty ones, each group stable in tid order.
    let ordered = refs
        .iter()
        .filter(|rs| rs.num_reads > 0)
        .chain(refs.iter().filter(|rs| rs.num_reads == 0));
    for rs in ordered {
        // Single pass mirroring coverage.c: a position contributes to
        // covbases / summed_coverage / summed_baseQ only when its (capped)
        // depth clears `min_depth`.
        let mut covbases = 0u64;
        let mut summed_coverage = 0u64;
        let mut summed_baseq = 0u64;
        let mut quality_bases = 0u64;
        for (i, &raw) in rs.depths.iter().enumerate() {
            let d = cap_depth(raw, max_depth);
            if d >= min_depth {
                covbases += 1;
                summed_coverage += d as u64;
                summed_baseq += rs.baseq_pos[i];
                quality_bases += rs.baseqn_pos[i] as u64;
            }
        }
        let coverage_pct = if rs.length > 0 {
            covbases as f64 / rs.length as f64 * 100.0
        } else {
            0.0
        };
        let meandepth = if rs.length > 0 {
            summed_coverage as f64 / rs.length as f64
        } else {
            0.0
        };
        let meanmapq = if rs.num_reads > 0 {
            rs.mapq_sum as f64 / rs.num_reads as f64
        } else {
            0.0
        };
        let meanbaseq = if quality_bases > 0 {
            summed_baseq as f64 / quality_bases as f64
        } else {
            0.0
        };
        // Upstream coverage.c:211 — `\t...\t%g\t%g\t%.3g\t%.3g`.
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rs.name,
            rs.output_start,
            rs.output_end,
            rs.num_reads,
            covbases,
            c_printf_g(coverage_pct, 6),
            c_printf_g(meandepth, 6),
            c_printf_g(meanbaseq, 3),
            c_printf_g(meanmapq, 3),
        )?;
    }
    Ok(())
}

/// Emits the upstream-shaped histogram/depth-plot view used by `coverage -m`,
/// `coverage -A`, and `coverage -D`.
fn emit_histogram(
    out: &mut dyn Write,
    refs: &[RefStats],
    min_depth: u32,
    max_depth: Option<u32>,
    n_bins: usize,
    full_utf: bool,
    plot_depth: bool,
) -> io::Result<()> {
    const N_ROWS: usize = 10;

    for (idx, rs) in refs.iter().enumerate() {
        if idx > 0 {
            writeln!(out)?;
        }
        let region_len = rs.length.max(1);
        let hist_size = n_bins.min(region_len).max(1);
        let bin_width = (region_len / hist_size).max(1);

        let mut hist = vec![0u64; hist_size];
        for (i, &depth) in rs.depths.iter().enumerate() {
            let depth = cap_depth(depth, max_depth);
            let bin = i / bin_width;
            if bin >= hist_size {
                continue;
            }
            if plot_depth {
                hist[bin] = hist[bin].saturating_add(depth as u64);
            } else if depth >= min_depth {
                hist[bin] += 1;
            }
        }
        let mut hist_data = vec![0.0f64; hist_size];
        let mut max_val = 0.0f64;
        for i in 0..hist_size {
            hist_data[i] = if plot_depth { 1.0 } else { 100.0 } * hist[i] as f64 / bin_width as f64;
            if hist_data[i] > max_val {
                max_val = hist_data[i];
            }
        }

        writeln!(
            out,
            "{} ({}bp)",
            rs.name,
            readable_bp_c(rs.reference_length as f64).trim_end()
        )?;

        let row_size = max_val / N_ROWS as f64;
        let render = HistRenderContext {
            min_depth,
            max_depth,
            bin_width,
            max_val,
            plot_depth,
        };
        for row in (0..N_ROWS).rev() {
            let current = row as f64 * row_size;
            if plot_depth {
                write!(out, ">{:8.1} ", current)?;
            } else {
                write!(out, ">{:7.2}% ", current)?;
            }
            write!(out, "{}", if full_utf { "│" } else { "|" })?;
            for &value in hist_data.iter().take(hist_size) {
                let diff = if row_size > 0.0 {
                    ((blockchar_len(full_utf) as f64 * (value - current) / row_size).round()
                        as isize)
                        - 1
                } else {
                    -1
                };
                if diff < 0 {
                    write!(out, " ")?;
                } else {
                    write!(out, "{}", block_char(diff as usize, full_utf))?;
                }
            }
            write!(out, "{}", if full_utf { "│" } else { "|" })?;
            write_hist_sidebar(out, row, rs, render)?;
            writeln!(out)?;
        }
        write_hist_axis(out, rs, bin_width, hist_size)?;
    }
    Ok(())
}

fn blockchar_len(full_utf: bool) -> usize {
    if full_utf { 8 } else { 2 }
}

fn block_char(diff: usize, full_utf: bool) -> &'static str {
    if full_utf {
        const BLOCKS: [&str; 8] = ["▁", "▂", "▃", "▄", "▅", "▆", "▇", "█"];
        BLOCKS[diff.min(BLOCKS.len() - 1)]
    } else {
        const BLOCKS: [&str; 2] = [".", ":"];
        BLOCKS[diff.min(BLOCKS.len() - 1)]
    }
}

#[derive(Clone, Copy)]
struct HistRenderContext {
    min_depth: u32,
    max_depth: Option<u32>,
    bin_width: usize,
    max_val: f64,
    plot_depth: bool,
}

fn write_hist_sidebar(
    out: &mut dyn Write,
    row: usize,
    rs: &RefStats,
    render: HistRenderContext,
) -> io::Result<()> {
    let summary = coverage_summary(rs, render.min_depth, render.max_depth);
    write!(out, " ")?;
    match row {
        9 => write!(out, "Number of reads: {}", rs.num_reads)?,
        8 => {
            let filtered = rs.total_reads.saturating_sub(rs.num_reads);
            if filtered > 0 {
                write!(out, "    ({} filtered)", filtered)?;
            }
        }
        7 => write!(
            out,
            "Covered bases:   {}bp",
            readable_bp_c(summary.covbases as f64)
        )?,
        6 => write!(
            out,
            "Percent covered: {}%",
            c_printf_g(summary.coverage_pct, 4)
        )?,
        5 => write!(
            out,
            "Mean coverage:   {}x",
            c_printf_g(summary.meandepth, 3)
        )?,
        4 => write!(out, "Mean baseQ:      {}", c_printf_g(summary.meanbaseq, 3))?,
        3 => write!(out, "Mean mapQ:       {}", c_printf_g(summary.meanmapq, 3))?,
        1 => write!(
            out,
            "Histo bin width: {}bp",
            readable_bp_c(render.bin_width as f64)
        )?,
        0 => {
            if render.plot_depth {
                write!(out, "Histo max cov:   {}", c_printf_g(render.max_val, 5))?;
            } else {
                write!(out, "Histo max bin:   {}%", c_printf_g(render.max_val, 5))?;
            }
        }
        _ => {}
    }
    Ok(())
}

fn write_hist_axis(
    out: &mut dyn Write,
    rs: &RefStats,
    bin_width: usize,
    hist_size: usize,
) -> io::Result<()> {
    write!(
        out,
        "     {}",
        center_text(&readable_bp_c(rs.output_start as f64), 10)
    )?;
    for rest in (10..10 * (hist_size / 10)).step_by(10) {
        write!(
            out,
            "{}",
            center_text(
                &readable_bp_c((rs.output_start - 1 + bin_width * rest) as f64),
                10
            )
        )?;
    }
    let last_padding = hist_size % 10;
    write!(
        out,
        "{:>width$}{}",
        " ",
        center_text(&readable_bp_c(rs.output_end as f64), 10),
        width = last_padding
    )?;
    writeln!(out)
}

fn center_text(text: &str, width: usize) -> String {
    let len = text.len();
    if len >= width {
        return text.to_string();
    }
    let padding = (width - len) / 2;
    let padding_ex = (width - len) % 2;
    if padding >= 1 {
        format!(
            " {:>text_width$}{:>right_width$}",
            text,
            "",
            text_width = len + padding,
            right_width = padding - 1 + padding_ex
        )
    } else {
        text.to_string()
    }
}

struct CoverageSummary {
    covbases: u64,
    coverage_pct: f64,
    meandepth: f64,
    meanbaseq: f64,
    meanmapq: f64,
}

fn coverage_summary(rs: &RefStats, min_depth: u32, max_depth: Option<u32>) -> CoverageSummary {
    let mut covbases = 0u64;
    let mut summed_coverage = 0u64;
    let mut summed_baseq = 0u64;
    let mut quality_bases = 0u64;
    for (i, &raw) in rs.depths.iter().enumerate() {
        let d = cap_depth(raw, max_depth);
        if d >= min_depth {
            covbases += 1;
            summed_coverage += d as u64;
            summed_baseq += rs.baseq_pos[i];
            quality_bases += rs.baseqn_pos[i] as u64;
        }
    }
    CoverageSummary {
        covbases,
        coverage_pct: if rs.length > 0 {
            covbases as f64 / rs.length as f64 * 100.0
        } else {
            0.0
        },
        meandepth: if rs.length > 0 {
            summed_coverage as f64 / rs.length as f64
        } else {
            0.0
        },
        meanbaseq: if quality_bases > 0 {
            summed_baseq as f64 / quality_bases as f64
        } else {
            0.0
        },
        meanmapq: if rs.num_reads > 0 {
            rs.mapq_sum as f64 / rs.num_reads as f64
        } else {
            0.0
        },
    }
}

fn cap_depth(depth: u32, max_depth: Option<u32>) -> u32 {
    max_depth.map_or(depth, |max| depth.min(max))
}

fn readable_bp_c(mut bp: f64) -> String {
    let units = ["", "K", "M", "G", "T"];
    let mut unit = 0usize;
    while bp >= 1000.0 && unit < units.len() - 1 {
        bp /= 1000.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{bp:.0}{}", units[unit])
    } else {
        format!("{bp:.unit$}{}", units[unit])
    }
}

struct RefStats {
    tid: usize,
    name: String,
    reference_length: usize,
    length: usize,
    output_start: usize,
    output_end: usize,
    start0: usize,
    end0: usize,
    depths: Vec<u32>,
    total_reads: u64,
    num_reads: u64,
    aligned_bases: u64,
    /// Per-position summed base quality / count of quality bases, so the
    /// `meandepth`/`meanbaseq` accumulators can be gated by `min_depth`
    /// exactly like upstream `coverage.c` (only positions whose depth
    /// clears `min_depth` contribute).
    baseq_pos: Vec<u64>,
    baseqn_pos: Vec<u32>,
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
                reference_length: ref_len,
                length,
                output_start,
                output_end,
                start0,
                end0,
                depths: vec![0u32; length],
                total_reads: 0,
                num_reads: 0,
                aligned_bases: 0u64,
                baseq_pos: vec![0u64; length],
                baseqn_pos: vec![0u32; length],
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
                    reference_length: length,
                    length,
                    output_start: 1,
                    output_end: length,
                    start0: 0,
                    end0: length,
                    depths: vec![0u32; length],
                    total_reads: 0,
                    num_reads: 0,
                    aligned_bases: 0u64,
                    baseq_pos: vec![0u64; length],
                    baseqn_pos: vec![0u32; length],
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
    let tid = match record.reference_sequence_id(header).transpose() {
        Ok(Some(t)) => t,
        _ => return,
    };
    let Some(rs) = refs.iter_mut().find(|rs| rs.tid == tid) else {
        return;
    };
    rs.total_reads += 1;
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
    if record
        .reference_sequence_id(header)
        .transpose()
        .unwrap_or_default()
        != Some(rs.tid)
    {
        return;
    }
    rs.total_reads += 1;
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
                                    rs.baseq_pos[offset] += quality as u64;
                                    rs.baseqn_pos[offset] += 1;
                                    true
                                }
                                Some(_) => false,
                                None => min_baseq == 0,
                            }
                        } else {
                            const MISSING_QUALITY: u8 = 0xff;
                            rs.baseq_pos[offset] += MISSING_QUALITY as u64;
                            rs.baseqn_pos[offset] += 1;
                            true
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
    writeln!(w, "  -m           show coverage histogram")?;
    writeln!(w, "  -D           plot per-bin mean depth")?;
    writeln!(w, "  -A           use ASCII histogram characters")?;
    writeln!(w, "  -H           no header")?;
    writeln!(w, "  -o FILE      output FILE")?;
    writeln!(w, "  -r REGION    restrict to REGION")?;
    Ok(())
}
