//! `samtools view` — SAM/BAM/CRAM conversion and filtering.
//!
//! This is the anchor subcommand and the upstream `sam_view.c` is the
//! largest file in samtools (68k LOC). The current Rust port covers the
//! common conversion and counting paths required by the basic
//! `test_view` cases in `samtools/test/test.pl`. Filters, region queries,
//! BED files, aux-tag manipulation, and the long tail of flags are still
//! TODO.
//!
//! Supported flags so far:
//!  - `-h` — include header in SAM output (default).
//!  - `-H` — print header only.
//!  - `-b` — write BAM output.
//!  - `-C` — write CRAM output (requires `-T`/`--reference`).
//!  - `-S` — input is SAM (ignored; format is auto-detected).
//!  - `-P` — accepted for multi-region BAM output; duplicate records are preserved.
//!  - `--no-PG` — do not add a `@PG` line.
//!  - `-c` — count matching records and print the count.
//!  - `-o FILE` — write output to FILE (default stdout).
//!  - `-T FILE` / `--reference FILE` — reference for CRAM I/O.
//!  - `-u` — write uncompressed BAM (accepted; treated as `-b -1` for now).
//!  - `-1` — fast compression level (accepted; treated as `-b` default for now).
//!
//! Anything else returns a "not yet supported" error so that test failures
//! are loud rather than silent.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{Category, Exact, detect_path};

use crate::diagnostics::{print_error, print_error_errno};
use crate::header_text::read_raw_header_text_with_format;

/// Entry point for `samtools view`.
pub fn main(args: &[OsString]) -> ExitCode {
    let opts = match parse_args(args) {
        Ok(o) => o,
        Err(ParseError::Usage) => {
            let _ = write_usage(&mut io::stdout());
            return ExitCode::SUCCESS;
        }
        Err(ParseError::Err(msg)) => {
            print_error("view", msg);
            let _ = write_usage(&mut io::stderr());
            return ExitCode::from(1);
        }
    };

    let Some(input) = opts.input.clone() else {
        print_error("view", "no input file (stdin not yet supported)");
        return ExitCode::from(1);
    };

    let format = match detect_path(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error(
                "view",
                format!("failed to open \"{}\": {}", input.display(), e),
            );
            return ExitCode::from(1);
        }
    };
    if format.category != Category::SequenceData {
        print_error("view", format!("{} is not sequence data", input.display()));
        return ExitCode::from(1);
    }

    // If a -L BED file was given, expand its entries into region strings.
    let mut opts = opts;
    if let Some(bed) = opts.bed_path.clone() {
        match load_bed_regions(&bed) {
            Ok(more) => opts.regions.extend(more),
            Err(e) => {
                print_error_errno("view", format!("failed to read \"{}\"", bed.display()), &e);
                return ExitCode::from(1);
            }
        }
    }

    match run(&opts, &input, format.exact) {
        Ok(code) => code,
        // Broken pipe from a downstream consumer (e.g. `samtools view | head`)
        // is a clean exit, not an error — matches upstream behavior.
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("view", "I/O error during view", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Default)]
struct Opts {
    input: Option<PathBuf>,
    output: Option<PathBuf>,
    output_fmt: OutputFmt,
    header: HeaderMode,
    count: bool,
    no_pg: bool,
    reference: Option<PathBuf>,
    /// `-f INT` — require ALL these flag bits to be set on the record.
    require_flags: u32,
    /// `-F INT` — exclude records with ANY of these flag bits set.
    exclude_flags: u32,
    /// `-G INT` — exclude records with ALL of these flag bits set (alternative
    /// "exclude" path; upstream uses it for excluding a particular flag combo).
    exclude_all_flags: u32,
    /// `-q INT` — minimum mapping quality (records with MAPQ < threshold are skipped).
    min_mapq: u8,
    /// Positional region strings after the input file (`samtools view file.bam ref1 ref2:100-200`).
    regions: Vec<String>,
    /// `-e EXPR` — HTSlib filter expression evaluated per record.
    filter_expr: Option<String>,
    /// `-L FILE` — BED file; only records overlapping a BED interval pass.
    bed_path: Option<PathBuf>,
    /// `-x TAG` (repeatable) — aux tags to strip from SAM-output records.
    remove_tags: Vec<String>,
    /// `--keep-tag TAG` (repeatable) — only these aux tags are kept.
    keep_tags: Vec<String>,
}

#[derive(Default, Clone, Copy, PartialEq, Eq)]
enum OutputFmt {
    /// Inherit a sensible default — SAM when writing to stdout/SAM extension,
    /// BAM otherwise. The current port always emits SAM unless explicitly told.
    #[default]
    Auto,
    Sam,
    Bam,
    Cram,
}

#[derive(Default, Clone, Copy, PartialEq, Eq)]
enum HeaderMode {
    /// Default for SAM output: no header. Default for binary output:
    /// header is always included.
    #[default]
    Default,
    /// `-h`: force-include header in SAM output.
    Include,
    /// `-H`: header only, no records.
    HeaderOnly,
}

enum ParseError {
    Usage,
    Err(String),
}

fn parse_flag_arg(args: &[OsString], i: usize, name: &str) -> Result<u32, ParseError> {
    let v = args
        .get(i)
        .and_then(|a| a.to_str())
        .ok_or_else(|| ParseError::Err(format!("missing value for {}", name)))?;
    crate::bam_flag::str_to_flag(v)
        .map(|x| x as u32)
        .ok_or_else(|| ParseError::Err(format!("Could not parse \"{}\"", v)))
}

/// Returns true iff the record's `flag` and `mapq` pass the filter mix.
fn record_passes(flag: u32, mapq: u8, opts: &Opts) -> bool {
    if opts.require_flags != 0 && (flag & opts.require_flags) != opts.require_flags {
        return false;
    }
    if opts.exclude_flags != 0 && (flag & opts.exclude_flags) != 0 {
        return false;
    }
    if opts.exclude_all_flags != 0 && (flag & opts.exclude_all_flags) == opts.exclude_all_flags {
        return false;
    }
    if mapq < opts.min_mapq {
        return false;
    }
    true
}

fn parse_args(args: &[OsString]) -> Result<Opts, ParseError> {
    let mut opts = Opts::default();
    let mut i = 1;
    while i < args.len() {
        let arg = &args[i];
        let Some(s) = arg.to_str() else {
            if opts.input.is_none() {
                opts.input = Some(PathBuf::from(arg));
                i += 1;
                continue;
            }
            return Err(ParseError::Err("too many positional arguments".to_string()));
        };

        match s {
            "-h" => {
                opts.header = HeaderMode::Include;
                i += 1;
            }
            "-H" => {
                opts.header = HeaderMode::HeaderOnly;
                i += 1;
            }
            "-S" => {
                i += 1;
            }
            "-P" => {
                i += 1;
            }
            "-b" => {
                opts.output_fmt = OutputFmt::Bam;
                i += 1;
            }
            "-C" => {
                opts.output_fmt = OutputFmt::Cram;
                i += 1;
            }
            "--sam" => {
                opts.output_fmt = OutputFmt::Sam;
                i += 1;
            }
            "--bam" => {
                opts.output_fmt = OutputFmt::Bam;
                i += 1;
            }
            "--cram" => {
                opts.output_fmt = OutputFmt::Cram;
                i += 1;
            }
            "-u" | "-1" => {
                if matches!(opts.output_fmt, OutputFmt::Auto | OutputFmt::Sam) {
                    opts.output_fmt = OutputFmt::Bam;
                }
                i += 1;
            }
            "-c" => {
                opts.count = true;
                i += 1;
            }
            "--no-PG" => {
                opts.no_pg = true;
                i += 1;
            }
            "-o" | "--output" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| ParseError::Err("missing value for -o".into()))?;
                opts.output = Some(PathBuf::from(v));
                i += 1;
            }
            "-f" => {
                i += 1;
                opts.require_flags = parse_flag_arg(args, i, "-f")?;
                i += 1;
            }
            "-F" => {
                i += 1;
                opts.exclude_flags = parse_flag_arg(args, i, "-F")?;
                i += 1;
            }
            "-G" => {
                i += 1;
                opts.exclude_all_flags = parse_flag_arg(args, i, "-G")?;
                i += 1;
            }
            "-q" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -q".into()))?;
                opts.min_mapq = v
                    .parse()
                    .map_err(|_| ParseError::Err(format!("invalid -q value \"{}\"", v)))?;
                i += 1;
            }
            "-e" | "--expr" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -e".into()))?;
                opts.filter_expr = Some(v.to_string());
                i += 1;
            }
            "-x" | "--remove-tag" | "--remove-tags" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -x".into()))?;
                if let Some(rest) = v.strip_prefix('^') {
                    for tag in rest.split(',') {
                        opts.keep_tags.push(tag.to_string());
                    }
                } else {
                    for tag in v.split(',') {
                        opts.remove_tags.push(tag.to_string());
                    }
                }
                i += 1;
            }
            "--keep-tag" | "--keep-tags" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for --keep-tag".into()))?;
                for tag in v.split(',') {
                    opts.keep_tags.push(tag.to_string());
                }
                i += 1;
            }
            "-L" | "--target-file" => {
                i += 1;
                opts.bed_path = args
                    .get(i)
                    .map(PathBuf::from)
                    .or_else(|| Some(PathBuf::new()))
                    .filter(|p| !p.as_os_str().is_empty());
                if opts.bed_path.is_none() {
                    return Err(ParseError::Err("missing value for -L".into()));
                }
                i += 1;
            }
            "-T" | "--reference" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| ParseError::Err("missing value for -T".into()))?;
                opts.reference = Some(PathBuf::from(v));
                i += 1;
            }
            "-O" | "--output-fmt" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -O".into()))?;
                // Accept formats like "cram", "bam", "sam", optionally with
                // ",opt=value" suffixes; we only honor the top-level format.
                let head = v.split(',').next().unwrap_or("").to_lowercase();
                opts.output_fmt = match head.as_str() {
                    "sam" => OutputFmt::Sam,
                    "bam" => OutputFmt::Bam,
                    "cram" => OutputFmt::Cram,
                    _ => {
                        return Err(ParseError::Err(format!(
                            "unsupported --output-fmt \"{}\"",
                            v
                        )));
                    }
                };
                i += 1;
            }
            "--help" => return Err(ParseError::Usage),
            _ if s.starts_with('-') && s != "-" => {
                return Err(ParseError::Err(format!(
                    "option `{}` is not yet supported in samtools-rs view",
                    s
                )));
            }
            _ => {
                if opts.input.is_none() {
                    opts.input = Some(PathBuf::from(arg));
                } else {
                    opts.regions.push(s.to_string());
                }
                i += 1;
            }
        }
    }
    Ok(opts)
}

fn run(opts: &Opts, input: &Path, input_exact: Exact) -> io::Result<ExitCode> {
    let effective_out_fmt = resolved_output_fmt(opts);

    // Count-only mode.
    if opts.count {
        let n = if let Some(expr) = opts.filter_expr.as_ref() {
            if input_exact != Exact::Sam {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "`-e EXPR` is currently wired up for SAM input only",
                ));
            }
            htslib_rs::alignment_compat::count_sam_records_matching_filter_from_path(input, expr)?
        } else if !opts.regions.is_empty() {
            count_region_records(input, input_exact, opts)?
        } else if has_filters(opts) {
            count_filtered_records(input, input_exact, opts)?
        } else {
            count_records(input, input_exact)?
        };
        let mut out = open_text_output(opts)?;
        writeln!(out, "{}", n)?;
        return Ok(ExitCode::SUCCESS);
    }

    // Header-only mode.
    if opts.header == HeaderMode::HeaderOnly {
        let header_text = read_raw_header_text_with_format(input, input_exact)?;
        let mut out = open_text_output(opts)?;
        out.write_all(header_text.as_bytes())?;
        return Ok(ExitCode::SUCCESS);
    }

    // SAM output paths.
    if effective_out_fmt == OutputFmt::Sam {
        let mut out = open_text_output(opts)?;
        // Upstream default for SAM output is records-only; `-h` opts in
        // to including the header.
        let include_header = match opts.header {
            HeaderMode::Include => true,
            HeaderMode::Default => false,
            HeaderMode::HeaderOnly => true, // handled above
        };
        if include_header {
            let header_text = read_raw_header_text_with_format(input, input_exact)?;
            out.write_all(header_text.as_bytes())?;
        }
        write_records_as_sam(&mut out, input, input_exact, opts)?;
        return Ok(ExitCode::SUCCESS);
    }

    // BAM output.
    if effective_out_fmt == OutputFmt::Bam {
        let dst = opts
            .output
            .clone()
            .ok_or_else(|| io::Error::other("BAM output to stdout requires -o file (TODO)"))?;
        let dst_file = File::create(&dst)?;
        match input_exact {
            Exact::Sam => {
                htslib_rs::alignment_compat::write_bam_from_sam_path(input, dst_file)?;
            }
            Exact::Bam => {
                if opts.regions.is_empty() {
                    htslib_rs::alignment_compat::write_bam_from_path(input, dst_file)?;
                } else {
                    let regions = parse_region_strings(input, &opts.regions)?;
                    htslib_rs::alignment_compat::write_bam_regions_from_path(
                        input, &regions, dst_file,
                    )?;
                }
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "CRAM->BAM via samtools-rs view is not yet wired up",
                ));
            }
        }
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Cram {
        let dst = opts
            .output
            .clone()
            .ok_or_else(|| io::Error::other("CRAM output to stdout requires -o file (TODO)"))?;
        let reference = opts.reference.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "CRAM output requires --reference / -T",
            )
        })?;
        let dst_file = File::create(&dst)?;
        match input_exact {
            Exact::Sam => {
                htslib_rs::alignment_compat::write_cram_from_sam_path_with_reference(
                    input, reference, dst_file,
                )?;
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "BAM/CRAM -> CRAM via samtools-rs view is not yet wired up",
                ));
            }
        }
        return Ok(ExitCode::SUCCESS);
    }

    Err(io::Error::other("unsupported output combination"))
}

fn resolved_output_fmt(opts: &Opts) -> OutputFmt {
    if opts.output_fmt != OutputFmt::Auto {
        return opts.output_fmt;
    }
    // Auto: infer from output extension if any, else SAM (stdout default).
    if let Some(p) = opts.output.as_ref()
        && let Some(ext) = p.extension().and_then(|e| e.to_str())
    {
        return match ext {
            "sam" => OutputFmt::Sam,
            "bam" => OutputFmt::Bam,
            "cram" => OutputFmt::Cram,
            _ => OutputFmt::Sam,
        };
    }
    OutputFmt::Sam
}

fn open_text_output(opts: &Opts) -> io::Result<Box<dyn Write>> {
    match opts.output.as_ref() {
        Some(p) => Ok(Box::new(File::create(p)?)),
        None => Ok(Box::new(io::stdout().lock())),
    }
}

fn load_bed_regions(path: &Path) -> io::Result<Vec<String>> {
    let file = File::open(path)?;
    let mut out = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        let s = line.trim_end();
        if s.is_empty()
            || s.starts_with('#')
            || s.starts_with("track ")
            || s.starts_with("browser ")
        {
            continue;
        }
        let mut fields = s.split('\t');
        let chrom = fields.next().unwrap_or("");
        let beg: u64 = fields.next().and_then(|t| t.parse().ok()).unwrap_or(0);
        let end: u64 = fields.next().and_then(|t| t.parse().ok()).unwrap_or(0);
        if chrom.is_empty() || end <= beg {
            continue;
        }
        // HTSlib region format is 1-based inclusive.
        out.push(format!("{}:{}-{}", chrom, beg + 1, end));
    }
    Ok(out)
}

fn parse_region_strings(
    _path: &Path,
    regions: &[String],
) -> io::Result<Vec<htslib_rs::core::Region>> {
    regions
        .iter()
        .map(|s| {
            s.parse::<htslib_rs::core::Region>().map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("region \"{}\": {}", s, e),
                )
            })
        })
        .collect()
}

fn count_region_records(path: &Path, exact: Exact, opts: &Opts) -> io::Result<usize> {
    if exact != Exact::Bam {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "region count is only wired up for BAM input",
        ));
    }
    let regions = parse_region_strings(path, &opts.regions)?;
    let mut n = 0usize;
    for region in &regions {
        n += htslib_rs::alignment_compat::count_bam_records_in_region_from_path(path, region)?;
    }
    Ok(n)
}

fn count_records(path: &Path, exact: Exact) -> io::Result<usize> {
    match exact {
        Exact::Sam => htslib_rs::alignment_compat::count_sam_records_from_path(path),
        Exact::Bam => htslib_rs::alignment_compat::count_bam_records_from_path(path),
        Exact::Cram => Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "CRAM counting via samtools-rs view is not yet wired up",
        )),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported format",
        )),
    }
}

fn write_records_as_sam<W: Write>(
    out: &mut W,
    path: &Path,
    exact: Exact,
    opts: &Opts,
) -> io::Result<()> {
    match exact {
        Exact::Sam => {
            if !opts.regions.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "region queries on SAM input require an index (BAM/CRAM only)",
                ));
            }
            stream_sam_records(out, path, opts)
        }
        Exact::Bam => {
            let text = if opts.regions.is_empty() {
                htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(path, None)?
            } else {
                let regions = parse_region_strings(path, &opts.regions)?;
                htslib_rs::alignment_compat::view_bam_regions_as_sam_text_from_path(path, &regions)?
            };
            let tail = strip_header_lines(text.as_bytes());
            // For BAM input we already have SAM text. Apply filters
            // line-by-line if any are set.
            if has_filters(opts) {
                for line in tail.split(|&b| b == b'\n') {
                    if line.is_empty() {
                        continue;
                    }
                    if line_passes(line, opts) {
                        out.write_all(line)?;
                        out.write_all(b"\n")?;
                    }
                }
                Ok(())
            } else {
                out.write_all(tail)
            }
        }
        Exact::Cram => Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "CRAM -> SAM via samtools-rs view is not yet wired up",
        )),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported format",
        )),
    }
}

fn has_filters(opts: &Opts) -> bool {
    opts.require_flags != 0
        || opts.exclude_flags != 0
        || opts.exclude_all_flags != 0
        || opts.min_mapq != 0
}

fn has_tag_filter(opts: &Opts) -> bool {
    !opts.remove_tags.is_empty() || !opts.keep_tags.is_empty()
}

/// Strip / keep aux tags from a SAM record line (without trailing newline).
fn apply_tag_filter(line: &[u8], opts: &Opts) -> Vec<u8> {
    if !has_tag_filter(opts) {
        return line.to_vec();
    }
    let mut out: Vec<u8> = Vec::with_capacity(line.len());
    for (i, field) in line.split(|&b| b == b'\t').enumerate() {
        if i > 0 {
            out.push(b'\t');
        }
        if i < 11 {
            out.extend_from_slice(field);
            continue;
        }
        if field.len() < 5 || field[2] != b':' {
            out.extend_from_slice(field);
            continue;
        }
        let tag = std::str::from_utf8(&field[..2]).unwrap_or("");
        let keep = if !opts.keep_tags.is_empty() {
            opts.keep_tags.iter().any(|t| t == tag)
        } else {
            !opts.remove_tags.iter().any(|t| t == tag)
        };
        if keep {
            out.extend_from_slice(field);
        } else {
            // drop the trailing tab we already wrote
            if out.last() == Some(&b'\t') {
                out.pop();
            }
        }
    }
    out
}

/// Apply view filters to a SAM record line, returning whether it should
/// be emitted. Parses the flag (column 2) and MAPQ (column 5).
fn line_passes(line: &[u8], opts: &Opts) -> bool {
    let mut fields = line.split(|&b| b == b'\t');
    let _qname = fields.next();
    let flag = fields
        .next()
        .and_then(|f| std::str::from_utf8(f).ok())
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(0);
    let _rname = fields.next();
    let _pos = fields.next();
    let mapq = fields
        .next()
        .and_then(|f| std::str::from_utf8(f).ok())
        .and_then(|s| s.parse::<u8>().ok())
        .unwrap_or(0);
    record_passes(flag, mapq, opts)
}

fn count_filtered_records(path: &Path, exact: Exact, opts: &Opts) -> io::Result<usize> {
    let mut count = 0usize;
    match exact {
        Exact::Sam => {
            let file = File::open(path)?;
            let reader: Box<dyn BufRead> = if is_bgzf_path(path)? {
                Box::new(BufReader::new(MultiGzDecoder::new(file)))
            } else {
                Box::new(BufReader::new(file))
            };
            let mut reader = reader;
            let mut line: Vec<u8> = Vec::with_capacity(1024);
            let mut in_records = false;
            loop {
                line.clear();
                let n = reader.read_until(b'\n', &mut line)?;
                if n == 0 {
                    break;
                }
                if !in_records {
                    if line.starts_with(b"@") {
                        continue;
                    }
                    in_records = true;
                }
                let line_no_nl = if line.last() == Some(&b'\n') {
                    &line[..line.len() - 1]
                } else {
                    &line[..]
                };
                if line_passes(line_no_nl, opts) {
                    count += 1;
                }
            }
        }
        Exact::Bam => {
            let text =
                htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(path, None)?;
            let tail = strip_header_lines(text.as_bytes());
            for line in tail.split(|&b| b == b'\n') {
                if line.is_empty() {
                    continue;
                }
                if line_passes(line, opts) {
                    count += 1;
                }
            }
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "filtered counts for this format are not yet wired up",
            ));
        }
    }
    Ok(count)
}

fn stream_sam_records<W: Write>(out: &mut W, path: &Path, opts: &Opts) -> io::Result<()> {
    let file = File::open(path)?;
    let bgzf = is_bgzf_path(path)?;
    let reader: Box<dyn BufRead> = if bgzf {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut reader = reader;
    let mut line: Vec<u8> = Vec::with_capacity(1024);
    let mut in_records = false;
    let filtering = has_filters(opts);
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        if !in_records {
            if line.starts_with(b"@") {
                continue;
            }
            in_records = true;
        }
        let body = if line.last() == Some(&b'\n') {
            &line[..line.len() - 1]
        } else {
            &line[..]
        };
        if filtering && !line_passes(body, opts) {
            continue;
        }
        if has_tag_filter(opts) {
            let filtered = apply_tag_filter(body, opts);
            out.write_all(&filtered)?;
            out.write_all(b"\n")?;
        } else {
            out.write_all(&line)?;
        }
    }
    Ok(())
}

fn is_bgzf_path(path: &Path) -> io::Result<bool> {
    let mut file = File::open(path)?;
    let mut hdr = [0u8; 2];
    let n = file.read(&mut hdr)?;
    Ok(n >= 2 && hdr[0] == 0x1f && hdr[1] == 0x8b)
}

fn strip_header_lines(bytes: &[u8]) -> &[u8] {
    let mut tail = bytes;
    while let Some(pos) = memchr::memchr(b'\n', tail) {
        if tail.starts_with(b"@") {
            tail = &tail[pos + 1..];
        } else {
            break;
        }
    }
    tail
}

fn write_usage<W: Write>(w: &mut W) -> io::Result<()> {
    writeln!(
        w,
        "Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram>"
    )?;
    writeln!(w, "  -h           include header in SAM output")?;
    writeln!(w, "  -H           print header only")?;
    writeln!(w, "  -b           output BAM")?;
    writeln!(w, "  -C           output CRAM (requires -T)")?;
    writeln!(w, "  -c           count records and print the count")?;
    writeln!(w, "  -o FILE      output FILE")?;
    writeln!(w, "  -T FILE      reference (for CRAM)")?;
    writeln!(w, "  --no-PG      do not add a @PG line")?;
    writeln!(w)?;
    writeln!(
        w,
        "Note: this samtools-rs port currently implements only the most common"
    )?;
    writeln!(
        w,
        "      conversion paths. Filters, region queries, BED files, aux-tag"
    )?;
    writeln!(w, "      manipulation, etc. are not yet supported.")?;
    Ok(())
}
