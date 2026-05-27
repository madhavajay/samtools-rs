//! `samtools view` — SAM/BAM/CRAM conversion and filtering.
//!
//! This is the anchor subcommand and the upstream `sam_view.c` is the
//! largest file in samtools (68k LOC). The current Rust port covers the
//! common conversion and counting paths required by the basic
//! `test_view` cases in `repos/samtools/test/test.pl`. Filters, region queries,
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
//!  - `-M` — accepted for multi-region iterator mode; current region paths
//!    already intersect `-L` and positional regions without duplicate collapsing.
//!  - `--no-PG` — do not add a `@PG` line.
//!  - `-c` — count matching records and print the count.
//!  - `-o FILE` — write output to FILE (default stdout).
//!  - `-U FILE` — for SAM output, write records not selected by flag/MAPQ filters to FILE.
//!  - `-p` — for SAM output, set UNMAP on records not selected by flag/MAPQ filters.
//!  - `-T FILE` / `--reference FILE` — reference for CRAM I/O.
//!  - `-t FILE` — reference `.fai` index used to add missing `@SQ` header lines.
//!  - `-u` — write uncompressed BAM (accepted; treated as `-b -1` for now).
//!  - `-1` — fast compression level (accepted; treated as `-b` default for now).
//!
//! Anything else returns a "not yet supported" error so that test failures
//! are loud rather than silent.

use std::collections::HashSet;
use std::error::Error;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::bgzf;
use htslib_rs::format::{Category, Exact};
use htslib_rs::sam::{self, alignment::RecordBuf};
use md5::{Digest, Md5};

use crate::aux_list::{AuxTag, parse_aux_list};
use crate::bedidx::load_bed_index;
use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
use crate::header_text::read_raw_header_text_with_format;
use crate::io as sam_io;
use crate::sam_global::current_global_args;
use crate::sanitize::{SanitizeFlags, parse_sanitize_options, sanitize_record};

#[derive(Debug)]
struct SamParseNoSqError {
    path: PathBuf,
    line: usize,
}

impl std::fmt::Display for SamParseNoSqError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("no SQ lines present in the header")
    }
}

impl Error for SamParseNoSqError {}

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

    let stdin_input = opts
        .input
        .as_deref()
        .is_none_or(|path| path == Path::new("-"));

    let mut opts = opts;
    opts.argv = Some(args.to_vec());
    if let Some(bed) = opts.bed_path.clone() {
        match load_bed_regions(&bed) {
            Ok(regions) => opts.bed_regions = regions,
            Err(e) => {
                print_error_errno("view", format!("failed to read \"{}\"", bed.display()), &e);
                return ExitCode::from(1);
            }
        }
    }

    if stdin_input {
        let mut data = Vec::new();
        if let Err(e) = io::stdin().lock().read_to_end(&mut data) {
            print_error_errno("view", "failed to read stdin", &e);
            return ExitCode::from(1);
        }

        if let Some(lib) = opts.library.clone() {
            let header_text: io::Result<String> = match stdin_format(&data) {
                StdinFormat::Sam => {
                    Ok(String::from_utf8_lossy(sam_header_lines(&data)).into_owned())
                }
                StdinFormat::Bam => htslib_rs::alignment_compat::view_bam_as_sam_text(
                    io::Cursor::new(&data),
                    Some(0),
                ),
                StdinFormat::Cram => cram_reference(&opts).and_then(|reference| {
                    htslib_rs::alignment_compat::view_cram_as_sam_text_with_reference(
                        io::Cursor::new(&data),
                        reference.as_path(),
                        Some(0),
                    )
                }),
            };
            match header_text {
                Ok(text) => {
                    opts.library_rg_ids = library_rg_ids_from_header(&text, &lib);
                }
                Err(e) => {
                    print_error_errno("view", "failed to read header for -l", &e);
                    return ExitCode::from(1);
                }
            }
        }

        let result = match stdin_format(&data) {
            StdinFormat::Sam => run_sam_stdin(&opts, &data),
            StdinFormat::Bam => run_bam_stdin(&opts, &data),
            StdinFormat::Cram => run_cram_stdin(&opts, &data),
        };

        return match result {
            Ok(code) => code,
            Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
            Err(e) => {
                print_error_errno("view", "I/O error during view", &e);
                ExitCode::from(1)
            }
        };
    }

    let input = opts.input.clone().expect("non-stdin input exists");
    if !input.exists() {
        print_hts_open_missing(&input);
        print_error(
            "view",
            format!(
                "failed to open \"{}\" for reading: No such file or directory",
                input.display()
            ),
        );
        return ExitCode::from(1);
    }
    let format = match sam_io::sam_open_format(&input) {
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

    if let Some(lib) = opts.library.clone() {
        match read_raw_header_text_with_format(&input, format.exact) {
            Ok(header_text) => {
                opts.library_rg_ids = library_rg_ids_from_header(&header_text, &lib);
            }
            Err(e) => {
                print_error_errno(
                    "view",
                    format!("failed to read header of \"{}\"", input.display()),
                    &e,
                );
                return ExitCode::from(1);
            }
        }
    }

    match run(&opts, &input, format.exact) {
        Ok(code) => {
            if opts.write_index
                && let Err(e) = write_output_index(&opts)
            {
                print_error_errno("view", "failed to write index", &e);
                return ExitCode::from(1);
            }
            if let Err(e) = write_saved_counts(&opts, &input, format.exact) {
                print_error_errno("view", "failed to write save-counts", &e);
                return ExitCode::from(1);
            }
            code
        }
        // Broken pipe from a downstream consumer (e.g. `samtools view | head`)
        // is a clean exit, not an error — matches upstream behavior.
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) if sam_parse_no_sq_error(&e).is_some() => {
            if let Some(err) = sam_parse_no_sq_error(&e) {
                print_sam_parse_no_sq_error(err);
            }
            ExitCode::from(1)
        }
        Err(e) => {
            print_error_errno("view", "I/O error during view", &e);
            ExitCode::from(1)
        }
    }
}

fn sam_parse_no_sq_error(err: &io::Error) -> Option<&SamParseNoSqError> {
    err.get_ref()?.downcast_ref::<SamParseNoSqError>()
}

fn print_sam_parse_no_sq_error(err: &SamParseNoSqError) {
    let _ = writeln!(
        io::stderr(),
        "[E::sam_parse1] no SQ lines present in the header"
    );
    let _ = writeln!(
        io::stderr(),
        "[W::sam_read1_sam] Parse error at line {}",
        err.line
    );
    print_error(
        "view",
        format!("error reading file \"{}\"", err.path.display()),
    );
}

fn write_saved_counts(opts: &Opts, input: &Path, exact: Exact) -> io::Result<()> {
    let Some(path) = opts.save_counts.as_deref() else {
        return Ok(());
    };
    let (processed, accepted) = if opts.fetch_pairs {
        let (_, processed, accepted) = fetch_pairs_sam_text(input, exact, opts)?;
        (processed, accepted)
    } else {
        (
            count_processed_records(input, exact, opts)?,
            count_selected_records(input, exact, opts)?,
        )
    };
    let mut out = File::create(path)?;
    writeln!(out, "{{")?;
    writeln!(out, "    \"records_processed\" : {processed},")?;
    writeln!(out, "    \"records_filter_accepted\" : {accepted},")?;
    writeln!(
        out,
        "    \"records_filter_rejected\" : {}",
        processed.saturating_sub(accepted)
    )?;
    writeln!(out, "}}")?;
    out.flush()
}

fn count_processed_records(input: &Path, exact: Exact, opts: &Opts) -> io::Result<usize> {
    let mut processed_opts = opts.clone();
    clear_selection_filters(&mut processed_opts);
    if has_requested_regions(&processed_opts) {
        count_region_records(input, exact, &processed_opts)
    } else {
        count_records(input, exact, &processed_opts)
    }
}

fn count_selected_records(input: &Path, exact: Exact, opts: &Opts) -> io::Result<usize> {
    if has_requested_regions(opts) {
        count_region_records(input, exact, opts)
    } else if let Some(expr) = opts.filter_expr.as_deref() {
        count_expr_records(input, exact, opts, expr)
    } else if has_filters(opts) {
        count_filtered_records(input, exact, opts)
    } else {
        count_records(input, exact, opts)
    }
}

fn clear_selection_filters(opts: &mut Opts) {
    opts.require_flags = 0;
    opts.exclude_flags = 0;
    opts.exclude_all_flags = 0;
    opts.min_mapq = 0;
    opts.min_query_len = 0;
    opts.qname_filter = None;
    opts.read_groups.clear();
    opts.exclude_no_rg = false;
    opts.library = None;
    opts.library_rg_ids.clear();
    opts.aux_tag_filter = None;
    opts.filter_expr = None;
    opts.only_unplaced = false;
}

fn has_requested_regions(opts: &Opts) -> bool {
    !opts.regions.is_empty() || !opts.bed_regions.is_empty()
}

/// Post-write index pass for `--write-index`.
fn write_output_index(opts: &Opts) -> io::Result<()> {
    let Some(out) = opts.output.as_deref() else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--write-index requires -o FILE",
        ));
    };
    match resolved_output_fmt(opts)? {
        OutputFmt::Bam => {
            let index = htslib_rs::index_compat::build_bam_csi(out)?;
            let mut idx = out.as_os_str().to_os_string();
            idx.push(".csi");
            htslib_rs::index_compat::write_csi(PathBuf::from(idx), &index)
        }
        OutputFmt::Sam => {
            let index = htslib_rs::index_compat::build_sam_csi(out)?;
            let mut idx = out.as_os_str().to_os_string();
            idx.push(".csi");
            htslib_rs::index_compat::write_csi(PathBuf::from(idx), &index)
        }
        OutputFmt::Cram => {
            let index = htslib_rs::index_compat::build_cram_crai(out)?;
            let mut idx = out.as_os_str().to_os_string();
            idx.push(".crai");
            htslib_rs::index_compat::write_cram_crai(PathBuf::from(idx), &index)
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--write-index is only supported for SAM/BAM/CRAM file output in samtools-rs view",
        )),
    }
}

#[derive(Clone, Default)]
struct Opts {
    input: Option<PathBuf>,
    output: Option<PathBuf>,
    unselected_output: Option<PathBuf>,
    output_fmt: OutputFmt,
    /// `-O cram,embed_ref=1`: embed the reference in CRAM output.
    embed_reference: bool,
    /// `-O cram,seqs_per_slice=N`: max records per CRAM slice
    /// (`None` keeps the encoder default).
    records_per_slice: Option<usize>,
    /// `-O cram,slices_per_slice=N`: max slices per CRAM container
    /// (`None` keeps the encoder default).
    slices_per_container: Option<usize>,
    header: HeaderMode,
    count: bool,
    no_pg: bool,
    /// Argv captured for `@PG` insertion when `--no-PG` is not set.
    /// `None` means the caller didn't supply an argv (e.g. internal tests).
    argv: Option<Vec<OsString>>,
    unmap_unselected: bool,
    /// `--write-index` — build an index next to a BAM file output.
    write_index: bool,
    /// `--save-counts FILE` — write processed/accepted/rejected filter counts.
    save_counts: Option<PathBuf>,
    /// `--fetch-pairs` — include mates of records fetched by region.
    fetch_pairs: bool,
    /// Region `*` — emit only unplaced (RNAME `*`) records.
    only_unplaced: bool,
    reference: Option<PathBuf>,
    /// `-t FILE` — FASTA index used to add missing `@SQ` lines.
    reference_index: Option<PathBuf>,
    /// `-f INT` — require ALL these flag bits to be set on the record.
    require_flags: u32,
    /// `-F INT` — exclude records with ANY of these flag bits set.
    exclude_flags: u32,
    /// `-G INT` — exclude records with ALL of these flag bits set (alternative
    /// "exclude" path; upstream uses it for excluding a particular flag combo).
    exclude_all_flags: u32,
    /// `-q INT` — minimum mapping quality (records with MAPQ < threshold are skipped).
    min_mapq: u8,
    /// `-m INT` — minimum query length from CIGAR query-consuming ops.
    min_query_len: usize,
    /// Positional region strings after the input file (`samtools view file.bam ref1 ref2:100-200`).
    regions: Vec<String>,
    /// `-e EXPR` — HTSlib filter expression evaluated per record.
    filter_expr: Option<String>,
    /// `-L FILE` — BED file; only records overlapping a BED interval pass.
    bed_path: Option<PathBuf>,
    /// Parsed `-L FILE` BED intervals, kept separate from positional regions
    /// so `-L bed chr:start-end` intersects them instead of unioning them.
    bed_regions: Vec<String>,
    /// `-x TAG` (repeatable) — aux tags to strip from SAM-output records.
    remove_tags: Vec<AuxTag>,
    /// `--keep-tag TAG` (repeatable) — only these aux tags are kept.
    keep_tags: Vec<AuxTag>,
    /// `-z FLAGS` / `--sanitize FLAGS` — upstream-style record sanitizer.
    sanitize_flags: SanitizeFlags,
    /// `-B` / `--remove-B` — collapse backward CIGAR operations.
    remove_b: bool,
    /// `--remove-flags FLAG` — clear these flag bits before writing.
    remove_flags: u32,
    /// `-s INT.FRAC` — deterministic subsampling by read name.
    subsample: Option<Subsample>,
    /// `-N FILE` / `--qname-file FILE` — read names listed in FILE (or
    /// `^FILE` to negate). Records whose qname appears in the set pass;
    /// `^FILE` flips to exclude. `None` means the filter is disabled.
    qname_filter: Option<QnameFilter>,
    /// `-r STR` / `-R FILE` — accumulated read-group IDs. Records with a
    /// matching `RG:Z:` aux value pass. Records with no `RG` aux tag also
    /// pass unless `-n` is set, matching upstream `samtools view`.
    read_groups: HashSet<Vec<u8>>,
    /// `-n` — exclude records that have no `RG:Z:` aux tag at all.
    exclude_no_rg: bool,
    /// `-l STR` / `--library STR` — only output records whose read group's
    /// `@RG LB:` value equals STR. `None` means the filter is off.
    library: Option<String>,
    /// Resolved from the header once the input is known: the set of `@RG`
    /// IDs whose `LB:` matches [`Opts::library`]. A record passes the
    /// library filter iff its `RG:Z:` value is in this set (so a record
    /// with no `RG` is excluded, matching upstream `bam_get_library`).
    library_rg_ids: HashSet<Vec<u8>>,
    /// `-X` / `--customized-index` — legacy synopsis where the second
    /// positional is an explicit index path
    /// (`view -X in.bam in.bam.bai [region…]`).
    customized_index: bool,
    /// Explicit index path captured under `-X`. Accepted for synopsis
    /// compatibility; our region queries build/find the index
    /// themselves, so this is currently informational (a no-op, matching
    /// `idxstats -X`).
    index_path: Option<PathBuf>,
    /// `-d TAG[:VAL]` / `-D TAG:FILE` — aux-tag presence (or value-set)
    /// filter. All `-d` / `-D` invocations must share the same TAG, and
    /// values accumulate into the same `AuxTagFilter`. `None` means the
    /// filter is off.
    aux_tag_filter: Option<AuxTagFilter>,
}

#[derive(Clone, Debug)]
struct AuxTagFilter {
    tag: [u8; 2],
    /// If `None`, any record carrying `tag` passes (presence-only).
    /// Otherwise the record's value for that tag must appear in the set.
    values: Option<HashSet<Vec<u8>>>,
}

#[derive(Clone, Copy, Debug)]
struct Subsample {
    seed: u32,
    fraction: f64,
}

#[derive(Clone, Debug)]
struct QnameFilter {
    negate: bool,
    names: HashSet<Vec<u8>>,
}

impl QnameFilter {
    fn matches(&self, qname: &[u8]) -> bool {
        let contained = self.names.contains(qname);
        if self.negate { !contained } else { contained }
    }
}

/// Scans the `@RG` lines of a SAM header text and returns the set of
/// read-group IDs whose `LB:` field equals `library` (upstream
/// `bam_get_library`/`-l` semantics). Only the leading `@`-prefixed
/// header block is inspected.
fn library_rg_ids_from_header(header_text: &str, library: &str) -> HashSet<Vec<u8>> {
    let mut ids = HashSet::new();
    for line in header_text.lines() {
        if !line.starts_with('@') {
            break;
        }
        if !line.starts_with("@RG\t") {
            continue;
        }
        let mut id = None;
        let mut lb = None;
        for field in line.split('\t').skip(1) {
            if let Some(v) = field.strip_prefix("ID:") {
                id = Some(v);
            } else if let Some(v) = field.strip_prefix("LB:") {
                lb = Some(v);
            }
        }
        if lb == Some(library)
            && let Some(id) = id
        {
            ids.insert(id.as_bytes().to_vec());
        }
    }
    ids
}

/// Walk a SAM record line and return the `RG:Z:` value, or `None` if not
/// present. Skips the first 11 mandatory fields.
fn extract_rg_value(line: &[u8]) -> Option<&[u8]> {
    for (i, field) in line.split(|&b| b == b'\t').enumerate() {
        if i < 11 {
            continue;
        }
        if field.len() >= 5 && field.starts_with(b"RG:Z:") {
            return Some(&field[5..]);
        }
    }
    None
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

#[derive(Debug)]
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

fn simple_filter_expr(opts: &Opts) -> Option<String> {
    let mut parts = Vec::new();
    if opts.require_flags != 0 {
        parts.push(format!(
            "(flag & {}) == {}",
            opts.require_flags, opts.require_flags
        ));
    }
    if opts.exclude_flags != 0 {
        parts.push(format!("(flag & {}) == 0", opts.exclude_flags));
    }
    if opts.exclude_all_flags != 0 {
        parts.push(format!(
            "(flag & {}) != {}",
            opts.exclude_all_flags, opts.exclude_all_flags
        ));
    }
    if opts.min_mapq != 0 {
        parts.push(format!("mapq >= {}", opts.min_mapq));
    }

    (!parts.is_empty()).then(|| parts.join(" && "))
}

fn combined_filter_expr(opts: &Opts) -> Option<String> {
    match (opts.filter_expr.as_deref(), simple_filter_expr(opts)) {
        (Some(expr), Some(simple)) => Some(format!("({expr}) && ({simple})")),
        (Some(expr), None) => Some(expr.to_string()),
        (None, Some(simple)) => Some(simple),
        (None, None) => None,
    }
}

fn prefilter_expr_for_sam_output(opts: &Opts) -> Option<String> {
    if opts.unselected_output.is_some() || opts.unmap_unselected {
        None
    } else {
        combined_filter_expr(opts)
    }
}

fn opts_after_prefiltered_expr(opts: &Opts) -> Opts {
    let mut opts = opts.clone();
    opts.filter_expr = None;
    opts
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
            if opts.customized_index && opts.index_path.is_none() {
                opts.index_path = Some(PathBuf::from(arg));
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
            "-P" | "-M" => {
                i += 1;
            }
            "--fetch-pairs" => {
                opts.fetch_pairs = true;
                i += 1;
            }
            "-B" | "--remove-B" => {
                opts.remove_b = true;
                i += 1;
            }
            "-s" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -s".into()))?;
                opts.subsample = Some(parse_subsample(v)?);
                i += 1;
            }
            _ if s.starts_with("-s") && s.len() > 2 => {
                opts.subsample = Some(parse_subsample(&s[2..])?);
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
            "-p" | "--unmap" => {
                opts.unmap_unselected = true;
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
            "-ho" => {
                opts.header = HeaderMode::Include;
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| ParseError::Err("missing value for -ho".into()))?;
                opts.output = Some(PathBuf::from(v));
                i += 1;
            }
            "-U" | "--unoutput" | "--output-unselected" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| ParseError::Err("missing value for -U".into()))?;
                opts.unselected_output = Some(PathBuf::from(v));
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
            "--exclude-flags" => {
                i += 1;
                opts.exclude_flags = parse_flag_arg(args, i, "--exclude-flags")?;
                i += 1;
            }
            "--remove-flags" => {
                i += 1;
                opts.remove_flags = parse_flag_arg(args, i, "--remove-flags")?;
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
            "-m" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -m".into()))?;
                opts.min_query_len = v
                    .parse()
                    .map_err(|_| ParseError::Err(format!("invalid -m value \"{}\"", v)))?;
                i += 1;
            }
            _ if s.starts_with("-m") && s.len() > 2 => {
                let v = &s[2..];
                opts.min_query_len = v
                    .parse()
                    .map_err(|_| ParseError::Err(format!("invalid -m value \"{}\"", v)))?;
                i += 1;
            }
            "-N" | "--qname-file" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -N".into()))?;
                let (path, negate) = match v.strip_prefix('^') {
                    Some(rest) => (rest, true),
                    None => (v, false),
                };
                let text = std::fs::read_to_string(path).map_err(|e| {
                    ParseError::Err(format!("failed to read -N value \"{}\": {}", path, e))
                })?;
                let names: HashSet<Vec<u8>> = text
                    .lines()
                    .map(|line| line.split_ascii_whitespace().next().unwrap_or(line))
                    .filter(|s| !s.is_empty())
                    .map(|s| s.as_bytes().to_vec())
                    .collect();
                opts.qname_filter = Some(QnameFilter { negate, names });
                i += 1;
            }
            "-r" | "--read-group" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -r".into()))?;
                opts.read_groups.insert(v.as_bytes().to_vec());
                i += 1;
            }
            "-R" | "--read-group-file" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -R".into()))?;
                let text = std::fs::read_to_string(v).map_err(|e| {
                    ParseError::Err(format!("failed to read -R value \"{}\": {}", v, e))
                })?;
                for line in text.lines() {
                    let id = line.split_ascii_whitespace().next().unwrap_or(line);
                    if !id.is_empty() {
                        opts.read_groups.insert(id.as_bytes().to_vec());
                    }
                }
                i += 1;
            }
            "-n" => {
                opts.exclude_no_rg = true;
                i += 1;
            }
            "-l" | "--library" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -l".into()))?;
                opts.library = Some(v.to_string());
                i += 1;
            }
            "-d" | "--tag" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -d".into()))?;
                if v.len() < 2 {
                    return Err(ParseError::Err(format!(
                        "invalid -d value \"{}\": expected TAG[:VALUE]",
                        v
                    )));
                }
                let tag = [v.as_bytes()[0], v.as_bytes()[1]];
                let value = if v.len() > 2 {
                    if v.as_bytes()[2] != b':' {
                        return Err(ParseError::Err(format!(
                            "invalid -d value \"{}\": expected TAG:VALUE separator",
                            v
                        )));
                    }
                    Some(v.as_bytes()[3..].to_vec())
                } else {
                    None
                };
                if let Some(existing) = opts.aux_tag_filter.as_mut() {
                    if existing.tag != tag {
                        return Err(ParseError::Err(format!(
                            "different tag \"{}{}\" specified after \"{}{}\"",
                            char::from(tag[0]),
                            char::from(tag[1]),
                            char::from(existing.tag[0]),
                            char::from(existing.tag[1]),
                        )));
                    }
                    if let Some(value) = value {
                        existing
                            .values
                            .get_or_insert_with(HashSet::new)
                            .insert(value);
                    }
                } else {
                    opts.aux_tag_filter = Some(AuxTagFilter {
                        tag,
                        values: value.map(|v| {
                            let mut set = HashSet::new();
                            set.insert(v);
                            set
                        }),
                    });
                }
                i += 1;
            }
            "-D" | "--tag-file" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -D".into()))?;
                if v.len() < 4 || (v.as_bytes()[2] != b':' && v.as_bytes()[2] != b';') {
                    return Err(ParseError::Err(format!(
                        "invalid -D value \"{}\": expected TAG:FILE",
                        v
                    )));
                }
                let tag = [v.as_bytes()[0], v.as_bytes()[1]];
                let path = &v[3..];
                let text = std::fs::read_to_string(path).map_err(|e| {
                    ParseError::Err(format!("failed to read -D value \"{}\": {}", path, e))
                })?;
                let new_values: HashSet<Vec<u8>> = text
                    .lines()
                    .filter(|line| !line.is_empty())
                    .map(|line| line.as_bytes().to_vec())
                    .collect();
                if let Some(existing) = opts.aux_tag_filter.as_mut() {
                    if existing.tag != tag {
                        return Err(ParseError::Err(format!(
                            "different tag \"{}{}\" specified after \"{}{}\"",
                            char::from(tag[0]),
                            char::from(tag[1]),
                            char::from(existing.tag[0]),
                            char::from(existing.tag[1]),
                        )));
                    }
                    existing
                        .values
                        .get_or_insert_with(HashSet::new)
                        .extend(new_values);
                } else {
                    opts.aux_tag_filter = Some(AuxTagFilter {
                        tag,
                        values: Some(new_values),
                    });
                }
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
                    extend_aux_tags(&mut opts.keep_tags, rest, "-x")?;
                } else {
                    extend_aux_tags(&mut opts.remove_tags, v, "-x")?;
                }
                i += 1;
            }
            "--keep-tag" | "--keep-tags" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for --keep-tag".into()))?;
                extend_aux_tags(&mut opts.keep_tags, v, "--keep-tag")?;
                i += 1;
            }
            "-z" | "--sanitize" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err(format!("missing value for {}", s)))?;
                opts.sanitize_flags = parse_sanitize_options(v).map_err(ParseError::Err)?;
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
            "-t" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| ParseError::Err("missing value for -t".into()))?;
                opts.reference_index = Some(PathBuf::from(v));
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
                let mut parts = v.split(',');
                let head = parts.next().unwrap_or("").to_lowercase();
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
                for opt in parts {
                    let opt = opt.trim();
                    if opt.is_empty() {
                        continue;
                    }
                    // Honor the options we support (embed_ref incl.
                    // embed_ref=1/2, seqs_per_slice, slices_per_slice,
                    // version, level); silently ignore other suffixes,
                    // matching the previous lenient `-O` behavior.
                    let _ = apply_output_fmt_option(&mut opts, opt);
                }
                i += 1;
            }
            "--output-fmt-option" => {
                i += 1;
                let v = args.get(i).and_then(|a| a.to_str()).ok_or_else(|| {
                    ParseError::Err("missing value for --output-fmt-option".into())
                })?;
                apply_output_fmt_option(&mut opts, v)?;
                i += 1;
            }
            _ if s.starts_with("--output-fmt-option=") => {
                let v = s
                    .split_once('=')
                    .map(|(_, value)| value)
                    .unwrap_or_default();
                apply_output_fmt_option(&mut opts, v)?;
                i += 1;
            }
            "-X" | "--customized-index" => {
                opts.customized_index = true;
                i += 1;
            }
            "--write-index" => {
                opts.write_index = true;
                i += 1;
            }
            "--save-counts" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| ParseError::Err("missing value for --save-counts".into()))?;
                opts.save_counts = Some(PathBuf::from(v));
                i += 1;
            }
            _ if s.starts_with("--save-counts=") => {
                let v = s
                    .split_once('=')
                    .map(|(_, value)| value)
                    .unwrap_or_default();
                if v.is_empty() {
                    return Err(ParseError::Err("missing value for --save-counts".into()));
                }
                opts.save_counts = Some(PathBuf::from(v));
                i += 1;
            }
            "--help" => return Err(ParseError::Usage),
            // Thread count: accepted and recorded. Output is byte-identical
            // regardless of the value (worker-pool wiring is a perf-only
            // follow-up — completed library batch #8); `-@ N`, `-@N`, `--threads N`.
            "-@" | "--threads" => {
                i += 1;
                let _ = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -@".into()))?;
                i += 1;
            }
            _ if s.starts_with("-@") && s.len() > 2 => {
                // Attached form `-@N`.
                i += 1;
            }
            _ if s.starts_with("--threads=") => {
                i += 1;
            }
            _ if s.starts_with('-') && s != "-" => {
                return Err(ParseError::Err(format!(
                    "option `{}` is not yet supported in samtools-rs view",
                    s
                )));
            }
            _ => {
                if opts.input.is_none() {
                    opts.input = Some(PathBuf::from(arg));
                } else if opts.customized_index && opts.index_path.is_none() {
                    // Legacy `-X` synopsis: the second positional is the
                    // explicit index path. Accepted as a no-op.
                    opts.index_path = Some(PathBuf::from(arg));
                } else {
                    opts.regions.push(s.to_string());
                }
                i += 1;
            }
        }
    }

    // HTSlib region grammar: `.` means "everything" — equivalent to no
    // region restriction (a whole-file pass). Drop it so the no-region
    // code paths handle it.
    opts.regions.retain(|r| r != ".");

    // HTSlib region grammar: `*` selects only unplaced ("no coordinate")
    // reads (RNAME `*`). Treat it as a whole-file pass with an
    // unplaced-only filter rather than a noodles region query.
    if opts.regions.iter().any(|r| r == "*") {
        opts.only_unplaced = true;
        opts.regions.retain(|r| r != "*");
    }

    Ok(opts)
}

fn apply_output_fmt_option(opts: &mut Opts, raw: &str) -> Result<(), ParseError> {
    let (key, value) = match raw.split_once('=') {
        Some((key, value)) => (key.trim(), Some(value.trim())),
        None => (raw.trim(), None),
    };

    fn parse_count(key: &str, value: Option<&str>) -> Result<usize, ParseError> {
        let value = value.ok_or_else(|| {
            ParseError::Err(format!("missing value for --output-fmt-option {key}"))
        })?;
        value
            .parse::<usize>()
            .ok()
            .filter(|n| *n > 0)
            .ok_or_else(|| ParseError::Err(format!("invalid {key} value \"{value}\"")))
    }

    match key {
        "embed_ref" => {
            opts.embed_reference = true;
            Ok(())
        }
        // CRAM slice/container sizing (htslib `seqs_per_slice` /
        // `slices_per_slice`), forwarded to the noodles CRAM encoder.
        "seqs_per_slice" => {
            opts.records_per_slice = Some(parse_count(key, value)?);
            Ok(())
        }
        "slices_per_slice" => {
            opts.slices_per_container = Some(parse_count(key, value)?);
            Ok(())
        }
        // Accepted for command-line parity. The current CRAM writer uses
        // noodles' default version/compression settings.
        "version" | "level" => Ok(()),
        "" => Err(ParseError::Err(
            "missing value for --output-fmt-option".into(),
        )),
        _ => Err(ParseError::Err(format!(
            "unsupported --output-fmt-option \"{}\"",
            raw
        ))),
    }
}

fn parse_subsample(raw: &str) -> Result<Subsample, ParseError> {
    let (seed, fraction_text) = match raw.split_once('.') {
        Some((seed, frac)) => {
            let seed = if seed.is_empty() {
                0
            } else {
                seed.parse::<u32>()
                    .map_err(|_| ParseError::Err(format!("invalid -s seed in \"{}\"", raw)))?
            };
            (seed, format!("0.{frac}"))
        }
        None => (0, raw.to_string()),
    };
    let fraction = fraction_text
        .parse::<f64>()
        .map_err(|_| ParseError::Err(format!("invalid -s fraction \"{}\"", raw)))?;
    if !(0.0..=1.0).contains(&fraction) {
        return Err(ParseError::Err(format!(
            "invalid -s fraction \"{}\": expected 0..1",
            raw
        )));
    }
    Ok(Subsample { seed, fraction })
}

/// Temp directory holding a `-X` data file plus its explicit index under the
/// default name; removed recursively on drop.
struct StagedIndexDir {
    dir: PathBuf,
    bam: PathBuf,
}

impl StagedIndexDir {
    fn bam_path(&self) -> &Path {
        &self.bam
    }
}

impl Drop for StagedIndexDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.dir);
    }
}

/// Picks the index suffix from the provided file: sniff the BAI/CSI magic,
/// then fall back to the path extension, then to `bai`.
fn customized_index_extension(index: &Path) -> &'static str {
    let mut magic = [0u8; 4];
    if let Ok(mut f) = File::open(index)
        && f.read_exact(&mut magic).is_ok()
    {
        match &magic {
            b"CSI\x01" => return "csi",
            b"BAI\x01" => return "bai",
            _ => {}
        }
    }
    match index.extension().and_then(|e| e.to_str()) {
        Some("csi") => "csi",
        Some("crai") => "crai",
        _ => "bai",
    }
}

/// `-X`/`--customized-index`: the region helpers otherwise resolve the index
/// next to the data file, ignoring the explicit path. Stage a temp directory
/// holding the data file plus the provided index under the default-probe name
/// so those lookups succeed without polluting the source tree.
fn stage_customized_index(input: &Path, index: &Path) -> io::Result<StagedIndexDir> {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let dir_name = format!("samtools-rs-xidx-{}-{}", std::process::id(), nanos);
    let dir = std::env::temp_dir().join(dir_name);
    fs::create_dir_all(&dir)?;

    let file_name = input
        .file_name()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "input has no file name"))?;
    let staged_bam = dir.join(file_name);
    // Hard-link when possible (same filesystem); fall back to a copy.
    if fs::hard_link(input, &staged_bam).is_err() {
        fs::copy(input, &staged_bam)?;
    }

    let mut staged_index = staged_bam.clone().into_os_string();
    staged_index.push(".");
    staged_index.push(customized_index_extension(index));
    fs::copy(index, PathBuf::from(staged_index))?;

    Ok(StagedIndexDir {
        dir,
        bam: staged_bam,
    })
}

fn run(opts: &Opts, input: &Path, input_exact: Exact) -> io::Result<ExitCode> {
    let staged_index_dir = if opts.customized_index
        && let Some(index) = opts.index_path.as_deref()
        && index.exists()
    {
        Some(stage_customized_index(input, index)?)
    } else {
        None
    };
    let input: &Path = staged_index_dir
        .as_ref()
        .map(StagedIndexDir::bam_path)
        .unwrap_or(input);

    let effective_out_fmt = resolved_output_fmt(opts)?;

    // Count-only mode.
    if opts.count {
        reject_unselected_for_count(opts)?;
        let n = if has_requested_regions(opts) {
            count_region_records(input, input_exact, opts)?
        } else if let Some(expr) = opts.filter_expr.as_deref() {
            let expr = combined_filter_expr(opts).unwrap_or_else(|| expr.to_string());
            count_expr_records(input, input_exact, opts, &expr)?
        } else if has_filters(opts) {
            count_filtered_records(input, input_exact, opts)?
        } else {
            count_records(input, input_exact, opts)?
        };
        let mut out = open_text_output(opts)?;
        writeln!(out, "{}", n)?;
        sam_io::check_sam_close(&mut out)?;
        return Ok(ExitCode::SUCCESS);
    }

    // Header-only mode.
    if opts.header == HeaderMode::HeaderOnly {
        let header_text =
            output_header_text(&read_raw_header_text_with_format(input, input_exact)?, opts)?;
        let mut out = open_text_output(opts)?;
        out.write_all(header_text.as_bytes())?;
        sam_io::check_sam_close(&mut out)?;
        return Ok(ExitCode::SUCCESS);
    }

    if opts.fetch_pairs {
        return run_fetch_pairs(opts, input, input_exact, effective_out_fmt);
    }

    // SAM output paths.
    if effective_out_fmt == OutputFmt::Sam {
        validate_unselected_sam_output(opts)?;
        let mut out = open_text_output(opts)?;
        let mut unselected = open_unselected_text_output(opts)?;
        // Upstream default for SAM output is records-only; `-h` opts in
        // to including the header.
        let include_header = match opts.header {
            HeaderMode::Include => true,
            HeaderMode::Default => false,
            HeaderMode::HeaderOnly => true, // handled above
        };
        if include_header {
            let header_text =
                output_header_text(&read_raw_header_text_with_format(input, input_exact)?, opts)?;
            out.write_all(header_text.as_bytes())?;
            if let Some(unselected) = unselected.as_mut() {
                unselected.write_all(header_text.as_bytes())?;
            }
        }
        write_records_as_sam(&mut out, &mut unselected, input, input_exact, opts)?;
        sam_io::check_sam_close(&mut out)?;
        if let Some(unselected) = unselected.as_mut() {
            sam_io::check_sam_close(unselected)?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    // BAM output.
    if effective_out_fmt == OutputFmt::Bam {
        let filter = combined_filter_expr(opts);
        let dst_file = open_binary_output(opts)?;
        if opts.unselected_output.is_some() || opts.unmap_unselected {
            validate_unselected_sam_output(opts)?;
            let text = sam_text_for_binary_split(input, input_exact, opts)?;
            let (selected, unselected) = build_split_sam_text(&text, opts)?;
            htslib_rs::alignment_compat::write_bam_from_sam_reader(
                BufReader::new(io::Cursor::new(sam_bytes_with_pg(&selected, opts)?)),
                dst_file,
            )?;
            if let Some(unselected_path) = opts.unselected_output.as_deref() {
                let unselected_dst = File::create(unselected_path)?;
                htslib_rs::alignment_compat::write_bam_from_sam_reader(
                    BufReader::new(io::Cursor::new(sam_bytes_with_pg(&unselected, opts)?)),
                    unselected_dst,
                )?;
            }
            return Ok(ExitCode::SUCCESS);
        }
        match input_exact {
            Exact::Sam => {
                if let Some(expr) = filter.as_deref() {
                    if has_record_rewrite(opts) {
                        let text =
                            htslib_rs::alignment_compat::view_sam_text_matching_filter_from_path(
                                input, expr,
                            )?;
                        let filtered = filtered_sam_text(text.as_bytes(), opts)?;
                        htslib_rs::alignment_compat::write_bam_from_sam_reader(
                            BufReader::new(io::Cursor::new(sam_bytes_with_pg(&filtered, opts)?)),
                            dst_file,
                        )?;
                    } else {
                        let text =
                            htslib_rs::alignment_compat::view_sam_text_matching_filter_from_path(
                                input, expr,
                            )?;
                        htslib_rs::alignment_compat::write_bam_from_sam_reader(
                            BufReader::new(io::Cursor::new(sam_bytes_with_pg(
                                text.as_bytes(),
                                opts,
                            )?)),
                            dst_file,
                        )?;
                    }
                } else if has_filters(opts) || has_record_rewrite(opts) {
                    let filtered = filtered_sam_text_from_path(input, opts)?;
                    htslib_rs::alignment_compat::write_bam_from_sam_reader(
                        BufReader::new(io::Cursor::new(sam_bytes_with_pg(&filtered, opts)?)),
                        dst_file,
                    )?;
                } else {
                    let raw = read_sam_path_bytes(input)?;
                    htslib_rs::alignment_compat::write_bam_from_sam_reader(
                        BufReader::new(io::Cursor::new(sam_bytes_with_pg(&raw, opts)?)),
                        dst_file,
                    )?;
                }
            }
            Exact::Bam => {
                if has_record_rewrite(opts) {
                    let text = if opts.regions.is_empty() {
                        if let Some(expr) = filter.as_deref() {
                            htslib_rs::alignment_compat::view_bam_as_sam_text_matching_filter_from_path(
                                input, expr,
                            )?
                        } else {
                            htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(
                                input, None,
                            )?
                        }
                    } else {
                        let regions = parse_region_strings(input, &opts.regions)?;
                        if let Some(expr) = filter.as_deref() {
                            htslib_rs::alignment_compat::view_bam_regions_as_sam_text_matching_filter_from_path(
                                input, &regions, expr,
                            )?
                        } else {
                            htslib_rs::alignment_compat::view_bam_regions_as_sam_text_from_path(
                                input, &regions,
                            )?
                        }
                    };
                    let filtered = filtered_binary_rewrite_text(
                        input,
                        input_exact,
                        text,
                        opts,
                        filter.is_some(),
                    )?;
                    htslib_rs::alignment_compat::write_bam_from_sam_reader(
                        BufReader::new(io::Cursor::new(sam_bytes_with_pg(&filtered, opts)?)),
                        dst_file,
                    )?;
                } else {
                    // Whether a binary @PG must be injected. When not,
                    // keep the fast direct binary-copy paths.
                    let want_pg = !opts.no_pg && opts.argv.is_some();
                    if opts.regions.is_empty() {
                        if let Some(expr) = filter.as_deref() {
                            if want_pg {
                                let text = htslib_rs::alignment_compat::view_bam_as_sam_text_matching_filter_from_path(
                                    input, expr,
                                )?;
                                htslib_rs::alignment_compat::write_bam_from_sam_reader(
                                    BufReader::new(io::Cursor::new(sam_bytes_with_pg(
                                        text.as_bytes(),
                                        opts,
                                    )?)),
                                    dst_file,
                                )?;
                            } else {
                                htslib_rs::alignment_compat::write_bam_matching_filter_from_path(
                                    input, expr, dst_file,
                                )?;
                            }
                        } else if !want_pg {
                            htslib_rs::alignment_compat::write_bam_from_path(input, dst_file)?;
                        } else {
                            // Rewrite only the header text; stream
                            // records unchanged.
                            htslib_rs::alignment_compat::write_bam_from_path_transforming_header(
                                input,
                                dst_file,
                                |header_text| apply_pg_to_header(header_text, opts),
                            )?;
                        }
                    } else {
                        let regions = parse_region_strings(input, &opts.regions)?;
                        if !want_pg {
                            if let Some(expr) = filter.as_deref() {
                                htslib_rs::alignment_compat::write_bam_regions_matching_filter_from_path(
                                    input, &regions, expr, dst_file,
                                )?;
                            } else {
                                htslib_rs::alignment_compat::write_bam_regions_from_path(
                                    input, &regions, dst_file,
                                )?;
                            }
                        } else {
                            let text = if let Some(expr) = filter.as_deref() {
                                htslib_rs::alignment_compat::view_bam_regions_as_sam_text_matching_filter_from_path(
                                    input, &regions, expr,
                                )?
                            } else {
                                htslib_rs::alignment_compat::view_bam_regions_as_sam_text_from_path(
                                    input, &regions,
                                )?
                            };
                            htslib_rs::alignment_compat::write_bam_from_sam_reader(
                                BufReader::new(io::Cursor::new(sam_bytes_with_pg(
                                    text.as_bytes(),
                                    opts,
                                )?)),
                                dst_file,
                            )?;
                        }
                    }
                }
            }
            Exact::Cram => {
                let reference_guard = optional_cram_input_reference_for_path(opts, input)?;
                if reference_guard.is_none() {
                    if !opts.regions.is_empty() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "CRAM region output requires --reference / -T, @SQ UR tags, or REF_PATH entries matching @SQ M5 tags",
                        ));
                    }
                    let text = cram_sam_text_from_path_maybe_synthesizing_reference(
                        input,
                        opts,
                        filter.as_deref(),
                        None,
                    )?;
                    let bytes = if has_record_rewrite(opts) {
                        filtered_binary_rewrite_text(
                            input,
                            input_exact,
                            text,
                            opts,
                            filter.is_some(),
                        )?
                    } else {
                        crate::sam_render::fix_sam_text(&text).into_bytes()
                    };
                    htslib_rs::alignment_compat::write_bam_from_sam_reader(
                        BufReader::new(io::Cursor::new(sam_bytes_with_pg(&bytes, opts)?)),
                        dst_file,
                    )?;
                    return Ok(ExitCode::SUCCESS);
                }

                let reference_guard = reference_guard.unwrap();
                let reference = reference_guard.path();
                if has_record_rewrite(opts) {
                    let text = if opts.regions.is_empty() {
                        if let Some(expr) = filter.as_deref() {
                            htslib_rs::alignment_compat::view_cram_as_sam_text_matching_filter_from_path_with_reference(
                                input, reference, expr,
                            )?
                        } else {
                            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                                input, reference, None,
                            )?
                        }
                    } else {
                        let regions = parse_region_strings(input, &opts.regions)?;
                        if let Some(expr) = filter.as_deref() {
                            htslib_rs::alignment_compat::view_cram_regions_as_sam_text_matching_filter_from_path_with_reference(
                                input, reference, &regions, expr,
                            )?
                        } else {
                            htslib_rs::alignment_compat::view_cram_regions_as_sam_text_from_path_with_reference(
                                input, reference, &regions, false,
                            )?
                        }
                    };
                    let filtered = filtered_binary_rewrite_text(
                        input,
                        input_exact,
                        text,
                        opts,
                        filter.is_some(),
                    )?;
                    htslib_rs::alignment_compat::write_bam_from_sam_reader(
                        BufReader::new(io::Cursor::new(sam_bytes_with_pg(&filtered, opts)?)),
                        dst_file,
                    )?;
                } else if opts.regions.is_empty() {
                    // CRAM-input binary @PG not yet wired: routing
                    // through the SAM-text reader regresses on CRAMs
                    // noodles can't re-decode with an external
                    // reference (the same limitation behind #2/#3).
                    // Keep the faithful direct binary copy.
                    if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::write_cram_records_matching_filter_as_bam_from_path_with_reference(
                            input, reference, expr, dst_file,
                        )?;
                    } else {
                        htslib_rs::alignment_compat::write_cram_records_with_required_flags_as_bam_from_path_with_reference(
                            input, reference, 0, dst_file,
                        )?;
                    }
                } else {
                    let regions = parse_region_strings(input, &opts.regions)?;
                    if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::write_cram_regions_matching_filter_as_bam_from_path_with_reference(
                            input, reference, &regions, expr, dst_file,
                        )?;
                    } else {
                        htslib_rs::alignment_compat::write_cram_regions_as_bam_from_path_with_reference(
                            input, reference, &regions, dst_file,
                        )?;
                    }
                }
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "this input format cannot be written as BAM yet",
                ));
            }
        }
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Cram {
        let filter = combined_filter_expr(opts);
        let reference_guard = cram_output_reference_for_input(opts, input, input_exact)?;
        let reference = reference_guard.path();
        let dst_file = open_binary_output(opts)?;
        if opts.unselected_output.is_some() || opts.unmap_unselected {
            validate_unselected_sam_output(opts)?;
            let text = sam_text_for_binary_split(input, input_exact, opts)?;
            let (selected, unselected) = build_split_sam_text(&text, opts)?;
            write_cram_from_sam_text_via_bam(opts, &selected, reference, dst_file)?;
            if let Some(unselected_path) = opts.unselected_output.as_deref() {
                let unselected_dst = File::create(unselected_path)?;
                write_cram_from_sam_text_via_bam(opts, &unselected, reference, unselected_dst)?;
            }
            return Ok(ExitCode::SUCCESS);
        }
        match input_exact {
            Exact::Sam => {
                if let Some(expr) = filter.as_deref() {
                    let text = if has_record_rewrite(opts) {
                        let t =
                            htslib_rs::alignment_compat::view_sam_text_matching_filter_from_path(
                                input, expr,
                            )?;
                        filtered_sam_text(t.as_bytes(), opts)?
                    } else {
                        htslib_rs::alignment_compat::view_sam_text_matching_filter_from_path(
                            input, expr,
                        )?
                        .into_bytes()
                    };
                    let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(
                        io::Cursor::new(sam_bytes_with_pg(&text, opts)?),
                    ));
                    cram_sam_reader_writer(opts, &mut reader, reference, dst_file)?;
                } else if has_filters(opts) || has_record_rewrite(opts) {
                    let filtered = filtered_sam_text_from_path(input, opts)?;
                    let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(
                        io::Cursor::new(sam_bytes_with_pg(&filtered, opts)?),
                    ));
                    cram_sam_reader_writer(opts, &mut reader, reference, dst_file)?;
                } else {
                    let raw = read_sam_path_bytes(input)?;
                    let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(
                        io::Cursor::new(sam_bytes_with_pg(&raw, opts)?),
                    ));
                    cram_sam_reader_writer(opts, &mut reader, reference, dst_file)?;
                }
            }
            Exact::Bam if opts.regions.is_empty() => {
                if has_record_rewrite(opts) {
                    let text = if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::view_bam_as_sam_text_matching_filter_from_path(
                            input, expr,
                        )?
                    } else {
                        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(
                            input, None,
                        )?
                    };
                    let filtered = filtered_binary_rewrite_text(
                        input,
                        input_exact,
                        text,
                        opts,
                        filter.is_some(),
                    )?;
                    write_cram_from_sam_text_via_bam(opts, &filtered, reference, dst_file)?;
                } else if (!opts.no_pg && opts.argv.is_some())
                    || opts.reference.is_some()
                    || opts.reference_index.is_some()
                {
                    let text = if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::view_bam_as_sam_text_matching_filter_from_path(
                            input, expr,
                        )?
                    } else {
                        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(
                            input, None,
                        )?
                    };
                    let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(
                        io::Cursor::new(sam_bytes_with_pg(text.as_bytes(), opts)?),
                    ));
                    cram_sam_reader_writer(opts, &mut reader, reference, dst_file)?;
                } else if let Some(expr) = filter.as_deref() {
                    htslib_rs::alignment_compat::write_cram_matching_filter_from_bam_path_with_reference(
                        input,
                        reference,
                        expr,
                        cram_write_options(opts),
                        dst_file,
                    )?;
                } else {
                    htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference_and_options(
                        input,
                        reference,
                        cram_write_options(opts),
                        dst_file,
                    )?;
                }
            }
            Exact::Bam => {
                let regions = parse_region_strings(input, &opts.regions)?;
                if has_record_rewrite(opts) {
                    let text = if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::view_bam_regions_as_sam_text_matching_filter_from_path(
                            input, &regions, expr,
                        )?
                    } else {
                        htslib_rs::alignment_compat::view_bam_regions_as_sam_text_from_path(
                            input, &regions,
                        )?
                    };
                    let filtered = filtered_binary_rewrite_text(
                        input,
                        input_exact,
                        text,
                        opts,
                        filter.is_some(),
                    )?;
                    write_cram_from_sam_text_via_bam(opts, &filtered, reference, dst_file)?;
                } else if !opts.no_pg && opts.argv.is_some() {
                    let text = if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::view_bam_regions_as_sam_text_matching_filter_from_path(
                            input, &regions, expr,
                        )?
                    } else {
                        htslib_rs::alignment_compat::view_bam_regions_as_sam_text_from_path(
                            input, &regions,
                        )?
                    };
                    let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(
                        io::Cursor::new(sam_bytes_with_pg(text.as_bytes(), opts)?),
                    ));
                    cram_sam_reader_writer(opts, &mut reader, reference, dst_file)?;
                } else if let Some(expr) = filter.as_deref() {
                    htslib_rs::alignment_compat::write_bam_regions_matching_filter_as_cram_from_path_with_reference(
                        input, reference, &regions, expr, dst_file,
                    )?;
                } else {
                    htslib_rs::alignment_compat::write_bam_regions_as_cram_from_path_with_reference(
                        input, reference, &regions, dst_file,
                    )?;
                }
            }
            Exact::Cram if opts.regions.is_empty() => {
                if has_record_rewrite(opts) {
                    let text = if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::view_cram_as_sam_text_matching_filter_from_path_with_reference(
                            input, reference, expr,
                        )?
                    } else {
                        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                            input, reference, None,
                        )?
                    };
                    let filtered = filtered_binary_rewrite_text(
                        input,
                        input_exact,
                        text,
                        opts,
                        filter.is_some(),
                    )?;
                    write_cram_from_sam_text_via_bam(opts, &filtered, reference, dst_file)?;
                } else if let Some(expr) = filter.as_deref() {
                    htslib_rs::alignment_compat::write_cram_matching_filter_from_path_with_reference(
                        input,
                        reference,
                        expr,
                        cram_write_options(opts),
                        dst_file,
                    )?;
                } else {
                    htslib_rs::alignment_compat::write_cram_from_path_with_reference_and_options(
                        input,
                        reference,
                        cram_write_options(opts),
                        dst_file,
                    )?;
                }
            }
            Exact::Cram => {
                let regions = parse_region_strings(input, &opts.regions)?;
                if has_record_rewrite(opts) {
                    let text = if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::view_cram_regions_as_sam_text_matching_filter_from_path_with_reference(
                            input, reference, &regions, expr,
                        )?
                    } else {
                        htslib_rs::alignment_compat::view_cram_regions_as_sam_text_from_path_with_reference(
                            input, reference, &regions, false,
                        )?
                    };
                    let filtered = filtered_binary_rewrite_text(
                        input,
                        input_exact,
                        text,
                        opts,
                        filter.is_some(),
                    )?;
                    write_cram_from_sam_text_via_bam(opts, &filtered, reference, dst_file)?;
                } else if let Some(expr) = filter.as_deref() {
                    htslib_rs::alignment_compat::write_cram_regions_matching_filter_from_path_with_reference(
                        input,
                        reference,
                        &regions,
                        expr,
                        cram_write_options(opts),
                        dst_file,
                    )?;
                } else {
                    htslib_rs::alignment_compat::write_cram_regions_from_path_with_reference(
                        input,
                        reference,
                        &regions,
                        cram_write_options(opts),
                        dst_file,
                    )?;
                }
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "this input/output CRAM conversion is not yet wired up",
                ));
            }
        }
        return Ok(ExitCode::SUCCESS);
    }

    Err(io::Error::other("unsupported output combination"))
}

fn run_fetch_pairs(
    opts: &Opts,
    input: &Path,
    exact: Exact,
    output_fmt: OutputFmt,
) -> io::Result<ExitCode> {
    let (sam_text, _, _) = fetch_pairs_sam_text(input, exact, opts)?;
    match output_fmt {
        OutputFmt::Sam => {
            let mut out = open_text_output(opts)?;
            out.write_all(sam_text.as_bytes())?;
            sam_io::check_sam_close(&mut out)?;
        }
        OutputFmt::Bam => {
            let dst_file = open_binary_output(opts)?;
            htslib_rs::alignment_compat::write_bam_from_sam_reader(
                BufReader::new(io::Cursor::new(sam_bytes_with_pg(
                    sam_text.as_bytes(),
                    opts,
                )?)),
                dst_file,
            )?;
        }
        OutputFmt::Cram => {
            let reference_guard = cram_output_reference_for_input(opts, input, exact)?;
            let reference = reference_guard.path();
            let dst_file = open_binary_output(opts)?;
            let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(
                sam_bytes_with_pg(sam_text.as_bytes(), opts)?,
            )));
            cram_sam_reader_writer(opts, &mut reader, reference, dst_file)?;
        }
        OutputFmt::Auto => unreachable!("resolved output format is never Auto"),
    }
    Ok(ExitCode::SUCCESS)
}

fn fetch_pairs_sam_text(
    input: &Path,
    exact: Exact,
    opts: &Opts,
) -> io::Result<(String, usize, usize)> {
    if opts.regions.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--fetch-pairs requires at least one region",
        ));
    }

    let text = match exact {
        Exact::Sam => String::from_utf8(read_sam_path_bytes(input)?)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        Exact::Bam => {
            htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(input, None)?
        }
        Exact::Cram => {
            let reference_guard = cram_input_reference_for_path(opts, input)?;
            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                input,
                reference_guard.path(),
                None,
            )?
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "--fetch-pairs is only supported for SAM/BAM/CRAM input",
            ));
        }
    };

    let header_len = sam_header_lines(text.as_bytes()).len();
    let header = &text[..header_len];
    let regions = parse_simple_regions(&opts.regions)?;
    let mut filter_opts = opts.clone();
    filter_opts.regions.clear();
    filter_opts.fetch_pairs = false;

    let mut qnames = HashSet::new();
    for line in text[header_len..].lines().filter(|line| !line.is_empty()) {
        let bytes = line.as_bytes();
        if sam_line_overlaps_regions(bytes, &regions)
            && let Some(qname) = bytes.split(|&b| b == b'\t').next()
        {
            qnames.insert(qname.to_vec());
        }
    }

    let mut out = String::with_capacity(text.len());
    out.push_str(&output_header_text(header, opts)?);
    let mut processed = 0usize;
    let mut accepted = 0usize;
    for line in text[header_len..].lines().filter(|line| !line.is_empty()) {
        let bytes = line.as_bytes();
        let Some(qname) = bytes.split(|&b| b == b'\t').next() else {
            continue;
        };
        if !qnames.contains(qname) {
            continue;
        }
        if line_passes(bytes, &filter_opts) {
            processed += 1;
            accepted += 1;
            let rendered = render_sam_record_line(bytes, opts);
            out.push_str(
                std::str::from_utf8(&rendered)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
            );
            out.push('\n');
        } else if sam_line_overlaps_regions(bytes, &regions) {
            processed += 1;
        }
    }

    Ok((out, processed, accepted))
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum StdinFormat {
    Sam,
    Bam,
    Cram,
}

fn stdin_format(input: &[u8]) -> StdinFormat {
    if input.starts_with(b"CRAM") {
        StdinFormat::Cram
    } else if input.starts_with(b"\x1f\x8b") || input.starts_with(b"BAM\x01") {
        StdinFormat::Bam
    } else {
        StdinFormat::Sam
    }
}

fn run_sam_stdin(opts: &Opts, input: &[u8]) -> io::Result<ExitCode> {
    if !opts.regions.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "region queries on stdin require an index (BAM/CRAM file input only)",
        ));
    }

    let effective_out_fmt = resolved_output_fmt(opts)?;

    if opts.count {
        reject_unselected_for_count(opts)?;
        let filter = combined_filter_expr(opts);
        let n = if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::count_sam_records_matching_filter(
                BufReader::new(io::Cursor::new(input)),
                expr,
            )?
        } else {
            count_sam_text_records(input, opts)
        };
        let mut out = open_text_output(opts)?;
        writeln!(out, "{}", n)?;
        sam_io::check_sam_close(&mut out)?;
        return Ok(ExitCode::SUCCESS);
    }

    if opts.header == HeaderMode::HeaderOnly {
        let mut out = open_text_output(opts)?;
        let header_text = output_header_text(
            std::str::from_utf8(sam_header_lines(input)).unwrap_or(""),
            opts,
        )?;
        out.write_all(header_text.as_bytes())?;
        sam_io::check_sam_close(&mut out)?;
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Sam {
        validate_unselected_sam_output(opts)?;
        let mut out = open_text_output(opts)?;
        let mut unselected = open_unselected_text_output(opts)?;
        if opts.header == HeaderMode::Include {
            let header_text = output_header_text(
                std::str::from_utf8(sam_header_lines(input)).unwrap_or(""),
                opts,
            )?;
            out.write_all(header_text.as_bytes())?;
            if let Some(unselected) = unselected.as_mut() {
                unselected.write_all(header_text.as_bytes())?;
            }
        }
        let filter = prefilter_expr_for_sam_output(opts);
        if let Some(expr) = filter.as_deref() {
            let text = htslib_rs::alignment_compat::view_sam_text_matching_filter(
                BufReader::new(io::Cursor::new(input)),
                expr,
            )?;
            let text = crate::sam_render::fix_sam_text(&text);
            let post_filter_opts = opts_after_prefiltered_expr(opts);
            write_sam_text_records_split(
                &mut out,
                &mut unselected,
                text.as_bytes(),
                &post_filter_opts,
            )?;
        } else {
            write_sam_text_records_split(&mut out, &mut unselected, input, opts)?;
        }
        sam_io::check_sam_close(&mut out)?;
        if let Some(unselected) = unselected.as_mut() {
            sam_io::check_sam_close(unselected)?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Bam {
        reject_unselected_binary_output(opts)?;
        let dst_file = open_binary_output(opts)?;
        let filter = combined_filter_expr(opts);
        if let Some(expr) = filter.as_deref() {
            if has_record_rewrite(opts) {
                let text = htslib_rs::alignment_compat::view_sam_text_matching_filter(
                    BufReader::new(io::Cursor::new(input)),
                    expr,
                )?;
                let filtered = filtered_sam_text(text.as_bytes(), opts)?;
                htslib_rs::alignment_compat::write_bam_from_sam_reader(
                    BufReader::new(io::Cursor::new(sam_bytes_with_pg(&filtered, opts)?)),
                    dst_file,
                )?;
            } else {
                let mut reader =
                    htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(input)));
                htslib_rs::alignment_compat::write_bam_matching_filter_from_sam_reader(
                    &mut reader,
                    expr,
                    dst_file,
                )?;
            }
        } else {
            let reader_input = if has_filters(opts) || has_record_rewrite(opts) {
                filtered_sam_text(input, opts)?
            } else {
                input.to_vec()
            };
            htslib_rs::alignment_compat::write_bam_from_sam_reader(
                BufReader::new(io::Cursor::new(sam_bytes_with_pg(&reader_input, opts)?)),
                dst_file,
            )?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Cram {
        reject_unselected_binary_output(opts)?;
        let reference = cram_reference(opts)?;
        let dst_file = open_binary_output(opts)?;
        let filter = combined_filter_expr(opts);
        if let Some(expr) = filter.as_deref() {
            if has_record_rewrite(opts) {
                let text = htslib_rs::alignment_compat::view_sam_text_matching_filter(
                    BufReader::new(io::Cursor::new(input)),
                    expr,
                )?;
                let filtered = filtered_sam_text(text.as_bytes(), opts)?;
                let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(
                    sam_bytes_with_pg(&filtered, opts)?,
                )));
                cram_sam_reader_writer(opts, &mut reader, reference, dst_file)?;
            } else {
                let mut reader =
                    htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(input)));
                htslib_rs::alignment_compat::write_cram_matching_filter_from_sam_reader_with_reference(
                    &mut reader,
                    reference,
                    expr,
                    cram_write_options(opts),
                    dst_file,
                )?;
            }
        } else {
            let reader_input = if has_filters(opts) || has_record_rewrite(opts) {
                filtered_sam_text(input, opts)?
            } else {
                input.to_vec()
            };
            let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(
                sam_bytes_with_pg(&reader_input, opts)?,
            )));
            cram_sam_reader_writer(opts, &mut reader, reference, dst_file)?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    Err(io::Error::other("unsupported stdin output combination"))
}

fn run_bam_stdin(opts: &Opts, input: &[u8]) -> io::Result<ExitCode> {
    if !opts.regions.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "region queries on stdin require an index (BAM/CRAM file input only)",
        ));
    }

    let effective_out_fmt = resolved_output_fmt(opts)?;

    if opts.count {
        reject_unselected_for_count(opts)?;
        let filter = combined_filter_expr(opts);
        let n = if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::count_bam_records_matching_filter(
                io::Cursor::new(input),
                expr,
            )?
        } else {
            let text =
                htslib_rs::alignment_compat::view_bam_as_sam_text(io::Cursor::new(input), None)?;
            count_sam_text_records(text.as_bytes(), opts)
        };
        let mut out = open_text_output(opts)?;
        writeln!(out, "{}", n)?;
        sam_io::check_sam_close(&mut out)?;
        return Ok(ExitCode::SUCCESS);
    }

    if opts.header == HeaderMode::HeaderOnly {
        let text =
            htslib_rs::alignment_compat::view_bam_as_sam_text(io::Cursor::new(input), Some(0))?;
        let header_text =
            output_header_text(std::str::from_utf8(text.as_bytes()).unwrap_or(""), opts)?;
        let mut out = open_text_output(opts)?;
        out.write_all(header_text.as_bytes())?;
        sam_io::check_sam_close(&mut out)?;
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Sam {
        validate_unselected_sam_output(opts)?;
        let filter = prefilter_expr_for_sam_output(opts);
        let text = if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::view_bam_as_sam_text_matching_filter(
                io::Cursor::new(input),
                expr,
            )?
        } else {
            htslib_rs::alignment_compat::view_bam_as_sam_text(io::Cursor::new(input), None)?
        };
        // Records came from binary; fix noodles' plain-decimal float
        // spelling to htslib's `%g` form (header lines pass through).
        let text = crate::sam_render::fix_sam_text(&text);
        let mut out = open_text_output(opts)?;
        let mut unselected = open_unselected_text_output(opts)?;
        if opts.header == HeaderMode::Include {
            let header_text = output_header_text(
                std::str::from_utf8(sam_header_lines(text.as_bytes())).unwrap_or(""),
                opts,
            )?;
            out.write_all(header_text.as_bytes())?;
            if let Some(unselected) = unselected.as_mut() {
                unselected.write_all(header_text.as_bytes())?;
            }
        }
        let post_filter_opts;
        let split_opts = if filter.is_some() {
            post_filter_opts = opts_after_prefiltered_expr(opts);
            &post_filter_opts
        } else {
            opts
        };
        write_sam_text_records_split(&mut out, &mut unselected, text.as_bytes(), split_opts)?;
        sam_io::check_sam_close(&mut out)?;
        if let Some(unselected) = unselected.as_mut() {
            sam_io::check_sam_close(unselected)?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Bam {
        reject_unselected_binary_output(opts)?;
        let filter = combined_filter_expr(opts);
        let dst_file = open_binary_output(opts)?;
        if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::write_bam_matching_filter(
                io::Cursor::new(input),
                expr,
                dst_file,
            )?;
        } else {
            htslib_rs::alignment_compat::write_bam(io::Cursor::new(input), dst_file)?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Cram {
        reject_unselected_binary_output(opts)?;
        let filter = combined_filter_expr(opts);
        let reference = cram_reference(opts)?;
        let dst_file = open_binary_output(opts)?;
        if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::write_cram_matching_filter_from_bam_reader_with_reference(
                io::Cursor::new(input),
                &reference,
                expr,
                cram_write_options(opts),
                dst_file,
            )?;
        } else if opts.reference.is_some() || opts.reference_index.is_some() {
            let text =
                htslib_rs::alignment_compat::view_bam_as_sam_text(io::Cursor::new(input), None)?;
            let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(
                sam_bytes_with_pg(text.as_bytes(), opts)?,
            )));
            cram_sam_reader_writer(opts, &mut reader, reference, dst_file)?;
        } else {
            htslib_rs::alignment_compat::write_cram_from_bam_reader_with_reference(
                io::Cursor::new(input),
                reference,
                dst_file,
            )?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    Err(io::Error::other("unsupported stdin output combination"))
}

fn run_cram_stdin(opts: &Opts, input: &[u8]) -> io::Result<ExitCode> {
    if !opts.regions.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "region queries on stdin require an index (BAM/CRAM file input only)",
        ));
    }

    let reference_guard = cram_stdin_reference(opts, input)?;
    let reference = reference_guard.path();
    let effective_out_fmt = resolved_output_fmt(opts)?;

    if opts.count {
        reject_unselected_for_count(opts)?;
        let filter = combined_filter_expr(opts);
        let n = if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::count_cram_records_matching_filter_with_reference(
                io::Cursor::new(input),
                reference,
                expr,
            )?
        } else {
            let text = htslib_rs::alignment_compat::view_cram_as_sam_text_with_reference(
                io::Cursor::new(input),
                reference,
                None,
            )?;
            count_sam_text_records(text.as_bytes(), opts)
        };
        let mut out = open_text_output(opts)?;
        writeln!(out, "{}", n)?;
        sam_io::check_sam_close(&mut out)?;
        return Ok(ExitCode::SUCCESS);
    }

    if opts.header == HeaderMode::HeaderOnly {
        let text = htslib_rs::alignment_compat::view_cram_as_sam_text_with_reference(
            io::Cursor::new(input),
            reference,
            Some(0),
        )?;
        let header_text =
            output_header_text(std::str::from_utf8(text.as_bytes()).unwrap_or(""), opts)?;
        let mut out = open_text_output(opts)?;
        out.write_all(header_text.as_bytes())?;
        sam_io::check_sam_close(&mut out)?;
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Sam {
        validate_unselected_sam_output(opts)?;
        let filter = prefilter_expr_for_sam_output(opts);
        let mut out = open_text_output(opts)?;
        let mut unselected = open_unselected_text_output(opts)?;
        if opts.header == HeaderMode::Include {
            let header = htslib_rs::alignment_compat::read_cram_header(io::Cursor::new(input))?;
            let mut header_text = Vec::new();
            crate::sam_render::write_header(&mut header_text, &header)?;
            let header_text = String::from_utf8(header_text)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let header_text = output_header_text(&header_text, opts)?;
            out.write_all(header_text.as_bytes())?;
            if let Some(unselected) = unselected.as_mut() {
                unselected.write_all(header_text.as_bytes())?;
            }
        }
        let text = if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::view_cram_as_sam_text_matching_filter_with_reference(
                io::Cursor::new(input),
                reference,
                expr,
            )?
        } else {
            htslib_rs::alignment_compat::view_cram_as_sam_text_with_reference(
                io::Cursor::new(input),
                reference,
                None,
            )?
        };
        let text = crate::sam_render::fix_sam_text(&text);
        let post_filter_opts;
        let split_opts = if filter.is_some() {
            post_filter_opts = opts_after_prefiltered_expr(opts);
            &post_filter_opts
        } else {
            opts
        };
        write_sam_text_records_split(&mut out, &mut unselected, text.as_bytes(), split_opts)?;
        sam_io::check_sam_close(&mut out)?;
        if let Some(unselected) = unselected.as_mut() {
            sam_io::check_sam_close(unselected)?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Bam {
        reject_unselected_binary_output(opts)?;
        let filter = combined_filter_expr(opts);
        let dst_file = open_binary_output(opts)?;
        if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::write_cram_records_matching_filter_as_bam_with_reference(
                io::Cursor::new(input),
                reference,
                expr,
                dst_file,
            )?;
        } else {
            htslib_rs::alignment_compat::write_cram_records_with_required_flags_as_bam_with_reference(
                io::Cursor::new(input),
                reference,
                0,
                dst_file,
            )?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Cram {
        reject_unselected_binary_output(opts)?;
        let filter = combined_filter_expr(opts);
        let dst_file = open_binary_output(opts)?;
        if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::write_cram_matching_filter_from_reader_with_reference(
                io::Cursor::new(input),
                reference,
                expr,
                cram_write_options(opts),
                dst_file,
            )?;
        } else {
            htslib_rs::alignment_compat::write_cram_from_reader_with_reference(
                io::Cursor::new(input),
                reference,
                dst_file,
            )?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    Err(io::Error::other("unsupported stdin output combination"))
}

fn resolved_output_fmt(opts: &Opts) -> io::Result<OutputFmt> {
    let explicit = match opts.output_fmt {
        OutputFmt::Auto => None,
        OutputFmt::Sam => Some(Exact::Sam),
        OutputFmt::Bam => Some(Exact::Bam),
        OutputFmt::Cram => Some(Exact::Cram),
    };
    let mode = sam_io::sam_open_mode(opts.output.as_deref(), explicit, Exact::Sam)?;
    match mode.exact {
        Exact::Sam => Ok(OutputFmt::Sam),
        Exact::Bam => Ok(OutputFmt::Bam),
        Exact::Cram => Ok(OutputFmt::Cram),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("unsupported output format inferred from {:?}", opts.output),
        )),
    }
}

/// Applies samtools' standard `@PG` chain entry to a SAM header text,
/// unless `--no-PG` was passed or the caller hasn't supplied argv.
fn apply_pg_to_header(header_text: &str, opts: &Opts) -> io::Result<String> {
    if opts.no_pg {
        return Ok(header_text.to_owned());
    }
    let Some(argv) = opts.argv.as_deref() else {
        return Ok(header_text.to_owned());
    };
    crate::pg::add_samtools_pg(header_text, argv)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn output_header_text(header_text: &str, opts: &Opts) -> io::Result<String> {
    let header_text = inject_reference_dictionary_if_missing(header_text, opts, false)?;
    let header_text = apply_pg_to_header(&header_text, opts)?;
    Ok(prune_rg_header_lines(&header_text, opts))
}

fn inject_reference_dictionary_if_missing(
    header_text: &str,
    opts: &Opts,
    include_uri: bool,
) -> io::Result<String> {
    let has_sq = header_text
        .lines()
        .any(|line| line.starts_with("@SQ\t") || line == "@SQ");
    if has_sq && include_uri {
        return augment_sq_lines_with_reference_metadata(header_text, opts);
    }
    if has_sq {
        return Ok(header_text.to_owned());
    }

    let Some(dictionary) = reference_dictionary_lines(opts, include_uri)? else {
        return Ok(header_text.to_owned());
    };
    if dictionary.is_empty() {
        return Ok(header_text.to_owned());
    }

    let mut out = String::with_capacity(header_text.len() + dictionary.len() + 1);
    out.push_str(header_text);
    if !out.is_empty() && !out.ends_with('\n') {
        out.push('\n');
    }
    out.push_str(&dictionary);
    Ok(out)
}

fn reference_dictionary_lines(opts: &Opts, include_uri: bool) -> io::Result<Option<String>> {
    if let Some(index) = opts.reference_index.as_deref() {
        let uri = include_uri
            .then(|| fasta_path_from_fai_path(index))
            .flatten();
        if include_uri
            && let Some(uri) = uri.as_deref()
            && uri.is_file()
        {
            return read_fasta_dictionary(uri, Some(uri), true).map(Some);
        }
        return read_fai_dictionary(index, uri.as_deref()).map(Some);
    }

    let Some(reference) = explicit_reference_path(opts) else {
        return Ok(None);
    };
    if include_uri {
        return read_fasta_dictionary(&reference, Some(&reference), true).map(Some);
    }

    let fai = reference_fai_path(&reference);
    if fai.is_file() {
        return read_fai_dictionary(&fai, None).map(Some);
    }

    read_fasta_dictionary(&reference, None, false).map(Some)
}

fn read_fai_dictionary(path: &Path, uri: Option<&Path>) -> io::Result<String> {
    let text = fs::read_to_string(path)?;
    let mut out = String::new();
    for (line_no, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let name = fields.next().unwrap_or_default();
        let len = fields.next().unwrap_or_default();
        if name.is_empty() || len.parse::<u64>().is_err() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid FASTA index line {} in {}",
                    line_no + 1,
                    path.display()
                ),
            ));
        }
        out.push_str("@SQ\tSN:");
        out.push_str(name);
        out.push_str("\tLN:");
        out.push_str(len);
        if let Some(uri) = uri {
            out.push_str("\tUR:");
            out.push_str(&uri.display().to_string());
        }
        out.push('\n');
    }
    Ok(out)
}

fn read_fasta_dictionary(path: &Path, uri: Option<&Path>, include_md5: bool) -> io::Result<String> {
    let file = File::open(path)?;
    let mut reader: Box<dyn BufRead> = if is_bgzf_path(path)? {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut out = String::new();
    let mut line = String::new();
    let mut name: Option<String> = None;
    let mut len = 0usize;
    let mut sequence = Vec::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(name) = name.take() {
                let md5 = include_md5.then(|| md5_hex(&sequence));
                push_sq_line(&mut out, &name, len, md5.as_deref(), uri);
            }
            name = rest.split_ascii_whitespace().next().map(str::to_owned);
            len = 0;
            sequence.clear();
        } else {
            for b in trimmed.bytes().filter(|b| !b.is_ascii_whitespace()) {
                len += 1;
                if include_md5 {
                    sequence.push(b.to_ascii_uppercase());
                }
            }
        }
    }
    if let Some(name) = name {
        let md5 = include_md5.then(|| md5_hex(&sequence));
        push_sq_line(&mut out, &name, len, md5.as_deref(), uri);
    }
    Ok(out)
}

fn push_sq_line(out: &mut String, name: &str, len: usize, md5: Option<&str>, uri: Option<&Path>) {
    if name.is_empty() {
        return;
    }
    out.push_str("@SQ\tSN:");
    out.push_str(name);
    out.push_str("\tLN:");
    out.push_str(&len.to_string());
    if let Some(md5) = md5 {
        out.push_str("\tM5:");
        out.push_str(md5);
    }
    if let Some(uri) = uri {
        out.push_str("\tUR:");
        out.push_str(&uri.display().to_string());
    }
    out.push('\n');
}

fn augment_sq_lines_with_reference_metadata(header_text: &str, opts: &Opts) -> io::Result<String> {
    let Some(reference) = explicit_reference_path(opts) else {
        return Ok(header_text.to_owned());
    };
    let md5s = fasta_md5_by_name(&reference)?;
    let mut out = String::with_capacity(header_text.len() + 96);
    for line in header_text.split_inclusive('\n') {
        let (body, newline) = match line.strip_suffix('\n') {
            Some(body) => (body, "\n"),
            None => (line, ""),
        };
        if !body.starts_with("@SQ\t") {
            out.push_str(line);
            continue;
        }

        let name = body
            .split('\t')
            .skip(1)
            .find_map(|field| field.strip_prefix("SN:"));
        let has_m5 = body
            .split('\t')
            .skip(1)
            .any(|field| field.starts_with("M5:"));
        let has_uri = body
            .split('\t')
            .skip(1)
            .any(|field| field.starts_with("UR:"));

        out.push_str(body);
        if !has_m5
            && let Some(name) = name
            && let Some(md5) = md5s.get(name)
        {
            out.push_str("\tM5:");
            out.push_str(md5);
        }
        if !has_uri {
            out.push_str("\tUR:");
            out.push_str(&reference.display().to_string());
        }
        out.push_str(newline);
    }
    Ok(out)
}

fn fasta_md5_by_name(path: &Path) -> io::Result<std::collections::HashMap<String, String>> {
    let file = File::open(path)?;
    let mut reader: Box<dyn BufRead> = if is_bgzf_path(path)? {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut map = std::collections::HashMap::new();
    let mut line = String::new();
    let mut name: Option<String> = None;
    let mut sequence = Vec::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(name) = name.take() {
                map.insert(name, md5_hex(&sequence));
            }
            name = rest.split_ascii_whitespace().next().map(str::to_owned);
            sequence.clear();
        } else {
            sequence.extend(
                trimmed
                    .bytes()
                    .filter(|b| !b.is_ascii_whitespace())
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
    }
    if let Some(name) = name {
        map.insert(name, md5_hex(&sequence));
    }
    Ok(map)
}

fn md5_hex(sequence: &[u8]) -> String {
    let mut hasher = Md5::new();
    hasher.update(sequence);
    format!("{:x}", hasher.finalize())
}

fn explicit_reference_path(opts: &Opts) -> Option<PathBuf> {
    opts.reference
        .clone()
        .or_else(|| current_global_args().reference)
        .or_else(|| {
            opts.reference_index
                .as_deref()
                .and_then(fasta_path_from_fai_path)
                .filter(|path| path.is_file())
        })
}

fn reference_fai_path(reference: &Path) -> PathBuf {
    let mut path = reference.as_os_str().to_os_string();
    path.push(".fai");
    PathBuf::from(path)
}

fn fasta_path_from_fai_path(path: &Path) -> Option<PathBuf> {
    let text = path.as_os_str().to_str()?;
    text.strip_suffix(".fai").map(PathBuf::from)
}

fn prune_rg_header_lines(header_text: &str, opts: &Opts) -> String {
    let Some(keep_ids) = output_header_rg_ids(opts) else {
        return header_text.to_owned();
    };

    let mut out = String::with_capacity(header_text.len());
    for line in header_text.split_inclusive('\n') {
        let body = line.strip_suffix('\n').unwrap_or(line);
        if body.starts_with("@RG\t") {
            if header_rg_id(body)
                .map(|id| keep_ids.contains(id.as_bytes()))
                .unwrap_or(false)
            {
                out.push_str(line);
            }
        } else {
            out.push_str(line);
        }
    }
    out
}

fn output_header_rg_ids(opts: &Opts) -> Option<HashSet<Vec<u8>>> {
    if opts.read_groups.is_empty() {
        return None;
    }

    let mut keep_ids = opts.read_groups.clone();
    if opts.library.is_some() {
        if opts.library_rg_ids.is_empty() {
            keep_ids.clear();
        } else {
            keep_ids.retain(|id| opts.library_rg_ids.contains(id));
        }
    }

    Some(keep_ids)
}

fn header_rg_id(line: &str) -> Option<&str> {
    line.split('\t')
        .skip(1)
        .find_map(|field| field.strip_prefix("ID:"))
}

/// Builds the CRAM encoder options from `-O cram,...` flags
/// (`embed_ref`, `seqs_per_slice`, `slices_per_slice`).
fn cram_write_options(opts: &Opts) -> htslib_rs::alignment_compat::CramWriteOptions {
    htslib_rs::alignment_compat::CramWriteOptions {
        embed_reference: opts.embed_reference,
        records_per_slice: opts.records_per_slice,
        slices_per_container: opts.slices_per_container,
    }
}

/// SAM-reader → CRAM writer, honoring `-O cram,embed_ref=1`
/// (`opts.embed_reference`) and the `seqs_per_slice` /
/// `slices_per_slice` slice/container sizing options.
fn cram_sam_reader_writer<R, Q, W>(
    opts: &Opts,
    reader: &mut htslib_rs::sam::io::Reader<R>,
    reference: Q,
    writer: W,
) -> io::Result<W>
where
    R: std::io::BufRead,
    Q: AsRef<Path>,
    W: io::Write,
{
    htslib_rs::alignment_compat::write_cram_from_sam_reader_with_reference_and_options(
        reader,
        reference,
        cram_write_options(opts),
        writer,
    )
}

fn write_cram_from_sam_text_via_bam<Q, W>(
    opts: &Opts,
    text: &[u8],
    reference: Q,
    writer: W,
) -> io::Result<W>
where
    Q: AsRef<Path>,
    W: io::Write,
{
    if opts.embed_reference {
        let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(
            sam_bytes_with_pg(text, opts)?,
        )));
        return cram_sam_reader_writer(opts, &mut reader, reference, writer);
    }

    let (mut tmp_bam, tmp_bam_path) =
        crate::tmp_file::create_temp_file("samtools-rs-view-rewrite", Some("bam"))?;
    htslib_rs::alignment_compat::write_bam_from_sam_reader(
        BufReader::new(io::Cursor::new(sam_bytes_with_pg(text, opts)?)),
        &mut tmp_bam,
    )?;
    drop(tmp_bam);

    htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference_and_options(
        tmp_bam_path.path(),
        reference,
        cram_write_options(opts),
        writer,
    )
}

/// Injects samtools' `@PG` chain entry into a SAM-text blob (header
/// lines split from the body), so binary BAM/CRAM output produced by
/// converting SAM text carries the `@PG` like upstream. A no-op under
/// `--no-PG` or when no argv was captured.
fn sam_bytes_with_pg(text: &[u8], opts: &Opts) -> io::Result<Vec<u8>> {
    let s = match std::str::from_utf8(text) {
        Ok(s) => s,
        Err(_) => return Ok(text.to_vec()),
    };
    let mut header_end = 0;
    for line in s.split_inclusive('\n') {
        if line.starts_with('@') {
            header_end += line.len();
        } else {
            break;
        }
    }
    let include_uri = resolved_output_fmt(opts)? == OutputFmt::Cram;
    let new_header = inject_reference_dictionary_if_missing(&s[..header_end], opts, include_uri)?;
    let new_header = apply_pg_to_header(&new_header, opts)?;
    let mut out = Vec::with_capacity(new_header.len() + (s.len() - header_end));
    out.extend_from_slice(new_header.as_bytes());
    out.extend_from_slice(&text[header_end..]);
    Ok(out)
}

fn open_text_output(opts: &Opts) -> io::Result<Box<dyn Write>> {
    if opts.write_index
        && resolved_output_fmt(opts)? == OutputFmt::Sam
        && let Some(path) = opts.output.as_deref()
    {
        let file = File::create(path)?;
        return Ok(Box::new(bgzf::io::Writer::new(file)));
    }
    sam_io::open_text_output(opts.output.as_deref())
}

fn open_binary_output(opts: &Opts) -> io::Result<Box<dyn Write>> {
    match opts.output.as_deref() {
        Some(path) => File::create(path).map(|file| Box::new(file) as Box<dyn Write>),
        None => Ok(Box::new(io::stdout())),
    }
}

fn open_unselected_text_output(opts: &Opts) -> io::Result<Option<Box<dyn Write>>> {
    opts.unselected_output
        .as_deref()
        .map(|path| sam_io::open_text_output(Some(path)))
        .transpose()
}

fn validate_unselected_sam_output(opts: &Opts) -> io::Result<()> {
    if opts.unselected_output.is_none() && !opts.unmap_unselected {
        return Ok(());
    }
    if opts.unselected_output.is_some() && opts.unmap_unselected {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "`-U/--output-unselected` and `-p/--unmap` are mutually exclusive",
        ));
    }
    Ok(())
}

fn reject_unselected_for_count(opts: &Opts) -> io::Result<()> {
    if opts.unselected_output.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "`-U/--output-unselected` is not supported with `-c` count output yet",
        ));
    }
    Ok(())
}

fn reject_unselected_binary_output(opts: &Opts) -> io::Result<()> {
    if opts.unselected_output.is_some() || opts.unmap_unselected {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "`-U/--output-unselected` and `-p/--unmap` for binary output are only supported for SAM input so far",
        ));
    }
    Ok(())
}

/// Builds two SAM text buffers — selected records and unselected records —
/// from a raw SAM byte buffer using the same filtering / `-p` / `-U`
/// logic as the SAM-output path. Returns `(selected, unselected)`.
/// `unselected` is empty when the caller did not pass `-U`.
///
/// Used by the SAM-input binary output paths to support `-p`/`-U` via a
/// text → BAM/CRAM roundtrip without needing to touch the binary record
/// representation directly.
fn build_split_sam_text(bytes: &[u8], opts: &Opts) -> io::Result<(Vec<u8>, Vec<u8>)> {
    let mut selected = Vec::with_capacity(bytes.len());
    let mut unselected = Vec::new();
    let header = sam_header_lines(bytes);
    selected.extend_from_slice(header);
    let want_unselected = opts.unselected_output.is_some();
    if want_unselected {
        unselected.extend_from_slice(header);
    }
    let tail = strip_header_lines(bytes);
    let expr_filter = opts
        .filter_expr
        .as_deref()
        .map(htslib_rs::expr::Filter::new);
    let regions = parse_simple_regions(&opts.regions)?;
    let bed_regions = parse_simple_regions(&opts.bed_regions)?;
    for line in tail.split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }
        if !sam_line_overlaps_requested_regions(line, &regions, &bed_regions) {
            continue;
        }
        if !line_selected(line, opts, expr_filter.as_ref())? {
            if want_unselected {
                write_sam_record_line(&mut unselected, line, opts)?;
            } else if opts.unmap_unselected {
                write_sam_unmapped_record_line(&mut selected, line, opts)?;
            }
            continue;
        }
        write_sam_record_line(&mut selected, line, opts)?;
    }
    Ok((selected, unselected))
}

/// Builds a SAM-text view suitable for binary `-U` / `-p` splitting.
///
/// The split logic is shared with SAM input, so BAM/CRAM inputs are first
/// rendered as SAM text, filtered/split, then encoded back to BAM or CRAM.
fn sam_text_for_binary_split(path: &Path, exact: Exact, opts: &Opts) -> io::Result<Vec<u8>> {
    let text = match exact {
        Exact::Sam => return read_sam_path_bytes(path),
        Exact::Bam => {
            let text =
                htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(path, None)?;
            sam_text_with_input_header(path, exact, text)?
        }
        Exact::Cram => {
            let reference_guard = cram_input_reference_for_path(opts, path)?;
            let text =
                htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                    path,
                    reference_guard.path(),
                    None,
                )?;
            sam_text_with_input_header(path, exact, text)?
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "binary unselected output is only wired up for SAM/BAM/CRAM input",
            ));
        }
    };
    Ok(crate::sam_render::fix_sam_text(&text).into_bytes())
}

fn filtered_binary_rewrite_text(
    path: &Path,
    exact: Exact,
    text: String,
    opts: &Opts,
    prefiltered: bool,
) -> io::Result<Vec<u8>> {
    let text = sam_text_with_input_header(path, exact, text)?;
    let text = crate::sam_render::fix_sam_text(&text);
    if prefiltered {
        let post_filter_opts = opts_after_prefiltered_expr(opts);
        filtered_sam_text(text.as_bytes(), &post_filter_opts)
    } else {
        filtered_sam_text(text.as_bytes(), opts)
    }
}

fn sam_text_with_input_header(path: &Path, exact: Exact, text: String) -> io::Result<String> {
    if text.as_bytes().starts_with(b"@") {
        return Ok(text);
    }
    let mut header = read_raw_header_text_with_format(path, exact)?;
    if !header.is_empty() && !header.ends_with('\n') {
        header.push('\n');
    }
    header.push_str(&text);
    Ok(header)
}

/// Returns the raw bytes of `path`, transparently decompressing BGZF input
/// the same way `filtered_sam_text_from_path` does. Used by binary output
/// paths that need a textual view of a SAM input for `-p`/`-U` handling.
fn read_sam_path_bytes(path: &Path) -> io::Result<Vec<u8>> {
    let file = File::open(path)?;
    let mut bytes = Vec::new();
    if is_bgzf_path(path)? {
        MultiGzDecoder::new(file).read_to_end(&mut bytes)?;
    } else {
        BufReader::new(file).read_to_end(&mut bytes)?;
    }
    Ok(bytes)
}

fn sam_text_from_path_for_region_filter(
    path: &Path,
    exact: Exact,
    opts: &Opts,
) -> io::Result<String> {
    match exact {
        Exact::Sam => String::from_utf8(read_sam_path_bytes(path)?)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        Exact::Bam => {
            htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(path, None)
        }
        Exact::Cram => {
            let reference_guard = cram_input_reference_for_path(opts, path)?;
            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                path,
                reference_guard.path(),
                None,
            )
        }
        _ => Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "region filtering is only wired up for SAM/BAM/CRAM input",
        )),
    }
}

fn load_bed_regions(path: &Path) -> io::Result<Vec<String>> {
    Ok(load_bed_index(path)?.to_region_strings())
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
    if !opts.bed_regions.is_empty() || has_filters(opts) {
        let text = sam_text_from_path_for_region_filter(path, exact, opts)?;
        return count_sam_text_region_records(text.as_bytes(), opts);
    }

    match exact {
        Exact::Sam => count_sam_region_records_from_path(path, opts),
        Exact::Bam => {
            let regions = parse_region_strings(path, &opts.regions)?;
            let mut n = 0usize;
            for region in &regions {
                n += htslib_rs::alignment_compat::count_bam_records_in_region_from_path(
                    path, region,
                )?;
            }
            Ok(n)
        }
        Exact::Cram => {
            let regions = parse_region_strings(path, &opts.regions)?;
            if cram_summary_count_supported(opts) {
                let mut n = 0usize;
                for region in &regions {
                    let records = if let Some(reference_guard) =
                        optional_cram_input_reference_for_path(opts, path)?
                    {
                        htslib_rs::alignment_compat::query_cram_records_from_path_with_reference(
                            path,
                            region,
                            reference_guard.path(),
                        )?
                    } else {
                        htslib_rs::alignment_compat::query_cram_records_from_path_synthesizing_reference(
                            path, region,
                        )?
                    };
                    n = n.saturating_add(records.len());
                }
                Ok(n)
            } else {
                let reference_guard = cram_input_reference_for_path(opts, path)?;
                let reference = reference_guard.path();
                let text =
                    htslib_rs::alignment_compat::view_cram_regions_as_sam_text_from_path_with_reference_and_limit(
                        path,
                        reference,
                        &regions,
                        None,
                    )?;
                count_sam_text_region_records(text.as_bytes(), opts)
            }
        }
        _ => Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "region count is only wired up for SAM/BAM/CRAM input",
        )),
    }
}

#[derive(Clone, Debug)]
struct SimpleRegion {
    reference: String,
    start: u64,
    end: u64,
}

fn parse_simple_regions(regions: &[String]) -> io::Result<Vec<SimpleRegion>> {
    regions
        .iter()
        .map(|region| parse_simple_region(region))
        .collect()
}

fn parse_simple_region(region: &str) -> io::Result<SimpleRegion> {
    let Some((reference, span)) = region.rsplit_once(':') else {
        return Ok(SimpleRegion {
            reference: region.to_string(),
            start: 1,
            end: u64::MAX,
        });
    };
    let (start, end) = match span.split_once('-') {
        Some((start, end)) => (
            parse_region_pos(start).unwrap_or(1),
            parse_region_pos(end).unwrap_or(u64::MAX),
        ),
        None => {
            let start = parse_region_pos(span).unwrap_or(1);
            (start, u64::MAX)
        }
    };
    Ok(SimpleRegion {
        reference: reference.to_string(),
        start,
        end,
    })
}

fn parse_region_pos(raw: &str) -> Option<u64> {
    raw.chars()
        .filter(|c| *c != ',')
        .collect::<String>()
        .parse()
        .ok()
}

fn count_sam_region_records_from_path(path: &Path, opts: &Opts) -> io::Result<usize> {
    let bytes = read_sam_path_bytes(path)?;
    count_sam_text_region_records(&bytes, opts)
}

fn count_sam_text_region_records(bytes: &[u8], opts: &Opts) -> io::Result<usize> {
    let regions = parse_simple_regions(&opts.regions)?;
    let bed_regions = parse_simple_regions(&opts.bed_regions)?;
    let expr_filter = opts
        .filter_expr
        .as_deref()
        .map(htslib_rs::expr::Filter::new);
    let mut count = 0usize;
    for line in strip_header_lines(bytes).split(|&b| b == b'\n') {
        if line.is_empty() || !sam_line_overlaps_requested_regions(line, &regions, &bed_regions) {
            continue;
        }
        if line_selected(line, opts, expr_filter.as_ref())? {
            count += 1;
        }
    }
    Ok(count)
}

fn count_records(path: &Path, exact: Exact, opts: &Opts) -> io::Result<usize> {
    match exact {
        Exact::Sam => htslib_rs::alignment_compat::count_sam_records_from_path(path),
        Exact::Bam => htslib_rs::alignment_compat::count_bam_records_from_path(path),
        Exact::Cram => Ok(summarize_cram_records_for_count(path, opts)?.len()),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported format",
        )),
    }
}

fn summarize_cram_records_for_count(
    path: &Path,
    opts: &Opts,
) -> io::Result<Vec<htslib_rs::alignment_compat::AlignmentRecordSummary>> {
    if let Some(reference_guard) = optional_cram_input_reference_for_path(opts, path)? {
        htslib_rs::alignment_compat::summarize_cram_records_from_path_with_reference(
            path,
            reference_guard.path(),
        )
    } else {
        htslib_rs::alignment_compat::summarize_cram_records_from_path_synthesizing_reference(path)
    }
}

fn write_records_as_sam<W: Write>(
    out: &mut W,
    unselected: &mut Option<Box<dyn Write>>,
    path: &Path,
    exact: Exact,
    opts: &Opts,
) -> io::Result<()> {
    let filter = prefilter_expr_for_sam_output(opts);
    match exact {
        Exact::Sam => {
            if has_sanitizer(opts) {
                let raw = read_sam_path_bytes(path)?;
                return write_sam_text_records_split(out, unselected, &raw, opts);
            }
            if let Some(expr) = filter.as_deref() {
                let text = htslib_rs::alignment_compat::view_sam_text_matching_filter_from_path(
                    path, expr,
                )?;
                let text = crate::sam_render::fix_sam_text(&text);
                let post_filter_opts = opts_after_prefiltered_expr(opts);
                return write_sam_text_records_split(
                    out,
                    unselected,
                    text.as_bytes(),
                    &post_filter_opts,
                );
            }
            stream_sam_records(out, unselected, path, opts)
        }
        Exact::Bam => {
            let text = if !opts.bed_regions.is_empty() || opts.regions.is_empty() {
                if let Some(expr) = filter.as_deref() {
                    htslib_rs::alignment_compat::view_bam_as_sam_text_matching_filter_from_path(
                        path, expr,
                    )?
                } else {
                    htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(
                        path, None,
                    )?
                }
            } else {
                let regions = parse_region_strings(path, &opts.regions)?;
                if let Some(expr) = filter.as_deref() {
                    htslib_rs::alignment_compat::view_bam_regions_as_sam_text_matching_filter_from_path(
                        path,
                        &regions,
                        expr,
                    )?
                } else {
                    htslib_rs::alignment_compat::view_bam_regions_as_sam_text_from_path(
                        path, &regions,
                    )?
                }
            };
            // For BAM input we already have SAM text. Records came from
            // binary, so fix noodles' plain-decimal float spelling to
            // htslib's `%g` form, then apply filters line-by-line.
            let text = crate::sam_render::fix_sam_text(&text);
            let post_filter_opts;
            let split_opts = if filter.is_some() {
                post_filter_opts = opts_after_prefiltered_expr(opts);
                &post_filter_opts
            } else {
                opts
            };
            write_sam_text_records_split(out, unselected, text.as_bytes(), split_opts)
        }
        Exact::Cram => {
            let text = if !opts.bed_regions.is_empty() || opts.regions.is_empty() {
                cram_sam_text_from_path_maybe_synthesizing_reference(
                    path,
                    opts,
                    filter.as_deref(),
                    None,
                )?
            } else {
                let reference_guard = cram_input_reference_for_path(opts, path)?;
                let reference = reference_guard.path();
                let regions = parse_region_strings(path, &opts.regions)?;
                if let Some(expr) = filter.as_deref() {
                    htslib_rs::alignment_compat::view_cram_regions_as_sam_text_matching_filter_from_path_with_reference(
                        path,
                        reference,
                        &regions,
                        expr,
                    )?
                } else {
                    htslib_rs::alignment_compat::view_cram_regions_as_sam_text_from_path_with_reference(
                        path, reference, &regions, false,
                    )?
                }
            };
            // CRAM records came from binary; fix float spelling.
            let text = crate::sam_render::fix_sam_text(&text);
            let post_filter_opts;
            let split_opts = if filter.is_some() {
                post_filter_opts = opts_after_prefiltered_expr(opts);
                &post_filter_opts
            } else {
                opts
            };
            write_sam_text_records_split(out, unselected, text.as_bytes(), split_opts)
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported format",
        )),
    }
}

fn cram_reference(opts: &Opts) -> io::Result<PathBuf> {
    explicit_reference_path(opts).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM input requires --reference / -T",
        )
    })
}

struct ReferenceGuard {
    path: PathBuf,
    cleanup: Vec<PathBuf>,
}

impl ReferenceGuard {
    fn new(path: PathBuf) -> Self {
        Self {
            path,
            cleanup: Vec::new(),
        }
    }

    fn temporary(path: PathBuf, cleanup: Vec<PathBuf>) -> Self {
        Self { path, cleanup }
    }

    fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for ReferenceGuard {
    fn drop(&mut self) {
        for path in self.cleanup.iter().rev() {
            let _ = fs::remove_file(path);
        }
    }
}

fn cram_output_reference_for_input(
    opts: &Opts,
    input: &Path,
    exact: Exact,
) -> io::Result<ReferenceGuard> {
    if let Some(reference) = explicit_reference_path(opts) {
        return Ok(ReferenceGuard::new(reference));
    }

    if let Some(reference) = reference_from_header_uri(input, exact)? {
        return Ok(reference);
    }

    reference_from_ref_path(input, exact)?.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM output requires --reference / -T, @SQ UR tags, or REF_PATH entries matching @SQ M5 tags",
        )
    })
}

fn cram_input_reference_for_path(opts: &Opts, input: &Path) -> io::Result<ReferenceGuard> {
    optional_cram_input_reference_for_path(opts, input)?.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM input requires --reference / -T, @SQ UR tags, or REF_PATH entries matching @SQ M5 tags",
        )
    })
}

fn optional_cram_input_reference_for_path(
    opts: &Opts,
    input: &Path,
) -> io::Result<Option<ReferenceGuard>> {
    if let Some(reference) = explicit_reference_path(opts) {
        return Ok(Some(ReferenceGuard::new(reference)));
    }

    if let Some(reference) = reference_from_header_uri(input, Exact::Cram)? {
        return Ok(Some(reference));
    }

    reference_from_ref_path(input, Exact::Cram)
}

fn reference_from_ref_path(input: &Path, exact: Exact) -> io::Result<Option<ReferenceGuard>> {
    let Some(ref_path) = std::env::var_os("REF_PATH") else {
        return Ok(None);
    };
    let ref_path = ref_path.to_string_lossy();
    let header_text = read_raw_header_text_with_format(input, exact)?;
    reference_from_ref_path_header(&header_text, &ref_path)
}

fn reference_from_header_uri(input: &Path, exact: Exact) -> io::Result<Option<ReferenceGuard>> {
    let header_text = read_raw_header_text_with_format(input, exact)?;
    reference_from_header_uri_text(&header_text)
}

fn reference_from_header_uri_text(header_text: &str) -> io::Result<Option<ReferenceGuard>> {
    for line in header_text.lines().filter(|line| line.starts_with("@SQ\t")) {
        for field in line.split('\t').skip(1) {
            let Some(uri) = field.strip_prefix("UR:") else {
                continue;
            };
            let path = uri.strip_prefix("file://").unwrap_or(uri);
            let path = PathBuf::from(path);
            if path.is_file() {
                return Ok(Some(ReferenceGuard::new(path)));
            }
        }
    }
    Ok(None)
}

fn reference_from_ref_path_header(
    header_text: &str,
    ref_path: &str,
) -> io::Result<Option<ReferenceGuard>> {
    let mut sequences = Vec::new();

    for line in header_text.lines().filter(|line| line.starts_with("@SQ\t")) {
        let mut name = None;
        let mut md5 = None;
        let mut len = None;
        for field in line.split('\t').skip(1) {
            if let Some(value) = field.strip_prefix("SN:") {
                name = Some(value);
            } else if let Some(value) = field.strip_prefix("LN:") {
                len = value.parse::<usize>().ok();
            } else if let Some(value) = field.strip_prefix("M5:") {
                md5 = Some(value);
            }
        }

        let (Some(name), Some(md5)) = (name, md5) else {
            continue;
        };
        if let Some(sequence) = read_ref_path_md5_sequence(ref_path, md5)? {
            sequences.push((name.to_string(), sequence));
        } else if let Some(len) = len {
            sequences.push((name.to_string(), "N".repeat(len)));
        }
    }

    if sequences.is_empty() {
        return Ok(None);
    }

    let fasta = temporary_reference_path("ref-path", "fa");
    {
        let mut out = File::create(&fasta)?;
        for (name, sequence) in &sequences {
            writeln!(out, ">{name}")?;
            out.write_all(sequence.as_bytes())?;
            out.write_all(b"\n")?;
        }
    }
    let fai = crate::reference::ensure_fai_index(&fasta, None)?;
    Ok(Some(ReferenceGuard::temporary(
        fasta.clone(),
        vec![fai, fasta],
    )))
}

fn cram_stdin_reference(opts: &Opts, input: &[u8]) -> io::Result<ReferenceGuard> {
    if let Some(reference) = explicit_reference_path(opts) {
        return Ok(ReferenceGuard::new(reference));
    }

    let header = htslib_rs::alignment_compat::read_cram_header(io::Cursor::new(input))?;
    let mut header_text = Vec::new();
    crate::sam_render::write_header(&mut header_text, &header)?;
    let header_text = String::from_utf8(header_text)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    if let Some(reference) = reference_from_header_uri_text(&header_text)? {
        return Ok(reference);
    }

    let Some(ref_path) = std::env::var_os("REF_PATH") else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM input requires --reference / -T, @SQ UR tags, or REF_PATH entries matching @SQ M5 tags",
        ));
    };

    reference_from_ref_path_header(&header_text, &ref_path.to_string_lossy())?.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM input requires --reference / -T, @SQ UR tags, or REF_PATH entries matching @SQ M5 tags",
        )
    })
}

fn read_ref_path_md5_sequence(ref_path: &str, md5: &str) -> io::Result<Option<String>> {
    for template in ref_path.split(':').filter(|part| !part.is_empty()) {
        let candidate = if template.contains("%s") {
            PathBuf::from(template.replace("%s", md5))
        } else {
            Path::new(template).join(md5)
        };
        match fs::read_to_string(&candidate) {
            Ok(text) => {
                let sequence: String = text.chars().filter(|c| !c.is_whitespace()).collect();
                if !sequence.is_empty() {
                    return Ok(Some(sequence));
                }
            }
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => return Err(e),
        }
    }
    Ok(None)
}

fn temporary_reference_path(stem: &str, ext: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        "samtools-rs-view-{stem}-{}-{nanos}.{ext}",
        std::process::id()
    ))
}

fn has_filters(opts: &Opts) -> bool {
    opts.require_flags != 0
        || opts.exclude_flags != 0
        || opts.exclude_all_flags != 0
        || opts.min_mapq != 0
        || opts.min_query_len != 0
        || opts.qname_filter.is_some()
        || !opts.read_groups.is_empty()
        || opts.exclude_no_rg
        || opts.library.is_some()
        || opts.aux_tag_filter.is_some()
        || opts.only_unplaced
        || opts.subsample.is_some()
}

fn has_sanitizer(opts: &Opts) -> bool {
    opts.sanitize_flags != SanitizeFlags::empty()
}

fn has_record_rewrite(opts: &Opts) -> bool {
    has_tag_filter(opts) || has_sanitizer(opts) || opts.remove_b || opts.remove_flags != 0
}

fn count_sam_text_records(bytes: &[u8], opts: &Opts) -> usize {
    strip_header_lines(bytes)
        .split(|&b| b == b'\n')
        .filter(|line| !line.is_empty())
        .filter(|line| !has_filters(opts) || line_passes(line, opts))
        .count()
}

fn filtered_sam_text_from_path(path: &Path, opts: &Opts) -> io::Result<Vec<u8>> {
    let file = File::open(path)?;
    let mut bytes = Vec::new();
    if is_bgzf_path(path)? {
        MultiGzDecoder::new(file).read_to_end(&mut bytes)?;
    } else {
        BufReader::new(file).read_to_end(&mut bytes)?;
    }
    filtered_sam_text(&bytes, opts)
}

fn filtered_sam_text(bytes: &[u8], opts: &Opts) -> io::Result<Vec<u8>> {
    let mut out = Vec::with_capacity(bytes.len());
    out.extend_from_slice(sam_header_lines(bytes));
    write_sam_text_records_split(&mut out, &mut None, bytes, opts)?;
    Ok(out)
}

fn write_sam_text_records_split<W: Write>(
    out: &mut W,
    unselected: &mut Option<Box<dyn Write>>,
    bytes: &[u8],
    opts: &Opts,
) -> io::Result<()> {
    if has_sanitizer(opts) {
        return write_sanitized_sam_text_records_split(out, unselected, bytes, opts);
    }

    let tail = strip_header_lines(bytes);
    let header_has_sq = sam_header_has_sq(bytes) || sam_reference_dictionary_present(opts);
    if has_filters(opts)
        || opts.filter_expr.is_some()
        || !opts.regions.is_empty()
        || !opts.bed_regions.is_empty()
        || has_tag_filter(opts)
        || opts.remove_b
        || opts.remove_flags != 0
        || unselected.is_some()
        || opts.unmap_unselected
    {
        let expr_filter = opts
            .filter_expr
            .as_deref()
            .map(htslib_rs::expr::Filter::new);
        let regions = parse_simple_regions(&opts.regions)?;
        let bed_regions = parse_simple_regions(&opts.bed_regions)?;
        for line in tail.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            validate_sam_record_line_with_header(line, header_has_sq)?;
            if !sam_line_overlaps_requested_regions(line, &regions, &bed_regions) {
                continue;
            }
            if !line_selected(line, opts, expr_filter.as_ref())? {
                if let Some(unselected) = unselected.as_mut() {
                    write_sam_record_line(unselected.as_mut(), line, opts)?;
                } else if opts.unmap_unselected {
                    write_sam_unmapped_record_line(out, line, opts)?;
                }
                continue;
            }
            write_sam_record_line(out, line, opts)?;
        }
        Ok(())
    } else {
        validate_sam_text_records(tail, header_has_sq)?;
        out.write_all(tail)
    }
}

fn write_sanitized_sam_text_records_split<W: Write>(
    out: &mut W,
    unselected: &mut Option<Box<dyn Write>>,
    bytes: &[u8],
    opts: &Opts,
) -> io::Result<()> {
    let mut reader = sam::io::Reader::new(BufReader::new(io::Cursor::new(bytes)));
    let header = reader.read_header()?;
    let raw_lines: Vec<&[u8]> = strip_header_lines(bytes)
        .split(|&b| b == b'\n')
        .filter(|line| !line.is_empty())
        .collect();
    let expr_filter = opts
        .filter_expr
        .as_deref()
        .map(htslib_rs::expr::Filter::new);
    let regions = parse_simple_regions(&opts.regions)?;
    let bed_regions = parse_simple_regions(&opts.bed_regions)?;

    for raw_line in raw_lines {
        if !sam_line_overlaps_requested_regions(raw_line, &regions, &bed_regions) {
            let mut skipped = RecordBuf::default();
            if reader.read_record_buf(&header, &mut skipped)? == 0 {
                break;
            }
            continue;
        }
        let selected = line_selected(raw_line, opts, expr_filter.as_ref())?;
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }

        if !selected && opts.unmap_unselected {
            set_record_unmapped_for_view(&mut record);
        }

        sanitize_record(&header, &mut record, opts.sanitize_flags);
        let line = record_to_sam_line(&header, &record)?;
        if !selected {
            if let Some(unselected) = unselected.as_mut() {
                write_prepared_sam_record_line(unselected.as_mut(), &line, opts)?;
            } else if opts.unmap_unselected {
                write_prepared_sam_record_line(out, &line, opts)?;
            }
            continue;
        }
        write_prepared_sam_record_line(out, &line, opts)?;
    }

    Ok(())
}

fn set_record_unmapped_for_view(record: &mut RecordBuf) {
    use sam::alignment::record::{Flags, MappingQuality};

    let mut flags = record.flags();
    flags.insert(Flags::UNMAPPED);
    *record.flags_mut() = flags;
    *record.mapping_quality_mut() = Some(MappingQuality::MIN);
    *record.cigar_mut() = Default::default();
    *record.template_length_mut() = 0;
}

fn record_to_sam_line(header: &sam::Header, record: &RecordBuf) -> io::Result<Vec<u8>> {
    use sam::alignment::io::Write as _;

    let mut buf = Vec::new();
    let mut writer = sam::io::Writer::new(&mut buf);
    writer.write_alignment_record(header, record)?;
    if buf.last() == Some(&b'\n') {
        buf.pop();
    }
    // noodles spells `f32` aux values as plain decimals; rewrite `:f:`
    // scalars and `B:f,` arrays to htslib's `%g` form for byte parity.
    let line =
        std::str::from_utf8(&buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    Ok(crate::sam_render::fix_sam_aux_floats(line).into_bytes())
}

fn write_prepared_sam_record_line<W: Write + ?Sized>(
    out: &mut W,
    line: &[u8],
    opts: &Opts,
) -> io::Result<()> {
    out.write_all(&render_sam_record_line(line, opts))?;
    out.write_all(b"\n")
}

fn write_sam_record_line<W: Write + ?Sized>(
    out: &mut W,
    line: &[u8],
    opts: &Opts,
) -> io::Result<()> {
    out.write_all(&render_sam_record_line(line, opts))?;
    out.write_all(b"\n")
}

fn render_sam_record_line(line: &[u8], opts: &Opts) -> Vec<u8> {
    let line = if opts.remove_b {
        remove_cigar_b_operator(line)
    } else {
        line.to_vec()
    };
    let line = if opts.remove_flags != 0 {
        remove_sam_record_flags(&line, opts.remove_flags)
    } else {
        line
    };

    if has_tag_filter(opts) {
        let filtered = apply_tag_filter(&line, opts);
        fix_sam_record_aux_floats(&filtered)
    } else {
        fix_sam_record_aux_floats(&line)
    }
}

fn remove_sam_record_flags(line: &[u8], remove_flags: u32) -> Vec<u8> {
    let mut fields: Vec<Vec<u8>> = line.split(|&b| b == b'\t').map(Vec::from).collect();
    if fields.len() < 2 {
        return line.to_vec();
    }
    let Some(flag) = std::str::from_utf8(&fields[1])
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
    else {
        return line.to_vec();
    };
    fields[1] = (flag & !remove_flags).to_string().into_bytes();
    fields.join(&b'\t')
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct SamCigarOp {
    len: usize,
    op: u8,
}

fn remove_cigar_b_operator(line: &[u8]) -> Vec<u8> {
    let mut fields: Vec<Vec<u8>> = line.split(|&b| b == b'\t').map(Vec::from).collect();
    if fields.len() < 11 || fields[5] == b"*" || !fields[5].contains(&b'B') {
        return line.to_vec();
    }

    let flag = std::str::from_utf8(&fields[1])
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(0);
    if flag & 0x4 != 0 {
        return line.to_vec();
    }

    let Some(cigar) = parse_sam_cigar_ops(&fields[5]) else {
        return line.to_vec();
    };
    if cigar.first().is_some_and(|op| op.op == b'B') {
        fields[1] = (flag | 0x4).to_string().into_bytes();
        return fields.join(&b'\t');
    }

    let mut seq = fields[9].clone();
    let mut qual = fields[10].clone();
    let has_qual = qual.len() == seq.len() && qual != b"*";
    let mut new_cigar: Vec<SamCigarOp> = Vec::with_capacity(cigar.len());
    let mut i = 0usize;
    let mut j = 0usize;
    let mut end_j: Option<usize> = None;
    let mut failed = false;

    for (idx, op) in cigar.iter().enumerate() {
        if op.op == b'B' {
            if idx == cigar.len() - 1 {
                break;
            }
            if op.len > j {
                failed = true;
                break;
            }

            let mut removed = 0usize;
            let mut cut = None;
            for t in (0..new_cigar.len()).rev() {
                if !sam_cigar_consumes_query(new_cigar[t].op) {
                    continue;
                }
                let len = new_cigar[t].len;
                if removed + len >= op.len {
                    let trim = op.len - removed;
                    new_cigar[t].len = new_cigar[t].len.saturating_sub(trim);
                    cut = Some(if new_cigar[t].len == 0 { t } else { t + 1 });
                    break;
                }
                removed += len;
            }

            let Some(cut) = cut else {
                failed = true;
                break;
            };
            new_cigar.truncate(cut);
            end_j = Some(j);
            j -= op.len;
            continue;
        }

        new_cigar.push(*op);
        if sam_cigar_consumes_query(op.op) {
            if i != j {
                for u in 0..op.len {
                    if i + u >= seq.len() || j + u >= seq.len() {
                        failed = true;
                        break;
                    }
                    let c = seq[i + u];
                    if end_j.is_some_and(|end| j + u < end) {
                        let c0 = seq[j + u];
                        if c != c0 {
                            if has_qual && phred(qual[j + u]) < phred(qual[i + u]) {
                                seq[j + u] = c;
                                qual[j + u] = phred_to_sam(phred(qual[i + u]) - phred(qual[j + u]));
                            } else if has_qual {
                                qual[j + u] = phred_to_sam(phred(qual[j + u]) - phred(qual[i + u]));
                            }
                        } else if has_qual {
                            qual[j + u] = qual[j + u].max(qual[i + u]);
                        }
                    } else {
                        seq[j + u] = c;
                        if has_qual {
                            qual[j + u] = qual[i + u];
                        }
                    }
                }
                if failed {
                    break;
                }
            }
            i += op.len;
            j += op.len;
        }
    }

    if failed {
        fields[1] = (flag | 0x4).to_string().into_bytes();
        return fields.join(&b'\t');
    }

    merge_adjacent_sam_cigar_ops(&mut new_cigar);
    let new_cigar = new_cigar
        .into_iter()
        .filter(|op| op.len != 0)
        .collect::<Vec<_>>();
    fields[5] = format_sam_cigar_ops(&new_cigar);
    seq.truncate(j);
    fields[9] = seq;
    if has_qual {
        qual.truncate(j);
        fields[10] = qual;
    }
    fields.join(&b'\t')
}

fn parse_sam_cigar_ops(cigar: &[u8]) -> Option<Vec<SamCigarOp>> {
    let mut ops = Vec::new();
    let mut len = 0usize;
    for &b in cigar {
        if b.is_ascii_digit() {
            len = len.checked_mul(10)?.checked_add(usize::from(b - b'0'))?;
            continue;
        }
        if len == 0
            || !matches!(
                b,
                b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'=' | b'X' | b'B'
            )
        {
            return None;
        }
        ops.push(SamCigarOp { len, op: b });
        len = 0;
    }
    (len == 0).then_some(ops)
}

fn validate_sam_text_records(bytes: &[u8], header_has_sq: bool) -> io::Result<()> {
    for line in bytes.split(|&b| b == b'\n') {
        if !line.is_empty() {
            validate_sam_record_line_with_header(line, header_has_sq)?;
        }
    }
    Ok(())
}

fn validate_sam_record_line_with_header(line: &[u8], header_has_sq: bool) -> io::Result<()> {
    let fields: Vec<_> = line.split(|&b| b == b'\t').collect();
    if fields.len() < 11 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "SAM record has fewer than 11 fields",
        ));
    }
    let cigar = fields[5];
    if cigar != b"*" && parse_sam_cigar_ops(cigar).is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid CIGAR {}", String::from_utf8_lossy(cigar)),
        ));
    }
    if !header_has_sq && fields[2] != b"*" && !record_flag_is_unmapped(fields[1]) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no SQ lines present in the header",
        ));
    }
    Ok(())
}

fn record_flag_is_unmapped(flag: &[u8]) -> bool {
    std::str::from_utf8(flag)
        .ok()
        .and_then(|s| s.parse::<u16>().ok())
        .is_some_and(|flag| flag & 0x4 != 0)
}

fn sam_record_line_is_mapped(line: &[u8]) -> bool {
    let mut fields = line.split(|&b| b == b'\t');
    let _qname = fields.next();
    let flag = fields.next().unwrap_or(b"");
    let reference = fields.next().unwrap_or(b"*");
    reference != b"*" && !record_flag_is_unmapped(flag)
}

fn sam_header_has_sq(bytes: &[u8]) -> bool {
    sam_header_lines(bytes)
        .split(|&b| b == b'\n')
        .any(|line| line.starts_with(b"@SQ\t") || line.starts_with(b"@SQ "))
}

/// `-t FILE` / `-T FILE` supply a reference dictionary that is injected as
/// `@SQ` lines into the output header (see
/// [`inject_reference_dictionary_if_missing`]). When either is given, mapped
/// records are valid even though the input SAM itself carries no `@SQ`, so
/// record validation must treat the effective header as having `@SQ`.
fn sam_reference_dictionary_present(opts: &Opts) -> bool {
    opts.reference_index.is_some() || explicit_reference_path(opts).is_some()
}

fn format_sam_cigar_ops(ops: &[SamCigarOp]) -> Vec<u8> {
    if ops.is_empty() {
        return b"*".to_vec();
    }
    let mut out = Vec::new();
    for op in ops {
        out.extend_from_slice(op.len.to_string().as_bytes());
        out.push(op.op);
    }
    out
}

fn merge_adjacent_sam_cigar_ops(ops: &mut [SamCigarOp]) {
    for i in 1..ops.len() {
        if ops[i].op == ops[i - 1].op {
            ops[i].len = ops[i].len.saturating_add(ops[i - 1].len);
            ops[i - 1].len = 0;
        }
    }
}

fn sam_cigar_consumes_query(op: u8) -> bool {
    matches!(op, b'M' | b'I' | b'S' | b'=' | b'X')
}

fn phred(q: u8) -> u8 {
    q.saturating_sub(33)
}

fn phred_to_sam(q: u8) -> u8 {
    q.saturating_add(33)
}

fn fix_sam_record_aux_floats(line: &[u8]) -> Vec<u8> {
    let Ok(text) = std::str::from_utf8(line) else {
        return line.to_vec();
    };
    crate::sam_render::fix_sam_aux_floats(text).into_bytes()
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
        let tag = [field[0], field[1]];
        let keep = if !opts.keep_tags.is_empty() {
            opts.keep_tags.contains(&tag)
        } else {
            !opts.remove_tags.contains(&tag)
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

fn extend_aux_tags(dst: &mut Vec<AuxTag>, raw: &str, option: &str) -> Result<(), ParseError> {
    let tags = parse_aux_list(raw)
        .map_err(|e| ParseError::Err(format!("invalid {option} value \"{raw}\": {e}")))?;
    dst.extend(tags);
    Ok(())
}

/// Apply view filters to a SAM record line, returning whether it should
/// be emitted. Parses the flag (column 2), MAPQ (column 5), and CIGAR
/// (column 6) when needed.
fn line_passes(line: &[u8], opts: &Opts) -> bool {
    let mut fields = line.split(|&b| b == b'\t');
    let qname = fields.next().unwrap_or(b"");
    let flag = fields
        .next()
        .and_then(|f| std::str::from_utf8(f).ok())
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(0);
    let rname = fields.next().unwrap_or(b"");
    if opts.only_unplaced && rname != b"*" {
        return false;
    }
    let _pos = fields.next();
    let mapq = fields
        .next()
        .and_then(|f| std::str::from_utf8(f).ok())
        .and_then(|s| s.parse::<u8>().ok())
        .unwrap_or(0);
    let cigar = fields.next().unwrap_or(b"");
    if opts.min_query_len != 0 && sam_cigar_query_len(cigar) < opts.min_query_len {
        return false;
    }
    if let Some(qfilter) = opts.qname_filter.as_ref()
        && !qfilter.matches(qname)
    {
        return false;
    }
    if let Some(subsample) = opts.subsample
        && !subsample_qname_passes(qname, subsample)
    {
        return false;
    }
    if !opts.read_groups.is_empty() || opts.exclude_no_rg {
        let rg = extract_rg_value(line);
        match rg {
            None if opts.exclude_no_rg => return false,
            None => {}
            Some(value) if !opts.read_groups.is_empty() && !opts.read_groups.contains(value) => {
                return false;
            }
            Some(_) => {}
        }
    }
    if opts.library.is_some() {
        // A record passes iff its read group's LB matches; no RG (or an
        // RG not in the resolved set) fails, like upstream's
        // bam_get_library check.
        match extract_rg_value(line) {
            Some(rg) if opts.library_rg_ids.contains(rg) => {}
            _ => return false,
        }
    }
    if let Some(filter) = opts.aux_tag_filter.as_ref() {
        let aux_value = extract_aux_value(line, filter.tag);
        match (&filter.values, aux_value) {
            (_, None) => return false,
            (Some(values), Some(value)) => {
                if !values.contains(value) {
                    return false;
                }
            }
            (None, Some(_)) => {}
        }
    }
    record_passes(flag, mapq, opts)
}

fn sam_line_overlaps_regions(line: &[u8], regions: &[SimpleRegion]) -> bool {
    if regions.is_empty() {
        return true;
    }
    let fields: Vec<_> = line.split(|&b| b == b'\t').collect();
    if fields.len() < 6 || fields[2] == b"*" {
        return false;
    }
    let Ok(reference) = std::str::from_utf8(fields[2]) else {
        return false;
    };
    let Some(start) = std::str::from_utf8(fields[3])
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
    else {
        return false;
    };
    let ref_len = sam_cigar_ref_len(fields[5]).max(1);
    let end = start.saturating_add(ref_len).saturating_sub(1);
    regions
        .iter()
        .any(|region| region.reference == reference && start <= region.end && region.start <= end)
}

fn sam_line_overlaps_requested_regions(
    line: &[u8],
    regions: &[SimpleRegion],
    bed_regions: &[SimpleRegion],
) -> bool {
    sam_line_overlaps_regions(line, regions) && sam_line_overlaps_regions(line, bed_regions)
}

fn sam_cigar_query_len(cigar: &[u8]) -> usize {
    let mut len = 0usize;
    let mut n = 0usize;
    for &b in cigar {
        if b.is_ascii_digit() {
            n = n.saturating_mul(10).saturating_add(usize::from(b - b'0'));
            continue;
        }
        if matches!(b, b'M' | b'I' | b'S' | b'=' | b'X') {
            len = len.saturating_add(n);
        }
        n = 0;
    }
    len
}

fn sam_cigar_ref_len(cigar: &[u8]) -> u64 {
    let mut len = 0u64;
    let mut n = 0u64;
    for &b in cigar {
        if b.is_ascii_digit() {
            n = n.saturating_mul(10).saturating_add(u64::from(b - b'0'));
            continue;
        }
        if matches!(b, b'M' | b'D' | b'N' | b'=' | b'X') {
            len = len.saturating_add(n);
        }
        n = 0;
    }
    len
}

fn subsample_qname_passes(qname: &[u8], subsample: Subsample) -> bool {
    if subsample.fraction >= 1.0 {
        return true;
    }
    if subsample.fraction <= 0.0 {
        return false;
    }
    let mut hash = 14695981039346656037u64 ^ u64::from(subsample.seed);
    for &b in qname {
        hash ^= u64::from(b);
        hash = hash.wrapping_mul(1099511628211);
    }
    hash ^= hash >> 30;
    hash = hash.wrapping_mul(0xbf58476d1ce4e5b9);
    hash ^= hash >> 27;
    hash = hash.wrapping_mul(0x94d049bb133111eb);
    hash ^= hash >> 31;
    (hash as f64 / (u64::MAX as f64 + 1.0)) < subsample.fraction
}

/// Walk a SAM record line and return the value of the aux tag `tag`, or
/// `None` if not present. Skips the first 11 mandatory fields. The raw
/// value bytes follow the `TAG:T:` prefix (5 bytes); for `B` array tags
/// the entire post-prefix portion (subtype prefix included) is returned.
fn extract_aux_value(line: &[u8], tag: [u8; 2]) -> Option<&[u8]> {
    for (i, field) in line.split(|&b| b == b'\t').enumerate() {
        if i < 11 {
            continue;
        }
        if field.len() < 5 || field[0] != tag[0] || field[1] != tag[1] || field[2] != b':' {
            continue;
        }
        return Some(&field[5..]);
    }
    None
}

fn line_selected(
    line: &[u8],
    opts: &Opts,
    expr_filter: Option<&htslib_rs::expr::Filter>,
) -> io::Result<bool> {
    if has_filters(opts) && !line_passes(line, opts) {
        return Ok(false);
    }
    if let Some(filter) = expr_filter {
        return line_expr_passes(line, filter);
    }
    Ok(true)
}

fn line_expr_passes(line: &[u8], filter: &htslib_rs::expr::Filter) -> io::Result<bool> {
    let context = SamLineFilterContext::new(line);
    let value = filter
        .eval_with(|symbol| context.lookup(symbol))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    Ok(value.truth())
}

struct SamLineFilterContext<'a> {
    fields: Vec<&'a [u8]>,
}

impl<'a> SamLineFilterContext<'a> {
    fn new(line: &'a [u8]) -> Self {
        Self {
            fields: line.split(|&b| b == b'\t').collect(),
        }
    }

    fn lookup(&self, symbol: &str) -> Option<(htslib_rs::expr::Value, usize)> {
        let (value, len) = if symbol.starts_with("qname") {
            (self.string_value(0), 5)
        } else if symbol.starts_with("rname") {
            (self.string_value(2), 5)
        } else if symbol.starts_with("rnext") {
            (self.rnext_value(), 5)
        } else if symbol.starts_with("mrname") {
            (self.rnext_value(), 6)
        } else if symbol.starts_with("cigar") {
            (self.optional_string_value(5, "*"), 5)
        } else if symbol.starts_with("seq") {
            (self.optional_string_value(9, "*"), 3)
        } else if symbol.starts_with("qual") {
            (self.optional_string_value(10, "*"), 4)
        } else if symbol.starts_with("pos") {
            (self.number_value(3), 3)
        } else if symbol.starts_with("pnext") {
            (self.number_value(7), 5)
        } else if symbol.starts_with("mpos") {
            (self.number_value(7), 4)
        } else if let Some((value, len)) = line_flag_expr_value(symbol, self.flag()) {
            (value, len)
        } else if symbol.starts_with("flag") {
            (self.number_value(1), 4)
        } else if symbol.starts_with("mapq") {
            (self.number_value(4), 4)
        } else if symbol.starts_with("qlen") {
            (htslib_rs::expr::Value::number(self.qlen() as f64), 4)
        } else if symbol.starts_with("rlen") {
            (htslib_rs::expr::Value::number(self.rlen() as f64), 4)
        } else if symbol.starts_with("sclen") {
            (htslib_rs::expr::Value::number(self.sclen() as f64), 5)
        } else if symbol.starts_with("hclen") {
            (htslib_rs::expr::Value::number(self.hclen() as f64), 5)
        } else if symbol.starts_with("endpos") {
            (self.endpos_value(), 6)
        } else if symbol.starts_with("ncigar") {
            (
                htslib_rs::expr::Value::number(self.cigar_metrics().op_count as f64),
                6,
            )
        } else if symbol.starts_with("tlen") {
            (self.number_value(8), 4)
        } else if let Some((tag, len)) = parse_expr_tag(symbol) {
            (self.aux_value(&tag), len)
        } else {
            return None;
        };
        Some((value, len))
    }

    fn string_value(&self, index: usize) -> htslib_rs::expr::Value {
        self.fields
            .get(index)
            .map(|field| htslib_rs::expr::Value::string(String::from_utf8_lossy(field)))
            .unwrap_or_else(htslib_rs::expr::Value::undefined)
    }

    fn rnext_value(&self) -> htslib_rs::expr::Value {
        self.fields
            .get(6)
            .map(|field| {
                let value = if *field == b"=" {
                    self.fields.get(2).copied().unwrap_or(b"*")
                } else {
                    *field
                };
                htslib_rs::expr::Value::string(String::from_utf8_lossy(value))
            })
            .unwrap_or_else(htslib_rs::expr::Value::undefined)
    }

    fn optional_string_value(&self, index: usize, missing: &str) -> htslib_rs::expr::Value {
        self.fields
            .get(index)
            .filter(|field| **field != missing.as_bytes())
            .map(|field| htslib_rs::expr::Value::string(String::from_utf8_lossy(field)))
            .unwrap_or_else(htslib_rs::expr::Value::undefined)
    }

    fn number_value(&self, index: usize) -> htslib_rs::expr::Value {
        self.fields
            .get(index)
            .and_then(|field| std::str::from_utf8(field).ok())
            .and_then(|s| s.parse::<f64>().ok())
            .map(htslib_rs::expr::Value::number)
            .unwrap_or_else(htslib_rs::expr::Value::undefined)
    }

    fn flag(&self) -> u16 {
        self.fields
            .get(1)
            .and_then(|field| std::str::from_utf8(field).ok())
            .and_then(|s| s.parse::<u16>().ok())
            .unwrap_or(0)
    }

    fn qlen(&self) -> usize {
        self.fields
            .get(9)
            .filter(|field| **field != b"*")
            .map_or(0, |field| field.len())
    }

    fn cigar_metrics(&self) -> CigarMetrics {
        self.fields
            .get(5)
            .map_or_else(CigarMetrics::default, |cigar| cigar_metrics(cigar))
    }

    fn rlen(&self) -> usize {
        self.cigar_metrics().reference_len
    }

    fn sclen(&self) -> usize {
        self.cigar_metrics().soft_clip_len
    }

    fn hclen(&self) -> usize {
        self.cigar_metrics().hard_clip_len
    }

    fn endpos_value(&self) -> htslib_rs::expr::Value {
        let pos = self
            .fields
            .get(3)
            .and_then(|field| std::str::from_utf8(field).ok())
            .and_then(|s| s.parse::<usize>().ok());
        match pos {
            Some(pos) if pos > 0 => {
                let rlen = self.rlen();
                htslib_rs::expr::Value::number((pos + rlen.saturating_sub(1)) as f64)
            }
            _ => htslib_rs::expr::Value::undefined(),
        }
    }

    fn aux_value(&self, tag: &[u8; 2]) -> htslib_rs::expr::Value {
        self.fields
            .iter()
            .skip(11)
            .find_map(|field| parse_aux_expr_value(field, tag))
            .unwrap_or_else(htslib_rs::expr::Value::undefined)
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct CigarMetrics {
    op_count: usize,
    reference_len: usize,
    soft_clip_len: usize,
    hard_clip_len: usize,
}

fn cigar_metrics(cigar: &[u8]) -> CigarMetrics {
    if cigar == b"*" {
        return CigarMetrics::default();
    }

    let mut ops = Vec::new();
    let mut len = 0usize;
    for &b in cigar {
        if b.is_ascii_digit() {
            len = len.saturating_mul(10).saturating_add(usize::from(b - b'0'));
            continue;
        }
        ops.push((b, len));
        len = 0;
    }

    let reference_len = ops
        .iter()
        .filter(|(op, _)| matches!(*op, b'M' | b'D' | b'N' | b'=' | b'X'))
        .map(|(_, len)| *len)
        .sum();

    let hard_clip_len = ops
        .iter()
        .enumerate()
        .filter(|(i, (op, _))| *op == b'H' && (*i == 0 || *i + 1 == ops.len()))
        .map(|(_, (_, len))| *len)
        .sum();

    let mut soft_clip_len = 0usize;
    let left = match ops.as_slice() {
        [(b'S', len), ..] => {
            soft_clip_len += *len;
            0
        }
        [(b'H', _), (b'S', len), ..] => {
            soft_clip_len += *len;
            1
        }
        _ => 0,
    };

    if ops.len().saturating_sub(1) > left && ops.last().is_some_and(|(op, _)| *op == b'S') {
        soft_clip_len += ops.last().map_or(0, |(_, len)| *len);
    } else if ops.len().saturating_sub(2) > left
        && ops.last().is_some_and(|(op, _)| *op == b'H')
        && ops
            .get(ops.len().saturating_sub(2))
            .is_some_and(|(op, _)| *op == b'S')
    {
        soft_clip_len += ops
            .get(ops.len().saturating_sub(2))
            .map_or(0, |(_, len)| *len);
    }

    CigarMetrics {
        op_count: ops.len(),
        reference_len,
        soft_clip_len,
        hard_clip_len,
    }
}

fn line_flag_expr_value(symbol: &str, flag: u16) -> Option<(htslib_rs::expr::Value, usize)> {
    let suffix = symbol.strip_prefix("flag.")?;
    let (mask, len): (u16, usize) = if suffix.starts_with("paired") {
        (0x1, "paired".len())
    } else if suffix.starts_with("proper_pair") {
        (0x2, "proper_pair".len())
    } else if suffix.starts_with("unmap") {
        (0x4, "unmap".len())
    } else if suffix.starts_with("munmap") {
        (0x8, "munmap".len())
    } else if suffix.starts_with("reverse") {
        (0x10, "reverse".len())
    } else if suffix.starts_with("mreverse") {
        (0x20, "mreverse".len())
    } else if suffix.starts_with("read1") {
        (0x40, "read1".len())
    } else if suffix.starts_with("read2") {
        (0x80, "read2".len())
    } else if suffix.starts_with("secondary") {
        (0x100, "secondary".len())
    } else if suffix.starts_with("qcfail") {
        (0x200, "qcfail".len())
    } else if suffix.starts_with("dup") {
        (0x400, "dup".len())
    } else if suffix.starts_with("supplementary") {
        (0x800, "supplementary".len())
    } else {
        return None;
    };

    Some((
        htslib_rs::expr::Value::number(f64::from(flag & mask)),
        "flag.".len() + len,
    ))
}

fn parse_expr_tag(symbol: &str) -> Option<([u8; 2], usize)> {
    let bytes = symbol.as_bytes();
    if bytes.len() < 4 || bytes[0] != b'[' || bytes[3] != b']' {
        return None;
    }
    Some(([bytes[1], bytes[2]], 4))
}

fn parse_aux_expr_value(field: &[u8], tag: &[u8; 2]) -> Option<htslib_rs::expr::Value> {
    if field.len() < 5 || field[0] != tag[0] || field[1] != tag[1] || field[2] != b':' {
        return None;
    }
    let value = std::str::from_utf8(&field[5..]).ok()?;
    match field[3] {
        b'A' | b'Z' | b'H' => Some(htslib_rs::expr::Value::string(value)),
        b'c' | b'C' | b's' | b'S' | b'i' | b'I' | b'f' => value
            .parse::<f64>()
            .ok()
            .map(htslib_rs::expr::Value::number),
        _ => None,
    }
}

fn count_expr_records(path: &Path, exact: Exact, opts: &Opts, expr: &str) -> io::Result<usize> {
    match exact {
        Exact::Sam if !opts.regions.is_empty() => Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "region queries on SAM input require an index (BAM/CRAM only)",
        )),
        Exact::Sam => {
            htslib_rs::alignment_compat::count_sam_records_matching_filter_from_path(path, expr)
        }
        Exact::Bam if !opts.regions.is_empty() => {
            let regions = parse_region_strings(path, &opts.regions)?;
            htslib_rs::alignment_compat::count_bam_records_in_regions_matching_filter_from_path(
                path, &regions, expr,
            )
        }
        Exact::Bam => {
            htslib_rs::alignment_compat::count_bam_records_matching_filter_from_path(path, expr)
        }
        Exact::Cram if !opts.regions.is_empty() => {
            let reference_guard = cram_input_reference_for_path(opts, path)?;
            let reference = reference_guard.path();
            let regions = parse_region_strings(path, &opts.regions)?;
            htslib_rs::alignment_compat::count_cram_records_in_regions_matching_filter_from_path_with_reference(
                path, reference, &regions, expr,
            )
        }
        Exact::Cram => {
            if let Some(reference_guard) = optional_cram_input_reference_for_path(opts, path)? {
                let reference = reference_guard.path();
                htslib_rs::alignment_compat::count_cram_records_matching_filter_from_path_with_reference(
                    path, reference, expr,
                )
            } else {
                count_cram_summary_expr_records(path, opts, expr)
            }
        }
        _ => Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "`-e EXPR` count is only wired up for SAM/BAM/CRAM input",
        )),
    }
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
        Exact::Cram => {
            if cram_summary_count_supported(opts) {
                count = summarize_cram_records_for_count(path, opts)?
                    .iter()
                    .filter(|record| cram_summary_passes(record, opts))
                    .count();
            } else {
                let reference_guard = cram_input_reference_for_path(opts, path)?;
                let reference = reference_guard.path();
                let text =
                    htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                        path, reference, None,
                    )?;
                count = count_sam_text_records(text.as_bytes(), opts);
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

fn count_cram_summary_expr_records(path: &Path, opts: &Opts, expr: &str) -> io::Result<usize> {
    let filter = htslib_rs::expr::Filter::new(expr);
    summarize_cram_records_for_count(path, opts)?
        .iter()
        .try_fold(0usize, |count, record| {
            summary_expr_passes(record, &filter).map(|passes| count + usize::from(passes))
        })
}

fn cram_sam_text_from_path_maybe_synthesizing_reference(
    path: &Path,
    opts: &Opts,
    filter: Option<&str>,
    limit: Option<usize>,
) -> io::Result<String> {
    if let Some(reference_guard) = optional_cram_input_reference_for_path(opts, path)? {
        let reference = reference_guard.path();
        if let Some(expr) = filter {
            htslib_rs::alignment_compat::view_cram_as_sam_text_matching_filter_from_path_with_reference(
                path, reference, expr,
            )
        } else {
            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                path, reference, limit,
            )
        }
    } else if opts.regions.is_empty() {
        if let Some(expr) = filter {
            htslib_rs::alignment_compat::view_cram_as_sam_text_matching_filter_from_path_synthesizing_reference(
                path, expr,
            )
        } else {
            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_synthesizing_reference_and_limit(
                path, limit,
            )
        }
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM region output requires --reference / -T, @SQ UR tags, or REF_PATH entries matching @SQ M5 tags",
        ))
    }
}

fn summary_expr_passes(
    record: &htslib_rs::alignment_compat::AlignmentRecordSummary,
    filter: &htslib_rs::expr::Filter,
) -> io::Result<bool> {
    let value = filter
        .eval_with(|symbol| summary_expr_value(record, symbol))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    Ok(value.truth())
}

fn summary_expr_value(
    record: &htslib_rs::alignment_compat::AlignmentRecordSummary,
    symbol: &str,
) -> Option<(htslib_rs::expr::Value, usize)> {
    let (value, len) = if symbol.starts_with("qname") {
        (
            record
                .name_bytes()
                .map(|name| htslib_rs::expr::Value::string(String::from_utf8_lossy(name)))
                .unwrap_or_else(htslib_rs::expr::Value::undefined),
            5,
        )
    } else if symbol.starts_with("pos") {
        (
            htslib_rs::expr::Value::number(record.alignment_start().unwrap_or(0) as f64),
            3,
        )
    } else if symbol.starts_with("pnext") {
        (
            htslib_rs::expr::Value::number(record.mate_alignment_start().unwrap_or(0) as f64),
            5,
        )
    } else if symbol.starts_with("mpos") {
        (
            htslib_rs::expr::Value::number(record.mate_alignment_start().unwrap_or(0) as f64),
            4,
        )
    } else if let Some((value, len)) = line_flag_expr_value(symbol, record.flags_u16()) {
        (value, len)
    } else if symbol.starts_with("flag") {
        (
            htslib_rs::expr::Value::number(f64::from(record.flags_u16())),
            4,
        )
    } else if symbol.starts_with("mapq") {
        (
            record
                .mapping_quality()
                .map(|mapq| htslib_rs::expr::Value::number(f64::from(mapq)))
                .unwrap_or_else(htslib_rs::expr::Value::undefined),
            4,
        )
    } else if symbol.starts_with("mrefid") {
        (
            htslib_rs::expr::Value::number(summary_ref_id_value(
                record.mate_reference_sequence_id(),
            )),
            6,
        )
    } else if symbol.starts_with("refid") {
        (
            htslib_rs::expr::Value::number(summary_ref_id_value(record.reference_sequence_id())),
            5,
        )
    } else if symbol.starts_with("qlen") {
        (
            htslib_rs::expr::Value::number(record.sequence_bytes().len() as f64),
            4,
        )
    } else if symbol.starts_with("tlen") {
        (
            htslib_rs::expr::Value::number(f64::from(record.template_length())),
            4,
        )
    } else {
        return None;
    };
    Some((value, len))
}

fn summary_ref_id_value(id: Option<usize>) -> f64 {
    id.map_or(-1.0, |id| id as f64)
}

fn cram_summary_count_supported(opts: &Opts) -> bool {
    opts.filter_expr.is_none() && opts.bed_regions.is_empty()
}

fn cram_summary_passes(
    record: &htslib_rs::alignment_compat::AlignmentRecordSummary,
    opts: &Opts,
) -> bool {
    if opts.only_unplaced && record.reference_sequence_id().is_some() {
        return false;
    }
    if opts.min_query_len != 0 && record.cigar_query_len() < opts.min_query_len {
        return false;
    }
    if let Some(qfilter) = opts.qname_filter.as_ref() {
        let Some(name) = record.name_bytes() else {
            return false;
        };
        if !qfilter.matches(name) {
            return false;
        }
    }
    if let Some(subsample) = opts.subsample {
        let Some(name) = record.name_bytes() else {
            return false;
        };
        if !subsample_qname_passes(name, subsample) {
            return false;
        }
    }
    if !opts.read_groups.is_empty() || opts.exclude_no_rg {
        match record.read_group_id() {
            None if opts.exclude_no_rg => return false,
            None => {}
            Some(value) if !opts.read_groups.is_empty() && !opts.read_groups.contains(value) => {
                return false;
            }
            Some(_) => {}
        }
    }
    if opts.library.is_some() {
        match record.read_group_id() {
            Some(rg) if opts.library_rg_ids.contains(rg) => {}
            _ => return false,
        }
    }
    if let Some(filter) = opts.aux_tag_filter.as_ref() {
        let aux_value = record.aux_value(filter.tag);
        match (&filter.values, aux_value) {
            (_, None) => return false,
            (Some(values), Some(value)) => {
                if !values.contains(value) {
                    return false;
                }
            }
            (None, Some(_)) => {}
        }
    }
    record_passes(
        u32::from(record.flags_u16()),
        record.mapping_quality().unwrap_or(0),
        opts,
    )
}

fn stream_sam_records<W: Write>(
    out: &mut W,
    unselected: &mut Option<Box<dyn Write>>,
    path: &Path,
    opts: &Opts,
) -> io::Result<()> {
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
    let mut header_has_sq = sam_reference_dictionary_present(opts);
    let mut line_no = 0usize;
    let expr_filter = opts
        .filter_expr
        .as_deref()
        .map(htslib_rs::expr::Filter::new);
    let regions = parse_simple_regions(&opts.regions)?;
    let bed_regions = parse_simple_regions(&opts.bed_regions)?;
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        line_no += 1;
        if !in_records {
            if line.starts_with(b"@") {
                if line.starts_with(b"@SQ\t") || line.starts_with(b"@SQ ") {
                    header_has_sq = true;
                }
                continue;
            }
            in_records = true;
        }
        let body = if line.last() == Some(&b'\n') {
            &line[..line.len() - 1]
        } else {
            &line[..]
        };
        if !header_has_sq && sam_record_line_is_mapped(body) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                SamParseNoSqError {
                    path: path.to_path_buf(),
                    line: line_no,
                },
            ));
        }
        validate_sam_record_line_with_header(body, header_has_sq)?;
        if !sam_line_overlaps_requested_regions(body, &regions, &bed_regions) {
            continue;
        }
        if !line_selected(body, opts, expr_filter.as_ref())? {
            if let Some(unselected) = unselected.as_mut() {
                write_sam_record_line(unselected.as_mut(), body, opts)?;
            } else if opts.unmap_unselected {
                write_sam_unmapped_record_line(out, body, opts)?;
            }
            continue;
        }
        if has_tag_filter(opts) || opts.remove_b || unselected.is_some() || opts.unmap_unselected {
            write_sam_record_line(out, body, opts)?;
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

fn write_sam_unmapped_record_line<W: Write + ?Sized>(
    out: &mut W,
    line: &[u8],
    opts: &Opts,
) -> io::Result<()> {
    let mut fields: Vec<Vec<u8>> = line.split(|&b| b == b'\t').map(Vec::from).collect();
    if fields.len() >= 11 {
        let flag = std::str::from_utf8(&fields[1])
            .ok()
            .and_then(|s| s.parse::<u32>().ok())
            .unwrap_or(0)
            | 0x4;
        fields[1] = flag.to_string().into_bytes();
        fields[4] = b"0".to_vec();
        fields[5] = b"*".to_vec();
        fields[8] = b"0".to_vec();
    }
    let unmapped = fields.join(&b'\t');
    write_sam_record_line(out, &unmapped, opts)
}

fn sam_header_lines(bytes: &[u8]) -> &[u8] {
    let tail = strip_header_lines(bytes);
    let len = bytes.len() - tail.len();
    &bytes[..len]
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
    writeln!(w, "  -m INT       minimum query length")?;
    writeln!(w, "  -o FILE      output FILE")?;
    writeln!(w, "  -U FILE      output records not selected by filters")?;
    writeln!(
        w,
        "  -p           set UNMAP on records not selected by filters"
    )?;
    writeln!(w, "  -z FLAGS     sanitize records before output")?;
    writeln!(w, "  -B           collapse backward CIGAR operations")?;
    writeln!(w, "  --remove-flags FLAGS")?;
    writeln!(w, "  -s FLOAT     subsample reads")?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn tmp_path(name: &str, ext: &str) -> PathBuf {
        static NEXT_TMP_ID: AtomicUsize = AtomicUsize::new(0);

        let id = NEXT_TMP_ID.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "samtools-rs-view-stdin-{}-{}-{}.{}",
            name,
            std::process::id(),
            id,
            ext
        ))
    }

    fn htslib_fixtures_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("repos")
            .join("htslib-rs")
            .join("repos")
            .join("htslib")
            .join("test")
    }

    fn stdin_sam() -> &'static [u8] {
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\nr2\t4\t*\t0\t0\t*\t*\t0\t0\tTG\t##\n"
    }

    fn stdin_bam() -> Vec<u8> {
        htslib_rs::alignment_compat::write_bam_from_sam_reader(
            BufReader::new(io::Cursor::new(stdin_sam())),
            Vec::new(),
        )
        .unwrap()
    }

    fn write_reference() -> (PathBuf, PathBuf) {
        let reference = tmp_path("ref", "fa");
        std::fs::write(
            &reference,
            b">ref\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        )
        .unwrap();
        let reference_index = crate::reference::ensure_fai_index(&reference, None).unwrap();
        (reference, reference_index)
    }

    fn assert_all_record_mapqs_at_least(text: &str, min_mapq: u8) {
        let mapqs: Vec<u8> = text
            .lines()
            .filter(|line| !line.starts_with('@'))
            .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
            .collect();
        assert!(!mapqs.is_empty());
        assert!(mapqs.iter().all(|mapq| *mapq >= min_mapq));
    }

    #[test]
    fn parses_output_fmt_option_version_and_embed_ref() {
        let args = vec![
            OsString::from("view"),
            OsString::from("-C"),
            OsString::from("--output-fmt-option"),
            OsString::from("version=3.0"),
            OsString::from("--output-fmt-option=embed_ref=1"),
            OsString::from("in.sam"),
        ];

        let opts = parse_args(&args).unwrap();

        assert!(matches!(opts.output_fmt, OutputFmt::Cram));
        assert!(opts.embed_reference);
        assert_eq!(opts.input.as_deref(), Some(Path::new("in.sam")));
    }

    #[test]
    fn stdin_count_expr_succeeds() {
        let out = tmp_path("count-expr", "txt");
        let opts = Opts {
            output: Some(out.clone()),
            count: true,
            filter_expr: Some("mapq >= 10".to_string()),
            ..Opts::default()
        };

        assert_eq!(
            run_sam_stdin(&opts, stdin_sam()).unwrap(),
            ExitCode::SUCCESS
        );
        assert_eq!(std::fs::read_to_string(&out).unwrap(), "1\n");
        let _ = std::fs::remove_file(out);
    }

    #[test]
    fn stdin_bam_output_succeeds() {
        let out = tmp_path("bam", "bam");
        let opts = Opts {
            output: Some(out.clone()),
            output_fmt: OutputFmt::Bam,
            ..Opts::default()
        };

        assert_eq!(
            run_sam_stdin(&opts, stdin_sam()).unwrap(),
            ExitCode::SUCCESS
        );
        assert_eq!(
            htslib_rs::alignment_compat::count_bam_records_from_path(&out).unwrap(),
            2
        );
        let _ = std::fs::remove_file(out);
    }

    #[test]
    fn stdin_cram_output_honors_mapq_filter() {
        let out = tmp_path("cram-filter", "cram");
        let (reference, reference_index) = write_reference();
        let opts = Opts {
            output: Some(out.clone()),
            output_fmt: OutputFmt::Cram,
            reference: Some(reference.clone()),
            min_mapq: 10,
            ..Opts::default()
        };

        assert_eq!(
            run_sam_stdin(&opts, stdin_sam()).unwrap(),
            ExitCode::SUCCESS
        );
        let text =
            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                &out, &reference, None,
            )
            .unwrap();
        let records: Vec<&str> = text.lines().filter(|line| !line.starts_with('@')).collect();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].split('\t').next(), Some("r1"));
        assert_eq!(records[0].split('\t').nth(4), Some("20"));
        let _ = std::fs::remove_file(out);
        let _ = std::fs::remove_file(reference);
        let _ = std::fs::remove_file(reference_index);
    }

    #[test]
    fn stdin_bam_bam_output_honors_mapq_filter() {
        let out = tmp_path("bam-bam-filter", "bam");
        let opts = Opts {
            output: Some(out.clone()),
            output_fmt: OutputFmt::Bam,
            min_mapq: 10,
            ..Opts::default()
        };

        assert_eq!(
            run_bam_stdin(&opts, &stdin_bam()).unwrap(),
            ExitCode::SUCCESS
        );

        let text =
            htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None)
                .unwrap();
        assert_all_record_mapqs_at_least(&text, 10);
        let _ = std::fs::remove_file(out);
    }

    #[test]
    fn stdin_bam_cram_output_honors_mapq_filter() {
        let out = tmp_path("bam-cram-filter", "cram");
        let (reference, reference_index) = write_reference();
        let opts = Opts {
            output: Some(out.clone()),
            output_fmt: OutputFmt::Cram,
            reference: Some(reference.clone()),
            min_mapq: 10,
            ..Opts::default()
        };

        assert_eq!(
            run_bam_stdin(&opts, &stdin_bam()).unwrap(),
            ExitCode::SUCCESS
        );

        let text =
            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                &out, &reference, None,
            )
            .unwrap();
        assert_all_record_mapqs_at_least(&text, 10);
        let _ = std::fs::remove_file(out);
        let _ = std::fs::remove_file(reference);
        let _ = std::fs::remove_file(reference_index);
    }

    #[test]
    fn stdin_bam_count_expr_succeeds() {
        let out = tmp_path("bam-count-expr", "txt");
        let opts = Opts {
            output: Some(out.clone()),
            count: true,
            filter_expr: Some("mapq >= 10".to_string()),
            ..Opts::default()
        };

        assert_eq!(
            run_bam_stdin(&opts, &stdin_bam()).unwrap(),
            ExitCode::SUCCESS
        );
        assert_eq!(std::fs::read_to_string(&out).unwrap(), "1\n");
        let _ = std::fs::remove_file(out);
    }

    #[test]
    fn stdin_bam_sam_output_succeeds() {
        let out = tmp_path("bam-sam", "sam");
        let opts = Opts {
            output: Some(out.clone()),
            header: HeaderMode::Include,
            min_mapq: 10,
            ..Opts::default()
        };

        assert_eq!(
            run_bam_stdin(&opts, &stdin_bam()).unwrap(),
            ExitCode::SUCCESS
        );

        let text = std::fs::read_to_string(&out).unwrap();
        assert!(text.starts_with("@HD\t"));
        assert!(text.contains("\tr1\t") || text.contains("\nr1\t"));
        assert!(!text.contains("\nr2\t"));
        let _ = std::fs::remove_file(out);
    }

    #[test]
    fn stdin_cram_count_expr_succeeds() {
        let fixtures = htslib_fixtures_dir();
        let reference = fixtures.join("ce.fa");
        let cram = std::fs::read(fixtures.join("range.cram")).unwrap();
        let out = tmp_path("cram-count-expr", "txt");
        let opts = Opts {
            output: Some(out.clone()),
            count: true,
            reference: Some(reference),
            filter_expr: Some("mapq >= 20".to_string()),
            ..Opts::default()
        };

        assert_eq!(run_cram_stdin(&opts, &cram).unwrap(), ExitCode::SUCCESS);
        let count = std::fs::read_to_string(&out)
            .unwrap()
            .trim()
            .parse::<usize>()
            .unwrap();
        assert!(count > 0);
        let _ = std::fs::remove_file(out);
    }

    #[test]
    fn stdin_cram_bam_output_succeeds() {
        let fixtures = htslib_fixtures_dir();
        let reference = fixtures.join("ce.fa");
        let cram = std::fs::read(fixtures.join("range.cram")).unwrap();
        let out = tmp_path("cram-bam", "bam");
        let opts = Opts {
            output: Some(out.clone()),
            output_fmt: OutputFmt::Bam,
            reference: Some(reference),
            filter_expr: Some("mapq >= 20".to_string()),
            ..Opts::default()
        };

        assert_eq!(run_cram_stdin(&opts, &cram).unwrap(), ExitCode::SUCCESS);
        assert!(htslib_rs::alignment_compat::count_bam_records_from_path(&out).unwrap() > 0);
        let _ = std::fs::remove_file(out);
    }

    #[test]
    fn stdin_cram_uses_header_ur_reference_without_dash_t() {
        let (reference, reference_index) = write_reference();
        let sam = tmp_path("cram-header-ur", "sam");
        let cram_path = tmp_path("cram-header-ur", "cram");
        let out = tmp_path("cram-header-ur", "sam.out");
        std::fs::write(
            &sam,
            format!(
                "@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\tUR:file://{}\n\
                 r1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\n",
                reference.display()
            ),
        )
        .unwrap();

        let mut cram_file = File::create(&cram_path).unwrap();
        htslib_rs::alignment_compat::write_cram_from_sam_path_with_reference(
            &sam,
            &reference,
            &mut cram_file,
        )
        .unwrap();
        drop(cram_file);

        let cram = std::fs::read(&cram_path).unwrap();
        let opts = Opts {
            output: Some(out.clone()),
            header: HeaderMode::Include,
            ..Opts::default()
        };

        assert_eq!(run_cram_stdin(&opts, &cram).unwrap(), ExitCode::SUCCESS);
        let text = std::fs::read_to_string(&out).unwrap();
        assert!(text.contains("\nr1\t"));

        let _ = std::fs::remove_file(out);
        let _ = std::fs::remove_file(cram_path);
        let _ = std::fs::remove_file(sam);
        let _ = std::fs::remove_file(reference);
        let _ = std::fs::remove_file(reference_index);
    }

    #[test]
    fn stdin_cram_bam_output_honors_mapq_filter() {
        let fixtures = htslib_fixtures_dir();
        let reference = fixtures.join("ce.fa");
        let cram = std::fs::read(fixtures.join("range.cram")).unwrap();
        let out = tmp_path("cram-bam-mapq-filter", "bam");
        let opts = Opts {
            output: Some(out.clone()),
            output_fmt: OutputFmt::Bam,
            reference: Some(reference),
            min_mapq: 20,
            ..Opts::default()
        };

        assert_eq!(run_cram_stdin(&opts, &cram).unwrap(), ExitCode::SUCCESS);
        let text =
            htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None)
                .unwrap();
        assert_all_record_mapqs_at_least(&text, 20);
        let _ = std::fs::remove_file(out);
    }

    #[test]
    fn stdin_cram_cram_output_honors_mapq_filter() {
        let fixtures = htslib_fixtures_dir();
        let reference = fixtures.join("ce.fa");
        let cram = std::fs::read(fixtures.join("range.cram")).unwrap();
        let out = tmp_path("cram-cram-mapq-filter", "cram");
        let opts = Opts {
            output: Some(out.clone()),
            output_fmt: OutputFmt::Cram,
            reference: Some(reference.clone()),
            min_mapq: 20,
            ..Opts::default()
        };

        assert_eq!(run_cram_stdin(&opts, &cram).unwrap(), ExitCode::SUCCESS);
        let text =
            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                &out, &reference, None,
            )
            .unwrap();
        assert_all_record_mapqs_at_least(&text, 20);
        let _ = std::fs::remove_file(out);
    }

    #[test]
    fn simple_region_single_position_spans_to_reference_end() {
        let regions = parse_simple_regions(&["chr1:15".to_string()]).unwrap();

        assert_eq!(regions[0].reference, "chr1");
        assert_eq!(regions[0].start, 15);
        assert_eq!(regions[0].end, u64::MAX);
    }
}
