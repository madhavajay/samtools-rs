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
//!  - `-U FILE` — for SAM output, write records not selected by flag/MAPQ filters to FILE.
//!  - `-p` — for SAM output, set UNMAP on records not selected by flag/MAPQ filters.
//!  - `-T FILE` / `--reference FILE` — reference for CRAM I/O.
//!  - `-u` — write uncompressed BAM (accepted; treated as `-b -1` for now).
//!  - `-1` — fast compression level (accepted; treated as `-b` default for now).
//!
//! Anything else returns a "not yet supported" error so that test failures
//! are loud rather than silent.

use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{Category, Exact};
use htslib_rs::sam::{self, alignment::RecordBuf};

use crate::aux_list::{AuxTag, parse_aux_list};
use crate::bedidx::load_bed_index;
use crate::diagnostics::{print_error, print_error_errno};
use crate::header_text::read_raw_header_text_with_format;
use crate::io as sam_io;
use crate::sam_global::current_global_args;
use crate::sanitize::{SanitizeFlags, parse_sanitize_options, sanitize_record};

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
            Ok(more) => opts.regions.extend(more),
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
                        &reference,
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
            code
        }
        // Broken pipe from a downstream consumer (e.g. `samtools view | head`)
        // is a clean exit, not an error — matches upstream behavior.
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("view", "I/O error during view", &e);
            ExitCode::from(1)
        }
    }
}

/// Post-write index pass for `--write-index` (BAM file output only).
fn write_output_index(opts: &Opts) -> io::Result<()> {
    let Some(out) = opts.output.as_deref() else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--write-index requires -o FILE",
        ));
    };
    match resolved_output_fmt(opts)? {
        OutputFmt::Bam => {
            let index = htslib_rs::index_compat::build_bai(out)?;
            let mut idx = out.as_os_str().to_os_string();
            idx.push(".bai");
            htslib_rs::index_compat::write_bai(PathBuf::from(idx), &index)
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--write-index is only supported for BAM file output in samtools-rs view",
        )),
    }
}

#[derive(Default)]
struct Opts {
    input: Option<PathBuf>,
    output: Option<PathBuf>,
    unselected_output: Option<PathBuf>,
    output_fmt: OutputFmt,
    header: HeaderMode,
    count: bool,
    no_pg: bool,
    /// Argv captured for `@PG` insertion when `--no-PG` is not set.
    /// `None` means the caller didn't supply an argv (e.g. internal tests).
    argv: Option<Vec<OsString>>,
    unmap_unselected: bool,
    /// `--write-index` — build an index next to a BAM file output.
    write_index: bool,
    /// Region `*` — emit only unplaced (RNAME `*`) records.
    only_unplaced: bool,
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
    remove_tags: Vec<AuxTag>,
    /// `--keep-tag TAG` (repeatable) — only these aux tags are kept.
    keep_tags: Vec<AuxTag>,
    /// `-z FLAGS` / `--sanitize FLAGS` — upstream-style record sanitizer.
    sanitize_flags: SanitizeFlags,
    /// `-N FILE` / `--qname-file FILE` — read names listed in FILE (or
    /// `^FILE` to negate). Records whose qname appears in the set pass;
    /// `^FILE` flips to exclude. `None` means the filter is disabled.
    qname_filter: Option<QnameFilter>,
    /// `-r STR` / `-R FILE` — accumulated read-group IDs. Records whose
    /// `RG:Z:` aux value is in the set pass. Empty means the filter is off.
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
            "-X" | "--customized-index" => {
                opts.customized_index = true;
                i += 1;
            }
            "--write-index" => {
                opts.write_index = true;
                i += 1;
            }
            "--help" => return Err(ParseError::Usage),
            // Thread count: accepted and recorded. Output is byte-identical
            // regardless of the value (worker-pool wiring is a perf-only
            // follow-up — TODO-NEXT #8); `-@ N`, `-@N`, `--threads N`.
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

fn run(opts: &Opts, input: &Path, input_exact: Exact) -> io::Result<ExitCode> {
    let effective_out_fmt = resolved_output_fmt(opts)?;

    // Count-only mode.
    if opts.count {
        reject_unselected_for_count(opts)?;
        let filter = combined_filter_expr(opts);
        let n = if let Some(expr) = filter.as_deref() {
            count_expr_records(input, input_exact, opts, expr)?
        } else if !opts.regions.is_empty() {
            count_region_records(input, input_exact, opts)?
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
            apply_pg_to_header(&read_raw_header_text_with_format(input, input_exact)?, opts)?;
        let mut out = open_text_output(opts)?;
        out.write_all(header_text.as_bytes())?;
        sam_io::check_sam_close(&mut out)?;
        return Ok(ExitCode::SUCCESS);
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
                apply_pg_to_header(&read_raw_header_text_with_format(input, input_exact)?, opts)?;
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
        if input_exact != Exact::Sam {
            reject_unselected_binary_output(opts)?;
        }
        let filter = combined_filter_expr(opts);
        let dst_file = open_binary_output(opts)?;
        match input_exact {
            Exact::Sam => {
                let needs_split = opts.unselected_output.is_some() || opts.unmap_unselected;
                if needs_split {
                    // Build SAM text with `-p`/`-U` semantics applied, then
                    // pipe each side into a BAM writer. This avoids
                    // touching binary records directly while still
                    // honoring the splitting modes for SAM input.
                    let raw = read_sam_path_bytes(input)?;
                    let (selected, unselected) = build_split_sam_text(&raw, opts)?;
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
                } else if let Some(expr) = filter.as_deref() {
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
                if has_sanitizer(opts) {
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
                    let filtered = filtered_sam_text(text.as_bytes(), opts)?;
                    htslib_rs::alignment_compat::write_bam_from_sam_reader(
                        BufReader::new(io::Cursor::new(filtered)),
                        dst_file,
                    )?;
                } else if opts.regions.is_empty() {
                    if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::write_bam_matching_filter_from_path(
                            input, expr, dst_file,
                        )?;
                    } else {
                        htslib_rs::alignment_compat::write_bam_from_path(input, dst_file)?;
                    }
                } else {
                    let regions = parse_region_strings(input, &opts.regions)?;
                    if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::write_bam_regions_matching_filter_from_path(
                            input, &regions, expr, dst_file,
                        )?;
                    } else {
                        htslib_rs::alignment_compat::write_bam_regions_from_path(
                            input, &regions, dst_file,
                        )?;
                    }
                }
            }
            Exact::Cram => {
                let reference = cram_reference(opts)?;
                if has_sanitizer(opts) {
                    let text = if opts.regions.is_empty() {
                        if let Some(expr) = filter.as_deref() {
                            htslib_rs::alignment_compat::view_cram_as_sam_text_matching_filter_from_path_with_reference(
                                input, &reference, expr,
                            )?
                        } else {
                            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                                input, &reference, None,
                            )?
                        }
                    } else {
                        let regions = parse_region_strings(input, &opts.regions)?;
                        if let Some(expr) = filter.as_deref() {
                            htslib_rs::alignment_compat::view_cram_regions_as_sam_text_matching_filter_from_path_with_reference(
                                input, &reference, &regions, expr,
                            )?
                        } else {
                            htslib_rs::alignment_compat::view_cram_regions_as_sam_text_from_path_with_reference(
                                input, &reference, &regions, false,
                            )?
                        }
                    };
                    let filtered = filtered_sam_text(text.as_bytes(), opts)?;
                    htslib_rs::alignment_compat::write_bam_from_sam_reader(
                        BufReader::new(io::Cursor::new(filtered)),
                        dst_file,
                    )?;
                } else if opts.regions.is_empty() {
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
        if input_exact != Exact::Sam {
            reject_unselected_binary_output(opts)?;
        }
        let filter = combined_filter_expr(opts);
        let dst = opts
            .output
            .clone()
            .ok_or_else(|| io::Error::other("CRAM output to stdout requires -o file (TODO)"))?;
        let reference = cram_reference(opts)?;
        let dst_file = File::create(&dst)?;
        match input_exact {
            Exact::Sam => {
                let needs_split = opts.unselected_output.is_some() || opts.unmap_unselected;
                if needs_split {
                    let raw = read_sam_path_bytes(input)?;
                    let (selected, unselected) = build_split_sam_text(&raw, opts)?;
                    let mut selected_reader =
                        htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(selected)));
                    htslib_rs::alignment_compat::write_cram_from_sam_reader_with_reference(
                        &mut selected_reader,
                        &reference,
                        dst_file,
                    )?;
                    if let Some(unselected_path) = opts.unselected_output.as_deref() {
                        let mut unselected_reader = htslib_rs::sam::io::Reader::new(
                            BufReader::new(io::Cursor::new(unselected)),
                        );
                        let unselected_dst = File::create(unselected_path)?;
                        htslib_rs::alignment_compat::write_cram_from_sam_reader_with_reference(
                            &mut unselected_reader,
                            &reference,
                            unselected_dst,
                        )?;
                    }
                } else if let Some(expr) = filter.as_deref() {
                    if has_record_rewrite(opts) {
                        let text =
                            htslib_rs::alignment_compat::view_sam_text_matching_filter_from_path(
                                input, expr,
                            )?;
                        let filtered = filtered_sam_text(text.as_bytes(), opts)?;
                        let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(
                            io::Cursor::new(filtered),
                        ));
                        htslib_rs::alignment_compat::write_cram_from_sam_reader_with_reference(
                            &mut reader,
                            reference,
                            dst_file,
                        )?;
                    } else {
                        htslib_rs::alignment_compat::write_cram_matching_filter_from_sam_path_with_reference(
                            input, reference, expr, dst_file,
                        )?;
                    }
                } else if has_filters(opts) || has_record_rewrite(opts) {
                    let filtered = filtered_sam_text_from_path(input, opts)?;
                    let mut reader =
                        htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(filtered)));
                    htslib_rs::alignment_compat::write_cram_from_sam_reader_with_reference(
                        &mut reader,
                        reference,
                        dst_file,
                    )?;
                } else {
                    htslib_rs::alignment_compat::write_cram_from_sam_path_with_reference(
                        input, reference, dst_file,
                    )?;
                }
            }
            Exact::Bam if opts.regions.is_empty() => {
                if has_sanitizer(opts) {
                    let text = if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::view_bam_as_sam_text_matching_filter_from_path(
                            input, expr,
                        )?
                    } else {
                        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(
                            input, None,
                        )?
                    };
                    let filtered = filtered_sam_text(text.as_bytes(), opts)?;
                    let mut reader =
                        htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(filtered)));
                    htslib_rs::alignment_compat::write_cram_from_sam_reader_with_reference(
                        &mut reader,
                        reference,
                        dst_file,
                    )?;
                } else if let Some(expr) = filter.as_deref() {
                    htslib_rs::alignment_compat::write_cram_matching_filter_from_bam_path_with_reference(
                        input, reference, expr, dst_file,
                    )?;
                } else {
                    htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference(
                        input, reference, dst_file,
                    )?;
                }
            }
            Exact::Bam => {
                let regions = parse_region_strings(input, &opts.regions)?;
                if has_sanitizer(opts) {
                    let text = if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::view_bam_regions_as_sam_text_matching_filter_from_path(
                            input, &regions, expr,
                        )?
                    } else {
                        htslib_rs::alignment_compat::view_bam_regions_as_sam_text_from_path(
                            input, &regions,
                        )?
                    };
                    let filtered = filtered_sam_text(text.as_bytes(), opts)?;
                    let mut reader =
                        htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(filtered)));
                    htslib_rs::alignment_compat::write_cram_from_sam_reader_with_reference(
                        &mut reader,
                        reference,
                        dst_file,
                    )?;
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
                if has_sanitizer(opts) {
                    let text = if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::view_cram_as_sam_text_matching_filter_from_path_with_reference(
                            input, &reference, expr,
                        )?
                    } else {
                        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                            input, &reference, None,
                        )?
                    };
                    let filtered = filtered_sam_text(text.as_bytes(), opts)?;
                    let mut reader =
                        htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(filtered)));
                    htslib_rs::alignment_compat::write_cram_from_sam_reader_with_reference(
                        &mut reader,
                        reference,
                        dst_file,
                    )?;
                } else if let Some(expr) = filter.as_deref() {
                    htslib_rs::alignment_compat::write_cram_matching_filter_from_path_with_reference(
                        input, reference, expr, dst_file,
                    )?;
                } else {
                    htslib_rs::alignment_compat::write_cram_from_path_with_reference(
                        input, reference, dst_file,
                    )?;
                }
            }
            Exact::Cram => {
                let regions = parse_region_strings(input, &opts.regions)?;
                if has_sanitizer(opts) {
                    let text = if let Some(expr) = filter.as_deref() {
                        htslib_rs::alignment_compat::view_cram_regions_as_sam_text_matching_filter_from_path_with_reference(
                            input, &reference, &regions, expr,
                        )?
                    } else {
                        htslib_rs::alignment_compat::view_cram_regions_as_sam_text_from_path_with_reference(
                            input, &reference, &regions, false,
                        )?
                    };
                    let filtered = filtered_sam_text(text.as_bytes(), opts)?;
                    let mut reader =
                        htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(filtered)));
                    htslib_rs::alignment_compat::write_cram_from_sam_reader_with_reference(
                        &mut reader,
                        reference,
                        dst_file,
                    )?;
                } else if let Some(expr) = filter.as_deref() {
                    htslib_rs::alignment_compat::write_cram_regions_matching_filter_from_path_with_reference(
                        input, reference, &regions, expr, dst_file,
                    )?;
                } else {
                    htslib_rs::alignment_compat::write_cram_regions_from_path_with_reference(
                        input, reference, &regions, dst_file,
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
        let header_text = apply_pg_to_header(
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
            let header_text = apply_pg_to_header(
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
            write_sam_text_records_split(&mut out, &mut unselected, text.as_bytes(), opts)?;
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
                    BufReader::new(io::Cursor::new(filtered)),
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
                BufReader::new(io::Cursor::new(reader_input)),
                dst_file,
            )?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Cram {
        reject_unselected_binary_output(opts)?;
        let dst = opts
            .output
            .clone()
            .ok_or_else(|| io::Error::other("CRAM output to stdout requires -o file (TODO)"))?;
        let reference = cram_reference(opts)?;
        let dst_file = File::create(&dst)?;
        let filter = combined_filter_expr(opts);
        if let Some(expr) = filter.as_deref() {
            if has_record_rewrite(opts) {
                let text = htslib_rs::alignment_compat::view_sam_text_matching_filter(
                    BufReader::new(io::Cursor::new(input)),
                    expr,
                )?;
                let filtered = filtered_sam_text(text.as_bytes(), opts)?;
                let mut reader =
                    htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(filtered)));
                htslib_rs::alignment_compat::write_cram_from_sam_reader_with_reference(
                    &mut reader,
                    reference,
                    dst_file,
                )?;
            } else {
                let mut reader =
                    htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(input)));
                htslib_rs::alignment_compat::write_cram_matching_filter_from_sam_reader_with_reference(
                    &mut reader,
                    reference,
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
            let mut reader =
                htslib_rs::sam::io::Reader::new(BufReader::new(io::Cursor::new(reader_input)));
            htslib_rs::alignment_compat::write_cram_from_sam_reader_with_reference(
                &mut reader,
                reference,
                dst_file,
            )?;
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
            apply_pg_to_header(std::str::from_utf8(text.as_bytes()).unwrap_or(""), opts)?;
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
            let header_text = apply_pg_to_header(
                std::str::from_utf8(sam_header_lines(text.as_bytes())).unwrap_or(""),
                opts,
            )?;
            out.write_all(header_text.as_bytes())?;
            if let Some(unselected) = unselected.as_mut() {
                unselected.write_all(header_text.as_bytes())?;
            }
        }
        write_sam_text_records_split(&mut out, &mut unselected, text.as_bytes(), opts)?;
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
        let dst = opts
            .output
            .clone()
            .ok_or_else(|| io::Error::other("CRAM output to stdout requires -o file (TODO)"))?;
        let reference = cram_reference(opts)?;
        let dst_file = File::create(&dst)?;
        if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::write_cram_matching_filter_from_bam_reader_with_reference(
                io::Cursor::new(input),
                &reference,
                expr,
                dst_file,
            )?;
        } else {
            htslib_rs::alignment_compat::write_cram_from_bam_reader_with_reference(
                io::Cursor::new(input),
                &reference,
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

    let reference = cram_reference(opts)?;
    let effective_out_fmt = resolved_output_fmt(opts)?;

    if opts.count {
        reject_unselected_for_count(opts)?;
        let filter = combined_filter_expr(opts);
        let n = if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::count_cram_records_matching_filter_with_reference(
                io::Cursor::new(input),
                &reference,
                expr,
            )?
        } else {
            let text = htslib_rs::alignment_compat::view_cram_as_sam_text_with_reference(
                io::Cursor::new(input),
                &reference,
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
            &reference,
            Some(0),
        )?;
        let header_text =
            apply_pg_to_header(std::str::from_utf8(text.as_bytes()).unwrap_or(""), opts)?;
        let mut out = open_text_output(opts)?;
        out.write_all(header_text.as_bytes())?;
        sam_io::check_sam_close(&mut out)?;
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Sam {
        validate_unselected_sam_output(opts)?;
        let filter = prefilter_expr_for_sam_output(opts);
        let text = if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::view_cram_as_sam_text_matching_filter_with_reference(
                io::Cursor::new(input),
                &reference,
                expr,
            )?
        } else {
            htslib_rs::alignment_compat::view_cram_as_sam_text_with_reference(
                io::Cursor::new(input),
                &reference,
                None,
            )?
        };
        let mut out = open_text_output(opts)?;
        let mut unselected = open_unselected_text_output(opts)?;
        if opts.header == HeaderMode::Include {
            let header_text = apply_pg_to_header(
                std::str::from_utf8(sam_header_lines(text.as_bytes())).unwrap_or(""),
                opts,
            )?;
            out.write_all(header_text.as_bytes())?;
            if let Some(unselected) = unselected.as_mut() {
                unselected.write_all(header_text.as_bytes())?;
            }
        }
        write_sam_text_records_split(&mut out, &mut unselected, text.as_bytes(), opts)?;
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
                &reference,
                expr,
                dst_file,
            )?;
        } else {
            htslib_rs::alignment_compat::write_cram_records_with_required_flags_as_bam_with_reference(
                io::Cursor::new(input),
                &reference,
                0,
                dst_file,
            )?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    if effective_out_fmt == OutputFmt::Cram {
        reject_unselected_binary_output(opts)?;
        let filter = combined_filter_expr(opts);
        let dst = opts
            .output
            .clone()
            .ok_or_else(|| io::Error::other("CRAM output to stdout requires -o file (TODO)"))?;
        let dst_file = File::create(&dst)?;
        if let Some(expr) = filter.as_deref() {
            htslib_rs::alignment_compat::write_cram_matching_filter_from_reader_with_reference(
                io::Cursor::new(input),
                &reference,
                expr,
                dst_file,
            )?;
        } else {
            htslib_rs::alignment_compat::write_cram_from_reader_with_reference(
                io::Cursor::new(input),
                &reference,
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

/// Injects samtools' `@PG` chain entry into a SAM-text blob (header
/// lines split from the body), so binary BAM/CRAM output produced by
/// converting SAM text carries the `@PG` like upstream. A no-op under
/// `--no-PG` or when no argv was captured.
fn sam_bytes_with_pg(text: &[u8], opts: &Opts) -> io::Result<Vec<u8>> {
    if opts.no_pg || opts.argv.is_none() {
        return Ok(text.to_vec());
    }
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
    let new_header = apply_pg_to_header(&s[..header_end], opts)?;
    let mut out = Vec::with_capacity(new_header.len() + (s.len() - header_end));
    out.extend_from_slice(new_header.as_bytes());
    out.extend_from_slice(&text[header_end..]);
    Ok(out)
}

fn open_text_output(opts: &Opts) -> io::Result<Box<dyn Write>> {
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
    if opts.unselected_output.is_some() || opts.unmap_unselected {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "`-U/--output-unselected` and `-p/--unmap` are not supported with `-c` count output yet",
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
    for line in tail.split(|&b| b == b'\n') {
        if line.is_empty() {
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
    let regions = parse_region_strings(path, &opts.regions)?;
    match exact {
        Exact::Bam => {
            let mut n = 0usize;
            for region in &regions {
                n += htslib_rs::alignment_compat::count_bam_records_in_region_from_path(
                    path, region,
                )?;
            }
            Ok(n)
        }
        Exact::Cram => {
            let reference = cram_reference(opts)?;
            let text =
                htslib_rs::alignment_compat::view_cram_regions_as_sam_text_from_path_with_reference_and_limit(
                    path,
                    reference,
                    &regions,
                    None,
                )?;
            Ok(count_sam_text_records(text.as_bytes(), opts))
        }
        _ => Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "region count is only wired up for BAM/CRAM input",
        )),
    }
}

fn count_records(path: &Path, exact: Exact, opts: &Opts) -> io::Result<usize> {
    match exact {
        Exact::Sam => htslib_rs::alignment_compat::count_sam_records_from_path(path),
        Exact::Bam => htslib_rs::alignment_compat::count_bam_records_from_path(path),
        Exact::Cram => {
            let reference = cram_reference(opts)?;
            let text =
                htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                    path, reference, None,
                )?;
            Ok(count_sam_text_records(text.as_bytes(), opts))
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported format",
        )),
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
            if !opts.regions.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "region queries on SAM input require an index (BAM/CRAM only)",
                ));
            }
            if has_sanitizer(opts) {
                let raw = read_sam_path_bytes(path)?;
                return write_sam_text_records_split(out, unselected, &raw, opts);
            }
            if let Some(expr) = filter.as_deref() {
                let text = htslib_rs::alignment_compat::view_sam_text_matching_filter_from_path(
                    path, expr,
                )?;
                return write_sam_text_records_split(out, unselected, text.as_bytes(), opts);
            }
            stream_sam_records(out, unselected, path, opts)
        }
        Exact::Bam => {
            let text = if opts.regions.is_empty() {
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
            write_sam_text_records_split(out, unselected, text.as_bytes(), opts)
        }
        Exact::Cram => {
            let reference = cram_reference(opts)?;
            let text = if opts.regions.is_empty() {
                if let Some(expr) = filter.as_deref() {
                    htslib_rs::alignment_compat::view_cram_as_sam_text_matching_filter_from_path_with_reference(
                        path, reference, expr,
                    )?
                } else {
                    htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                        path,
                        reference,
                        None,
                    )?
                }
            } else {
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
            write_sam_text_records_split(out, unselected, text.as_bytes(), opts)
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported format",
        )),
    }
}

fn cram_reference(opts: &Opts) -> io::Result<PathBuf> {
    opts.reference
        .clone()
        .or_else(|| current_global_args().reference)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "CRAM input requires --reference / -T",
            )
        })
}

fn has_filters(opts: &Opts) -> bool {
    opts.require_flags != 0
        || opts.exclude_flags != 0
        || opts.exclude_all_flags != 0
        || opts.min_mapq != 0
        || opts.qname_filter.is_some()
        || !opts.read_groups.is_empty()
        || opts.exclude_no_rg
        || opts.library.is_some()
        || opts.aux_tag_filter.is_some()
        || opts.only_unplaced
}

fn has_sanitizer(opts: &Opts) -> bool {
    opts.sanitize_flags != SanitizeFlags::empty()
}

fn has_record_rewrite(opts: &Opts) -> bool {
    has_tag_filter(opts) || has_sanitizer(opts)
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
    if has_filters(opts)
        || opts.filter_expr.is_some()
        || has_tag_filter(opts)
        || unselected.is_some()
        || opts.unmap_unselected
    {
        let expr_filter = opts
            .filter_expr
            .as_deref()
            .map(htslib_rs::expr::Filter::new);
        for line in tail.split(|&b| b == b'\n') {
            if line.is_empty() {
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

    for raw_line in raw_lines {
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
    if has_tag_filter(opts) {
        out.write_all(&apply_tag_filter(line, opts))?;
    } else {
        out.write_all(line)?;
    }
    out.write_all(b"\n")
}

fn write_sam_record_line<W: Write + ?Sized>(
    out: &mut W,
    line: &[u8],
    opts: &Opts,
) -> io::Result<()> {
    if has_tag_filter(opts) {
        let filtered = apply_tag_filter(line, opts);
        out.write_all(&filtered)?;
    } else {
        out.write_all(line)?;
    }
    out.write_all(b"\n")
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
/// be emitted. Parses the flag (column 2) and MAPQ (column 5).
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
    if let Some(qfilter) = opts.qname_filter.as_ref()
        && !qfilter.matches(qname)
    {
        return false;
    }
    if !opts.read_groups.is_empty() || opts.exclude_no_rg {
        let rg = extract_rg_value(line);
        match rg {
            None if opts.exclude_no_rg => return false,
            None => {
                if !opts.read_groups.is_empty() {
                    return false;
                }
            }
            Some(value) => {
                if !opts.read_groups.is_empty() && !opts.read_groups.contains(value) {
                    return false;
                }
            }
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
        } else if symbol.starts_with("cigar") {
            (self.optional_string_value(5, "*"), 5)
        } else if symbol.starts_with("seq") {
            (self.optional_string_value(9, "*"), 3)
        } else if symbol.starts_with("qual") {
            (self.optional_string_value(10, "*"), 4)
        } else if symbol.starts_with("pos") {
            (self.number_value(3), 3)
        } else if symbol.starts_with("flag") {
            (self.number_value(1), 4)
        } else if symbol.starts_with("mapq") {
            (self.number_value(4), 4)
        } else if symbol.starts_with("qlen") {
            (htslib_rs::expr::Value::number(self.qlen() as f64), 4)
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

    fn qlen(&self) -> usize {
        self.fields
            .get(9)
            .filter(|field| **field != b"*")
            .map_or(0, |field| field.len())
    }

    fn aux_value(&self, tag: &[u8; 2]) -> htslib_rs::expr::Value {
        self.fields
            .iter()
            .skip(11)
            .find_map(|field| parse_aux_expr_value(field, tag))
            .unwrap_or_else(htslib_rs::expr::Value::undefined)
    }
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
            let reference = cram_reference(opts)?;
            let regions = parse_region_strings(path, &opts.regions)?;
            htslib_rs::alignment_compat::count_cram_records_in_regions_matching_filter_from_path_with_reference(
                path, reference, &regions, expr,
            )
        }
        Exact::Cram => {
            let reference = cram_reference(opts)?;
            htslib_rs::alignment_compat::count_cram_records_matching_filter_from_path_with_reference(
                path, reference, expr,
            )
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
            let reference = cram_reference(opts)?;
            let text =
                htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                    path, reference, None,
                )?;
            count = count_sam_text_records(text.as_bytes(), opts);
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
    let expr_filter = opts
        .filter_expr
        .as_deref()
        .map(htslib_rs::expr::Filter::new);
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
        if !line_selected(body, opts, expr_filter.as_ref())? {
            if let Some(unselected) = unselected.as_mut() {
                write_sam_record_line(unselected.as_mut(), body, opts)?;
            } else if opts.unmap_unselected {
                write_sam_unmapped_record_line(out, body, opts)?;
            }
            continue;
        }
        if has_tag_filter(opts) || unselected.is_some() || opts.unmap_unselected {
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
    writeln!(w, "  -o FILE      output FILE")?;
    writeln!(w, "  -U FILE      output records not selected by filters")?;
    writeln!(
        w,
        "  -p           set UNMAP on records not selected by filters"
    )?;
    writeln!(w, "  -z FLAGS     sanitize records before output")?;
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
            .join("htslib-rs")
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
}
