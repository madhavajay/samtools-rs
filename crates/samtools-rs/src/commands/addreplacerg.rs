//! `samtools addreplacerg` — add or replace `@RG` lines and `RG:Z:` aux tags.
//!
//! Mirrors `main_addreplacerg` in `bam_addrprg.c`. Initial Rust port works on
//! **SAM/BAM input → SAM/BAM output**:
//!  - `-r '@RG\tID:foo\tSM:bar'` — full `@RG` line spec; merged into header.
//!  - `-r 'ID:foo'` — incremental tag form (one tag per `-r`); combined into
//!    a single `@RG` line.
//!  - `-R ID` — set every record's `RG:Z` to this existing ID. The ID
//!    must already be present in the input header (upstream rejects an
//!    unknown ID with "RG ID supplied does not exist in header").
//!  - no `-r` / `-R` — default to the first `@RG` ID in the input header
//!    (matching upstream); error only if the input has no `@RG` line.
//!  - `-m overwrite_all|orphan_only` — how to handle existing `RG:Z` tags.
//!    Matches upstream's two modes; `overwrite_all` is the default. In
//!    `overwrite_all` with `-r`, all other `@RG` header lines are removed
//!    so only the new one remains (mirrors `sam_hdr_remove_except`).
//!  - `-O sam|bam` — output format (default: sam).
//!  - `--no-PG` — silently accepted (no `@PG` is added by this port).
//!  - `-w` — overwrite an existing `@RG` header line with the same ID
//!    instead of erroring.
//!
//! **Pending:** CRAM input/output, paired-end mate update, full `orphan_first`
//! semantics.

use bstr::BString;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use flate2::read::MultiGzDecoder;
use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf};

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

#[derive(Clone, Copy)]
enum Mode {
    OrphanOnly,
    OverwriteAll,
}

/// Bundle of the resolved read-group rewrite parameters, threaded through
/// the SAM and binary rewrite paths.
struct RgRewrite<'a> {
    rg_line: Option<&'a str>,
    rg_id: &'a str,
    mode: Mode,
    overwrite_hdr_rg: bool,
    pg_argv: Option<&'a [OsString]>,
}

/// Entry point for `samtools addreplacerg`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut tag_pieces: Vec<String> = Vec::new();
    let mut replace_id: Option<String> = None;
    let mut mode = Mode::OverwriteAll;
    let mut output: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let mut output_fmt = OutFmt::Sam;
    let mut no_pg = false;
    let mut overwrite_hdr_rg = false;
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-r" => {
                let v = iter.next().and_then(|a| a.to_str()).unwrap_or("");
                tag_pieces.push(v.to_string());
            }
            "-R" => {
                replace_id = iter.next().and_then(|a| a.to_str()).map(|s| s.to_string());
            }
            "-m" => {
                let v = iter.next().and_then(|a| a.to_str()).unwrap_or("");
                mode = match v {
                    "overwrite_all" => Mode::OverwriteAll,
                    "orphan_only" => Mode::OrphanOnly,
                    _ => {
                        print_error(
                            "addreplacerg",
                            format!(
                                "unknown -m value \"{}\" (expected overwrite_all or orphan_only)",
                                v
                            ),
                        );
                        return ExitCode::from(1);
                    }
                };
            }
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "-O" => {
                let value = iter.next().and_then(|a| a.to_str()).unwrap_or("sam");
                output_fmt = match parse_output_format(value) {
                    Ok(fmt) => fmt,
                    Err(e) => {
                        print_error("addreplacerg", e);
                        return ExitCode::from(1);
                    }
                };
            }
            "--no-PG" => {
                no_pg = true;
            }
            "-w" => {
                overwrite_hdr_rg = true;
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("addreplacerg", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                }
            }
        }
    }

    let Some(input) = input else {
        let _ = print_usage();
        return ExitCode::from(1);
    };

    let format = match sam_io::sam_open_format(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error("addreplacerg", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam) {
        print_error(
            "addreplacerg",
            "only SAM and BAM input are currently supported (CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    let rg_line = match build_rg_line(&tag_pieces) {
        Ok(line) => line,
        Err(msg) => {
            print_error("addreplacerg", msg);
            return ExitCode::from(1);
        }
    };

    let rg_id = match (replace_id.as_ref(), rg_line.as_ref()) {
        (Some(id), _) => Some(id.to_string()),
        (None, Some(line)) => extract_id(line),
        _ => match read_first_rg_id_from_header(&input, format.exact) {
            Ok(id) => id,
            Err(e) => {
                print_error_errno(
                    "addreplacerg",
                    format!("failed to read header from \"{}\"", input.display()),
                    &e,
                );
                return ExitCode::from(1);
            }
        },
    };

    let Some(rg_id) = rg_id else {
        print_error(
            "addreplacerg",
            "an @RG ID is required (use `-r 'ID:...'` or `-R ID`) when the input has no @RG header",
        );
        return ExitCode::from(1);
    };

    // `-R ID` (no `-r`) requires the ID to already exist in the header,
    // matching upstream's "RG ID supplied does not exist in header" check.
    if replace_id.is_some() && rg_line.is_none() {
        match header_has_rg_id(&input, format.exact, &rg_id) {
            Ok(true) => {}
            Ok(false) => {
                print_error(
                    "addreplacerg",
                    "RG ID supplied does not exist in header. Supply full @RG line with -r instead?",
                );
                return ExitCode::from(1);
            }
            Err(e) => {
                print_error_errno(
                    "addreplacerg",
                    format!("failed to read header from \"{}\"", input.display()),
                    &e,
                );
                return ExitCode::from(1);
            }
        }
    }

    let pg_argv = if no_pg { None } else { Some(args) };
    let rewrite = RgRewrite {
        rg_line: rg_line.as_deref(),
        rg_id: &rg_id,
        mode,
        overwrite_hdr_rg,
        pg_argv,
    };
    let result = match (format.exact, output_fmt) {
        (Exact::Sam, OutFmt::Sam) => {
            let mut writer = match sam_io::open_text_output(output.as_deref()) {
                Ok(writer) => writer,
                Err(e) => {
                    print_error_errno("addreplacerg", "open -o output", &e);
                    return ExitCode::from(1);
                }
            };
            let result = rewrite_sam(&input, &mut writer, &rewrite);
            result.and_then(|()| sam_io::check_sam_close(&mut writer))
        }
        _ => rewrite_records(&input, output.as_deref(), output_fmt, &rewrite),
    };

    if let Err(e) = result {
        if e.kind() == io::ErrorKind::BrokenPipe {
            return ExitCode::SUCCESS;
        }
        print_error_errno("addreplacerg", "rewrite failed", &e);
        return ExitCode::from(1);
    }
    ExitCode::SUCCESS
}

#[derive(Clone, Copy)]
enum OutFmt {
    Sam,
    Bam,
}

fn parse_output_format(raw: &str) -> Result<OutFmt, String> {
    match raw.to_ascii_lowercase().as_str() {
        "sam" => Ok(OutFmt::Sam),
        "bam" => Ok(OutFmt::Bam),
        _ => Err(format!("unsupported output format \"{}\"", raw)),
    }
}

fn build_rg_line(pieces: &[String]) -> Result<Option<String>, String> {
    if pieces.is_empty() {
        return Ok(None);
    }
    // If a single piece already begins with @RG\t, use it as-is.
    if pieces.len() == 1 && pieces[0].starts_with("@RG\t") {
        return Ok(Some(pieces[0].clone()));
    }
    // Otherwise treat each piece as `KEY:VALUE` and assemble a tab-separated
    // @RG line. An ID must be present somewhere.
    let mut out = String::from("@RG");
    let mut have_id = false;
    for p in pieces {
        out.push('\t');
        out.push_str(p);
        if p.starts_with("ID:") {
            have_id = true;
        }
    }
    if !have_id {
        return Err("missing ID: field in -r tag pieces".into());
    }
    Ok(Some(out))
}

/// Returns whether the input header contains an `@RG` line with `rg_id`.
fn header_has_rg_id(input: &Path, exact: Exact, rg_id: &str) -> io::Result<bool> {
    let header_text = crate::header_text::read_raw_header_text_with_format(input, exact)?;
    Ok(header_text.lines().any(|line| {
        line.starts_with("@RG\t")
            && line
                .split('\t')
                .skip(1)
                .any(|field| field.strip_prefix("ID:") == Some(rg_id))
    }))
}

/// Read the input header and return the ID of the first `@RG` line, if any.
/// Matches upstream's default behavior when neither `-r` nor `-R` is supplied.
fn read_first_rg_id_from_header(input: &Path, exact: Exact) -> io::Result<Option<String>> {
    let header_text = crate::header_text::read_raw_header_text_with_format(input, exact)?;
    for line in header_text.lines() {
        if !line.starts_with("@RG\t") {
            continue;
        }
        for field in line.split('\t').skip(1) {
            if let Some(id) = field.strip_prefix("ID:") {
                return Ok(Some(id.to_string()));
            }
        }
    }
    Ok(None)
}

fn extract_id(rg_line: &str) -> Option<String> {
    for field in rg_line.split('\t') {
        if let Some(rest) = field.strip_prefix("ID:") {
            return Some(rest.to_string());
        }
    }
    None
}

fn rewrite_sam(path: &Path, out: &mut dyn Write, rw: &RgRewrite<'_>) -> io::Result<()> {
    let RgRewrite {
        rg_line,
        rg_id,
        mode,
        overwrite_hdr_rg,
        pg_argv,
    } = *rw;
    let file = File::open(path)?;
    let mut probe = File::open(path)?;
    let mut hdr = [0u8; 2];
    let n = io::Read::read(&mut probe, &mut hdr)?;
    let bgzf = n >= 2 && hdr[0] == 0x1f && hdr[1] == 0x8b;
    let reader: Box<dyn BufRead> = if bgzf {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut reader = reader;
    let rg_tag = format!("RG:Z:{}", rg_id);

    // Buffer header lines so we can apply upstream's @RG header semantics
    // (mirrors `bam_addrprg.c` init_state) before emitting any record bytes.
    let mut header_lines: Vec<String> = Vec::new();
    let mut first_record: Option<String> = None;
    let mut line = String::new();
    loop {
        line.clear();
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            break;
        }
        if line.starts_with('@') {
            header_lines.push(line.trim_end_matches(['\r', '\n']).to_string());
        } else {
            first_record = Some(line.clone());
            break;
        }
    }

    if let Some(rg) = rg_line {
        let existing_rg_id = header_lines.iter().any(|l| {
            l.starts_with("@RG") && extract_id(l).as_deref().is_some_and(|id| id == rg_id)
        });
        if existing_rg_id {
            if !overwrite_hdr_rg {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "RG line with ID:{rg_id} already present in the header. Use -w to overwrite."
                    ),
                ));
            }
            header_lines.retain(|l| {
                !(l.starts_with("@RG") && extract_id(l).as_deref().is_some_and(|id| id == rg_id))
            });
        }
        header_lines.push(rg.trim_end_matches(['\r', '\n']).to_string());
        if matches!(mode, Mode::OverwriteAll) {
            header_lines
                .retain(|l| !l.starts_with("@RG") || extract_id(l).as_deref() == Some(rg_id));
        }
    }

    let mut header = String::new();
    if header_lines.is_empty() && rg_line.is_some() {
        header.push_str("@HD\tVN:1.6\n");
    }
    for l in &header_lines {
        header.push_str(l);
        header.push('\n');
    }

    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg(&header, argv)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    out.write_all(header.as_bytes())?;

    if let Some(first) = first_record {
        let stripped = first.trim_end_matches(&['\r', '\n'][..]);
        let new = rewrite_record(stripped, &rg_tag, mode);
        out.write_all(new.as_bytes())?;
        out.write_all(b"\n")?;
    }

    loop {
        line.clear();
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            break;
        }
        let stripped = line.trim_end_matches(&['\r', '\n'][..]);
        let new = rewrite_record(stripped, &rg_tag, mode);
        out.write_all(new.as_bytes())?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

fn rewrite_record(line: &str, new_rg: &str, mode: Mode) -> String {
    let mut fields: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
    let mut existing_rg_idx: Option<usize> = None;
    for (i, f) in fields.iter().enumerate().skip(11) {
        if f.starts_with("RG:Z:") {
            existing_rg_idx = Some(i);
            break;
        }
    }
    match (existing_rg_idx, mode) {
        (Some(idx), Mode::OverwriteAll) => {
            fields[idx] = new_rg.to_string();
        }
        (Some(_), Mode::OrphanOnly) => {
            // Leave existing RG untouched.
        }
        (None, _) => {
            fields.push(new_rg.to_string());
        }
    }
    fields.join("\t")
}

fn rewrite_records(
    input: &Path,
    output: Option<&Path>,
    output_fmt: OutFmt,
    rw: &RgRewrite<'_>,
) -> io::Result<()> {
    let format = sam_io::sam_open_format(input)?;
    let (mut header, mut records) = match format.exact {
        Exact::Sam => read_sam_records(input)?,
        Exact::Bam => read_bam_records(input)?,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "only SAM and BAM input are currently supported (CRAM TODO)",
            ));
        }
    };

    header = update_header(
        header,
        rw.rg_line,
        rw.rg_id,
        rw.mode,
        rw.overwrite_hdr_rg,
        rw.pg_argv,
    )?;
    let rg_id = rw.rg_id;
    let mode = rw.mode;
    for record in &mut records {
        rewrite_record_buf(record, rg_id, mode);
    }

    let mut writer = open_record_output(output, output_fmt, &header)?;
    for record in &records {
        writer.write_record(&header, record)?;
    }
    Ok(())
}

fn read_sam_records(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(input)?));
    let header = reader.read_header()?;
    let mut records = Vec::new();
    loop {
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        records.push(record);
    }
    Ok((header, records))
}

fn read_bam_records(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;
    let mut records = Vec::new();
    loop {
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        records.push(record);
    }
    Ok((header, records))
}

fn update_header(
    mut header: sam::Header,
    rg_line: Option<&str>,
    rg_id: &str,
    mode: Mode,
    overwrite_hdr_rg: bool,
    pg_argv: Option<&[OsString]>,
) -> io::Result<sam::Header> {
    if let Some(rg) = rg_line {
        let rg_id_key = BString::from(rg_id);
        let already_present = header.read_groups().contains_key(rg_id.as_bytes());
        if already_present {
            if !overwrite_hdr_rg {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "RG line with ID:{rg_id} already present in the header. Use -w to overwrite."
                    ),
                ));
            }
            header.read_groups_mut().shift_remove(&rg_id_key);
        }

        let mut parser = sam::header::Parser::default();
        parser
            .parse_partial(rg.as_bytes())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let parsed = parser.finish();
        for (id, read_group) in parsed.read_groups() {
            header
                .read_groups_mut()
                .insert(id.clone(), read_group.clone());
        }

        if matches!(mode, Mode::OverwriteAll) {
            header.read_groups_mut().retain(|k, _| *k == rg_id_key);
        }
    }

    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

    Ok(header)
}

fn rewrite_record_buf(record: &mut RecordBuf, rg_id: &str, mode: Mode) {
    use sam::alignment::record::data::field::Tag;
    use sam::alignment::record_buf::data::field::Value;

    let tag = Tag::from([b'R', b'G']);
    let has_rg = record.data().get(&tag).is_some();
    match (has_rg, mode) {
        (true, Mode::OrphanOnly) => {}
        _ => {
            record
                .data_mut()
                .insert(tag, Value::String(BString::from(rg_id)));
        }
    }
}

trait RecordSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(File);
struct SamStdout(io::Stdout);

impl RecordSink for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}

impl RecordSink for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}

impl RecordSink for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        // Shared renderer: htslib `%g` float aux spelling.
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}

impl RecordSink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}

fn open_record_output(
    output: Option<&Path>,
    fmt: OutFmt,
    header: &sam::Header,
) -> io::Result<Box<dyn RecordSink>> {
    match (output, fmt) {
        (Some(path), OutFmt::Sam) => {
            let mut file = File::create(path)?;
            crate::sam_render::write_header(&mut file, header)?;
            Ok(Box::new(SamFile(file)))
        }
        (Some(path), OutFmt::Bam) => {
            let file = File::create(path)?;
            let mut writer = bam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(BamFile(writer)))
        }
        (None, OutFmt::Sam) => {
            let mut stdout = io::stdout();
            crate::sam_render::write_header(&mut stdout, header)?;
            Ok(Box::new(SamStdout(stdout)))
        }
        (None, OutFmt::Bam) => {
            let mut writer = bam::io::Writer::new(io::stdout());
            writer.write_header(header)?;
            Ok(Box::new(BamStdout(writer)))
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(
        w,
        "Usage: samtools addreplacerg [options] -r 'tag-spec' <in.sam|in.bam>"
    )?;
    writeln!(
        w,
        "  -r SPEC       @RG line or 'KEY:VALUE' tag (repeatable)"
    )?;
    writeln!(w, "  -R ID         existing @RG ID to apply")?;
    writeln!(
        w,
        "  -m MODE       overwrite_all | orphan_only [overwrite_all]"
    )?;
    writeln!(w, "  -o FILE       output FILE (default stdout)")?;
    writeln!(w, "  -O FMT        output format: sam|bam")?;
    Ok(())
}
