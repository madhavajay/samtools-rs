//! `samtools reheader` — replace the header of an alignment file.
//!
//! Mirrors `main_reheader` in `bam_reheader.c`. The basic mode is
//! `samtools reheader <new.hdr.sam> <in.bam>` → write a new BAM to stdout
//! with the records from `<in.bam>` and the header from `<new.hdr.sam>`.
//!
//! Not yet supported: `--in-place` (CRAM rewrite) and BGZF block-level fast
//! paths.

use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufReader, Cursor, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

use htslib_rs::bam;
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf};

use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools reheader`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut in_place = false;
    let mut no_pg = false;
    let mut external_cmd: Option<String> = None;
    let mut positional: Vec<PathBuf> = Vec::new();

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-i" | "--in-place" => {
                in_place = true;
            }
            "-P" | "--no-PG" => {
                no_pg = true;
            }
            "-c" | "--command" => {
                external_cmd = iter.next().and_then(|a| a.to_str()).map(|s| s.to_string());
            }
            "-h" | "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("reheader", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => positional.push(PathBuf::from(arg)),
        }
    }

    let input_path = if external_cmd.is_some() {
        if positional.len() != 1 {
            let _ = print_usage();
            return ExitCode::from(1);
        }
        &positional[0]
    } else {
        if positional.len() != 2 {
            let _ = print_usage();
            return ExitCode::from(1);
        }
        &positional[1]
    };

    if args.iter().any(|a| a == "-c" || a == "--command") && external_cmd.is_none() {
        print_error("reheader", "missing value for -c/--command");
        let _ = print_usage();
        return ExitCode::from(1);
    }

    if external_cmd.is_none() {
        let header_path = &positional[0];
        if header_path.as_os_str() != "-" && !header_path.exists() {
            print_hts_open_missing(header_path);
            print_error(
                "reheader",
                format!(
                    "fail to read the header from '{}': No such file or directory",
                    header_path.display()
                ),
            );
            return ExitCode::from(1);
        }
    }

    if input_path.as_os_str() != "-" && !input_path.exists() {
        print_hts_open_missing(input_path);
        print_error(
            "reheader",
            format!(
                "fail to open file '{}': No such file or directory",
                input_path.display()
            ),
        );
        return ExitCode::from(1);
    }

    let format = match sam_io::sam_open_format(input_path) {
        Ok(f) => f,
        Err(e) => {
            print_error("reheader", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Bam | Exact::Cram) {
        print_error(
            "reheader",
            format!("input file '{}' must be BAM or CRAM", input_path.display()),
        );
        return ExitCode::from(1);
    }

    let result = if in_place {
        if format.exact != Exact::Cram {
            print_error("reheader", "--in-place is only supported for CRAM input");
            return ExitCode::from(1);
        }
        if let Some(command) = external_cmd.as_deref() {
            run_reheader_command_in_place(command, input_path, !no_pg, args)
        } else {
            run_reheader_in_place(&positional[0], input_path, !no_pg, args)
        }
    } else if let Some(command) = external_cmd.as_deref() {
        run_reheader_command(command, input_path, !no_pg, args)
    } else {
        run_reheader(&positional[0], input_path, !no_pg, args)
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("reheader", "reheader failed", &e);
            ExitCode::from(1)
        }
    }
}

fn run_reheader(
    new_header_path: &Path,
    input_bam: &Path,
    add_pg: bool,
    argv: &[OsString],
) -> io::Result<()> {
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    run_reheader_to_writer(new_header_path, input_bam, add_pg, argv, &mut handle)
}

fn run_reheader_command(
    command: &str,
    input_bam: &Path,
    add_pg: bool,
    argv: &[OsString],
) -> io::Result<()> {
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    run_reheader_command_to_writer(command, input_bam, add_pg, argv, &mut handle)
}

fn run_reheader_in_place(
    new_header_path: &Path,
    input: &Path,
    add_pg: bool,
    argv: &[OsString],
) -> io::Result<()> {
    if sam_io::sam_open_format(input)?.exact == Exact::Cram {
        let header_text = read_header_text_from_path(new_header_path, add_pg, argv)?;
        return rewrite_cram_sam_text_in_place(&header_text, input);
    }
    let new_header = read_header_from_path(new_header_path)?;
    rewrite_cram_in_place(new_header, input, add_pg, argv)
}

fn run_reheader_command_in_place(
    command: &str,
    input: &Path,
    add_pg: bool,
    argv: &[OsString],
) -> io::Result<()> {
    let header_text = crate::header_text::read_raw_header_text_with_format(input, Exact::Cram)?;
    let new_header_text = read_header_text_from_command(command, &header_text, add_pg, argv)?;
    rewrite_cram_sam_text_in_place(&new_header_text, input)
}

fn rewrite_cram_sam_text_in_place(new_header_text: &str, input: &Path) -> io::Result<()> {
    let tmp = temporary_sibling_path(input, "reheader", "sam");
    {
        let output = File::create(&tmp)?;
        rewrite_cram_as_sam_with_header_text(new_header_text, input, output)?;
    }
    fs::rename(&tmp, input).inspect_err(|_| {
        let _ = fs::remove_file(&tmp);
    })
}

fn read_header_text_from_path(path: &Path, add_pg: bool, argv: &[OsString]) -> io::Result<String> {
    let text = fs::read_to_string(path)?;
    if add_pg {
        crate::pg::add_samtools_pg(&text, argv)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    } else {
        Ok(text)
    }
}

fn read_header_text_from_command(
    command: &str,
    header_text: &str,
    add_pg: bool,
    argv: &[OsString],
) -> io::Result<String> {
    let text = run_header_command(command, header_text)?;
    if add_pg {
        crate::pg::add_samtools_pg(&text, argv)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    } else {
        Ok(text)
    }
}

fn run_header_command(command: &str, header_text: &str) -> io::Result<String> {
    let output = run_header_command_bytes(command, header_text)?;
    String::from_utf8(output).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn run_header_command_bytes(command: &str, header_text: &str) -> io::Result<Vec<u8>> {
    let mut child = Command::new("sh")
        .arg("-c")
        .arg(command)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;

    {
        let mut stdin = child.stdin.take().ok_or_else(|| {
            io::Error::new(io::ErrorKind::BrokenPipe, "failed to open command stdin")
        })?;
        stdin.write_all(header_text.as_bytes())?;
    }

    let output = child.wait_with_output()?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let message = stderr.trim();
        return Err(io::Error::other(if message.is_empty() {
            format!("header command exited with status {}", output.status)
        } else {
            format!(
                "header command exited with status {}: {message}",
                output.status
            )
        }));
    }

    Ok(output.stdout)
}

pub(crate) fn run_reheader_to_writer<W: Write>(
    new_header_path: &Path,
    input: &Path,
    add_pg: bool,
    argv: &[OsString],
    output: W,
) -> io::Result<()> {
    if sam_io::sam_open_format(input)?.exact == Exact::Cram {
        let header_text = read_header_text_from_path(new_header_path, add_pg, argv)?;
        return rewrite_cram_as_sam_with_header_text(&header_text, input, output);
    }
    let new_header = read_header_from_path(new_header_path)?;
    rewrite_with_header(new_header, input, add_pg, argv, output)
}

pub(crate) fn run_reheader_command_to_writer<W: Write>(
    command: &str,
    input: &Path,
    add_pg: bool,
    argv: &[OsString],
    output: W,
) -> io::Result<()> {
    let exact = sam_io::sam_open_format(input)?.exact;
    let header_text = crate::header_text::read_raw_header_text_with_format(input, exact)?;
    if exact == Exact::Cram {
        let new_header_text = read_header_text_from_command(command, &header_text, add_pg, argv)?;
        return rewrite_cram_as_sam_with_header_text(&new_header_text, input, output);
    }
    let new_header = read_header_from_command(command, &header_text)?;
    rewrite_with_header(new_header, input, add_pg, argv, output)
}

fn read_header_from_path(path: &Path) -> io::Result<sam::Header> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(path)?));
    reader.read_header()
}

fn read_header_from_command(command: &str, header_text: &str) -> io::Result<sam::Header> {
    let mut reader =
        sam::io::Reader::new(Cursor::new(run_header_command_bytes(command, header_text)?));
    reader.read_header()
}

fn rewrite_bam_with_header<W: Write>(
    new_header: sam::Header,
    input_bam: &Path,
    add_pg: bool,
    argv: &[OsString],
    output: W,
) -> io::Result<()> {
    let new_header = if add_pg {
        crate::pg::add_samtools_pg_to_header(&new_header, argv)?
    } else {
        new_header
    };

    let mut input = bam::io::Reader::new(File::open(input_bam)?);
    let input_header = input.read_header()?;

    let mut writer = bam::io::Writer::new(output);
    writer.write_header(&new_header)?;

    let mut record = bam::Record::default();
    loop {
        let n = input.read_record(&mut record)?;
        if n == 0 {
            break;
        }
        writer.write_record(&input_header, &record)?;
    }
    Ok(())
}

fn rewrite_cram_in_place(
    new_header: sam::Header,
    input: &Path,
    add_pg: bool,
    argv: &[OsString],
) -> io::Result<()> {
    let tmp = temporary_sibling_path(input, "reheader", "cram");
    {
        let output = File::create(&tmp)?;
        rewrite_cram_with_header(new_header, input, add_pg, argv, output)?;
    }
    fs::rename(&tmp, input).inspect_err(|_| {
        let _ = fs::remove_file(&tmp);
    })
}

fn rewrite_with_header<W: Write>(
    new_header: sam::Header,
    input: &Path,
    add_pg: bool,
    argv: &[OsString],
    output: W,
) -> io::Result<()> {
    match sam_io::sam_open_format(input)?.exact {
        Exact::Sam => rewrite_sam_with_header(new_header, input, add_pg, argv, output),
        Exact::Bam => rewrite_bam_with_header(new_header, input, add_pg, argv, output),
        Exact::Cram => rewrite_cram_with_header(new_header, input, add_pg, argv, output),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "only SAM, BAM, and CRAM input are currently supported",
        )),
    }
}

fn rewrite_cram_with_header<W: Write>(
    new_header: sam::Header,
    input_cram: &Path,
    add_pg: bool,
    argv: &[OsString],
    output: W,
) -> io::Result<()> {
    let new_header = if add_pg {
        crate::pg::add_samtools_pg_to_header(&new_header, argv)?
    } else {
        new_header
    };

    let reference = cram_reference_for_input(input_cram)?;
    let mut merged = Vec::new();
    crate::sam_render::write_header(&mut merged, &new_header)?;
    let records = htslib_rs::alignment_compat::query_cram_records_all_from_path_with_reference(
        input_cram,
        reference.path(),
    )?;
    for record in records {
        crate::sam_render::write_record(&mut merged, &new_header, &record)?;
    }

    let mut reader = sam::io::Reader::new(BufReader::new(Cursor::new(merged)));
    htslib_rs::alignment_compat::write_cram_from_sam_reader_with_reference(
        &mut reader,
        reference.path(),
        output,
    )?;
    Ok(())
}

fn rewrite_cram_as_sam_with_header_text<W: Write>(
    new_header_text: &str,
    input_cram: &Path,
    mut output: W,
) -> io::Result<()> {
    let reference = cram_reference_for_input(input_cram)?;
    let old_text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            input_cram,
            reference.path(),
            None,
        )?;
    output.write_all(new_header_text.as_bytes())?;
    for line in strip_header_lines(old_text.as_bytes()).split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }
        let normalized = normalize_reheader_cram_sam_line(line);
        output.write_all(normalized.as_bytes())?;
        output.write_all(b"\n")?;
    }
    Ok(())
}

fn normalize_reheader_cram_sam_line(line: &[u8]) -> String {
    let mut fields: Vec<String> = line
        .split(|&b| b == b'\t')
        .map(|field| String::from_utf8_lossy(field).into_owned())
        .collect();
    if fields.len() < 11 {
        return String::from_utf8_lossy(line).into_owned();
    }

    if fields[2] == "*" {
        fields[4] = "0".to_string();
    }
    if fields[6] != "=" && fields[6] != "*" && fields[2] != fields[6] {
        fields[8] = "0".to_string();
    }

    let aux = fields.split_off(11);
    let mut md = None;
    let mut nm = None;
    let mut rg = None;
    let mut rest = Vec::with_capacity(aux.len());
    for field in aux {
        if field.starts_with("MD:") {
            md = Some(field);
        } else if field.starts_with("NM:") {
            nm = Some(field);
        } else if field.starts_with("RG:") {
            rg = Some(field);
        } else {
            rest.push(field);
        }
    }
    fields.extend(rest);
    if let Some(field) = md {
        fields.push(field);
    }
    if let Some(field) = nm {
        fields.push(field);
    }
    if let Some(field) = rg {
        fields.push(field);
    }

    crate::sam_render::fix_sam_aux_floats(&fields.join("\t"))
}

fn rewrite_sam_with_header<W: Write>(
    new_header: sam::Header,
    input_sam: &Path,
    add_pg: bool,
    argv: &[OsString],
    mut output: W,
) -> io::Result<()> {
    let new_header = if add_pg {
        crate::pg::add_samtools_pg_to_header(&new_header, argv)?
    } else {
        new_header
    };

    let mut input = sam::io::Reader::new(BufReader::new(File::open(input_sam)?));
    let input_header = input.read_header()?;
    // Shared renderer: htslib `%g` float aux spelling.
    crate::sam_render::write_header(&mut output, &new_header)?;

    loop {
        let mut record = RecordBuf::default();
        if input.read_record_buf(&input_header, &mut record)? == 0 {
            break;
        }
        crate::sam_render::write_record(&mut output, &new_header, &record)?;
    }
    Ok(())
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

fn cram_reference_for_input(input: &Path) -> io::Result<ReferenceGuard> {
    if let Some(reference) = current_global_args().reference {
        return Ok(ReferenceGuard::new(reference));
    }

    reference_from_ref_path(input)?.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM input requires --reference or REF_PATH entries matching @SQ M5 tags",
        )
    })
}

fn reference_from_ref_path(input: &Path) -> io::Result<Option<ReferenceGuard>> {
    let Some(ref_path) = std::env::var_os("REF_PATH") else {
        return Ok(None);
    };
    let ref_path = ref_path.to_string_lossy();
    let header_text = crate::header_text::read_raw_header_text_with_format(input, Exact::Cram)?;
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
        if let Some(sequence) = read_ref_path_md5_sequence(&ref_path, md5)? {
            sequences.push((name.to_string(), sequence));
        } else if let Some(len) = len {
            sequences.push((name.to_string(), "N".repeat(len)));
        }
    }

    if sequences.is_empty() {
        return Ok(None);
    }

    let fasta = temporary_sibling_path(input, "ref-path", "fa");
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

fn temporary_sibling_path(input: &Path, stem: &str, ext: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let name = format!("samtools-rs-{stem}-{}-{nanos}.{ext}", std::process::id());
    input
        .parent()
        .map(|parent| parent.join(&name))
        .unwrap_or_else(|| PathBuf::from(name))
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools reheader [options] in.header.sam in.bam")?;
    writeln!(w, "Options:")?;
    writeln!(w, "  -P, --no-PG       do not add a @PG line")?;
    writeln!(
        w,
        "  -i, --in-place    edit the file in-place (CRAM only; TODO)"
    )?;
    writeln!(w, "  -c, --command CMD pipe existing header through CMD")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixtures_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("samtools")
            .join("test")
    }

    fn tmp_path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "samtools-rs-reheader-{}-{}",
            name,
            std::process::id()
        ))
    }

    #[test]
    fn reheader_adds_pg_by_default() {
        let fixtures = fixtures_dir();
        let header = fixtures.join("reheader").join("hdr.sam");
        let input = fixtures.join("checksum").join("chk1.bam");
        let output = tmp_path("pg.bam");
        let argv = [
            OsString::from("reheader"),
            header.clone().into(),
            input.clone().into(),
        ];

        {
            let out = File::create(&output).unwrap();
            run_reheader_to_writer(&header, &input, true, &argv, out).unwrap();
        }

        let header_text = crate::header_text::read_raw_header_text(&output).unwrap();
        // Upstream `@PG` field order is ID, PN, PP, VN, CL — `hdr.sam`
        // has `@PG ID:prog1`, so the new samtools entry chains `PP:prog1`
        // before VN/CL.
        assert!(header_text.contains("\tPN:samtools\tPP:prog1\tVN:"));
        assert!(header_text.contains("\tCL:reheader "));
        let _ = std::fs::remove_file(output);
    }

    #[test]
    fn reheader_no_pg_suppresses_pg() {
        let fixtures = fixtures_dir();
        let header = fixtures.join("reheader").join("hdr.sam");
        let input = fixtures.join("checksum").join("chk1.bam");
        let output = tmp_path("no-pg.bam");
        let argv = [
            OsString::from("reheader"),
            header.clone().into(),
            input.clone().into(),
        ];

        {
            let out = File::create(&output).unwrap();
            run_reheader_to_writer(&header, &input, false, &argv, out).unwrap();
        }

        let header_text = crate::header_text::read_raw_header_text(&output).unwrap();
        assert!(!header_text.contains("\tCL:reheader "));
        let _ = std::fs::remove_file(output);
    }

    #[test]
    fn reheader_command_filters_existing_header() {
        let fixtures = fixtures_dir();
        let input = fixtures.join("checksum").join("chk1.bam");
        let output = tmp_path("command.bam");
        let command = "cat; printf '@CO\\treheader command\\n'";
        let argv = [
            OsString::from("reheader"),
            OsString::from("-c"),
            OsString::from(command),
            input.clone().into(),
        ];

        {
            let out = File::create(&output).unwrap();
            run_reheader_command_to_writer(command, &input, true, &argv, out).unwrap();
        }

        let header_text = crate::header_text::read_raw_header_text(&output).unwrap();
        assert!(header_text.contains("@CO\treheader command"));
        assert!(header_text.contains("\tCL:reheader "));
        let _ = std::fs::remove_file(output);
    }

    #[test]
    fn reheader_sam_input_writes_sam_with_replacement_header() {
        let input = tmp_path("input.sam");
        let header = tmp_path("header.sam");
        let argv = [
            OsString::from("reheader"),
            header.clone().into(),
            input.clone().into(),
        ];
        std::fs::write(
            &input,
            "@HD\tVN:1.6\n@SQ\tSN:old\tLN:20\nr1\t0\told\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        )
        .unwrap();
        std::fs::write(
            &header,
            "@HD\tVN:1.6\n@SQ\tSN:old\tLN:20\n@CO\tnew header\n",
        )
        .unwrap();

        let mut output = Vec::new();
        run_reheader_to_writer(&header, &input, false, &argv, &mut output).unwrap();
        let text = String::from_utf8(output).unwrap();
        assert!(text.contains("@CO\tnew header\n"));
        assert!(text.contains("r1\t0\told\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n"));
        assert!(!text.contains("\tCL:reheader "));

        let _ = std::fs::remove_file(input);
        let _ = std::fs::remove_file(header);
    }
}
