//! `samtools reheader` — replace the header of an alignment file.
//!
//! Mirrors `main_reheader` in `bam_reheader.c`. The basic mode is
//! `samtools reheader <new.hdr.sam> <in.bam>` → write a new BAM to stdout
//! with the records from `<in.bam>` and the header from `<new.hdr.sam>`.
//! SAM input is also supported and writes SAM output.
//!
//! Not yet supported: `--in-place` (CRAM rewrite) and BGZF block-level fast
//! paths.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Cursor, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Stdio};

use htslib_rs::bam;
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf};

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

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

    if in_place {
        print_error(
            "reheader",
            "the `--in-place` mode is not yet supported (CRAM only)",
        );
        return ExitCode::from(1);
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

    let format = match sam_io::sam_open_format(input_path) {
        Ok(f) => f,
        Err(e) => {
            print_error("reheader", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam) {
        print_error(
            "reheader",
            "only SAM and BAM input are currently supported (CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    let result = if let Some(command) = external_cmd.as_deref() {
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

pub(crate) fn run_reheader_to_writer<W: Write>(
    new_header_path: &Path,
    input: &Path,
    add_pg: bool,
    argv: &[OsString],
    output: W,
) -> io::Result<()> {
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
    let new_header = read_header_from_command(command, &header_text)?;
    rewrite_with_header(new_header, input, add_pg, argv, output)
}

fn read_header_from_path(path: &Path) -> io::Result<sam::Header> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(path)?));
    reader.read_header()
}

fn read_header_from_command(command: &str, header_text: &str) -> io::Result<sam::Header> {
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

    let mut reader = sam::io::Reader::new(Cursor::new(output.stdout));
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
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "only SAM and BAM input are currently supported",
        )),
    }
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
