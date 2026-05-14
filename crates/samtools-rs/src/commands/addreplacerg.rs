//! `samtools addreplacerg` — add or replace `@RG` lines and `RG:Z:` aux tags.
//!
//! Mirrors `main_addreplacerg` in `bam_addrprg.c`. Initial Rust port works on
//! **SAM input → SAM output** by streaming lines and rewriting them as text:
//!  - `-r '@RG\tID:foo\tSM:bar'` — full `@RG` line spec; merged into header.
//!  - `-r 'ID:foo'` — incremental tag form (one tag per `-r`); combined into
//!    a single `@RG` line.
//!  - `-R ID` — set every record's `RG:Z` to this existing ID.
//!  - `-m overwrite_all|orphan_only|orphan_first` — how to handle existing
//!    `RG:Z` tags. `orphan_only` is the default; `overwrite_all` always
//!    replaces.
//!  - `-O sam` — only SAM output is currently supported.
//!  - `--no-PG` — silently accepted (no `@PG` is added by this port).
//!  - `-w` — accepted but currently a no-op (editing-only mode).
//!
//! **Pending:** BAM/CRAM input/output, paired-end mate update, full
//! `orphan_first` semantics.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use flate2::read::MultiGzDecoder;
use htslib_rs::format::Exact;

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

#[derive(Clone, Copy)]
enum Mode {
    OrphanOnly,
    OverwriteAll,
}

/// Entry point for `samtools addreplacerg`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut tag_pieces: Vec<String> = Vec::new();
    let mut replace_id: Option<String> = None;
    let mut mode = Mode::OrphanOnly;
    let mut output: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let mut output_fmt = "sam".to_string();
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
                    _ => Mode::OrphanOnly,
                };
            }
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "-O" => {
                output_fmt = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .unwrap_or("sam")
                    .to_lowercase();
            }
            "-w" | "--no-PG" | "-@" | "--threads" => {
                if matches!(s, "-@" | "--threads") {
                    let _ = iter.next();
                }
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

    if output_fmt != "sam" {
        print_error(
            "addreplacerg",
            "only -O sam is currently supported in samtools-rs",
        );
        return ExitCode::from(1);
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
    if format.exact != Exact::Sam {
        print_error(
            "addreplacerg",
            "only SAM input is currently supported (BAM/CRAM TODO)",
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
        _ => None,
    };

    let Some(rg_id) = rg_id else {
        print_error(
            "addreplacerg",
            "an @RG ID is required (use `-r 'ID:...'` or `-R ID`)",
        );
        return ExitCode::from(1);
    };

    let mut writer = match sam_io::open_text_output(output.as_deref()) {
        Ok(writer) => writer,
        Err(e) => {
            print_error_errno("addreplacerg", "open -o output", &e);
            return ExitCode::from(1);
        }
    };

    if let Err(e) = rewrite_sam(&input, &mut writer, rg_line.as_deref(), &rg_id, mode) {
        if e.kind() == io::ErrorKind::BrokenPipe {
            return ExitCode::SUCCESS;
        }
        print_error_errno("addreplacerg", "rewrite failed", &e);
        return ExitCode::from(1);
    }
    match sam_io::check_sam_close(&mut writer) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("addreplacerg", "close output", &e);
            ExitCode::from(1)
        }
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

fn extract_id(rg_line: &str) -> Option<String> {
    for field in rg_line.split('\t') {
        if let Some(rest) = field.strip_prefix("ID:") {
            return Some(rest.to_string());
        }
    }
    None
}

fn rewrite_sam(
    path: &Path,
    out: &mut dyn Write,
    rg_line: Option<&str>,
    rg_id: &str,
    mode: Mode,
) -> io::Result<()> {
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
    let mut header_emitted = false;
    let mut wrote_rg_line = false;
    let rg_tag = format!("RG:Z:{}", rg_id);

    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with('@') {
            // Suppress duplicate @RG lines with the same ID; otherwise pass through.
            if line.starts_with("@RG")
                && let Some(existing_id) = extract_id(line.trim_end())
                && existing_id == rg_id
            {
                continue;
            }
            out.write_all(line.as_bytes())?;
            header_emitted = true;
        } else {
            // First record line — flush our new @RG line ahead of it if needed.
            if !wrote_rg_line {
                if let Some(rg) = rg_line {
                    if !header_emitted {
                        // Nothing was emitted yet; add a minimal HD line so
                        // downstream SAM parsers are happy.
                        out.write_all(b"@HD\tVN:1.6\n")?;
                    }
                    out.write_all(rg.as_bytes())?;
                    if !rg.ends_with('\n') {
                        out.write_all(b"\n")?;
                    }
                }
                wrote_rg_line = true;
            }
            let stripped = line.trim_end_matches(&['\r', '\n'][..]);
            let new = rewrite_record(stripped, &rg_tag, mode);
            out.write_all(new.as_bytes())?;
            out.write_all(b"\n")?;
        }
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

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(
        w,
        "Usage: samtools addreplacerg [options] -r 'tag-spec' <in.sam>"
    )?;
    writeln!(
        w,
        "  -r SPEC       @RG line or 'KEY:VALUE' tag (repeatable)"
    )?;
    writeln!(w, "  -R ID         existing @RG ID to apply")?;
    writeln!(
        w,
        "  -m MODE       overwrite_all | orphan_only [orphan_only]"
    )?;
    writeln!(w, "  -o FILE       output FILE (default stdout)")?;
    writeln!(w, "  -O FMT        output format (only 'sam' supported)")?;
    Ok(())
}
