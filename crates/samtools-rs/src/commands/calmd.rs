//! `samtools calmd` / `samtools fillmd` — recalculate MD/NM tags and BAQ.
//!
//! Mirrors `bam_fillmd` in `bam_md.c`. The full upstream version recomputes
//! MD/NM tags against a reference and optionally applies BAQ (Base Alignment
//! Quality). The initial Rust port wraps the BAQ paths already available in
//! `htslib_rs::alignment_compat`:
//!  - default (`-r` not set): recomputes MD/NM tags when a reference is
//!    supplied, otherwise emits the input as SAM. SAM, BAM, and reference-
//!    backed CRAM input are supported for this record-text path.
//!  - `-r`: recalculate BAQ (`recalculate_baq_from_sam_path`).
//!  - `-r -e`: apply existing BAQ to the quality strings.
//!  - `-E`: extended BAQ (`recalculate_extended_baq_from_sam_path`).
//!  - `-d`: drop existing `BQ` tags from the SAM-text output.
//!
//! **Pending:** `-A` (always-apply), `-C cap`, BAM/CRAM output, and
//! BAM/CRAM BAQ helper paths.

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::format::Exact;

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

/// Entry point for `samtools calmd` / `samtools fillmd`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut realn = false;
    let mut extended = false;
    let mut apply_existing = false;
    let mut drop_baq = false;
    let mut reference: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut no_pg = false;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-r" => realn = true,
            "-E" => extended = true,
            "-e" => apply_existing = true,
            "--no-PG" => no_pg = true,
            "-d" => {
                drop_baq = true;
            }
            "-A" | "-C" | "-n" | "-S" | "-b" | "-u" | "-N" | "-h" | "-q" | "-Q" => {
                if matches!(s, "-C" | "-n") {
                    let _ = iter.next();
                }
            }
            "-T" | "--reference" => {
                reference = iter.next().map(PathBuf::from);
            }
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("calmd", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                } else if reference.is_none() {
                    reference = Some(PathBuf::from(arg));
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
            print_error("calmd", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam | Exact::Cram) {
        print_error("calmd", "only SAM, BAM, and CRAM input are supported");
        return ExitCode::from(1);
    }
    if format.exact == Exact::Cram && reference.is_none() {
        print_error("calmd", "CRAM input requires -T/--reference or ref.fa");
        return ExitCode::from(1);
    }
    if format.exact != Exact::Sam && (realn || apply_existing || extended) {
        print_error(
            "calmd",
            "BAQ recalculation/application is currently supported for SAM input only",
        );
        return ExitCode::from(1);
    }

    let text = if realn {
        let Some(reference) = reference.as_ref() else {
            print_error("calmd", "-r/--reference required for BAQ recalculation");
            return ExitCode::from(1);
        };
        if extended {
            htslib_rs::alignment_compat::recalculate_extended_baq_from_sam_path(&input, reference)
        } else if apply_existing {
            htslib_rs::alignment_compat::recalculate_and_apply_baq_from_sam_path(&input, reference)
        } else {
            htslib_rs::alignment_compat::recalculate_baq_from_sam_path(&input, reference)
        }
    } else if apply_existing {
        htslib_rs::alignment_compat::apply_existing_baq_from_sam_path(&input)
    } else {
        input_as_sam_text(&input, format.exact, reference.as_deref())
    };

    let mut text = match text {
        Ok(t) => t,
        Err(e) => {
            print_error_errno(
                "calmd",
                format!("calmd failed for \"{}\"", input.display()),
                &e,
            );
            return ExitCode::from(1);
        }
    };

    if let Some(reference) = reference.as_ref() {
        match recalculate_md_nm_sam_text(&text, reference) {
            Ok(modified) => text = modified,
            Err(e) => {
                print_error_errno("calmd", "recalculate MD/NM", &e);
                return ExitCode::from(1);
            }
        }
    }
    if drop_baq {
        text = remove_aux_tag_sam_text(&text, "BQ");
    }

    let text = if no_pg {
        text
    } else {
        match inject_pg_into_sam_text(&text, args) {
            Ok(modified) => modified,
            Err(e) => {
                print_error_errno("calmd", "inject @PG line", &e);
                return ExitCode::from(1);
            }
        }
    };

    let mut out = match sam_io::open_text_output(output.as_deref()) {
        Ok(out) => out,
        Err(e) => {
            print_error_errno("calmd", "open -o output", &e);
            return ExitCode::from(1);
        }
    };
    if let Err(e) = sam_io::write_all_and_close(&mut out, text.as_bytes())
        && e.kind() != io::ErrorKind::BrokenPipe
    {
        print_error_errno("calmd", "write output", &e);
        return ExitCode::from(1);
    }
    ExitCode::SUCCESS
}

fn recalculate_md_nm_sam_text(text: &str, reference: &Path) -> io::Result<String> {
    let references = read_fasta(reference)?;
    let mut out = String::with_capacity(text.len());

    for line in text.split_inclusive('\n') {
        let (line_body, newline) = line
            .strip_suffix('\n')
            .map(|body| (body, "\n"))
            .unwrap_or((line, ""));

        if line_body.starts_with('@') || line_body.is_empty() {
            out.push_str(line_body);
            out.push_str(newline);
            continue;
        }

        let fields: Vec<&str> = line_body.split('\t').collect();
        if fields.len() < 11 || fields[2] == "*" || fields[5] == "*" {
            out.push_str(line_body);
            out.push_str(newline);
            continue;
        }

        let Some(reference_sequence) = references.get(fields[2]) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("reference sequence {} not found", fields[2]),
            ));
        };

        let start = fields[3].parse::<usize>().map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid POS {}: {e}", fields[3]),
            )
        })?;
        if start == 0 {
            out.push_str(line_body);
            out.push_str(newline);
            continue;
        }

        let (md, nm) = calculate_md_nm(fields[5], fields[9], reference_sequence, start - 1)?;
        write_record_with_md_nm(&mut out, &fields, &md, nm);
        out.push_str(newline);
    }

    Ok(out)
}

fn input_as_sam_text(input: &Path, exact: Exact, reference: Option<&Path>) -> io::Result<String> {
    match exact {
        Exact::Sam => htslib_rs::alignment_compat::view_sam_text_from_path_with_limit(input, None),
        Exact::Bam => {
            htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(input, None)
        }
        Exact::Cram => {
            let reference = reference.ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "CRAM input requires reference")
            })?;
            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                input, reference, None,
            )
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "unsupported alignment format",
        )),
    }
}

fn remove_aux_tag_sam_text(text: &str, tag: &str) -> String {
    let prefix = format!("{tag}:");
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        let (line_body, newline) = line
            .strip_suffix('\n')
            .map(|body| (body, "\n"))
            .unwrap_or((line, ""));
        if line_body.starts_with('@') || line_body.is_empty() {
            out.push_str(line_body);
            out.push_str(newline);
            continue;
        }

        let mut first = true;
        for field in line_body.split('\t') {
            if field.starts_with(&prefix) {
                continue;
            }
            if !first {
                out.push('\t');
            }
            first = false;
            out.push_str(field);
        }
        out.push_str(newline);
    }
    out
}

fn read_fasta(path: &Path) -> io::Result<HashMap<String, Vec<u8>>> {
    let text = fs::read_to_string(path)?;
    let mut references = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = current_name.replace(
                rest.split_whitespace()
                    .next()
                    .unwrap_or_default()
                    .to_string(),
            ) {
                references.insert(name, std::mem::take(&mut current_seq));
            }
        } else if !line.starts_with(';') {
            current_seq.extend(line.bytes().filter(|b| !b.is_ascii_whitespace()).map(|b| {
                if b.is_ascii_lowercase() {
                    b.to_ascii_uppercase()
                } else {
                    b
                }
            }));
        }
    }

    if let Some(name) = current_name {
        references.insert(name, current_seq);
    }

    Ok(references)
}

fn calculate_md_nm(
    cigar: &str,
    sequence: &str,
    reference: &[u8],
    start: usize,
) -> io::Result<(String, usize)> {
    let read = sequence.as_bytes();
    let mut read_i = 0usize;
    let mut ref_i = start;
    let mut nm = 0usize;
    let mut md = String::new();
    let mut match_count = 0usize;
    let mut n = 0usize;

    for b in cigar.bytes() {
        if b.is_ascii_digit() {
            n = n * 10 + usize::from(b - b'0');
            continue;
        }
        if n == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid CIGAR operation in {cigar}"),
            ));
        }

        match b {
            b'M' | b'=' | b'X' => {
                for _ in 0..n {
                    let read_base = read.get(read_i).copied().ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "CIGAR consumes past read")
                    })?;
                    let ref_base = reference.get(ref_i).copied().ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "CIGAR consumes past reference")
                    })?;
                    if bases_match(read_base, ref_base) {
                        match_count += 1;
                    } else {
                        md.push_str(&match_count.to_string());
                        match_count = 0;
                        md.push(char::from(ref_base.to_ascii_uppercase()));
                        nm += 1;
                    }
                    read_i += 1;
                    ref_i += 1;
                }
            }
            b'I' => {
                read_i += n;
                nm += n;
                if read_i > read.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CIGAR consumes past read",
                    ));
                }
            }
            b'D' => {
                md.push_str(&match_count.to_string());
                match_count = 0;
                md.push('^');
                for _ in 0..n {
                    let ref_base = reference.get(ref_i).copied().ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "CIGAR consumes past reference")
                    })?;
                    md.push(char::from(ref_base.to_ascii_uppercase()));
                    ref_i += 1;
                }
                nm += n;
            }
            b'N' => {
                ref_i += n;
                if ref_i > reference.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CIGAR consumes past reference",
                    ));
                }
            }
            b'S' => {
                read_i += n;
                if read_i > read.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CIGAR consumes past read",
                    ));
                }
            }
            b'H' | b'P' => {}
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("unsupported CIGAR operation {}", char::from(b)),
                ));
            }
        }
        n = 0;
    }

    if n != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("trailing CIGAR length in {cigar}"),
        ));
    }

    md.push_str(&match_count.to_string());
    Ok((md, nm))
}

fn bases_match(read_base: u8, ref_base: u8) -> bool {
    read_base.eq_ignore_ascii_case(&ref_base)
}

fn write_record_with_md_nm(out: &mut String, fields: &[&str], md: &str, nm: usize) {
    for (i, field) in fields.iter().enumerate() {
        if i >= 11 && (field.starts_with("MD:Z:") || field.starts_with("NM:i:")) {
            continue;
        }
        if !out.is_empty() && !out.ends_with('\n') {
            out.push('\t');
        }
        out.push_str(field);
    }
    out.push('\t');
    out.push_str("NM:i:");
    out.push_str(&nm.to_string());
    out.push('\t');
    out.push_str("MD:Z:");
    out.push_str(md);
}

/// Inserts samtools' `@PG` chain entry into a SAM text blob. Splits the
/// header (lines starting with `@`) from the body, applies the shared
/// `pg::add_samtools_pg` helper to the header, then rejoins.
fn inject_pg_into_sam_text(text: &str, argv: &[OsString]) -> io::Result<String> {
    let mut header_end = 0;
    for line in text.split_inclusive('\n') {
        if line.starts_with('@') {
            header_end += line.len();
        } else {
            break;
        }
    }
    let header_slice = &text[..header_end];
    let body_slice = &text[header_end..];
    let new_header = crate::pg::add_samtools_pg(header_slice, argv)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut combined = String::with_capacity(new_header.len() + body_slice.len());
    combined.push_str(&new_header);
    combined.push_str(body_slice);
    Ok(combined)
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools calmd [options] <in.sam> [ref.fa]")?;
    writeln!(w, "  -r          recalculate BAQ (requires --reference)")?;
    writeln!(w, "  -E          extended BAQ (with -r)")?;
    writeln!(w, "  -e          apply existing BAQ (BQ:Z) to qualities")?;
    writeln!(w, "  -T FILE     reference FASTA")?;
    writeln!(w, "  -o FILE     output FILE")?;
    Ok(())
}
