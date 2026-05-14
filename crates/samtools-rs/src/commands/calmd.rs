//! `samtools calmd` / `samtools fillmd` — recalculate MD/NM tags and BAQ.
//!
//! Mirrors `bam_fillmd` in `bam_md.c`. The full upstream version recomputes
//! MD/NM tags against a reference and optionally applies BAQ (Base Alignment
//! Quality). The initial Rust port wraps the BAQ paths already available in
//! `htslib_rs::alignment_compat`:
//!  - default (`-r` not set): no-op pass-through (emits the input as SAM).
//!  - `-r`: recalculate BAQ (`recalculate_baq_from_sam_path`).
//!  - `-r -e`: apply existing BAQ to the quality strings.
//!  - `-E`: extended BAQ (`recalculate_extended_baq_from_sam_path`).
//!
//! **Pending:** MD/NM tag recomputation from reference, `-A` (always-apply),
//! `-d` (drop BAQ), `-C cap`, BAM/CRAM I/O.

use std::ffi::OsString;
use std::io::{self, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use htslib_rs::format::Exact;

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

/// Entry point for `samtools calmd` / `samtools fillmd`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut realn = false;
    let mut extended = false;
    let mut apply_existing = false;
    let mut reference: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-r" => realn = true,
            "-E" => extended = true,
            "-e" => apply_existing = true,
            "-A" | "-C" | "-n" | "-S" | "-b" | "-u" | "-N" | "-h" | "-q" | "-Q" | "--no-PG"
            | "-d" => {
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
    if format.exact != Exact::Sam {
        print_error(
            "calmd",
            "only SAM input is currently supported (BAM/CRAM TODO)",
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
        // No-op mode: just copy the input through.
        htslib_rs::alignment_compat::view_sam_text_from_path_with_limit(&input, None)
    };

    let text = match text {
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
