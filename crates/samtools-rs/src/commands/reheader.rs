//! `samtools reheader` — replace the header of a BAM file.
//!
//! Mirrors `main_reheader` in `bam_reheader.c`. The basic mode is
//! `samtools reheader <new.hdr.sam> <in.bam>` → write a new BAM to stdout
//! with the records from `<in.bam>` and the header from `<new.hdr.sam>`.
//!
//! Not yet supported: `--in-place` (CRAM rewrite), `-c <command>` (filter
//! existing header through a shell pipe), and BGZF block-level fast paths.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::format::{Exact, detect_path};
use htslib_rs::sam;

use crate::diagnostics::{print_error, print_error_errno};

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

    let _ = no_pg;

    if external_cmd.is_some() {
        print_error(
            "reheader",
            "the `-c <command>` external-header mode is not yet supported",
        );
        return ExitCode::from(1);
    }

    if in_place {
        print_error(
            "reheader",
            "the `--in-place` mode is not yet supported (CRAM only)",
        );
        return ExitCode::from(1);
    }

    if positional.len() != 2 {
        let _ = print_usage();
        return ExitCode::from(1);
    }
    let new_header_path = &positional[0];
    let input_bam_path = &positional[1];

    let format = match detect_path(input_bam_path) {
        Ok(f) => f,
        Err(e) => {
            print_error(
                "reheader",
                format!(
                    "failed to detect format of \"{}\": {}",
                    input_bam_path.display(),
                    e
                ),
            );
            return ExitCode::from(1);
        }
    };
    if format.exact != Exact::Bam {
        print_error(
            "reheader",
            "only BAM input is currently supported (CRAM/SAM TODO)",
        );
        return ExitCode::from(1);
    }

    match run_reheader(new_header_path, input_bam_path) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("reheader", "reheader failed", &e);
            ExitCode::from(1)
        }
    }
}

fn run_reheader(new_header_path: &Path, input_bam: &Path) -> io::Result<()> {
    let new_header: sam::Header = {
        let mut reader = sam::io::Reader::new(BufReader::new(File::open(new_header_path)?));
        reader.read_header()?
    };

    let mut input = bam::io::Reader::new(File::open(input_bam)?);
    let input_header = input.read_header()?;

    let mut writer = bam::io::Writer::new(io::stdout());
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

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools reheader [options] in.header.sam in.bam")?;
    writeln!(w, "Options:")?;
    writeln!(w, "  -P, --no-PG       do not add a @PG line")?;
    writeln!(
        w,
        "  -i, --in-place    edit the file in-place (CRAM only; TODO)"
    )?;
    writeln!(
        w,
        "  -c, --command CMD pipe existing header through CMD (TODO)"
    )?;
    Ok(())
}
