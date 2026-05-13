//! `samtools fastq` / `samtools fasta` / `samtools bam2fq` — convert
//! SAM/BAM records to FASTQ or FASTA text.
//!
//! Mirrors `main_bam2fq` in `bam_fastq.c`. The initial Rust port supports
//! single-output mode (all reads written to stdout, `-o FILE`, or `-0 FILE`) and
//! basic paired-output split (`-1`/`-2`/`-s`).
//!
//! **Not yet supported:** exact name-grouped singleton/other routing, barcode
//! tag handling (`-T`, `-i`), index files (`--i1`/`--i2`).

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use htslib_rs::format::{Exact, detect_path};

use crate::diagnostics::{print_error, print_error_errno};

/// Entry point for `samtools fastq` / `samtools fasta` / `samtools bam2fq`.
pub fn main(args: &[OsString]) -> ExitCode {
    let sub_name = args.first().and_then(|a| a.to_str()).unwrap_or("fastq");
    let fasta_mode = sub_name == "fasta";

    let mut output: Option<PathBuf> = None;
    let mut other_output: Option<PathBuf> = None;
    let mut read1_output: Option<PathBuf> = None;
    let mut read2_output: Option<PathBuf> = None;
    let mut singleton_output: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let mut require_flags = 0u16;
    let mut exclude_flags = 0u16;
    let mut exclude_all_flags = 0u16;
    let mut aux_tags: Vec<[u8; 2]> = Vec::new();
    let mut append_read_number_override: Option<bool> = None;
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-o" => {
                output = iter.next().map(PathBuf::from);
            }
            "-0" => {
                other_output = iter.next().map(PathBuf::from);
            }
            "-1" => {
                read1_output = iter.next().map(PathBuf::from);
            }
            "-2" => {
                read2_output = iter.next().map(PathBuf::from);
            }
            "-s" => {
                singleton_output = iter.next().map(PathBuf::from);
            }
            "-f" => {
                require_flags = match parse_flag_arg(iter.next(), "-f", sub_name) {
                    Ok(flag) => flag,
                    Err(code) => return code,
                };
            }
            "-F" => {
                exclude_flags = match parse_flag_arg(iter.next(), "-F", sub_name) {
                    Ok(flag) => flag,
                    Err(code) => return code,
                };
            }
            "-G" => {
                exclude_all_flags = match parse_flag_arg(iter.next(), "-G", sub_name) {
                    Ok(flag) => flag,
                    Err(code) => return code,
                };
            }
            "-T" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error(sub_name, "missing value for -T");
                    return ExitCode::from(1);
                };
                aux_tags = parse_aux_tag_list(raw);
            }
            "-n" => {
                append_read_number_override = Some(false);
            }
            "-N" => {
                append_read_number_override = Some(true);
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "--help" => {
                let _ = print_usage(sub_name);
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error(
                    sub_name,
                    format!(
                        "option `{}` is not yet supported in samtools-rs {}",
                        s, sub_name
                    ),
                );
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                }
            }
        }
    }

    let stdin_input = input.as_ref().is_none_or(|path| path.as_os_str() == "-");

    let format = if stdin_input {
        None
    } else {
        let input = input.as_ref().expect("non-stdin input exists");
        match detect_path(input) {
            Ok(f) => Some(f),
            Err(e) => {
                print_error(
                    sub_name,
                    format!("failed to detect format of \"{}\": {}", input.display(), e),
                );
                return ExitCode::from(1);
            }
        }
    };

    let filtering = require_flags != 0 || exclude_flags != 0 || exclude_all_flags != 0;
    let split_mode = read1_output.is_some()
        || read2_output.is_some()
        || singleton_output.is_some()
        || other_output.is_some();
    let singleton_only = singleton_output.is_some()
        && read1_output.is_none()
        && read2_output.is_none()
        && other_output.is_none();
    let append_read_number = append_read_number_override.unwrap_or(!split_mode || singleton_only);

    if split_mode {
        let split = if stdin_input {
            let stdin = io::stdin().lock();
            let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(stdin));
            if fasta_mode {
                htslib_rs::alignment_compat::view_sam_as_fasta_split_text_from_reader_with_flag_filter_and_suffix(
                    &mut reader,
                    require_flags,
                    exclude_flags,
                    exclude_all_flags,
                    append_read_number,
                )
            } else {
                htslib_rs::alignment_compat::view_sam_as_fastq_split_text_from_reader_with_flag_filter_suffix_and_aux(
                    &mut reader,
                    require_flags,
                    exclude_flags,
                    exclude_all_flags,
                    append_read_number,
                    (!aux_tags.is_empty()).then_some(aux_tags.as_slice()),
                )
            }
        } else {
            let input = input.as_ref().expect("non-stdin input exists");
            match (format.expect("non-stdin format exists").exact, fasta_mode) {
            (Exact::Sam, false) => {
                htslib_rs::alignment_compat::view_sam_as_fastq_split_text_from_path_with_flag_filter_and_suffix(
                    input,
                    require_flags,
                    exclude_flags,
                    exclude_all_flags,
                    append_read_number,
                )
            }
            (Exact::Sam, true) => {
                htslib_rs::alignment_compat::view_sam_as_fasta_split_text_from_path_with_flag_filter_and_suffix(
                    input,
                    require_flags,
                    exclude_flags,
                    exclude_all_flags,
                    append_read_number,
                )
            }
            (Exact::Bam, false) => {
                htslib_rs::alignment_compat::view_bam_as_fastq_split_text_from_path_with_flag_filter_and_suffix(
                    input,
                    require_flags,
                    exclude_flags,
                    exclude_all_flags,
                    append_read_number,
                )
            }
            (Exact::Bam, true) => {
                htslib_rs::alignment_compat::view_bam_as_fasta_split_text_from_path_with_flag_filter_and_suffix(
                    input,
                    require_flags,
                    exclude_flags,
                    exclude_all_flags,
                    append_read_number,
                )
            }
            _ => {
                print_error(
                    sub_name,
                    "only SAM and BAM input are currently supported (CRAM TODO)",
                );
                return ExitCode::from(1);
            }
            }
        };

        let split = match split {
            Ok(split) => split,
            Err(e) => {
                print_error_errno(
                    sub_name,
                    format!(
                        "conversion failed for \"{}\"",
                        input
                            .as_ref()
                            .map(|path| path.display().to_string())
                            .unwrap_or_else(|| "-".to_string())
                    ),
                    &e,
                );
                return ExitCode::from(1);
            }
        };

        if let Some(path) = read1_output.as_ref()
            && let Err(e) = write_text_file(path, split.read1.as_bytes())
        {
            print_error_errno(sub_name, format!("open/write {}", path.display()), &e);
            return ExitCode::from(1);
        }
        if let Some(path) = read2_output.as_ref()
            && let Err(e) = write_text_file(path, split.read2.as_bytes())
        {
            print_error_errno(sub_name, format!("open/write {}", path.display()), &e);
            return ExitCode::from(1);
        }
        if let Some(path) = singleton_output.as_ref()
            && singleton_only
        {
            let mut all_singletons = String::new();
            all_singletons.push_str(&split.read1);
            all_singletons.push_str(&split.read2);
            all_singletons.push_str(&split.singleton);
            if let Err(e) = write_text_file(path, all_singletons.as_bytes()) {
                print_error_errno(sub_name, format!("open/write {}", path.display()), &e);
                return ExitCode::from(1);
            }
        } else if let Some(path) = singleton_output.as_ref()
            && let Err(e) = write_text_file(path, split.singleton.as_bytes())
        {
            print_error_errno(sub_name, format!("open/write {}", path.display()), &e);
            return ExitCode::from(1);
        }
        if let Some(path) = other_output.as_ref()
            && let Err(e) = write_text_file(path, split.singleton.as_bytes())
        {
            print_error_errno(sub_name, format!("open/write {}", path.display()), &e);
            return ExitCode::from(1);
        }

        return ExitCode::SUCCESS;
    }

    let text = if stdin_input {
        let stdin = io::stdin().lock();
        let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(stdin));
        if fasta_mode {
            htslib_rs::alignment_compat::view_sam_as_fasta_text_from_reader_with_flag_filter_and_suffix(
                &mut reader,
                require_flags,
                exclude_flags,
                exclude_all_flags,
                append_read_number,
            )
        } else {
            htslib_rs::alignment_compat::view_sam_as_fastq_text_from_reader_with_flag_filter_and_suffix(
                &mut reader,
                require_flags,
                exclude_flags,
                exclude_all_flags,
                append_read_number,
            )
        }
    } else {
        let input = input.as_ref().expect("non-stdin input exists");
        match (format.expect("non-stdin format exists").exact, fasta_mode, filtering) {
        (Exact::Sam, false, false) => {
            htslib_rs::alignment_compat::view_sam_as_fastq_text_from_path_with_limit_and_suffix(
                input,
                None,
                append_read_number,
            )
        }
        (Exact::Sam, false, true) => {
            htslib_rs::alignment_compat::view_sam_as_fastq_text_from_path_with_flag_filter_and_suffix(
                input,
                require_flags,
                exclude_flags,
                exclude_all_flags,
                append_read_number,
            )
        }
        (Exact::Sam, true, false) => {
            htslib_rs::alignment_compat::view_sam_as_fasta_text_from_path_with_limit_and_suffix(
                input,
                None,
                append_read_number,
            )
        }
        (Exact::Sam, true, true) => {
            htslib_rs::alignment_compat::view_sam_as_fasta_text_from_path_with_flag_filter_and_suffix(
                input,
                require_flags,
                exclude_flags,
                exclude_all_flags,
                append_read_number,
            )
        }
        (Exact::Bam, false, false) => {
            htslib_rs::alignment_compat::view_bam_as_fastq_text_from_path_with_limit_and_suffix(
                input,
                None,
                append_read_number,
            )
        }
        (Exact::Bam, false, true) => {
            htslib_rs::alignment_compat::view_bam_as_fastq_text_from_path_with_flag_filter_and_suffix(
                input,
                require_flags,
                exclude_flags,
                exclude_all_flags,
                append_read_number,
            )
        }
        (Exact::Bam, true, false) => {
            htslib_rs::alignment_compat::view_bam_as_fasta_text_from_path_with_limit_and_suffix(
                input,
                None,
                append_read_number,
            )
        }
        (Exact::Bam, true, true) => {
            htslib_rs::alignment_compat::view_bam_as_fasta_text_from_path_with_flag_filter_and_suffix(
                input,
                require_flags,
                exclude_flags,
                exclude_all_flags,
                append_read_number,
            )
        }
        _ => {
            print_error(
                sub_name,
                "only SAM and BAM input are currently supported (CRAM TODO)",
            );
            return ExitCode::from(1);
        }
        }
    };

    let text = match text {
        Ok(t) => t,
        Err(e) => {
            print_error_errno(
                sub_name,
                format!(
                    "conversion failed for \"{}\"",
                    input
                        .as_ref()
                        .map(|path| path.display().to_string())
                        .unwrap_or_else(|| "-".to_string())
                ),
                &e,
            );
            return ExitCode::from(1);
        }
    };

    let output = output.as_ref().or(other_output.as_ref());
    let mut out: Box<dyn Write> = match output {
        Some(p) => match File::create(p) {
            Ok(f) => Box::new(f),
            Err(e) => {
                print_error_errno(sub_name, "open -o output", &e);
                return ExitCode::from(1);
            }
        },
        None => Box::new(io::stdout().lock()),
    };
    if let Err(e) = out.write_all(text.as_bytes())
        && e.kind() != io::ErrorKind::BrokenPipe
    {
        print_error_errno(sub_name, "write output", &e);
        return ExitCode::from(1);
    }
    ExitCode::SUCCESS
}

fn write_text_file(path: &std::path::Path, text: &[u8]) -> io::Result<()> {
    let mut file = File::create(path)?;
    file.write_all(text)
}

fn parse_flag_arg(arg: Option<&OsString>, opt: &str, sub_name: &str) -> Result<u16, ExitCode> {
    let Some(raw) = arg.and_then(|a| a.to_str()) else {
        print_error(sub_name, format!("missing value for {}", opt));
        return Err(ExitCode::from(1));
    };
    match crate::bam_flag::str_to_flag(raw) {
        Some(flag) => u16::try_from(flag).map_err(|_| {
            print_error(sub_name, format!("Could not parse \"{}\"", raw));
            ExitCode::from(1)
        }),
        None => {
            print_error(sub_name, format!("Could not parse \"{}\"", raw));
            Err(ExitCode::from(1))
        }
    }
}

fn parse_aux_tag_list(raw: &str) -> Vec<[u8; 2]> {
    raw.split(',')
        .filter_map(|tag| {
            let bytes = tag.as_bytes();
            (bytes.len() == 2).then_some([bytes[0], bytes[1]])
        })
        .collect()
}

fn print_usage(sub: &str) -> io::Result<()> {
    let mut w = io::stderr().lock();
    let suffix = if sub == "fasta" { "fasta" } else { "fastq" };
    writeln!(
        w,
        "Usage: samtools {} [options] <in.bam>  > out.{}",
        sub, suffix
    )?;
    writeln!(w, "  -o FILE      write output to FILE (default stdout)")?;
    writeln!(
        w,
        "  -0 FILE      write all reads to FILE in single-output mode"
    )?;
    writeln!(w, "  -n           do not append /1 or /2 to read names")?;
    writeln!(w, "  -N           append /1 or /2 to read names")?;
    writeln!(w, "  -T TAGLIST   copy aux tags to FASTQ comments")?;
    writeln!(
        w,
        "  -f FLAG      only include reads with all FLAG bits set"
    )?;
    writeln!(w, "  -F FLAG      exclude reads with any FLAG bits set")?;
    writeln!(w, "  -G FLAG      exclude reads with all FLAG bits set")?;
    Ok(())
}
