//! `samtools split` — split a BAM file by `@RG` ID.
//!
//! Mirrors `main_split` in `bam_split.c`. The initial Rust port supports:
//!  - Input BAM only (SAM/CRAM TODO).
//!  - Output template via `-f` with `%!` (RG ID), `%#` (RG index), and
//!    `%.` (extension); other placeholders are passed through verbatim.
//!  - `-u <path>` — write records with unknown/missing RG to this file.
//!  - `--output-fmt sam|bam` — only `bam` (default) and `sam` are honored.
//!  - `-p N` — width for `%#` padding.

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::{Exact, detect_path};
use htslib_rs::sam;

use crate::diagnostics::{print_error, print_error_errno};

/// Entry point for `samtools split`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut template = String::from("%*_%!.%.");
    let mut unaccounted: Option<PathBuf> = None;
    let mut output_fmt = OutFmt::Bam;
    let mut pad_width: usize = 0;
    let mut input: Option<PathBuf> = None;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-f" => {
                template = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .map(|s| s.to_string())
                    .unwrap_or(template);
            }
            "-u" => {
                unaccounted = iter.next().map(PathBuf::from);
            }
            "-p" => {
                pad_width = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
            }
            "--output-fmt" => {
                let v = iter.next().and_then(|a| a.to_str()).unwrap_or("bam");
                output_fmt = match v.to_lowercase().as_str() {
                    "sam" => OutFmt::Sam,
                    "bam" => OutFmt::Bam,
                    _ => OutFmt::Bam,
                };
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error(
                    "split",
                    format!("option `{}` is not yet supported in samtools-rs split", s),
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

    let Some(input) = input else {
        let _ = print_usage();
        return ExitCode::from(1);
    };

    let format = match detect_path(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error(
                "split",
                format!("failed to detect format of \"{}\": {}", input.display(), e),
            );
            return ExitCode::from(1);
        }
    };
    if format.exact != Exact::Bam {
        print_error(
            "split",
            "only BAM input is currently supported (SAM/CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    match run_split(
        &input,
        &template,
        unaccounted.as_deref(),
        output_fmt,
        pad_width,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("split", "split failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Copy)]
enum OutFmt {
    Sam,
    Bam,
}

fn run_split(
    input: &Path,
    template: &str,
    unaccounted: Option<&Path>,
    fmt: OutFmt,
    pad_width: usize,
) -> io::Result<()> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;

    let rg_ids: Vec<String> = header
        .read_groups()
        .iter()
        .map(|(id, _)| String::from_utf8_lossy(id).into_owned())
        .collect();

    let mut outputs: Vec<Box<dyn BamLike>> = Vec::with_capacity(rg_ids.len());
    let mut id_to_index: HashMap<String, usize> = HashMap::new();
    for (idx, rg) in rg_ids.iter().enumerate() {
        let path = render_template(template, rg, idx, fmt, pad_width);
        outputs.push(open_output(Path::new(&path), fmt, &header)?);
        id_to_index.insert(rg.clone(), idx);
    }

    let mut unaccounted_out: Option<Box<dyn BamLike>> = match unaccounted {
        Some(p) => Some(open_output(p, fmt, &header)?),
        None => None,
    };

    let mut record = bam::Record::default();
    loop {
        let n = reader.read_record(&mut record)?;
        if n == 0 {
            break;
        }
        let rg = record_rg_tag(&record);
        let sink = match rg {
            Some(rg_id) => match id_to_index.get(&rg_id) {
                Some(i) => Some(outputs[*i].as_mut()),
                None => unaccounted_out.as_deref_mut(),
            },
            None => unaccounted_out.as_deref_mut(),
        };
        if let Some(sink) = sink {
            sink.write_record(&header, &record)?;
        }
    }
    Ok(())
}

fn record_rg_tag(record: &bam::Record) -> Option<String> {
    record
        .data()
        .get(b"RG")
        .and_then(|r| r.ok())
        .and_then(|value| match value {
            sam::alignment::record::data::field::Value::String(s) => {
                Some(String::from_utf8_lossy(s.as_ref()).into_owned())
            }
            _ => None,
        })
}

fn render_template(template: &str, rg_id: &str, idx: usize, fmt: OutFmt, pad: usize) -> String {
    let ext = match fmt {
        OutFmt::Sam => "sam",
        OutFmt::Bam => "bam",
    };
    let idx_str = if pad > 0 {
        format!("{:0width$}", idx, width = pad)
    } else {
        idx.to_string()
    };

    let mut out = String::with_capacity(template.len() + 16);
    let mut chars = template.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '%' {
            match chars.next() {
                Some('!') => out.push_str(rg_id),
                Some('#') => out.push_str(&idx_str),
                Some('.') => out.push_str(ext),
                Some('*') => {
                    // No SM lookup yet; substitute the RG ID as a placeholder.
                    out.push_str(rg_id);
                }
                Some(other) => {
                    out.push('%');
                    out.push(other);
                }
                None => out.push('%'),
            }
        } else {
            out.push(c);
        }
    }
    out
}

trait BamLike {
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
impl BamLike for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        self.0.write_record(header, record)
    }
}

struct SamFile(sam::io::Writer<File>);
impl BamLike for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}

fn open_output(path: &Path, fmt: OutFmt, header: &sam::Header) -> io::Result<Box<dyn BamLike>> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)?;
    }
    let file = File::create(path)?;
    match fmt {
        OutFmt::Sam => {
            let mut writer = sam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(SamFile(writer)))
        }
        OutFmt::Bam => {
            let mut writer = bam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(BamFile(writer)))
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools split [options] <in.bam>")?;
    writeln!(w, "Options:")?;
    writeln!(
        w,
        "  -f STR     output filename template (default: %*_%!.%.)"
    )?;
    writeln!(w, "  -u FILE    write records with unknown @RG to FILE")?;
    writeln!(w, "  -p N       pad the %# index to N digits")?;
    writeln!(w, "  --output-fmt sam|bam")?;
    Ok(())
}
