//! `samtools collate` / `samtools bamshuf` — group reads with the same name
//! together, with a (pseudo-)randomised order otherwise.
//!
//! Mirrors `main_bamshuf` in `bamshuf.c`. The upstream implementation uses
//! name-hash bucketing with on-disk temp files for memory bounding. This
//! initial Rust port performs an in-memory name sort for BAM/SAM inputs, which
//! gives the same per-name grouping result but does not scale to inputs larger
//! than memory.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf};

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

/// Entry point for `samtools collate` (alias `bamshuf`).
pub fn main(args: &[OsString]) -> ExitCode {
    let mut output_prefix: Option<String> = None;
    let mut to_stdout = false;
    let mut output_fmt = OutFmt::Bam;
    let mut input: Option<PathBuf> = None;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-o" => {
                output_prefix = iter.next().and_then(|a| a.to_str()).map(|s| s.to_string());
            }
            "-O" => {
                to_stdout = true;
            }
            "--output-fmt" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("collate", "missing value for --output-fmt");
                    return ExitCode::from(1);
                };
                output_fmt = match parse_output_format(v) {
                    Ok(fmt) => fmt,
                    Err(e) => {
                        print_error("collate", e);
                        return ExitCode::from(1);
                    }
                };
            }
            "--no-PG" => {}
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error(
                    "collate",
                    format!("option `{}` is not yet supported in samtools-rs collate", s),
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

    let format = match sam_io::sam_open_format(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error("collate", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam) {
        print_error(
            "collate",
            "only SAM and BAM input are currently supported (CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    let output = if to_stdout {
        OutputTarget::Stdout
    } else {
        let prefix = output_prefix.unwrap_or_else(|| "collated".to_string());
        let ext = match output_fmt {
            OutFmt::Sam => "sam",
            OutFmt::Bam => "bam",
        };
        let path = if prefix.ends_with(&format!(".{}", ext)) {
            PathBuf::from(prefix)
        } else {
            PathBuf::from(format!("{}.{}", prefix, ext))
        };
        OutputTarget::File(path)
    };

    match run_collate(&input, output, output_fmt) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("collate", "collate failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Copy)]
enum OutFmt {
    Sam,
    Bam,
}

enum OutputTarget {
    Stdout,
    File(PathBuf),
}

fn parse_output_format(raw: &str) -> Result<OutFmt, String> {
    match raw.to_ascii_lowercase().as_str() {
        "sam" => Ok(OutFmt::Sam),
        "bam" => Ok(OutFmt::Bam),
        _ => Err(format!("unsupported output format \"{}\"", raw)),
    }
}

fn run_collate(input: &Path, output: OutputTarget, fmt: OutFmt) -> io::Result<()> {
    let format = sam_io::sam_open_format(input)?;
    let (header, mut records) = match format.exact {
        Exact::Sam => read_sam_records(input)?,
        Exact::Bam => read_bam_records(input)?,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "only SAM and BAM input are currently supported (CRAM TODO)",
            ));
        }
    };

    records.sort_by(|a, b| {
        let an = a.name().map(|s| s.to_vec()).unwrap_or_default();
        let bn = b.name().map(|s| s.to_vec()).unwrap_or_default();
        an.cmp(&bn)
    });

    let mut writer = open_output(&output, fmt, &header)?;
    for rec in &records {
        writer.write_record(&header, rec)?;
    }
    Ok(())
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

trait CollateSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(sam::io::Writer<File>);
struct SamStdout(sam::io::Writer<io::Stdout>);

impl CollateSink for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl CollateSink for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl CollateSink for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl CollateSink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}

fn open_output(
    out: &OutputTarget,
    fmt: OutFmt,
    header: &sam::Header,
) -> io::Result<Box<dyn CollateSink>> {
    match (out, fmt) {
        (OutputTarget::Stdout, OutFmt::Sam) => {
            let mut writer = sam::io::Writer::new(io::stdout());
            writer.write_header(header)?;
            Ok(Box::new(SamStdout(writer)))
        }
        (OutputTarget::Stdout, OutFmt::Bam) => {
            let mut writer = bam::io::Writer::new(io::stdout());
            writer.write_header(header)?;
            Ok(Box::new(BamStdout(writer)))
        }
        (OutputTarget::File(p), OutFmt::Sam) => {
            let file = File::create(p)?;
            let mut writer = sam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(SamFile(writer)))
        }
        (OutputTarget::File(p), OutFmt::Bam) => {
            let file = File::create(p)?;
            let mut writer = bam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(BamFile(writer)))
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(
        w,
        "Usage: samtools collate [options] <in.bam|in.sam> [<out.prefix>]"
    )?;
    writeln!(w, "  -o PREFIX   output prefix or path")?;
    writeln!(w, "  -O          write to stdout")?;
    writeln!(w, "  --output-fmt sam|bam")?;
    Ok(())
}
