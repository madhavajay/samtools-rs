//! `samtools cat` — concatenate alignment files of the same format.
//!
//! Mirrors `main_cat` in `bam_cat.c`. The upstream implementation
//! concatenates BAM files at the BGZF block level for speed and supports
//! `-h <hdr>` (replace header), `-o <out>` (output), and `-p N/M` for CRAM.
//!
//! This Rust port implements a record-level concatenation (decompress +
//! re-encode) for BAM. CRAM concatenation and `-p` are not yet supported.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam;

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

/// Entry point for `samtools cat`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut header_file: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut no_pg = false;
    let mut inputs: Vec<PathBuf> = Vec::new();
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-h" => {
                let v = iter.next().map(PathBuf::from);
                header_file = v;
            }
            "-o" => {
                let v = iter.next().map(PathBuf::from);
                output = v;
            }
            "--no-PG" => {
                no_pg = true;
            }
            "-r" | "-p" | "-q" | "-f" | "-b" => {
                // Reserved upstream flags not yet supported.
                print_error(
                    "cat",
                    format!("option `{}` is not yet supported in samtools-rs cat", s),
                );
                return ExitCode::from(1);
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("cat", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => {
                inputs.push(PathBuf::from(arg));
            }
        }
    }

    let _ = no_pg; // currently we never add a @PG line ourselves.

    if inputs.is_empty() {
        let _ = writeln!(
            io::stderr(),
            "Usage: samtools cat [-h hdr] [-o out] in1.bam ..."
        );
        return ExitCode::from(1);
    }

    // Determine input format from the first file.
    let format = match sam_io::sam_open_format(&inputs[0]) {
        Ok(f) => f,
        Err(e) => {
            print_error("cat", e.to_string());
            return ExitCode::from(1);
        }
    };

    if format.exact != Exact::Bam {
        print_error(
            "cat",
            "only BAM input is currently supported (CRAM/SAM TODO)",
        );
        return ExitCode::from(1);
    }

    match run_bam_cat(&inputs, header_file.as_deref(), output.as_deref()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("cat", "concatenation failed", &e);
            ExitCode::from(1)
        }
    }
}

fn run_bam_cat(
    inputs: &[PathBuf],
    header_file: Option<&Path>,
    output: Option<&Path>,
) -> io::Result<()> {
    // Pick the header. If -h <hdr.sam>, parse that; otherwise use the first
    // input file's header.
    let header: sam::Header = match header_file {
        Some(p) => {
            let mut reader = sam::io::Reader::new(BufReader::new(File::open(p)?));
            reader.read_header()?
        }
        None => {
            let mut reader = bam::io::Reader::new(File::open(&inputs[0])?);
            reader.read_header()?
        }
    };

    // Open output writer.
    let mut writer: Box<dyn BamSink> = match output {
        Some(p) => Box::new(FileSink::new(p)?),
        None => Box::new(StdoutSink::new()?),
    };
    writer.write_header(&header)?;

    for input in inputs {
        let mut reader = bam::io::Reader::new(File::open(input)?);
        let input_header = reader.read_header()?;
        let mut record = bam::Record::default();
        loop {
            let n = reader.read_record(&mut record)?;
            if n == 0 {
                break;
            }
            writer.write_record(&input_header, &record)?;
        }
    }
    writer.finish()
}

trait BamSink {
    fn write_header(&mut self, header: &sam::Header) -> io::Result<()>;
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()>;
    fn finish(self: Box<Self>) -> io::Result<()>;
}

struct FileSink {
    writer: bam::io::Writer<bgzf::io::Writer<File>>,
}

impl FileSink {
    fn new(path: &Path) -> io::Result<Self> {
        Ok(Self {
            writer: bam::io::Writer::new(File::create(path)?),
        })
    }
}

impl BamSink for FileSink {
    fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.writer.write_header(header)
    }
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        self.writer.write_record(header, record)
    }
    fn finish(self: Box<Self>) -> io::Result<()> {
        Ok(())
    }
}

struct StdoutSink {
    writer: bam::io::Writer<bgzf::io::Writer<io::Stdout>>,
}

impl StdoutSink {
    fn new() -> io::Result<Self> {
        Ok(Self {
            writer: bam::io::Writer::new(io::stdout()),
        })
    }
}

impl BamSink for StdoutSink {
    fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.writer.write_header(header)
    }
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        self.writer.write_record(header, record)
    }
    fn finish(self: Box<Self>) -> io::Result<()> {
        Ok(())
    }
}
