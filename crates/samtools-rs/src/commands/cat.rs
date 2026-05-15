//! `samtools cat` — concatenate alignment files of the same format.
//!
//! Mirrors `main_cat` in `bam_cat.c`. The upstream implementation
//! concatenates BAM files at the BGZF block level for speed and supports
//! `-b <fofn>` (input file list), `-h <hdr>` (replace header), `-o <out>`
//! (output), and `-p N/M` for CRAM.
//!
//! This Rust port implements record-level concatenation (decompress +
//! re-encode) for SAM and BAM. CRAM concatenation and `-p` are not yet
//! supported.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf};

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

/// Entry point for `samtools cat`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut header_file: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut no_pg = false;
    let mut region: Option<String> = None;
    let mut input_lists: Vec<PathBuf> = Vec::new();
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
            "-b" => {
                let Some(v) = iter.next() else {
                    print_error("cat", "missing value for -b");
                    return ExitCode::from(1);
                };
                input_lists.push(PathBuf::from(v));
            }
            "--no-PG" => {
                no_pg = true;
            }
            "-r" => {
                region = iter.next().and_then(|a| a.to_str().map(str::to_owned));
            }
            "-p" | "-q" | "-f" => {
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

    if !input_lists.is_empty() {
        let mut expanded = Vec::new();
        for list in &input_lists {
            match read_input_list(list) {
                Ok(mut list_inputs) => expanded.append(&mut list_inputs),
                Err(e) => {
                    print_error("cat", e.to_string());
                    return ExitCode::from(1);
                }
            }
        }
        expanded.extend(inputs);
        inputs = expanded;
    }

    if inputs.is_empty() {
        let _ = writeln!(
            io::stderr(),
            "Usage: samtools cat [-b list] [-h hdr] [-o out] in1.bam ..."
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

    let result = match format.exact {
        Exact::Sam => run_sam_cat(
            &inputs,
            header_file.as_deref(),
            output.as_deref(),
            !no_pg,
            args,
            region.as_deref(),
        ),
        Exact::Bam => run_bam_cat(
            &inputs,
            header_file.as_deref(),
            output.as_deref(),
            !no_pg,
            args,
            region.as_deref(),
        ),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "only SAM and BAM input are currently supported (CRAM TODO)",
        )),
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("cat", "concatenation failed", &e);
            ExitCode::from(1)
        }
    }
}

fn read_input_list(path: &Path) -> io::Result<Vec<PathBuf>> {
    let file = File::open(path)?;
    let mut inputs = Vec::new();

    for line in BufReader::new(file).lines() {
        let line = line?;
        let line = line.trim();
        if !line.is_empty() {
            inputs.push(PathBuf::from(line));
        }
    }

    Ok(inputs)
}

fn run_bam_cat(
    inputs: &[PathBuf],
    header_file: Option<&Path>,
    output: Option<&Path>,
    add_pg: bool,
    argv: &[OsString],
    region: Option<&str>,
) -> io::Result<()> {
    let parsed_region = match region {
        Some(r) => Some(r.parse::<htslib_rs::core::Region>().map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("invalid -r region \"{r}\": {e}"),
            )
        })?),
        None => None,
    };

    // Pick the header. If -h <hdr.sam>, parse that; otherwise use the first
    // input file's header.
    let mut header: sam::Header = match header_file {
        Some(p) => {
            let mut reader = sam::io::Reader::new(BufReader::new(File::open(p)?));
            reader.read_header()?
        }
        None => {
            let mut reader = bam::io::Reader::new(File::open(&inputs[0])?);
            reader.read_header()?
        }
    };
    if add_pg {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

    // Open output writer.
    let mut writer: Box<dyn BamSink> = match output {
        Some(p) => Box::new(FileSink::new(p)?),
        None => Box::new(StdoutSink::new()?),
    };
    writer.write_header(&header)?;

    for input in inputs {
        let mut reader = bam::io::Reader::new(File::open(input)?);
        let input_header = reader.read_header()?;
        match parsed_region.as_ref() {
            Some(region) => {
                for record in
                    htslib_rs::alignment_compat::query_bam_records_from_path(input, region)?
                {
                    writer.write_record(&input_header, &record)?;
                }
            }
            None => {
                let mut record = bam::Record::default();
                loop {
                    let n = reader.read_record(&mut record)?;
                    if n == 0 {
                        break;
                    }
                    writer.write_record(&input_header, &record)?;
                }
            }
        }
    }
    writer.finish()
}

fn run_sam_cat(
    inputs: &[PathBuf],
    header_file: Option<&Path>,
    output: Option<&Path>,
    add_pg: bool,
    argv: &[OsString],
    region: Option<&str>,
) -> io::Result<()> {
    if let Some(region) = region {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("-r region \"{region}\" is currently supported for indexed BAM input only"),
        ));
    }

    for input in inputs {
        let format = sam_io::sam_open_format(input)?;
        if format.exact != Exact::Sam {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "all cat inputs must use the same format",
            ));
        }
    }

    let mut header: sam::Header = match header_file {
        Some(p) => {
            let mut reader = sam::io::Reader::new(BufReader::new(File::open(p)?));
            reader.read_header()?
        }
        None => {
            let mut reader = sam::io::Reader::new(BufReader::new(File::open(&inputs[0])?));
            reader.read_header()?
        }
    };
    if add_pg {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

    let mut writer: Box<dyn SamSink> = match output {
        Some(p) => Box::new(SamFileSink::new(p, &header)?),
        None => Box::new(SamStdoutSink::new(&header)?),
    };

    for input in inputs {
        let mut reader = sam::io::Reader::new(BufReader::new(File::open(input)?));
        let input_header = reader.read_header()?;
        loop {
            let mut record = RecordBuf::default();
            if reader.read_record_buf(&input_header, &mut record)? == 0 {
                break;
            }
            writer.write_record(&header, &record)?;
        }
    }
    writer.finish()
}

trait SamSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
    fn finish(self: Box<Self>) -> io::Result<()>;
}

struct SamFileSink {
    writer: File,
}

impl SamFileSink {
    fn new(path: &Path, header: &sam::Header) -> io::Result<Self> {
        let mut writer = File::create(path)?;
        crate::sam_render::write_header(&mut writer, header)?;
        Ok(Self { writer })
    }
}

impl SamSink for SamFileSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        // Shared renderer: htslib `%g` float aux spelling.
        crate::sam_render::write_record(&mut self.writer, header, record)
    }

    fn finish(self: Box<Self>) -> io::Result<()> {
        Ok(())
    }
}

struct SamStdoutSink {
    writer: io::Stdout,
}

impl SamStdoutSink {
    fn new(header: &sam::Header) -> io::Result<Self> {
        let mut writer = io::stdout();
        crate::sam_render::write_header(&mut writer, header)?;
        Ok(Self { writer })
    }
}

impl SamSink for SamStdoutSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        crate::sam_render::write_record(&mut self.writer, header, record)
    }

    fn finish(self: Box<Self>) -> io::Result<()> {
        Ok(())
    }
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
