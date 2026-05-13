//! `samtools fixmate` — fix mate-related flags and positions on paired records.
//!
//! Mirrors `bam_mate.c` in upstream samtools. Initial Rust port handles
//! **name-sorted BAM input**: adjacent records with the same `qname` are
//! paired up and their `FMUNMAP`/`FMREVERSE` flags + `mate_reference_sequence_id`
//! + `mate_alignment_start` are made consistent.
//!
//! **Not yet supported:** MC/MQ aux-tag updates, `-r` (rescore secondary
//! alignments), `-c` (calculate CT), `-m` (add ms score), CRAM input/output.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::{Exact, detect_path};
use htslib_rs::sam::{
    self,
    alignment::{RecordBuf, record::Flags},
};

use crate::diagnostics::{print_error, print_error_errno};

/// Entry point for `samtools fixmate`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut output: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let mut output_fmt = OutFmt::Bam;
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-O" | "--output-fmt" => {
                let v = iter.next().and_then(|a| a.to_str()).unwrap_or("bam");
                output_fmt = match v.to_lowercase().as_str() {
                    "sam" => OutFmt::Sam,
                    "bam" => OutFmt::Bam,
                    _ => OutFmt::Bam,
                };
            }
            "-@" | "--threads" | "-l" => {
                let _ = iter.next();
            }
            "-r" | "-c" | "-m" | "-p" | "-z" | "--no-PG" => {
                // Accepted but not yet implemented.
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("fixmate", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                } else if output.is_none() {
                    output = Some(PathBuf::from(arg));
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
                "fixmate",
                format!("failed to detect format of \"{}\": {}", input.display(), e),
            );
            return ExitCode::from(1);
        }
    };
    if format.exact != Exact::Bam {
        print_error(
            "fixmate",
            "only BAM input is currently supported (SAM/CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    match run_fixmate(&input, output.as_deref(), output_fmt) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("fixmate", "fixmate failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Copy)]
enum OutFmt {
    Sam,
    Bam,
}

fn run_fixmate(input: &Path, output: Option<&Path>, fmt: OutFmt) -> io::Result<()> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;

    let mut pending: Option<RecordBuf> = None;
    let mut sink = open_output(output, fmt, &header)?;

    let mut next = RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut next)?;
        if n == 0 {
            break;
        }
        match pending.take() {
            None => pending = Some(next.clone()),
            Some(prev) => {
                let prev_name = prev.name().map(|n| n.to_vec());
                let next_name = next.name().map(|n| n.to_vec());
                if prev_name == next_name && next_name.is_some() {
                    let (a, b) = pair_fixmate(prev, next.clone());
                    sink.write_record(&header, &a)?;
                    sink.write_record(&header, &b)?;
                    pending = None;
                } else {
                    // Singleton: emit prev as-is.
                    sink.write_record(&header, &prev)?;
                    pending = Some(next.clone());
                }
            }
        }
    }
    if let Some(rec) = pending {
        sink.write_record(&header, &rec)?;
    }
    Ok(())
}

fn pair_fixmate(mut a: RecordBuf, mut b: RecordBuf) -> (RecordBuf, RecordBuf) {
    let a_tid = a.reference_sequence_id();
    let b_tid = b.reference_sequence_id();
    let a_pos = a.alignment_start();
    let b_pos = b.alignment_start();
    apply_mate_flags(&mut a, &b);
    apply_mate_flags(&mut b, &a);
    *a.mate_reference_sequence_id_mut() = b_tid;
    *b.mate_reference_sequence_id_mut() = a_tid;
    *a.mate_alignment_start_mut() = b_pos;
    *b.mate_alignment_start_mut() = a_pos;
    (a, b)
}

fn apply_mate_flags(target: &mut RecordBuf, mate: &RecordBuf) {
    let mut flags = target.flags();
    flags.insert(Flags::SEGMENTED);
    if mate.flags().is_unmapped() {
        flags.insert(Flags::MATE_UNMAPPED);
    } else {
        flags.remove(Flags::MATE_UNMAPPED);
    }
    if mate.flags().is_reverse_complemented() {
        flags.insert(Flags::MATE_REVERSE_COMPLEMENTED);
    } else {
        flags.remove(Flags::MATE_REVERSE_COMPLEMENTED);
    }
    *target.flags_mut() = flags;
}

trait Sink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(sam::io::Writer<File>);
struct SamStdout(sam::io::Writer<io::Stdout>);

impl Sink for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl Sink for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl Sink for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl Sink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}

fn open_output(out: Option<&Path>, fmt: OutFmt, header: &sam::Header) -> io::Result<Box<dyn Sink>> {
    match (out, fmt) {
        (Some(p), OutFmt::Sam) => {
            let mut w = sam::io::Writer::new(File::create(p)?);
            w.write_header(header)?;
            Ok(Box::new(SamFile(w)))
        }
        (Some(p), OutFmt::Bam) => {
            let mut w = bam::io::Writer::new(File::create(p)?);
            w.write_header(header)?;
            Ok(Box::new(BamFile(w)))
        }
        (None, OutFmt::Sam) => {
            let mut w = sam::io::Writer::new(io::stdout());
            w.write_header(header)?;
            Ok(Box::new(SamStdout(w)))
        }
        (None, OutFmt::Bam) => {
            let mut w = bam::io::Writer::new(io::stdout());
            w.write_header(header)?;
            Ok(Box::new(BamStdout(w)))
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools fixmate [options] <in.bam> [<out.bam>]")?;
    writeln!(w, "  -O sam|bam   output format (default: bam)")?;
    Ok(())
}
