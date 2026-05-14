//! `samtools fixmate` — fix mate-related flags and positions on paired records.
//!
//! Mirrors `bam_mate.c` in upstream samtools. Initial Rust port handles
//! **name-sorted BAM/SAM input**: adjacent records with the same `qname` are
//! paired up and their `FMUNMAP`/`FMREVERSE` flags + `mate_reference_sequence_id`
//! + `mate_alignment_start` are made consistent.
//!
//! **Not yet supported:** MC/MQ aux-tag updates, `-r` (rescore secondary
//! alignments), `-c` (calculate CT), `-m` (add ms score), CRAM input/output.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{
    self,
    alignment::{RecordBuf, record::Flags},
};

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sanitize::{SanitizeFlags, parse_sanitize_options};

/// Entry point for `samtools fixmate`.
pub fn main(args: &[OsString]) -> ExitCode {
    let opts = match parse_args(args) {
        Ok(opts) => opts,
        Err(ParseError::Usage) => {
            let _ = print_usage();
            return ExitCode::SUCCESS;
        }
        Err(ParseError::Err(e)) => {
            print_error("fixmate", e);
            return ExitCode::from(1);
        }
    };

    let _ = opts.sanitize_flags; // Parsed for parity; record mutation is still TODO.

    let Some(input) = opts.input else {
        let _ = print_usage();
        return ExitCode::from(1);
    };

    let format = match sam_io::sam_open_format(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error("fixmate", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam) {
        print_error(
            "fixmate",
            "only SAM and BAM input are currently supported (CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    match run_fixmate(&input, opts.output.as_deref(), opts.output_fmt) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("fixmate", "fixmate failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct Opts {
    output: Option<PathBuf>,
    input: Option<PathBuf>,
    output_fmt: OutFmt,
    sanitize_flags: Option<SanitizeFlags>,
}

impl Default for Opts {
    fn default() -> Self {
        Self {
            output: None,
            input: None,
            output_fmt: OutFmt::Bam,
            sanitize_flags: None,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum ParseError {
    Usage,
    Err(String),
}

fn parse_args(args: &[OsString]) -> Result<Opts, ParseError> {
    let mut opts = Opts::default();
    let mut iter = args.iter().skip(1);

    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-O" | "--output-fmt" => {
                let v = next_value(&mut iter, s)?;
                opts.output_fmt = match v.to_lowercase().as_str() {
                    "sam" => OutFmt::Sam,
                    "bam" => OutFmt::Bam,
                    _ => OutFmt::Bam,
                };
            }
            "-@" | "--threads" | "-l" => {
                let _ = next_value(&mut iter, s)?;
            }
            "-z" | "--sanitize" => {
                let v = next_value(&mut iter, s)?;
                opts.sanitize_flags = Some(parse_sanitize_options(&v).map_err(ParseError::Err)?);
            }
            "-r" | "-c" | "-m" | "-p" | "--no-PG" => {
                // Accepted but not yet implemented.
            }
            "--help" => return Err(ParseError::Usage),
            _ if s.starts_with('-') && s != "-" => {
                return Err(ParseError::Err(format!("unknown option {}", s)));
            }
            _ => {
                if opts.input.is_none() {
                    opts.input = Some(PathBuf::from(arg));
                } else if opts.output.is_none() {
                    opts.output = Some(PathBuf::from(arg));
                }
            }
        }
    }

    Ok(opts)
}

fn next_value<'a, I>(iter: &mut I, option: &str) -> Result<String, ParseError>
where
    I: Iterator<Item = &'a OsString>,
{
    iter.next()
        .and_then(|a| a.to_str().map(str::to_owned))
        .ok_or_else(|| ParseError::Err(format!("missing value for {option}")))
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum OutFmt {
    Sam,
    Bam,
}

fn run_fixmate(input: &Path, output: Option<&Path>, fmt: OutFmt) -> io::Result<()> {
    let format = sam_io::sam_open_format(input)?;
    match format.exact {
        Exact::Sam => run_fixmate_sam(input, output, fmt),
        Exact::Bam => run_fixmate_bam(input, output, fmt),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "only SAM and BAM input are currently supported (CRAM TODO)",
        )),
    }
}

fn run_fixmate_bam(input: &Path, output: Option<&Path>, fmt: OutFmt) -> io::Result<()> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;
    let mut sink = open_output(output, fmt, &header)?;
    let mut pending: Option<RecordBuf> = None;
    let mut next = RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut next)?;
        if n == 0 {
            break;
        }
        write_fixed_record(&header, sink.as_mut(), &mut pending, next.clone())?;
    }
    if let Some(rec) = pending {
        sink.write_record(&header, &rec)?;
    }
    Ok(())
}

fn run_fixmate_sam(input: &Path, output: Option<&Path>, fmt: OutFmt) -> io::Result<()> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(input)?));
    let header = reader.read_header()?;
    let mut sink = open_output(output, fmt, &header)?;
    let mut pending: Option<RecordBuf> = None;
    let mut next = RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut next)?;
        if n == 0 {
            break;
        }
        write_fixed_record(&header, sink.as_mut(), &mut pending, next.clone())?;
    }
    if let Some(rec) = pending {
        sink.write_record(&header, &rec)?;
    }
    Ok(())
}

fn write_fixed_record(
    header: &sam::Header,
    sink: &mut dyn Sink,
    pending: &mut Option<RecordBuf>,
    next: RecordBuf,
) -> io::Result<()> {
    match pending.take() {
        None => *pending = Some(next),
        Some(prev) => {
            let prev_name = prev.name().map(|n| n.to_vec());
            let next_name = next.name().map(|n| n.to_vec());
            if prev_name == next_name && next_name.is_some() {
                let (a, b) = pair_fixmate(prev, next);
                sink.write_record(header, &a)?;
                sink.write_record(header, &b)?;
                *pending = None;
            } else {
                sink.write_record(header, &prev)?;
                *pending = Some(next);
            }
        }
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
    writeln!(w, "  -z, --sanitize FLAG[,FLAG]")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(rest: &[&str]) -> Vec<OsString> {
        std::iter::once(OsString::from("fixmate"))
            .chain(rest.iter().map(OsString::from))
            .collect()
    }

    #[test]
    fn parses_sanitize_option_without_treating_value_as_input() {
        let opts = parse_args(&argv(&["-z", "on", "in.bam", "out.bam"])).unwrap();

        assert_eq!(opts.input.as_deref(), Some(Path::new("in.bam")));
        assert_eq!(opts.output.as_deref(), Some(Path::new("out.bam")));
        assert!(opts.sanitize_flags.unwrap().contains(SanitizeFlags::CIGAR));
    }

    #[test]
    fn rejects_missing_sanitize_value() {
        assert_eq!(
            parse_args(&argv(&["--sanitize"])).unwrap_err(),
            ParseError::Err(String::from("missing value for --sanitize"))
        );
    }

    #[test]
    fn rejects_invalid_sanitize_value() {
        assert_eq!(
            parse_args(&argv(&["--sanitize", "nope"])).unwrap_err(),
            ParseError::Err(String::from("unrecognised sanitize keyword \"nope\""))
        );
    }
}
