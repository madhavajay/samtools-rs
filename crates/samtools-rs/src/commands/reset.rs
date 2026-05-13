//! `samtools reset` — revert aligner changes in reads.
//!
//! Mirrors `main_reset` in `reset.c`. Initial Rust port operates on BAM
//! input and writes BAM output via `RecordBuf` mutation:
//!  - clears `reference_sequence_id`, `alignment_start`, `cigar`,
//!    `mapping_quality`, `mate_reference_sequence_id`,
//!    `mate_alignment_start`, `template_length`
//!  - clears flag bits that depend on alignment (FUNMAP set to 1,
//!    FSECONDARY/FSUPPLEMENTARY/FPROPER_PAIR/FMUNMAP/FREVERSE/FMREVERSE
//!    cleared)
//!  - drops common aligner-added aux tags (NM, MD, AS, XS, SA, MC, MQ,
//!    NH, HI) by default
//!
//! Reverse-strand sequence/quality re-reversal is **not yet implemented**.
//! `--no-RG`, `--reject-PG`, `--no-PG`, and `--dupflag` are accepted but
//! ignored for now. `-x`/`--keep-tag` aux-tag filtering is honored.

use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{
    self,
    alignment::{RecordBuf, record::Flags},
};

use crate::aux_list::{AuxTag, parse_aux_list};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

const DEFAULT_DROP_TAGS: &[&[u8; 2]] = &[
    b"NM", b"MD", b"AS", b"XS", b"SA", b"MC", b"MQ", b"NH", b"HI", b"ms",
];

/// Entry point for `samtools reset`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut output: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let mut output_fmt = OutFmt::Bam;
    let mut extra_drop: Vec<AuxTag> = Vec::new();
    let mut keep_only: Option<HashSet<AuxTag>> = None;

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
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "-x" | "--remove-tag" | "--remove-tags" => {
                if let Some(v) = iter.next().and_then(|a| a.to_str()) {
                    if let Some(rest) = v.strip_prefix('^') {
                        match parse_aux_list(rest) {
                            Ok(tags) => keep_only = Some(tags),
                            Err(e) => {
                                print_error("reset", format!("invalid -x value \"{rest}\": {e}"));
                                return ExitCode::from(1);
                            }
                        }
                    } else {
                        match parse_aux_list(v) {
                            Ok(tags) => extra_drop.extend(tags),
                            Err(e) => {
                                print_error("reset", format!("invalid -x value \"{v}\": {e}"));
                                return ExitCode::from(1);
                            }
                        }
                    }
                }
            }
            "--keep-tag" | "--keep-tags" => {
                if let Some(v) = iter.next().and_then(|a| a.to_str()) {
                    match parse_aux_list(v) {
                        Ok(tags) => keep_only = Some(tags),
                        Err(e) => {
                            print_error("reset", format!("invalid --keep-tag value \"{v}\": {e}"));
                            return ExitCode::from(1);
                        }
                    }
                }
            }
            "--no-RG" | "--reject-PG" | "--no-PG" | "--dupflag" | "-T" => {
                if matches!(s, "--reject-PG" | "-T") {
                    let _ = iter.next();
                }
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("reset", format!("unknown option {}", s));
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
            print_error("reset", e.to_string());
            return ExitCode::from(1);
        }
    };
    if format.exact != Exact::Bam {
        print_error(
            "reset",
            "only BAM input is currently supported (SAM/CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    match run_reset(
        &input,
        output.as_deref(),
        output_fmt,
        &extra_drop,
        keep_only.as_ref(),
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("reset", "reset failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Copy)]
enum OutFmt {
    Sam,
    Bam,
}

fn run_reset(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    extra_drop: &[[u8; 2]],
    keep_only: Option<&HashSet<[u8; 2]>>,
) -> io::Result<()> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;
    let mut sink = open_output(output, fmt, &header)?;

    let mut record = RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut record)?;
        if n == 0 {
            break;
        }
        reset_record(&mut record, extra_drop, keep_only);
        sink.write_record(&header, &record)?;
    }
    Ok(())
}

fn reset_record(
    record: &mut RecordBuf,
    extra_drop: &[[u8; 2]],
    keep_only: Option<&HashSet<[u8; 2]>>,
) {
    // Reset alignment fields.
    *record.reference_sequence_id_mut() = None;
    *record.alignment_start_mut() = None;
    *record.cigar_mut() = sam::alignment::record_buf::Cigar::default();
    *record.mapping_quality_mut() = None;
    *record.mate_reference_sequence_id_mut() = None;
    *record.mate_alignment_start_mut() = None;
    *record.template_length_mut() = 0;

    // Reset flag bits.
    let mut flags = record.flags();
    flags.remove(Flags::PROPERLY_SEGMENTED);
    flags.remove(Flags::SECONDARY);
    flags.remove(Flags::SUPPLEMENTARY);
    flags.remove(Flags::DUPLICATE);
    flags.remove(Flags::MATE_UNMAPPED);
    flags.remove(Flags::REVERSE_COMPLEMENTED);
    flags.remove(Flags::MATE_REVERSE_COMPLEMENTED);
    flags.insert(Flags::UNMAPPED);
    *record.flags_mut() = flags;

    // Drop aligner-added aux tags.
    let data = record.data_mut();
    let mut to_drop: HashSet<[u8; 2]> = HashSet::new();
    for tag in DEFAULT_DROP_TAGS {
        to_drop.insert(**tag);
    }
    for tag in extra_drop {
        to_drop.insert(*tag);
    }
    let mut keys: Vec<sam::alignment::record::data::field::Tag> =
        data.iter().map(|(t, _)| t).collect();
    for k in keys.drain(..) {
        let bytes: [u8; 2] = k.into();
        let should_drop = match keep_only {
            Some(keep) => !keep.contains(&bytes),
            None => to_drop.contains(&bytes),
        };
        if should_drop {
            data.remove(&k);
        }
    }
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
    writeln!(w, "Usage: samtools reset [options] <in.bam>")?;
    writeln!(w, "  -o FILE                 output FILE")?;
    writeln!(w, "  -O sam|bam              output format")?;
    writeln!(
        w,
        "  -x/--remove-tag TAG     drop the listed aux tags (comma-separated, ^ for keep)"
    )?;
    writeln!(w, "  --keep-tag TAG          only keep the listed aux tags")?;
    Ok(())
}
