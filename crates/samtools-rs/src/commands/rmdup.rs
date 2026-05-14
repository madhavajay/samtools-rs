//! `samtools rmdup` — remove PCR duplicates (deprecated; `markdup` is preferred).
//!
//! Mirrors `bam_rmdup` / `bam_rmdupse` in upstream samtools. Upstream's
//! implementation is paired-aware and works on coordinate-sorted BAMs.
//!
//! This Rust port implements a single-end variant: for each
//! `(reference_sequence_id, alignment_start, reverse-flag)` group, keep
//! the record with the highest mapping quality and drop the rest.
//! Requires coordinate-sorted BAM input.

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf, alignment::io::Write as _};

use crate::bam_flag::{BAM_FREVERSE, BAM_FUNMAP};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

/// Entry point for `samtools rmdup`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut output: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let iter = args.iter().skip(1);
    for arg in iter {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-S" | "-s" => {
                // upstream: -s treats paired-end as single-end; here single-end
                // is the only mode.
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("rmdup", format!("unknown option {}", s));
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

    let format = match sam_io::sam_open_format(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error("rmdup", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam) {
        print_error(
            "rmdup",
            "only SAM and BAM input are currently supported (CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    let result = match format.exact {
        Exact::Sam => run_sam_rmdup(&input, output.as_deref()),
        Exact::Bam => run_bam_rmdup(&input, output.as_deref()),
        _ => unreachable!("format checked above"),
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("rmdup", "rmdup failed", &e);
            ExitCode::from(1)
        }
    }
}

fn run_bam_rmdup(input: &Path, output: Option<&Path>) -> io::Result<()> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;

    let mut records: Vec<bam::Record> = Vec::new();
    let mut record = bam::Record::default();
    loop {
        let n = reader.read_record(&mut record)?;
        if n == 0 {
            break;
        }
        records.push(record.clone());
    }

    // Group by (tid, pos, reverse-flag) and keep the record with the
    // highest MAPQ in each group. Unmapped records pass through.
    let mut best_per_group: HashMap<(i32, i64, bool), usize> = HashMap::new();
    let mut keep = vec![false; records.len()];
    for (i, rec) in records.iter().enumerate() {
        let flag = u16::from(rec.flags()) as u32;
        if flag & BAM_FUNMAP != 0 {
            keep[i] = true;
            continue;
        }
        let tid = rec
            .reference_sequence_id()
            .and_then(|res| res.ok())
            .map(|t| t as i32)
            .unwrap_or(-1);
        let pos = rec
            .alignment_start()
            .and_then(|res| res.ok())
            .map(|p| usize::from(p) as i64)
            .unwrap_or(0);
        let rev = flag & BAM_FREVERSE != 0;
        let mapq = rec.mapping_quality().map(u8::from).unwrap_or(0);
        let key = (tid, pos, rev);
        match best_per_group.get(&key) {
            Some(&idx) => {
                let prev_mapq = records[idx].mapping_quality().map(u8::from).unwrap_or(0);
                if mapq > prev_mapq {
                    keep[idx] = false;
                    keep[i] = true;
                    best_per_group.insert(key, i);
                }
            }
            None => {
                keep[i] = true;
                best_per_group.insert(key, i);
            }
        }
    }

    let mut writer = open_bam_output(output, &header)?;
    for (i, rec) in records.iter().enumerate() {
        if keep[i] {
            writer.write_record(&header, rec)?;
        }
    }
    Ok(())
}

fn run_sam_rmdup(input: &Path, output: Option<&Path>) -> io::Result<()> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(input)?));
    let header = reader.read_header()?;

    let mut records: Vec<RecordBuf> = Vec::new();
    loop {
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        records.push(record);
    }

    let keep = duplicate_keep_mask_for_sam_records(&records);
    let mut writer = open_sam_output(output, &header)?;
    for (i, rec) in records.iter().enumerate() {
        if keep[i] {
            writer.write_record(&header, rec)?;
        }
    }
    Ok(())
}

fn duplicate_keep_mask_for_sam_records(records: &[RecordBuf]) -> Vec<bool> {
    let mut best_per_group: HashMap<(i32, i64, bool), usize> = HashMap::new();
    let mut keep = vec![false; records.len()];
    for (i, rec) in records.iter().enumerate() {
        let flag = rec.flags().bits() as u32;
        if flag & BAM_FUNMAP != 0 {
            keep[i] = true;
            continue;
        }
        let tid = rec.reference_sequence_id().map(|t| t as i32).unwrap_or(-1);
        let pos = rec.alignment_start().map(usize::from).unwrap_or(0) as i64;
        let rev = flag & BAM_FREVERSE != 0;
        let mapq = rec.mapping_quality().map(u8::from).unwrap_or(0);
        let key = (tid, pos, rev);
        match best_per_group.get(&key) {
            Some(&idx) => {
                let prev_mapq = records[idx].mapping_quality().map(u8::from).unwrap_or(0);
                if mapq > prev_mapq {
                    keep[idx] = false;
                    keep[i] = true;
                    best_per_group.insert(key, i);
                }
            }
            None => {
                keep[i] = true;
                best_per_group.insert(key, i);
            }
        }
    }
    keep
}

trait BamLike {
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);

impl BamLike for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        self.0.write_record(header, record)
    }
}
impl BamLike for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        self.0.write_record(header, record)
    }
}

fn open_bam_output(out: Option<&Path>, header: &sam::Header) -> io::Result<Box<dyn BamLike>> {
    match out {
        Some(p) => {
            let mut writer = bam::io::Writer::new(File::create(p)?);
            writer.write_header(header)?;
            Ok(Box::new(BamFile(writer)))
        }
        None => {
            let mut writer = bam::io::Writer::new(io::stdout());
            writer.write_header(header)?;
            Ok(Box::new(BamStdout(writer)))
        }
    }
}

trait SamLike {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct SamFile(sam::io::Writer<File>);
struct SamStdout(sam::io::Writer<io::Stdout>);

impl SamLike for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        self.0.write_alignment_record(header, record)
    }
}
impl SamLike for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        self.0.write_alignment_record(header, record)
    }
}

fn open_sam_output(out: Option<&Path>, header: &sam::Header) -> io::Result<Box<dyn SamLike>> {
    match out {
        Some(p) => {
            let mut writer = sam::io::Writer::new(File::create(p)?);
            writer.write_header(header)?;
            Ok(Box::new(SamFile(writer)))
        }
        None => {
            let mut writer = sam::io::Writer::new(io::stdout());
            writer.write_header(header)?;
            Ok(Box::new(SamStdout(writer)))
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools rmdup [-sS] <in.bam|in.sam> [<out>]")?;
    writeln!(
        w,
        "  -s    treat reads as single-end (this port: single-end only)"
    )?;
    writeln!(w, "  -S    treat paired-end as single-end (alias of -s)")?;
    writeln!(w)?;
    writeln!(w, "NOTE: rmdup is deprecated; prefer `samtools markdup`.")?;
    Ok(())
}
