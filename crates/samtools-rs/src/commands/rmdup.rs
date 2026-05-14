//! `samtools rmdup` — remove PCR duplicates (deprecated; `markdup` is preferred).
//!
//! Mirrors `bam_rmdup` / `bam_rmdupse` in upstream samtools. Upstream's
//! implementation is paired-aware and works on coordinate-sorted BAMs.
//!
//! This Rust port implements single-end and adjacent paired-end duplicate
//! removal: SE records are keyed by `(reference_sequence_id, alignment_start,
//! reverse-flag)`, while PE records are paired by qname and keyed by the
//! canonical pair of end coordinates. The record or pair with the highest
//! mapping quality score is retained.

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

use crate::bam_flag::{BAM_FMUNMAP, BAM_FPAIRED, BAM_FREVERSE, BAM_FUNMAP};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

/// Entry point for `samtools rmdup`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut output: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let mut no_pg = false;
    let mut single_end = false;
    let iter = args.iter().skip(1);
    for arg in iter {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-S" | "-s" => {
                single_end = true;
            }
            "--no-PG" => {
                no_pg = true;
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

    let pg_argv = if no_pg { None } else { Some(args) };
    let result = match format.exact {
        Exact::Sam => run_sam_rmdup(&input, output.as_deref(), pg_argv, single_end),
        Exact::Bam => run_bam_rmdup(&input, output.as_deref(), pg_argv, single_end),
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

fn run_bam_rmdup(
    input: &Path,
    output: Option<&Path>,
    pg_argv: Option<&[OsString]>,
    single_end: bool,
) -> io::Result<()> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let mut header = reader.read_header()?;
    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

    let mut records: Vec<RecordBuf> = Vec::new();
    let mut record = RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut record)?;
        if n == 0 {
            break;
        }
        records.push(record.clone());
    }

    let keep = duplicate_keep_mask_for_records(&records, single_end);
    let mut writer = open_bam_output(output, &header)?;
    for (i, rec) in records.iter().enumerate() {
        if keep[i] {
            writer.write_record(&header, rec)?;
        }
    }
    Ok(())
}

fn run_sam_rmdup(
    input: &Path,
    output: Option<&Path>,
    pg_argv: Option<&[OsString]>,
    single_end: bool,
) -> io::Result<()> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(input)?));
    let mut header = reader.read_header()?;
    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

    let mut records: Vec<RecordBuf> = Vec::new();
    loop {
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        records.push(record);
    }

    let keep = duplicate_keep_mask_for_records(&records, single_end);
    let mut writer = open_sam_output(output, &header)?;
    for (i, rec) in records.iter().enumerate() {
        if keep[i] {
            writer.write_record(&header, rec)?;
        }
    }
    Ok(())
}

fn duplicate_keep_mask_for_records(records: &[RecordBuf], single_end: bool) -> Vec<bool> {
    type PosKey = (i32, i64, bool);
    type PairKey = (PosKey, PosKey);
    type PairIdx = (usize, usize);
    let mut se_best: HashMap<PosKey, usize> = HashMap::new();
    let mut pair_pending: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut pair_best: HashMap<PairKey, PairIdx> = HashMap::new();
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
        let me = (tid, pos, rev);

        let paired_both_mapped = !single_end && flag & BAM_FPAIRED != 0 && flag & BAM_FMUNMAP == 0;
        if paired_both_mapped {
            let name = rec.name().map(|n| n.to_vec()).unwrap_or_default();
            match pair_pending.remove(&name) {
                None => {
                    pair_pending.insert(name, i);
                }
                Some(first_idx) => {
                    let first = pos_key(&records[first_idx]);
                    let key = if first <= me {
                        (first, me)
                    } else {
                        (me, first)
                    };
                    let score = pair_score(records, first_idx, i);
                    match pair_best.get(&key).copied() {
                        Some((prev_first, prev_second)) => {
                            let prev_score = pair_score(records, prev_first, prev_second);
                            if score > prev_score {
                                keep[prev_first] = false;
                                keep[prev_second] = false;
                                keep[first_idx] = true;
                                keep[i] = true;
                                pair_best.insert(key, (first_idx, i));
                            }
                        }
                        None => {
                            keep[first_idx] = true;
                            keep[i] = true;
                            pair_best.insert(key, (first_idx, i));
                        }
                    }
                }
            }
            continue;
        }

        match se_best.get(&me) {
            Some(&idx) => {
                let prev_mapq = records[idx].mapping_quality().map(u8::from).unwrap_or(0);
                if mapq > prev_mapq {
                    keep[idx] = false;
                    keep[i] = true;
                    se_best.insert(me, i);
                }
            }
            None => {
                keep[i] = true;
                se_best.insert(me, i);
            }
        }
    }

    for idx in pair_pending.into_values() {
        keep[idx] = true;
    }
    keep
}

trait BamLike {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);

impl BamLike for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        self.0.write_alignment_record(header, record)
    }
}
impl BamLike for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        self.0.write_alignment_record(header, record)
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

fn pos_key(record: &RecordBuf) -> (i32, i64, bool) {
    let flag = record.flags().bits() as u32;
    (
        record
            .reference_sequence_id()
            .map(|t| t as i32)
            .unwrap_or(-1),
        record.alignment_start().map(usize::from).unwrap_or(0) as i64,
        flag & BAM_FREVERSE != 0,
    )
}

fn pair_score(records: &[RecordBuf], first: usize, second: usize) -> u32 {
    u32::from(records[first].mapping_quality().map(u8::from).unwrap_or(0))
        + u32::from(records[second].mapping_quality().map(u8::from).unwrap_or(0))
}
