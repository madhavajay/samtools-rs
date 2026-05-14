//! `samtools markdup` — mark duplicate alignments.
//!
//! Mirrors `bam_markdup.c`. The upstream implementation is paired-aware,
//! barcode-aware, and supports optical-duplicate clustering. This initial
//! Rust port handles single-end and paired-end primary alignments for SAM
//! and BAM inputs. SE records are grouped by `(tid, alignment_start,
//! reverse-flag)`; PE records pair by qname and group on the canonical
//! combined position key. In each group the entry with the highest
//! (combined) MAPQ stays primary and the rest receive the `BAM_FDUP`
//! flag. Secondary and supplementary records inherit the duplicate flag from
//! duplicate primary records with the same query name.
//!
//! Supported flags:
//!  - `-r` — remove duplicates from the output (rather than just flagging).
//!  - `-s` — emit basic counts to stderr.
//!  - `-b TAG` / `--barcode-tag TAG` — include a string aux tag in the
//!    duplicate key.
//!  - `-O sam|bam` / `--output-fmt sam|bam` — output format (default `bam`).
//!  - `-o FILE` — output file (default stdout).
//!  - `--no-PG` — suppress the default samtools `@PG` line.
//!
//! Not yet supported: optical-duplicate distance (`-d`), full upstream stats
//! output, CRAM.

use std::collections::{HashMap, HashSet};
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
    alignment::{
        RecordBuf, io::Write as _, record::data::field::Tag, record_buf::data::field::Value,
    },
};

use crate::bam_flag::{
    BAM_FDUP, BAM_FMUNMAP, BAM_FPAIRED, BAM_FREVERSE, BAM_FSECONDARY, BAM_FSUPPLEMENTARY,
    BAM_FUNMAP,
};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

#[derive(Clone, Copy)]
enum OutFmt {
    Sam,
    Bam,
}

/// Entry point for `samtools markdup`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut output_fmt = OutFmt::Bam;
    let mut remove_dups = false;
    let mut emit_stats = false;
    let mut no_pg = false;
    let mut barcode_tag: Option<Tag> = None;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-r" => remove_dups = true,
            "-s" => emit_stats = true,
            "-O" | "--output-fmt" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("markdup", format!("missing value for {}", s));
                    return ExitCode::from(1);
                };
                output_fmt = match v.to_ascii_lowercase().as_str() {
                    "sam" => OutFmt::Sam,
                    "bam" => OutFmt::Bam,
                    _ => {
                        print_error("markdup", format!("unsupported output format \"{}\"", v));
                        return ExitCode::from(1);
                    }
                };
            }
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "-b" | "--barcode-tag" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("markdup", format!("missing value for {}", s));
                    return ExitCode::from(1);
                };
                barcode_tag = match parse_tag(v) {
                    Ok(tag) => Some(tag),
                    Err(e) => {
                        print_error("markdup", e);
                        return ExitCode::from(1);
                    }
                };
            }
            "--no-PG" => no_pg = true,
            "-@" | "--threads" | "-l" | "-m" | "-c" | "-d" | "-t" | "-T" => {
                // Accepted-but-ignored for compatibility.
                let _ = iter.next();
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error(
                    "markdup",
                    format!("option `{}` is not yet supported in samtools-rs markdup", s),
                );
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
            print_error("markdup", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam) {
        print_error(
            "markdup",
            "only SAM and BAM input are currently supported (CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    let pg_argv = if no_pg { None } else { Some(args) };
    let result = match format.exact {
        Exact::Sam => run_sam_markdup(
            &input,
            output.as_deref(),
            output_fmt,
            pg_argv,
            remove_dups,
            emit_stats,
            barcode_tag,
        ),
        Exact::Bam => run_bam_markdup(
            &input,
            output.as_deref(),
            output_fmt,
            pg_argv,
            remove_dups,
            emit_stats,
            barcode_tag,
        ),
        _ => unreachable!("format checked above"),
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("markdup", "markdup failed", &e);
            ExitCode::from(1)
        }
    }
}

struct MarkdupStats {
    duplicates: u64,
    examined: u64,
}

fn run_bam_markdup(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    remove_dups: bool,
    emit_stats: bool,
    barcode_tag: Option<Tag>,
) -> io::Result<()> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let mut header = reader.read_header()?;
    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }
    let mut records: Vec<RecordBuf> = Vec::new();
    let mut record = RecordBuf::default();
    loop {
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        records.push(record.clone());
    }
    let stats = mark_duplicates(&mut records, barcode_tag);
    let mut sink = open_output(output, fmt, &header)?;
    for rec in &records {
        if remove_dups && rec.flags().bits() as u32 & BAM_FDUP != 0 {
            continue;
        }
        sink.write_record(&header, rec)?;
    }
    if emit_stats {
        emit_basic_stats(&stats);
    }
    Ok(())
}

fn run_sam_markdup(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    remove_dups: bool,
    emit_stats: bool,
    barcode_tag: Option<Tag>,
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
    let stats = mark_duplicates(&mut records, barcode_tag);
    let mut sink = open_output(output, fmt, &header)?;
    for rec in &records {
        if remove_dups && rec.flags().bits() as u32 & BAM_FDUP != 0 {
            continue;
        }
        sink.write_record(&header, rec)?;
    }
    if emit_stats {
        emit_basic_stats(&stats);
    }
    Ok(())
}

/// Marks duplicates across single-end and paired-end primary alignments.
///
/// SE records (unpaired, or paired with the mate unmapped) are grouped by
/// `(tid, alignment_start, reverse_flag)` and the record with the highest
/// MAPQ in each group is kept as the primary.
///
/// PE records (paired + both mates mapped) are grouped into pairs by
/// qname. Each pair gets a canonical key combining both reads' `(tid,
/// pos, strand)` triples (sorted so the order does not matter), and the
/// pair with the highest summed MAPQ within a key group is kept; all
/// other records of duplicate pairs are flagged `BAM_FDUP`.
///
/// Secondary and supplementary alignments are not assessed for duplicates
/// themselves but inherit the dup flag from their primary if it gets one.
/// (Upstream's full algorithm matches this with extra qname-based linkage;
/// here we approximate by carrying the dup flag forward across records
/// with the same qname when at least one primary in that qname is dup.)
fn mark_duplicates(records: &mut [RecordBuf], barcode_tag: Option<Tag>) -> MarkdupStats {
    type BarcodeKey = Option<Vec<u8>>;
    type PosKey = (i32, i64, bool);
    type SeKey = (PosKey, BarcodeKey);
    type PairBarcodeKey = (BarcodeKey, BarcodeKey);
    type PairKey = (PosKey, PosKey, PairBarcodeKey);
    type PairIdx = (usize, usize);
    let mut se_best: HashMap<SeKey, usize> = HashMap::new();
    let mut pair_pending: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut pair_best: HashMap<PairKey, PairIdx> = HashMap::new();
    let mut duplicates = 0u64;
    let mut examined = 0u64;

    for i in 0..records.len() {
        let flag = records[i].flags().bits() as u32;
        if flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY) != 0 {
            continue;
        }
        examined += 1;
        let tid = records[i]
            .reference_sequence_id()
            .map(|t| t as i32)
            .unwrap_or(-1);
        let pos = records[i].alignment_start().map(usize::from).unwrap_or(0) as i64;
        let rev = flag & BAM_FREVERSE != 0;
        let mapq = records[i].mapping_quality().map(u8::from).unwrap_or(0);
        let me = (tid, pos, rev);

        let paired_both_mapped = flag & BAM_FPAIRED != 0 && flag & BAM_FMUNMAP == 0;
        if paired_both_mapped {
            // Look for the partner record by qname; pair them when both
            // ends have been seen. The first read of the pair sits in
            // `pair_pending`; the second read computes the pair key and
            // resolves dedup against any prior pair sharing the same key.
            let name = records[i].name().map(|n| n.to_vec()).unwrap_or_default();
            match pair_pending.remove(&name) {
                None => {
                    pair_pending.insert(name, i);
                }
                Some(first_idx) => {
                    let first_flag = records[first_idx].flags().bits() as u32;
                    let first_tid = records[first_idx]
                        .reference_sequence_id()
                        .map(|t| t as i32)
                        .unwrap_or(-1);
                    let first_pos = records[first_idx]
                        .alignment_start()
                        .map(usize::from)
                        .unwrap_or(0) as i64;
                    let first_rev = first_flag & BAM_FREVERSE != 0;
                    let first_mapq = records[first_idx]
                        .mapping_quality()
                        .map(u8::from)
                        .unwrap_or(0);
                    let first_barcode = barcode_value(&records[first_idx], barcode_tag);
                    let barcode = barcode_value(&records[i], barcode_tag);

                    let first = (first_tid, first_pos, first_rev);
                    let key = if first <= me {
                        (first, me, (first_barcode, barcode))
                    } else {
                        (me, first, (barcode, first_barcode))
                    };
                    let score = first_mapq as u32 + mapq as u32;

                    match pair_best.get(&key).copied() {
                        Some((prev_first, prev_second)) => {
                            let prev_score = records[prev_first]
                                .mapping_quality()
                                .map(u8::from)
                                .unwrap_or(0) as u32
                                + records[prev_second]
                                    .mapping_quality()
                                    .map(u8::from)
                                    .unwrap_or(0) as u32;
                            if score > prev_score {
                                set_dup(&mut records[prev_first]);
                                set_dup(&mut records[prev_second]);
                                pair_best.insert(key, (first_idx, i));
                            } else {
                                set_dup(&mut records[first_idx]);
                                set_dup(&mut records[i]);
                            }
                            duplicates += 2;
                        }
                        None => {
                            pair_best.insert(key, (first_idx, i));
                        }
                    }
                }
            }
            continue;
        }

        // SE path: mate not mapped (singleton) or unpaired primary.
        let se_key = (me, barcode_value(&records[i], barcode_tag));
        match se_best.get(&se_key).copied() {
            Some(idx) => {
                let prev_mapq = records[idx].mapping_quality().map(u8::from).unwrap_or(0);
                if mapq > prev_mapq {
                    set_dup(&mut records[idx]);
                    se_best.insert(se_key, i);
                } else {
                    set_dup(&mut records[i]);
                }
                duplicates += 1;
            }
            None => {
                se_best.insert(se_key, i);
            }
        }
    }

    let duplicate_primary_names: HashSet<Vec<u8>> = records
        .iter()
        .filter_map(|record| {
            let flag = record.flags().bits() as u32;
            (flag & BAM_FDUP != 0
                && flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY) == 0
                && flag & BAM_FUNMAP == 0)
                .then(|| record.name().map(|name| name.to_vec()))
                .flatten()
        })
        .collect();

    for record in records {
        let flag = record.flags().bits() as u32;
        if flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY) == 0 {
            continue;
        }
        let Some(name) = record.name() else {
            continue;
        };
        if duplicate_primary_names.contains(name) && set_dup(record) {
            duplicates += 1;
        }
    }

    MarkdupStats {
        duplicates,
        examined,
    }
}

fn set_dup(record: &mut RecordBuf) -> bool {
    let mut flags = record.flags();
    let was_duplicate = flags.contains(sam::alignment::record::Flags::DUPLICATE);
    flags.insert(sam::alignment::record::Flags::DUPLICATE);
    *record.flags_mut() = flags;
    !was_duplicate
}

fn parse_tag(s: &str) -> Result<Tag, String> {
    let bytes = s.as_bytes();
    if bytes.len() == 2 && bytes.iter().all(|b| b.is_ascii_alphanumeric()) {
        Ok(Tag::from([bytes[0], bytes[1]]))
    } else {
        Err(format!(
            "invalid aux tag \"{}\"; expected two alphanumeric characters",
            s
        ))
    }
}

fn barcode_value(record: &RecordBuf, tag: Option<Tag>) -> Option<Vec<u8>> {
    let tag = tag?;
    match record.data().get(&tag) {
        Some(Value::String(s)) | Some(Value::Hex(s)) => Some(s.to_vec()),
        Some(Value::Character(c)) => Some(vec![*c]),
        Some(value) => value.as_int().map(|n| n.to_string().into_bytes()),
        None => None,
    }
}

fn emit_basic_stats(stats: &MarkdupStats) {
    let _ = writeln!(
        io::stderr(),
        "READ: {}\nDUPLICATE TOTAL: {}",
        stats.examined,
        stats.duplicates
    );
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
        self.0.write_alignment_record(header, record)
    }
}
impl Sink for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        self.0.write_alignment_record(header, record)
    }
}
impl Sink for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        self.0.write_alignment_record(header, record)
    }
}
impl Sink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
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
    writeln!(w, "Usage: samtools markdup [options] <in.bam> [out.bam]")?;
    writeln!(w, "  -r            remove duplicate records")?;
    writeln!(w, "  -s            emit basic counts to stderr")?;
    writeln!(
        w,
        "  -b TAG        include barcode aux tag in duplicate key"
    )?;
    writeln!(w, "  -O sam|bam    output format (default: bam)")?;
    writeln!(w, "  -o FILE       output file (default stdout)")?;
    writeln!(w, "  --no-PG       do not add a @PG line")?;
    Ok(())
}
