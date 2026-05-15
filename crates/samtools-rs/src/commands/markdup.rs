//! `samtools markdup` — mark duplicate alignments.
//!
//! Faithful port of `bam_markdup.c`'s primary duplicate detection. Reads
//! are processed in coordinate order; paired reads build the upstream
//! `make_pair_key` (template / `--mode s` sequence) and a shared
//! `make_single_key`, and the kept read of a colliding key is the one
//! with the higher `calc_score (+ ms)` (sum of base quals ≥ 15 plus the
//! mate-score tag), ties broken by qname. The left/right (`R_LE`/`R_RI`)
//! component keeps the two ends of a template distinct so only
//! corresponding mates of duplicate templates collide. With `-S`,
//! duplicate reads carrying `SA`/`XA` or an unmapped mate seed a
//! qname `dup_hash` that flags matching supplementary/secondary/unmapped
//! records. Byte-exact vs upstream `markdup/{5,6,7,13}` (template,
//! sequence, supplementary, barcode-tag). **Not yet:** optical-duplicate
//! chain re-tagging (`find_duplicate_chains`), `--read-coords` /
//! `--coords-order` / `--barcode-rgx` / `--barcode-name` /
//! `--use-read-groups` / `--duplicate-count`, exact `-s` stats, CRAM.
//!
//! Supported flags:
//!  - `-r` — remove duplicates from the output (rather than just flagging).
//!  - `-s` — emit upstream-shaped summary counts to stderr.
//!  - `-S` — accepted; supplementary propagation is always performed.
//!  - `-c` — clear existing duplicate flags and duplicate metadata tags before
//!    marking.
//!  - `-t` — add `do:Z:<original>` duplicate-origin tags for duplicates.
//!  - `-d DISTANCE` — add `dt:Z:SQ` / `dt:Z:LB` duplicate-type tags using
//!    Illumina-style read-name tile/x/y optical-distance checks.
//!  - `--include-fails` — include QCFAIL reads in duplicate marking.
//!  - `-m t|s` / `--mode t|s` — accepted duplicate-decision mode selector.
//!  - `-b TAG` / `--barcode-tag TAG` — include a string aux tag in the
//!    duplicate key.
//!  - `-O sam|bam` / `--output-fmt sam|bam` — output format (default `bam`).
//!  - `-o FILE` — output file (default stdout).
//!  - `--no-PG` — suppress the default samtools `@PG` line.
//!
//! Not yet supported: exact upstream stats output, CRAM.

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use bstr::BString;
use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{
    self,
    alignment::{
        RecordBuf, io::Write as _, record::cigar::op::Kind, record::data::field::Tag,
        record_buf::data::field::Value,
    },
};

use crate::bam_flag::{
    BAM_FDUP, BAM_FMUNMAP, BAM_FPAIRED, BAM_FQCFAIL, BAM_FREVERSE, BAM_FSECONDARY,
    BAM_FSUPPLEMENTARY, BAM_FUNMAP,
};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

#[derive(Clone, Copy)]
enum OutFmt {
    Sam,
    Bam,
}

/// Duplicate-decision mode (`-m`/`--mode`). Upstream default is
/// `MD_MODE_TEMPLATE`; `s` selects `MD_MODE_SEQUENCE`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum DupMode {
    Template,
    Sequence,
}

#[derive(Clone, Copy)]
struct MarkdupOptions {
    remove_dups: bool,
    emit_stats: bool,
    clear_existing_dups: bool,
    duplicate_origin_tag: bool,
    optical_distance: Option<u32>,
    include_fails: bool,
    barcode_tag: Option<Tag>,
    mode: DupMode,
    supp: bool,
}

/// Entry point for `samtools markdup`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut output_fmt = OutFmt::Bam;
    let mut remove_dups = false;
    let mut emit_stats = false;
    let mut clear_existing_dups = false;
    let mut duplicate_origin_tag = false;
    let mut optical_distance = None;
    let mut include_fails = false;
    let mut no_pg = false;
    let mut barcode_tag: Option<Tag> = None;
    let mut mode = DupMode::Template;
    let mut supp = false;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-r" => remove_dups = true,
            "-s" => emit_stats = true,
            "-S" => supp = true,
            "--include-fails" => include_fails = true,
            "-c" => clear_existing_dups = true,
            "-t" => duplicate_origin_tag = true,
            "-m" | "--mode" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("markdup", format!("missing value for {}", s));
                    return ExitCode::from(1);
                };
                mode = match v {
                    "t" => DupMode::Template,
                    "s" => DupMode::Sequence,
                    _ => {
                        print_error("markdup", format!("unknown mode '{}'", v));
                        return ExitCode::from(1);
                    }
                };
            }
            "-d" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("markdup", "missing value for -d");
                    return ExitCode::from(1);
                };
                match v.parse::<u32>() {
                    Ok(distance) => optical_distance = Some(distance),
                    Err(_) => {
                        print_error("markdup", format!("invalid optical distance \"{}\"", v));
                        return ExitCode::from(1);
                    }
                }
            }
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
            "-@" | "--threads" | "-l" | "-T" => {
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
                } else if output.is_none() && s != "-" {
                    // A `-` output operand means stdout (output stays None).
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
    let options = MarkdupOptions {
        remove_dups,
        emit_stats,
        clear_existing_dups,
        duplicate_origin_tag,
        optical_distance,
        include_fails,
        barcode_tag,
        mode,
        supp,
    };
    let result = match format.exact {
        Exact::Sam => run_sam_markdup(&input, output.as_deref(), output_fmt, pg_argv, options),
        Exact::Bam => run_bam_markdup(&input, output.as_deref(), output_fmt, pg_argv, options),
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
    read: u64,
    written: u64,
    excluded: u64,
    examined: u64,
    paired: u64,
    single: u64,
    duplicate_pair: u64,
    duplicate_single: u64,
    duplicate_pair_optical: u64,
    duplicate_single_optical: u64,
    duplicate_non_primary: u64,
    duplicate_non_primary_optical: u64,
}

#[derive(Clone, Copy)]
enum DuplicateType {
    Library,
    Optical,
}

impl DuplicateType {
    fn tag_value(self) -> &'static [u8] {
        match self {
            Self::Library => b"LB",
            Self::Optical => b"SQ",
        }
    }

    fn is_optical(self) -> bool {
        matches!(self, Self::Optical)
    }
}

struct DuplicateMetadata {
    origin: Vec<u8>,
    duplicate_type: Option<DuplicateType>,
}

fn run_bam_markdup(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    options: MarkdupOptions,
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
    if options.clear_existing_dups {
        clear_duplicate_marks(&mut records);
    }
    let mut stats = mark_duplicates(&mut records, options);
    stats.written = output_record_count(&records, options.remove_dups);
    let mut sink = open_output(output, fmt, &header)?;
    for rec in &records {
        if options.remove_dups && rec.flags().bits() as u32 & BAM_FDUP != 0 {
            continue;
        }
        sink.write_record(&header, rec)?;
    }
    if options.emit_stats {
        write_markdup_stats(&mut io::stderr().lock(), &stats)?;
    }
    Ok(())
}

fn run_sam_markdup(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    options: MarkdupOptions,
) -> io::Result<()> {
    let mut reader = crate::sam_compat::open_sam_reader_tolerant(input)?;
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
    if options.clear_existing_dups {
        clear_duplicate_marks(&mut records);
    }
    let mut stats = mark_duplicates(&mut records, options);
    stats.written = output_record_count(&records, options.remove_dups);
    let mut sink = open_output(output, fmt, &header)?;
    for rec in &records {
        if options.remove_dups && rec.flags().bits() as u32 & BAM_FDUP != 0 {
            continue;
        }
        sink.write_record(&header, rec)?;
    }
    if options.emit_stats {
        write_markdup_stats(&mut io::stderr().lock(), &stats)?;
    }
    Ok(())
}

// ---- Upstream-faithful duplicate detection (bam_markdup.c port) ----

const MD_MIN_QUALITY: u8 = 15;
const O_FF: i8 = 2;
const O_RR: i8 = 3;
const O_FR: i8 = 5;
const O_RF: i8 = 7;
const R_LE: i8 = 11;
const R_RI: i8 = 13;
const BAM_FREAD1: u32 = 0x40;
const BAM_FMREVERSE: u32 = 0x20;

/// Mirror of upstream `key_data_t`. For single keys the `other_*`/`leftmost`
/// fields are forced to 0 so derived `Eq`/`Hash` ignore them, matching
/// `key_equal`'s `!a.single` guard.
#[derive(Clone, PartialEq, Eq, Hash)]
struct KeyData {
    single: bool,
    this_ref: i32,
    this_coord: i64,
    other_ref: i32,
    other_coord: i64,
    leftmost: i8,
    orientation: i8,
    barcode: Vec<u8>,
    read_group: i32,
}

fn rec_flags(r: &RecordBuf) -> u32 {
    r.flags().bits() as u32
}
fn is_rev(r: &RecordBuf) -> bool {
    rec_flags(r) & BAM_FREVERSE != 0
}
fn is_mrev(r: &RecordBuf) -> bool {
    rec_flags(r) & BAM_FMREVERSE != 0
}
fn is_read1(r: &RecordBuf) -> bool {
    rec_flags(r) & BAM_FREAD1 != 0
}
fn rec_tid(r: &RecordBuf) -> i32 {
    r.reference_sequence_id().map(|t| t as i32).unwrap_or(-1)
}
fn rec_mtid(r: &RecordBuf) -> i32 {
    r.mate_reference_sequence_id()
        .map(|t| t as i32)
        .unwrap_or(-1)
}
/// 1-based alignment start (upstream `core.pos + 1`); 0 if unmapped.
fn pos1(r: &RecordBuf) -> i64 {
    r.alignment_start().map(|p| p.get() as i64).unwrap_or(0)
}
/// 0-based mate position (upstream `core.mpos`); -1 if absent.
fn mpos0(r: &RecordBuf) -> i64 {
    r.mate_alignment_start()
        .map(|p| p.get() as i64 - 1)
        .unwrap_or(-1)
}
fn rec_name(r: &RecordBuf) -> Vec<u8> {
    r.name().map(|n| n.to_vec()).unwrap_or_default()
}

fn clip_lead(cigar: &sam::alignment::record_buf::Cigar) -> i64 {
    let mut n = 0i64;
    for op in cigar.as_ref() {
        match op.kind() {
            Kind::SoftClip | Kind::HardClip => n += op.len() as i64,
            _ => break,
        }
    }
    n
}
fn clip_trail(cigar: &sam::alignment::record_buf::Cigar) -> i64 {
    let mut n = 0i64;
    for op in cigar.as_ref().iter().rev() {
        match op.kind() {
            Kind::SoftClip | Kind::HardClip => n += op.len() as i64,
            _ => break,
        }
    }
    n
}
fn ref_len(cigar: &sam::alignment::record_buf::Cigar) -> i64 {
    let mut n = 0i64;
    for op in cigar.as_ref() {
        match op.kind() {
            Kind::Match
            | Kind::Deletion
            | Kind::Skip
            | Kind::SequenceMatch
            | Kind::SequenceMismatch => n += op.len() as i64,
            _ => {}
        }
    }
    n
}

/// `unclipped_start`: `core.pos - leading_clip + 1`.
fn unclipped_start(r: &RecordBuf) -> i64 {
    pos1(r) - clip_lead(r.cigar())
}
/// `unclipped_end`: `bam_endpos + trailing_clip`, where
/// `bam_endpos = core.pos + ref_len`.
fn unclipped_end(r: &RecordBuf) -> i64 {
    pos1(r) - 1 + ref_len(r.cigar()) + clip_trail(r.cigar())
}

/// Iterate `<num><op>` tokens of an MC-style CIGAR string. A non-digit run
/// yields `num = 1` (mirrors upstream `strtol` fallback).
fn mc_tokens(mc: &str) -> Vec<(i64, u8)> {
    let bytes = mc.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() && bytes[i] != b'*' {
        let mut num: i64 = 0;
        if bytes[i].is_ascii_digit() {
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                num = num * 10 + (bytes[i] - b'0') as i64;
                i += 1;
            }
        } else {
            num = 1;
        }
        if i >= bytes.len() {
            break;
        }
        let op = bytes[i];
        out.push((num, op));
        i += 1;
    }
    out
}
/// `unclipped_other_start(op, cigar)` = `op - leading_clip + 1`.
fn unclipped_other_start(op: i64, mc: &str) -> i64 {
    let mut clipped = 0i64;
    for (num, c) in mc_tokens(mc) {
        if c == b'S' || c == b'H' {
            clipped += num;
        } else {
            break;
        }
    }
    op - clipped + 1
}
/// `unclipped_other_end(op, cigar)` = `op + ref_consumed` where reference
/// consumption counts M/D/N/=/X and post-leading-clip S/H.
fn unclipped_other_end(op: i64, mc: &str) -> i64 {
    let mut refpos = 0i64;
    let mut skip = true;
    for (num, c) in mc_tokens(mc) {
        match c {
            b'M' | b'D' | b'N' | b'=' | b'X' => {
                refpos += num;
                skip = false;
            }
            b'S' | b'H' if !skip => refpos += num,
            _ => {}
        }
    }
    op + refpos
}

/// Sum of base qualities `>= MD_MIN_QUALITY` (upstream `calc_score`).
fn calc_score(r: &RecordBuf) -> i64 {
    r.quality_scores()
        .as_ref()
        .iter()
        .filter(|&&q| q >= MD_MIN_QUALITY)
        .map(|&q| q as i64)
        .sum()
}
/// Mate score from the `ms` aux tag (set by `fixmate -m`).
fn mate_score(r: &RecordBuf) -> i64 {
    r.data()
        .get(&Tag::from([b'm', b's']))
        .and_then(|v| v.as_int())
        .unwrap_or(0)
}

fn barcode_bytes(r: &RecordBuf, tag: Option<Tag>) -> Vec<u8> {
    barcode_value(r, tag).unwrap_or_default()
}

fn mc_string(r: &RecordBuf) -> Option<String> {
    match r.data().get(&Tag::from([b'M', b'C']))? {
        Value::String(s) => Some(String::from_utf8_lossy(s).into_owned()),
        _ => None,
    }
}

/// has_mate: paired, mate mapped, and mate coordinates present.
fn has_mate(r: &RecordBuf) -> bool {
    let f = rec_flags(r);
    f & BAM_FPAIRED != 0 && f & BAM_FMUNMAP == 0 && !(rec_mtid(r) == -1 && mpos0(r) == -1)
}

#[allow(clippy::collapsible_else_if)]
fn make_pair_key(r: &RecordBuf, mode: DupMode, barcode_tag: Option<Tag>, rg: i32) -> KeyData {
    let this_ref = rec_tid(r) + 1;
    let other_ref = rec_mtid(r) + 1;
    let mut this_coord = unclipped_start(r);
    let this_end = unclipped_end(r);
    let mc = mc_string(r).unwrap_or_default();
    let mp = mpos0(r);
    let mut other_coord = unclipped_other_start(mp, &mc);
    let other_end = unclipped_other_end(mp, &mc);
    let rev = is_rev(r);
    let mrev = is_mrev(r);
    let read1 = is_read1(r);
    let orientation: i8;
    let leftmost: bool;

    if mode == DupMode::Template {
        let lm = if this_ref != other_ref {
            this_ref < other_ref
        } else if rev == mrev {
            if !rev {
                this_coord <= other_coord
            } else {
                this_end <= other_end
            }
        } else if rev {
            this_end <= other_coord
        } else {
            this_coord <= other_end
        };
        leftmost = lm;

        if lm {
            if rev == mrev {
                other_coord = other_end;
                orientation = if !rev {
                    if read1 { O_FF } else { O_RR }
                } else {
                    if read1 { O_RR } else { O_FF }
                };
            } else if !rev {
                orientation = O_FR;
                other_coord = other_end;
            } else {
                orientation = O_RF;
                this_coord = this_end;
            }
        } else {
            if rev == mrev {
                this_coord = this_end;
                orientation = if !rev {
                    if read1 { O_RR } else { O_FF }
                } else {
                    if read1 { O_FF } else { O_RR }
                };
            } else if !rev {
                orientation = O_RF;
                other_coord = other_end;
            } else {
                orientation = O_FR;
                this_coord = this_end;
            }
        }
    } else {
        // MD_MODE_SEQUENCE
        let diff: i64 = if this_ref != other_ref {
            (this_ref - other_ref) as i64
        } else if rev == mrev {
            if !rev {
                this_coord - other_coord
            } else {
                this_end - other_end
            }
        } else if rev {
            this_end - other_coord
        } else {
            this_coord - other_end
        };
        let pos0 = pos1(r) - 1;
        leftmost = if diff < 0 {
            true
        } else if diff > 0 {
            false
        } else if pos0 == mp {
            read1
        } else {
            pos0 < mp
        };

        orientation = if leftmost {
            if rev == mrev {
                if !rev { O_FF } else { O_RR }
            } else if !rev {
                O_FR
            } else {
                O_RF
            }
        } else if rev == mrev {
            if !rev { O_RR } else { O_FF }
        } else if !rev {
            O_RF
        } else {
            O_FR
        };

        this_coord = if !rev {
            unclipped_start(r)
        } else {
            unclipped_end(r)
        };
        other_coord = if !mrev {
            unclipped_other_start(mp, &mc)
        } else {
            unclipped_other_end(mp, &mc)
        };
    }

    let left_read = if !leftmost { R_RI } else { R_LE };

    KeyData {
        single: false,
        this_ref,
        this_coord,
        other_ref,
        other_coord,
        leftmost: left_read,
        orientation,
        barcode: barcode_bytes(r, barcode_tag),
        read_group: rg,
    }
}

fn make_single_key(r: &RecordBuf, barcode_tag: Option<Tag>, rg: i32) -> KeyData {
    let this_ref = rec_tid(r) + 1;
    let (this_coord, orientation) = if is_rev(r) {
        (unclipped_end(r), O_RR)
    } else {
        (unclipped_start(r), O_FF)
    };
    KeyData {
        single: true,
        this_ref,
        this_coord,
        other_ref: 0,
        other_coord: 0,
        leftmost: 0,
        orientation,
        barcode: barcode_bytes(r, barcode_tag),
        read_group: rg,
    }
}

/// Paired score = `calc_score + ms`, with QCFAIL-asymmetry override.
fn pair_scores(stored: &RecordBuf, incoming: &RecordBuf) -> (i64, i64) {
    let qf_s = rec_flags(stored) & BAM_FQCFAIL != 0;
    let qf_i = rec_flags(incoming) & BAM_FQCFAIL != 0;
    if qf_s != qf_i {
        if qf_s { (0, 1) } else { (1, 0) }
    } else {
        (
            calc_score(stored) + mate_score(stored),
            calc_score(incoming) + mate_score(incoming),
        )
    }
}

/// Faithful port of `bam_markdup.c`'s primary duplicate decision.
///
/// Reads are processed in input (coordinate) order. Paired reads build an
/// upstream `make_pair_key` (template/sequence mode) plus a `make_single_key`
/// kept in a shared single hash so a true singleton always loses to a pair.
/// The kept read of a colliding key is the one with the higher
/// `calc_score (+ mate ms)`, breaking ties by qname (`strcmp`). Each
/// mate's key encodes left/right (`R_LE`/`R_RI`) so the two ends of one
/// template get distinct keys and only corresponding mates of duplicate
/// templates collide. Secondary/supplementary records inherit the flag
/// from their duplicate primary by qname.
fn mark_duplicates(records: &mut [RecordBuf], options: MarkdupOptions) -> MarkdupStats {
    let mut single_hash: HashMap<KeyData, usize> = HashMap::new();
    let mut pair_hash: HashMap<KeyData, usize> = HashMap::new();
    let mut duplicate_primary_metadata: HashMap<Vec<u8>, DuplicateMetadata> = HashMap::new();
    // Upstream `dup_hash`: qnames of marked-duplicate reads that carry an
    // `SA`/`XA` tag or have an unmapped mate; only built with `-S`.
    let mut supp_dups: HashMap<Vec<u8>, DuplicateMetadata> = HashMap::new();
    let mut stats = MarkdupStats {
        read: records.len() as u64,
        written: records.len() as u64,
        excluded: 0,
        examined: 0,
        paired: 0,
        single: 0,
        duplicate_pair: 0,
        duplicate_single: 0,
        duplicate_pair_optical: 0,
        duplicate_single_optical: 0,
        duplicate_non_primary: 0,
        duplicate_non_primary_optical: 0,
    };

    let opt = options.optical_distance;
    let do_tag = options.duplicate_origin_tag;

    let supp = options.supp;
    // Marks `dup_idx` as a duplicate of `orig_idx`, recording origin/type.
    // Mirrors `mark_duplicates`: with `-S`, a dup read carrying `SA`/`XA`
    // or an unmapped mate seeds `dup_hash` (first qname wins).
    let mark = |records: &mut [RecordBuf],
                duplicate_primary_metadata: &mut HashMap<Vec<u8>, DuplicateMetadata>,
                supp_dups: &mut HashMap<Vec<u8>, DuplicateMetadata>,
                dup_idx: usize,
                orig_idx: usize|
     -> Option<DuplicateType> {
        let origin = rec_name(&records[orig_idx]);
        let dtype = duplicate_type(&records[dup_idx], &records[orig_idx], opt);
        mark_duplicate(&mut records[dup_idx], Some(&origin), do_tag, dtype);
        remember_duplicate_metadata(
            duplicate_primary_metadata,
            &records[dup_idx],
            &origin,
            dtype,
        );
        if supp {
            let dup = &records[dup_idx];
            let mate_unmapped = rec_flags(dup) & BAM_FMUNMAP != 0;
            let has_sa = dup.data().get(&Tag::from([b'S', b'A'])).is_some();
            let has_xa = dup.data().get(&Tag::from([b'X', b'A'])).is_some();
            if mate_unmapped || has_sa || has_xa {
                supp_dups.entry(rec_name(dup)).or_insert(DuplicateMetadata {
                    origin: origin.clone(),
                    duplicate_type: dtype,
                });
            }
        }
        dtype
    };

    for i in 0..records.len() {
        let flag = rec_flags(&records[i]);
        if flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY) != 0 {
            continue;
        }
        if flag & BAM_FQCFAIL != 0 && !options.include_fails {
            stats.excluded += 1;
            continue;
        }
        stats.examined += 1;

        if has_mate(&records[i]) {
            let pair_key = make_pair_key(&records[i], options.mode, options.barcode_tag, 0);
            let single_key = make_single_key(&records[i], options.barcode_tag, 0);
            stats.paired += 1;

            // Single hash: a true singleton already stored loses to this pair.
            match single_hash.get(&single_key).copied() {
                None => {
                    single_hash.insert(single_key.clone(), i);
                }
                Some(j) => {
                    if !has_mate(&records[j]) {
                        let dtype = mark(
                            records,
                            &mut duplicate_primary_metadata,
                            &mut supp_dups,
                            j,
                            i,
                        );
                        stats.duplicate_single += 1;
                        if dtype.is_some_and(DuplicateType::is_optical) {
                            stats.duplicate_single_optical += 1;
                        }
                        single_hash.insert(single_key.clone(), i);
                    }
                }
            }

            // Pair hash: corresponding mates of duplicate templates collide.
            match pair_hash.get(&pair_key).copied() {
                None => {
                    pair_hash.insert(pair_key, i);
                }
                Some(j) => {
                    let (old_s, new_s) = pair_scores(&records[j], &records[i]);
                    let tie = if new_s == old_s {
                        if rec_name(&records[i]) < rec_name(&records[j]) {
                            1
                        } else {
                            -1
                        }
                    } else {
                        0
                    };
                    let (dup_idx, orig_idx, swap) = if new_s + tie > old_s {
                        (j, i, true)
                    } else {
                        (i, j, false)
                    };
                    let dtype = mark(
                        records,
                        &mut duplicate_primary_metadata,
                        &mut supp_dups,
                        dup_idx,
                        orig_idx,
                    );
                    if swap {
                        pair_hash.insert(pair_key, i);
                    }
                    stats.duplicate_pair += 1;
                    if dtype.is_some_and(DuplicateType::is_optical) {
                        stats.duplicate_pair_optical += 1;
                    }
                }
            }
        } else {
            // Single (or effectively single) reads.
            let single_key = make_single_key(&records[i], options.barcode_tag, 0);
            stats.single += 1;
            match single_hash.get(&single_key).copied() {
                None => {
                    single_hash.insert(single_key, i);
                }
                Some(j) => {
                    if has_mate(&records[j]) {
                        let dtype = mark(
                            records,
                            &mut duplicate_primary_metadata,
                            &mut supp_dups,
                            i,
                            j,
                        );
                        stats.duplicate_single += 1;
                        if dtype.is_some_and(DuplicateType::is_optical) {
                            stats.duplicate_single_optical += 1;
                        }
                    } else {
                        let old_s = calc_score(&records[j]);
                        let new_s = calc_score(&records[i]);
                        let (dup_idx, orig_idx, swap) = if new_s > old_s {
                            (j, i, true)
                        } else {
                            (i, j, false)
                        };
                        let dtype = mark(
                            records,
                            &mut duplicate_primary_metadata,
                            &mut supp_dups,
                            dup_idx,
                            orig_idx,
                        );
                        if swap {
                            single_hash.insert(single_key, i);
                        }
                        stats.duplicate_single += 1;
                        if dtype.is_some_and(DuplicateType::is_optical) {
                            stats.duplicate_single_optical += 1;
                        }
                    }
                }
            }
        }
    }

    // Upstream supplementary pass (`-S` only): supplementary, secondary, and
    // unmapped records whose qname seeded `dup_hash` inherit the flag.
    if options.supp {
        for record in records {
            let flag = record.flags().bits() as u32;
            if flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP) == 0 {
                continue;
            }
            let Some(name) = record.name() else {
                continue;
            };
            if let Some(metadata) = supp_dups.get(name)
                && mark_duplicate(
                    record,
                    Some(metadata.origin.as_slice()),
                    options.duplicate_origin_tag,
                    metadata.duplicate_type,
                )
            {
                stats.duplicate_non_primary += 1;
                if metadata
                    .duplicate_type
                    .is_some_and(DuplicateType::is_optical)
                {
                    stats.duplicate_non_primary_optical += 1;
                }
            }
        }
    }

    stats
}

fn output_record_count(records: &[RecordBuf], remove_dups: bool) -> u64 {
    records
        .iter()
        .filter(|record| !remove_dups || record.flags().bits() as u32 & BAM_FDUP == 0)
        .count() as u64
}

fn clear_duplicate_marks(records: &mut [RecordBuf]) {
    let duplicate_origin_tag = Tag::from([b'd', b'o']);
    let duplicate_type_tag = Tag::from([b'd', b't']);
    for record in records {
        let mut flags = record.flags();
        flags.remove(sam::alignment::record::Flags::DUPLICATE);
        *record.flags_mut() = flags;
        record.data_mut().remove(&duplicate_origin_tag);
        record.data_mut().remove(&duplicate_type_tag);
    }
}

fn remember_duplicate_metadata(
    metadata: &mut HashMap<Vec<u8>, DuplicateMetadata>,
    record: &RecordBuf,
    origin: &[u8],
    duplicate_type: Option<DuplicateType>,
) {
    if let Some(name) = record.name() {
        metadata.insert(
            name.to_vec(),
            DuplicateMetadata {
                origin: origin.to_vec(),
                duplicate_type,
            },
        );
    }
}

fn mark_duplicate(
    record: &mut RecordBuf,
    origin: Option<&[u8]>,
    add_origin_tag: bool,
    duplicate_type: Option<DuplicateType>,
) -> bool {
    let was_new = set_dup(record);
    if add_origin_tag && let Some(origin) = origin {
        record.data_mut().insert(
            Tag::from([b'd', b'o']),
            Value::String(BString::from(origin.to_vec())),
        );
    }
    if let Some(duplicate_type) = duplicate_type {
        record.data_mut().insert(
            Tag::from([b'd', b't']),
            Value::String(BString::from(duplicate_type.tag_value().to_vec())),
        );
    }
    was_new
}

fn duplicate_type(
    duplicate: &RecordBuf,
    original: &RecordBuf,
    optical_distance: Option<u32>,
) -> Option<DuplicateType> {
    let distance = optical_distance?;
    Some(if is_optical_duplicate(duplicate, original, distance) {
        DuplicateType::Optical
    } else {
        DuplicateType::Library
    })
}

fn is_optical_duplicate(duplicate: &RecordBuf, original: &RecordBuf, distance: u32) -> bool {
    let (Some(dn), Some(on)) = (duplicate.name(), original.name()) else {
        return false;
    };
    let dn = dn.as_ref();
    let on = on.as_ref();
    let Some((d_beg, d_end, dx, dy)) = get_coordinates_colons(dn) else {
        return false;
    };
    let Some((o_beg, o_end, ox, oy)) = get_coordinates_colons(on) else {
        return false;
    };
    let o_len = o_end - o_beg;
    let d_len = d_end - d_beg;
    o_len == d_len
        && on[o_beg..o_end] == dn[d_beg..d_end]
        && ox.abs_diff(dx) <= distance as u64
        && oy.abs_diff(dy) <= distance as u64
}

/// Port of `get_coordinates_colons`: from an Illumina-style read name,
/// pick x/y by colon-separator count and return
/// `(tile_beg, tile_end, x, y)` where `name[tile_beg..tile_end]` is the
/// prefix compared for string equality. `None` if undecipherable.
fn get_coordinates_colons(qname: &[u8]) -> Option<(usize, usize, i64, i64)> {
    let mut sep = 0;
    let mut xpos = 0usize;
    let mut ypos = 0usize;
    for (pos, &c) in qname.iter().enumerate() {
        if c == b':' {
            sep += 1;
            match sep {
                2 => xpos = pos + 1,
                3 => ypos = pos + 1,
                4 => {
                    xpos = ypos;
                    ypos = pos + 1;
                }
                5 => xpos = pos + 1,
                6 => ypos = pos + 1,
                _ => {}
            }
        }
    }
    if !(sep == 3 || sep == 4 || sep == 6 || sep == 7) {
        return None;
    }
    let x = parse_strtol(qname, xpos)?;
    let y = parse_strtol(qname, ypos)?;
    Some((0, xpos, x, y))
}

/// `strtol`-style: parse a base-10 integer at `start`; `None` if no
/// digits are consumed (mirrors `(qname+pos) == end`).
fn parse_strtol(bytes: &[u8], start: usize) -> Option<i64> {
    let mut i = start;
    let mut neg = false;
    if i < bytes.len() && (bytes[i] == b'-' || bytes[i] == b'+') {
        neg = bytes[i] == b'-';
        i += 1;
    }
    let digit_start = i;
    let mut v: i64 = 0;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        v = v * 10 + (bytes[i] - b'0') as i64;
        i += 1;
    }
    if i == digit_start {
        return None;
    }
    Some(if neg { -v } else { v })
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

fn write_markdup_stats(mut w: impl Write, stats: &MarkdupStats) -> io::Result<()> {
    let duplicate_primary_total = stats.duplicate_pair + stats.duplicate_single;
    let duplicate_total = duplicate_primary_total + stats.duplicate_non_primary;
    let estimated_library_size = estimate_library_size(
        stats.paired,
        stats.duplicate_pair,
        stats.duplicate_pair_optical,
    );
    writeln!(w, "READ: {}", stats.read)?;
    writeln!(w, "WRITTEN: {}", stats.written)?;
    writeln!(w, "EXCLUDED: {}", stats.excluded)?;
    writeln!(w, "EXAMINED: {}", stats.examined)?;
    writeln!(w, "PAIRED: {}", stats.paired)?;
    writeln!(w, "SINGLE: {}", stats.single)?;
    writeln!(w, "DUPLICATE PAIR: {}", stats.duplicate_pair)?;
    writeln!(w, "DUPLICATE SINGLE: {}", stats.duplicate_single)?;
    writeln!(
        w,
        "DUPLICATE PAIR OPTICAL: {}",
        stats.duplicate_pair_optical
    )?;
    writeln!(
        w,
        "DUPLICATE SINGLE OPTICAL: {}",
        stats.duplicate_single_optical
    )?;
    writeln!(w, "DUPLICATE NON PRIMARY: {}", stats.duplicate_non_primary)?;
    writeln!(
        w,
        "DUPLICATE NON PRIMARY OPTICAL: {}",
        stats.duplicate_non_primary_optical
    )?;
    writeln!(w, "DUPLICATE PRIMARY TOTAL: {duplicate_primary_total}")?;
    writeln!(w, "DUPLICATE TOTAL: {duplicate_total}")?;
    writeln!(w, "ESTIMATED_LIBRARY_SIZE: {estimated_library_size}")?;
    Ok(())
}

fn estimate_library_size(paired_reads: u64, paired_duplicate_reads: u64, optical: u64) -> u64 {
    let non_optical_pairs = paired_reads.saturating_sub(optical) / 2;
    let unique_pairs = paired_reads.saturating_sub(paired_duplicate_reads) / 2;
    let duplicate_pairs = paired_duplicate_reads.saturating_sub(optical) / 2;

    if non_optical_pairs == 0
        || duplicate_pairs == 0
        || unique_pairs == 0
        || non_optical_pairs <= duplicate_pairs
    {
        return 0;
    }

    let unique_pairs_f = unique_pairs as f64;
    let non_optical_pairs_f = non_optical_pairs as f64;
    let mut lower = 1.0;
    let mut upper = 100.0;

    if coverage_equation(lower * unique_pairs_f, unique_pairs_f, non_optical_pairs_f) < 0.0 {
        return 0;
    }

    while coverage_equation(upper * unique_pairs_f, unique_pairs_f, non_optical_pairs_f) > 0.0 {
        upper *= 10.0;
    }

    for _ in 0..40 {
        let midpoint = (lower + upper) / 2.0;
        let value = coverage_equation(
            midpoint * unique_pairs_f,
            unique_pairs_f,
            non_optical_pairs_f,
        );
        if value > 0.0 {
            lower = midpoint;
        } else if value < 0.0 {
            upper = midpoint;
        } else {
            break;
        }
    }

    (unique_pairs_f * (lower + upper) / 2.0) as u64
}

fn coverage_equation(x: f64, unique_pairs: f64, non_optical_pairs: f64) -> f64 {
    unique_pairs / x - 1.0 + (-non_optical_pairs / x).exp()
}

trait Sink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(File);
struct SamStdout(io::Stdout);

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
        // Shared renderer: htslib `%g` float aux spelling.
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}
impl Sink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}

fn open_output(out: Option<&Path>, fmt: OutFmt, header: &sam::Header) -> io::Result<Box<dyn Sink>> {
    match (out, fmt) {
        (Some(p), OutFmt::Sam) => {
            let mut w = File::create(p)?;
            crate::sam_render::write_header(&mut w, header)?;
            Ok(Box::new(SamFile(w)))
        }
        (Some(p), OutFmt::Bam) => {
            let mut w = bam::io::Writer::new(File::create(p)?);
            w.write_header(header)?;
            Ok(Box::new(BamFile(w)))
        }
        (None, OutFmt::Sam) => {
            let mut w = io::stdout();
            crate::sam_render::write_header(&mut w, header)?;
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
    writeln!(w, "  -s            emit summary counts to stderr")?;
    writeln!(w, "  -S            mark supplementary duplicates (default)")?;
    writeln!(
        w,
        "  -c            clear existing duplicate flags/tags first"
    )?;
    writeln!(w, "  -t            add duplicate-origin do tags")?;
    writeln!(w, "  -d DISTANCE   add duplicate-type dt tags")?;
    writeln!(w, "  --include-fails include QCFAIL reads")?;
    writeln!(w, "  -m t|s        duplicate decision mode")?;
    writeln!(
        w,
        "  -b TAG        include barcode aux tag in duplicate key"
    )?;
    writeln!(w, "  -O sam|bam    output format (default: bam)")?;
    writeln!(w, "  -o FILE       output file (default stdout)")?;
    writeln!(w, "  --no-PG       do not add a @PG line")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn markdup_stats_text_uses_upstream_field_names() {
        let stats = MarkdupStats {
            read: 200,
            written: 180,
            excluded: 0,
            examined: 200,
            paired: 200,
            single: 2,
            duplicate_pair: 20,
            duplicate_single: 1,
            duplicate_pair_optical: 0,
            duplicate_single_optical: 0,
            duplicate_non_primary: 1,
            duplicate_non_primary_optical: 0,
        };
        let mut out = Vec::new();
        write_markdup_stats(&mut out, &stats).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("READ: 200\n"));
        assert!(text.contains("WRITTEN: 180\n"));
        assert!(text.contains("DUPLICATE PAIR: 20\n"));
        assert!(text.contains("DUPLICATE SINGLE: 1\n"));
        assert!(text.contains("DUPLICATE NON PRIMARY: 1\n"));
        assert!(text.contains("DUPLICATE PRIMARY TOTAL: 21\n"));
        assert!(text.contains("DUPLICATE TOTAL: 22\n"));
        assert!(text.contains("ESTIMATED_LIBRARY_SIZE: 466\n"));
    }

    #[test]
    fn estimate_library_size_subtracts_optical_pairs() {
        assert_eq!(estimate_library_size(200, 20, 0), 466);
        assert_eq!(estimate_library_size(200, 20, 2), 510);
        assert_eq!(estimate_library_size(20, 0, 0), 0);
    }
}
