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

use std::collections::{HashMap, HashSet};
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use bstr::BString;
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

#[derive(Clone, Copy)]
struct MarkdupOptions {
    remove_dups: bool,
    emit_stats: bool,
    clear_existing_dups: bool,
    duplicate_origin_tag: bool,
    optical_distance: Option<u32>,
    include_fails: bool,
    barcode_tag: Option<Tag>,
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

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-r" => remove_dups = true,
            "-s" => emit_stats = true,
            "-S" => {
                // Supplementary duplicate propagation is always performed.
            }
            "--include-fails" => include_fails = true,
            "-c" => clear_existing_dups = true,
            "-t" => duplicate_origin_tag = true,
            "-m" | "--mode" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("markdup", format!("missing value for {}", s));
                    return ExitCode::from(1);
                };
                if !matches!(v, "t" | "s") {
                    print_error("markdup", format!("unknown mode '{}'", v));
                    return ExitCode::from(1);
                }
                // Current PE grouping uses the implemented coordinate key for both modes.
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
    let options = MarkdupOptions {
        remove_dups,
        emit_stats,
        clear_existing_dups,
        duplicate_origin_tag,
        optical_distance,
        include_fails,
        barcode_tag,
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
fn mark_duplicates(records: &mut [RecordBuf], options: MarkdupOptions) -> MarkdupStats {
    type BarcodeKey = Option<Vec<u8>>;
    type PosKey = (i32, i64, bool);
    type SeKey = (PosKey, BarcodeKey);
    type PairBarcodeKey = (BarcodeKey, BarcodeKey);
    type PairKey = (PosKey, PosKey, PairBarcodeKey);
    type PairIdx = (usize, usize);
    let mut se_best: HashMap<SeKey, usize> = HashMap::new();
    let mut pair_pending: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut pair_best: HashMap<PairKey, PairIdx> = HashMap::new();
    let mut duplicate_primary_metadata: HashMap<Vec<u8>, DuplicateMetadata> = HashMap::new();
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

    for i in 0..records.len() {
        let flag = records[i].flags().bits() as u32;
        if flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY) != 0 {
            continue;
        }
        if flag & BAM_FQCFAIL != 0 && !options.include_fails {
            stats.excluded += 1;
            continue;
        }
        stats.examined += 1;
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
            stats.paired += 1;
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
                    let first_barcode = barcode_value(&records[first_idx], options.barcode_tag);
                    let barcode = barcode_value(&records[i], options.barcode_tag);

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
                                let origin = records[first_idx]
                                    .name()
                                    .map(|name| name.to_vec())
                                    .unwrap_or_default();
                                let prev_duplicate_type = duplicate_type(
                                    &records[prev_first],
                                    &records[first_idx],
                                    options.optical_distance,
                                );
                                mark_duplicate(
                                    &mut records[prev_first],
                                    Some(&origin),
                                    options.duplicate_origin_tag,
                                    prev_duplicate_type,
                                );
                                mark_duplicate(
                                    &mut records[prev_second],
                                    Some(&origin),
                                    options.duplicate_origin_tag,
                                    prev_duplicate_type,
                                );
                                remember_duplicate_metadata(
                                    &mut duplicate_primary_metadata,
                                    &records[prev_first],
                                    &origin,
                                    prev_duplicate_type,
                                );
                                remember_duplicate_metadata(
                                    &mut duplicate_primary_metadata,
                                    &records[prev_second],
                                    &origin,
                                    prev_duplicate_type,
                                );
                                if prev_duplicate_type.is_some_and(DuplicateType::is_optical) {
                                    stats.duplicate_pair_optical += 2;
                                }
                                pair_best.insert(key, (first_idx, i));
                            } else {
                                let origin = records[prev_first]
                                    .name()
                                    .map(|name| name.to_vec())
                                    .unwrap_or_default();
                                let current_duplicate_type = duplicate_type(
                                    &records[first_idx],
                                    &records[prev_first],
                                    options.optical_distance,
                                );
                                mark_duplicate(
                                    &mut records[first_idx],
                                    Some(&origin),
                                    options.duplicate_origin_tag,
                                    current_duplicate_type,
                                );
                                mark_duplicate(
                                    &mut records[i],
                                    Some(&origin),
                                    options.duplicate_origin_tag,
                                    current_duplicate_type,
                                );
                                remember_duplicate_metadata(
                                    &mut duplicate_primary_metadata,
                                    &records[first_idx],
                                    &origin,
                                    current_duplicate_type,
                                );
                                remember_duplicate_metadata(
                                    &mut duplicate_primary_metadata,
                                    &records[i],
                                    &origin,
                                    current_duplicate_type,
                                );
                                if current_duplicate_type.is_some_and(DuplicateType::is_optical) {
                                    stats.duplicate_pair_optical += 2;
                                }
                            }
                            stats.duplicate_pair += 2;
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
        stats.single += 1;
        let se_key = (me, barcode_value(&records[i], options.barcode_tag));
        match se_best.get(&se_key).copied() {
            Some(idx) => {
                let prev_mapq = records[idx].mapping_quality().map(u8::from).unwrap_or(0);
                if mapq > prev_mapq {
                    let origin = records[i]
                        .name()
                        .map(|name| name.to_vec())
                        .unwrap_or_default();
                    let prev_duplicate_type =
                        duplicate_type(&records[idx], &records[i], options.optical_distance);
                    mark_duplicate(
                        &mut records[idx],
                        Some(&origin),
                        options.duplicate_origin_tag,
                        prev_duplicate_type,
                    );
                    remember_duplicate_metadata(
                        &mut duplicate_primary_metadata,
                        &records[idx],
                        &origin,
                        prev_duplicate_type,
                    );
                    if prev_duplicate_type.is_some_and(DuplicateType::is_optical) {
                        stats.duplicate_single_optical += 1;
                    }
                    se_best.insert(se_key, i);
                } else {
                    let origin = records[idx]
                        .name()
                        .map(|name| name.to_vec())
                        .unwrap_or_default();
                    let current_duplicate_type =
                        duplicate_type(&records[i], &records[idx], options.optical_distance);
                    mark_duplicate(
                        &mut records[i],
                        Some(&origin),
                        options.duplicate_origin_tag,
                        current_duplicate_type,
                    );
                    remember_duplicate_metadata(
                        &mut duplicate_primary_metadata,
                        &records[i],
                        &origin,
                        current_duplicate_type,
                    );
                    if current_duplicate_type.is_some_and(DuplicateType::is_optical) {
                        stats.duplicate_single_optical += 1;
                    }
                }
                stats.duplicate_single += 1;
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
        if duplicate_primary_names.contains(name) {
            let metadata = duplicate_primary_metadata.get(name);
            if mark_duplicate(
                record,
                metadata.map(|metadata| metadata.origin.as_slice()),
                options.duplicate_origin_tag,
                metadata.and_then(|metadata| metadata.duplicate_type),
            ) {
                stats.duplicate_non_primary += 1;
                if metadata
                    .and_then(|metadata| metadata.duplicate_type)
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
    let Some(duplicate_location) = duplicate
        .name()
        .and_then(|name| optical_location(name.as_ref()))
    else {
        return false;
    };
    let Some(original_location) = original
        .name()
        .and_then(|name| optical_location(name.as_ref()))
    else {
        return false;
    };
    duplicate_location.is_within_distance(original_location, distance)
}

#[derive(Clone, Copy)]
struct OpticalLocation {
    tile: i64,
    x: i64,
    y: i64,
}

impl OpticalLocation {
    fn is_within_distance(self, other: Self, distance: u32) -> bool {
        self.tile == other.tile
            && self.x.abs_diff(other.x) <= distance as u64
            && self.y.abs_diff(other.y) <= distance as u64
    }
}

fn optical_location(name: &[u8]) -> Option<OpticalLocation> {
    let mut fields = name.rsplit(|&b| b == b':');
    let y = parse_i64_ascii(fields.next()?)?;
    let x = parse_i64_ascii(fields.next()?)?;
    let tile = parse_i64_ascii(fields.next()?)?;
    Some(OpticalLocation { tile, x, y })
}

fn parse_i64_ascii(bytes: &[u8]) -> Option<i64> {
    std::str::from_utf8(bytes).ok()?.parse().ok()
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
