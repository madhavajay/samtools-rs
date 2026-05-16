//! `samtools ampliconclip` — clip amplicon primers from read ends.
//!
//! Faithful port of `bam_ampliconclip.c`. A BED file of primer sites is
//! loaded per reference (sorted by `right`); each primary mapped read is
//! matched to an overlapping site (`matching_clip_site`) and the
//! overlapping primer bases are soft- or hard-clipped from the read's 5'
//! end (forward / `bam_trim_left`) or 3' end (reverse / `bam_trim_right`),
//! optionally on `--both-ends`. Supports `--strand`, `--original` (`OA`
//! tag), `--keep-tag` (default deletes `NM`/`MD`), `--filter-len`,
//! `--fail-len`, `--unmap-len`, `--clipped`, `--no-excluded`,
//! `--rejects-file`, `--primer-counts`, `--tolerance`, `-f` stats,
//! `-o`/`-O sam|bam`, `-b`, default `@PG`, `--no-PG`. Coordinate-sorted
//! input headers have `@HD SO:` rewritten to `unknown`.
//!
//! SAM output is byte-exact vs the upstream `test/ampliconclip` fixtures
//! (raw-header preserving). **Pending:** CRAM, BGZF fast paths.

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use bstr::BString;
use htslib_rs::core::Position;
use htslib_rs::sam::{
    self,
    alignment::{
        RecordBuf,
        record::{
            MappingQuality,
            cigar::op::{Kind, Op},
            data::field::Tag,
        },
        record_buf::{Cigar, data::field::Value},
    },
};

use crate::diagnostics::{print_error, print_error_errno};

#[derive(Clone, Copy, PartialEq, Eq)]
enum Clip {
    Soft,
    Hard,
}

struct Params {
    add_pg: bool,
    use_strand: bool,
    write_clipped: bool,
    mark_fail: bool,
    both: bool,
    fail_len: i64,
    filter_len: i64,
    unmapped: bool,
    oa_tag: bool,
    del_tag: bool,
    tol: i64,
    unmap_len: i64,
    stats_file: Option<PathBuf>,
    primer_counts_file: Option<PathBuf>,
    rejects_file: Option<PathBuf>,
}

impl Default for Params {
    fn default() -> Self {
        // Mirrors upstream `cl_param_t param = {1,0,0,0,0,-1,-1,0,0,1,5,0,...}`.
        Params {
            add_pg: true,
            use_strand: false,
            write_clipped: false,
            mark_fail: false,
            both: false,
            fail_len: -1,
            filter_len: -1,
            unmapped: false,
            oa_tag: false,
            del_tag: true,
            tol: 5,
            unmap_len: 0,
            stats_file: None,
            primer_counts_file: None,
            rejects_file: None,
        }
    }
}

#[derive(Clone)]
struct BedEntry {
    left: i64,
    right: i64,
    name: String,
    score: String,
    rev: bool,
    num_reads: i64,
}

struct BedList {
    bp: Vec<BedEntry>,
    longest: i64,
}

const BAM_FUNMAP: u32 = 0x4;
const BAM_FREVERSE: u32 = 0x10;
const BAM_FQCFAIL: u32 = 0x200;

fn consumes_query(k: Kind) -> bool {
    matches!(
        k,
        Kind::Match
            | Kind::Insertion
            | Kind::SoftClip
            | Kind::SequenceMatch
            | Kind::SequenceMismatch
    )
}
fn consumes_ref(k: Kind) -> bool {
    matches!(
        k,
        Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
    )
}

fn cigar_ops(r: &RecordBuf) -> Vec<(u32, Kind)> {
    r.cigar()
        .as_ref()
        .iter()
        .map(|op| (op.len() as u32, op.kind()))
        .collect()
}

fn cigar_string(ops: &[(u32, Kind)]) -> String {
    if ops.is_empty() {
        return "*".to_string();
    }
    let mut s = String::new();
    for (len, k) in ops {
        s.push_str(&len.to_string());
        s.push(match k {
            Kind::Match => 'M',
            Kind::Insertion => 'I',
            Kind::Deletion => 'D',
            Kind::Skip => 'N',
            Kind::SoftClip => 'S',
            Kind::HardClip => 'H',
            Kind::Pad => 'P',
            Kind::SequenceMatch => '=',
            Kind::SequenceMismatch => 'X',
        });
    }
    s
}

fn set_cigar(r: &mut RecordBuf, ops: &[(u32, Kind)]) {
    let v: Vec<Op> = ops.iter().map(|(l, k)| Op::new(*k, *l as usize)).collect();
    *r.cigar_mut() = Cigar::from(v);
}

fn pos0(r: &RecordBuf) -> i64 {
    r.alignment_start()
        .map(|p| p.get() as i64 - 1)
        .unwrap_or(-1)
}

/// 0-based exclusive reference end (htslib `bam_endpos`).
fn bam_endpos(r: &RecordBuf) -> i64 {
    let p = pos0(r);
    let rlen: i64 = r
        .cigar()
        .as_ref()
        .iter()
        .filter(|op| consumes_ref(op.kind()))
        .map(|op| op.len() as i64)
        .sum();
    if rlen == 0 { p + 1 } else { p + rlen }
}

/// `active_query_len`: query-consuming bases excluding soft clips.
fn active_query_len(ops: &[(u32, Kind)]) -> i64 {
    ops.iter()
        .filter(|(_, k)| consumes_query(*k) && *k != Kind::SoftClip)
        .map(|(l, _)| *l as i64)
        .sum()
}

fn load_bed(
    path: &Path,
    get_strand: bool,
) -> Result<(HashMap<String, BedList>, Vec<String>), String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("unable to open file {}: {e}", path.display()))?;
    let mut map: HashMap<String, BedList> = HashMap::new();
    let mut order: Vec<String> = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with("track ") || line.starts_with("browser ") {
            continue;
        }
        let cols: Vec<&str> = line.split_whitespace().collect();
        let need = if get_strand { 6 } else { 3 };
        if cols.len() < need {
            return Err(format!(
                "invalid bed file format in line {} of {}",
                lineno + 1,
                path.display()
            ));
        }
        let refn = cols[0].to_string();
        let left: i64 = cols[1]
            .parse()
            .map_err(|_| format!("bad left coord line {}", lineno + 1))?;
        let right: i64 = cols[2]
            .parse()
            .map_err(|_| format!("bad right coord line {}", lineno + 1))?;
        let name = cols.get(3).map(|s| s.to_string()).unwrap_or_default();
        let score = cols.get(4).map(|s| s.to_string()).unwrap_or_default();
        let rev = if get_strand {
            match cols[5] {
                "+" => false,
                "-" => true,
                other => {
                    return Err(format!(
                        "bad strand value in line {}, expecting '+' or '-', found '{}'.",
                        lineno + 1,
                        other
                    ));
                }
            }
        } else {
            false
        };
        let list = map.entry(refn.clone()).or_insert_with(|| {
            order.push(refn.clone());
            BedList {
                bp: Vec::new(),
                longest: 0,
            }
        });
        if right - left > list.longest {
            list.longest = right - left;
        }
        list.bp.push(BedEntry {
            left,
            right,
            name,
            score,
            rev,
            num_reads: 0,
        });
    }
    if map.is_empty() {
        return Err("no usable bed entries".to_string());
    }
    // sort_by_pos: stable sort by `right` (mirrors qsort key bed_entry_sort).
    for list in map.values_mut() {
        list.bp.sort_by_key(|e| e.right);
    }
    Ok((map, order))
}

/// Port of `matching_clip_site`: returns the clip size (reference bases to
/// trim) and bumps the chosen primer's `num_reads`.
fn matching_clip_site(
    sites: &mut BedList,
    pos: i64,
    is_rev: bool,
    use_strand: bool,
    longest: i64,
    tol: i64,
) -> i64 {
    let len = sites.bp.len();
    if len == 0 {
        return 0;
    }
    let mut l = 0usize;
    let mut r = len;
    let mut mid = len / 2;
    let pos_tol = if is_rev {
        if pos > tol { pos - tol } else { 0 }
    } else {
        pos
    };
    while r - l > 1 {
        if sites.bp[mid].right <= pos_tol {
            l = mid;
        } else {
            r = mid;
        }
        mid = (l + r) / 2;
    }

    let mut size = 0i64;
    let mut used_i: isize = -1;
    let mut i = l;
    while i < len {
        let b = &sites.bp[i];
        if use_strand && is_rev != b.rev {
            i += 1;
            continue;
        }
        let (mod_left, mod_right) = if is_rev {
            (b.left, b.right + tol)
        } else {
            (if b.left > tol { b.left - tol } else { 0 }, b.right)
        };
        if pos + longest + tol < mod_right {
            break;
        }
        if pos >= mod_left && pos <= mod_right {
            if is_rev {
                if size < pos - b.left {
                    size = pos - b.left;
                    used_i = i as isize;
                }
            } else if size < b.right - pos {
                size = b.right - pos;
                used_i = i as isize;
            }
        }
        i += 1;
    }
    if used_i >= 0 {
        sites.bp[used_i as usize].num_reads += 1;
    }
    size
}

/// Result of a trim: new cigar, new 0-based pos, query bases physically
/// removed (hard-clip only), and whether the record became empty.
struct Trim {
    cigar: Vec<(u32, Kind)>,
    new_pos0: i64,
    qry_removed: u32,
    empty: bool,
}

fn trim_left(ops: &[(u32, Kind)], pos0: i64, l_qseq: u32, bases: u32, clip: Clip) -> Trim {
    let mut ref_remove = bases;
    let mut qry_removed: u32 = 0;
    let mut hardclip: u32 = 0;
    let mut new_pos = pos0;
    let n = ops.len();
    let mut i = 0;
    while i < n {
        let (len, k) = ops[i];
        if k == Kind::HardClip {
            hardclip += len;
            i += 1;
            continue;
        }
        if consumes_ref(k) {
            if len <= ref_remove {
                ref_remove -= len;
                new_pos += len as i64;
            } else {
                break;
            }
        }
        if consumes_query(k) {
            qry_removed += len;
        }
        i += 1;
    }

    if i < n {
        let (_, k) = ops[i];
        if consumes_ref(k) {
            new_pos += ref_remove as i64;
        }
        if consumes_query(k) {
            qry_removed += ref_remove;
        }
    } else {
        if clip == Clip::Hard {
            return Trim {
                cigar: Vec::new(),
                new_pos0: pos0,
                qry_removed: l_qseq,
                empty: true,
            };
        }
        qry_removed = l_qseq;
    }

    let mut out: Vec<(u32, Kind)> = Vec::new();
    if clip == Clip::Hard && hardclip + qry_removed > 0 {
        out.push((hardclip + qry_removed, Kind::HardClip));
    }
    if clip == Clip::Soft {
        if hardclip > 0 {
            out.push((hardclip, Kind::HardClip));
        }
        if qry_removed > 0 {
            out.push((qry_removed, Kind::SoftClip));
        }
    }
    if i < n && ops[i].0 > ref_remove {
        out.push((ops[i].0 - ref_remove, ops[i].1));
        for op in &ops[i + 1..] {
            out.push(*op);
        }
    }
    let removed = if clip == Clip::Soft { 0 } else { qry_removed };
    Trim {
        cigar: out,
        new_pos0: new_pos,
        qry_removed: removed,
        empty: false,
    }
}

fn trim_right(ops: &[(u32, Kind)], pos0: i64, l_qseq: u32, bases: u32, clip: Clip) -> Trim {
    let mut ref_remove = bases;
    let mut qry_removed: u32 = 0;
    let mut hardclip: u32 = 0;
    let n = ops.len() as i64;
    let mut i = n - 1;
    while i >= 0 {
        let (len, k) = ops[i as usize];
        if k == Kind::HardClip {
            hardclip += len;
            i -= 1;
            continue;
        }
        if consumes_ref(k) {
            if len <= ref_remove {
                ref_remove -= len;
            } else {
                break;
            }
        }
        if consumes_query(k) {
            qry_removed += len;
        }
        i -= 1;
    }

    if i >= 0 {
        let (_, k) = ops[i as usize];
        if consumes_query(k) {
            qry_removed += ref_remove;
        }
    } else if clip == Clip::Hard {
        return Trim {
            cigar: Vec::new(),
            new_pos0: pos0,
            qry_removed: l_qseq,
            empty: true,
        };
    } else {
        qry_removed = l_qseq;
    }

    let mut out: Vec<(u32, Kind)> = Vec::new();
    if i >= 0 {
        for op in &ops[..i as usize] {
            out.push(*op);
        }
        let keep = ops[i as usize].0 - ref_remove;
        if keep > 0 {
            out.push((keep, ops[i as usize].1));
        }
    }
    if clip == Clip::Hard {
        if hardclip + qry_removed > 0 {
            out.push((hardclip + qry_removed, Kind::HardClip));
        }
    } else {
        if qry_removed > 0 {
            out.push((qry_removed, Kind::SoftClip));
        }
        if hardclip > 0 {
            out.push((hardclip, Kind::HardClip));
        }
    }
    let removed = if clip == Clip::Soft { 0 } else { qry_removed };
    Trim {
        cigar: out,
        new_pos0: pos0,
        qry_removed: removed,
        empty: false,
    }
}

fn apply_trim(r: &mut RecordBuf, t: &Trim, left: bool) {
    if t.empty {
        *r.cigar_mut() = Cigar::default();
        r.sequence_mut().as_mut().clear();
        r.quality_scores_mut().as_mut().clear();
        return;
    }
    set_cigar(r, &t.cigar);
    if t.qry_removed > 0 {
        let q = t.qry_removed as usize;
        {
            let seq = r.sequence_mut().as_mut();
            if left {
                seq.drain(..q.min(seq.len()));
            } else {
                let keep = seq.len().saturating_sub(q);
                seq.truncate(keep);
            }
        }
        {
            let qual = r.quality_scores_mut().as_mut();
            if left {
                qual.drain(..q.min(qual.len()));
            } else {
                let keep = qual.len().saturating_sub(q);
                qual.truncate(keep);
            }
        }
    }
    if left {
        *r.alignment_start_mut() = Position::new((t.new_pos0 + 1).max(1) as usize);
    }
}

fn aux_del(r: &mut RecordBuf, tag: Tag) {
    let kept: Vec<_> = r
        .data()
        .iter()
        .filter(|(t, _)| *t != tag)
        .map(|(t, v)| (t, v.clone()))
        .collect();
    *r.data_mut() = kept.into_iter().collect();
}

fn aux_append_str(r: &mut RecordBuf, tag: Tag, value: String) {
    let mut fields: Vec<_> = r
        .data()
        .iter()
        .filter(|(t, _)| *t != tag)
        .map(|(t, v)| (t, v.clone()))
        .collect();
    fields.push((tag, Value::String(BString::from(value))));
    *r.data_mut() = fields.into_iter().collect();
}

fn aux_get_int(r: &RecordBuf, tag: Tag) -> Option<i64> {
    r.data().get(&tag).and_then(|v| v.as_int())
}
fn aux_get_str(r: &RecordBuf, tag: Tag) -> Option<String> {
    match r.data().get(&tag)? {
        Value::String(s) => Some(String::from_utf8_lossy(s).into_owned()),
        _ => None,
    }
}

/// Format `OA:Z:[old]qname,pos+1,strand,CIGAR,MAPQ,NM;`.
fn tag_original_data(orig: &RecordBuf) -> String {
    let mut s = String::new();
    if let Some(old) = aux_get_str(orig, Tag::from([b'O', b'A'])) {
        s.push_str(&old);
    }
    let strand = if rec_flags(orig) & BAM_FREVERSE != 0 {
        '-'
    } else {
        '+'
    };
    let nm = aux_get_int(orig, Tag::from([b'N', b'M']));
    let qname = orig.name().map(|n| n.to_vec()).unwrap_or_default();
    let mapq = orig.mapping_quality().map(u8::from).unwrap_or(255);
    s.push_str(&format!(
        "{},{},{},",
        String::from_utf8_lossy(&qname),
        pos0(orig) + 1,
        strand
    ));
    s.push_str(&cigar_string(&cigar_ops(orig)));
    match nm {
        Some(nm) => s.push_str(&format!(",{},{};", mapq, nm)),
        None => s.push_str(&format!("{},;", mapq)),
    }
    s
}

fn rec_flags(r: &RecordBuf) -> u32 {
    r.flags().bits() as u32
}
fn set_flags(r: &mut RecordBuf, bits: u32) {
    *r.flags_mut() = sam::alignment::record::Flags::from(bits as u16);
}
fn or_flag(r: &mut RecordBuf, bits: u32) {
    let f = rec_flags(r) | bits;
    set_flags(r, f);
}

/// Rewrites `@HD ... SO:coordinate` to `SO:unknown` in raw header text,
/// preserving field order and the rest verbatim.
fn header_so_unknown(raw: &str) -> String {
    let mut out: Vec<String> = Vec::new();
    for line in raw.lines() {
        if line.starts_with("@HD") {
            let fields: Vec<String> = line
                .split('\t')
                .map(|f| {
                    if f == "SO:coordinate" {
                        "SO:unknown".to_string()
                    } else {
                        f.to_string()
                    }
                })
                .collect();
            out.push(fields.join("\t"));
        } else {
            out.push(line.to_string());
        }
    }
    let mut s = out.join("\n");
    s.push('\n');
    s
}

/// Entry point for `samtools ampliconclip`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut bedfile: Option<PathBuf> = None;
    let mut fnout: Option<PathBuf> = None;
    let mut out_bam = true; // default wmode "wb"
    let mut clip = Clip::Soft;
    let mut p = Params::default();
    let mut input: Option<PathBuf> = None;

    let mut it = args.iter().skip(1).peekable();
    while let Some(arg) = it.next() {
        let s = arg.to_str().unwrap_or("");
        if let Some(v) = s.strip_prefix("--output-fmt=") {
            out_bam = !v.eq_ignore_ascii_case("sam");
            continue;
        }
        if let Some(v) = s.strip_prefix("--filter-len=") {
            p.filter_len = v.parse().unwrap_or(-1);
            continue;
        }
        if let Some(v) = s.strip_prefix("--fail-len=") {
            p.fail_len = v.parse().unwrap_or(-1);
            continue;
        }
        match s {
            "-b" => bedfile = it.next().map(PathBuf::from),
            "-o" => fnout = it.next().map(PathBuf::from),
            "-f" => p.stats_file = it.next().map(PathBuf::from),
            "-u" => out_bam = false,
            "-O" | "--output-fmt" => {
                out_bam = !it
                    .next()
                    .and_then(|a| a.to_str())
                    .is_some_and(|v| v.eq_ignore_ascii_case("sam"));
            }
            "-@" | "--threads" => {
                let _ = it.next();
            }
            "--no-PG" => p.add_pg = false,
            "--soft-clip" => clip = Clip::Soft,
            "--hard-clip" => clip = Clip::Hard,
            "--strand" => p.use_strand = true,
            "--clipped" => p.write_clipped = true,
            "--fail" => p.mark_fail = true,
            "--both-ends" => p.both = true,
            "--filter-len" => {
                p.filter_len = it
                    .next()
                    .and_then(|v| v.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(-1)
            }
            "--fail-len" => {
                p.fail_len = it
                    .next()
                    .and_then(|v| v.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(-1)
            }
            "--unmap-len" => {
                p.unmap_len = it
                    .next()
                    .and_then(|v| v.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0)
            }
            "--no-excluded" => p.unmapped = true,
            "--rejects-file" => p.rejects_file = it.next().map(PathBuf::from),
            "--primer-counts" => p.primer_counts_file = it.next().map(PathBuf::from),
            "--original" => p.oa_tag = true,
            "--keep-tag" => p.del_tag = false,
            "--tolerance" => {
                p.tol = it
                    .next()
                    .and_then(|v| v.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(5)
            }
            "--help" => return ExitCode::SUCCESS,
            _ if s.starts_with('-') && s != "-" => {
                print_error(
                    "ampliconclip",
                    format!(
                        "option `{}` is not supported in samtools-rs ampliconclip",
                        s
                    ),
                );
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                }
            }
        }
    }

    let Some(bedfile) = bedfile else {
        print_error("ampliconclip", "a BED file is required (-b)");
        return ExitCode::from(1);
    };
    let Some(input) = input else {
        print_error("ampliconclip", "an input file is required");
        return ExitCode::from(1);
    };
    if p.tol < 0 {
        eprintln!(
            "[ampliconclip] warning: invalid tolerance of {}, resetting tolerance to default of 5.",
            p.tol
        );
        p.tol = 5;
    }

    match run(
        &input,
        fnout.as_deref(),
        out_bam,
        clip,
        &bedfile,
        &mut p,
        args,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("ampliconclip", "ampliconclip failed", &e);
            ExitCode::from(1)
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn run(
    input: &Path,
    output: Option<&Path>,
    out_bam: bool,
    clip: Clip,
    bedfile: &Path,
    p: &mut Params,
    args: &[OsString],
) -> io::Result<()> {
    let (mut bed, ref_order) = load_bed(bedfile, p.use_strand).map_err(io::Error::other)?;
    let raw = crate::header_text::read_raw_header_text(input)?;
    let (header, records) = crate::sam_compat::read_sam_records_tolerant(input)?;

    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect();

    let mut header_text = header_so_unknown(&raw);
    if p.add_pg {
        header_text = crate::pg::add_samtools_pg(&header_text, args).map_err(io::Error::other)?;
    }

    let mut counts = Stats::default();
    let mut kept: Vec<RecordBuf> = Vec::new();
    let mut rejected: Vec<RecordBuf> = Vec::new();

    for rec in records {
        let mut b = rec;
        counts.total += 1;
        let flag = rec_flags(&b);
        let excluded = flag & (BAM_FUNMAP | BAM_FQCFAIL) != 0;
        let ref_name = b
            .reference_sequence_id()
            .and_then(|i| ref_names.get(i))
            .cloned();
        let site_key = ref_name.filter(|n| bed.contains_key(n));

        let mut filter = false;

        if !excluded && let Some(key) = site_key {
            let oat = if p.oa_tag {
                Some(tag_original_data(&b))
            } else {
                None
            };
            let longest = bed[&key].longest;
            let tol = p.tol;
            let mut been_clipped = false;

            if !p.both {
                let rev = flag & BAM_FREVERSE != 0;
                let (pos, is_rev) = if rev {
                    (bam_endpos(&b), true)
                } else {
                    (pos0(&b), false)
                };
                let sites = bed.get_mut(&key).unwrap();
                let psize = matching_clip_site(sites, pos, is_rev, p.use_strand, longest, tol);
                if psize > 0 {
                    let ops = cigar_ops(&b);
                    let l_qseq = b.sequence().len() as u32;
                    let t = if is_rev {
                        counts.rev += 1;
                        trim_right(&ops, pos0(&b), l_qseq, psize as u32, clip)
                    } else {
                        counts.fwd += 1;
                        trim_left(&ops, pos0(&b), l_qseq, psize as u32, clip)
                    };
                    apply_trim(&mut b, &t, !is_rev);
                    finalize_clip(&mut b, p, oat.as_deref());
                    been_clipped = true;
                } else {
                    if p.mark_fail {
                        or_flag(&mut b, BAM_FQCFAIL);
                    }
                    counts.not_clipped += 1;
                }
            } else {
                let mut left_done = false;
                let mut right_done = false;
                {
                    let pos = pos0(&b);
                    let sites = bed.get_mut(&key).unwrap();
                    let psize = matching_clip_site(sites, pos, false, p.use_strand, longest, tol);
                    if psize > 0 {
                        let ops = cigar_ops(&b);
                        let l_qseq = b.sequence().len() as u32;
                        let t = trim_left(&ops, pos0(&b), l_qseq, psize as u32, clip);
                        apply_trim(&mut b, &t, true);
                        counts.fwd += 1;
                        left_done = true;
                        been_clipped = true;
                    }
                }
                {
                    let pos = bam_endpos(&b);
                    let sites = bed.get_mut(&key).unwrap();
                    let psize = matching_clip_site(sites, pos, true, p.use_strand, longest, tol);
                    if psize > 0 {
                        let ops = cigar_ops(&b);
                        let l_qseq = b.sequence().len() as u32;
                        let t = trim_right(&ops, pos0(&b), l_qseq, psize as u32, clip);
                        apply_trim(&mut b, &t, false);
                        counts.rev += 1;
                        right_done = true;
                        been_clipped = true;
                    }
                }
                if left_done || right_done {
                    finalize_clip(&mut b, p, oat.as_deref());
                }
                if left_done && right_done {
                    counts.both += 1;
                } else if !left_done && !right_done {
                    if p.mark_fail {
                        or_flag(&mut b, BAM_FQCFAIL);
                    }
                    counts.not_clipped += 1;
                }
            }

            if p.fail_len >= 0 || p.filter_len >= 0 || p.unmap_len >= 0 {
                let aql = active_query_len(&cigar_ops(&b));
                if p.fail_len >= 0 && aql <= p.fail_len {
                    or_flag(&mut b, BAM_FQCFAIL);
                }
                if p.filter_len >= 0 && aql <= p.filter_len {
                    filter = true;
                }
                if p.unmap_len >= 0 && aql <= p.unmap_len {
                    or_flag(&mut b, BAM_FUNMAP);
                    *b.mapping_quality_mut() = MappingQuality::new(0);
                    *b.cigar_mut() = Cigar::default();
                }
            }

            if rec_flags(&b) & BAM_FQCFAIL != 0 {
                counts.failed += 1;
            }
            if p.write_clipped && !been_clipped {
                filter = true;
            }
        } else {
            counts.excluded += 1;
            if p.unmapped {
                filter = true;
            }
        }

        if !filter {
            kept.push(b);
            counts.written += 1;
        } else {
            if p.rejects_file.is_some() {
                rejected.push(b);
            }
            counts.filtered += 1;
        }
    }

    write_sam_or_bam(output, out_bam, &header, &header_text, &kept)?;
    if let Some(rf) = &p.rejects_file {
        write_sam_or_bam(Some(rf), out_bam, &header, &header_text, &rejected)?;
    }

    let cmd = crate::pg::stringify_argv(args);
    let stats_text = format!(
        "COMMAND: {cmd}\nTOTAL READS: {}\nTOTAL CLIPPED: {}\nFORWARD CLIPPED: {}\nREVERSE CLIPPED: {}\nBOTH CLIPPED: {}\nNOT CLIPPED: {}\nEXCLUDED: {}\nFILTERED: {}\nFAILED: {}\nWRITTEN: {}\n",
        counts.total,
        counts.fwd + counts.rev,
        counts.fwd,
        counts.rev,
        counts.both,
        counts.not_clipped,
        counts.excluded,
        counts.filtered,
        counts.failed,
        counts.written,
    );
    if let Some(sf) = &p.stats_file {
        std::fs::write(sf, stats_text)?;
    } else {
        eprint!("{stats_text}");
    }

    if let Some(pcf) = &p.primer_counts_file {
        let mut w = BufWriter::new(File::create(pcf)?);
        writeln!(w, "#CHR\tLEFT\tRIGHT\tNAME\tSCORE\tSTRAND\tNUM_CLIPPED")?;
        for refn in &ref_order {
            if let Some(list) = bed.get(refn) {
                for e in &list.bp {
                    let strand = if p.use_strand {
                        if e.rev { "-" } else { "+" }
                    } else {
                        "."
                    };
                    writeln!(
                        w,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        refn, e.left, e.right, e.name, e.score, strand, e.num_reads
                    )?;
                }
            }
        }
        w.flush()?;
    }

    Ok(())
}

/// Post-clip housekeeping: add `OA`, delete `NM`/`MD` unless `--keep-tag`.
fn finalize_clip(b: &mut RecordBuf, p: &Params, oat: Option<&str>) {
    if let Some(oa) = oat {
        aux_append_str(b, Tag::from([b'O', b'A']), oa.to_string());
    }
    if p.del_tag {
        aux_del(b, Tag::from([b'N', b'M']));
        aux_del(b, Tag::from([b'M', b'D']));
    }
}

fn write_sam_or_bam(
    output: Option<&Path>,
    out_bam: bool,
    header: &sam::Header,
    header_text: &str,
    records: &[RecordBuf],
) -> io::Result<()> {
    if out_bam {
        use sam::alignment::io::Write as _;
        let mut hdr = header.clone();
        if let Ok(parsed) = header_text.parse::<sam::Header>() {
            hdr = parsed;
        }
        match output {
            Some(p) => {
                let mut w = htslib_rs::bam::io::Writer::new(File::create(p)?);
                w.write_header(&hdr)?;
                for r in records {
                    w.write_alignment_record(&hdr, r)?;
                }
            }
            None => {
                let mut w = htslib_rs::bam::io::Writer::new(io::stdout().lock());
                w.write_header(&hdr)?;
                for r in records {
                    w.write_alignment_record(&hdr, r)?;
                }
            }
        }
        return Ok(());
    }
    let mut out: Box<dyn Write> = match output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };
    out.write_all(header_text.as_bytes())?;
    for r in records {
        crate::sam_render::write_record(&mut out, header, r)?;
    }
    out.flush()?;
    Ok(())
}

#[derive(Default)]
struct Stats {
    total: i64,
    fwd: i64,
    rev: i64,
    both: i64,
    not_clipped: i64,
    excluded: i64,
    filtered: i64,
    failed: i64,
    written: i64,
}
