//! `samtools ampliconstats` — amplicon-level coverage / usage statistics.
//!
//! Faithful port of `amplicon_stats.c`. Loads a primer BED (per reference,
//! file order, with strand), collapses LEFT/RIGHT primer pairs into
//! amplicons (`bed2amplicon`), builds a position→amplicon lookup
//! (`initialise_amp_pos_lookup`, ±`pos-margin`), accumulates per-read
//! stats with read-pair overlap removal (`accumulate_stats`), and emits
//! the multi-section `SS`/`AMPLICON`/`F*`/`C*` report (`dump_stats`).
//! SAM/BAM input; one row per file plus a `COMBINED` mean/stddev block.
//!
//! Byte-exact (modulo the harness-stripped `Samtools version`/`Command
//! line` lines) vs upstream `test/ampliconstats/stats{,_mixed,_partial}`.
//! **Pending:** `--tcoord-bin` aggregation, CRAM, `--use-sample-name`.

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::sam::alignment::{RecordBuf, record::cigar::op::Kind};

use crate::diagnostics::{print_error, print_error_errno};

const MAX_DEPTH: usize = 5;

struct Args {
    flag_require: u32,
    flag_filter: u32,
    max_delta: i64,
    min_depth: [i64; MAX_DEPTH],
    max_amp: i64,
    max_amp_len: i64,
    tlen_adj: i64,
    tcoord_min_count: i64,
    tcoord_bin: i64,
    depth_bin: f64,
    multi_ref: bool,
    argv: String,
    out: Option<PathBuf>,
}

impl Default for Args {
    fn default() -> Self {
        Args {
            flag_require: 0,
            flag_filter: 0x10B04,
            max_delta: 30,
            min_depth: [1, 0, 0, 0, 0],
            max_amp: 1000 + 1,
            max_amp_len: 1000 + 1,
            tlen_adj: 0,
            tcoord_min_count: 10,
            tcoord_bin: 1,
            depth_bin: 0.01,
            multi_ref: true,
            argv: String::new(),
            out: None,
        }
    }
}

#[derive(Clone)]
struct BedEntry {
    left: i64,
    right: i64,
    rev: bool,
}

const BAM_FPAIRED: u32 = 0x1;
const BAM_FUNMAP: u32 = 0x4;
const BAM_FMUNMAP: u32 = 0x8;
const BAM_FREVERSE: u32 = 0x10;
const BAM_FSECONDARY: u32 = 0x100;
const BAM_FSUPPLEMENTARY: u32 = 0x800;

/// One amplicon: alt-primer LEFT/RIGHT inner+outer coordinates.
#[derive(Default, Clone)]
struct Amplicon {
    left: Vec<i64>,
    right: Vec<i64>,
    max_left: i64,
    min_right: i64,
    min_left: i64,
    max_right: i64,
}

/// Per-(file|combined) accumulators for one reference.
struct Stats {
    nseq: i64,
    nfiltered: i64,
    nfailprimer: i64,
    max_amp_len: i64,
    nreads: Vec<i64>,
    nreads2: Vec<i64>,
    nfull_reads: Vec<f64>,
    nrperc: Vec<f64>,
    nrperc2: Vec<f64>,
    nbases: Vec<i64>,
    nbases2: Vec<i64>,
    coverage: Vec<i64>,
    covered_perc: Vec<[f64; MAX_DEPTH]>,
    covered_perc2: Vec<[f64; MAX_DEPTH]>,
    amp_dist: Vec<[i64; 3]>,
    depth_all: Vec<i64>,
    depth_valid: Vec<i64>,
    // tcoord[anum+1]: map (start|end<<32) -> (freq | status<<32)
    tcoord: Vec<HashMap<u64, u64>>,
    qend: HashMap<Vec<u8>, u64>,
}

impl Stats {
    fn new(len: i64, max_amp: i64, max_amp_len: i64) -> Self {
        let na = max_amp as usize;
        Stats {
            nseq: 0,
            nfiltered: 0,
            nfailprimer: 0,
            max_amp_len,
            nreads: vec![0; na],
            nreads2: vec![0; na],
            nfull_reads: vec![0.0; na],
            nrperc: vec![0.0; na],
            nrperc2: vec![0.0; na],
            nbases: vec![0; na],
            nbases2: vec![0; na],
            coverage: vec![0; na * max_amp_len as usize],
            covered_perc: vec![[0.0; MAX_DEPTH]; na],
            covered_perc2: vec![[0.0; MAX_DEPTH]; na],
            amp_dist: vec![[0; 3]; na],
            depth_all: vec![0; len.max(1) as usize],
            depth_valid: vec![0; len.max(1) as usize],
            tcoord: (0..na + 1).map(|_| HashMap::new()).collect(),
            qend: HashMap::new(),
        }
    }
    fn reset(&mut self) {
        self.nseq = 0;
        self.nfiltered = 0;
        self.nfailprimer = 0;
        self.nreads.iter_mut().for_each(|x| *x = 0);
        self.nfull_reads.iter_mut().for_each(|x| *x = 0.0);
        self.nbases.iter_mut().for_each(|x| *x = 0);
        self.coverage.iter_mut().for_each(|x| *x = 0);
        self.amp_dist.iter_mut().for_each(|x| *x = [0; 3]);
        self.depth_all.iter_mut().for_each(|x| *x = 0);
        self.depth_valid.iter_mut().for_each(|x| *x = 0);
        self.tcoord.iter_mut().for_each(|m| m.clear());
        self.qend.clear();
    }
}

struct RefAmp {
    namp: usize,
    len: i64,
    sites: Vec<BedEntry>,
    amp: Vec<Amplicon>,
    refname: String,
    first_amp: usize,
    lstats: Stats,
    gstats: Stats,
}

type BedData = (Vec<String>, HashMap<String, Vec<BedEntry>>);

fn load_bed(path: &Path) -> Result<BedData, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("unable to open file {}: {e}", path.display()))?;
    let mut order: Vec<String> = Vec::new();
    let mut map: HashMap<String, Vec<BedEntry>> = HashMap::new();
    for (lineno, line) in text.lines().enumerate() {
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track ")
            || line.starts_with("browser ")
        {
            continue;
        }
        let c: Vec<&str> = line.split_whitespace().collect();
        if c.len() < 6 {
            return Err(format!(
                "invalid bed file format in line {} of {}",
                lineno + 1,
                path.display()
            ));
        }
        let rev = match c[5] {
            "+" => false,
            "-" => true,
            o => {
                return Err(format!(
                    "bad strand value in line {}, expecting '+' or '-', found '{}'.",
                    lineno + 1,
                    o
                ));
            }
        };
        let refn = c[0].to_string();
        let e = BedEntry {
            left: c[1].parse().map_err(|_| "bad left".to_string())?,
            right: c[2].parse().map_err(|_| "bad right".to_string())?,
            rev,
        };
        map.entry(refn.clone())
            .or_insert_with(|| {
                order.push(refn.clone());
                Vec::new()
            })
            .push(e);
    }
    if map.is_empty() {
        return Err("no usable bed entries".to_string());
    }
    Ok((order, map))
}

/// Count amplicons: a `+` primer following a `-` primer starts a new one.
fn count_amplicon(sites: &[BedEntry]) -> usize {
    let mut namp = 0;
    let mut last_rev = false;
    for s in sites {
        if !s.rev && last_rev {
            namp += 1;
        }
        last_rev = s.rev;
    }
    namp + 1
}

/// Collapse LEFT end / RIGHT start primers into amplicons, emitting the
/// `AMPLICON` block when `do_title`.
fn bed2amplicon(
    args: &Args,
    sites: &[BedEntry],
    amp: &mut Vec<Amplicon>,
    do_title: bool,
    refname: &str,
    first_amp: usize,
    out: &mut dyn Write,
) -> Result<usize, String> {
    amp.clear();
    amp.push(Amplicon {
        max_left: 0,
        min_right: i64::MAX,
        min_left: i64::MAX,
        max_right: 0,
        ..Default::default()
    });
    if do_title {
        writeln!(out, "# Amplicon locations from BED file.").ok();
        writeln!(
            out,
            "# LEFT/RIGHT are <start>-<end> format and comma-separated for alt-primers."
        )
        .ok();
        if args.multi_ref {
            writeln!(out, "#\n# AMPLICON\tREF\tNUMBER\tLEFT\tRIGHT").ok();
        } else {
            writeln!(out, "#\n# AMPLICON\tNUMBER\tLEFT\tRIGHT").ok();
        }
    }
    let mut j = 0usize;
    let mut last_rev = false;
    let mut line = String::new();
    for (i, s) in sites.iter().enumerate() {
        if i == 0 && s.rev {
            return Err("BED file should start with the + strand primer".into());
        }
        if !s.rev && last_rev {
            j += 1;
            if j as i64 >= args.max_amp {
                return Err(format!("too many amplicons ({})", j));
            }
            amp.push(Amplicon {
                max_left: 0,
                min_right: i64::MAX,
                min_left: i64::MAX,
                max_right: 0,
                ..Default::default()
            });
        }
        if !s.rev {
            if i == 0 || last_rev {
                if j > 0 {
                    line.push('\n');
                }
                if args.multi_ref {
                    line.push_str(&format!("AMPLICON\t{}\t{}", refname, j + 1 + first_amp));
                } else {
                    line.push_str(&format!("AMPLICON\t{}", j + 1));
                }
            }
            let a = &mut amp[j];
            a.left.push(s.right);
            if a.max_left < s.right + 1 {
                a.max_left = s.right + 1;
            }
            if a.min_left > s.right + 1 {
                a.min_left = s.right + 1;
            }
            line.push(if a.left.len() > 1 { ',' } else { '\t' });
            line.push_str(&format!("{}-{}", s.left + 1, s.right));
        } else {
            let a = &mut amp[j];
            a.right.push(s.left);
            if a.min_right > s.left - 1 {
                a.min_right = s.left - 1;
            }
            if a.max_right < s.left - 1 {
                a.max_right = s.left - 1;
                if a.max_right - a.min_left + 1 >= args.max_amp_len {
                    return Err("amplicon longer than max_amp_len".into());
                }
            }
            line.push(if a.right.len() > 1 { ',' } else { '\t' });
            line.push_str(&format!("{}-{}", s.left + 1, s.right));
        }
        last_rev = s.rev;
    }
    if !last_rev {
        return Err("bed file does not end on a reverse strand primer.".into());
    }
    let namp = j + 1;
    if !line.is_empty() {
        line.push('\n');
    }
    if do_title || !line.is_empty() {
        out.write_all(line.as_bytes()).ok();
    }
    Ok(namp)
}

fn op_consumes_ref(k: Kind) -> bool {
    matches!(
        k,
        Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
    )
}

fn pos0(r: &RecordBuf) -> i64 {
    r.alignment_start()
        .map(|p| p.get() as i64 - 1)
        .unwrap_or(-1)
}
fn bam_endpos(r: &RecordBuf) -> i64 {
    let p = pos0(r);
    let rlen: i64 = r
        .cigar()
        .as_ref()
        .iter()
        .filter(|op| op_consumes_ref(op.kind()))
        .map(|op| op.len() as i64)
        .sum();
    if rlen == 0 { p + 1 } else { p + rlen }
}
fn rec_flags(r: &RecordBuf) -> u32 {
    r.flags().bits() as u32
}
fn isize_of(r: &RecordBuf) -> i64 {
    r.template_length() as i64
}

/// Position→amplicon lookup (±`max_delta`); -1 = no amplicon.
fn build_pos_lookup(args: &Args, ra: &RefAmp) -> (Vec<i64>, Vec<i64>) {
    let max_len = ra.len.max(0) as usize;
    let mut p2s = vec![-1i64; max_len + 1];
    let mut p2e = vec![-1i64; max_len + 1];
    for (i, a) in ra.amp.iter().enumerate().take(ra.namp) {
        for &l in &a.left {
            let mut p = l - args.max_delta;
            while p <= l + args.max_delta {
                if p >= 1 && p <= ra.len {
                    p2s[(p - 1) as usize] = i as i64;
                }
                p += 1;
            }
        }
        for &r in &a.right {
            let mut p = r - args.max_delta;
            while p <= r + args.max_delta {
                if p >= 1 && p <= ra.len {
                    p2e[(p - 1) as usize] = i as i64;
                }
                p += 1;
            }
        }
    }
    (p2s, p2e)
}

#[allow(clippy::too_many_arguments)]
fn accumulate_stats(args: &Args, ra: &mut RefAmp, b: &RecordBuf, p2s: &[i64], p2e: &[i64]) {
    let stats = &mut ra.lstats;
    let len = ra.len;
    stats.nseq += 1;
    let flag = rec_flags(b);
    if (flag & args.flag_require) != args.flag_require || (flag & args.flag_filter) != 0 {
        stats.nfiltered += 1;
        return;
    }
    let start = pos0(b);
    let mut mstart = start;
    let mut end = bam_endpos(b);
    let mut prev_start = 0i64;
    let mut prev_end = 0i64;
    if flag & BAM_FPAIRED != 0 && flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY) == 0 {
        let name = b.name().map(|n| n.to_vec()).unwrap_or_default();
        if let Some(v) = stats.qend.remove(&name) {
            prev_start = (v & 0xffff_ffff) as i64;
            prev_end = (v >> 32) as i64;
            mstart = mstart.max(prev_end);
        } else {
            stats
                .qend
                .insert(name, (start as u64 & 0xffff_ffff) | ((end as u64) << 32));
        }
    }
    let mut i = mstart;
    while i < end && i < len {
        if i >= 0 {
            stats.depth_all[i as usize] += 1;
        }
        i += 1;
    }

    let anum = if flag & BAM_FREVERSE != 0 || flag & BAM_FPAIRED == 0 {
        if end > 0 && end - 1 < len {
            p2e[(end - 1) as usize]
        } else {
            -1
        }
    } else if start >= 0 && start < len {
        p2s[start as usize]
    } else {
        -1
    };

    if end == start && (args.flag_filter & BAM_FUNMAP) != 0 {
        stats.nfiltered += 1;
        return;
    }
    if anum == -1 {
        stats.nfailprimer += 1;
    }
    let mut start = start;
    if anum >= 0 {
        let ai = anum as usize;
        let amp = &ra.amp[ai];
        let c = end.min(amp.min_right + 1) - start.max(amp.max_left);
        if c > 0 {
            stats.nreads[ai] += 1;
            stats.nbases[ai] += c;
            if start < 0 {
                start = 0;
            }
            if end > len {
                end = len;
            }
            let ostart = start.max(amp.min_left - 1);
            let oend = end.min(amp.max_right);
            let offset = amp.min_left - 1;
            let mut k = ostart;
            while k < oend {
                let idx = ai * stats.max_amp_len as usize + (k - offset) as usize;
                stats.coverage[idx] += 1;
                k += 1;
            }
        } else {
            stats.nfailprimer += 1;
        }
    }

    let mut t_end;
    let mut oth_anum: i64 = -1;
    if flag & BAM_FPAIRED != 0 {
        let isz = isize_of(b);
        t_end = if flag & BAM_FREVERSE != 0 { end } else { start } + isz;
        t_end += if isz > 0 {
            -args.tlen_adj
        } else {
            args.tlen_adj
        };
        if t_end > 0 && t_end < len && isz != 0 {
            oth_anum = if flag & BAM_FREVERSE != 0 {
                p2s[t_end as usize]
            } else {
                p2e[t_end as usize]
            };
        }
    } else {
        oth_anum = if start >= 0 && start < len {
            p2s[start as usize]
        } else {
            -1
        };
        t_end = end;
    }

    let mut astatus = 2;
    if anum != -1 && oth_anum != -1 {
        astatus = if oth_anum == anum { 0 } else { 1 };
        if start <= t_end {
            ra.lstats.amp_dist[anum as usize][astatus] += 1;
        }
    } else if anum >= 0 {
        astatus = 2;
        ra.lstats.amp_dist[anum as usize][2] += 1;
    }

    if astatus == 0 && flag & (BAM_FUNMAP | BAM_FMUNMAP) == 0 {
        let half = if flag & BAM_FPAIRED != 0 { 0.5 } else { 1.0 };
        if prev_end != 0 && mstart > prev_end {
            let mut k = prev_start;
            while k < prev_end {
                if k >= 0 && (k as usize) < ra.lstats.depth_valid.len() {
                    ra.lstats.depth_valid[k as usize] -= 1;
                }
                k += 1;
            }
            ra.lstats.nfull_reads[anum as usize] -= half;
        } else {
            let mut k = mstart;
            while k < end {
                if k >= 0 && (k as usize) < ra.lstats.depth_valid.len() {
                    ra.lstats.depth_valid[k as usize] += 1;
                }
                k += 1;
            }
            ra.lstats.nfull_reads[anum as usize] += half;
        }
    }

    if flag & BAM_FPAIRED != 0 && isize_of(b) <= 0 {
        return;
    }
    let start2 = pos0(b);
    let t_end2 = if flag & BAM_FPAIRED != 0 {
        start2 + isize_of(b) - 1
    } else {
        bam_endpos(b)
    };
    let tk = ((start2 + 1).min(u32::MAX as i64) as u64)
        | (((t_end2 + 1).min(u32::MAX as i64) as u64) << 32);
    let m = &mut ra.lstats.tcoord[(anum + 1) as usize];
    let freq = (m.get(&tk).map(|v| v & 0xffff_ffff).unwrap_or(0)) + 1;
    m.insert(tk, (freq & 0xffff_ffff) | ((astatus as u64) << 32));
}

fn append_lstats(l: &Stats, g: &mut Stats, namp: usize, all_nseq: i64) {
    g.nseq += l.nseq;
    g.nfiltered += l.nfiltered;
    g.nfailprimer += l.nfailprimer;
    for a in 0..=namp {
        // a here is anum+1 index space: a==0 is the unmatched bucket.
        for (&k, &v) in &l.tcoord[a] {
            if v & 0xffff_ffff == 0 {
                continue;
            }
            let cur = g.tcoord[a].get(&k).map(|x| x & 0xffff_ffff).unwrap_or(0);
            let add = l.tcoord[a][&k] & 0xffff_ffff;
            let status = v & !0xffff_ffff;
            g.tcoord[a].insert(k, ((cur + add) & 0xffff_ffff) | status);
        }
        if a == 0 {
            continue;
        }
        let ai = a - 1;
        g.nreads[ai] += l.nreads[ai];
        g.nreads2[ai] += l.nreads[ai] * l.nreads[ai];
        g.nfull_reads[ai] += l.nfull_reads[ai];
        let nrperc = if all_nseq != 0 {
            100.0 * l.nreads[ai] as f64 / all_nseq as f64
        } else {
            0.0
        };
        g.nrperc[ai] += nrperc;
        g.nrperc2[ai] += nrperc * nrperc;
        g.nbases[ai] += l.nbases[ai];
        g.nbases2[ai] += l.nbases[ai] * l.nbases[ai];
        for d in 0..MAX_DEPTH {
            g.covered_perc[ai][d] += l.covered_perc[ai][d];
            g.covered_perc2[ai][d] += l.covered_perc[ai][d] * l.covered_perc[ai][d];
        }
        for d in 0..3 {
            g.amp_dist[ai][d] += l.amp_dist[ai][d];
        }
    }
    for k in 0..l.depth_all.len().min(g.depth_all.len()) {
        g.depth_valid[k] += l.depth_valid[k];
        g.depth_all[k] += l.depth_all[k];
    }
}

#[allow(clippy::too_many_arguments)]
fn dump_stats(
    args: &Args,
    type_ch: char,
    name: &str,
    nfile: i64,
    refs: &mut [RefAmp],
    nref: usize,
    local: bool,
    out: &mut dyn Write,
) {
    let combined = type_ch == 'C';
    fn pick(ra: &RefAmp, local: bool) -> &Stats {
        if local { &ra.lstats } else { &ra.gstats }
    }
    macro_rules! S {
        ($ra:expr) => {
            pick($ra, local)
        };
    }

    writeln!(out, "# Summary stats.").ok();
    writeln!(
        out,
        "# Use 'grep ^{}SS | cut -f 2-' to extract this part.",
        type_ch
    )
    .ok();
    for ra in refs.iter().take(nref) {
        if ra.sites.is_empty() {
            continue;
        }
        let stats = S!(ra);
        let nmatch = stats.nseq - stats.nfiltered - stats.nfailprimer;
        let name_ref = if args.multi_ref {
            format!("{}\t{}", name, ra.refname)
        } else {
            name.to_string()
        };
        writeln!(
            out,
            "{}SS\t{}\traw total sequences:\t{}",
            type_ch, name_ref, stats.nseq
        )
        .ok();
        writeln!(
            out,
            "{}SS\t{}\tfiltered sequences:\t{}",
            type_ch, name_ref, stats.nfiltered
        )
        .ok();
        writeln!(
            out,
            "{}SS\t{}\tfailed primer match:\t{}",
            type_ch, name_ref, stats.nfailprimer
        )
        .ok();
        writeln!(
            out,
            "{}SS\t{}\tmatching sequences:\t{}",
            type_ch, name_ref, nmatch
        )
        .ok();
        let mut d = 0;
        loop {
            let mut start = 0i64;
            let mut covered = 0i64;
            let mut total = 0i64;
            for i in 0..ra.namp {
                let amp = &ra.amp[i];
                let offset = amp.min_left - 1;
                let mut j = start.max(amp.max_left - 1);
                let jend = start.max(amp.min_right);
                while j < jend {
                    let idx = i * stats.max_amp_len as usize + (j - offset) as usize;
                    if stats.coverage[idx] >= args.min_depth[d] {
                        covered += 1;
                    }
                    total += 1;
                    j += 1;
                }
                start = start.max(amp.min_right);
            }
            writeln!(
                out,
                "{}SS\t{}\tconsensus depth count < {} and >= {}:\t{}\t{}",
                type_ch,
                name_ref,
                args.min_depth[d],
                args.min_depth[d],
                total - covered,
                covered
            )
            .ok();
            d += 1;
            if d >= MAX_DEPTH || args.min_depth[d] == 0 {
                break;
            }
        }
    }

    // READS
    writeln!(out, "# Absolute matching read counts per amplicon.").ok();
    writeln!(
        out,
        "# Use 'grep ^{}READS | cut -f 2-' to extract this part.",
        type_ch
    )
    .ok();
    let mut buf = format!("{}READS\t{}", type_ch, name);
    for ra in refs.iter().take(nref) {
        if ra.sites.is_empty() {
            continue;
        }
        let s = S!(ra);
        for i in 0..ra.namp {
            buf.push_str(&format!("\t{}", s.nreads[i]));
        }
    }
    writeln!(out, "{}", buf).ok();

    let mut buf = format!("{}VDEPTH\t{}", type_ch, name);
    for ra in refs.iter().take(nref) {
        if ra.sites.is_empty() {
            continue;
        }
        let s = S!(ra);
        for i in 0..ra.namp {
            buf.push_str(&format!("\t{}", s.nfull_reads[i] as i64));
        }
    }
    writeln!(out, "{}", buf).ok();

    if combined {
        let mut m = String::from("CREADS\tMEAN");
        for ra in refs.iter().take(nref) {
            if ra.sites.is_empty() {
                continue;
            }
            let s = S!(ra);
            for i in 0..ra.namp {
                m.push_str(&format!("\t{:.1}", s.nreads[i] as f64 / nfile as f64));
            }
        }
        writeln!(out, "{}", m).ok();
        let mut m = String::from("CREADS\tSTDDEV");
        for ra in refs.iter().take(nref) {
            if ra.sites.is_empty() {
                continue;
            }
            let s = S!(ra);
            for i in 0..ra.namp {
                let n1 = s.nreads[i] as f64;
                let v = if nfile > 1 && s.nreads2[i] > 0 {
                    (s.nreads2[i] as f64 / nfile as f64 - (n1 / nfile as f64).powi(2)).max(0.0)
                } else {
                    0.0
                };
                m.push_str(&format!("\t{:.1}", v.sqrt()));
            }
        }
        writeln!(out, "{}", m).ok();
    }

    // RPERC
    writeln!(out, "# Read percentage of distribution between amplicons.").ok();
    writeln!(
        out,
        "# Use 'grep ^{}RPERC | cut -f 2-' to extract this part.",
        type_ch
    )
    .ok();
    let mut all_nseq = 0i64;
    for ra in refs.iter().take(nref) {
        if ra.sites.is_empty() {
            continue;
        }
        let s = S!(ra);
        all_nseq += s.nseq - s.nfiltered - s.nfailprimer;
    }
    let mut buf = format!("{}RPERC\t{}", type_ch, name);
    for ra in refs.iter().take(nref) {
        if ra.sites.is_empty() {
            continue;
        }
        let s = S!(ra);
        for i in 0..ra.namp {
            if combined {
                buf.push_str(&format!("\t{:.3}", s.nrperc[i] / nfile as f64));
            } else {
                let v = if all_nseq != 0 {
                    100.0 * s.nreads[i] as f64 / all_nseq as f64
                } else {
                    0.0
                };
                buf.push_str(&format!("\t{:.3}", v));
            }
        }
    }
    writeln!(out, "{}", buf).ok();
    if combined {
        let mut m = String::from("CRPERC\tMEAN");
        for ra in refs.iter().take(nref) {
            if ra.sites.is_empty() {
                continue;
            }
            let s = S!(ra);
            for i in 0..ra.namp {
                m.push_str(&format!("\t{:.3}", s.nrperc[i] / nfile as f64));
            }
        }
        writeln!(out, "{}", m).ok();
        let mut m = String::from("CRPERC\tSTDDEV");
        for ra in refs.iter().take(nref) {
            if ra.sites.is_empty() {
                continue;
            }
            let s = S!(ra);
            for i in 0..ra.namp {
                let n1 = s.nrperc[i];
                let v = s.nrperc2[i] / nfile as f64 - (n1 / nfile as f64).powi(2);
                m.push_str(&format!("\t{:.3}", if v > 0.0 { v.sqrt() } else { 0.0 }));
            }
        }
        writeln!(out, "{}", m).ok();
    }

    // DEPTH
    writeln!(out, "# Read depth per amplicon.").ok();
    writeln!(
        out,
        "# Use 'grep ^{}DEPTH | cut -f 2-' to extract this part.",
        type_ch
    )
    .ok();
    let mut buf = format!("{}DEPTH\t{}", type_ch, name);
    for ra in refs.iter().take(nref) {
        if ra.sites.is_empty() {
            continue;
        }
        let s = S!(ra);
        let nseq = s.nseq - s.nfiltered - s.nfailprimer;
        for i in 0..ra.namp {
            let amp = &ra.amp[i];
            let alen = amp.min_right - amp.max_left + 1;
            let v = if nseq != 0 {
                s.nbases[i] as f64 / alen as f64
            } else {
                0.0
            };
            buf.push_str(&format!("\t{:.1}", v));
        }
    }
    writeln!(out, "{}", buf).ok();
    if combined {
        let mut m = String::from("CDEPTH\tMEAN");
        for ra in refs.iter().take(nref) {
            if ra.sites.is_empty() {
                continue;
            }
            let s = S!(ra);
            let nseq = s.nseq - s.nfiltered - s.nfailprimer;
            for i in 0..ra.namp {
                let amp = &ra.amp[i];
                let alen = amp.min_right - amp.max_left + 1;
                let v = if nseq != 0 {
                    s.nbases[i] as f64 / alen as f64 / nfile as f64
                } else {
                    0.0
                };
                m.push_str(&format!("\t{:.1}", v));
            }
        }
        writeln!(out, "{}", m).ok();
        let mut m = String::from("CDEPTH\tSTDDEV");
        for ra in refs.iter().take(nref) {
            if ra.sites.is_empty() {
                continue;
            }
            let s = S!(ra);
            for i in 0..ra.namp {
                let amp = &ra.amp[i];
                let alen = (amp.min_right - amp.max_left + 1) as f64;
                let n1 = s.nbases[i] as f64 / alen;
                let v = s.nbases2[i] as f64 / (alen * alen) / nfile as f64
                    - (n1 / nfile as f64).powi(2);
                m.push_str(&format!("\t{:.1}", if v > 0.0 { v.sqrt() } else { 0.0 }));
            }
        }
        writeln!(out, "{}", m).ok();
    }

    // PCOV
    if type_ch == 'F' {
        writeln!(out, "# Percentage coverage per amplicon").ok();
        writeln!(
            out,
            "# Use 'grep ^{}PCOV | cut -f 2-' to extract this part.",
            type_ch
        )
        .ok();
        let mut d = 0;
        loop {
            let mut buf = format!("FPCOV-{}\t{}", args.min_depth[d], name);
            for ra in refs.iter_mut().take(nref) {
                if ra.sites.is_empty() {
                    continue;
                }
                for i in 0..ra.namp {
                    let amp = ra.amp[i].clone();
                    let offset = amp.min_left - 1;
                    let mut covered = 0i64;
                    let mut j = amp.max_left - 1;
                    while j < amp.min_right {
                        let idx = i * ra.lstats.max_amp_len as usize + (j - offset) as usize;
                        let cv = if local {
                            ra.lstats.coverage[idx]
                        } else {
                            ra.gstats.coverage[idx]
                        };
                        if cv >= args.min_depth[d] {
                            covered += 1;
                        }
                        j += 1;
                    }
                    let alen = amp.min_right - amp.max_left + 1;
                    let pc = 100.0 * covered as f64 / alen as f64;
                    if local {
                        ra.lstats.covered_perc[i][d] = pc;
                    } else {
                        ra.gstats.covered_perc[i][d] = pc;
                    }
                    buf.push_str(&format!("\t{:.2}", pc));
                }
            }
            writeln!(out, "{}", buf).ok();
            d += 1;
            if d >= MAX_DEPTH || args.min_depth[d] == 0 {
                break;
            }
        }
    } else if combined {
        let mut d = 0;
        loop {
            let mut m = format!("CPCOV-{}\tMEAN", args.min_depth[d]);
            for ra in refs.iter().take(nref) {
                if ra.sites.is_empty() {
                    continue;
                }
                let s = S!(ra);
                for i in 0..ra.namp {
                    m.push_str(&format!("\t{:.1}", s.covered_perc[i][d] / nfile as f64));
                }
            }
            writeln!(out, "{}", m).ok();
            let mut m = format!("CPCOV-{}\tSTDDEV", args.min_depth[d]);
            for ra in refs.iter().take(nref) {
                if ra.sites.is_empty() {
                    continue;
                }
                let s = S!(ra);
                for i in 0..ra.namp {
                    let n1 = s.covered_perc[i][d] / nfile as f64;
                    let v = s.covered_perc2[i][d] / nfile as f64 - n1 * n1;
                    m.push_str(&format!("\t{:.1}", if v > 0.0 { v.sqrt() } else { 0.0 }));
                }
            }
            writeln!(out, "{}", m).ok();
            d += 1;
            if d >= MAX_DEPTH || args.min_depth[d] == 0 {
                break;
            }
        }
    }

    // DP_ALL / DP_VALID (run-length encoded).
    let rle = |depth: &[i64], len: i64, out: &mut dyn Write| {
        let mut i = 0i64;
        while i < len {
            let di = depth[i as usize];
            let (mut dmin, mut dmax) = (di, di);
            let mut dmid = (dmin + dmax) as f64 / 2.0;
            let mut low = dmid * (1.0 - args.depth_bin);
            let mut high = dmid * (1.0 + args.depth_bin);
            let mut j = i + 1;
            while j < len {
                let d = depth[j as usize];
                if (d as f64) < low || (d as f64) > high {
                    break;
                }
                if dmin > d {
                    dmin = d;
                    dmid = (dmin + dmax) as f64 / 2.0;
                    low = dmid * (1.0 - args.depth_bin);
                    high = dmid * (1.0 + args.depth_bin);
                } else if dmax < d {
                    dmax = d;
                    dmid = (dmin + dmax) as f64 / 2.0;
                    low = dmid * (1.0 - args.depth_bin);
                    high = dmid * (1.0 + args.depth_bin);
                }
                j += 1;
            }
            write!(out, "\t{},{}", dmid as i64, j - i).ok();
            i = j;
        }
    };
    writeln!(out, "# Depth per reference base for ALL data.").ok();
    writeln!(
        out,
        "# Use 'grep ^{}DP_ALL | cut -f 2-' to extract this part.",
        type_ch
    )
    .ok();
    for ra in refs.iter().take(nref) {
        if ra.sites.is_empty() {
            continue;
        }
        let s = S!(ra);
        if args.multi_ref {
            write!(out, "{}DP_ALL\t{}\t{}", type_ch, name, ra.refname).ok();
        } else {
            write!(out, "{}DP_ALL\t{}", type_ch, name).ok();
        }
        rle(&s.depth_all, ra.len, out);
        writeln!(out).ok();
    }
    writeln!(
        out,
        "# Depth per reference base for full-length valid amplicon data."
    )
    .ok();
    writeln!(
        out,
        "# Use 'grep ^{}DP_VALID | cut -f 2-' to extract this part.",
        type_ch
    )
    .ok();
    for ra in refs.iter().take(nref) {
        if ra.sites.is_empty() {
            continue;
        }
        let s = S!(ra);
        if args.multi_ref {
            write!(out, "{}DP_VALID\t{}\t{}", type_ch, name, ra.refname).ok();
        } else {
            write!(out, "{}DP_VALID\t{}", type_ch, name).ok();
        }
        rle(&s.depth_valid, ra.len, out);
        writeln!(out).ok();
    }

    // TCOORD
    writeln!(out, "# Distribution of aligned template coordinates.").ok();
    writeln!(
        out,
        "# Use 'grep ^{}TCOORD | cut -f 2-' to extract this part.",
        type_ch
    )
    .ok();
    for ra in refs.iter().take(nref) {
        if ra.sites.is_empty() {
            continue;
        }
        let s = S!(ra);
        let start_i: i64 = if nref == 1 { -1 } else { 0 };
        for i in start_i..ra.namp as i64 {
            let map = &s.tcoord[(i + 1) as usize];
            let mut tp: Vec<(u32, u32, u32, u32)> = Vec::new();
            for (&k, &v) in map {
                if v & 0xffff_ffff == 0 {
                    continue;
                }
                tp.push((
                    (k & 0xffff_ffff) as u32,
                    (k >> 32) as u32,
                    (v & 0xffff_ffff) as u32,
                    (v >> 32) as u32,
                ));
            }
            // Upstream emits in hash-iteration order (unspecified) but
            // fixtures have a single entry per amplicon here.
            tp.sort_by(|a, b| b.2.cmp(&a.2).then(a.0.cmp(&b.0)).then(a.1.cmp(&b.1)));
            write!(
                out,
                "{}TCOORD\t{}\t{}",
                type_ch,
                name,
                i + 1 + ra.first_amp as i64
            )
            .ok();
            for (s0, e0, f0, st0) in &tp {
                if (*f0 as i64) < args.tcoord_min_count {
                    continue;
                }
                write!(out, "\t{},{},{},{}", s0, e0, f0, st0).ok();
            }
            writeln!(out).ok();
        }
    }

    // AMP classification.
    writeln!(out, "# Classification of amplicon status.  Columns are").ok();
    writeln!(
        out,
        "# number with both primers from this amplicon, number with"
    )
    .ok();
    writeln!(
        out,
        "# primers from different amplicon, and number with a position"
    )
    .ok();
    writeln!(out, "# not matching any valid amplicon primer site").ok();
    writeln!(
        out,
        "# Use 'grep ^{}AMP | cut -f 2-' to extract this part.",
        type_ch
    )
    .ok();
    let mut ad = [0i64; 3];
    for ra in refs.iter().take(nref) {
        if ra.sites.is_empty() {
            continue;
        }
        let s = S!(ra);
        for i in 0..ra.namp {
            ad[0] += s.amp_dist[i][0];
            ad[1] += s.amp_dist[i][1];
            ad[2] += s.amp_dist[i][2];
        }
    }
    writeln!(
        out,
        "{}AMP\t{}\t0\t{}\t{}\t{}",
        type_ch, name, ad[0], ad[1], ad[2]
    )
    .ok();
    for ra in refs.iter().take(nref) {
        if ra.sites.is_empty() {
            continue;
        }
        let s = S!(ra);
        for i in 0..ra.namp {
            writeln!(
                out,
                "{}AMP\t{}\t{}\t{}\t{}\t{}",
                type_ch,
                name,
                i + 1 + ra.first_amp,
                s.amp_dist[i][0],
                s.amp_dist[i][1],
                s.amp_dist[i][2]
            )
            .ok();
        }
    }
}

/// Entry point for `samtools ampliconstats`.
pub fn main(argv: &[OsString]) -> ExitCode {
    let mut args = Args::default();
    let mut bedfile: Option<PathBuf> = None;
    let mut files: Vec<PathBuf> = Vec::new();

    let mut it = argv.iter().skip(1).peekable();
    while let Some(a) = it.next() {
        let s = a.to_str().unwrap_or("");
        let mut val = || it.next().and_then(|v| v.to_str()).map(|x| x.to_string());
        match s {
            "-S" | "--single-ref" => args.multi_ref = false,
            "-s" | "--use-sample-name" => {}
            "-o" | "--output" => args.out = val().map(PathBuf::from),
            "-f" | "--flag-require" => {
                args.flag_require = val().and_then(|v| parse_flag(&v)).unwrap_or(0)
            }
            "-F" | "--flag-filter" => {
                args.flag_filter = val().and_then(|v| parse_flag(&v)).unwrap_or(0)
            }
            "-m" | "--pos-margin" => {
                args.max_delta = val().and_then(|v| v.parse().ok()).unwrap_or(30)
            }
            "-D" | "--depth-bin" => {
                args.depth_bin = val().and_then(|v| v.parse().ok()).unwrap_or(0.01)
            }
            "-t" | "--tlen-adjust" => {
                args.tlen_adj = val().and_then(|v| v.parse().ok()).unwrap_or(0)
            }
            "-c" | "--tcoord-min-count" => {
                args.tcoord_min_count = val().and_then(|v| v.parse().ok()).unwrap_or(10)
            }
            "-b" | "--tcoord-bin" => {
                args.tcoord_bin = val().and_then(|v| v.parse().ok()).unwrap_or(1).max(1)
            }
            "-a" | "--max-amplicons" => {
                args.max_amp = val().and_then(|v| v.parse::<i64>().ok()).unwrap_or(1000) + 1
            }
            "-l" | "--max-amplicon-length" => {
                args.max_amp_len = val().and_then(|v| v.parse::<i64>().ok()).unwrap_or(1000) + 1
            }
            "-d" | "--min-depth" => {
                if let Some(v) = val() {
                    args.min_depth = [0; MAX_DEPTH];
                    for (k, part) in v.split(',').take(MAX_DEPTH).enumerate() {
                        args.min_depth[k] = part.parse().unwrap_or(0);
                    }
                }
            }
            "-@" | "--threads" => {
                let _ = val();
            }
            "-h" | "--help" => return ExitCode::SUCCESS,
            _ if s.starts_with('-') && s != "-" => {
                print_error("ampliconstats", format!("unsupported option `{}`", s));
                return ExitCode::from(1);
            }
            _ => {
                if bedfile.is_none() {
                    bedfile = Some(PathBuf::from(a));
                } else {
                    files.push(PathBuf::from(a));
                }
            }
        }
    }
    let Some(bedfile) = bedfile else {
        print_error("ampliconstats", "a BED file is required");
        return ExitCode::from(1);
    };
    args.argv = crate::pg::stringify_argv(argv);

    match run(&args, &bedfile, &files) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("ampliconstats", "ampliconstats failed", &e);
            ExitCode::from(1)
        }
    }
}

fn parse_flag(s: &str) -> Option<u32> {
    if let Some(h) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        u32::from_str_radix(h, 16).ok()
    } else {
        s.parse().ok()
    }
}

fn run(args: &Args, bedfile: &Path, files: &[PathBuf]) -> io::Result<()> {
    let (bed_order, bed_map) = load_bed(bedfile).map_err(io::Error::other)?;
    let _ = bed_order;
    if files.is_empty() {
        return Ok(());
    }

    let mut out: Box<dyn Write> = match &args.out {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };

    // Header / SS from the first file.
    let (header0, _) = crate::sam_compat::read_sam_records_tolerant(&files[0])?;
    let sq: Vec<(String, i64)> = header0
        .reference_sequences()
        .iter()
        .map(|(k, v)| {
            (
                String::from_utf8_lossy(k.as_ref()).into_owned(),
                usize::from(v.length()) as i64,
            )
        })
        .collect();
    let nref = sq.len();

    writeln!(out, "# Summary statistics, used for scaling the plots.").ok();
    writeln!(out, "SS\tSamtools version: samtools-rs").ok();
    writeln!(out, "SS\tCommand line: {}", args.argv).ok();
    writeln!(out, "SS\tNumber of files:\t{}", files.len()).ok();

    let mut refs: Vec<RefAmp> = Vec::with_capacity(nref);
    for (name, len) in &sq {
        if let Some(sites) = bed_map.get(name) {
            let namp = count_amplicon(sites);
            if args.multi_ref {
                writeln!(out, "SS\tNumber of amplicons:\t{}\t{}", name, namp).ok();
            } else {
                writeln!(out, "SS\tNumber of amplicons:\t{}", namp).ok();
            }
            if args.multi_ref {
                writeln!(out, "SS\tReference length:\t{}\t{}", name, len).ok();
            } else {
                writeln!(out, "SS\tReference length:\t{}", len).ok();
            }
            refs.push(RefAmp {
                namp,
                len: *len,
                sites: sites.clone(),
                amp: Vec::new(),
                refname: name.clone(),
                first_amp: 0,
                lstats: Stats::new(*len, args.max_amp, args.max_amp_len),
                gstats: Stats::new(*len, args.max_amp, args.max_amp_len),
            });
        } else {
            refs.push(RefAmp {
                namp: 0,
                len: *len,
                sites: Vec::new(),
                amp: Vec::new(),
                refname: name.clone(),
                first_amp: 0,
                lstats: Stats::new(1, args.max_amp, args.max_amp_len),
                gstats: Stats::new(1, args.max_amp, args.max_amp_len),
            });
        }
    }
    writeln!(out, "SS\tEnd of summary").ok();

    // bed2amplicon for each ref (cumulative first_amp).
    let mut offset = 0usize;
    let mut first = true;
    for ra in refs.iter_mut() {
        if ra.sites.is_empty() {
            continue;
        }
        ra.first_amp = offset;
        let mut amp = Vec::new();
        let namp = bed2amplicon(
            args,
            &ra.sites,
            &mut amp,
            first,
            &ra.refname,
            offset,
            out.as_mut(),
        )
        .map_err(io::Error::other)?;
        ra.amp = amp;
        ra.namp = namp;
        offset += namp;
        first = false;
    }

    // Per-file.
    for f in files {
        let (header, records) = crate::sam_compat::read_sam_records_tolerant(f)?;
        let fhsq: Vec<String> = header
            .reference_sequences()
            .keys()
            .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
            .collect();
        if fhsq.len() != nref {
            return Err(io::Error::other("SAM headers are not consistent"));
        }
        let sname = sample_name(f);

        for ra in refs.iter_mut() {
            ra.lstats.reset();
        }

        let mut last_ref: i64 = -9;
        let mut p2s: Vec<i64> = Vec::new();
        let mut p2e: Vec<i64> = Vec::new();
        for b in &records {
            let tid = match b.reference_sequence_id() {
                Some(t) => t as i64,
                None => continue,
            };
            if tid < 0 || tid as usize >= refs.len() {
                continue;
            }
            if last_ref != tid {
                last_ref = tid;
                let (a, c) = build_pos_lookup(args, &refs[tid as usize]);
                p2s = a;
                p2e = c;
            }
            if refs[tid as usize].sites.is_empty() {
                // still counts toward nseq for that ref
                refs[tid as usize].lstats.nseq += 1;
                continue;
            }
            accumulate_stats(args, &mut refs[tid as usize], b, &p2s, &p2e);
        }

        dump_stats(
            args,
            'F',
            &sname,
            files.len() as i64,
            &mut refs,
            nref,
            true,
            out.as_mut(),
        );

        // append_stats
        let mut all_nseq = 0i64;
        for ra in refs.iter() {
            if ra.sites.is_empty() {
                continue;
            }
            all_nseq += ra.lstats.nseq - ra.lstats.nfiltered - ra.lstats.nfailprimer;
        }
        for ra in refs.iter_mut() {
            if ra.sites.is_empty() {
                continue;
            }
            let RefAmp {
                lstats,
                gstats,
                namp,
                ..
            } = ra;
            append_lstats(lstats, gstats, *namp, all_nseq);
        }
    }

    dump_stats(
        args,
        'C',
        "COMBINED",
        files.len() as i64,
        &mut refs,
        nref,
        false,
        out.as_mut(),
    );
    out.flush()?;
    Ok(())
}

fn sample_name(path: &Path) -> String {
    let n = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
    for ext in [".bam", ".sam", ".cram"] {
        if let Some(stem) = n.strip_suffix(ext) {
            return stem.to_string();
        }
    }
    n.to_string()
}
