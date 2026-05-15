//! `samtools consensus` — simple-mode consensus (FASTA / FASTQ / pileup).
//!
//! Ports the `MODE_SIMPLE` frequency-counting path of `bam_consensus.c`
//! (`calculate_consensus_simple`) on top of the htslib-rs pileup engine.
//! Per reference position the base with the greatest summed weight wins;
//! a candidate must clear `call_fract * total_score` (`-c`, default 0.75)
//! and `min_depth` (`-d`, default 1) or it becomes `N` (or `*` for an
//! all-gap column). `-q`/`--use-qual` weights by base quality. Insertions
//! are folded in as sub-columns when `--show-ins` (default yes); a `*`
//! consensus column is emitted only with `--show-del` (default no).
//! FASTQ/pileup quality is `100 * used_score / total_score` (capped at 93
//! for the FASTQ ASCII char).
//!
//! **Scope:** `--mode simple` only (the default `recall`/Bayesian modes
//! are not ported). No reference is required (simple mode is freq-based).

use std::collections::BTreeMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use htslib_rs::alignment_compat::{PileupRead, pileup_from_alignment_paths_with_options};

use crate::diagnostics::{print_error, print_error_errno};

#[derive(Clone, Copy, PartialEq, Eq)]
enum Format {
    Fasta,
    Fastq,
    Pileup,
}

struct Config {
    input: Option<PathBuf>,
    output: Option<PathBuf>,
    format: Format,
    call_fract: f64,
    het_fract: f64,
    min_depth: usize,
    use_qual: bool,
    min_qual: u8,
    ambig: bool,
    show_del: bool,
    show_ins: bool,
    line_len: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            input: None,
            output: None,
            format: Format::Fasta,
            call_fract: 0.75,
            het_fract: 0.5,
            min_depth: 1,
            use_qual: false,
            min_qual: 0,
            ambig: false,
            show_del: false,
            show_ins: true,
            line_len: 70,
        }
    }
}

fn yes(v: Option<&str>) -> bool {
    matches!(v, Some(s) if s.starts_with('y') || s.starts_with('Y'))
}

/// Entry point for `samtools consensus`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut cfg = Config::default();
    let mut mode_simple = true;

    let mut iter = args.iter().skip(1);
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-o" | "--output" => cfg.output = iter.next().map(PathBuf::from),
            "-f" | "--format" => {
                cfg.format = match iter.next().and_then(|a| a.to_str()) {
                    Some(v) if v.eq_ignore_ascii_case("fasta") => Format::Fasta,
                    Some(v) if v.eq_ignore_ascii_case("fastq") => Format::Fastq,
                    Some(v) if v.eq_ignore_ascii_case("pileup") => Format::Pileup,
                    other => {
                        print_error("consensus", format!("unknown format {other:?}"));
                        return ExitCode::from(1);
                    }
                };
            }
            "-m" | "--mode" => {
                mode_simple = matches!(iter.next().and_then(|a| a.to_str()), Some("simple"));
            }
            "-c" | "--call-fract" => {
                cfg.call_fract = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0.75);
            }
            "-H" | "--het-fract" => {
                cfg.het_fract = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0.5);
            }
            "-d" | "--min-depth" => {
                cfg.min_depth = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(1);
            }
            "-q" | "--use-qual" => cfg.use_qual = true,
            "--no-use-qual" => cfg.use_qual = false,
            "-A" | "--ambig" => cfg.ambig = true,
            "--show-del" => cfg.show_del = yes(iter.next().and_then(|a| a.to_str())),
            "--show-ins" => cfg.show_ins = yes(iter.next().and_then(|a| a.to_str())),
            "-l" | "--line-len" => {
                cfg.line_len = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .filter(|&n| n > 0)
                    .unwrap_or(usize::MAX);
            }
            "--min-BQ" => {
                cfg.min_qual = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0);
            }
            "--no-PG" | "-a" | "-aa" => {}
            "-@" | "--threads" | "-r" | "--region" | "-T" | "--reference" => {
                let _ = iter.next();
            }
            _ if s.starts_with('-') && s != "-" => { /* tolerate */ }
            _ => cfg.input = Some(PathBuf::from(arg)),
        }
    }

    if !mode_simple {
        print_error(
            "consensus",
            "only --mode simple is supported in samtools-rs consensus",
        );
        return ExitCode::from(1);
    }

    let Some(input) = cfg.input.clone() else {
        print_error("consensus", "no input file");
        return ExitCode::from(1);
    };

    match run(&cfg, &input) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("consensus", "consensus failed", &e);
            ExitCode::from(1)
        }
    }
}

/// IUPAC table indexed by the nt16-bit combination (`A=1 C=2 G=4 T=8`)
/// plus the gap slot 16. Mirrors HTSlib's `het[]`.
const HET: &[u8; 32] = b"NACMGRSVTWYHKDBN*ac?g???t???????";

fn base4(b: u8) -> usize {
    match b.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        _ => 15, // N / ambiguous: contributes to no pure base
    }
}

// seqi -> per-pure-base weights (HTSlib seqi2{A,C,G,T}).
const SEQI2A: [u64; 16] = [0, 8, 0, 4, 0, 4, 0, 2, 0, 4, 0, 2, 0, 2, 0, 1];
const SEQI2C: [u64; 16] = [0, 0, 8, 4, 0, 0, 4, 2, 0, 0, 4, 2, 0, 0, 2, 1];
const SEQI2G: [u64; 16] = [0, 0, 0, 0, 8, 4, 4, 1, 0, 0, 0, 0, 4, 2, 2, 1];
const SEQI2T: [u64; 16] = [0, 0, 0, 0, 0, 0, 0, 0, 8, 4, 4, 2, 8, 2, 2, 1];

/// One pileup entry for the simple consensus: `code` is the nt16 value, or
/// 16 for a gap (`*`); `qual` is the base quality.
struct Entry {
    code: usize,
    qual: u8,
}

/// Ports `calculate_consensus_simple`. Returns `(base_char, quality)`.
fn consensus_simple(entries: &[Entry], cfg: &Config) -> (u8, u32) {
    let mut freq = [0u64; 17];
    let mut score = [0u64; 17];
    let mut tot_depth = 0u64;

    for e in entries {
        if e.qual < cfg.min_qual {
            continue;
        }
        let qw = if cfg.use_qual { u64::from(e.qual) } else { 1 };
        if e.code < 16 {
            for (slot, table) in [(1usize, &SEQI2A), (2, &SEQI2C), (4, &SEQI2G), (8, &SEQI2T)] {
                let q = table[e.code] * qw;
                if q != 0 {
                    freq[slot] += 1;
                    score[slot] += q;
                }
            }
        } else {
            freq[16] += 1;
            score[16] += 8 * qw;
        }
        tot_depth += 1;
    }

    let tscore: u64 = [1usize, 2, 4, 8, 16].iter().map(|&c| score[c]).sum();

    let (mut call1, mut call2) = (15usize, 15usize);
    let (mut score1, mut score2) = (0u64, 0u64);
    for &c in &[1usize, 2, 4, 8, 16] {
        if score1 < score[c] {
            score2 = score1;
            call2 = call1;
            score1 = score[c];
            call1 = c;
        } else if score2 < score[c] {
            score2 = score[c];
            call2 = c;
        }
    }

    let mut used_base = call1;
    let mut used_score = score1;
    if cfg.ambig && score1 > 0 && (score2 as f64) >= cfg.het_fract * (score1 as f64) {
        used_base |= call2;
        used_score += score2;
    }

    if (tot_depth as usize) < cfg.min_depth
        || (used_score as f64) < cfg.call_fract * (tscore as f64)
    {
        used_base = if call1 == 16 { 16 } else { 0 };
    }

    let ch = HET[used_base & 31];
    let qual = if used_base != 0 && tscore > 0 {
        (100.0 * used_score as f64 / tscore as f64) as u32
    } else {
        0
    };
    (ch, qual)
}

fn ref_entries(reads: &[PileupRead]) -> Vec<Entry> {
    reads
        .iter()
        .filter(|r| !r.is_refskip)
        .map(|r| Entry {
            code: match r.base {
                Some(b) => base4(b),
                None => 16, // deletion → gap
            },
            qual: r.qpos_quality,
        })
        .collect()
}

fn insertion_entries(reads: &[PileupRead], nth: usize) -> Vec<Entry> {
    reads
        .iter()
        .filter(|r| !r.is_refskip)
        .map(|r| {
            let code = if r.indel > 0 && r.insertion.len() >= nth {
                base4(r.insertion[nth - 1])
            } else {
                16 // no insertion here → gap
            };
            Entry {
                code,
                qual: r.qpos_quality,
            }
        })
        .collect()
}

fn fastq_qual_char(q: u32) -> u8 {
    if q > 93 { 126 } else { (q as u8) + 33 }
}

struct RefSeq {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u32>,
}

fn run(cfg: &Config, input: &PathBuf) -> io::Result<()> {
    let columns =
        pileup_from_alignment_paths_with_options(std::slice::from_ref(input), &Default::default())?;

    let mut writer: Box<dyn Write> = match cfg.output.as_ref() {
        Some(p) => Box::new(io::BufWriter::new(File::create(p)?)),
        None => Box::new(io::BufWriter::new(io::stdout().lock())),
    };

    // Group columns by reference, preserving first-seen order.
    let mut order: Vec<String> = Vec::new();
    let mut by_ref: BTreeMap<String, RefSeq> = BTreeMap::new();
    let mut pileup_rows: Vec<u8> = Vec::new();

    for col in &columns {
        let reads = &col.reads_by_input[0];
        if !by_ref.contains_key(&col.reference_name) {
            order.push(col.reference_name.clone());
            by_ref.insert(
                col.reference_name.clone(),
                RefSeq {
                    name: col.reference_name.clone(),
                    seq: Vec::new(),
                    qual: Vec::new(),
                },
            );
        }
        let rs = by_ref.get_mut(&col.reference_name).unwrap();

        // Reference position (nth = 0).
        let entries = ref_entries(reads);
        let depth = entries.iter().filter(|e| e.qual >= cfg.min_qual).count();
        let (cb, cq) = consensus_simple(&entries, cfg);

        if cfg.format == Format::Pileup {
            if cb != b'*' || cfg.show_del {
                emit_pileup_row(
                    &mut pileup_rows,
                    &col.reference_name,
                    col.position,
                    0,
                    depth,
                    cb,
                    cq,
                    reads,
                    None,
                );
            }
        } else if cb != b'*' || cfg.show_del {
            rs.seq.push(cb);
            rs.qual.push(cq);
        }

        if !cfg.show_ins {
            continue;
        }

        // Insertion sub-columns (nth = 1..max insertion length).
        let max_ins = reads
            .iter()
            .filter(|r| !r.is_refskip && r.indel > 0)
            .map(|r| r.insertion.len())
            .max()
            .unwrap_or(0);
        for nth in 1..=max_ins {
            let ins = insertion_entries(reads, nth);
            let idepth = ins.iter().filter(|e| e.qual >= cfg.min_qual).count();
            let (ib, iq) = consensus_simple(&ins, cfg);
            if cfg.format == Format::Pileup {
                emit_pileup_row(
                    &mut pileup_rows,
                    &col.reference_name,
                    col.position,
                    nth,
                    idepth,
                    ib,
                    iq,
                    reads,
                    Some(nth),
                );
            } else if ib != b'*' {
                rs.seq.push(ib);
                rs.qual.push(iq);
            }
        }
    }

    if cfg.format == Format::Pileup {
        writer.write_all(&pileup_rows)?;
        return writer.flush();
    }

    for name in &order {
        let rs = &by_ref[name];
        // Trim leading/trailing all-N as HTSlib does for un-padded output.
        let start = rs
            .seq
            .iter()
            .position(|&b| b != b'N')
            .unwrap_or(rs.seq.len());
        let end = rs.seq.iter().rposition(|&b| b != b'N').map_or(0, |i| i + 1);
        if start >= end {
            continue;
        }
        let seq = &rs.seq[start..end];
        let qual = &rs.qual[start..end];
        let lead = if cfg.format == Format::Fastq {
            b'@'
        } else {
            b'>'
        };
        writer.write_all(&[lead])?;
        writer.write_all(rs.name.as_bytes())?;
        writer.write_all(b"\n")?;
        for chunk in seq.chunks(cfg.line_len) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
        if cfg.format == Format::Fastq {
            writer.write_all(b"+\n")?;
            let q: Vec<u8> = qual.iter().map(|&q| fastq_qual_char(q)).collect();
            for chunk in q.chunks(cfg.line_len) {
                writer.write_all(chunk)?;
                writer.write_all(b"\n")?;
            }
        }
    }

    writer.flush()
}

#[allow(clippy::too_many_arguments)]
fn emit_pileup_row(
    out: &mut Vec<u8>,
    name: &str,
    pos: usize,
    nth: usize,
    depth: usize,
    cons: u8,
    score: u32,
    reads: &[PileupRead],
    insertion_nth: Option<usize>,
) {
    let mut bases = Vec::new();
    let mut quals = Vec::new();
    for r in reads {
        if r.is_refskip {
            continue;
        }
        let c = match insertion_nth {
            None => match r.base {
                Some(b) => b.to_ascii_uppercase(),
                None => b'*',
            },
            Some(n) => {
                if r.indel > 0 && r.insertion.len() >= n {
                    r.insertion[n - 1].to_ascii_uppercase()
                } else {
                    b'*'
                }
            }
        };
        bases.push(c);
        quals.push(fastq_qual_char(u32::from(r.qpos_quality)));
    }
    out.extend_from_slice(name.as_bytes());
    let _ = write!(out, "\t{pos}\t{nth}\t{depth}\t{}\t{score}\t", cons as char);
    out.extend_from_slice(&bases);
    out.push(b'\t');
    out.extend_from_slice(&quals);
    out.push(b'\n');
}
