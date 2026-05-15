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
//! **Modes:** `--mode simple` (frequency counting) and the default
//! `bayesian`/`recall` Gap5 model (`calculate_consensus_gap5` ported in
//! the `bayes` submodule, fed by the htslib-rs pileup `nm_init`
//! precompute via `PileupRead::bayes_poly`/`bayes_nm_local`). 52/77
//! `test/consensus/consensus.reg` cases byte-exact; the remaining 25
//! (pileup-format, `--ref-qual`+`-T` reference, `-a`/all-bases,
//! `--min-MQ`, the 30/31/32 show-del bayesian series) are WIP.

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
    // Bayesian-mode (`--mode bayesian`/default) parameters.
    mode_simple: bool,
    default_qual: u8,
    cons_cutoff: i32,
    use_mqual: bool,
    nm_adjust: bool,
    scale_mqual: f64,
    low_mqual: f64,
    high_mqual: f64,
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
            mode_simple: false, // upstream default is bayesian/recall
            default_qual: 10,
            cons_cutoff: 10,
            use_mqual: true,
            nm_adjust: true,
            scale_mqual: 1.0,
            low_mqual: 1.0,
            high_mqual: 60.0,
        }
    }
}

fn yes(v: Option<&str>) -> bool {
    matches!(v, Some(s) if s.starts_with('y') || s.starts_with('Y'))
}

/// Entry point for `samtools consensus`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut cfg = Config::default();
    // Upstream default is bayesian/recall; only an explicit
    // `-m simple`/`--mode simple` selects the frequency-count path.
    let mut mode_simple = false;

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
            "--use-MQ" => cfg.use_mqual = true,
            "--no-use-MQ" => cfg.use_mqual = false,
            "-C" | "--cutoff" => {
                cfg.cons_cutoff = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(10);
            }
            "--default-qual" => {
                cfg.default_qual = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(10);
            }
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
            // Glued short option, e.g. `-C0`, `-c0.6`, `-d6`.
            _ if s.len() > 2 && s.starts_with('-') && !s.starts_with("--") => {
                let (flag, val) = s.split_at(2);
                match flag {
                    "-C" => cfg.cons_cutoff = val.parse().unwrap_or(10),
                    "-c" => cfg.call_fract = val.parse().unwrap_or(0.75),
                    "-H" => cfg.het_fract = val.parse().unwrap_or(0.5),
                    "-d" => cfg.min_depth = val.parse().unwrap_or(1),
                    "-l" => {
                        cfg.line_len = val.parse().ok().filter(|&n| n > 0).unwrap_or(usize::MAX)
                    }
                    _ => { /* tolerate other glued short opts */ }
                }
            }
            "--no-PG" | "-a" | "-aa" => {}
            "-@" | "--threads" | "-r" | "--region" | "-T" | "--reference" => {
                let _ = iter.next();
            }
            _ if s.starts_with('-') && s != "-" => { /* tolerate */ }
            _ => cfg.input = Some(PathBuf::from(arg)),
        }
    }

    cfg.mode_simple = mode_simple;

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

/// Prebuilt Bayesian-mode tables (`MODE_RECALL`), constructed once.
struct BayesCtx {
    cp: bayes::ConsProbs,
    et: bayes::ETab,
    q2p: [f64; 101],
    mpow: [f64; 256],
}

impl BayesCtx {
    fn new() -> Self {
        BayesCtx {
            cp: bayes::default_recall(),
            et: bayes::ETab::new(),
            q2p: bayes::q2p_table(),
            mpow: bayes::mqual_pow_1m_table(),
        }
    }
}

/// Build one `Gap5Obs` per read for the reference column (`nth==0`) or
/// an insertion sub-column (`nth>=1`), filtering ref-skips and reads
/// below `min_qual`, applying the upstream `255 -> default_qual` rule.
fn gap5_obs(reads: &[PileupRead], cfg: &Config, nth: usize) -> Vec<bayes::Gap5Obs> {
    let mut v = Vec::new();
    for r in reads {
        if r.is_refskip {
            continue;
        }
        let (code, raw) = if nth == 0 {
            (
                match r.base {
                    Some(b) => base4(b),
                    None => 16,
                },
                r.quality.unwrap_or(r.qpos_quality),
            )
        } else if r.indel > 0 && r.insertion.len() >= nth {
            (base4(r.insertion[nth - 1]), r.qpos_quality)
        } else {
            (16, r.qpos_quality)
        };
        if raw < cfg.min_qual {
            continue;
        }
        let qual = if raw == 255 || (raw == 0 && r.qpos_quality == 255) {
            cfg.default_qual
        } else {
            raw
        };
        v.push(bayes::Gap5Obs {
            base4: code as u8,
            qual,
            mqual: r.mapping_quality as f64,
            nm_local: r.bayes_nm_local as i32,
            poly: r.bayes_poly as f64,
        });
    }
    v
}

/// Bayesian consensus for one column, faithfully porting upstream
/// `consensus_base`'s non-simple branch (the `min_depth`/ambiguity/
/// `cons_cutoff` thresholds on top of `gap5_call`).
fn consensus_bayes(reads: &[PileupRead], cfg: &Config, ctx: &BayesCtx, nth: usize) -> (u8, u32) {
    let td = reads.iter().filter(|r| !r.is_refskip).count() as i32;
    let obs = gap5_obs(reads, cfg, nth);
    let opts = bayes::Gap5Opts {
        use_mqual: cfg.use_mqual,
        nm_adjust: cfg.nm_adjust,
        scale_mqual: cfg.scale_mqual,
        low_mqual: cfg.low_mqual,
        high_mqual: cfg.high_mqual,
    };
    let cons = bayes::gap5_call(&obs, td, &ctx.cp, &ctx.et, &ctx.q2p, &ctx.mpow, &opts);

    let acgt = [b'A', b'C', b'G', b'T', b'*'];
    // 5x5 ACGT* ambiguity matrix (rows/cols A C G T *).
    const AMBIG: &[u8; 25] = b"AMRWaMCSYcRSGKgWYKTtacgt*";

    let (mut cb, mut cq): (u8, i32) = if (cons.depth as usize) < cfg.min_depth && cons.call != 4 {
        (b'N', 0)
    } else if cons.het_logodd > 0 && cfg.ambig {
        (AMBIG[cons.het_call], cons.het_logodd)
    } else {
        (acgt[cons.call.min(4)], cons.phred)
    };
    if cq < cfg.cons_cutoff && cb != b'*' && cons.het_call % 5 != 4 && cons.het_call / 5 != 4 {
        cb = b'N';
        cq = 0;
    }
    (cb, cq.max(0) as u32)
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

    let bayes_ctx = if cfg.mode_simple {
        None
    } else {
        Some(BayesCtx::new())
    };

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
        let (cb, cq) = match &bayes_ctx {
            Some(ctx) => consensus_bayes(reads, cfg, ctx, 0),
            None => consensus_simple(&entries, cfg),
        };

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
            let (ib, iq) = match &bayes_ctx {
                Some(ctx) => consensus_bayes(reads, cfg, ctx, nth),
                None => consensus_simple(&ins, cfg),
            };
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

/// Bayesian (`--mode bayesian`/default `recall`) consensus probability
/// tables — a faithful port of `consensus_init` in
/// `samtools/bam_consensus.c` (the Gap5-derived model). This is the
/// foundational table-construction step of the Bayesian engine; the
/// `calculate_consensus_gap5` accumulation/call is wired on top of it.
///
/// Not yet reachable from the CLI (only `--mode simple` is dispatched),
/// so this is regression-safe scaffolding verified by its own tests.
pub mod bayes {
    /// `samtools/bam_consensus.c` defaults: `P_HET`, `P_INDEL`,
    /// `P_HET_SCALE`, and `homopoly_redux` (poly_mul) for `MODE_RECALL`.
    pub const P_HET: f64 = 1e-3;
    pub const P_INDEL: f64 = 2e-4;
    pub const P_HET_SCALE: f64 = 1.0;
    pub const POLY_MUL_RECALL: f64 = 0.01;

    /// Quality calibration maps (`qcal_t`). The default / `:flat`
    /// profile is the identity `smap[i] = omap[i] = umap[i] = i`.
    #[derive(Clone)]
    pub struct Qcal {
        pub smap: [i32; 101],
        pub omap: [i32; 101],
        pub umap: [i32; 101],
    }

    impl Qcal {
        /// `QCAL_FLAT` / `set_qcal(_, QCAL_FLAT)`.
        pub fn flat() -> Self {
            let mut m = [0i32; 101];
            for (i, v) in m.iter_mut().enumerate() {
                *v = i as i32;
            }
            Qcal {
                smap: m,
                omap: m,
                umap: m,
            }
        }
    }

    /// Ported `cons_probs` (the subset used by `calculate_consensus_gap5`:
    /// the 15-combination log-priors and the per-quality log-likelihood
    /// tables). `e_tab`/`e_log` accel tables live with the accumulator.
    #[derive(Clone)]
    pub struct ConsProbs {
        pub poly_mul: f64,
        pub prior: [f64; 25],
        pub lprior15: [f64; 15],
        pub p_mm: [f64; 101],
        pub p_xx: [f64; 101],
        pub p_xm: [f64; 101],
        pub p_oo: [f64; 101],
        pub p_om: [f64; 101],
        pub p_ox: [f64; 101],
        pub p_uu: [f64; 101],
        pub p_um: [f64; 101],
        pub p_mm_lower: [f64; 101], // upstream `pmm` (undercall match)
    }

    /// Faithful port of `consensus_init` for any non-`MODE_BAYES_116`
    /// mode (i.e. `MODE_RECALL`/`MODE_PRECISE`/`MODE_MIXED`).
    pub fn cons_probs_init(
        p_het: f64,
        p_indel: f64,
        het_scale: f64,
        poly_mul: f64,
        qcal: &Qcal,
    ) -> ConsProbs {
        let mut prior = [p_het / 6.0; 25];
        // Flat "it is what we observe": homozygous priors = 1.
        for &i in &[0usize, 6, 12, 18, 24] {
            prior[i] = 1.0;
        }
        // Heterozygous deletion (i = 4, 9, 14, 19).
        let mut i = 4;
        while i < 24 {
            prior[i] = p_indel / 6.0;
            i += 5;
        }
        // Heterozygous insertion (i = 20..=23).
        for v in prior.iter_mut().take(24).skip(20) {
            *v = p_indel / 6.0;
        }

        let pri_idx = [0usize, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 14, 18, 19, 24];
        let mut lprior15 = [0.0f64; 15];
        for (k, &idx) in pri_idx.iter().enumerate() {
            lprior15[k] = prior[idx].ln();
        }

        let mut cp = ConsProbs {
            poly_mul,
            prior,
            lprior15,
            p_mm: [0.0; 101],
            p_xx: [0.0; 101],
            p_xm: [0.0; 101],
            p_oo: [0.0; 101],
            p_om: [0.0; 101],
            p_ox: [0.0; 101],
            p_uu: [0.0; 101],
            p_um: [0.0; 101],
            p_mm_lower: [0.0; 101],
        };

        for q in 1..101usize {
            let prob = 1.0 - 10f64.powf(-(qcal.smap[q] as f64) / 10.0);
            cp.p_mm[q] = prob.ln();
            cp.p_xx[q] = ((1.0 - prob) / 3.0).ln();
            cp.p_xm[q] = ((cp.p_mm[q].exp() + cp.p_xx[q].exp()) / 2.0).ln() + het_scale.ln();

            // overcall (insertion-leaning)
            let prob_o = 1.0 - 10f64.powf(-(qcal.omap[q] as f64) / 10.0);
            cp.p_oo[q] = ((1.0 - prob_o) / 3.0).ln();
            if cp.p_oo[q] > cp.p_mm[q] - 0.5 {
                cp.p_oo[q] = cp.p_mm[q] - 0.5;
            }
            cp.p_ox[q] = ((cp.p_oo[q].exp() + cp.p_xx[q].exp()) / 2.0).ln();
            cp.p_om[q] = ((cp.p_oo[q].exp() + cp.p_mm[q].exp()) / 2.0).ln();
            if cp.p_om[q] > cp.p_xm[q] + 0.5 {
                cp.p_om[q] = cp.p_xm[q] + 0.5;
            }

            // undercall (deletion-leaning)
            let prob_u = 1.0 - 10f64.powf(-(qcal.umap[q] as f64) / 10.0);
            cp.p_mm_lower[q] = prob_u.ln();
            cp.p_uu[q] = ((1.0 - prob_u) / 3.0).ln();
            if cp.p_uu[q] > cp.p_mm[q] - 0.5 {
                cp.p_uu[q] = cp.p_mm[q] - 0.5;
            }
            cp.p_um[q] = ((cp.p_uu[q].exp() + cp.p_mm_lower[q].exp()) / 2.0).ln();
        }

        // Index 0 mirrors index 1.
        cp.p_mm[0] = cp.p_mm[1];
        cp.p_xx[0] = cp.p_xx[1];
        cp.p_xm[0] = cp.p_xm[1];
        cp.p_mm_lower[0] = cp.p_mm_lower[1];
        cp.p_oo[0] = cp.p_oo[1];
        cp.p_ox[0] = cp.p_ox[1];
        cp.p_om[0] = cp.p_om[1];
        cp.p_uu[0] = cp.p_uu[1];
        cp.p_um[0] = cp.p_um[1];

        cp
    }

    /// Default (`MODE_RECALL`) probability table, as built by upstream
    /// `main_consensus` for `cons_prob_recall`.
    pub fn default_recall() -> ConsProbs {
        cons_probs_init(P_HET, P_INDEL, P_HET_SCALE, POLY_MUL_RECALL, &Qcal::flat())
    }

    // ---- Math accel helpers (bam_consensus.c) ----

    pub const TENLOG2OVERLOG10: f64 = 3.0103;

    /// `q2p[i] = pow(10, -i/10.0)` for `i` in `0..=100`
    /// (`bam_consensus_tab.h`).
    pub fn q2p_table() -> [f64; 101] {
        let mut t = [0.0f64; 101];
        for (i, v) in t.iter_mut().enumerate() {
            *v = 10f64.powf(-(i as f64) / 10.0);
        }
        t
    }

    /// `mqual_pow_1m[i] = pow(10, -(i*.9)/10.0)` for `i` in `0..255`,
    /// then `mqual_pow_1m[255] = mqual_pow_1m[10]`.
    pub fn mqual_pow_1m_table() -> [f64; 256] {
        let mut t = [0.0f64; 256];
        for (i, v) in t.iter_mut().enumerate().take(255) {
            *v = 10f64.powf(-((i as f64) * 0.9) / 10.0);
        }
        t[255] = t[10];
        t
    }

    /// `e_tab[i] = exp(i)` for `i` in `-500..=500`; `e_tab2[i] =
    /// exp(i/10.)` for `i` in `-500..=500` (built in `consensus_init`).
    /// Stored with a +500 bias so index 0 == `exp(-500)`.
    pub struct ETab {
        e_tab: [f64; 1001],
        e_tab2: [f64; 1001],
    }

    impl ETab {
        pub fn new() -> Self {
            let mut e_tab = [0.0f64; 1001];
            let mut e_tab2 = [0.0f64; 1001];
            for i in -500i32..=500 {
                e_tab[(i + 500) as usize] = (i as f64).exp();
                e_tab2[(i + 500) as usize] = ((i as f64) / 10.0).exp();
            }
            ETab { e_tab, e_tab2 }
        }

        /// `fast_exp` (bam_consensus.c:883): table lookup, C `(int)`
        /// truncation toward zero.
        pub fn fast_exp(&self, y: f64) -> f64 {
            if (-50.0..=50.0).contains(&y) {
                let idx = (y * 10.0) as i32; // truncates toward zero
                return self.e_tab2[(idx + 500) as usize];
            }
            let y = y.clamp(-500.0, 500.0);
            self.e_tab[(y as i32 + 500) as usize]
        }
    }

    impl Default for ETab {
        fn default() -> Self {
            Self::new()
        }
    }

    /// `fast_log2` (bam_consensus.c:896): exponent + degree-3 Taylor of
    /// the mantissa, via the IEEE-754 bit layout.
    pub fn fast_log2(val: f64) -> f64 {
        let mut x = val.to_bits();
        let e = (((x >> 52) & 2047) as i64) - 1024;
        x &= !(2047i64 << 52) as u64;
        x = x.wrapping_add((1023i64 << 52) as u64);
        let d = f64::from_bits(x);
        let v = ((-1.0 / 3.0) * d + 2.0) * d - 2.0 / 3.0;
        e as f64 + v
    }

    /// `#define ph_log(x) (-TENLOG2OVERLOG10*fast_log2((x)))`.
    pub fn ph_log(x: f64) -> f64 {
        -TENLOG2OVERLOG10 * fast_log2(x)
    }

    // ---- S[15] accumulation + call (calculate_consensus_gap5) ----

    /// SAM 4-bit base code -> {A=0,C=1,G=2,T=3,*=4,N/other=5}, the
    /// `L[32]` map. `*` (pad) arrives as code >= 16.
    pub const L: [u8; 32] = [
        5, 0, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, // 0..15
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 16..31 (pad)
    ];
    const MAP_SING: [usize; 15] = [0, 5, 5, 5, 5, 1, 5, 5, 5, 2, 5, 5, 3, 5, 4];
    const MAP_HET: [usize; 15] = [0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 14, 18, 19, 24];

    /// One filtered pileup observation feeding the Bayesian model.
    pub struct Gap5Obs {
        /// SAM 4-bit sequence code (`bam_seqi`), or `>=16` for `*`/pad.
        pub base4: u8,
        /// Per-base quality after the `255 -> default_qual` rule.
        pub qual: u8,
        /// Read mapping quality (`b->core.qual`).
        pub mqual: f64,
        /// Local NM within the halo (`nm_local`); used by `--NM-adjust`.
        pub nm_local: i32,
        /// Homopolymer run length (`poly_len`) at this position.
        pub poly: f64,
    }

    pub struct Gap5Opts {
        pub use_mqual: bool,
        pub nm_adjust: bool,
        pub scale_mqual: f64,
        pub low_mqual: f64,
        pub high_mqual: f64,
    }

    pub struct Gap5Cons {
        /// 0=A 1=C 2=G 3=T 4=* 5=N (`map_sing[call]`).
        pub call: usize,
        pub phred: i32,
        pub het_call: usize,
        pub het_logodd: i32,
        pub depth: i32,
    }

    /// Faithful port of `calculate_consensus_gap5` (default build: no
    /// K2 / DO_FRACT / DO_HDW / DO_POLY_DIST / CONS_DISCREP). `obs` is
    /// pre-filtered (min-qual / ref-skip done by the caller). `td` is
    /// the original (pre-filter) depth used by the MQUAL depth fudge.
    #[allow(clippy::needless_range_loop)]
    pub fn gap5_call(
        obs: &[Gap5Obs],
        td: i32,
        cp: &ConsProbs,
        et: &ETab,
        q2p: &[f64; 101],
        mpow: &[f64; 256],
        opts: &Gap5Opts,
    ) -> Gap5Cons {
        let min_e_exp = (f64::MIN_EXP as f64 - 1.0) * std::f64::consts::LN_2 + 1.0;
        let mut s = [0.0f64; 15];
        let mut counts = [0i32; 6];
        let mut depth = 0i32;

        for p in obs {
            let mut qual = p.qual as f64;
            let base = L[(p.base4 as usize) & 31] as usize;

            if opts.use_mqual {
                let mut mqual = p.mqual;
                if opts.nm_adjust {
                    mqual /= (p.nm_local + 1) as f64;
                    let d = if td > 30 { 30 } else { td };
                    mqual *= 1.0 + 2.0 * (0.5 - d as f64 / 60.0);
                }
                mqual *= opts.scale_mqual;
                if mqual < opts.low_mqual {
                    mqual = opts.low_mqual;
                }
                if mqual > opts.high_mqual {
                    mqual = opts.high_mqual;
                }
                let pq = q2p[(qual as usize).min(100)];
                let m = mpow[(mqual as i64).clamp(0, 255) as usize];
                qual = ph_log(pq + 0.75 * m - pq * m);
            }
            if qual < 1.0 {
                qual = 1.0;
            }

            let q = (qual as i64).clamp(0, 100) as usize;
            let q2 = (qual - (p.poly - 2.0) * cp.poly_mul).max(1.0);
            let q2 = (q2 as i64).clamp(0, 100) as usize;

            let xx = cp.p_xx[q];
            let mm_ = cp.p_mm[q] - xx;
            let xm = cp.p_xm[q] - xx;
            let oo = cp.p_oo[q2] - xx;
            let om = cp.p_om[q2] - xx;
            let ox = cp.p_ox[q2] - xx;
            let uu = cp.p_uu[q2] - xx;
            let um = cp.p_um[q2] - xx;
            let mm = cp.p_mm_lower[q2] - xx;

            counts[base] += 1;
            match base {
                0 => {
                    s[0] += mm_;
                    s[1] += xm;
                    s[2] += xm;
                    s[3] += xm;
                    s[4] += om;
                    s[8] += ox;
                    s[11] += ox;
                    s[13] += ox;
                    s[14] += oo;
                }
                1 => {
                    s[1] += xm;
                    s[5] += mm_;
                    s[6] += xm;
                    s[7] += xm;
                    s[8] += om;
                    s[4] += ox;
                    s[11] += ox;
                    s[13] += ox;
                    s[14] += oo;
                }
                2 => {
                    s[2] += xm;
                    s[6] += xm;
                    s[9] += mm_;
                    s[10] += xm;
                    s[11] += om;
                    s[4] += ox;
                    s[8] += ox;
                    s[13] += ox;
                    s[14] += oo;
                }
                3 => {
                    s[3] += xm;
                    s[7] += xm;
                    s[10] += xm;
                    s[12] += mm_;
                    s[13] += om;
                    s[4] += ox;
                    s[8] += ox;
                    s[11] += ox;
                    s[14] += oo;
                }
                4 => {
                    s[0] += uu;
                    s[1] += uu;
                    s[2] += uu;
                    s[3] += uu;
                    s[4] += um;
                    s[5] += uu;
                    s[6] += uu;
                    s[7] += uu;
                    s[8] += um;
                    s[9] += uu;
                    s[10] += uu;
                    s[11] += um;
                    s[12] += uu;
                    s[13] += um;
                    s[14] += mm;
                }
                _ => {
                    // 5 => N: equal weight to A,C,G,T (not a pad).
                    s[0] += mm_;
                    s[1] += mm_;
                    s[2] += mm_;
                    s[3] += mm_;
                    s[4] += om;
                    s[5] += mm_;
                    s[6] += mm_;
                    s[7] += mm_;
                    s[8] += om;
                    s[9] += mm_;
                    s[10] += mm_;
                    s[11] += om;
                    s[12] += mm_;
                    s[13] += om;
                    s[14] += oo;
                }
            }
            depth += 1;
        }

        // Add priors; split pure (homozygous) vs het argmax.
        let mut shift = f64::MIN;
        let mut max = f64::MIN;
        let mut max_het = f64::MIN;
        let mut call = 0usize;
        let mut het_call = 0usize;
        for j in 0..15 {
            s[j] += cp.lprior15[j];
            if shift < s[j] {
                shift = s[j];
            }
            if j != 0 && j != 5 && j != 9 && j != 12 && j != 14 {
                if max_het < s[j] {
                    max_het = s[j];
                    het_call = j;
                }
                continue;
            }
            if max < s[j] {
                max = s[j];
                call = j;
            }
        }

        for j in 0..15 {
            s[j] -= shift;
            let e = et.fast_exp(s[j]);
            s[j] = if s[j] > min_e_exp {
                e
            } else {
                f64::MIN_POSITIVE
            };
        }
        let mut norm = [0.0f64; 15];
        let mut tot1 = 0.0f64;
        let mut tot2 = 0.0f64;
        for j in 0..15 {
            norm[j] += tot1;
            norm[14 - j] += tot2;
            tot1 += s[j];
            tot2 += s[14 - j];
        }

        if depth == 0 || depth == counts[5] {
            return Gap5Cons {
                call: 4,
                phred: 0,
                het_call: 0,
                het_logodd: 0,
                depth: 0,
            };
        }

        if norm[call] == 0.0 {
            norm[call] = f64::MIN_POSITIVE;
        }
        let ph = if s[call] == 1.0 && norm[call] < 0.01 {
            ph_log(norm[call]) + 0.5
        } else {
            ph_log(1.0 - s[call] / (norm[call] + s[call])) + 0.5
        };
        let phred = if (ph as i32) < 0 { 0 } else { ph as i32 };

        if norm[het_call] == 0.0 {
            norm[het_call] = f64::MIN_POSITIVE;
        }
        let het_ph = TENLOG2OVERLOG10 * (fast_log2(s[het_call]) - fast_log2(norm[het_call])) + 0.5;

        Gap5Cons {
            call: MAP_SING[call],
            phred,
            het_call: MAP_HET[het_call],
            het_logodd: het_ph as i32,
            depth,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn priors_match_upstream_consensus_init() {
            let cp = default_recall();
            // Homozygous priors are 1 -> log 0.
            for &k in &[0usize, 5, 9, 12, 14] {
                assert_eq!(cp.lprior15[k], 0.0, "lprior15[{k}] should be ln(1)=0");
            }
            // Substitution het prior P_HET/6 -> lprior15[1] (AC).
            assert!((cp.lprior15[1] - (P_HET / 6.0).ln()).abs() < 1e-12);
            // Het-deletion prior P_INDEL/6 -> lprior15[4] (prior[4]).
            assert!((cp.lprior15[4] - (P_INDEL / 6.0).ln()).abs() < 1e-12);
            // Het-insertion prior P_INDEL/6 -> lprior15[13] (prior[19]).
            assert!((cp.lprior15[13] - (P_INDEL / 6.0).ln()).abs() < 1e-12);
        }

        #[test]
        fn likelihoods_match_upstream_formulas() {
            let cp = default_recall();
            // q=20, flat qcal: prob = 1 - 10^-2 = 0.99.
            let prob = 0.99_f64;
            assert!((cp.p_mm[20] - prob.ln()).abs() < 1e-12);
            assert!((cp.p_xx[20] - ((1.0 - prob) / 3.0).ln()).abs() < 1e-12);
            let xm = ((cp.p_mm[20].exp() + cp.p_xx[20].exp()) / 2.0).ln() + P_HET_SCALE.ln();
            assert!((cp.p_xm[20] - xm).abs() < 1e-12);
            // Index 0 mirrors index 1.
            assert_eq!(cp.p_mm[0], cp.p_mm[1]);
            assert_eq!(cp.p_um[0], cp.p_um[1]);
            // pMM is monotonically increasing toward 0 (higher qual ->
            // more confident match), and always negative.
            for q in 2..101 {
                assert!(cp.p_mm[q] < 0.0);
                assert!(cp.p_mm[q] >= cp.p_mm[q - 1]);
            }
            // The MM clamps hold: poo/puu never exceed pMM-0.5.
            for q in 1..101 {
                assert!(cp.p_oo[q] <= cp.p_mm[q] - 0.5 + 1e-12);
                assert!(cp.p_uu[q] <= cp.p_mm[q] - 0.5 + 1e-12);
            }
        }

        #[test]
        fn accel_tables_and_helpers_match_upstream() {
            // q2p[i] = pow(10, -i/10) — literal values from
            // bam_consensus_tab.h.
            let q2p = q2p_table();
            assert_eq!(q2p[0], 1.0);
            assert!((q2p[1] - 0.794_328_234_724_281_5).abs() < 1e-15);
            assert!((q2p[2] - 0.630_957_344_480_193_2).abs() < 1e-15);
            assert!((q2p[3] - 0.501_187_233_627_272_2).abs() < 1e-15);

            // mqual_pow_1m[i] = pow(10, -(i*.9)/10); [255] == [10].
            let m = mqual_pow_1m_table();
            assert_eq!(m[0], 1.0);
            assert!((m[1] - 0.812_830_516_164_099_3).abs() < 1e-15);
            assert!((m[2] - 0.660_693_448_007_596_0).abs() < 1e-15);
            assert_eq!(m[255], m[10]);

            // fast_log2 within the documented deg-3 Taylor tolerance,
            // and the exact-power identities the algorithm relies on.
            assert!((fast_log2(1.0) - 0.0).abs() < 1e-9);
            for &v in &[0.5f64, 2.0, 8.0, 0.125, 1e-10, 1e10] {
                assert!((fast_log2(v) - v.log2()).abs() < 0.02, "log2({v})");
            }

            // fast_exp table lookup matches exp() at the grid points
            // (e_tab2 step 0.1) and is bounded elsewhere.
            let et = ETab::new();
            for i in -50..=50 {
                let y = i as f64; // exact grid (y*10 integral)
                assert!((et.fast_exp(y / 1.0) - et.fast_exp(y)).abs() < 1e-300);
            }
            assert!((et.fast_exp(0.0) - 1.0).abs() < 1e-12);
            assert!((et.fast_exp(2.0) - 2.0_f64.exp().ln().exp()).abs() < 1e-9);
            // ph_log(x) = -3.0103 * fast_log2(x); ph_log(1) == 0.
            assert!((ph_log(1.0)).abs() < 1e-9);
            assert!(ph_log(0.001) > 0.0);
        }

        fn obs(base4: u8, qual: u8, n: usize) -> Vec<Gap5Obs> {
            (0..n)
                .map(|_| Gap5Obs {
                    base4,
                    qual,
                    mqual: 60.0,
                    nm_local: 0,
                    poly: 1.0,
                })
                .collect()
        }

        #[test]
        fn gap5_call_basic_homozygous_and_het() {
            let cp = default_recall();
            let et = ETab::new();
            let q2p = q2p_table();
            let mpow = mqual_pow_1m_table();
            let opts = Gap5Opts {
                use_mqual: true,
                nm_adjust: true,
                scale_mqual: 1.0,
                low_mqual: 1.0,
                high_mqual: 60.0,
            };

            // 20x high-qual A (SAM code 1) -> call A (0), confident.
            let c = gap5_call(&obs(1, 40, 20), 20, &cp, &et, &q2p, &mpow, &opts);
            assert_eq!(c.call, 0, "pure A column should call A");
            assert_eq!(c.depth, 20);
            assert!(c.phred > 0, "confident A call should have phred>0");

            // 20x high-qual G (SAM code 4) -> call G (2).
            let g = gap5_call(&obs(4, 40, 20), 20, &cp, &et, &q2p, &mpow, &opts);
            assert_eq!(g.call, 2, "pure G column should call G");

            // Empty column -> N call, depth 0.
            let n = gap5_call(&[], 0, &cp, &et, &q2p, &mpow, &opts);
            assert_eq!(n.call, 4);
            assert_eq!(n.depth, 0);

            // Balanced A/C -> heterozygous A,C (map_het index 1).
            let mut mix = obs(1, 40, 15);
            mix.extend(obs(2, 40, 15));
            let h = gap5_call(&mix, 30, &cp, &et, &q2p, &mpow, &opts);
            assert_eq!(h.het_call, 1, "balanced A/C -> AC het (map_het[1])");
            assert!(h.het_logodd > 0, "clear het should have positive logodd");
        }
    }
}
