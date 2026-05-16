//! `samtools phase` — heterozygote phasing.
//!
//! This ports the standalone algorithm in upstream `phase.c`: discover
//! heterozygous markers from pileup genotype likelihoods, phase adjacent
//! read-backed markers with the local haplotype dynamic program, mask suspect
//! regions, and optionally split reads into haplotype/chimera BAMs. Upstream
//! ships no dedicated `test_phase` fixtures, so coverage is focused unit tests.

use std::collections::{HashMap, HashSet};
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::alignment_compat::{
    PileupColumn, PileupOptions, PileupRead, pileup_from_alignment_paths_with_options,
    pileup_from_alignment_paths_with_reference_and_options,
};
use htslib_rs::format::Exact;
use htslib_rs::math::kf_lgamma;

use crate::bam_flag::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FUNMAP};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io::sam_open_format;

const MAX_VARS: usize = 256;
const FLIP_PENALTY: i32 = 2;
const FLIP_THRES: i32 = 4;
const MASK_THRES: i32 = 3;
const FLAG_FIX_CHIMERA: u8 = 0x1;
const FLAG_LIST_EXCL: u8 = 0x4;
const FLAG_DROP_AMBIG: u8 = 0x8;
const ERR_DEP: f64 = 0.83;

#[derive(Clone, Debug)]
struct Config {
    input: PathBuf,
    reference: Option<PathBuf>,
    output_prefix: Option<String>,
    site_list: Option<PathBuf>,
    flags: u8,
    k: usize,
    min_base_q: u8,
    min_var_lod: i32,
    max_depth: usize,
    no_pg: bool,
    argv: Vec<OsString>,
}

impl Config {
    fn default_with_argv(argv: &[OsString]) -> Self {
        Self {
            input: PathBuf::new(),
            reference: crate::sam_global::current_global_args().reference,
            output_prefix: None,
            site_list: None,
            flags: FLAG_FIX_CHIMERA,
            k: 13,
            min_base_q: 13,
            min_var_lod: 37,
            max_depth: 256,
            no_pg: false,
            argv: argv.to_vec(),
        }
    }
}

/// Entry point for `samtools phase`.
pub fn main(args: &[OsString]) -> ExitCode {
    let cfg = match parse_args(args) {
        Ok(cfg) => cfg,
        Err(ParseError::Usage) => {
            let _ = usage(&mut io::stderr().lock());
            return ExitCode::from(1);
        }
        Err(ParseError::Message(msg)) => {
            print_error("phase", msg);
            let _ = usage(&mut io::stderr().lock());
            return ExitCode::from(1);
        }
    };

    let mut report = Vec::new();
    match run(&cfg, &mut report) {
        Ok(result) => {
            if let Err(e) = io::stdout().lock().write_all(&report) {
                print_error_errno("phase", "write phase report", &e);
                return ExitCode::from(1);
            }
            if let Some(prefix) = cfg.output_prefix.as_deref()
                && let Err(e) = write_split_bams(&cfg, prefix, &result.assignments)
            {
                print_error_errno("phase", "write split BAMs", &e);
                return ExitCode::from(1);
            }
            ExitCode::SUCCESS
        }
        Err(e) => {
            print_error_errno("phase", "phasing failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Debug)]
enum ParseError {
    Usage,
    Message(String),
}

fn normalize_args(args: &[OsString]) -> Vec<OsString> {
    let mut out = Vec::with_capacity(args.len().saturating_sub(1));
    for arg in args.iter().skip(1) {
        let s = arg.to_string_lossy();
        if let Some(rest) = s.strip_prefix("--") {
            if let Some((name, val)) = rest.split_once('=') {
                out.push(OsString::from(format!("--{name}")));
                out.push(OsString::from(val));
            } else {
                out.push(arg.clone());
            }
            continue;
        }

        if s.len() > 2 && s.starts_with('-') {
            let opt = s.as_bytes()[1];
            if matches!(opt, b'Q' | b'q' | b'k' | b'b' | b'l' | b'D' | b'f') {
                out.push(OsString::from(format!("-{}", opt as char)));
                out.push(OsString::from(&s[2..]));
                continue;
            }
        }

        out.push(arg.clone());
    }
    out
}

fn parse_args(args: &[OsString]) -> Result<Config, ParseError> {
    let mut cfg = Config::default_with_argv(args);
    let mut input = None;

    let normalized = normalize_args(args);
    let mut iter = normalized.iter();
    while let Some(arg) = iter.next() {
        let s = arg.to_string_lossy();
        match s.as_ref() {
            "-h" | "--help" => return Err(ParseError::Usage),
            "-D" => cfg.max_depth = parse_next(&mut iter, "-D")?.parse_num("-D")?,
            "-q" => cfg.min_var_lod = parse_next(&mut iter, "-q")?.parse_num("-q")?,
            "-Q" | "--min-BQ" | "--min-bq" => {
                cfg.min_base_q = parse_next(&mut iter, s.as_ref())?.parse_num(s.as_ref())?;
            }
            "-k" => cfg.k = parse_next(&mut iter, "-k")?.parse_num("-k")?,
            "-F" => cfg.flags &= !FLAG_FIX_CHIMERA,
            "-e" => cfg.flags |= FLAG_LIST_EXCL,
            "-A" => cfg.flags |= FLAG_DROP_AMBIG,
            "-b" => cfg.output_prefix = Some(parse_next(&mut iter, "-b")?),
            "-l" => cfg.site_list = Some(PathBuf::from(parse_next(&mut iter, "-l")?)),
            "-f" | "--reference" | "--fasta-ref" => {
                cfg.reference = Some(PathBuf::from(parse_next(&mut iter, s.as_ref())?));
            }
            "--no-PG" => cfg.no_pg = true,
            _ if s.starts_with('-') => {
                return Err(ParseError::Message(format!("unknown option {s}")));
            }
            _ => {
                if input.is_some() {
                    return Err(ParseError::Message(
                        "multiple input files are not supported".into(),
                    ));
                }
                input = Some(PathBuf::from(arg));
            }
        }
    }

    cfg.input = input.ok_or(ParseError::Usage)?;
    if cfg.k == 0 || cfg.k > 31 {
        return Err(ParseError::Message("-k must be between 1 and 31".into()));
    }
    if cfg.max_depth == 0 {
        return Err(ParseError::Message("-D must be greater than zero".into()));
    }
    if cfg.site_list.is_none() {
        cfg.flags &= !FLAG_LIST_EXCL;
    }

    Ok(cfg)
}

trait ParseNum {
    fn parse_num<T: std::str::FromStr>(self, option: &str) -> Result<T, ParseError>;
}

impl ParseNum for String {
    fn parse_num<T: std::str::FromStr>(self, option: &str) -> Result<T, ParseError> {
        self.parse()
            .map_err(|_| ParseError::Message(format!("invalid {option} value")))
    }
}

fn parse_next<'a>(
    iter: &mut std::slice::Iter<'a, OsString>,
    option: &str,
) -> Result<String, ParseError> {
    iter.next()
        .and_then(|a| a.to_str())
        .map(str::to_owned)
        .ok_or_else(|| ParseError::Message(format!("option {option} requires an argument")))
}

fn usage(mut w: impl Write) -> io::Result<()> {
    writeln!(w)?;
    writeln!(w, "Usage:   samtools phase [options] <in.bam>")?;
    writeln!(w)?;
    writeln!(w, "Options: -k INT    block length [13]")?;
    writeln!(w, "         -b STR    prefix of BAMs to output [null]")?;
    writeln!(w, "         -q INT    min het phred-LOD [37]")?;
    writeln!(w, "         -Q, --min-BQ INT")?;
    writeln!(w, "                   min base quality in het calling [13]")?;
    writeln!(w, "         -D INT    max read depth [256]")?;
    writeln!(w, "         -F        do not attempt to fix chimeras")?;
    writeln!(w, "         -A        drop reads with ambiguous phase")?;
    writeln!(w, "         --no-PG   do not add a PG line")?;
    writeln!(w)
}

#[derive(Clone, Debug, Default)]
struct PhaseRun {
    assignments: Vec<Assignment>,
}

fn run(cfg: &Config, out: &mut dyn Write) -> io::Result<PhaseRun> {
    let site_set = match cfg.site_list.as_ref() {
        Some(path) => load_positions(path)?,
        None => HashSet::new(),
    };
    let columns = pileup_columns(cfg)?;
    let errmod = ErrMod::new(1.0 - ERR_DEP)?;
    let mut state = PhaseState::new(cfg, site_set);

    write_report_header(&mut *out)?;

    for column in &columns {
        if state.last_ref.as_deref() != Some(column.reference_name.as_str()) {
            if let Some(ref_name) = state.last_ref.take() {
                state.phase_current_block(&ref_name, out)?;
                state.update_vpos(i32::MAX as usize);
            }
            state.last_ref = Some(column.reference_name.clone());
            state.vpos = 0;
            state.vpos_shift = 0;
            state.cns.clear();
        }

        state.process_column(column, &errmod, out)?;
    }

    if let Some(ref_name) = state.last_ref.take() {
        state.phase_current_block(&ref_name, out)?;
    }

    Ok(PhaseRun {
        assignments: state.assignments,
    })
}

fn pileup_columns(cfg: &Config) -> io::Result<Vec<PileupColumn>> {
    let options = PileupOptions {
        exclude_flags: (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) as u16,
        detect_overlaps: false,
        discard_orphans: false,
        ..Default::default()
    };

    if let Some(reference) = cfg.reference.as_ref() {
        pileup_from_alignment_paths_with_reference_and_options(
            std::slice::from_ref(&cfg.input),
            reference,
            &options,
        )
    } else {
        let format = sam_open_format(&cfg.input)?;
        if format.exact == Exact::Cram {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "CRAM input requires -f/--reference FILE",
            ));
        }
        pileup_from_alignment_paths_with_options(std::slice::from_ref(&cfg.input), &options)
    }
}

fn write_report_header(mut out: impl Write) -> io::Result<()> {
    writeln!(out, "CC")?;
    writeln!(out, "CC\tDescriptions:\nCC")?;
    writeln!(out, "CC\t  CC      comments")?;
    writeln!(out, "CC\t  PS      start of a phase set")?;
    writeln!(out, "CC\t  FL      filtered region")?;
    writeln!(
        out,
        "CC\t  M[012]  markers; 0 for singletons, 1 for phased and 2 for filtered"
    )?;
    writeln!(out, "CC\t  EV      supporting reads; SAM format")?;
    writeln!(out, "CC\t  //      end of a phase set\nCC")?;
    writeln!(
        out,
        "CC\tFormats of PS, FL and M[012] lines (1-based coordinates):\nCC"
    )?;
    writeln!(out, "CC\t  PS  chr  phaseSetStart  phaseSetEnd")?;
    writeln!(out, "CC\t  FL  chr  filterStart    filterEnd")?;
    writeln!(
        out,
        "CC\t  M?  chr  PS  pos  allele0  allele1  hetIndex  #supports0  #errors0  #supp1  #err1"
    )?;
    writeln!(out, "CC\nCC")
}

#[derive(Clone, Debug)]
struct PhaseState<'a> {
    cfg: &'a Config,
    site_set: HashSet<(String, usize)>,
    last_ref: Option<String>,
    vpos: usize,
    vpos_shift: usize,
    cns: Vec<u64>,
    frags: HashMap<u64, Fragment>,
    assignments: Vec<Assignment>,
}

impl<'a> PhaseState<'a> {
    fn new(cfg: &'a Config, site_set: HashSet<(String, usize)>) -> Self {
        Self {
            cfg,
            site_set,
            last_ref: None,
            vpos: 0,
            vpos_shift: 0,
            cns: Vec::new(),
            frags: HashMap::new(),
            assignments: Vec::new(),
        }
    }

    fn process_column(
        &mut self,
        column: &PileupColumn,
        errmod: &ErrMod,
        out: &mut dyn Write,
    ) -> io::Result<()> {
        let reads = column
            .reads_by_input
            .first()
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        if reads.len() > self.cfg.max_depth {
            return Ok(());
        }

        let mut bases = genotype_bases(reads, self.cfg.min_base_q);
        if bases.is_empty() {
            return Ok(());
        }

        let q = errmod.cal(&mut bases, 4);
        let cns = gl2cns(&q);
        let pos0 = column.position.saturating_sub(1);
        let in_set = self
            .site_set
            .contains(&(column.reference_name.clone(), pos0));
        if !in_set && (self.cfg.flags & FLAG_LIST_EXCL) != 0 {
            return Ok(());
        }
        if !in_set && consensus_lod(cns) < self.cfg.min_var_lod {
            return Ok(());
        }

        if self.vpos == self.cns.len() {
            self.cns.push(0);
        }
        self.cns[self.vpos] = ((pos0 as u64) << 32) | u64::from(cns);

        let mut dophase = true;
        for read in reads {
            if read.is_deletion || read.is_refskip || read.mapping_quality == 0 {
                continue;
            }
            let Some(name) = read.name.as_deref() else {
                continue;
            };
            let observed = read.base.and_then(base_index);
            let code = match observed {
                Some(base) if u32::from(base) == (cns & 3) => 1,
                Some(base) if u32::from(base) == ((cns >> 16) & 3) => 2,
                _ => 0,
            };
            let key = x31_hash_string(name);
            if let Some(frag) = self.frags.get_mut(&key) {
                let len = self.vpos.saturating_sub(frag.vpos) + 1;
                if len < MAX_VARS {
                    frag.vlen = len;
                    frag.seq[len - 1] = code;
                    frag.end = pos0 + 1;
                }
                dophase = false;
            } else {
                let mut seq = [0u8; MAX_VARS];
                seq[0] = code;
                self.frags.insert(
                    key,
                    Fragment {
                        key,
                        seq,
                        vpos: self.vpos,
                        beg: pos0,
                        end: pos0 + 1,
                        vlen: 1,
                        single: false,
                        flip: false,
                        phase: false,
                        phased: false,
                        ambig: false,
                        in_count: 0,
                        out_count: 0,
                    },
                );
            }
        }

        if dophase {
            let ref_name = column.reference_name.as_str();
            let next_cns = self.cns[self.vpos];
            self.phase_current_block(ref_name, out)?;
            self.update_vpos(self.vpos);
            self.cns.clear();
            self.cns.push(next_cns);
            self.vpos = 0;
        }
        self.vpos += 1;

        Ok(())
    }

    fn phase_current_block(&mut self, reference_name: &str, out: &mut dyn Write) -> io::Result<()> {
        if self.vpos == 0 {
            return Ok(());
        }

        self.clean_seqs();
        if self.vpos == 1 {
            let c = self.cns[0];
            writeln!(out, "PS\t{}\t{}\t{}", reference_name, pos1(c), pos1(c))?;
            writeln!(
                out,
                "M0\t{}\t{}\t{}\t{}\t{}\t{}\t0\t0\t0\t0",
                reference_name,
                pos1(c),
                pos1(c),
                allele_char((c & 3) as u8),
                allele_char(((c >> 16) & 3) as u8),
                self.vpos_shift + 1
            )?;
            writeln!(out, "//")?;

            for frag in self.frags.values_mut() {
                if frag.vpos == 0 {
                    frag.flip = false;
                    if frag.seq[0] == 0 {
                        frag.phased = false;
                    } else {
                        frag.phased = true;
                        frag.phase = frag.seq[0] == 2;
                    }
                    self.assignments
                        .push(Assignment::from_fragment(reference_name, frag));
                }
            }
            self.vpos_shift += 1;
            return Ok(());
        }

        writeln!(
            out,
            "PS\t{}\t{}\t{}",
            reference_name,
            pos1(self.cns[0]),
            pos1(self.cns[self.vpos - 1])
        )?;

        let counts = count_all(self.cfg.k, self.vpos, &mut self.frags);
        let path = dynaprog(self.cfg.k, self.vpos, &counts);
        let mut pcnt = fragphase(self.vpos, &path, &mut self.frags, false);
        let masks = genmask(self.vpos, &pcnt);
        let mut site_mask = vec![false; self.vpos];
        for &(beg, end) in &masks {
            writeln!(
                out,
                "FL\t{}\t{}\t{}",
                reference_name,
                pos1(self.cns[beg]),
                pos1(self.cns[end])
            )?;
            for masked in site_mask.iter_mut().take(end + 1).skip(beg) {
                *masked = true;
            }
        }
        if (self.cfg.flags & FLAG_FIX_CHIMERA) != 0 {
            pcnt = fragphase(self.vpos, &path, &mut self.frags, true);
        }

        for i in 0..self.vpos {
            let x = pcnt[i];
            let c0 = display_allele_low(self.cns[i]);
            let c1 = display_allele_high(self.cns[i]);
            let alleles = [c0, c1];
            writeln!(
                out,
                "M{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                if site_mask[i] { 2 } else { 1 },
                reference_name,
                pos1(self.cns[0]),
                pos1(self.cns[i]),
                allele_char(alleles[path[i] as usize]),
                allele_char(alleles[1 - path[i] as usize]),
                i + self.vpos_shift + 1,
                x & 0xffff,
                (x >> 16) & 0xffff,
                (x >> 32) & 0xffff,
                (x >> 48) & 0xffff
            )?;
        }

        let mut seqs: Vec<_> = self
            .frags
            .values()
            .filter(|f| f.vpos < self.vpos && !f.single)
            .cloned()
            .collect();
        seqs.sort_by_key(|f| f.vpos);
        for frag in &seqs {
            write!(
                out,
                "EV\t0\t{}\t{}\t40\t{}M\t*\t0\t0\t",
                reference_name,
                frag.vpos + 1 + self.vpos_shift,
                frag.vlen
            )?;
            for j in 0..frag.vlen {
                let c = self.cns[frag.vpos + j];
                let base = match frag.seq[j] {
                    0 => b'N',
                    1 => b"ACGT"[(c & 3) as usize],
                    2 => b"ACGT"[((c >> 16) & 3) as usize],
                    _ => b'N',
                };
                out.write_all(&[base])?;
            }
            writeln!(
                out,
                "\t*\tYP:i:{}\tYF:i:{}\tYI:i:{}\tYO:i:{}\tYS:i:{}",
                u8::from(frag.phase),
                u8::from(frag.flip),
                frag.in_count,
                frag.out_count,
                frag.beg + 1
            )?;
        }
        writeln!(out, "//")?;

        for frag in self.frags.values() {
            if frag.vpos < self.vpos {
                self.assignments
                    .push(Assignment::from_fragment(reference_name, frag));
            }
        }
        self.vpos_shift += self.vpos;
        Ok(())
    }

    fn clean_seqs(&mut self) -> bool {
        let mut had_future = false;
        self.frags.retain(|_, frag| {
            if frag.vpos >= self.vpos {
                had_future = true;
                return true;
            }
            let beg = frag.seq[..frag.vlen]
                .iter()
                .position(|&c| c != 0)
                .unwrap_or(frag.vlen);
            let end = frag.seq[..frag.vlen]
                .iter()
                .rposition(|&c| c != 0)
                .map(|i| i + 1)
                .unwrap_or(0);
            if end <= beg {
                return false;
            }
            if beg != 0 {
                frag.seq.copy_within(beg..end, 0);
            }
            frag.vpos += beg;
            frag.vlen = end - beg;
            frag.single = frag.vlen == 1;
            true
        });
        had_future
    }

    fn update_vpos(&mut self, vpos: usize) {
        self.frags.retain(|_, frag| {
            if frag.vpos < vpos {
                false
            } else {
                frag.vpos -= vpos;
                true
            }
        });
        if vpos < self.cns.len() {
            self.cns.drain(0..vpos);
        } else {
            self.cns.clear();
        }
    }
}

#[derive(Clone, Debug)]
struct Fragment {
    key: u64,
    seq: [u8; MAX_VARS],
    vpos: usize,
    beg: usize,
    end: usize,
    vlen: usize,
    single: bool,
    flip: bool,
    phase: bool,
    phased: bool,
    ambig: bool,
    in_count: u16,
    out_count: u16,
}

#[derive(Clone, Debug)]
struct Assignment {
    key: u64,
    reference_name: String,
    beg: usize,
    end: usize,
    phase: bool,
    phased: bool,
    flip: bool,
    ambig: bool,
}

impl Assignment {
    fn from_fragment(reference_name: &str, frag: &Fragment) -> Self {
        Self {
            key: frag.key,
            reference_name: reference_name.to_owned(),
            beg: frag.beg,
            end: frag.end,
            phase: frag.phase,
            phased: frag.phased,
            flip: frag.flip,
            ambig: frag.ambig,
        }
    }
}

fn genotype_bases(reads: &[PileupRead], min_base_q: u8) -> Vec<u16> {
    let mut bases = Vec::new();
    for read in reads {
        if read.is_refskip || read.is_deletion {
            continue;
        }
        let base_q = read.qpos_quality;
        if base_q < min_base_q {
            continue;
        }
        let Some(base) = read.base.and_then(base_index) else {
            continue;
        };
        let q = base_q.min(read.mapping_quality).clamp(4, 63);
        bases.push((u16::from(q) << 5) | (u16::from(read.is_reverse) << 4) | u16::from(base));
    }
    bases
}

fn gl2cns(q: &[f32; 16]) -> u32 {
    let mut min = 1e30f32;
    let mut min2 = 1e30f32;
    let mut min_ij = 0usize;
    for i in 0..4usize {
        for j in i..4usize {
            let v = q[i << 2 | j];
            if v < min {
                min_ij = i << 2 | j;
                min2 = min;
                min = v;
            } else if v < min2 {
                min2 = v;
            }
        }
    }
    if ((min_ij >> 2) & 3) == (min_ij & 3) {
        0
    } else {
        (1 << 18)
            | (((min_ij >> 2) as u32 & 3) << 16)
            | ((min_ij as u32) & 3)
            | (((min2 - min + 0.499) as u32) << 2)
    }
}

fn consensus_lod(cns: u32) -> i32 {
    ((cns & 0xffff) >> 2) as i32
}

fn count_all(l: usize, vpos: usize, frags: &mut HashMap<u64, Fragment>) -> Vec<Vec<i32>> {
    let cnt_len = 1usize << l;
    let mut counts = vec![vec![0i32; cnt_len]; vpos];
    let mut seq = vec![0u8; l];
    for frag in frags.values_mut() {
        if frag.vpos >= vpos || frag.single {
            continue;
        }
        if frag.vlen == 1 {
            frag.single = true;
            continue;
        }
        for j in 1..frag.vlen {
            for (i, slot) in seq.iter_mut().enumerate() {
                *slot = if j < l - 1 - i {
                    0
                } else {
                    frag.seq[j - (l - 1 - i)]
                };
            }
            count1(l, &seq, &mut counts[frag.vpos + j]);
        }
    }
    counts
}

fn count1(l: usize, seq: &[u8], cnt: &mut [i32]) {
    if seq[l - 1] == 0 {
        return;
    }
    let n_ambi = seq.iter().take(l).filter(|&&c| c == 0).count();
    if l - n_ambi <= 1 {
        return;
    }
    for x in 0..(1usize << n_ambi) {
        let mut z = 0usize;
        let mut j = 0usize;
        for &base in seq.iter().take(l) {
            let c = if base != 0 {
                usize::from(base - 1)
            } else {
                let c = (x >> j) & 1;
                j += 1;
                c
            };
            z = (z << 1) | c;
        }
        cnt[z] += 1;
    }
}

fn dynaprog(l: usize, vpos: usize, w: &[Vec<i32>]) -> Vec<u8> {
    let z = 1usize << (l - 1);
    let mask = (1usize << l) - 1;
    let mut prev = vec![0i32; z];
    let mut curr = vec![0i32; z];
    let mut back = vec![vec![0u8; z]; vpos];

    for i in 0..vpos {
        for x in 0..z {
            let xc = (!x) & mask;
            let y0 = x >> 1;
            let y1 = xc >> 1;
            let c0 = prev[y0] + w[i][x] + w[i][xc];
            let c1 = prev[y1] + w[i][x] + w[i][xc];
            if c0 > c1 {
                back[i][x] = 0;
                curr[x] = c0;
            } else {
                back[i][x] = 1;
                curr[x] = c1;
            }
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    let mut h = vec![0u8; vpos];
    let mut x = prev
        .iter()
        .enumerate()
        .max_by_key(|(_, v)| *v)
        .map(|(i, _)| i)
        .unwrap_or(0);
    let mut which = false;
    for i in (0..vpos).rev() {
        h[i] = if which {
            ((!x) & 1) as u8
        } else {
            (x & 1) as u8
        };
        if back[i][x] != 0 {
            which = !which;
            x = ((!x) & mask) >> 1;
        } else {
            x >>= 1;
        }
    }
    h
}

fn fragphase(
    vpos: usize,
    path: &[u8],
    frags: &mut HashMap<u64, Fragment>,
    fix_flip: bool,
) -> Vec<u64> {
    let mut pcnt = vec![0u64; vpos];
    for frag in frags.values_mut() {
        if frag.vpos >= vpos {
            continue;
        }
        let mut c = [0u16; 2];
        for i in 0..frag.vlen {
            if frag.seq[i] == 0 {
                continue;
            }
            let idx = usize::from(frag.seq[i] != path[frag.vpos + i] + 1);
            c[idx] += 1;
        }
        frag.phase = c[0] <= c[1];
        frag.in_count = c[usize::from(frag.phase)];
        frag.out_count = c[1 - usize::from(frag.phase)];
        frag.phased = frag.in_count != frag.out_count;
        frag.ambig = frag.in_count != 0
            && frag.out_count != 0
            && frag.out_count < 3
            && frag.in_count <= frag.out_count + 1;
        frag.flip = false;

        if fix_flip && c[0] >= 3 && c[1] >= 3 {
            maybe_fix_chimera(frag, path);
        }

        if !frag.single {
            for i in 0..frag.vlen {
                if frag.seq[i] == 0 {
                    continue;
                }
                let base = if frag.phase {
                    2 - frag.seq[i]
                } else {
                    frag.seq[i] - 1
                };
                let idx = frag.vpos + i;
                if base == path[idx] {
                    if !frag.phase {
                        pcnt[idx] += 1;
                    } else {
                        pcnt[idx] += 1u64 << 32;
                    }
                } else if !frag.phase {
                    pcnt[idx] += 1u64 << 16;
                } else {
                    pcnt[idx] += 1u64 << 48;
                }
            }
        }
    }
    pcnt
}

fn maybe_fix_chimera(frag: &mut Fragment, path: &[u8]) {
    let mut left = vec![0u32; frag.vlen];
    let mut right = vec![0u32; frag.vlen];
    let mut sum = [0u16; 2];
    for (i, slot) in left.iter_mut().enumerate().take(frag.vlen) {
        if frag.seq[i] != 0 {
            let c = if frag.phase {
                2 - frag.seq[i]
            } else {
                frag.seq[i] - 1
            };
            sum[usize::from(c != path[frag.vpos + i])] += 1;
        }
        *slot = (u32::from(sum[1]) << 16) | u32::from(sum[0]);
    }
    sum = [0; 2];
    for i in (0..frag.vlen).rev() {
        if frag.seq[i] != 0 {
            let c = if frag.phase {
                2 - frag.seq[i]
            } else {
                frag.seq[i] - 1
            };
            sum[usize::from(c != path[frag.vpos + i])] += 1;
        }
        right[i] = (u32::from(sum[1]) << 16) | u32::from(sum[0]);
    }

    let mut best = 0;
    let mut best_i = None;
    let mut best_dir = 0;
    for i in 0..frag.vlen.saturating_sub(1) {
        let a0 = low16(left[i]) + high16(right[i + 1]) - low16(right[i + 1]) * FLIP_PENALTY;
        let a1 = high16(left[i]) + low16(right[i + 1]) - high16(right[i + 1]) * FLIP_PENALTY;
        if a0 > a1 {
            if a0 > best {
                best = a0;
                best_i = Some(i);
                best_dir = 0;
            }
        } else if a1 > best {
            best = a1;
            best_i = Some(i);
            best_dir = 1;
        }
    }

    if let Some(i) = best_i
        && best - i32::from(frag.in_count) >= FLIP_THRES
        && best - i32::from(frag.out_count) >= FLIP_THRES
    {
        frag.flip = true;
        let range: Box<dyn Iterator<Item = usize>> = if best_dir == 0 {
            Box::new((i + 1)..frag.vlen)
        } else {
            Box::new(0..=i)
        };
        for j in range {
            if frag.seq[j] == 1 {
                frag.seq[j] = 2;
            } else if frag.seq[j] == 2 {
                frag.seq[j] = 1;
            }
        }
    }
}

fn low16(v: u32) -> i32 {
    (v & 0xffff) as i32
}

fn high16(v: u32) -> i32 {
    ((v >> 16) & 0xffff) as i32
}

fn genmask(vpos: usize, pcnt: &[u64]) -> Vec<(usize, usize)> {
    let mut max = 0;
    let mut max_i = 0usize;
    let mut beg = 0usize;
    let mut score = 0;
    let mut list = Vec::new();
    let mut i = 0usize;
    while i < vpos {
        let x = pcnt[i];
        let c0 = (x & 0xffff) as i32;
        let c1 = ((x >> 16) & 0xffff) as i32;
        let c2 = ((x >> 32) & 0xffff) as i32;
        let c3 = ((x >> 48) & 0xffff) as i32;
        let pre = score;
        let mut s = if c1 + c3 == 0 {
            -(c0 + c2)
        } else {
            c1 + c3 - 1
        };
        if c3 > c2 {
            s += c3 - c2;
        }
        if c1 > c0 {
            s += c1 - c0;
        }
        score += s;
        if score < 0 {
            score = 0;
        }
        if pre == 0 && score > 0 {
            beg = i;
        }
        if (i == vpos - 1 || score == 0) && max >= MASK_THRES {
            list.push((beg, max_i));
            i = max_i;
            score = 0;
        } else if score > max {
            max = score;
            max_i = i;
        }
        if score == 0 {
            max = 0;
        }
        i += 1;
    }
    list
}

fn write_split_bams(cfg: &Config, prefix: &str, assignments: &[Assignment]) -> io::Result<()> {
    let text = input_as_sam_text(cfg)?;
    let (header, records) = split_sam_text(&text);
    let header = if cfg.no_pg {
        header
    } else {
        crate::pg::add_samtools_pg(&header, &cfg.argv)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
    };

    let mut outs = [header.clone(), header.clone(), header];
    for line in records {
        let mut record = SamRecordLine::parse(line)?;
        let which = classify_record(&record, assignments, (cfg.flags & FLAG_DROP_AMBIG) != 0);
        if which.add_zp {
            record.append_aux("ZP:A:Y");
        }
        outs[which.index].push_str(record.line.as_ref());
        outs[which.index].push('\n');
    }

    for (middle, text) in [("0", &outs[0]), ("1", &outs[1]), ("chimera", &outs[2])] {
        let path = format!("{prefix}.{middle}.bam");
        let mut file = File::create(path)?;
        htslib_rs::alignment_compat::write_bam_from_sam_reader(
            Cursor::new(text.as_bytes()),
            &mut file,
        )?;
    }
    Ok(())
}

fn input_as_sam_text(cfg: &Config) -> io::Result<String> {
    match sam_open_format(&cfg.input)?.exact {
        Exact::Sam => {
            let mut s = String::new();
            File::open(&cfg.input)?.read_to_string(&mut s)?;
            Ok(s)
        }
        Exact::Bam => {
            htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&cfg.input, None)
        }
        Exact::Cram => {
            let reference = cfg.reference.as_ref().ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "CRAM input requires -f/--reference FILE",
                )
            })?;
            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                &cfg.input, reference, None,
            )
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "phase input must be SAM, BAM, or CRAM",
        )),
    }
}

fn split_sam_text(text: &str) -> (String, Vec<&str>) {
    let mut header = String::new();
    let mut records = Vec::new();
    for line in text.lines() {
        if line.starts_with('@') {
            header.push_str(line);
            header.push('\n');
        } else if !line.is_empty() {
            records.push(line);
        }
    }
    (header, records)
}

#[derive(Debug)]
struct SamRecordLine<'a> {
    line: String,
    qname: &'a str,
    flag: u16,
    rname: &'a str,
    start0: usize,
    end0: usize,
}

impl<'a> SamRecordLine<'a> {
    fn parse(line: &'a str) -> io::Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "malformed SAM record",
            ));
        }
        let flag = fields[1].parse().unwrap_or(0);
        let pos1 = fields[3].parse::<usize>().unwrap_or(0);
        let start0 = pos1.saturating_sub(1);
        let end0 = start0 + cigar_ref_len(fields[5]);
        Ok(Self {
            line: line.to_owned(),
            qname: fields[0],
            flag,
            rname: fields[2],
            start0,
            end0,
        })
    }

    fn append_aux(&mut self, aux: &str) {
        self.line.push('\t');
        self.line.push_str(aux);
    }
}

#[derive(Clone, Copy)]
struct SplitClass {
    index: usize,
    add_zp: bool,
}

fn classify_record(
    record: &SamRecordLine<'_>,
    assignments: &[Assignment],
    drop_ambig: bool,
) -> SplitClass {
    if record.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) as u16 != 0 {
        return SplitClass {
            index: hash_to_bucket(record.qname),
            add_zp: false,
        };
    }
    let key = x31_hash_string(record.qname.as_bytes());
    let assignment = assignments.iter().find(|assignment| {
        assignment.key == key
            && assignment.reference_name == record.rname
            && record.end0 > assignment.beg
            && record.start0 < assignment.end
    });
    match assignment {
        Some(a) if a.ambig && drop_ambig => SplitClass {
            index: 2,
            add_zp: false,
        },
        Some(a) if a.phased && a.flip => SplitClass {
            index: 2,
            add_zp: false,
        },
        Some(a) if a.phased => SplitClass {
            index: usize::from(a.phase),
            add_zp: true,
        },
        _ => SplitClass {
            index: hash_to_bucket(record.qname),
            add_zp: false,
        },
    }
}

fn hash_to_bucket(qname: &str) -> usize {
    (x31_hash_string(qname.as_bytes()) & 1) as usize
}

fn cigar_ref_len(cigar: &str) -> usize {
    let mut n = 0usize;
    let mut len = 0usize;
    for b in cigar.bytes() {
        if b.is_ascii_digit() {
            n = n * 10 + usize::from(b - b'0');
        } else {
            if matches!(b, b'M' | b'D' | b'N' | b'=' | b'X') {
                len += n;
            }
            n = 0;
        }
    }
    len
}

fn load_positions(path: &Path) -> io::Result<HashSet<(String, usize)>> {
    let text = std::fs::read_to_string(path)?;
    let mut set = HashSet::new();
    for line in text.lines() {
        let mut fields = line.split_whitespace();
        let Some(reference) = fields.next() else {
            continue;
        };
        let Some(pos) = fields.next() else {
            continue;
        };
        if let Ok(pos1) = pos.parse::<usize>() {
            set.insert((reference.to_owned(), pos1.saturating_sub(1)));
        }
    }
    Ok(set)
}

fn pos1(cns: u64) -> usize {
    ((cns >> 32) as usize) + 1
}

fn display_allele_low(cns: u64) -> u8 {
    if ((cns & 0xffff) >> 2) == 0 {
        4
    } else {
        (cns & 3) as u8
    }
}

fn display_allele_high(cns: u64) -> u8 {
    if (((cns >> 16) & 0xffff) >> 2) == 0 {
        4
    } else {
        ((cns >> 16) & 3) as u8
    }
}

fn allele_char(base: u8) -> char {
    match base {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => 'X',
    }
}

fn base_index(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn x31_hash_string(name: &[u8]) -> u64 {
    let Some((&first, rest)) = name.split_first() else {
        return 0;
    };
    let mut h = u64::from(first);
    for &b in rest {
        h = (h << 5).wrapping_sub(h).wrapping_add(u64::from(b));
    }
    h
}

#[derive(Clone, Debug)]
struct ErrMod {
    fk: Vec<f64>,
    beta: Vec<f64>,
    lhet: Vec<f64>,
}

impl ErrMod {
    fn new(depcorr: f64) -> io::Result<Self> {
        let eta = 0.03;
        let mut fk = vec![0.0; 256];
        fk[0] = 1.0;
        for (n, v) in fk.iter_mut().enumerate().skip(1) {
            *v = (1.0 - depcorr).powi(n as i32) * (1.0 - eta) + eta;
        }

        let logbinom = logbinomial_table();
        let mut beta = vec![0.0; 256 * 256 * 64];
        for q in 1..64usize {
            let e = 10.0_f64.powf(-(q as f64) / 10.0);
            let le = e.ln();
            let le1 = (1.0 - e).ln();
            for n in 1..=255usize {
                let offset = q << 16 | n << 8;
                let mut sum1 = logbinom[n << 8 | n] + n as f64 * le;
                beta[offset | n] = f64::INFINITY;
                for k in (0..n).rev() {
                    let sum = sum1
                        + (logbinom[n << 8 | k] + k as f64 * le + (n - k) as f64 * le1 - sum1)
                            .exp()
                            .ln_1p();
                    beta[offset | k] = -10.0 / std::f64::consts::LN_10 * (sum1 - sum);
                    sum1 = sum;
                }
            }
        }

        let mut lhet = vec![0.0; 256 * 256];
        for n in 0..256usize {
            for k in 0..256usize {
                lhet[n << 8 | k] = logbinom[n << 8 | k] - std::f64::consts::LN_2 * n as f64;
            }
        }

        if beta.iter().any(|v| v.is_nan()) || lhet.iter().any(|v| v.is_nan()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "failed to initialise phase error model",
            ));
        }

        Ok(Self { fk, beta, lhet })
    }

    fn cal(&self, bases: &mut [u16], m: usize) -> [f32; 16] {
        let mut q = [0f32; 16];
        if bases.is_empty() {
            return q;
        }

        let n = if bases.len() > 255 {
            downsample_to_255(bases)
        } else {
            bases.len()
        };
        bases[..n].sort_unstable();

        let mut bsum = [0f64; 16];
        let mut counts = [0usize; 16];
        let mut strand_counts = [0usize; 32];

        for &packed in bases[..n].iter().rev() {
            let qual = ((packed >> 5) as usize).clamp(4, 63);
            let basestrand = (packed & 0x1f) as usize;
            let base = (packed & 0x0f) as usize;
            bsum[base] +=
                self.fk[strand_counts[basestrand]] * self.beta[qual << 16 | n << 8 | counts[base]];
            counts[base] += 1;
            strand_counts[basestrand] += 1;
        }

        for j in 0..m {
            let mut tmp1 = 0.0;
            let mut tmp2 = 0usize;
            for k in 0..m {
                if k != j {
                    tmp1 += bsum[k];
                    tmp2 += counts[k];
                }
            }
            if tmp2 != 0 {
                q[j * m + j] = tmp1 as f32;
            }

            for k in (j + 1)..m {
                let cjk = counts[j] + counts[k];
                tmp1 = 0.0;
                tmp2 = 0;
                for i in 0..m {
                    if i != j && i != k {
                        tmp1 += bsum[i];
                        tmp2 += counts[i];
                    }
                }
                let val = -4.343 * self.lhet[cjk << 8 | counts[k]] + tmp1;
                q[j * m + k] = val as f32;
                q[k * m + j] = q[j * m + k];
                if tmp2 == 0 {
                    q[j * m + k] = (-4.343 * self.lhet[cjk << 8 | counts[k]]) as f32;
                    q[k * m + j] = q[j * m + k];
                }
            }

            for k in 0..m {
                if q[j * m + k] < 0.0 {
                    q[j * m + k] = 0.0;
                }
            }
        }

        q
    }
}

fn logbinomial_table() -> Vec<f64> {
    let mut logbinom = vec![0.0; 256 * 256];
    for n in 1..256usize {
        let lfn = kf_lgamma(n as f64 + 1.0);
        for k in 1..=n {
            logbinom[n << 8 | k] =
                lfn - kf_lgamma(k as f64 + 1.0) - kf_lgamma((n - k) as f64 + 1.0);
        }
    }
    logbinom
}

fn downsample_to_255(bases: &mut [u16]) -> usize {
    let len = bases.len();
    for i in 0..255usize {
        let src = i * len / 255;
        bases[i] = bases[src];
    }
    255
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_dir(name: &str) -> PathBuf {
        let mut path = std::env::temp_dir();
        path.push(format!("samtools-rs-phase-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        path
    }

    fn write_phase_sam(path: &Path) {
        let sam = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:12
r0\t0\tchr1\t1\t60\t12M\t*\t0\t0\tAAAAAACCCCCC\tFFFFFFFFFFFF
r1\t0\tchr1\t1\t60\t12M\t*\t0\t0\tAAAAAACCCCCC\tFFFFFFFFFFFF
r2\t0\tchr1\t1\t60\t12M\t*\t0\t0\tCCCCCCAAAAAA\tFFFFFFFFFFFF
r3\t0\tchr1\t1\t60\t12M\t*\t0\t0\tCCCCCCAAAAAA\tFFFFFFFFFFFF
";
        std::fs::write(path, sam).unwrap();
    }

    fn test_cfg(input: &Path) -> Config {
        parse_args(&[
            OsString::from("phase"),
            OsString::from("-q"),
            OsString::from("1"),
            OsString::from("-Q"),
            OsString::from("1"),
            OsString::from("-k"),
            OsString::from("3"),
            input.as_os_str().to_os_string(),
        ])
        .unwrap()
    }

    #[test]
    fn phase_emits_phase_set_markers_and_evidence() {
        let dir = tmp_dir("report");
        let input = dir.join("in.sam");
        write_phase_sam(&input);
        let cfg = test_cfg(&input);

        let mut out = Vec::new();
        let result = run(&cfg, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();

        assert!(text.starts_with("CC\n"));
        assert!(text.contains("\nPS\tchr1\t"));
        assert!(text.contains("\nM1\tchr1\t"));
        assert!(text.contains("\nEV\t0\tchr1\t"));
        assert!(!result.assignments.is_empty());
    }

    #[test]
    fn phase_respects_min_base_quality() {
        let dir = tmp_dir("baseq");
        let input = dir.join("in.sam");
        write_phase_sam(&input);
        let mut cfg = test_cfg(&input);
        cfg.min_base_q = 50;

        let mut out = Vec::new();
        run(&cfg, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();

        assert!(!text.contains("\nPS\tchr1\t"));
    }

    #[test]
    fn phase_split_bams_are_written() {
        let dir = tmp_dir("split");
        let input = dir.join("in.sam");
        write_phase_sam(&input);
        let mut cfg = test_cfg(&input);
        let prefix = dir.join("phase-out");
        cfg.output_prefix = Some(prefix.to_string_lossy().into_owned());
        cfg.no_pg = true;

        let mut out = Vec::new();
        let result = run(&cfg, &mut out).unwrap();
        write_split_bams(
            &cfg,
            cfg.output_prefix.as_deref().unwrap(),
            &result.assignments,
        )
        .unwrap();

        for middle in ["0", "1", "chimera"] {
            let path = dir.join(format!("phase-out.{middle}.bam"));
            assert!(path.exists(), "{} missing", path.display());
        }
    }
}
