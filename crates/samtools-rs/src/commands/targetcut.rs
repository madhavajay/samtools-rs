//! `samtools targetcut` — fosmid pool target cutting.
//!
//! This is a faithful Rust port of the standalone algorithm in
//! `cut_target.c`: per-position pileup consensus is scored with HTSlib's
//! revised MAQ error model, then a two-state dynamic program emits long
//! target intervals as SAM-like records. Upstream ships no dedicated
//! `test_targetcut` fixtures, so this module is covered by focused unit tests.

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Write};
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

const ERR_DEP: f64 = 0.83;

#[derive(Clone, Debug)]
struct ScoreParam {
    e: [[i32; 3]; 2],
    p: [[i32; 2]; 2],
}

impl Default for ScoreParam {
    fn default() -> Self {
        Self {
            e: [[0, 0, 0], [-4, 1, 6]],
            p: [[0, -14000], [0, 0]],
        }
    }
}

#[derive(Clone, Debug)]
struct Config {
    input: PathBuf,
    reference: Option<PathBuf>,
    output: Option<PathBuf>,
    min_base_q: u8,
    score: ScoreParam,
}

/// Entry point for `samtools targetcut`.
pub fn main(args: &[OsString]) -> ExitCode {
    let cfg = match parse_args(args) {
        Ok(cfg) => cfg,
        Err(ParseError::Usage) => {
            let _ = usage(&mut io::stderr().lock());
            return ExitCode::from(1);
        }
        Err(ParseError::Message(msg)) => {
            print_error("targetcut", msg);
            let _ = usage(&mut io::stderr().lock());
            return ExitCode::from(1);
        }
    };

    let mut writer: Box<dyn Write> = match cfg.output.as_ref() {
        Some(path) => match File::create(path) {
            Ok(file) => Box::new(io::BufWriter::new(file)),
            Err(e) => {
                print_error_errno("targetcut", "open -o output", &e);
                return ExitCode::from(1);
            }
        },
        None => Box::new(io::BufWriter::new(io::stdout().lock())),
    };

    match run(&cfg, &mut writer) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("targetcut", "target cutting failed", &e);
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
            if matches!(opt, b'Q' | b'i' | b'0' | b'1' | b'2' | b'f' | b'o') {
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
    let mut input = None;
    let mut reference = crate::sam_global::current_global_args().reference;
    let mut output = None;
    let mut min_base_q = 13u8;
    let mut score = ScoreParam::default();

    let normalized = normalize_args(args);
    let mut iter = normalized.iter();
    while let Some(arg) = iter.next() {
        let s = arg.to_string_lossy();
        match s.as_ref() {
            "-h" | "--help" => return Err(ParseError::Usage),
            "-Q" => {
                min_base_q = parse_next(&mut iter, "-Q")?
                    .parse()
                    .map_err(|_| ParseError::Message("invalid -Q value".into()))?;
            }
            "-i" => {
                let v: i32 = parse_next(&mut iter, "-i")?
                    .parse()
                    .map_err(|_| ParseError::Message("invalid -i value".into()))?;
                score.p[0][1] = -v;
            }
            "-0" => {
                score.e[1][0] = parse_next(&mut iter, "-0")?
                    .parse()
                    .map_err(|_| ParseError::Message("invalid -0 value".into()))?;
            }
            "-1" => {
                score.e[1][1] = parse_next(&mut iter, "-1")?
                    .parse()
                    .map_err(|_| ParseError::Message("invalid -1 value".into()))?;
            }
            "-2" => {
                score.e[1][2] = parse_next(&mut iter, "-2")?
                    .parse()
                    .map_err(|_| ParseError::Message("invalid -2 value".into()))?;
            }
            "-f" | "--reference" | "--fasta-ref" => {
                reference = Some(PathBuf::from(parse_next(&mut iter, s.as_ref())?));
            }
            "-o" | "--output" => {
                output = Some(PathBuf::from(parse_next(&mut iter, s.as_ref())?));
            }
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

    let input = input.ok_or(ParseError::Usage)?;
    Ok(Config {
        input,
        reference,
        output,
        min_base_q,
        score,
    })
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
    writeln!(
        w,
        "Usage: samtools targetcut [-Q minQ] [-i inPen] [-0 em0] [-1 em1] [-2 em2] <in.bam>"
    )
}

fn run(cfg: &Config, out: &mut dyn Write) -> io::Result<()> {
    let refs = read_reference_lengths(&cfg.input)?;
    let ref_index: HashMap<&str, usize> = refs
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.as_str(), i))
        .collect();

    let columns = pileup_columns(cfg)?;
    let errmod = ErrMod::new(1.0 - ERR_DEP)?;

    let mut last_ref: Option<usize> = None;
    let mut cns: Vec<u16> = Vec::new();

    for column in &columns {
        let Some(&idx) = ref_index.get(column.reference_name.as_str()) else {
            continue;
        };
        if last_ref != Some(idx) {
            if let Some(prev) = last_ref {
                process_cns(&refs[prev].0, &cns, &cfg.score, out)?;
            }
            cns.clear();
            cns.resize(refs[idx].1, 0);
            last_ref = Some(idx);
        }

        let pos0 = column.position.saturating_sub(1);
        if pos0 < cns.len() {
            cns[pos0] = gen_cns(&errmod, cfg.min_base_q, &column.reads_by_input[0]);
        }
    }

    if let Some(prev) = last_ref {
        process_cns(&refs[prev].0, &cns, &cfg.score, out)?;
    }

    Ok(())
}

fn read_reference_lengths(path: &Path) -> io::Result<Vec<(String, usize)>> {
    let format = sam_open_format(path)?;
    let header = match format.exact {
        Exact::Sam => htslib_rs::alignment_compat::read_sam_header_from_path(path)?,
        Exact::Bam => htslib_rs::alignment_compat::read_bam_header_from_path(path)?,
        Exact::Cram => htslib_rs::alignment_compat::read_cram_header_from_path(path)?,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "targetcut input must be SAM, BAM, or CRAM",
            ));
        }
    };

    Ok(header
        .reference_sequences()
        .iter()
        .map(|(name, desc)| {
            (
                String::from_utf8_lossy(name).into_owned(),
                usize::from(desc.length()),
            )
        })
        .collect())
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

fn gen_cns(errmod: &ErrMod, min_base_q: u8, reads: &[PileupRead]) -> u16 {
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
        let mut q = base_q.min(read.mapping_quality);
        q = q.clamp(4, 63);
        bases.push((u16::from(q) << 5) | (u16::from(read.is_reverse) << 4) | u16::from(base));
    }

    if bases.is_empty() {
        return 0;
    }

    let q = errmod.cal(&mut bases, 4);
    let mut sum = [0i32; 4];
    for i in 0..4 {
        sum[i] = ((q[i << 2 | i] + 0.499) as i32) << 2 | i as i32;
    }
    for i in 1..4 {
        let mut j = i;
        while j > 0 && sum[j] < sum[j - 1] {
            sum.swap(j, j - 1);
            j -= 1;
        }
    }

    let qual = ((sum[1] >> 2) - (sum[0] >> 2)).min(63);
    let depth = bases.len().min(255) as u16;
    let ret = ((qual as u16) << 2) | ((sum[0] & 3) as u16);
    (ret << 8) | depth
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

fn process_cns(
    reference_name: &str,
    cns: &[u16],
    score: &ScoreParam,
    out: &mut dyn Write,
) -> io::Result<()> {
    if cns.is_empty() {
        return Ok(());
    }

    let mut b = vec![0u8; cns.len()];
    let mut prev = [0i32, 0i32];
    let mut curr = [0i32, 0i32];

    for (i, &call) in cns.iter().enumerate() {
        let c = if call == 0 {
            0
        } else if call >> 8 == 0 {
            1
        } else {
            2
        };

        let tmp0 = prev[0] + score.e[0][c] + score.p[0][0];
        let tmp1 = prev[1] + score.e[0][c] + score.p[1][0];
        if tmp0 > tmp1 {
            curr[0] = tmp0;
            b[i] = 0;
        } else {
            curr[0] = tmp1;
            b[i] = 1;
        }

        let tmp0 = prev[0] + score.e[1][c] + score.p[0][1];
        let tmp1 = prev[1] + score.e[1][c] + score.p[1][1];
        if tmp0 > tmp1 {
            curr[1] = tmp0;
        } else {
            curr[1] = tmp1;
            b[i] |= 1 << 1;
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    let mut state = if prev[0] > prev[1] { 0u8 } else { 1u8 };
    for i in (1..cns.len()).rev() {
        b[i] |= state << 2;
        state = (b[i] >> state) & 1;
    }

    let mut start: Option<usize> = None;
    for (i, in_target) in b
        .iter()
        .map(|bits| ((bits >> 2) & 3) != 0)
        .chain(std::iter::once(false))
        .enumerate()
    {
        if !in_target {
            if let Some(s) = start.take() {
                write_segment(reference_name, s, i, cns, out)?;
            }
        } else if start.is_none() {
            start = Some(i);
        }
    }

    Ok(())
}

fn write_segment(
    reference_name: &str,
    start0: usize,
    end0: usize,
    cns: &[u16],
    out: &mut dyn Write,
) -> io::Result<()> {
    write!(
        out,
        "{}:{}-{}\t0\t{}\t{}\t60\t{}M\t*\t0\t0\t",
        reference_name,
        start0 + 1,
        end0,
        reference_name,
        start0 + 1,
        end0 - start0
    )?;

    for &call in &cns[start0..end0] {
        let c = call >> 8;
        if c == 0 {
            out.write_all(b"N")?;
        } else {
            out.write_all(&[b"ACGT"[(c & 3) as usize]])?;
        }
    }
    out.write_all(b"\t")?;
    for &call in &cns[start0..end0] {
        out.write_all(&[33 + ((call >> 8) >> 2) as u8])?;
    }
    out.write_all(b"\n")
}

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
                "failed to initialise targetcut error model",
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

        let mut fsum = [0f64; 16];
        let mut bsum = [0f64; 16];
        let mut counts = [0usize; 16];
        let mut strand_counts = [0usize; 32];

        for &packed in bases[..n].iter().rev() {
            let mut qual = (packed >> 5) as usize;
            qual = qual.clamp(4, 63);
            let basestrand = (packed & 0x1f) as usize;
            let base = (packed & 0x0f) as usize;
            fsum[base] += self.fk[strand_counts[basestrand]];
            bsum[base] +=
                self.fk[strand_counts[basestrand]] * self.beta[qual << 16 | n << 8 | counts[base]];
            counts[base] += 1;
            strand_counts[basestrand] += 1;
        }

        for j in 0..m {
            let mut tmp1;
            let mut tmp2;

            tmp1 = 0.0;
            tmp2 = 0usize;
            for k in 0..m {
                if k == j {
                    continue;
                }
                tmp1 += bsum[k];
                tmp2 += counts[k];
            }
            if tmp2 != 0 {
                q[j * m + j] = tmp1 as f32;
            }

            for k in (j + 1)..m {
                let cjk = counts[j] + counts[k];
                tmp1 = 0.0;
                tmp2 = 0;
                for i in 0..m {
                    if i == j || i == k {
                        continue;
                    }
                    tmp1 += bsum[i];
                    tmp2 += counts[i];
                }
                let val = -4.343 * self.lhet[cjk << 8 | counts[k]] + tmp1;
                q[j * m + k] = val as f32;
                q[k * m + j] = q[j * m + k];
                if tmp2 == 0 {
                    // Same expression as the C `else`; kept explicit to
                    // mirror the branch and document that `tmp1 == 0`.
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
        path.push(format!(
            "samtools-rs-targetcut-{name}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        path
    }

    #[test]
    fn targetcut_emits_long_supported_interval() {
        let tmp = tmp_dir("long-interval");
        let sam = tmp.join("in.sam");
        let len = 2600usize;
        let seq = "A".repeat(len);
        let qual = "I".repeat(len);
        std::fs::write(
            &sam,
            format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{len}\nr1\t0\tchr1\t1\t60\t{len}M\t*\t0\t0\t{seq}\t{qual}\n"),
        )
        .unwrap();

        let cfg = Config {
            input: sam,
            reference: None,
            output: None,
            min_base_q: 13,
            score: ScoreParam::default(),
        };
        let mut out = Vec::new();
        run(&cfg, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let fields: Vec<&str> = text.trim_end().split('\t').collect();
        assert_eq!(fields.len(), 11);
        assert!(fields[0].starts_with("chr1:"));
        assert_eq!(fields[1], "0");
        assert_eq!(fields[2], "chr1");
        assert_eq!(fields[4], "60");
        assert!(fields[5].ends_with('M'));
        assert!(fields[9].chars().all(|c| c == 'A'));
        assert_eq!(fields[9].len(), fields[10].len());
    }

    #[test]
    fn targetcut_respects_min_base_quality() {
        let tmp = tmp_dir("min-baseq");
        let sam = tmp.join("in.sam");
        let len = 2600usize;
        let seq = "A".repeat(len);
        let qual = "!".repeat(len);
        std::fs::write(
            &sam,
            format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{len}\nr1\t0\tchr1\t1\t60\t{len}M\t*\t0\t0\t{seq}\t{qual}\n"),
        )
        .unwrap();

        let cfg = Config {
            input: sam,
            reference: None,
            output: None,
            min_base_q: 13,
            score: ScoreParam::default(),
        };
        let mut out = Vec::new();
        run(&cfg, &mut out).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn targetcut_parses_attached_scoring_options() {
        let args = vec![
            OsString::from("targetcut"),
            OsString::from("-Q20"),
            OsString::from("-i123"),
            OsString::from("-010"),
            OsString::from("-111"),
            OsString::from("-212"),
            OsString::from("in.bam"),
        ];
        let cfg = parse_args(&args).unwrap();
        assert_eq!(cfg.min_base_q, 20);
        assert_eq!(cfg.score.p[0][1], -123);
        assert_eq!(cfg.score.e[1], [10, 11, 12]);
        assert_eq!(cfg.input, PathBuf::from("in.bam"));
    }
}
