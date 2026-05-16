//! `samtools mpileup` — multi-way pileup (default text output).
//!
//! Mirrors the text-pileup path of `bam_plcmd.c`. Each covered reference
//! position emits:
//!
//! ```text
//! <chrom>\t<1-based pos>\t<ref base>\t<depth1>\t<bases1>\t<quals1>[\t<depth2>...]
//! ```
//!
//! one depth/bases/quals triple per input file (one sample per file; `@RG`
//! sample grouping is not yet modelled). The read-base encoding follows
//! HTSlib's `pileup_seq` (`.`/`,` reference match, upper/lower mismatch,
//! `^`+mapq head, `$` tail, `*` deletion, `<`/`>` reference skip, `+`/`-`
//! indels). The per-read base-quality gate uses HTSlib's default
//! `--min-BQ 13` and `--ff UNMAP,SECONDARY,QCFAIL,DUP`.
//!
//! Smart overlap removal (`MPLP_SMART_OVERLAPS`) and the orphan filter
//! (`MPLP_NO_ORPHAN`, cleared by `-A`) are applied via the htslib-rs pileup
//! engine. Byte parity verified against upstream `mpileup.out.3`
//! (`-B --ff`) and `mpileup.out.5` (overlap).
//!
//! **Parity gap (tracked):** BAQ recomputation (HTSlib default, disabled
//! with `-B`; needs completed library batch #11) is not yet applied, so a handful of
//! base-quality characters differ on non-`-B` inputs (depths and read
//! bases match exactly — see `mpileup.out.1`).

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::alignment_compat::{
    PileupColumn, PileupOptions, PileupRead, pileup_from_alignment_paths_with_options,
    pileup_from_alignment_paths_with_reference_and_options,
};

use crate::bam_flag::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FUNMAP, str_to_flag};
use crate::diagnostics::{print_error, print_error_errno};

struct Config {
    inputs: Vec<PathBuf>,
    reference: Option<PathBuf>,
    region: Option<String>,
    output: Option<PathBuf>,
    min_base_q: u8,
    min_map_q: u8,
    exclude_flags: u16,
    require_flags: u16,
    detect_overlaps: bool,
    count_orphans: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            inputs: Vec::new(),
            reference: None,
            region: None,
            output: None,
            min_base_q: 13,
            min_map_q: 0,
            exclude_flags: (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) as u16,
            require_flags: 0,
            detect_overlaps: true,
            count_orphans: false,
        }
    }
}

/// Splits attached short-option values (`-r17:1-2`) and `--long=value` into
/// separate tokens so the parser only sees the canonical split form.
fn normalize_args(args: &[OsString]) -> Vec<OsString> {
    // Short options that take a value.
    const VALUE_SHORT: &[u8] = b"frbozQqGCdlT";
    let mut out = Vec::with_capacity(args.len());
    for arg in args.iter().skip(1) {
        let s = arg.to_str().unwrap_or("");
        if let Some(rest) = s.strip_prefix("--") {
            if let Some((name, val)) = rest.split_once('=') {
                out.push(OsString::from(format!("--{name}")));
                out.push(OsString::from(val));
                continue;
            }
            out.push(arg.clone());
        } else if s.len() > 2 && s.starts_with('-') && VALUE_SHORT.contains(&s.as_bytes()[1]) {
            out.push(OsString::from(format!("-{}", &s[1..2])));
            out.push(OsString::from(&s[2..]));
        } else {
            out.push(arg.clone());
        }
    }
    out
}

/// Entry point for `samtools mpileup`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut cfg = Config::default();
    let mut bam_list: Option<PathBuf> = None;

    let normalized = normalize_args(args);
    let mut iter = normalized.iter();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-f" | "--fasta-ref" => cfg.reference = iter.next().map(PathBuf::from),
            "-r" | "--region" => {
                cfg.region = iter.next().and_then(|a| a.to_str().map(str::to_owned));
            }
            "-b" | "--bam-list" => bam_list = iter.next().map(PathBuf::from),
            "-o" | "--output" => cfg.output = iter.next().map(PathBuf::from),
            "-Q" | "--min-BQ" => {
                cfg.min_base_q = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(13);
            }
            "-q" | "--min-MQ" => {
                cfg.min_map_q = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0);
            }
            "--ff" | "--excl-flags" => match iter.next().and_then(|a| a.to_str()) {
                Some(v) => match str_to_flag(v) {
                    Some(f) => cfg.exclude_flags = f as u16,
                    None => {
                        eprintln!("Could not parse --ff {v}");
                        return ExitCode::from(1);
                    }
                },
                None => return ExitCode::from(1),
            },
            "--rf" | "--incl-flags" | "--require-flags" => {
                match iter.next().and_then(|a| a.to_str()) {
                    Some(v) => match str_to_flag(v) {
                        Some(f) => cfg.require_flags = f as u16,
                        None => return ExitCode::from(1),
                    },
                    None => return ExitCode::from(1),
                }
            }
            // Accepted no-ops / not-yet-modelled boolean options.
            "-A" | "--count-orphans" => cfg.count_orphans = true,
            "-B" | "--no-BAQ" | "-E" | "--redo-BAQ" | "-a" | "-aa" | "-R" | "--ignore-RG"
            | "-s" | "--output-MQ" | "-O" | "--output-BP" | "-M" | "--output-mods" | "-6"
            | "--illumina1.3+" | "--no-PG" => {}
            "-x" | "--ignore-overlaps-removal" | "--ignore-overlaps" => {
                cfg.detect_overlaps = false;
            }
            // Accepted options whose value is consumed but not yet modelled.
            "-d" | "--max-depth" | "-l" | "--positions" | "-G" | "--exclude-RG" | "-C"
            | "--adjust-MQ" => {
                let _ = iter.next();
            }
            _ if s.starts_with('-') && s != "-" => {
                // Unknown flag: ignore conservatively.
            }
            _ => cfg.inputs.push(PathBuf::from(arg)),
        }
    }

    if let Some(list) = bam_list {
        match read_input_list(&list) {
            Ok(mut paths) => cfg.inputs.append(&mut paths),
            Err(e) => {
                let msg = format!("failed to read {}", list.display());
                print_error_errno("mpileup", &msg, &e);
                return ExitCode::from(1);
            }
        }
    }

    if cfg.inputs.is_empty() {
        print_error("mpileup", "no input files");
        return ExitCode::from(1);
    }

    match run(&cfg) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("mpileup", "pileup failed", &e);
            ExitCode::from(1)
        }
    }
}

fn read_input_list(path: &Path) -> io::Result<Vec<PathBuf>> {
    let file = File::open(path)?;
    let mut out = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let line = line.strip_prefix("file://").unwrap_or(line);
        out.push(PathBuf::from(line));
    }
    Ok(out)
}

/// Parses a region string into `(name, 1-based-begin, 1-based-end-inclusive)`.
fn parse_region(spec: &str) -> (String, usize, usize) {
    match spec.rsplit_once(':') {
        None => (spec.to_string(), 1, usize::MAX),
        Some((name, range)) => {
            let range = range.replace(',', "");
            let (b, e) = match range.split_once('-') {
                None => {
                    let b = range.parse().unwrap_or(1);
                    (b, b)
                }
                Some((bs, es)) => {
                    let b = if bs.is_empty() {
                        1
                    } else {
                        bs.parse().unwrap_or(1)
                    };
                    let e = if es.is_empty() {
                        usize::MAX
                    } else {
                        es.parse().unwrap_or(usize::MAX)
                    };
                    (b, e)
                }
            };
            (name.to_string(), b, e)
        }
    }
}

/// Reads a (possibly bgzipped) FASTA into a name → sequence map.
fn read_fasta(path: &Path) -> io::Result<HashMap<String, Vec<u8>>> {
    let mut bytes = Vec::new();
    let file = File::open(path)?;
    if path.extension().is_some_and(|e| e == "gz" || e == "bgz") {
        htslib_rs::bgzf::io::Reader::new(file).read_to_end(&mut bytes)?;
    } else {
        BufReader::new(file).read_to_end(&mut bytes)?;
    }

    let mut refs = HashMap::new();
    let mut name: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    for line in bytes.split(|&b| b == b'\n') {
        if line.first() == Some(&b'>') {
            if let Some(n) = name.take() {
                refs.insert(n, std::mem::take(&mut seq));
            }
            let header = &line[1..];
            let end = header
                .iter()
                .position(|b| b.is_ascii_whitespace())
                .unwrap_or(header.len());
            name = Some(String::from_utf8_lossy(&header[..end]).into_owned());
        } else if line.first() != Some(&b';') {
            seq.extend(line.iter().filter(|b| !b.is_ascii_whitespace()));
        }
    }
    if let Some(n) = name {
        refs.insert(n, seq);
    }
    Ok(refs)
}

fn run(cfg: &Config) -> io::Result<()> {
    let options = PileupOptions {
        exclude_flags: cfg.exclude_flags,
        require_flags: cfg.require_flags,
        min_mapping_quality: cfg.min_map_q,
        detect_overlaps: cfg.detect_overlaps,
        discard_orphans: !cfg.count_orphans,
    };

    let has_cram = cfg.inputs.iter().any(|p| is_cram(p));
    let columns: Vec<PileupColumn> = if has_cram || cfg.reference.is_some() {
        let reference = cfg.reference.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "CRAM input requires -f/--fasta-ref",
            )
        })?;
        match pileup_from_alignment_paths_with_reference_and_options(
            &cfg.inputs,
            reference,
            &options,
        ) {
            Ok(c) => c,
            Err(e) if !has_cram => {
                // Reference index missing but only used for the ref-base
                // column: fall back to a plain SAM/BAM pileup.
                let _ = e;
                pileup_from_alignment_paths_with_options(&cfg.inputs, &options)?
            }
            Err(e) => return Err(e),
        }
    } else {
        pileup_from_alignment_paths_with_options(&cfg.inputs, &options)?
    };

    let refs = match cfg.reference.as_ref() {
        Some(p) => Some(read_fasta(p)?),
        None => None,
    };
    let region = cfg.region.as_deref().map(parse_region);

    let mut writer: Box<dyn Write> = match cfg.output.as_ref() {
        Some(p) => Box::new(io::BufWriter::new(File::create(p)?)),
        None => Box::new(io::BufWriter::new(io::stdout().lock())),
    };

    eprintln!(
        "[mpileup] {} samples in {} input files",
        cfg.inputs.len(),
        cfg.inputs.len()
    );

    let mut line: Vec<u8> = Vec::new();
    for column in &columns {
        if let Some((name, beg, end)) = &region
            && (column.reference_name != *name || column.position < *beg || column.position > *end)
        {
            continue;
        }

        let ref_seq = refs
            .as_ref()
            .and_then(|m| m.get(&column.reference_name))
            .map(Vec::as_slice);
        let pos0 = column.position - 1;
        let ref_base = ref_seq
            .and_then(|s| s.get(pos0).copied())
            .unwrap_or(b'N')
            .to_ascii_uppercase();

        line.clear();
        line.extend_from_slice(column.reference_name.as_bytes());
        write!(line, "\t{}\t", column.position).unwrap();
        line.push(ref_base);

        for reads in &column.reads_by_input {
            let mut seq: Vec<u8> = Vec::new();
            let mut qual: Vec<u8> = Vec::new();
            let mut cnt = 0usize;
            for r in reads {
                let bq = r.qpos_quality;
                if bq < cfg.min_base_q {
                    continue;
                }
                encode_read(&mut seq, r, pos0, ref_seq);
                qual.push(if (bq as u16) + 33 < 126 { bq + 33 } else { 126 });
                cnt += 1;
            }
            write!(line, "\t{cnt}\t").unwrap();
            if seq.is_empty() {
                line.push(b'*');
            } else {
                line.extend_from_slice(&seq);
            }
            line.push(b'\t');
            if qual.is_empty() {
                line.push(b'*');
            } else {
                line.extend_from_slice(&qual);
            }
        }
        line.push(b'\n');
        writer.write_all(&line)?;
    }

    writer.flush()
}

fn is_cram(p: &Path) -> bool {
    p.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("cram"))
}

/// Ports HTSlib's `pileup_seq` for one read at a column.
fn encode_read(out: &mut Vec<u8>, r: &PileupRead, pos0: usize, ref_seq: Option<&[u8]>) {
    if r.is_head {
        out.push(b'^');
        let mq = r.mapping_quality;
        out.push(if mq > 93 { 126 } else { mq + 33 });
    }

    if !r.is_deletion && !r.is_refskip {
        let base = r.base.unwrap_or(b'N');
        let matches = ref_seq
            .and_then(|s| s.get(pos0).copied())
            .is_some_and(|rb| rb.eq_ignore_ascii_case(&base));
        if matches {
            out.push(if r.is_reverse { b',' } else { b'.' });
        } else if r.is_reverse {
            out.push(base.to_ascii_lowercase());
        } else {
            out.push(base.to_ascii_uppercase());
        }
    } else if r.is_refskip {
        out.push(if r.is_reverse { b'<' } else { b'>' });
    } else {
        out.push(b'*');
    }

    if r.indel > 0 {
        out.push(b'+');
        write!(out, "{}", r.indel).unwrap();
        for &b in &r.insertion {
            out.push(if r.is_reverse {
                b.to_ascii_lowercase()
            } else {
                b.to_ascii_uppercase()
            });
        }
    } else if r.indel < 0 {
        let del = (-r.indel) as usize;
        out.push(b'-');
        write!(out, "{del}").unwrap();
        for j in 1..=del {
            let c = ref_seq
                .and_then(|s| s.get(pos0 + j).copied())
                .unwrap_or(b'N');
            out.push(if r.is_reverse {
                c.to_ascii_lowercase()
            } else {
                c.to_ascii_uppercase()
            });
        }
    }

    if r.is_tail {
        out.push(b'$');
    }
}
