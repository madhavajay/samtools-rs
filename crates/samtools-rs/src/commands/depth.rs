//! `samtools depth` — per-position depth of alignments.
//!
//! Mirrors `main_depth` in `bam2depth.c`. Upstream uses pileup. This Rust
//! port builds a per-reference depth vector by walking each record's CIGAR,
//! then emits `<chr>\t<pos>\t<depth>` lines for positions with `depth >=
//! min-depth` (default 1).
//!
//! Output covers only positions with depth ≥ threshold; pass `-a` to emit
//! every reference position, or `-aa` to also include references with no
//! coverage.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::core::Region;
use htslib_rs::format::{Exact, detect_path};

use crate::bam_flag::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FUNMAP};
use crate::diagnostics::{print_error, print_error_errno};

#[derive(Default, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AMode {
    #[default]
    /// Only positions with depth > 0
    None,
    /// `-a` — all positions on references that have any coverage
    AllPositions,
    /// `-aa` — all positions on all references
    AllRefsAllPositions,
}

/// Entry point for `samtools depth`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut min_mapq: u8 = 0;
    let mut min_depth: u32 = 1;
    let mut a_mode = AMode::None;
    let mut output: Option<PathBuf> = None;
    let mut region: Option<String> = None;
    let mut bed: Option<PathBuf> = None;
    let mut inputs: Vec<PathBuf> = Vec::new();

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-a" => {
                a_mode = match a_mode {
                    AMode::None => AMode::AllPositions,
                    _ => AMode::AllRefsAllPositions,
                };
            }
            "-aa" => {
                a_mode = AMode::AllRefsAllPositions;
            }
            "-q" | "--min-MQ" => {
                min_mapq = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0);
            }
            "-d" | "--min-depth" => {
                min_depth = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(1);
            }
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "-r" | "--region" => {
                region = iter.next().and_then(|a| a.to_str().map(str::to_owned));
            }
            "-b" => {
                bed = iter.next().map(PathBuf::from);
            }
            "-l" | "--min-read-len" | "-Q" | "--min-BQ" | "--rf" | "--ff" | "-m"
            | "--max-depth" | "-G" | "-g" | "-f" => {
                let _ = iter.next();
            }
            "-H" | "-J" | "-s" | "-x" | "-X" => {
                // Header / scientific / strip / index variations not yet supported.
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("depth", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => inputs.push(PathBuf::from(arg)),
        }
    }

    if inputs.is_empty() {
        let _ = print_usage();
        return ExitCode::from(1);
    }

    for path in &inputs {
        let format = match detect_path(path) {
            Ok(f) => f,
            Err(e) => {
                print_error(
                    "depth",
                    format!("failed to detect format of \"{}\": {}", path.display(), e),
                );
                return ExitCode::from(1);
            }
        };
        if format.exact != Exact::Bam {
            print_error(
                "depth",
                "only BAM input is currently supported (SAM/CRAM TODO)",
            );
            return ExitCode::from(1);
        }
    }

    let mut writer: Box<dyn Write> = match output.as_ref() {
        Some(p) => match File::create(p) {
            Ok(f) => Box::new(f),
            Err(e) => {
                print_error_errno("depth", "open -o output", &e);
                return ExitCode::from(1);
            }
        },
        None => Box::new(io::stdout().lock()),
    };

    match run_depth(
        &inputs,
        &mut *writer,
        min_mapq,
        min_depth,
        a_mode,
        region.as_deref(),
        bed.as_deref(),
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("depth", "depth failed", &e);
            ExitCode::from(1)
        }
    }
}

pub(crate) fn run_depth(
    inputs: &[PathBuf],
    out: &mut dyn Write,
    min_mapq: u8,
    min_depth: u32,
    a_mode: AMode,
    region: Option<&str>,
    bed: Option<&Path>,
) -> io::Result<()> {
    let exclude_flags = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    let mut regions = Vec::new();
    if let Some(region) = region {
        regions.push(parse_region(region)?);
    }
    if let Some(bed) = bed {
        regions.extend(load_bed_regions(bed)?);
    }

    for path in inputs {
        let mut reader = bam::io::Reader::new(File::open(path)?);
        let header = reader.read_header()?;
        let mut targets = depth_targets(&header, &regions)?;

        if regions.is_empty() {
            let mut record = bam::Record::default();
            loop {
                let n = reader.read_record(&mut record)?;
                if n == 0 {
                    break;
                }
                update_targets(&mut targets, &record, exclude_flags, min_mapq);
            }
        } else {
            for (i, region) in regions.iter().enumerate() {
                for record in
                    htslib_rs::alignment_compat::query_bam_records_from_path(path, region)?
                {
                    update_target(&mut targets[i], &record, exclude_flags, min_mapq);
                }
            }
        }

        emit_depths(out, &targets, min_depth, a_mode)?;
    }
    Ok(())
}

fn parse_region(s: &str) -> io::Result<Region> {
    s.parse::<Region>().map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("region \"{}\": {}", s, e),
        )
    })
}

fn load_bed_regions(path: &Path) -> io::Result<Vec<Region>> {
    let file = File::open(path)?;
    let mut out = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        let s = line.trim_end();
        if s.is_empty()
            || s.starts_with('#')
            || s.starts_with("track ")
            || s.starts_with("browser ")
        {
            continue;
        }

        let mut fields = s.split('\t');
        let chrom = fields.next().unwrap_or("");
        let beg: u64 = fields.next().and_then(|t| t.parse().ok()).unwrap_or(0);
        let end: u64 = fields.next().and_then(|t| t.parse().ok()).unwrap_or(0);
        if chrom.is_empty() || end <= beg {
            continue;
        }

        out.push(parse_region(&format!("{}:{}-{}", chrom, beg + 1, end))?);
    }
    Ok(out)
}

fn emit_depths(
    out: &mut dyn Write,
    targets: &[DepthTarget],
    min_depth: u32,
    a_mode: AMode,
) -> io::Result<()> {
    for target in targets {
        let has_any = target.depths.iter().any(|&d| d > 0);
        if !has_any && !matches!(a_mode, AMode::AllRefsAllPositions) {
            continue;
        }
        for (i, &d) in target.depths.iter().enumerate() {
            if d == 0 && a_mode == AMode::None {
                continue;
            }
            if d < min_depth && a_mode == AMode::None {
                continue;
            }
            writeln!(out, "{}\t{}\t{}", target.name, target.output_start + i, d)?;
        }
    }
    Ok(())
}

struct DepthTarget {
    tid: usize,
    name: String,
    output_start: usize,
    start0: usize,
    end0: usize,
    depths: Vec<u32>,
}

fn depth_targets(
    header: &htslib_rs::sam::Header,
    regions: &[Region],
) -> io::Result<Vec<DepthTarget>> {
    if regions.is_empty() {
        Ok(header
            .reference_sequences()
            .iter()
            .enumerate()
            .map(|(tid, (name, def))| {
                let length = usize::from(def.length());
                DepthTarget {
                    tid,
                    name: String::from_utf8_lossy(name).into_owned(),
                    output_start: 1,
                    start0: 0,
                    end0: length,
                    depths: vec![0u32; length],
                }
            })
            .collect())
    } else {
        regions
            .iter()
            .map(|region| depth_target_for_region(header, region))
            .collect()
    }
}

fn depth_target_for_region(
    header: &htslib_rs::sam::Header,
    region: &Region,
) -> io::Result<DepthTarget> {
    let tid = header
        .reference_sequences()
        .get_index_of(region.name())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "region reference sequence does not exist: {}",
                    String::from_utf8_lossy(region.name())
                ),
            )
        })?;
    let (_, def) = header.reference_sequences().get_index(tid).ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence ID")
    })?;
    let ref_len = usize::from(def.length());
    let interval = region.interval();
    let output_start = interval.start().map(usize::from).unwrap_or(1);
    let output_end = interval.end().map(usize::from).unwrap_or(ref_len);
    if output_end < output_start {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid region interval: {}", region),
        ));
    }

    let start0 = output_start - 1;
    let end0 = output_end.min(ref_len);
    let length = output_end - output_start + 1;

    Ok(DepthTarget {
        tid,
        name: String::from_utf8_lossy(region.name()).into_owned(),
        output_start,
        start0,
        end0,
        depths: vec![0u32; length],
    })
}

fn update_targets(
    targets: &mut [DepthTarget],
    record: &bam::Record,
    exclude_flags: u32,
    min_mapq: u8,
) {
    let flag = u16::from(record.flags()) as u32;
    if flag & exclude_flags != 0 {
        return;
    }
    let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
    if mapq < min_mapq {
        return;
    }
    let tid = match record.reference_sequence_id().and_then(|r| r.ok()) {
        Some(t) => t,
        None => return,
    };
    let start = match record.alignment_start().and_then(|r| r.ok()) {
        Some(p) => usize::from(p) - 1,
        None => return,
    };
    for target in targets.iter_mut().filter(|target| target.tid == tid) {
        update_target_cigar(target, record, start);
    }
}

fn update_target(target: &mut DepthTarget, record: &bam::Record, exclude_flags: u32, min_mapq: u8) {
    let flag = u16::from(record.flags()) as u32;
    if flag & exclude_flags != 0 {
        return;
    }
    let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
    if mapq < min_mapq {
        return;
    }
    if record.reference_sequence_id().and_then(|r| r.ok()) != Some(target.tid) {
        return;
    }
    let start = match record.alignment_start().and_then(|r| r.ok()) {
        Some(p) => usize::from(p) - 1,
        None => return,
    };
    update_target_cigar(target, record, start);
}

fn update_target_cigar(target: &mut DepthTarget, record: &bam::Record, start: usize) {
    let mut ref_pos = start;
    for op in record.cigar().iter() {
        let op = match op {
            Ok(op) => op,
            Err(_) => break,
        };
        let len = op.len();
        use htslib_rs::sam::alignment::record::cigar::op::Kind;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let op_end = ref_pos.saturating_add(len);
                let lo = ref_pos.max(target.start0);
                let hi = op_end.min(target.end0);
                if hi > lo {
                    for p in lo..hi {
                        let offset = p - target.start0;
                        if offset < target.depths.len() {
                            target.depths[offset] = target.depths[offset].saturating_add(1);
                        }
                    }
                }
                ref_pos = op_end;
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos = ref_pos.saturating_add(len);
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools depth [options] <in.bam>")?;
    writeln!(w, "  -a          all positions (on covered refs)")?;
    writeln!(w, "  -aa         all positions on all refs")?;
    writeln!(w, "  -d INT      minimum depth threshold [1]")?;
    writeln!(w, "  -q INT      minimum mapq [0]")?;
    writeln!(w, "  -o FILE     output FILE")?;
    writeln!(w, "  -r REGION   restrict to REGION")?;
    writeln!(w, "  -b FILE     restrict to BED regions")?;
    Ok(())
}
