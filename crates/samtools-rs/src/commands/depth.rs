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
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::core::Region;
use htslib_rs::format::Exact;
use htslib_rs::sam;

use crate::bam_flag::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FUNMAP};
use crate::bedidx::load_bed_index;
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

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

pub(crate) struct DepthRunConfig<'a> {
    pub(crate) min_mapq: u8,
    pub(crate) min_depth: u32,
    pub(crate) a_mode: AMode,
    pub(crate) region: Option<&'a str>,
    pub(crate) bed: Option<&'a Path>,
    pub(crate) reference: Option<&'a Path>,
}

#[derive(Clone, Copy)]
struct DepthWalkConfig {
    exclude_flags: u32,
    min_mapq: u8,
    min_depth: u32,
    a_mode: AMode,
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

    let mut has_cram = false;
    for path in &inputs {
        let format = match sam_io::sam_open_format(path) {
            Ok(f) => f,
            Err(e) => {
                print_error("depth", e.to_string());
                return ExitCode::from(1);
            }
        };
        match format.exact {
            Exact::Sam | Exact::Bam => {}
            Exact::Cram => has_cram = true,
            _ => {
                print_error(
                    "depth",
                    "only SAM, BAM, and reference-backed CRAM input are currently supported",
                );
                return ExitCode::from(1);
            }
        }
    }

    let reference = if has_cram {
        match current_global_args().reference {
            Some(reference) => Some(reference),
            None => {
                print_error("depth", "CRAM input requires top-level --reference FILE");
                return ExitCode::from(1);
            }
        }
    } else {
        None
    };

    if has_cram && regions_need_index(region.as_deref(), bed.as_deref()).is_err() {
        print_error(
            "depth",
            "CRAM input requires indexed region-compatible input",
        );
        return ExitCode::from(1);
    }

    let mut writer = match sam_io::open_text_output(output.as_deref()) {
        Ok(writer) => writer,
        Err(e) => {
            print_error_errno("depth", "open -o output", &e);
            return ExitCode::from(1);
        }
    };

    match run_depth(
        &inputs,
        &mut *writer,
        DepthRunConfig {
            min_mapq,
            min_depth,
            a_mode,
            region: region.as_deref(),
            bed: bed.as_deref(),
            reference: reference.as_deref(),
        },
    ) {
        Ok(()) => match sam_io::check_sam_close(&mut writer) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
            Err(e) => {
                print_error_errno("depth", "close output", &e);
                ExitCode::from(1)
            }
        },
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("depth", "depth failed", &e);
            ExitCode::from(1)
        }
    }
}

fn regions_need_index(region: Option<&str>, bed: Option<&Path>) -> io::Result<()> {
    if let Some(region) = region {
        parse_region(region)?;
    }
    if let Some(bed) = bed {
        let _ = load_bed_regions(bed)?;
    }
    Ok(())
}

pub(crate) fn run_depth(
    inputs: &[PathBuf],
    out: &mut dyn Write,
    config: DepthRunConfig<'_>,
) -> io::Result<()> {
    let exclude_flags = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    let walk = DepthWalkConfig {
        exclude_flags,
        min_mapq: config.min_mapq,
        min_depth: config.min_depth,
        a_mode: config.a_mode,
    };
    let mut regions = Vec::new();
    if let Some(region) = config.region {
        regions.push(parse_region(region)?);
    }
    if let Some(bed) = config.bed {
        regions.extend(load_bed_regions(bed)?);
    }

    let mut per_input_targets = Vec::with_capacity(inputs.len());
    for path in inputs {
        let format = sam_io::sam_open_format(path)?;
        let targets = match format.exact {
            Exact::Sam => collect_sam_depth(path, walk, &regions)?,
            Exact::Bam => collect_bam_depth(path, walk, &regions)?,
            Exact::Cram => collect_cram_depth(
                path,
                config.reference.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "CRAM input requires top-level --reference FILE",
                    )
                })?,
                walk,
                &regions,
            )?,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "only SAM, BAM, and reference-backed CRAM input are currently supported",
                ));
            }
        };
        per_input_targets.push(targets);
    }

    emit_depths(out, &per_input_targets, walk.min_depth, walk.a_mode)
}

fn collect_sam_depth(
    path: &Path,
    config: DepthWalkConfig,
    regions: &[Region],
) -> io::Result<Vec<DepthTarget>> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(path)?));
    let header = reader.read_header()?;
    let mut targets = depth_targets(&header, regions)?;

    for result in reader.records() {
        let record = result?;
        update_targets(
            &header,
            &mut targets,
            &record,
            config.exclude_flags,
            config.min_mapq,
        );
    }

    Ok(targets)
}

fn collect_bam_depth(
    path: &Path,
    config: DepthWalkConfig,
    regions: &[Region],
) -> io::Result<Vec<DepthTarget>> {
    let mut reader = bam::io::Reader::new(File::open(path)?);
    let header = reader.read_header()?;
    let mut targets = depth_targets(&header, regions)?;

    if regions.is_empty() {
        let mut record = bam::Record::default();
        loop {
            let n = reader.read_record(&mut record)?;
            if n == 0 {
                break;
            }
            update_targets(
                &header,
                &mut targets,
                &record,
                config.exclude_flags,
                config.min_mapq,
            );
        }
    } else {
        for (i, region) in regions.iter().enumerate() {
            for record in htslib_rs::alignment_compat::query_bam_records_from_path(path, region)? {
                update_target(
                    &header,
                    &mut targets[i],
                    &record,
                    config.exclude_flags,
                    config.min_mapq,
                );
            }
        }
    }

    Ok(targets)
}

fn collect_cram_depth(
    path: &Path,
    reference: &Path,
    config: DepthWalkConfig,
    regions: &[Region],
) -> io::Result<Vec<DepthTarget>> {
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(path)?;
    let mut targets = depth_targets(&header, regions)?;

    if regions.is_empty() {
        for target in &mut targets {
            let region = target_region(target)?;
            for record in htslib_rs::alignment_compat::query_cram_records_from_path_with_reference(
                path, &region, reference,
            )? {
                update_target(
                    &header,
                    target,
                    &record,
                    config.exclude_flags,
                    config.min_mapq,
                );
            }
        }
    } else {
        for (i, region) in regions.iter().enumerate() {
            for record in htslib_rs::alignment_compat::query_cram_records_from_path_with_reference(
                path, region, reference,
            )? {
                update_target(
                    &header,
                    &mut targets[i],
                    &record,
                    config.exclude_flags,
                    config.min_mapq,
                );
            }
        }
    }

    Ok(targets)
}

fn target_region(target: &DepthTarget) -> io::Result<Region> {
    format!("{}:{}-{}", target.name, target.output_start, target.end0)
        .parse::<Region>()
        .map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "region \"{}:{}-{}\": {}",
                    target.name, target.output_start, target.end0, e
                ),
            )
        })
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
    load_bed_index(path)?.to_htslib_regions()
}

fn emit_depths(
    out: &mut dyn Write,
    per_input_targets: &[Vec<DepthTarget>],
    min_depth: u32,
    a_mode: AMode,
) -> io::Result<()> {
    let Some(first_targets) = per_input_targets.first() else {
        return Ok(());
    };

    for targets in &per_input_targets[1..] {
        ensure_compatible_targets(first_targets, targets)?;
    }

    for (target_index, target) in first_targets.iter().enumerate() {
        let has_any = per_input_targets
            .iter()
            .any(|targets| targets[target_index].depths.iter().any(|&d| d > 0));
        if !has_any && !matches!(a_mode, AMode::AllRefsAllPositions) {
            continue;
        }
        for i in 0..target.depths.len() {
            let depths = per_input_targets
                .iter()
                .map(|targets| targets[target_index].depths[i]);
            if a_mode == AMode::None && !depths.clone().any(|d| d > 0) {
                continue;
            }
            if a_mode == AMode::None && !depths.clone().any(|d| d >= min_depth) {
                continue;
            }
            write!(out, "{}\t{}", target.name, target.output_start + i)?;
            for d in depths {
                write!(out, "\t{}", d)?;
            }
            writeln!(out)?;
        }
    }
    Ok(())
}

fn ensure_compatible_targets(expected: &[DepthTarget], actual: &[DepthTarget]) -> io::Result<()> {
    if expected.len() != actual.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "depth inputs have incompatible reference dictionaries",
        ));
    }

    for (left, right) in expected.iter().zip(actual) {
        if left.name != right.name
            || left.output_start != right.output_start
            || left.start0 != right.start0
            || left.end0 != right.end0
            || left.depths.len() != right.depths.len()
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "depth inputs have incompatible reference dictionaries",
            ));
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
    header: &sam::Header,
    targets: &mut [DepthTarget],
    record: &(impl sam::alignment::Record + ?Sized),
    exclude_flags: u32,
    min_mapq: u8,
) {
    let flag = match record.flags() {
        Ok(flags) => u16::from(flags) as u32,
        Err(_) => return,
    };
    if flag & exclude_flags != 0 {
        return;
    }
    let mapq = match record.mapping_quality() {
        Some(Ok(q)) => u8::from(q),
        Some(Err(_)) => return,
        None => 0,
    };
    if mapq < min_mapq {
        return;
    }
    let tid = match record.reference_sequence_id(header).transpose() {
        Ok(Some(t)) => t,
        _ => return,
    };
    let start = match record.alignment_start().transpose() {
        Ok(Some(p)) => usize::from(p) - 1,
        _ => return,
    };
    for target in targets.iter_mut().filter(|target| target.tid == tid) {
        update_target_cigar(target, record, start);
    }
}

fn update_target(
    header: &sam::Header,
    target: &mut DepthTarget,
    record: &(impl sam::alignment::Record + ?Sized),
    exclude_flags: u32,
    min_mapq: u8,
) {
    let flag = match record.flags() {
        Ok(flags) => u16::from(flags) as u32,
        Err(_) => return,
    };
    if flag & exclude_flags != 0 {
        return;
    }
    let mapq = match record.mapping_quality() {
        Some(Ok(q)) => u8::from(q),
        Some(Err(_)) => return,
        None => 0,
    };
    if mapq < min_mapq {
        return;
    }
    if record
        .reference_sequence_id(header)
        .transpose()
        .unwrap_or_default()
        != Some(target.tid)
    {
        return;
    }
    let start = match record.alignment_start().transpose() {
        Ok(Some(p)) => usize::from(p) - 1,
        _ => return,
    };
    update_target_cigar(target, record, start);
}

fn update_target_cigar(
    target: &mut DepthTarget,
    record: &(impl sam::alignment::Record + ?Sized),
    start: usize,
) {
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
