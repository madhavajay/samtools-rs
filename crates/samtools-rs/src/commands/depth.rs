//! `samtools depth` — per-position depth of alignments.
//!
//! Mirrors `main_depth` in `bam2depth.c`. Upstream uses pileup. This Rust
//! port accumulates depth by walking each record's CIGAR into a **sparse**
//! per-reference map (only covered positions), so very large references
//! (e.g. the upstream `large_pos` fixture, `LN:10001009800`) do not OOM,
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
use htslib_rs::format::Exact;
use htslib_rs::sam;

use crate::bam_flag::{
    BAM_FDUP, BAM_FPAIRED, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FUNMAP, str_to_flag,
};
use crate::bedidx::parse_bed_line;
use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
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
    pub(crate) min_read_len: usize,
    pub(crate) a_mode: AMode,
    pub(crate) show_header: bool,
    pub(crate) exclude_flags: u32,
    pub(crate) include_any_flags: u32,
    pub(crate) require_flags: u32,
    pub(crate) include_deletions: bool,
    pub(crate) remove_overlaps: bool,
    pub(crate) region: Option<&'a str>,
    pub(crate) bed: Option<&'a Path>,
    pub(crate) reference: Option<&'a Path>,
}

#[derive(Clone, Copy)]
struct DepthWalkConfig {
    exclude_flags: u32,
    include_any_flags: u32,
    require_flags: u32,
    min_mapq: u8,
    min_read_len: usize,
    min_depth: u32,
    a_mode: AMode,
    include_deletions: bool,
    remove_overlaps: bool,
}

struct DepthRegion {
    region: Region,
    emit_empty: bool,
}

/// Entry point for `samtools depth`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut min_mapq: u8 = 0;
    let mut min_depth: u32 = 1;
    let mut min_read_len: usize = 0;
    let mut a_mode = AMode::None;
    let mut output: Option<PathBuf> = None;
    let mut region: Option<String> = None;
    let mut bed: Option<PathBuf> = None;
    let mut show_header = false;
    let mut exclude_flags = default_exclude_flags();
    let mut include_any_flags = 0;
    let mut require_flags = 0;
    let mut include_deletions = false;
    let mut remove_overlaps = false;
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
            "-l" | "--min-read-len" => {
                min_read_len = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0);
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
            "-g" => match parse_flag_value(iter.next(), "-g") {
                Ok(flags) => exclude_flags &= !flags,
                Err(()) => return ExitCode::from(1),
            },
            "-G" | "--excl-flags" => match parse_flag_value(iter.next(), s) {
                Ok(flags) => exclude_flags |= flags,
                Err(()) => return ExitCode::from(1),
            },
            "--incl-flags" => match parse_flag_value(iter.next(), s) {
                Ok(flags) => include_any_flags |= flags,
                Err(()) => return ExitCode::from(1),
            },
            "--require-flags" => match parse_flag_value(iter.next(), s) {
                Ok(flags) => require_flags |= flags,
                Err(()) => return ExitCode::from(1),
            },
            "-f" => {
                let Some(path) = iter.next().map(PathBuf::from) else {
                    print_error("depth", "option -f requires an argument");
                    return ExitCode::from(1);
                };
                match read_input_list(&path) {
                    Ok(listed_inputs) => inputs.extend(listed_inputs),
                    Err(e) => {
                        print_error_errno("depth", "read -f input list", &e);
                        return ExitCode::from(1);
                    }
                }
            }
            "-Q" | "--min-BQ" | "--rf" | "--ff" | "-m" | "--max-depth" => {
                let _ = iter.next();
            }
            "-H" => {
                show_header = true;
            }
            "-J" => {
                include_deletions = true;
            }
            "-s" => {
                remove_overlaps = true;
            }
            "-x" | "-X" => {
                // Scientific / custom-index variations not yet supported.
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
        if path.as_os_str() != "-" && !path.exists() {
            print_hts_open_missing(path);
            print_error(
                "depth",
                format!(
                    "Cannot open input file \"{}\": No such file or directory",
                    path.display()
                ),
            );
            return ExitCode::from(1);
        }
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
                    "only SAM, BAM, and CRAM input are currently supported",
                );
                return ExitCode::from(1);
            }
        }
    }
    let reference = has_cram.then(|| current_global_args().reference).flatten();

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
            min_read_len,
            a_mode,
            show_header,
            exclude_flags,
            include_any_flags,
            require_flags,
            include_deletions,
            remove_overlaps,
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

pub(crate) fn default_exclude_flags() -> u32 {
    BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP
}

fn parse_flag_value(value: Option<&OsString>, option: &str) -> Result<u32, ()> {
    let Some(raw) = value.and_then(|a| a.to_str()) else {
        print_error("depth", format!("option {option} requires an argument"));
        return Err(());
    };
    let Some(flags) = str_to_flag(raw) else {
        print_error("depth", format!("unknown flag '{}'", raw));
        return Err(());
    };
    Ok(flags as u32)
}

fn read_input_list(path: &Path) -> io::Result<Vec<PathBuf>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut inputs = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end_matches('\r').trim();
        if !trimmed.is_empty() {
            inputs.push(PathBuf::from(trimmed));
        }
    }
    Ok(inputs)
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
    let walk = DepthWalkConfig {
        exclude_flags: config.exclude_flags,
        include_any_flags: config.include_any_flags,
        require_flags: config.require_flags,
        min_mapq: config.min_mapq,
        min_read_len: config.min_read_len,
        min_depth: config.min_depth,
        a_mode: config.a_mode,
        include_deletions: config.include_deletions,
        remove_overlaps: config.remove_overlaps,
    };
    let mut regions = Vec::new();
    if let Some(region) = config.region {
        regions.push(DepthRegion {
            region: parse_region(region)?,
            emit_empty: true,
        });
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
            Exact::Cram => collect_cram_depth(path, config.reference, walk, &regions)?,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "only SAM, BAM, and CRAM input are currently supported",
                ));
            }
        };
        per_input_targets.push(targets);
    }

    emit_depths(
        out,
        inputs,
        &per_input_targets,
        walk.min_depth,
        walk.a_mode,
        config.show_header,
    )
}

fn collect_sam_depth(
    path: &Path,
    config: DepthWalkConfig,
    regions: &[DepthRegion],
) -> io::Result<Vec<DepthTarget>> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(path)?));
    let header = reader.read_header()?;
    let mut targets = depth_targets(&header, regions)?;
    let mut overlap_state = OverlapState::default();

    for result in reader.records() {
        let record = result?;
        update_targets(&header, &mut targets, &record, config, &mut overlap_state);
    }

    Ok(targets)
}

fn collect_bam_depth(
    path: &Path,
    config: DepthWalkConfig,
    regions: &[DepthRegion],
) -> io::Result<Vec<DepthTarget>> {
    let mut reader = bam::io::Reader::new(File::open(path)?);
    let header = reader.read_header()?;
    let mut targets = depth_targets(&header, regions)?;
    let mut overlap_state = OverlapState::default();

    if regions.is_empty() {
        let mut record = bam::Record::default();
        loop {
            let n = reader.read_record(&mut record)?;
            if n == 0 {
                break;
            }
            update_targets(&header, &mut targets, &record, config, &mut overlap_state);
        }
    } else {
        for (i, region) in regions.iter().enumerate() {
            let mut overlap_state = OverlapState::default();
            for record in
                htslib_rs::alignment_compat::query_bam_records_from_path(path, &region.region)?
            {
                update_target(
                    &header,
                    &mut targets[i],
                    i,
                    &record,
                    config,
                    &mut overlap_state,
                );
            }
        }
        if config.a_mode == AMode::AllPositions {
            mark_bam_reference_coverage(path, &header, &mut targets, config)?;
        }
    }

    Ok(targets)
}

fn collect_cram_depth(
    path: &Path,
    reference: Option<&Path>,
    config: DepthWalkConfig,
    regions: &[DepthRegion],
) -> io::Result<Vec<DepthTarget>> {
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(path)?;
    let mut targets = depth_targets(&header, regions)?;

    if regions.is_empty() {
        let mut overlap_state = OverlapState::default();
        for target in &mut targets {
            let region = target_region(target)?;
            for record in query_cram_depth_records(path, &region, reference)? {
                update_target(&header, target, 0, &record, config, &mut overlap_state);
            }
        }
    } else {
        for (i, region) in regions.iter().enumerate() {
            let mut overlap_state = OverlapState::default();
            for record in query_cram_depth_records(path, &region.region, reference)? {
                update_target(
                    &header,
                    &mut targets[i],
                    i,
                    &record,
                    config,
                    &mut overlap_state,
                );
            }
        }
    }

    Ok(targets)
}

fn query_cram_depth_records(
    path: &Path,
    region: &Region,
    reference: Option<&Path>,
) -> io::Result<Vec<sam::alignment::RecordBuf>> {
    if let Some(reference) = reference {
        htslib_rs::alignment_compat::query_cram_records_from_path_with_reference(
            path, region, reference,
        )
    } else {
        htslib_rs::alignment_compat::query_cram_records_from_path_synthesizing_reference(
            path, region,
        )
    }
}

fn mark_bam_reference_coverage(
    path: &Path,
    header: &sam::Header,
    targets: &mut [DepthTarget],
    config: DepthWalkConfig,
) -> io::Result<()> {
    let mut reader = bam::io::Reader::new(File::open(path)?);
    let _ = reader.read_header()?;
    let mut record = bam::Record::default();

    loop {
        let n = reader.read_record(&mut record)?;
        if n == 0 {
            break;
        }
        mark_reference_coverage(header, targets, &record, config);
    }

    Ok(())
}

fn mark_reference_coverage(
    header: &sam::Header,
    targets: &mut [DepthTarget],
    record: &(impl sam::alignment::Record + ?Sized),
    config: DepthWalkConfig,
) {
    let flag = match record.flags() {
        Ok(flags) => u16::from(flags) as u32,
        Err(_) => return,
    };
    if !flag_passes(flag, config) {
        return;
    }
    let mapq = match record.mapping_quality() {
        Some(Ok(q)) => u8::from(q),
        Some(Err(_)) => return,
        None => 0,
    };
    if mapq < config.min_mapq {
        return;
    }
    if config.min_read_len != 0
        && read_length_used(record.cigar().iter()).unwrap_or_default() < config.min_read_len
    {
        return;
    }
    let tid = match record.reference_sequence_id(header).transpose() {
        Ok(Some(t)) => t,
        _ => return,
    };
    if record
        .alignment_start()
        .transpose()
        .ok()
        .flatten()
        .is_none()
    {
        return;
    }
    for target in targets.iter_mut().filter(|target| target.tid == tid) {
        target.reference_has_coverage = true;
    }
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

fn load_bed_regions(path: &Path) -> io::Result<Vec<DepthRegion>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut regions = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let Some(interval) = parse_bed_line(&line) else {
            continue;
        };
        let region = interval.to_region_string();
        regions.push(DepthRegion {
            region: region.parse::<Region>().map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("region \"{}\": {}", region, e),
                )
            })?,
            emit_empty: false,
        });
    }

    Ok(regions)
}

fn emit_depths(
    out: &mut dyn Write,
    inputs: &[PathBuf],
    per_input_targets: &[Vec<DepthTarget>],
    min_depth: u32,
    a_mode: AMode,
    show_header: bool,
) -> io::Result<()> {
    let Some(first_targets) = per_input_targets.first() else {
        return Ok(());
    };

    for targets in &per_input_targets[1..] {
        ensure_compatible_targets(first_targets, targets)?;
    }

    if show_header {
        write!(out, "#CHROM\tPOS")?;
        for input in inputs {
            write!(out, "\t{}", input.display())?;
        }
        writeln!(out)?;
    }

    for (target_index, target) in first_targets.iter().enumerate() {
        let has_any = per_input_targets
            .iter()
            .any(|targets| !targets[target_index].depths.is_empty());
        let reference_has_coverage = per_input_targets
            .iter()
            .any(|targets| targets[target_index].reference_has_coverage);
        let should_emit_empty_span = target.emit_empty
            || matches!(a_mode, AMode::AllRefsAllPositions)
            || (a_mode == AMode::AllPositions && reference_has_coverage);
        if !has_any && !should_emit_empty_span {
            continue;
        }

        let mut emit_offset = |off: usize| -> io::Result<()> {
            let depths: Vec<u32> = per_input_targets
                .iter()
                .map(|targets| targets[target_index].depths.get(off))
                .collect();
            let has_zero_mark = per_input_targets
                .iter()
                .any(|targets| targets[target_index].depths.is_zero_marked(off));
            if a_mode == AMode::None && !depths.iter().any(|&d| d >= min_depth) && !has_zero_mark {
                return Ok(());
            }
            write!(out, "{}\t{}", target.name, target.output_start + off)?;
            for d in &depths {
                write!(out, "\t{d}")?;
            }
            writeln!(out)
        };

        if a_mode == AMode::None {
            // Only covered positions matter; iterate the sorted union of
            // covered offsets across inputs (sparse — safe for huge refs).
            let mut offsets: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
            for targets in per_input_targets {
                offsets.extend(targets[target_index].depths.covered());
            }
            for off in offsets {
                emit_offset(off)?;
            }
        } else {
            // `-a`/`-aa`: every position in the span.
            for off in 0..target.depths.span() {
                emit_offset(off)?;
            }
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
            || left.depths.span() != right.depths.span()
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "depth inputs have incompatible reference dictionaries",
            ));
        }
    }
    Ok(())
}

/// Sparse per-position depth. A dense `Vec<u32>` over the reference span
/// is infeasible for very large references (e.g. the upstream `large_pos`
/// fixture has `LN:10001009800`, which would need ~40 GB). Only covered
/// positions are stored.
#[derive(Default)]
struct Depths {
    span: usize,
    map: std::collections::BTreeMap<usize, u32>,
}

impl Depths {
    fn new(span: usize) -> Self {
        Self {
            span,
            map: std::collections::BTreeMap::new(),
        }
    }
    fn add(&mut self, offset: usize) {
        if offset < self.span {
            let e = self.map.entry(offset).or_insert(0);
            *e = e.saturating_add(1);
        }
    }
    fn mark(&mut self, offset: usize) {
        if offset < self.span {
            self.map.entry(offset).or_insert(0);
        }
    }
    fn get(&self, offset: usize) -> u32 {
        self.map.get(&offset).copied().unwrap_or(0)
    }
    fn is_zero_marked(&self, offset: usize) -> bool {
        self.map.get(&offset).copied() == Some(0)
    }
    fn span(&self) -> usize {
        self.span
    }
    fn is_empty(&self) -> bool {
        self.map.is_empty()
    }
    fn covered(&self) -> impl Iterator<Item = usize> + '_ {
        self.map.keys().copied()
    }
}

struct DepthTarget {
    tid: usize,
    name: String,
    output_start: usize,
    start0: usize,
    end0: usize,
    emit_empty: bool,
    reference_has_coverage: bool,
    depths: Depths,
}

fn depth_targets(
    header: &htslib_rs::sam::Header,
    regions: &[DepthRegion],
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
                    emit_empty: false,
                    reference_has_coverage: false,
                    depths: Depths::new(length),
                }
            })
            .collect())
    } else {
        regions
            .iter()
            .map(|region| depth_target_for_region(header, &region.region, region.emit_empty))
            .collect()
    }
}

fn depth_target_for_region(
    header: &htslib_rs::sam::Header,
    region: &Region,
    emit_empty: bool,
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
        emit_empty,
        reference_has_coverage: false,
        depths: Depths::new(length),
    })
}

fn update_targets(
    header: &sam::Header,
    targets: &mut [DepthTarget],
    record: &(impl sam::alignment::Record + ?Sized),
    config: DepthWalkConfig,
    overlap_state: &mut OverlapState,
) {
    let flag = match record.flags() {
        Ok(flags) => u16::from(flags) as u32,
        Err(_) => return,
    };
    if !flag_passes(flag, config) {
        return;
    }
    let mapq = match record.mapping_quality() {
        Some(Ok(q)) => u8::from(q),
        Some(Err(_)) => return,
        None => 0,
    };
    if mapq < config.min_mapq {
        return;
    }
    if config.min_read_len != 0
        && read_length_used(record.cigar().iter()).unwrap_or_default() < config.min_read_len
    {
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
        target.reference_has_coverage = true;
    }
    for (target_index, target) in targets
        .iter_mut()
        .enumerate()
        .filter(|(_, target)| target.tid == tid)
    {
        update_target_cigar(
            target,
            target_index,
            record,
            flag,
            start,
            config,
            overlap_state,
        );
    }
}

fn update_target(
    header: &sam::Header,
    target: &mut DepthTarget,
    target_index: usize,
    record: &(impl sam::alignment::Record + ?Sized),
    config: DepthWalkConfig,
    overlap_state: &mut OverlapState,
) {
    let flag = match record.flags() {
        Ok(flags) => u16::from(flags) as u32,
        Err(_) => return,
    };
    if !flag_passes(flag, config) {
        return;
    }
    let mapq = match record.mapping_quality() {
        Some(Ok(q)) => u8::from(q),
        Some(Err(_)) => return,
        None => 0,
    };
    if mapq < config.min_mapq {
        return;
    }
    if config.min_read_len != 0
        && read_length_used(record.cigar().iter()).unwrap_or_default() < config.min_read_len
    {
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
    target.reference_has_coverage = true;
    update_target_cigar(
        target,
        target_index,
        record,
        flag,
        start,
        config,
        overlap_state,
    );
}

fn flag_passes(flag: u32, config: DepthWalkConfig) -> bool {
    if flag & config.exclude_flags != 0 {
        return false;
    }
    if config.include_any_flags != 0 && flag & config.include_any_flags == 0 {
        return false;
    }
    if config.require_flags != 0 && flag & config.require_flags != config.require_flags {
        return false;
    }
    true
}

fn read_length_used(
    cigar: impl Iterator<Item = io::Result<htslib_rs::sam::alignment::record::cigar::Op>>,
) -> io::Result<usize> {
    use htslib_rs::sam::alignment::record::cigar::op::Kind;

    let mut len = 0usize;
    for op in cigar {
        let op = op?;
        match op.kind() {
            Kind::Match | Kind::Insertion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                len = len.saturating_add(op.len());
            }
            Kind::Deletion | Kind::Skip | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
    Ok(len)
}

#[derive(Default)]
struct OverlapState {
    seen: std::collections::HashMap<OverlapKey, std::collections::BTreeSet<usize>>,
}

type OverlapKey = (usize, usize, Vec<u8>);

#[derive(Default)]
struct DepthOffsets {
    count: Vec<usize>,
    touch: Vec<usize>,
}

fn update_target_cigar(
    target: &mut DepthTarget,
    target_index: usize,
    record: &(impl sam::alignment::Record + ?Sized),
    flag: u32,
    start: usize,
    config: DepthWalkConfig,
    overlap_state: &mut OverlapState,
) {
    let positions = target_depth_offsets(target, record, start, config.include_deletions);

    if config.remove_overlaps
        && flag & BAM_FPAIRED != 0
        && let Some(name) = record.name()
    {
        let key = (target_index, target.tid, name.to_vec());
        let seen = overlap_state.seen.entry(key).or_default();
        for abs_pos in positions.count {
            if seen.insert(abs_pos) {
                target.depths.add(abs_pos - target.start0);
            }
        }
        for abs_pos in positions.touch {
            if seen.insert(abs_pos) {
                target.depths.mark(abs_pos - target.start0);
            }
        }
    } else {
        for abs_pos in positions.count {
            target.depths.add(abs_pos - target.start0);
        }
        for abs_pos in positions.touch {
            target.depths.mark(abs_pos - target.start0);
        }
    }
}

fn target_depth_offsets(
    target: &DepthTarget,
    record: &(impl sam::alignment::Record + ?Sized),
    start: usize,
    include_deletions: bool,
) -> DepthOffsets {
    let mut positions = DepthOffsets::default();
    let mut ref_pos = start;
    for op in record.cigar().iter() {
        let op = match op {
            Ok(op) => op,
            Err(_) => return positions,
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
                        positions.count.push(p);
                    }
                }
                ref_pos = op_end;
            }
            Kind::Deletion => {
                let op_end = ref_pos.saturating_add(len);
                let lo = ref_pos.max(target.start0);
                let hi = op_end.min(target.end0);
                if hi > lo {
                    let offsets = if include_deletions {
                        &mut positions.count
                    } else {
                        &mut positions.touch
                    };
                    for p in lo..hi {
                        offsets.push(p);
                    }
                }
                ref_pos = op_end;
            }
            Kind::Skip => {
                ref_pos = ref_pos.saturating_add(len);
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
    positions
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
    writeln!(
        w,
        "  -g FLAGS    remove FLAGS from the default filter-out set"
    )?;
    writeln!(w, "  -G FLAGS    add FLAGS to the filter-out set")?;
    Ok(())
}
