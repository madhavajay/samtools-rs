//! `samtools stats` — alignment statistics summary.
//!
//! Mirrors `stats.c` (123k LOC). Upstream produces a large blob with many
//! `SN`/`FFQ`/`LFQ`/`COV`/`GCF`/etc. sections. This initial Rust port emits
//! only the basic `SN` summary numbers that can be computed from
//! `AlignmentRecordSummary` without per-base access or pileup.
//!
//! **Pending:** quality histograms (FFQ/LFQ), insert size distributions
//! (IS), GC content (GCF), coverage histograms (COV), per-cycle stats,
//! BAQ adjustments, reference-based mismatch counts.

use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use htslib_rs::alignment_compat::AlignmentRecordSummary;
use htslib_rs::core::Region;
use htslib_rs::format::Exact;
use htslib_rs::sam;

use crate::bam_flag::{
    BAM_FDUP, BAM_FMUNMAP, BAM_FPAIRED, BAM_FPROPER_PAIR, BAM_FQCFAIL, BAM_FREAD1, BAM_FREAD2,
    BAM_FSECONDARY, BAM_FUNMAP,
};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sam_global::current_global_args;
use crate::version::SAMTOOLS_VERSION;

#[derive(Clone, Copy, Debug, Default)]
struct StatsConfig {
    remove_dups: bool,
}

/// Entry point for `samtools stats`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut target_file: Option<PathBuf> = None;
    let mut config = StatsConfig::default();
    let mut regions: Vec<String> = Vec::new();
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "-t" | "--target-regions" => {
                target_file = iter.next().map(PathBuf::from);
            }
            "-@" | "--threads" | "-r" | "--reference" | "-l" | "--read-length" | "-I" | "--id"
            | "-S" | "--split" | "-P" | "--split-prefix" | "-g" | "--cov-threshold" | "-G" => {
                let _ = iter.next();
            }
            "-d" | "--remove-dups" => {
                config.remove_dups = true;
            }
            "-s" | "--sparse" | "-x" | "--sam" | "-p" | "--remove-overlaps" | "--no-PG" => {
                // Accepted but not yet implemented.
            }
            "--help" | "-h" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("stats", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                } else {
                    regions.push(s.to_owned());
                }
            }
        }
    }

    let Some(input) = input else {
        let _ = print_usage();
        return ExitCode::from(1);
    };

    let format = match sam_io::sam_open_format(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error("stats", e.to_string());
            return ExitCode::from(1);
        }
    };

    let target_regions = match target_file.as_deref() {
        Some(path) => match read_target_regions(path) {
            Ok(regions) => regions,
            Err(e) => {
                print_error_errno("stats", "invalid target file", &e);
                return ExitCode::from(1);
            }
        },
        None => Vec::new(),
    };

    let parsed_regions = match regions
        .iter()
        .map(|region| parse_region(region))
        .chain(target_regions.into_iter().map(Ok))
        .collect::<io::Result<Vec<_>>>()
    {
        Ok(regions) => regions,
        Err(e) => {
            print_error_errno("stats", "invalid region", &e);
            return ExitCode::from(1);
        }
    };

    enum StatsInput {
        Summaries(Vec<AlignmentRecordSummary>),
        Counts(StatsCounts),
    }

    let stats_input = match format.exact {
        Exact::Sam if parsed_regions.is_empty() => {
            htslib_rs::alignment_compat::summarize_sam_records_from_path(&input)
                .map(StatsInput::Summaries)
        }
        Exact::Sam => {
            collect_sam_region_stats(&input, &parsed_regions, config).map(StatsInput::Counts)
        }
        Exact::Bam if parsed_regions.is_empty() => {
            htslib_rs::alignment_compat::summarize_bam_records_from_path(&input)
                .map(StatsInput::Summaries)
        }
        Exact::Bam => {
            collect_bam_region_stats(&input, &parsed_regions, config).map(StatsInput::Counts)
        }
        Exact::Cram => {
            let Some(reference) = current_global_args().reference else {
                print_error("stats", "CRAM input requires top-level --reference FILE");
                return ExitCode::from(1);
            };
            if parsed_regions.is_empty() {
                htslib_rs::alignment_compat::summarize_cram_records_from_path_with_reference(
                    &input, reference,
                )
                .map(StatsInput::Summaries)
            } else {
                collect_cram_region_stats(&input, reference, &parsed_regions, config)
                    .map(StatsInput::Counts)
            }
        }
        _ => {
            print_error(
                "stats",
                "only SAM, BAM, and CRAM input are currently supported",
            );
            return ExitCode::from(1);
        }
    };
    let stats_input = match stats_input {
        Ok(v) => v,
        Err(e) => {
            print_error_errno(
                "stats",
                format!("error reading from \"{}\"", input.display()),
                &e,
            );
            return ExitCode::from(1);
        }
    };

    let mut writer = match sam_io::open_text_output(output.as_deref()) {
        Ok(writer) => writer,
        Err(e) => {
            print_error_errno("stats", "open -o output", &e);
            return ExitCode::from(1);
        }
    };
    let write_result = match stats_input {
        StatsInput::Summaries(summaries) => write_stats(&mut writer, &summaries, config),
        StatsInput::Counts(counts) => write_stats_counts(&mut writer, &counts),
    };
    if let Err(e) = write_result {
        if e.kind() == io::ErrorKind::BrokenPipe {
            return ExitCode::SUCCESS;
        }
        print_error_errno("stats", "write failed", &e);
        return ExitCode::from(1);
    }
    match sam_io::check_sam_close(&mut writer) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("stats", "close output", &e);
            ExitCode::from(1)
        }
    }
}

fn parse_region(s: &str) -> io::Result<Region> {
    s.parse::<Region>().map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("region \"{}\": {}", s, e),
        )
    })
}

fn read_target_regions(path: &std::path::Path) -> io::Result<Vec<Region>> {
    let file = File::open(path)?;
    let mut regions = Vec::new();
    for (line_no, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') {
            continue;
        }

        let fields: Vec<_> = s.split_whitespace().collect();
        if fields.len() < 3 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}:{}: expected CHROM START END",
                    path.display(),
                    line_no + 1
                ),
            ));
        }

        let start: u64 = fields[1].parse().map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}:{}: invalid target start \"{}\": {}",
                    path.display(),
                    line_no + 1,
                    fields[1],
                    e
                ),
            )
        })?;
        let end: u64 = fields[2].parse().map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}:{}: invalid target end \"{}\": {}",
                    path.display(),
                    line_no + 1,
                    fields[2],
                    e
                ),
            )
        })?;
        if end < start {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}:{}: target end precedes start",
                    path.display(),
                    line_no + 1
                ),
            ));
        }

        regions.push(parse_region(&format!("{}:{}-{}", fields[0], start, end))?);
    }
    Ok(regions)
}

#[derive(Clone, Debug)]
struct RegionTarget {
    tid: usize,
    start: usize,
    end: usize,
}

fn region_targets(header: &sam::Header, regions: &[Region]) -> io::Result<Vec<RegionTarget>> {
    regions
        .iter()
        .map(|region| {
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
            let start = interval.start().map(usize::from).unwrap_or(1);
            let end = interval.end().map(usize::from).unwrap_or(ref_len);
            if end < start {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid region interval: {}", region),
                ));
            }

            Ok(RegionTarget { tid, start, end })
        })
        .collect()
}

fn collect_sam_region_stats(
    input: &PathBuf,
    regions: &[Region],
    config: StatsConfig,
) -> io::Result<StatsCounts> {
    let mut reader = File::open(input)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader.read_header()?;
    let targets = region_targets(&header, regions)?;
    let mut counts = StatsCounts::default();
    let mut seen = HashSet::new();

    for result in reader.records() {
        let record = result?;
        if record_overlaps_targets(&header, &record, &targets)
            && seen.insert(record_identity(&header, &record))
        {
            counts.update_record(&header, &record, config);
        }
    }

    Ok(counts)
}

fn collect_bam_region_stats(
    input: &PathBuf,
    regions: &[Region],
    config: StatsConfig,
) -> io::Result<StatsCounts> {
    let header = htslib_rs::alignment_compat::read_bam_header_from_path(input)?;
    let mut counts = StatsCounts::default();
    let mut seen = HashSet::new();
    for region in regions {
        for record in htslib_rs::alignment_compat::query_bam_records_from_path(input, region)? {
            if seen.insert(record_identity(&header, &record)) {
                counts.update_record(&header, &record, config);
            }
        }
    }
    Ok(counts)
}

fn collect_cram_region_stats(
    input: &PathBuf,
    reference: PathBuf,
    regions: &[Region],
    config: StatsConfig,
) -> io::Result<StatsCounts> {
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(input)?;
    let mut counts = StatsCounts::default();
    let mut seen = HashSet::new();
    for region in regions {
        for record in htslib_rs::alignment_compat::query_cram_records_from_path_with_reference(
            input, region, &reference,
        )? {
            if seen.insert(record_identity(&header, &record)) {
                counts.update_record(&header, &record, config);
            }
        }
    }
    Ok(counts)
}

fn record_identity(
    header: &sam::Header,
    record: &(impl sam::alignment::Record + ?Sized),
) -> Vec<u8> {
    let name = record.name().map(|name| name.to_vec()).unwrap_or_default();
    let flags = record.flags().map(u16::from).unwrap_or_default();
    let tid = record
        .reference_sequence_id(header)
        .transpose()
        .unwrap_or_default();
    let start = record
        .alignment_start()
        .transpose()
        .unwrap_or_default()
        .map(usize::from);
    let mate_tid = record
        .mate_reference_sequence_id(header)
        .transpose()
        .unwrap_or_default();
    let mate_start = record
        .mate_alignment_start()
        .transpose()
        .unwrap_or_default()
        .map(usize::from);

    format!(
        "{}\t{}\t{:?}\t{:?}\t{:?}\t{:?}\t{:?}",
        String::from_utf8_lossy(&name),
        flags,
        tid,
        start,
        mate_tid,
        mate_start,
        record.template_length().ok()
    )
    .into_bytes()
}

fn record_overlaps_targets(
    header: &sam::Header,
    record: &(impl sam::alignment::Record + ?Sized),
    targets: &[RegionTarget],
) -> bool {
    let tid = match record.reference_sequence_id(header).transpose() {
        Ok(Some(tid)) => tid,
        _ => return false,
    };
    let start = match record.alignment_start().transpose() {
        Ok(Some(start)) => usize::from(start),
        _ => return false,
    };
    let end = match record.alignment_end().transpose() {
        Ok(Some(end)) => usize::from(end),
        _ => return false,
    };

    targets
        .iter()
        .any(|target| target.tid == tid && start <= target.end && target.start <= end)
}

#[derive(Default)]
struct StatsCounts {
    raw_total: u64,
    filtered: u64,
    total: u64,
    mapped: u64,
    unmapped: u64,
    paired: u64,
    proper_paired: u64,
    read1: u64,
    read2: u64,
    mapped_and_paired: u64,
    dup: u64,
    qc_fail: u64,
    secondary: u64,
    mq0: u64,
    singletons: u64,
    diffchr: u64,
}

impl StatsCounts {
    fn update_record(
        &mut self,
        header: &sam::Header,
        rec: &(impl sam::alignment::Record + ?Sized),
        config: StatsConfig,
    ) {
        let Ok(flags) = rec.flags() else {
            return;
        };
        let flag = u16::from(flags) as u32;
        let mapq = rec.mapping_quality().and_then(Result::ok).map(u8::from);
        let reference_sequence_id = rec
            .reference_sequence_id(header)
            .transpose()
            .unwrap_or_default();
        let mate_reference_sequence_id = rec
            .mate_reference_sequence_id(header)
            .transpose()
            .unwrap_or_default();

        self.update(
            flag,
            mapq,
            reference_sequence_id,
            mate_reference_sequence_id,
            config,
        );
    }

    fn update_summary(&mut self, rec: &AlignmentRecordSummary, config: StatsConfig) {
        self.update(
            rec.flags_u16() as u32,
            rec.mapping_quality(),
            rec.reference_sequence_id(),
            rec.mate_reference_sequence_id(),
            config,
        );
    }

    fn update(
        &mut self,
        flag: u32,
        mapq: Option<u8>,
        reference_sequence_id: Option<usize>,
        mate_reference_sequence_id: Option<usize>,
        config: StatsConfig,
    ) {
        if flag & BAM_FSECONDARY != 0 {
            self.secondary += 1;
            return;
        }
        self.raw_total += 1;
        if config.remove_dups && flag & BAM_FDUP != 0 {
            self.filtered += 1;
            return;
        }
        self.total += 1;

        if flag & BAM_FUNMAP == 0 {
            self.mapped += 1;
            if mapq == Some(0) {
                self.mq0 += 1;
            }
        } else {
            self.unmapped += 1;
        }
        if flag & BAM_FPAIRED != 0 {
            self.paired += 1;
            if flag & BAM_FREAD1 != 0 {
                self.read1 += 1;
            }
            if flag & BAM_FREAD2 != 0 {
                self.read2 += 1;
            }
            if flag & BAM_FPROPER_PAIR != 0 {
                self.proper_paired += 1;
            }
            if flag & BAM_FUNMAP == 0 && flag & BAM_FMUNMAP == 0 {
                self.mapped_and_paired += 1;
                if reference_sequence_id != mate_reference_sequence_id
                    && reference_sequence_id.is_some()
                    && mate_reference_sequence_id.is_some()
                {
                    self.diffchr += 1;
                }
            }
            if flag & BAM_FMUNMAP != 0 && flag & BAM_FUNMAP == 0 {
                self.singletons += 1;
            }
        }
        if flag & BAM_FDUP != 0 {
            self.dup += 1;
        }
        if flag & BAM_FQCFAIL != 0 {
            self.qc_fail += 1;
        }
    }
}

fn write_stats(
    out: &mut dyn Write,
    recs: &[AlignmentRecordSummary],
    config: StatsConfig,
) -> io::Result<()> {
    let mut counts = StatsCounts::default();
    for rec in recs {
        counts.update_summary(rec, config);
    }
    write_stats_counts(out, &counts)
}

fn write_stats_counts(out: &mut dyn Write, counts: &StatsCounts) -> io::Result<()> {
    writeln!(
        out,
        "# This file was produced by samtools-rs stats (samtools-{}+htslib-rs)",
        SAMTOOLS_VERSION
    )?;
    writeln!(out, "# This file contains statistics for all reads.")?;
    writeln!(out, "SN\traw total sequences:\t{}", counts.raw_total)?;
    writeln!(out, "SN\tfiltered sequences:\t{}", counts.filtered)?;
    writeln!(out, "SN\tsequences:\t{}", counts.total)?;
    writeln!(out, "SN\t1st fragments:\t{}", counts.read1)?;
    writeln!(out, "SN\tlast fragments:\t{}", counts.read2)?;
    writeln!(out, "SN\treads mapped:\t{}", counts.mapped)?;
    writeln!(
        out,
        "SN\treads mapped and paired:\t{}",
        counts.mapped_and_paired
    )?;
    writeln!(out, "SN\treads unmapped:\t{}", counts.unmapped)?;
    writeln!(out, "SN\treads properly paired:\t{}", counts.proper_paired)?;
    writeln!(out, "SN\treads paired:\t{}", counts.paired)?;
    writeln!(out, "SN\treads duplicated:\t{}", counts.dup)?;
    writeln!(out, "SN\treads MQ0:\t{}", counts.mq0)?;
    writeln!(out, "SN\treads QC failed:\t{}", counts.qc_fail)?;
    writeln!(out, "SN\tnon-primary alignments:\t{}", counts.secondary)?;
    writeln!(out, "SN\tsupplementary alignments:\t0")?;
    writeln!(out, "SN\tsingletons:\t{}", counts.singletons)?;
    writeln!(
        out,
        "SN\tpairs on different chromosomes:\t{}",
        counts.diffchr / 2
    )?;
    Ok(())
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools stats [options] <in.bam>")?;
    writeln!(w, "  -o FILE      output FILE")?;
    writeln!(w)?;
    writeln!(
        w,
        "Note: only the `SN` summary lines are currently produced."
    )?;
    Ok(())
}
