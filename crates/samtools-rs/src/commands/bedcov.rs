//! `samtools bedcov` — read depth per BED region.
//!
//! Mirrors `main_bedcov` in `bedcov.c`. The upstream implementation uses
//! pileup to compute exact per-base coverage with filtering. This Rust
//! port computes a simpler per-region **total aligned-base coverage**
//! by walking each record's CIGAR ops over the region — sufficient for
//! a "sum of bases mapped within region" metric. The output per BED line
//! is the original line followed by one coverage column per BAM input.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::format::Exact;
use htslib_rs::sam;

use crate::bam_flag::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FUNMAP, str_to_flag};
use crate::bedidx::parse_bed_line;
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools bedcov`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut min_mapq: u8 = 0;
    let mut print_header = false;
    let mut print_read_count = false;
    let mut min_depth: Option<u32> = None;
    let mut exclude_flags = default_exclude_flags();
    let mut skip_deletions = false;
    let mut positionals: Vec<PathBuf> = Vec::new();

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-Q" | "--min-MQ" | "--min-mq" => {
                min_mapq = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0);
            }
            "-H" => {
                print_header = true;
            }
            "-c" => {
                print_read_count = true;
            }
            "-d" => {
                min_depth = Some(
                    iter.next()
                        .and_then(|a| a.to_str())
                        .and_then(|v| v.parse().ok())
                        .unwrap_or(0),
                );
            }
            "-g" => match parse_flag_value(iter.next(), "-g") {
                Ok(flags) => exclude_flags &= !flags,
                Err(()) => return ExitCode::from(1),
            },
            "-G" => match parse_flag_value(iter.next(), "-G") {
                Ok(flags) => exclude_flags |= flags,
                Err(()) => return ExitCode::from(1),
            },
            "-j" => {
                skip_deletions = true;
            }
            "-X" | "--max-depth" => {
                // Filters/columns not yet supported; consume value-taking flags.
                if matches!(s, "--max-depth") {
                    let _ = iter.next();
                }
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            // Attached short-option values (`-g512`, `-G2048`, `-Q20`,
            // `-d5`) as upstream / option-grouping callers pass them.
            _ if s.starts_with("-g") && s.len() > 2 => {
                let v = OsString::from(&s[2..]);
                match parse_flag_value(Some(&v), "-g") {
                    Ok(flags) => exclude_flags &= !flags,
                    Err(()) => return ExitCode::from(1),
                }
            }
            _ if s.starts_with("-G") && s.len() > 2 => {
                let v = OsString::from(&s[2..]);
                match parse_flag_value(Some(&v), "-G") {
                    Ok(flags) => exclude_flags |= flags,
                    Err(()) => return ExitCode::from(1),
                }
            }
            _ if s.starts_with("-Q") && s.len() > 2 => {
                min_mapq = s[2..].parse().unwrap_or(0);
            }
            _ if s.starts_with("-d") && s.len() > 2 => {
                min_depth = Some(s[2..].parse().unwrap_or(0));
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("bedcov", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => positionals.push(PathBuf::from(arg)),
        }
    }

    if positionals.len() < 2 {
        let _ = print_usage();
        return ExitCode::from(1);
    }
    let bed = positionals[0].clone();
    let bams: Vec<PathBuf> = positionals[1..].to_vec();

    let mut has_cram = false;
    for path in &bams {
        let format = match sam_io::sam_open_format(path) {
            Ok(f) => f,
            Err(e) => {
                print_error("bedcov", e.to_string());
                return ExitCode::from(1);
            }
        };
        match format.exact {
            Exact::Sam | Exact::Bam => {}
            Exact::Cram => has_cram = true,
            _ => {
                print_error(
                    "bedcov",
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
                print_error("bedcov", "CRAM input requires top-level --reference FILE");
                return ExitCode::from(1);
            }
        }
    } else {
        None
    };

    let opts = BedcovOpts {
        min_mapq,
        exclude_flags,
        skip_deletions,
        print_header,
        print_read_count,
        min_depth,
        reference: reference.as_deref(),
    };
    let mut stdout = io::stdout().lock();
    match run_bedcov(&mut stdout, &bed, &bams, &opts) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("bedcov", "bedcov failed", &e);
            ExitCode::from(1)
        }
    }
}

fn default_exclude_flags() -> u32 {
    BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP
}

fn parse_flag_value(value: Option<&OsString>, option: &str) -> Result<u32, ()> {
    let Some(raw) = value.and_then(|a| a.to_str()) else {
        print_error("bedcov", format!("option {option} requires an argument"));
        return Err(());
    };
    let Some(flags) = str_to_flag(raw) else {
        print_error("bedcov", format!("flag value \"{}\" is not supported", raw));
        return Err(());
    };
    Ok(flags as u32)
}

#[derive(Clone, Copy, Debug, Default)]
struct BedcovOpts<'a> {
    min_mapq: u8,
    exclude_flags: u32,
    skip_deletions: bool,
    print_header: bool,
    print_read_count: bool,
    min_depth: Option<u32>,
    reference: Option<&'a Path>,
}

fn run_bedcov<W>(out: &mut W, bed: &Path, bams: &[PathBuf], opts: &BedcovOpts<'_>) -> io::Result<()>
where
    W: Write + ?Sized,
{
    let bed_file = File::open(bed)?;
    let reader = BufReader::new(bed_file);
    let mut wrote_header = !opts.print_header;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('#') {
            if !wrote_header && trimmed.starts_with("#chrom\t") {
                write_bedcov_header(out, Some(trimmed), 0, bams, opts)?;
                wrote_header = true;
            }
            continue;
        }
        if trimmed.starts_with("track ") || trimmed.starts_with("browser ") {
            continue;
        }

        if !wrote_header {
            let field_count = trimmed.split('\t').count();
            write_bedcov_header(out, None, field_count, bams, opts)?;
            wrote_header = true;
        }

        let Some(interval) = parse_bed_line(trimmed) else {
            continue;
        };
        let (Ok(start), Ok(end)) = (i64::try_from(interval.start), i64::try_from(interval.end))
        else {
            continue;
        };

        write!(out, "{}", trimmed)?;
        let mut depth_counts = Vec::new();
        let mut read_counts = Vec::new();
        for bam_path in bams {
            let metrics = compute_region_metrics(bam_path, &interval.chrom, start, end, opts)?;
            let cov = metrics.coverage;
            write!(out, "\t{}", cov)?;
            if opts.min_depth.is_some() {
                depth_counts.push(metrics.depth_bases);
            }
            if opts.print_read_count {
                read_counts.push(metrics.read_count);
            }
        }
        for depth_count in depth_counts {
            write!(out, "\t{}", depth_count)?;
        }
        for count in read_counts {
            write!(out, "\t{}", count)?;
        }
        writeln!(out)?;
    }
    Ok(())
}

fn write_bedcov_header<W>(
    out: &mut W,
    bed_header: Option<&str>,
    field_count: usize,
    bams: &[PathBuf],
    opts: &BedcovOpts<'_>,
) -> io::Result<()>
where
    W: Write + ?Sized,
{
    if let Some(header) = bed_header {
        write!(out, "{header}")?;
    } else {
        const BED_COLS: [&str; 12] = [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "blockCount",
            "blockSizes",
            "blockStarts",
        ];
        for i in 0..field_count {
            if i > 0 {
                write!(out, "\t")?;
            } else {
                write!(out, "#")?;
            }
            write!(out, "{}", BED_COLS.get(i).copied().unwrap_or("."))?;
        }
    }

    for bam in bams {
        write!(out, "\t{}_cov", bam.display())?;
    }
    if opts.min_depth.is_some() {
        for bam in bams {
            write!(out, "\t{}_depth", bam.display())?;
        }
    }
    if opts.print_read_count {
        for bam in bams {
            write!(out, "\t{}_count", bam.display())?;
        }
    }
    writeln!(out)?;
    Ok(())
}

/// Sum of bases (CIGAR `M`/`=`/`X`) that fall inside `[beg, end)` for each
/// record overlapping the region. Excludes records flagged unmapped,
/// secondary, QCFAIL, or duplicate, and those below `min_mapq`.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct RegionMetrics {
    coverage: i64,
    depth_bases: u64,
    read_count: u64,
}

#[derive(Clone, Copy)]
struct RegionMetricConfig {
    beg: i64,
    end: i64,
    exclude_flags: u32,
    min_mapq: u8,
    skip_deletions: bool,
}

fn compute_region_metrics(
    alignment_path: &Path,
    chrom: &str,
    beg: i64,
    end: i64,
    opts: &BedcovOpts<'_>,
) -> io::Result<RegionMetrics> {
    let region_str = format!("{}:{}-{}", chrom, beg + 1, end);
    let region: htslib_rs::core::Region = region_str.parse().map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid region \"{}\": {}", region_str, e),
        )
    })?;
    let mut depths = if opts.min_depth.is_some() && end > beg {
        Some(vec![0u32; (end - beg) as usize])
    } else {
        None
    };

    let mut metrics = RegionMetrics::default();
    match sam_io::sam_open_format(alignment_path)?.exact {
        Exact::Sam => {
            let mut reader = sam::io::Reader::new(BufReader::new(File::open(alignment_path)?));
            let header = reader.read_header()?;
            for result in reader.records() {
                let rec = result?;
                update_region_metrics(
                    &header,
                    &rec,
                    RegionMetricConfig {
                        beg,
                        end,
                        exclude_flags: opts.exclude_flags,
                        min_mapq: opts.min_mapq,
                        skip_deletions: opts.skip_deletions,
                    },
                    &mut metrics,
                    depths.as_mut(),
                );
            }
        }
        Exact::Bam => {
            let header = htslib_rs::alignment_compat::read_bam_header_from_path(alignment_path)?;
            for rec in
                htslib_rs::alignment_compat::query_bam_records_from_path(alignment_path, &region)?
            {
                update_region_metrics(
                    &header,
                    &rec,
                    RegionMetricConfig {
                        beg,
                        end,
                        exclude_flags: opts.exclude_flags,
                        min_mapq: opts.min_mapq,
                        skip_deletions: opts.skip_deletions,
                    },
                    &mut metrics,
                    depths.as_mut(),
                );
            }
        }
        Exact::Cram => {
            let Some(reference) = opts.reference else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "CRAM input requires top-level --reference FILE",
                ));
            };
            let header = htslib_rs::alignment_compat::read_cram_header_from_path(alignment_path)?;
            for rec in htslib_rs::alignment_compat::query_cram_records_from_path_with_reference(
                alignment_path,
                &region,
                reference,
            )? {
                update_region_metrics(
                    &header,
                    &rec,
                    RegionMetricConfig {
                        beg,
                        end,
                        exclude_flags: opts.exclude_flags,
                        min_mapq: opts.min_mapq,
                        skip_deletions: opts.skip_deletions,
                    },
                    &mut metrics,
                    depths.as_mut(),
                );
            }
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "only SAM, BAM, and reference-backed CRAM input are currently supported",
            ));
        }
    }

    if let Some(min_depth) = opts.min_depth {
        metrics.depth_bases = depths
            .as_deref()
            .unwrap_or(&[])
            .iter()
            .filter(|&&depth| depth >= min_depth)
            .count() as u64;
    }
    Ok(metrics)
}

fn update_region_metrics(
    header: &sam::Header,
    rec: &(impl sam::alignment::Record + ?Sized),
    config: RegionMetricConfig,
    metrics: &mut RegionMetrics,
    mut depths: Option<&mut Vec<u32>>,
) {
    let flag = match rec.flags() {
        Ok(flags) => u16::from(flags) as u32,
        Err(_) => return,
    };
    if flag & config.exclude_flags != 0 {
        return;
    }
    let mapq = match rec.mapping_quality() {
        Some(Ok(q)) => u8::from(q),
        Some(Err(_)) => return,
        None => 0,
    };
    if mapq < config.min_mapq {
        return;
    }
    if rec.reference_sequence_id(header).transpose().is_err() {
        return;
    }

    metrics.read_count += 1;
    let start = match rec.alignment_start().transpose() {
        Ok(Some(p)) => usize::from(p) as i64,
        _ => return,
    };
    let mut ref_pos = start - 1; // 0-based
    for op_result in rec.cigar().iter() {
        let op = match op_result {
            Ok(op) => op,
            Err(_) => break,
        };
        let kind = op.kind();
        let len = op.len() as i64;
        use htslib_rs::sam::alignment::record::cigar::op::Kind;
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // Aligned to reference; count overlap with [beg, end).
                let op_end = ref_pos + len;
                let lo = ref_pos.max(config.beg);
                let hi = op_end.min(config.end);
                if hi > lo {
                    metrics.coverage += hi - lo;
                    if let Some(depths) = depths.as_deref_mut() {
                        let offset_start = (lo - config.beg) as usize;
                        let offset_end = (hi - config.beg) as usize;
                        for depth in &mut depths[offset_start..offset_end] {
                            *depth = depth.saturating_add(1);
                        }
                    }
                }
                ref_pos = op_end;
            }
            Kind::Deletion | Kind::Skip => {
                let op_end = ref_pos + len;
                if !config.skip_deletions {
                    let lo = ref_pos.max(config.beg);
                    let hi = op_end.min(config.end);
                    if hi > lo {
                        metrics.coverage += hi - lo;
                        if let Some(depths) = depths.as_deref_mut() {
                            let offset_start = (lo - config.beg) as usize;
                            let offset_end = (hi - config.beg) as usize;
                            for depth in &mut depths[offset_start..offset_end] {
                                *depth = depth.saturating_add(1);
                            }
                        }
                    }
                }
                ref_pos = op_end;
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(
        w,
        "Usage: samtools bedcov [options] <in.bed> <in1.bam> [<in2.bam>...]"
    )?;
    writeln!(w, "Options:")?;
    writeln!(w, "  -Q INT   minimum mapping quality threshold [0]")?;
    writeln!(w, "  -H       print a header line")?;
    writeln!(w, "  -c       add read-count columns")?;
    writeln!(w, "  -d INT   add bases-at-depth-threshold columns")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::write_bedcov_header;
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn tmp_dir(name: &str) -> PathBuf {
        static NEXT_TMP_ID: AtomicUsize = AtomicUsize::new(0);

        let id = NEXT_TMP_ID.fetch_add(1, Ordering::Relaxed);
        let p = std::env::temp_dir().join(format!(
            "samtools-rs-bedcov-{}-{}-{}",
            name,
            std::process::id(),
            id
        ));
        std::fs::create_dir_all(&p).unwrap();
        p
    }

    fn write_bedcov_to_string(
        bed: &Path,
        sam: &Path,
        exclude_flags: u32,
    ) -> std::io::Result<String> {
        write_bedcov_to_string_with_opts(
            bed,
            sam,
            super::BedcovOpts {
                exclude_flags,
                ..super::BedcovOpts::default()
            },
        )
    }

    fn write_bedcov_to_string_with_opts(
        bed: &Path,
        sam: &Path,
        opts: super::BedcovOpts<'_>,
    ) -> std::io::Result<String> {
        let mut out = Vec::new();
        super::run_bedcov(&mut out, bed, &[sam.to_path_buf()], &opts)?;
        Ok(String::from_utf8(out).unwrap())
    }

    #[test]
    fn writes_bedcov_header_from_existing_bed_header() {
        let mut out = Vec::new();
        let bams = vec![PathBuf::from("a.bam")];
        let opts = super::BedcovOpts::default();

        write_bedcov_header(
            &mut out,
            Some("#chrom\tchromStart\tchromEnd\tT1"),
            0,
            &bams,
            &opts,
        )
        .unwrap();

        assert_eq!(
            String::from_utf8(out).unwrap(),
            "#chrom\tchromStart\tchromEnd\tT1\ta.bam_cov\n"
        );
    }

    #[test]
    fn writes_bedcov_header_from_bed_field_count() {
        let mut out = Vec::new();
        let bams = vec![PathBuf::from("a.bam"), PathBuf::from("b.bam")];
        let opts = super::BedcovOpts::default();

        write_bedcov_header(&mut out, None, 14, &bams, &opts).unwrap();

        assert_eq!(
            String::from_utf8(out).unwrap(),
            "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\t.\t.\ta.bam_cov\tb.bam_cov\n"
        );
    }

    #[test]
    fn writes_bedcov_header_with_read_count_columns() {
        let mut out = Vec::new();
        let bams = vec![PathBuf::from("a.bam")];
        let opts = super::BedcovOpts {
            print_read_count: true,
            ..Default::default()
        };

        write_bedcov_header(&mut out, None, 3, &bams, &opts).unwrap();

        assert_eq!(
            String::from_utf8(out).unwrap(),
            "#chrom\tchromStart\tchromEnd\ta.bam_cov\ta.bam_count\n"
        );
    }

    #[test]
    fn writes_bedcov_header_with_depth_and_count_columns() {
        let mut out = Vec::new();
        let bams = vec![PathBuf::from("a.bam")];
        let opts = super::BedcovOpts {
            print_read_count: true,
            min_depth: Some(3),
            ..Default::default()
        };

        write_bedcov_header(&mut out, None, 3, &bams, &opts).unwrap();

        assert_eq!(
            String::from_utf8(out).unwrap(),
            "#chrom\tchromStart\tchromEnd\ta.bam_cov\ta.bam_depth\ta.bam_count\n"
        );
    }

    #[test]
    fn flag_mask_controls_duplicate_record_filtering() {
        let tmp = tmp_dir("flag-mask");
        let bed = tmp.join("r.bed");
        let sam = tmp.join("in.sam");
        std::fs::write(&bed, "chr1\t0\t4\n").unwrap();
        std::fs::write(
            &sam,
            concat!(
                "@HD\tVN:1.6\n",
                "@SQ\tSN:chr1\tLN:8\n",
                "normal\t0\tchr1\t1\t60\t2M\t*\t0\t0\tAA\tII\n",
                "dup\t1024\tchr1\t2\t60\t2M\t*\t0\t0\tCC\tII\n",
            ),
        )
        .unwrap();

        assert_eq!(
            write_bedcov_to_string(&bed, &sam, super::default_exclude_flags()).unwrap(),
            "chr1\t0\t4\t2\n"
        );
        assert_eq!(
            write_bedcov_to_string(
                &bed,
                &sam,
                super::default_exclude_flags() & !super::BAM_FDUP
            )
            .unwrap(),
            "chr1\t0\t4\t4\n"
        );
    }

    #[test]
    fn deletion_and_refskip_bases_are_counted_unless_j_is_set() {
        let tmp = tmp_dir("skip-deletions");
        let bed = tmp.join("r.bed");
        let sam = tmp.join("in.sam");
        std::fs::write(&bed, "chr1\t0\t6\n").unwrap();
        std::fs::write(
            &sam,
            concat!(
                "@HD\tVN:1.6\n",
                "@SQ\tSN:chr1\tLN:8\n",
                "with_del\t0\tchr1\t1\t60\t2M2D2M\t*\t0\t0\tAAAA\tIIII\n",
                "with_skip\t0\tchr1\t1\t60\t2M2N2M\t*\t0\t0\tCCCC\tIIII\n",
            ),
        )
        .unwrap();

        assert_eq!(
            write_bedcov_to_string_with_opts(
                &bed,
                &sam,
                super::BedcovOpts {
                    exclude_flags: super::default_exclude_flags(),
                    ..super::BedcovOpts::default()
                }
            )
            .unwrap(),
            "chr1\t0\t6\t12\n"
        );
        assert_eq!(
            write_bedcov_to_string_with_opts(
                &bed,
                &sam,
                super::BedcovOpts {
                    exclude_flags: super::default_exclude_flags(),
                    skip_deletions: true,
                    ..super::BedcovOpts::default()
                }
            )
            .unwrap(),
            "chr1\t0\t6\t8\n"
        );
    }
}
