//! Native Rust API wrappers for samtools operations needed by BioScript.

use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::write::GzEncoder;
use htslib_rs::format::Exact;

use crate::bedidx::load_bed_index;
use crate::commands::{
    depth as depth_command, merge as merge_command, quickcheck, sort as sort_command,
};
use crate::io as sam_io;
use crate::tmp_file::{self, TempPath};

/// Per-position depth for one reference coordinate.
#[derive(Clone, Debug, PartialEq)]
pub struct PerBaseDepth {
    pub reference_name: String,
    pub position: usize,
    pub depth: u32,
}

/// Summary statistics for per-position depth values.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthSummary {
    pub mean: f64,
    pub median: f64,
    pub min: u32,
    pub max: u32,
    pub uncovered: usize,
}

/// Builds a BAI index for a BAM file and returns the written index path.
pub fn index<P, Q>(
    input_bam: P,
    output_bai: Option<Q>,
    threads: Option<usize>,
) -> io::Result<PathBuf>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    index_native(input_bam, output_bai, threads)
}

/// Builds a BAI index for a BAM file and returns the written index path.
pub fn index_native<P, Q>(
    input_bam: P,
    output_bai: Option<Q>,
    _threads: Option<usize>,
) -> io::Result<PathBuf>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let input_bam = input_bam.as_ref();
    let output_bai = output_bai
        .as_ref()
        .map(|p| p.as_ref().to_path_buf())
        .unwrap_or_else(|| append_extension(input_bam, "bai"));
    let index = htslib_rs::index_compat::build_bai(input_bam)?;
    htslib_rs::index_compat::write_bai(&output_bai, &index)?;
    Ok(output_bai)
}

/// Converts a BAM file to paired FASTQ outputs.
///
/// When `name_sort` is true, the BAM is first name-sorted through the existing
/// in-memory sorter. Output paths ending in `.gz` are gzip-compressed.
pub fn bam_to_fastq_pair<P, Q, R>(
    input_bam: P,
    fastq_1: Q,
    fastq_2: R,
    other_fastq: Option<&Path>,
    singleton_fastq: Option<&Path>,
    name_sort: bool,
    threads: Option<usize>,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    fastq_native(
        input_bam,
        fastq_1,
        fastq_2,
        other_fastq,
        singleton_fastq,
        name_sort,
        threads,
    )
}

/// Converts a BAM file to paired FASTQ outputs.
pub fn fastq_native<P, Q, R>(
    input_bam: P,
    fastq_1: Q,
    fastq_2: R,
    other_fastq: Option<&Path>,
    singleton_fastq: Option<&Path>,
    name_sort: bool,
    _threads: Option<usize>,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    let input_bam = input_bam.as_ref();
    let sorted_bam = if name_sort {
        let (_file, sorted_bam) =
            tmp_file::create_temp_file("samtools-rs-fastq-name-sort", Some("bam"))?;
        sort_command::run_sort(
            Some(input_bam),
            Some(sorted_bam.path()),
            true,
            None,
            sort_command::OutFmt::Bam,
            false,
            None,
            None,
            true,
            false,
        )?;
        Some(sorted_bam)
    } else {
        None
    };
    let source = sorted_bam.as_ref().map(TempPath::path).unwrap_or(input_bam);
    let split =
        htslib_rs::alignment_compat::view_bam_as_fastq_split_text_from_path_with_flag_filter_and_suffix(
            source,
            0,
            0,
            0,
            false,
        )?;

    write_text_or_gzip(fastq_1.as_ref(), split.read1.as_bytes())?;
    write_text_or_gzip(fastq_2.as_ref(), split.read2.as_bytes())?;

    if let Some(path) = other_fastq {
        write_text_or_gzip(path, b"")?;
    }
    if let Some(path) = singleton_fastq {
        write_text_or_gzip(path, split.singleton.as_bytes())?;
    }

    Ok(())
}

/// Returns per-base BAM depths for a region.
pub fn depth<P, S>(
    input_bam: P,
    region: S,
    include_zero: bool,
    threads: Option<usize>,
) -> io::Result<Vec<PerBaseDepth>>
where
    P: AsRef<Path>,
    S: AsRef<str>,
{
    depth_native(input_bam, region, include_zero, threads)
}

/// Returns per-base BAM depths for a region.
pub fn depth_native<P, S>(
    input_bam: P,
    region: S,
    include_zero: bool,
    _threads: Option<usize>,
) -> io::Result<Vec<PerBaseDepth>>
where
    P: AsRef<Path>,
    S: AsRef<str>,
{
    let mut out = Vec::new();
    let a_mode = if include_zero {
        depth_command::AMode::AllPositions
    } else {
        depth_command::AMode::None
    };
    depth_command::run_depth(
        &[input_bam.as_ref().to_path_buf()],
        &mut out,
        depth_command::DepthRunConfig {
            min_mapq: 0,
            min_depth: 1,
            min_read_len: 0,
            a_mode,
            show_header: false,
            exclude_flags: depth_command::default_exclude_flags(),
            include_any_flags: 0,
            require_flags: 0,
            region: Some(region.as_ref()),
            bed: None,
            reference: None,
        },
    )?;
    parse_depth_output(&out)
}

/// Returns summary depth statistics for a BAM region.
pub fn depth_summary<P, S>(
    input_bam: P,
    region: S,
    include_zero: bool,
    threads: Option<usize>,
) -> io::Result<DepthSummary>
where
    P: AsRef<Path>,
    S: AsRef<str>,
{
    let depths = depth_native(input_bam, region, include_zero, threads)?;
    Ok(summarize_depths(&depths))
}

/// Sorts a BAM file to BAM output.
pub fn sort<P, Q>(
    input_bam: P,
    output_bam: Q,
    by_name: bool,
    threads: Option<usize>,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    sort_native(input_bam, output_bam, by_name, threads)
}

/// Sorts a BAM file to BAM output.
pub fn sort_native<P, Q>(
    input_bam: P,
    output_bam: Q,
    by_name: bool,
    _threads: Option<usize>,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    sort_command::run_sort(
        Some(input_bam.as_ref()),
        Some(output_bam.as_ref()),
        by_name,
        None,
        sort_command::OutFmt::Bam,
        false,
        None,
        None,
        true,
        false,
    )
}

/// Merges BAM files to BAM output.
pub fn merge<P, Q>(
    output_bam: P,
    input_bams: &[Q],
    force: bool,
    threads: Option<usize>,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    merge_native(output_bam, input_bams, force, threads)
}

/// Merges BAM files to BAM output.
pub fn merge_native<P, Q>(
    output_bam: P,
    input_bams: &[Q],
    force: bool,
    _threads: Option<usize>,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let output_bam = output_bam.as_ref();
    if output_bam.exists() && !force {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!(
                "output file \"{}\" exists; pass force = true to overwrite",
                output_bam.display()
            ),
        ));
    }
    let inputs: Vec<PathBuf> = input_bams
        .iter()
        .map(|path| path.as_ref().to_path_buf())
        .collect();
    merge_command::run_merge(
        &inputs,
        Some(output_bam),
        merge_command::MergeOrder::Coordinate,
        merge_command::OutFmt::Bam,
        merge_command::MergeIdMode::default(),
        false,
        None,
        0,
        merge_command::MergeRestriction::None,
    )
}

/// Quickly validates a BAM/CRAM/SAM input.
pub fn quickcheck<P>(input_alignment: P, verbose: bool) -> io::Result<()>
where
    P: AsRef<Path>,
{
    quickcheck_native(input_alignment, verbose)
}

/// Quickly validates a BAM/CRAM/SAM input.
pub fn quickcheck_native<P>(input_alignment: P, verbose: bool) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let mut stderr = Vec::new();
    let state = quickcheck::check_file(
        input_alignment.as_ref(),
        i32::from(verbose),
        !verbose,
        false,
        &mut stderr,
    );
    if state == 0 {
        Ok(())
    } else {
        let message = String::from_utf8_lossy(&stderr);
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("quickcheck failed with status {state}: {message}"),
        ))
    }
}

/// Extracts records with all required flag bits set to BAM output.
pub fn extract_unmapped_pairs<P, Q>(
    input_alignment: P,
    output_bam: Q,
    flag: u16,
    threads: Option<usize>,
    reference_fasta: Option<&Path>,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    extract_unmapped_pairs_native(input_alignment, output_bam, flag, threads, reference_fasta)
}

/// Extracts records with all required flag bits set to BAM output.
pub fn extract_unmapped_pairs_native<P, Q>(
    input_alignment: P,
    output_bam: Q,
    flag: u16,
    _threads: Option<usize>,
    reference_fasta: Option<&Path>,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let input_alignment = input_alignment.as_ref();
    let output = File::create(output_bam)?;
    match sam_io::sam_open_format(input_alignment)?.exact {
        Exact::Bam => {
            htslib_rs::alignment_compat::write_bam_records_with_required_flags_from_path(
                input_alignment,
                flag,
                output,
            )?;
        }
        Exact::Cram => {
            let reference_fasta = reference_fasta.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "CRAM extraction requires reference_fasta",
                )
            })?;
            htslib_rs::alignment_compat::write_cram_records_with_required_flags_as_bam_from_path_with_reference(
                input_alignment,
                reference_fasta,
                flag,
                output,
            )?;
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "extract_unmapped_pairs currently supports BAM and CRAM input",
            ));
        }
    }

    Ok(())
}

/// Writes BAM records overlapping one indexed region to a BAM file.
///
/// `threads` and `reference_fasta` are accepted for API compatibility with
/// samtools-shaped callers. BAM region slicing does not use either yet.
pub fn view_region<P, Q, S>(
    input_bam: P,
    region: S,
    output_bam: Q,
    threads: Option<usize>,
    reference_fasta: Option<&Path>,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    S: AsRef<str>,
{
    view_region_native(input_bam, region, output_bam, threads, reference_fasta)
}

/// Writes BAM records overlapping one indexed region to a BAM file.
pub fn view_region_native<P, Q, S>(
    input_bam: P,
    region: S,
    output_bam: Q,
    _threads: Option<usize>,
    reference_fasta: Option<&Path>,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    S: AsRef<str>,
{
    let input = input_bam.as_ref();
    let region = parse_region(region.as_ref())?;
    let output = File::create(output_bam)?;
    match sam_io::sam_open_format(input)?.exact {
        Exact::Bam => {
            htslib_rs::alignment_compat::write_bam_regions_from_path(input, &[region], output)?;
        }
        Exact::Cram => {
            let reference_fasta = reference_fasta.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "CRAM region slicing requires reference_fasta",
                )
            })?;
            htslib_rs::alignment_compat::write_cram_regions_as_bam_from_path_with_reference(
                input,
                reference_fasta,
                &[region],
                output,
            )?;
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "view_region currently supports BAM and CRAM input",
            ));
        }
    }
    Ok(())
}

/// Writes BAM records overlapping BED intervals to a BAM file.
///
/// BED intervals are converted from 0-based half-open to samtools-style
/// 1-based inclusive regions.
pub fn view_bed<P, Q, R>(
    input_bam: P,
    bed_file: R,
    output_bam: Q,
    threads: Option<usize>,
    reference_fasta: Option<&Path>,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    view_bed_native(input_bam, bed_file, output_bam, threads, reference_fasta)
}

/// Writes BAM records overlapping BED intervals to a BAM file.
pub fn view_bed_native<P, Q, R>(
    input_bam: P,
    bed_file: R,
    output_bam: Q,
    _threads: Option<usize>,
    reference_fasta: Option<&Path>,
) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: AsRef<Path>,
{
    let input = input_bam.as_ref();
    let regions = load_bed_regions(bed_file.as_ref())?;
    let output = File::create(output_bam)?;
    match sam_io::sam_open_format(input)?.exact {
        Exact::Bam => {
            htslib_rs::alignment_compat::write_bam_regions_from_path(input, &regions, output)?;
        }
        Exact::Cram => {
            let reference_fasta = reference_fasta.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "CRAM region slicing requires reference_fasta",
                )
            })?;
            htslib_rs::alignment_compat::write_cram_regions_as_bam_from_path_with_reference(
                input,
                reference_fasta,
                &regions,
                output,
            )?;
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "view_bed currently supports BAM and CRAM input",
            ));
        }
    }
    Ok(())
}

fn parse_region(region: &str) -> io::Result<htslib_rs::core::Region> {
    region.parse().map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("region \"{}\": {}", region, e),
        )
    })
}

fn load_bed_regions(path: &Path) -> io::Result<Vec<htslib_rs::core::Region>> {
    load_bed_index(path)?.to_htslib_regions()
}

fn append_extension(path: &Path, ext: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

fn parse_depth_output(buf: &[u8]) -> io::Result<Vec<PerBaseDepth>> {
    let text =
        std::str::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut out = Vec::new();
    for line in text.lines() {
        let mut fields = line.split('\t');
        let reference_name = fields.next().unwrap_or("").to_string();
        let position = fields
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing depth position"))?
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let depth = fields
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing depth value"))?
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        out.push(PerBaseDepth {
            reference_name,
            position,
            depth,
        });
    }
    Ok(out)
}

fn summarize_depths(depths: &[PerBaseDepth]) -> DepthSummary {
    if depths.is_empty() {
        return DepthSummary {
            mean: 0.0,
            median: 0.0,
            min: 0,
            max: 0,
            uncovered: 0,
        };
    }

    let mut values: Vec<u32> = depths.iter().map(|d| d.depth).collect();
    values.sort_unstable();
    let sum: u64 = values.iter().map(|&d| u64::from(d)).sum();
    let mean = sum as f64 / values.len() as f64;
    let median = if values.len().is_multiple_of(2) {
        let hi = values.len() / 2;
        (f64::from(values[hi - 1]) + f64::from(values[hi])) / 2.0
    } else {
        f64::from(values[values.len() / 2])
    };
    let uncovered = values.iter().filter(|&&d| d == 0).count();

    DepthSummary {
        mean,
        median,
        min: values[0],
        max: values[values.len() - 1],
        uncovered,
    }
}

fn write_text_or_gzip(path: &Path, text: &[u8]) -> io::Result<()> {
    let file = File::create(path)?;
    if path.extension().and_then(|ext| ext.to_str()) == Some("gz") {
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(text)?;
        encoder.finish()?;
    } else {
        let mut file = file;
        file.write_all(text)?;
    }
    Ok(())
}
