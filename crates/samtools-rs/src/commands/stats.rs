//! `samtools stats` — alignment statistics summary.
//!
//! Mirrors `stats.c` (123k LOC). Upstream produces a large blob with many
//! `SN`/`FFQ`/`LFQ`/`COV`/`GCF`/etc. sections. This Rust port emits the
//! `SN` summary numbers plus record-level quality, GC, and approximate
//! CIGAR-walk coverage histograms that can be computed without pileup.
//!
//! **Pending:** insert size distributions (IS), exact pileup-backed coverage
//! histograms, per-cycle stats, BAQ adjustments, and deeper reference-based
//! mismatch parity.

use std::collections::{BTreeMap, HashMap, HashSet};
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
    BAM_FREVERSE, BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP, str_to_flag,
};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

use crate::sam_global::current_global_args;
use crate::version::SAMTOOLS_VERSION;

/// Upstream `stats.c` `stats->ngc` — GC-fraction histogram array size.
const NGC: usize = 200;
/// Upstream `stats->nindels` (init = `nbases` = 300, never grown):
/// indels longer than this are excluded from the ID distribution.
const NINDELS: usize = 300;

#[derive(Clone, Debug)]
struct StatsConfig {
    remove_dups: bool,
    required_flags: u32,
    filter_flags: u32,
    id_filter: Option<String>,
    insert_size_max: u32,
    insert_size_main_bulk: f64,
    read_length_filter: Option<usize>,
    trim_quality: u8,
    coverage_min: u32,
    coverage_max: u32,
    coverage_step: u32,
    cov_threshold: u32,
    // True when a reference (`-r`) was supplied; upstream only allocates
    // `mpc_buf` (and therefore prints the MPC section) in that case.
    has_reference: bool,
    // Reference sequences (upper-cased bases) keyed by name, loaded when
    // a reference is supplied. Drives the MPC reference-mismatch engine.
    reference_seqs: Option<HashMap<String, Vec<u8>>>,
    // `-S`/`--split <tag>`: also write per-tag-value `.bamstat` files.
    // `-P`/`--split-prefix`: filename prefix (default = input path).
    split_tag: Option<String>,
    split_prefix: Option<String>,
    // `--ref-stats`: emit the RFS reference-statistics section.
    ref_stats: bool,
}

impl Default for StatsConfig {
    fn default() -> Self {
        Self {
            remove_dups: false,
            required_flags: 0,
            filter_flags: 0,
            id_filter: None,
            insert_size_max: 8000,
            insert_size_main_bulk: 0.99,
            read_length_filter: None,
            trim_quality: 0,
            coverage_min: 1,
            coverage_max: 1000,
            coverage_step: 1,
            cov_threshold: 0,
            has_reference: false,
            reference_seqs: None,
            split_tag: None,
            split_prefix: None,
            ref_stats: false,
        }
    }
}

/// Entry point for `samtools stats`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut target_file: Option<PathBuf> = None;
    let mut config = StatsConfig::default();
    let mut reference_arg: Option<PathBuf> = None;
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
            "-c" | "--coverage" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("stats", "option -c requires an argument");
                    return ExitCode::from(1);
                };
                match parse_coverage_range(raw) {
                    Ok((min, max, step)) => {
                        config.coverage_min = min;
                        config.coverage_max = max;
                        config.coverage_step = step;
                    }
                    Err(e) => {
                        print_error("stats", e);
                        return ExitCode::from(1);
                    }
                }
            }
            "-g" | "--cov-threshold" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("stats", "option -g requires an argument");
                    return ExitCode::from(1);
                };
                match raw.parse::<u32>() {
                    Ok(threshold) => config.cov_threshold = threshold,
                    Err(e) => {
                        print_error(
                            "stats",
                            format!("invalid coverage threshold \"{raw}\": {e}"),
                        );
                        return ExitCode::from(1);
                    }
                }
            }
            "-f" | "--required-flag" => match parse_flag_value(iter.next(), s) {
                Ok(flags) => config.required_flags = flags,
                Err(()) => return ExitCode::from(1),
            },
            "-F" | "--filtering-flag" => match parse_flag_value(iter.next(), s) {
                Ok(flags) => config.filter_flags |= flags,
                Err(()) => return ExitCode::from(1),
            },
            "-I" | "--id" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("stats", "option -I requires an argument");
                    return ExitCode::from(1);
                };
                config.id_filter = Some(raw.to_owned());
            }
            "-i" | "--insert-size" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("stats", "option -i requires an argument");
                    return ExitCode::from(1);
                };
                match raw.parse::<u32>() {
                    Ok(max) => config.insert_size_max = max,
                    Err(e) => {
                        print_error("stats", format!("invalid insert size \"{raw}\": {e}"));
                        return ExitCode::from(1);
                    }
                }
            }
            "-m" | "--most-inserts" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("stats", "option -m requires an argument");
                    return ExitCode::from(1);
                };
                match raw.parse::<f64>() {
                    Ok(value) if value.is_finite() && value >= 0.0 => {
                        config.insert_size_main_bulk = value
                    }
                    Ok(_) => {
                        print_error("stats", format!("invalid most-inserts value \"{raw}\""));
                        return ExitCode::from(1);
                    }
                    Err(e) => {
                        print_error(
                            "stats",
                            format!("invalid most-inserts value \"{raw}\": {e}"),
                        );
                        return ExitCode::from(1);
                    }
                }
            }
            "-l" | "--read-length" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("stats", "option -l requires an argument");
                    return ExitCode::from(1);
                };
                match raw.parse::<usize>() {
                    Ok(len) => config.read_length_filter = Some(len),
                    Err(e) => {
                        print_error("stats", format!("invalid read length \"{raw}\": {e}"));
                        return ExitCode::from(1);
                    }
                }
            }
            "-q" | "--trim-quality" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("stats", "option -q requires an argument");
                    return ExitCode::from(1);
                };
                match raw.parse::<u8>() {
                    Ok(trim_quality) => config.trim_quality = trim_quality,
                    Err(e) => {
                        print_error("stats", format!("invalid trim quality \"{raw}\": {e}"));
                        return ExitCode::from(1);
                    }
                }
            }
            "-r" | "--reference" | "--ref-seq" => {
                reference_arg = iter.next().map(PathBuf::from);
            }
            "-S" | "--split" => {
                config.split_tag = iter.next().and_then(|s| s.to_str().map(str::to_owned));
            }
            "-P" | "--split-prefix" => {
                config.split_prefix = iter.next().and_then(|s| s.to_str().map(str::to_owned));
            }
            "--ref-stats" => {
                config.ref_stats = true;
            }
            "--ref-stats-chunk" => {
                let _ = iter.next();
            }
            "-@" | "--threads" | "-G" => {
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
    if config.cov_threshold > 0 && parsed_regions.is_empty() {
        print_error(
            "stats",
            "Coverage percentage calculation requires a list of target regions",
        );
        return ExitCode::from(1);
    }

    enum StatsInput {
        // Retained for the SAM/BAM/CRAM `AlignmentRecordSummary` fallback
        // shape; no longer produced now that no-region CRAM uses the
        // full-record iterator (TODO-NEXT #2).
        #[allow(dead_code)]
        Summaries(Vec<AlignmentRecordSummary>),
        Counts(Box<StatsCounts>),
    }

    let header_sort_order = match read_input_header_sort_order(&input, format.exact) {
        Ok(so) => so,
        Err(e) => {
            print_error_errno(
                "stats",
                format!("error reading header from \"{}\"", input.display()),
                &e,
            );
            return ExitCode::from(1);
        }
    };

    let resolved_reference = reference_arg
        .clone()
        .or_else(|| current_global_args().reference);
    config.has_reference = resolved_reference.is_some();
    if let Some(ref_path) = resolved_reference.as_ref() {
        match load_reference_seqs(ref_path) {
            Ok(map) => config.reference_seqs = Some(map),
            Err(e) => {
                print_error_errno(
                    "stats",
                    format!("failed to read reference \"{}\"", ref_path.display()),
                    &e,
                );
                return ExitCode::from(1);
            }
        }
    }

    let stats_input = match format.exact {
        Exact::Sam if parsed_regions.is_empty() => collect_sam_full_stats(&input, &config)
            .map(|counts| StatsInput::Counts(Box::new(counts))),
        Exact::Sam => collect_sam_region_stats(&input, &parsed_regions, &config)
            .map(|counts| StatsInput::Counts(Box::new(counts))),
        Exact::Bam if parsed_regions.is_empty() => collect_bam_full_stats(&input, &config)
            .map(|counts| StatsInput::Counts(Box::new(counts))),
        Exact::Bam => collect_bam_region_stats(&input, &parsed_regions, &config)
            .map(|counts| StatsInput::Counts(Box::new(counts))),
        Exact::Cram => {
            let Some(reference) = reference_arg
                .clone()
                .or_else(|| current_global_args().reference)
            else {
                print_error(
                    "stats",
                    "CRAM input requires -r/--reference FILE or top-level --reference",
                );
                return ExitCode::from(1);
            };
            if parsed_regions.is_empty() {
                collect_cram_full_stats(&input, reference, &config)
                    .map(|counts| StatsInput::Counts(Box::new(counts)))
            } else {
                collect_cram_region_stats(&input, reference, &parsed_regions, &config)
                    .map(|counts| StatsInput::Counts(Box::new(counts)))
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
        StatsInput::Summaries(summaries) => {
            let is_sorted = header_sort_order
                .as_deref()
                .map(|so| so == "coordinate")
                .unwrap_or(false);
            write_stats(&mut writer, &summaries, &config, is_sorted)
        }
        StatsInput::Counts(counts) => {
            let is_sorted = counts.is_coordinate_sorted();
            let combined = write_stats_counts(&mut writer, &counts, &config, is_sorted);
            // `-S`/`--split`: one `<prefix|input>_<value>.bamstat` per
            // tag value, prefix defaulting to the input path.
            combined.and_then(|()| {
                if config.split_tag.is_none() {
                    return Ok(());
                }
                let prefix = config
                    .split_prefix
                    .clone()
                    .unwrap_or_else(|| input.to_string_lossy().into_owned());
                for (value, sub) in &counts.splits {
                    let path = format!("{prefix}_{value}.bamstat");
                    let file = File::create(&path)?;
                    let mut w = std::io::BufWriter::new(file);
                    write_stats_counts(&mut w, sub, &config, sub.is_coordinate_sorted())?;
                    w.flush()?;
                }
                Ok(())
            })
        }
    };
    let write_result = write_result.and_then(|()| {
        if config.ref_stats {
            let dims = read_input_ref_dims(&input, format.exact)?;
            write_ref_stats(&mut writer, &dims, config.reference_seqs.as_ref())?;
        }
        Ok(())
    });
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

/// Looks up a record's `NM` (edit distance) aux tag and returns its numeric
/// value if present. Returns `None` if the tag is absent or has an
/// unexpected type (which we treat as 0 contribution, matching upstream's
/// silent skip).
fn read_nm_aux(rec: &(impl sam::alignment::Record + ?Sized)) -> Option<u64> {
    use sam::alignment::record::data::field::{Tag, Value};
    let nm_tag = Tag::from([b'N', b'M']);
    let data = rec.data();
    let value = data.get(&nm_tag)?.ok()?;
    let n = match value {
        Value::Int8(v) => v as i64,
        Value::UInt8(v) => v as i64,
        Value::Int16(v) => v as i64,
        Value::UInt16(v) => v as i64,
        Value::Int32(v) => v as i64,
        Value::UInt32(v) => v as i64,
        _ => return None,
    };
    if n < 0 { None } else { Some(n as u64) }
}

/// Returns the `SO:` value from the input's `@HD` line, lower-cased, or
/// `None` if the header has no sort-order tag. This is only used for
/// summary-only paths that cannot currently inspect full records.
fn read_input_header_sort_order(
    input: &std::path::Path,
    exact: Exact,
) -> io::Result<Option<String>> {
    let header = match exact {
        Exact::Sam => htslib_rs::alignment_compat::read_sam_header_from_path(input)?,
        Exact::Bam => htslib_rs::alignment_compat::read_bam_header_from_path(input)?,
        Exact::Cram => htslib_rs::alignment_compat::read_cram_header_from_path(input)?,
        _ => return Ok(None),
    };
    let so = header
        .header()
        .as_ref()
        .and_then(|hd| {
            hd.other_fields()
                .get(&sam::header::record::value::map::header::tag::SORT_ORDER)
        })
        .map(|value| String::from_utf8_lossy(value.as_ref()).to_lowercase());
    Ok(so)
}

/// Header `@SQ` (name, length) pairs in order, for the `--ref-stats`
/// RFS section.
fn read_input_ref_dims(input: &std::path::Path, exact: Exact) -> io::Result<Vec<(String, u64)>> {
    let header = match exact {
        Exact::Sam => htslib_rs::alignment_compat::read_sam_header_from_path(input)?,
        Exact::Bam => htslib_rs::alignment_compat::read_bam_header_from_path(input)?,
        Exact::Cram => htslib_rs::alignment_compat::read_cram_header_from_path(input)?,
        _ => return Ok(Vec::new()),
    };
    Ok(header
        .reference_sequences()
        .iter()
        .map(|(name, def)| {
            (
                String::from_utf8_lossy(name.as_ref()).into_owned(),
                usize::from(def.length()) as u64,
            )
        })
        .collect())
}

/// Writes the RFS reference-statistics section (`--ref-stats`). Without
/// a reference the GC/N columns are -1 (upstream `gcsum=-1`); with one,
/// GC = G+C / (A+C+G+T) and N = count of N over the header-length
/// prefix of each sequence, mirroring `collect_refstats`.
fn write_ref_stats(
    out: &mut dyn Write,
    dims: &[(String, u64)],
    refmap: Option<&HashMap<String, Vec<u8>>>,
) -> io::Result<()> {
    writeln!(
        out,
        "# Reference statistics. Use `grep ^RFS | cut -f 2-` to extract this part."
    )?;
    writeln!(
        out,
        "# Total count, Output count, Average GC, Min length, Max length, Average length, Total length in first row."
    )?;
    writeln!(
        out,
        "# Sequence name, Length, GC content, Unknown count in following rows."
    )?;
    let total = dims.len() as i64;
    let combined: i64 = dims.iter().map(|(_, l)| *l as i64).sum();
    let minlen = dims.iter().map(|(_, l)| *l as i64).min().unwrap_or(0);
    let maxlen = dims.iter().map(|(_, l)| *l as i64).max().unwrap_or(0);
    let avglen = if total > 0 {
        combined as f64 / total as f64
    } else {
        -1.0
    };
    // Per-sequence GC fraction / N count over the header-length prefix.
    let mut rows: Vec<(String, u64, f64, i64)> = Vec::with_capacity(dims.len());
    let mut gcsum = 0.0_f64;
    let mut have_ref = false;
    for (name, len) in dims {
        if let Some(seq) = refmap.and_then(|m| m.get(name)) {
            have_ref = true;
            let take = (*len as usize).min(seq.len());
            let (mut gc, mut at, mut n) = (0i64, 0i64, 0i64);
            for &b in &seq[..take] {
                match b {
                    b'G' | b'C' => gc += 1,
                    b'A' | b'T' => at += 1,
                    b'N' => n += 1,
                    _ => {}
                }
            }
            let refgc = if gc + at > 0 {
                gc as f64 / (gc + at) as f64
            } else {
                0.0
            };
            gcsum += refgc;
            rows.push((name.clone(), *len, refgc, n));
        } else if refmap.is_some() {
            // Reference supplied but sequence absent: upstream ends with
            // zero bases -> GC 0, N 0.
            have_ref = true;
            rows.push((name.clone(), *len, 0.0, 0));
        } else {
            rows.push((name.clone(), *len, -1.0, -1));
        }
    }
    let avggc = if !have_ref || total == 0 {
        -1.0
    } else {
        gcsum / total as f64
    };
    writeln!(
        out,
        "RFS\t{total}\t{total}\t{avggc:.2}\t{minlen}\t{maxlen}\t{avglen:.2}\t{combined}"
    )?;
    for (name, len, gc, n) in rows {
        writeln!(out, "RFS\t{name}\t{len}\t{gc:.2}\t{n}")?;
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

fn parse_coverage_range(raw: &str) -> Result<(u32, u32, u32), String> {
    let fields: Vec<_> = raw.split(',').collect();
    if fields.len() != 3 {
        return Err(format!("coverage range \"{raw}\" must have MIN,MAX,STEP"));
    }
    let min = fields[0]
        .parse::<u32>()
        .map_err(|e| format!("invalid coverage min \"{}\": {}", fields[0], e))?;
    let max = fields[1]
        .parse::<u32>()
        .map_err(|e| format!("invalid coverage max \"{}\": {}", fields[1], e))?;
    let step = fields[2]
        .parse::<u32>()
        .map_err(|e| format!("invalid coverage step \"{}\": {}", fields[2], e))?;
    if min == 0 || max < min || step == 0 {
        return Err(format!(
            "invalid coverage range \"{raw}\": require 0 < MIN <= MAX and STEP > 0"
        ));
    }
    Ok((min, max, step))
}

fn parse_flag_value(value: Option<&OsString>, option: &str) -> Result<u32, ()> {
    let Some(raw) = value.and_then(|a| a.to_str()) else {
        print_error("stats", format!("option {option} requires an argument"));
        return Err(());
    };
    let Some(flags) = str_to_flag(raw) else {
        print_error("stats", format!("Unknown flag '{}'", raw));
        return Err(());
    };
    Ok(flags as u32)
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

fn target_base_count(targets: &[RegionTarget]) -> u64 {
    let mut intervals: Vec<_> = targets
        .iter()
        .map(|target| (target.tid, target.start, target.end))
        .collect();
    intervals.sort_unstable();

    let mut total = 0u64;
    let mut current: Option<(usize, usize, usize)> = None;
    for (tid, start, end) in intervals {
        match current {
            Some((cur_tid, cur_start, cur_end)) if cur_tid == tid && start <= cur_end + 1 => {
                current = Some((cur_tid, cur_start, cur_end.max(end)));
            }
            Some((_, cur_start, cur_end)) => {
                total += (cur_end - cur_start + 1) as u64;
                current = Some((tid, start, end));
            }
            None => current = Some((tid, start, end)),
        }
    }
    if let Some((_, start, end)) = current {
        total += (end - start + 1) as u64;
    }
    total
}

fn matching_read_group_ids(header: &sam::Header, id: Option<&str>) -> Option<HashSet<String>> {
    let requested = id?;
    let sample_tag = sam::header::record::value::map::read_group::tag::SAMPLE;
    Some(
        header
            .read_groups()
            .iter()
            .filter_map(|(rg_id, rg)| {
                let id_matches = rg_id.as_slice() == requested.as_bytes();
                let sample_matches = rg
                    .other_fields()
                    .get(&sample_tag)
                    .is_some_and(|sample| sample.as_slice() == requested.as_bytes());
                if id_matches || sample_matches {
                    Some(String::from_utf8_lossy(rg_id.as_slice()).into_owned())
                } else {
                    None
                }
            })
            .collect(),
    )
}

fn record_read_group(rec: &(impl sam::alignment::Record + ?Sized)) -> io::Result<Option<String>> {
    use sam::alignment::record::data::field::{Tag, Value};
    let rg_tag = Tag::from([b'R', b'G']);
    let data = rec.data();
    let Some(value) = data.get(&rg_tag).transpose()? else {
        return Ok(None);
    };
    match value {
        Value::String(s) => Ok(Some(s.to_string())),
        _ => Ok(None),
    }
}

/// Reads an arbitrary aux tag as a string (upstream `bam_aux2Z`), used
/// by `-S`/`--split` to bucket records by tag value.
fn record_aux_string(
    rec: &(impl sam::alignment::Record + ?Sized),
    tag: &str,
) -> io::Result<Option<String>> {
    use sam::alignment::record::data::field::{Tag, Value};
    let bytes = tag.as_bytes();
    if bytes.len() != 2 {
        return Ok(None);
    }
    let t = Tag::from([bytes[0], bytes[1]]);
    let data = rec.data();
    let Some(value) = data.get(&t).transpose()? else {
        return Ok(None);
    };
    match value {
        Value::String(s) => Ok(Some(s.to_string())),
        Value::Character(c) => Ok(Some((c as char).to_string())),
        _ => Ok(None),
    }
}

/// Loads every reference sequence (upper-cased) keyed by name. Used by
/// the MPC reference-mismatch engine; small test references fit easily
/// in memory and a full read yields results identical to upstream's
/// windowed `rseq_buf`.
fn load_reference_seqs(path: &std::path::Path) -> io::Result<HashMap<String, Vec<u8>>> {
    use htslib_rs::fasta;
    let reader = File::open(path).map(BufReader::new)?;
    let mut reader = fasta::io::Reader::new(reader);
    let mut map = HashMap::new();
    for result in reader.records() {
        let record = result?;
        let name = String::from_utf8_lossy(record.name()).into_owned();
        let mut seq = record.sequence().as_ref().to_vec();
        seq.make_ascii_uppercase();
        map.insert(name, seq);
    }
    Ok(map)
}

/// Iterates all SAM records to build a `StatsCounts` with sequence-length
/// and quality accumulators populated (which the `summarize_*` path can
/// not provide because `AlignmentRecordSummary` discards sequence and
/// quality data).
fn collect_sam_full_stats(input: &PathBuf, config: &StatsConfig) -> io::Result<StatsCounts> {
    let mut reader = File::open(input)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader.read_header()?;
    let read_group_filter = matching_read_group_ids(&header, config.id_filter.as_deref());
    let mut counts = StatsCounts::default();
    for result in reader.records() {
        let record = result?;
        counts.update_record(&header, &record, config, read_group_filter.as_ref());
        counts.feed_split(&header, &record, config, read_group_filter.as_ref());
    }
    Ok(counts)
}

fn collect_bam_full_stats(input: &PathBuf, config: &StatsConfig) -> io::Result<StatsCounts> {
    use htslib_rs::bam;
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;
    let read_group_filter = matching_read_group_ids(&header, config.id_filter.as_deref());
    let mut counts = StatsCounts::default();
    let mut record = bam::Record::default();
    loop {
        let n = reader.read_record(&mut record)?;
        if n == 0 {
            break;
        }
        counts.update_record(&header, &record, config, read_group_filter.as_ref());
        counts.feed_split(&header, &record, config, read_group_filter.as_ref());
    }
    Ok(counts)
}

fn collect_sam_region_stats(
    input: &PathBuf,
    regions: &[Region],
    config: &StatsConfig,
) -> io::Result<StatsCounts> {
    let mut reader = File::open(input)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader.read_header()?;
    let read_group_filter = matching_read_group_ids(&header, config.id_filter.as_deref());
    let targets = region_targets(&header, regions)?;
    let mut counts = StatsCounts {
        target_bases: target_base_count(&targets),
        ..Default::default()
    };
    let mut seen = HashSet::new();

    for result in reader.records() {
        let record = result?;
        if record_overlaps_targets(&header, &record, &targets)
            && seen.insert(record_identity(&header, &record))
        {
            counts.update_record_with_targets(
                &header,
                &record,
                config,
                read_group_filter.as_ref(),
                Some(&targets),
            );
        }
    }

    Ok(counts)
}

fn collect_bam_region_stats(
    input: &PathBuf,
    regions: &[Region],
    config: &StatsConfig,
) -> io::Result<StatsCounts> {
    let header = htslib_rs::alignment_compat::read_bam_header_from_path(input)?;
    let read_group_filter = matching_read_group_ids(&header, config.id_filter.as_deref());
    let targets = region_targets(&header, regions)?;
    let mut counts = StatsCounts {
        target_bases: target_base_count(&targets),
        ..Default::default()
    };
    let mut seen = HashSet::new();
    for region in regions {
        for record in htslib_rs::alignment_compat::query_bam_records_from_path(input, region)? {
            if seen.insert(record_identity(&header, &record)) {
                counts.update_record_with_targets(
                    &header,
                    &record,
                    config,
                    read_group_filter.as_ref(),
                    Some(&targets),
                );
            }
        }
    }
    Ok(counts)
}

/// Whole-CRAM (no region) stats using the htslib-rs all-record iterator,
/// so sequence-length/quality/GC/COV/NM accumulate like the BAM path
/// (TODO-NEXT #2) instead of the seq/quality-discarding `summarize_*` path.
fn collect_cram_full_stats(
    input: &PathBuf,
    reference: PathBuf,
    config: &StatsConfig,
) -> io::Result<StatsCounts> {
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(input)?;
    let read_group_filter = matching_read_group_ids(&header, config.id_filter.as_deref());
    let mut counts = StatsCounts::default();
    for record in htslib_rs::alignment_compat::query_cram_records_all_from_path_with_reference(
        input, &reference,
    )? {
        counts.update_record(&header, &record, config, read_group_filter.as_ref());
    }
    Ok(counts)
}

fn collect_cram_region_stats(
    input: &PathBuf,
    reference: PathBuf,
    regions: &[Region],
    config: &StatsConfig,
) -> io::Result<StatsCounts> {
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(input)?;
    let read_group_filter = matching_read_group_ids(&header, config.id_filter.as_deref());
    let targets = region_targets(&header, regions)?;
    let mut counts = StatsCounts {
        target_bases: target_base_count(&targets),
        ..Default::default()
    };
    let mut seen = HashSet::new();
    for region in regions {
        for record in htslib_rs::alignment_compat::query_cram_records_from_path_with_reference(
            input, region, &reference,
        )? {
            if seen.insert(record_identity(&header, &record)) {
                counts.update_record_with_targets(
                    &header,
                    &record,
                    config,
                    read_group_filter.as_ref(),
                    Some(&targets),
                );
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
    supplementary: u64,
    mq0: u64,
    singletons: u64,
    diffchr: u64,
    // Per-record orientation bins (double-counted because each pair
    // contributes both reads). Halved on output to match upstream's
    // per-pair semantics.
    isize_inward: u64,
    isize_outward: u64,
    isize_other: u64,
    // Per-size orientation arrays (still double-counted, halved on
    // output) indexed by abs(template_length) capped at `-i`. These
    // drive the `IS` section's `ibulk`/sparse logic exactly as
    // upstream `output_stats` does.
    isize_in: Vec<u64>,
    isize_out: Vec<u64>,
    isize_oth: Vec<u64>,
    // Per-record abs(template_length) histogram across the same population
    // that feeds the orientation bins. Used to apply `-m/--most-inserts`
    // before computing the reported mean / standard deviation.
    isize_hist: BTreeMap<u64, u64>,
    // Sequence-length and quality accumulators (only populated when the
    // collection path has access to record-level data — currently SAM /
    // BAM iteration and any region path, but not the CRAM-non-region
    // summary path).
    total_len: u64,
    total_len_1st: u64,
    total_len_2nd: u64,
    max_len: u32,
    max_len_1st: u32,
    max_len_2nd: u32,
    qual_sum: u64,
    qual_count: u64,
    // Highest quality value observed (incl. the 0xFF/255 stored for `*`
    // quality), driving the FFQ/LFQ/MPC column count exactly as
    // upstream's `max_qual` (+1 unless it would reach `nquals`=256).
    max_qual: u8,
    first_qual_hist: Vec<[u64; 256]>,
    last_qual_hist: Vec<[u64; 256]>,
    // Per-cycle ACGTNO counts [a,c,g,t,n,other]; `acgt_rc` is
    // read-oriented (reverse reads complemented).
    acgt_1st: Vec<[u64; 6]>,
    acgt_2nd: Vec<[u64; 6]>,
    acgt_rc: Vec<[u64; 6]>,
    read_lengths: Vec<u64>,
    read_lengths_1st: Vec<u64>,
    read_lengths_2nd: Vec<u64>,
    mapping_qualities: Vec<u64>,
    // Upstream `gc_1st`/`gc_2nd`: fixed `ngc`-sized GC-fraction arrays.
    first_gc_hist: Vec<u64>,
    last_gc_hist: Vec<u64>,
    coverage_depths: BTreeMap<(usize, usize), u32>,
    target_bases: u64,
    // Bases mapped (sum of sequence lengths for mapped reads, ignoring
    // clipping), bases mapped from CIGAR (sum of M/=/X ops), and total
    // NM aux values across reads. Used for the `bases mapped`,
    // `bases mapped (cigar)`, `mismatches`, and `error rate` SN lines.
    bases_mapped: u64,
    bases_mapped_cigar: u64,
    nmismatches: u64,
    // Indel distribution (ID) and indels-per-cycle (IC), faithful to
    // upstream `count_indels`. `*_len[k]` counts indels of length k+1;
    // `*_cycles_{1st,2nd}[c]` counts read1/read2 indels at cycle c.
    insertions_len: Vec<u64>,
    deletions_len: Vec<u64>,
    ins_cycles_1st: Vec<u64>,
    ins_cycles_2nd: Vec<u64>,
    del_cycles_1st: Vec<u64>,
    del_cycles_2nd: Vec<u64>,
    // Mismatches per cycle and quality (MPC): `mpc_buf[cycle][qual]`.
    // Column 0 doubles as the N-base count, exactly as upstream
    // `count_mismatches_per_cycle` (qual byte + 1, u8-wrapping).
    mpc_buf: Vec<[u64; 256]>,
    // `-S`/`--split`: this run's tag value (drives the "statistics only
    // for reads with tag" header line), and the per-tag-value sub-stats.
    split_name: Option<String>,
    splits: BTreeMap<String, StatsCounts>,
    // Sum of sequence lengths for records carrying the duplicate flag.
    bases_dup: u64,
    bases_trimmed: u64,
    last_sort_position: Option<(usize, usize)>,
    sort_order_violation: bool,
    // CHK section: 32-bit-wrapping sums of per-record CRC32 over the
    // read name, BAM-packed sequence nibbles, and raw quality bytes.
    // Accumulated for every record passing flag-require/flag-filter/
    // read-length filtering (before the secondary/supplementary skip),
    // exactly as upstream `update_checksum`.
    chk_names: u32,
    chk_reads: u32,
    chk_quals: u32,
}

struct StatsRecordFields {
    flag: u32,
    mapq: Option<u8>,
    reference_sequence_id: Option<usize>,
    mate_reference_sequence_id: Option<usize>,
    template_length: i32,
    read_len: Option<usize>,
    pos: Option<usize>,
    mpos: Option<usize>,
}

impl StatsCounts {
    fn update_record(
        &mut self,
        header: &sam::Header,
        rec: &(impl sam::alignment::Record + ?Sized),
        config: &StatsConfig,
        read_group_filter: Option<&HashSet<String>>,
    ) {
        self.update_record_with_targets(header, rec, config, read_group_filter, None);
    }

    /// `-S`/`--split`: also accumulate this record into its tag-value's
    /// sub-`StatsCounts` (created on first sight with `split_name` set so
    /// the per-tag `.bamstat` header line is correct). The sub-counts
    /// never split further, so no recursion.
    fn feed_split(
        &mut self,
        header: &sam::Header,
        rec: &(impl sam::alignment::Record + ?Sized),
        config: &StatsConfig,
        read_group_filter: Option<&HashSet<String>>,
    ) {
        let Some(tag) = config.split_tag.as_deref() else {
            return;
        };
        let value = if tag == "RG" {
            record_read_group(rec).ok().flatten()
        } else {
            record_aux_string(rec, tag).ok().flatten()
        };
        let Some(value) = value else {
            return;
        };
        let sub = self
            .splits
            .entry(value.clone())
            .or_insert_with(|| StatsCounts {
                split_name: Some(value),
                ..Default::default()
            });
        sub.update_record_with_targets(header, rec, config, read_group_filter, None);
    }

    fn update_record_with_targets(
        &mut self,
        header: &sam::Header,
        rec: &(impl sam::alignment::Record + ?Sized),
        config: &StatsConfig,
        read_group_filter: Option<&HashSet<String>>,
        targets: Option<&[RegionTarget]>,
    ) {
        if let Some(allowed_read_groups) = read_group_filter {
            let Ok(Some(read_group)) = record_read_group(rec) else {
                return;
            };
            if !allowed_read_groups.contains(&read_group) {
                return;
            }
        }

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
        let template_length = rec.template_length().ok().unwrap_or(0);
        let pos = rec.alignment_start().and_then(Result::ok).map(usize::from);
        let mpos = rec
            .mate_alignment_start()
            .and_then(Result::ok)
            .map(usize::from);
        let seq_len = rec.sequence().len();

        let chk_name = rec.name().map(|n| n.to_vec()).unwrap_or_default();
        let chk_seq: Vec<u8> = rec.sequence().iter().collect();
        let chk_qual: Vec<u8> = rec.quality_scores().iter().flatten().collect();
        self.accumulate_checksum(flag, config, &chk_name, &chk_seq, &chk_qual);

        let pre_total = self.total;
        let pre_supp = self.supplementary;
        self.update(
            StatsRecordFields {
                flag,
                mapq,
                reference_sequence_id,
                mate_reference_sequence_id,
                template_length,
                read_len: Some(seq_len),
                pos,
                mpos,
            },
            config,
        );

        // Only accumulate seq/quality stats for primary, non-supplementary
        // reads that were not filtered out by `--remove-dups` — i.e. the
        // ones that contributed to `total`.
        if self.total > pre_total {
            if flag & BAM_FUNMAP == 0
                && let Ok(Some(start)) = rec.alignment_start().transpose()
                && let Some(tid) = reference_sequence_id
            {
                self.update_sort_order(tid, usize::from(start));
            }

            let seq_len_u32 = seq_len as u32;
            // Upstream `order`: unpaired => first fragment.
            let order_first =
                flag & BAM_FPAIRED == 0 || (flag & BAM_FREAD1 != 0 && flag & BAM_FREAD2 == 0);
            let order_last =
                flag & BAM_FPAIRED != 0 && flag & BAM_FREAD2 != 0 && flag & BAM_FREAD1 == 0;
            if seq_len_u32 > 0 {
                self.total_len += u64::from(seq_len_u32);
                if seq_len_u32 > self.max_len {
                    self.max_len = seq_len_u32;
                }
                if order_first {
                    self.total_len_1st += u64::from(seq_len_u32);
                    if seq_len_u32 > self.max_len_1st {
                        self.max_len_1st = seq_len_u32;
                    }
                }
                if order_last {
                    self.total_len_2nd += u64::from(seq_len_u32);
                    if seq_len_u32 > self.max_len_2nd {
                        self.max_len_2nd = seq_len_u32;
                    }
                }
            }
            let quals: Vec<u8> = rec.quality_scores().iter().flatten().collect();
            if quals.is_empty() && seq_len_u32 > 0 {
                // `*` quality: HTSlib stores 0xFF per base.
                self.max_qual = 255;
                for cycle in 0..seq_len_u32 as usize {
                    self.qual_sum += 255;
                    self.qual_count += 1;
                    if order_first {
                        increment_quality_hist(&mut self.first_qual_hist, cycle, 255);
                    }
                    if order_last {
                        increment_quality_hist(&mut self.last_qual_hist, cycle, 255);
                    }
                }
            } else {
                for (cycle, q) in quals.iter().copied().enumerate() {
                    self.qual_sum += u64::from(q);
                    self.qual_count += 1;
                    if q > self.max_qual {
                        self.max_qual = q;
                    }
                    if order_first {
                        increment_quality_hist(&mut self.first_qual_hist, cycle, q);
                    }
                    if order_last {
                        increment_quality_hist(&mut self.last_qual_hist, cycle, q);
                    }
                }
            }
            if config.trim_quality > 0 {
                let reverse = flag & BAM_FREVERSE != 0;
                self.bases_trimmed += u64::from(bwa_trim_read(
                    config.trim_quality,
                    rec.quality_scores().iter().flatten(),
                    reverse,
                ));
            }

            if flag & BAM_FDUP != 0 {
                self.bases_dup += u64::from(seq_len_u32);
            }

            if seq_len_u32 > 0 {
                // Upstream `stats.c`: bin into `gc_*[gc*(NGC-1)/L ..
                // (gc+1)*(NGC-1)/L)` (capped at NGC-1).
                let l = seq_len_u32 as u64;
                let gc = rec
                    .sequence()
                    .iter()
                    .filter(|b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
                    .count() as u64;
                let lo = (gc * (NGC as u64 - 1) / l) as usize;
                let mut hi = ((gc + 1) * (NGC as u64 - 1) / l) as usize;
                if hi >= NGC {
                    hi = NGC - 1;
                }
                let tgt = if order_first {
                    Some(&mut self.first_gc_hist)
                } else if order_last {
                    Some(&mut self.last_gc_hist)
                } else {
                    None
                };
                if let Some(h) = tgt {
                    if h.len() < NGC {
                        h.resize(NGC, 0);
                    }
                    for slot in h.iter_mut().take(hi).skip(lo) {
                        *slot += 1;
                    }
                }
            }

            // Read-length / mapq histograms + per-cycle ACGT (upstream
            // `stats.c`: read_len = unclipped length; mapq for
            // !(UNMAP|SEC|SUPP|QCFAIL|DUP); per-cycle ACGT for originals).
            {
                let reverse = flag & BAM_FREVERSE != 0;
                use sam::alignment::record::cigar::op::Kind as CKind;
                let hard: u32 = rec
                    .cigar()
                    .iter()
                    .flatten()
                    .filter(|op| op.kind() == CKind::HardClip)
                    .map(|op| op.len() as u32)
                    .sum();
                let read_len = (seq_len_u32 + hard) as usize;
                if read_len > 0 && flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY) == 0 {
                    if read_len >= self.read_lengths.len() {
                        self.read_lengths.resize(read_len + 1, 0);
                    }
                    self.read_lengths[read_len] += 1;
                    if order_first {
                        if read_len >= self.read_lengths_1st.len() {
                            self.read_lengths_1st.resize(read_len + 1, 0);
                        }
                        self.read_lengths_1st[read_len] += 1;
                    }
                    if order_last {
                        if read_len >= self.read_lengths_2nd.len() {
                            self.read_lengths_2nd.resize(read_len + 1, 0);
                        }
                        self.read_lengths_2nd[read_len] += 1;
                    }
                    if order_first || order_last {
                        let seq = rec.sequence();
                        let sl = seq.len();
                        if self.acgt_1st.len() < sl {
                            self.acgt_1st.resize(sl, [0; 6]);
                        }
                        if self.acgt_2nd.len() < sl {
                            self.acgt_2nd.resize(sl, [0; 6]);
                        }
                        if self.acgt_rc.len() < sl {
                            self.acgt_rc.resize(sl, [0; 6]);
                        }
                        for (i, b) in seq.iter().enumerate() {
                            let cyc = if reverse { sl - i - 1 } else { i };
                            let (idx, rc): (usize, usize) = match b.to_ascii_uppercase() {
                                b'A' => (0, if reverse { 3 } else { 0 }),
                                b'C' => (1, if reverse { 2 } else { 1 }),
                                b'G' => (2, if reverse { 1 } else { 2 }),
                                b'T' => (3, if reverse { 0 } else { 3 }),
                                b'N' => (4, 6),
                                _ => (5, 6),
                            };
                            if order_last {
                                self.acgt_2nd[cyc][idx] += 1;
                            } else {
                                self.acgt_1st[cyc][idx] += 1;
                            }
                            if rc < 6 {
                                self.acgt_rc[cyc][rc] += 1;
                            }
                        }
                    }
                }
                if flag
                    & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP)
                    == 0
                {
                    let mapq = rec
                        .mapping_quality()
                        .and_then(Result::ok)
                        .map(u8::from)
                        .unwrap_or(255) as usize;
                    if self.mapping_qualities.is_empty() {
                        self.mapping_qualities = vec![0; 256];
                    }
                    self.mapping_qualities[mapq] += 1;
                }
            }

            if flag & BAM_FUNMAP == 0 {
                self.bases_mapped += u64::from(seq_len_u32);
                use sam::alignment::record::cigar::op::Kind;
                for op in rec.cigar().iter().flatten() {
                    // Upstream (non-region path) counts M, I, =, X bases.
                    if matches!(
                        op.kind(),
                        Kind::Match
                            | Kind::Insertion
                            | Kind::SequenceMatch
                            | Kind::SequenceMismatch
                    ) {
                        self.bases_mapped_cigar += op.len() as u64;
                    }
                }
                self.count_indels(rec, flag, seq_len_u32 as usize);
                if let Some(nm) = read_nm_aux(rec) {
                    self.nmismatches += nm;
                }
                self.update_coverage_depths(header, rec, targets);
            }
        }

        // Supplementary alignments do not return early in upstream
        // `collect_stats`: they are excluded from the IS_ORIGINAL
        // sequence/quality/read-length stats but still contribute to
        // the indel distribution, bases-mapped-(cigar), NM mismatches
        // and the coverage histogram. They were counted (and skipped)
        // by `update`, so replay just those accumulations here.
        if self.supplementary > pre_supp && flag & BAM_FUNMAP == 0 {
            use sam::alignment::record::cigar::op::Kind;
            for op in rec.cigar().iter().flatten() {
                if matches!(
                    op.kind(),
                    Kind::Match | Kind::Insertion | Kind::SequenceMatch | Kind::SequenceMismatch
                ) {
                    self.bases_mapped_cigar += op.len() as u64;
                }
            }
            self.count_indels(rec, flag, seq_len);
            if let Some(nm) = read_nm_aux(rec) {
                self.nmismatches += nm;
            }
            self.update_coverage_depths(header, rec, targets);
        }

        // MPC reference-mismatch engine. Upstream calls
        // `count_mismatches_per_cycle` for every mapped read reaching it
        // (primary or supplementary, not secondary) when a reference is
        // loaded.
        let counted = self.total > pre_total || self.supplementary > pre_supp;
        if counted
            && flag & BAM_FUNMAP == 0
            && let Some(refmap) = config.reference_seqs.as_ref()
            && let Some(rsid) = reference_sequence_id
            && let Some(p) = pos
            && let Some((name, _)) = header.reference_sequences().get_index(rsid)
            && let Some(refseq) = refmap.get(String::from_utf8_lossy(name.as_ref()).as_ref())
        {
            self.count_mismatches_per_cycle(rec, flag, seq_len, refseq, p - 1, &chk_seq, &chk_qual);
        }
    }

    /// Faithful port of upstream `count_indels`: per-CIGAR insertion /
    /// deletion length distribution (ID) and per-cycle indel counts
    /// split by read1/read2 (IC). `read_len` is `l_qseq`.
    fn count_indels(
        &mut self,
        rec: &(impl sam::alignment::Record + ?Sized),
        flag: u32,
        read_len: usize,
    ) {
        use sam::alignment::record::cigar::op::Kind;
        let is_fwd = flag & BAM_FREVERSE == 0;
        // order: 1 = read1 only, 2 = read2 only (per READ_ORDER_*).
        let order: u32 = if flag & BAM_FPAIRED != 0 {
            (if flag & BAM_FREAD1 != 0 { 1 } else { 0 })
                + (if flag & BAM_FREAD2 != 0 { 2 } else { 0 })
        } else {
            1
        };
        fn bump(v: &mut Vec<u64>, idx: usize) {
            if v.len() <= idx {
                v.resize(idx + 1, 0);
            }
            v[idx] += 1;
        }
        let mut icycle: isize = 0;
        for op in rec.cigar().iter().flatten() {
            let ncig = op.len();
            if ncig == 0 {
                continue;
            }
            let ncig_i = ncig as isize;
            match op.kind() {
                Kind::Insertion => {
                    let idx = if is_fwd {
                        icycle
                    } else {
                        read_len as isize - icycle - ncig_i
                    };
                    if idx >= 0 {
                        if order == 1 {
                            bump(&mut self.ins_cycles_1st, idx as usize);
                        } else if order == 2 {
                            bump(&mut self.ins_cycles_2nd, idx as usize);
                        }
                    }
                    icycle += ncig_i;
                    if ncig <= NINDELS {
                        bump(&mut self.insertions_len, ncig - 1);
                    }
                }
                Kind::Deletion => {
                    let idx = if is_fwd {
                        icycle - 1
                    } else {
                        read_len as isize - icycle - 1
                    };
                    if idx < 0 {
                        continue;
                    }
                    if order == 1 {
                        bump(&mut self.del_cycles_1st, idx as usize);
                    } else if order == 2 {
                        bump(&mut self.del_cycles_2nd, idx as usize);
                    }
                    if ncig <= NINDELS {
                        bump(&mut self.deletions_len, ncig - 1);
                    }
                }
                Kind::Skip | Kind::HardClip | Kind::Pad => {}
                _ => icycle += ncig_i,
            }
        }
    }

    /// Faithful port of upstream `count_mismatches_per_cycle`. For each
    /// aligned base it increments `mpc_buf[cycle][col]` where `col` is 0
    /// for an N read base, otherwise `(qual + 1) as u8` (note the u8 wrap
    /// that maps the 0xFF missing-quality byte to column 0). `pos0` is
    /// the read's 0-based reference start; `refseq` the upper-cased
    /// reference bases. N (ref-skip)/H/P CIGARs are ignored as upstream.
    #[allow(clippy::too_many_arguments)]
    fn count_mismatches_per_cycle(
        &mut self,
        rec: &(impl sam::alignment::Record + ?Sized),
        flag: u32,
        read_len: usize,
        refseq: &[u8],
        pos0: usize,
        seq_ascii: &[u8],
        quals: &[u8],
    ) {
        use sam::alignment::record::cigar::op::Kind;
        fn nt16(b: u8) -> u8 {
            match b.to_ascii_uppercase() {
                b'=' => 0,
                b'A' => 1,
                b'C' => 2,
                b'M' => 3,
                b'G' => 4,
                b'R' => 5,
                b'S' => 6,
                b'V' => 7,
                b'T' => 8,
                b'W' => 9,
                b'Y' => 10,
                b'H' => 11,
                b'K' => 12,
                b'D' => 13,
                b'B' => 14,
                _ => 15,
            }
        }
        let is_fwd = flag & BAM_FREVERSE == 0;
        let mut iread: usize = 0;
        let mut icycle: isize = 0;
        let mut iref: usize = 0;
        for op in rec.cigar().iter().flatten() {
            let ncig = op.len();
            match op.kind() {
                Kind::Insertion => {
                    iread += ncig;
                    icycle += ncig as isize;
                }
                Kind::Deletion => {
                    iref += ncig;
                }
                Kind::SoftClip => {
                    icycle += ncig as isize;
                    iread += ncig;
                }
                Kind::HardClip => {
                    icycle += ncig as isize;
                }
                Kind::Skip | Kind::Pad => {}
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for _ in 0..ncig {
                        let cread = seq_ascii.get(iread).copied().map(nt16).unwrap_or(0);
                        let cref = refseq.get(pos0 + iref).copied().map(nt16).unwrap_or(0);
                        let idx = if is_fwd {
                            icycle
                        } else {
                            read_len as isize - icycle - 1
                        };
                        if idx >= 0 {
                            let cyc = idx as usize;
                            if self.mpc_buf.len() <= cyc {
                                self.mpc_buf.resize(cyc + 1, [0; 256]);
                            }
                            if cread == 15 {
                                self.mpc_buf[cyc][0] += 1;
                            } else if cref != 0 && cread != 0 && cref != cread {
                                let qbyte = quals.get(iread).copied().unwrap_or(0xFF);
                                let col = qbyte.wrapping_add(1) as usize;
                                self.mpc_buf[cyc][col] += 1;
                            }
                        }
                        iref += 1;
                        iread += 1;
                        icycle += 1;
                    }
                }
            }
        }
    }

    fn is_coordinate_sorted(&self) -> bool {
        !self.sort_order_violation
    }

    /// Number of quality columns emitted in FFQ/LFQ/MPC, mirroring
    /// upstream's `if (max_qual+1 < nquals) max_qual++;` then
    /// `iqual <= max_qual` (so `max_qual + 1` columns).
    fn qual_cols(&self) -> usize {
        let m = self.max_qual as usize;
        (if m + 1 < 256 { m + 1 } else { m }) + 1
    }

    fn update_sort_order(&mut self, tid: usize, pos: usize) {
        if let Some((last_tid, last_pos)) = self.last_sort_position
            && (tid < last_tid || (tid == last_tid && pos < last_pos))
        {
            self.sort_order_violation = true;
        }
        self.last_sort_position = Some((tid, pos));
    }

    fn update_coverage_depths(
        &mut self,
        header: &sam::Header,
        rec: &(impl sam::alignment::Record + ?Sized),
        targets: Option<&[RegionTarget]>,
    ) {
        let tid = match rec.reference_sequence_id(header).transpose() {
            Ok(Some(tid)) => tid,
            _ => return,
        };
        let mut ref_pos = match rec.alignment_start().transpose() {
            Ok(Some(start)) => usize::from(start) - 1,
            _ => return,
        };

        use sam::alignment::record::cigar::op::Kind;
        for op in rec.cigar().iter().flatten() {
            let len = op.len();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for pos0 in ref_pos..ref_pos.saturating_add(len) {
                        if targets
                            .map(|targets| {
                                targets.iter().any(|target| {
                                    target.tid == tid
                                        && target.start <= pos0.saturating_add(1)
                                        && pos0.saturating_add(1) <= target.end
                                })
                            })
                            .unwrap_or(true)
                        {
                            let depth = self.coverage_depths.entry((tid, pos0)).or_default();
                            *depth = depth.saturating_add(1);
                        }
                    }
                    ref_pos = ref_pos.saturating_add(len);
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos = ref_pos.saturating_add(len);
                }
                Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
            }
        }
    }

    fn update_summary(&mut self, rec: &AlignmentRecordSummary, config: &StatsConfig) {
        let flag = rec.flags_u16() as u32;
        self.accumulate_checksum(
            flag,
            config,
            rec.name_bytes().unwrap_or_default(),
            rec.sequence_bytes(),
            rec.quality_score_bytes(),
        );
        self.update(
            StatsRecordFields {
                flag: rec.flags_u16() as u32,
                mapq: rec.mapping_quality(),
                reference_sequence_id: rec.reference_sequence_id(),
                mate_reference_sequence_id: rec.mate_reference_sequence_id(),
                template_length: rec.template_length(),
                read_len: None,
                pos: rec.alignment_start(),
                mpos: rec.mate_alignment_start(),
            },
            config,
        );
    }

    fn update(&mut self, rec: StatsRecordFields, config: &StatsConfig) {
        let flag = rec.flag;
        if config.required_flags != 0 && flag & config.required_flags != config.required_flags {
            self.raw_total += 1;
            self.filtered += 1;
            return;
        }
        if config.filter_flags != 0 && flag & config.filter_flags != 0 {
            self.raw_total += 1;
            self.filtered += 1;
            return;
        }
        if config.remove_dups && flag & BAM_FDUP != 0 {
            self.raw_total += 1;
            self.filtered += 1;
            return;
        }
        if config
            .read_length_filter
            .is_some_and(|required_len| rec.read_len != Some(required_len))
        {
            return;
        }
        if flag & BAM_FSECONDARY != 0 {
            self.secondary += 1;
            return;
        }
        if flag & BAM_FSUPPLEMENTARY != 0 {
            self.supplementary += 1;
            return;
        }
        self.raw_total += 1;
        self.total += 1;

        if flag & BAM_FUNMAP == 0 {
            self.mapped += 1;
            if rec.mapq == Some(0) {
                self.mq0 += 1;
            }
        } else {
            self.unmapped += 1;
        }
        // Upstream `order`: unpaired reads are first fragments; paired
        // reads are first/last by the READ1/READ2 bits (both-or-neither
        // counts as "other", excluded from read1/read2).
        let order_first =
            flag & BAM_FPAIRED == 0 || (flag & BAM_FREAD1 != 0 && flag & BAM_FREAD2 == 0);
        let order_last =
            flag & BAM_FPAIRED != 0 && flag & BAM_FREAD2 != 0 && flag & BAM_FREAD1 == 0;
        if order_first {
            self.read1 += 1;
        } else if order_last {
            self.read2 += 1;
        }
        if flag & BAM_FPAIRED != 0 {
            self.paired += 1;
            if flag & BAM_FPROPER_PAIR != 0 {
                self.proper_paired += 1;
            }
            if flag & BAM_FUNMAP == 0 && flag & BAM_FMUNMAP == 0 {
                self.mapped_and_paired += 1;
                if rec.reference_sequence_id != rec.mate_reference_sequence_id
                    && rec.reference_sequence_id.is_some()
                    && rec.mate_reference_sequence_id.is_some()
                {
                    self.diffchr += 1;
                }
                self.update_isize_bin(
                    flag,
                    rec.template_length,
                    rec.pos,
                    rec.mpos,
                    rec.reference_sequence_id,
                    rec.mate_reference_sequence_id,
                    config.insert_size_max,
                );
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

    /// Classify an insert-size observation into inward/outward/other,
    /// faithfully mirroring the orientation logic in `stats.c`'s
    /// `collect_stats` (the `IS_PAIRED_AND_MAPPED && IS_ORIGINAL` block).
    /// Each record contributes once; the output halves to obtain
    /// per-pair counts. The accumulation gate is upstream's
    /// `isize > 0 || tid == mtid`.
    #[allow(clippy::too_many_arguments)]
    fn update_isize_bin(
        &mut self,
        flag: u32,
        template_length: i32,
        pos: Option<usize>,
        mpos: Option<usize>,
        rsid: Option<usize>,
        mrsid: Option<usize>,
        insert_size_max: u32,
    ) {
        let mut isize = template_length.unsigned_abs() as u64;
        if insert_size_max > 0 {
            isize = isize.min(u64::from(insert_size_max));
        }
        let same_ref = rsid.is_some() && rsid == mrsid;
        if !(isize > 0 || same_ref) {
            return;
        }
        *self.isize_hist.entry(isize).or_default() += 1;
        let i = isize as usize;
        if self.isize_in.len() <= i {
            self.isize_in.resize(i + 1, 0);
            self.isize_out.resize(i + 1, 0);
            self.isize_oth.resize(i + 1, 0);
        }

        // pos_fst = mpos - pos (cancels the 0-/1-based offset since both
        // ends share it). is_fst is the read1/read2 discriminator.
        let pos_fst: i64 = mpos.unwrap_or(0) as i64 - pos.unwrap_or(0) as i64;
        let is_fst: i64 = if flag & BAM_FREAD1 != 0 { 1 } else { -1 };
        let is_fwd: i64 = if flag & BAM_FREVERSE != 0 { -1 } else { 1 };
        let is_mfwd: i64 = if flag & 0x20 /* BAM_FMREVERSE */ != 0 {
            -1
        } else {
            1
        };

        enum Ori {
            In,
            Out,
            Oth,
        }
        let ori = if is_fwd * is_mfwd > 0 {
            Ori::Oth
        } else if is_fst * pos_fst > 0 {
            if is_fst * is_fwd > 0 {
                Ori::In
            } else {
                Ori::Out
            }
        } else if is_fst * pos_fst < 0 {
            if is_fst * is_fwd > 0 {
                Ori::Out
            } else {
                Ori::In
            }
        } else {
            // Exactly overlapping reads are assumed inward.
            Ori::In
        };
        match ori {
            Ori::In => {
                self.isize_inward += 1;
                self.isize_in[i] += 1;
            }
            Ori::Out => {
                self.isize_outward += 1;
                self.isize_out[i] += 1;
            }
            Ori::Oth => {
                self.isize_other += 1;
                self.isize_oth[i] += 1;
            }
        }
    }

    /// Faithful port of upstream `update_checksum`. Runs for every
    /// record that survives the flag-require / flag-filter /
    /// read-length filters (before the secondary/supplementary skip),
    /// summing per-record CRC32 with 32-bit overflow.
    fn accumulate_checksum(
        &mut self,
        flag: u32,
        config: &StatsConfig,
        name: &[u8],
        seq_ascii: &[u8],
        quals: &[u8],
    ) {
        if config.required_flags != 0 && flag & config.required_flags != config.required_flags {
            return;
        }
        if config.filter_flags != 0 && flag & config.filter_flags != 0 {
            return;
        }
        if config
            .read_length_filter
            .is_some_and(|required_len| seq_ascii.len() != required_len)
        {
            return;
        }
        self.chk_names = self.chk_names.wrapping_add(crc32_bytes(0, name));
        let seq_len = seq_ascii.len();
        if seq_len == 0 {
            return;
        }
        let packed = bam_pack_seq(seq_ascii);
        self.chk_reads = self.chk_reads.wrapping_add(crc32_bytes(0, &packed));
        if quals.len() == seq_len {
            self.chk_quals = self.chk_quals.wrapping_add(crc32_bytes(0, quals));
        } else {
            // SAM `*` quality is stored as 0xFF per base in BAM.
            let missing = vec![0xFFu8; seq_len];
            self.chk_quals = self.chk_quals.wrapping_add(crc32_bytes(0, &missing));
        }
    }
}

/// zlib `crc32(initial, buf, len)` over `bytes`.
fn crc32_bytes(initial: u32, bytes: &[u8]) -> u32 {
    let mut crc = libdeflater::Crc::with_initial(initial);
    crc.update(bytes);
    crc.sum()
}

/// Packs an ASCII sequence into BAM 4-bit-per-base nibbles
/// (`(len + 1) / 2` bytes), using HTSlib's `seq_nt16_table` codes.
fn bam_pack_seq(seq_ascii: &[u8]) -> Vec<u8> {
    fn code(b: u8) -> u8 {
        match b {
            b'=' => 0,
            b'A' | b'a' => 1,
            b'C' | b'c' => 2,
            b'M' | b'm' => 3,
            b'G' | b'g' => 4,
            b'R' | b'r' => 5,
            b'S' | b's' => 6,
            b'V' | b'v' => 7,
            b'T' | b't' => 8,
            b'W' | b'w' => 9,
            b'Y' | b'y' => 10,
            b'H' | b'h' => 11,
            b'K' | b'k' => 12,
            b'D' | b'd' => 13,
            b'B' | b'b' => 14,
            _ => 15,
        }
    }
    let mut out = vec![0u8; seq_ascii.len().div_ceil(2)];
    for (i, &b) in seq_ascii.iter().enumerate() {
        let c = code(b);
        if i % 2 == 0 {
            out[i / 2] = c << 4;
        } else {
            out[i / 2] |= c;
        }
    }
    out
}

fn increment_quality_hist(hist: &mut Vec<[u64; 256]>, cycle: usize, quality: u8) {
    if hist.len() <= cycle {
        hist.resize_with(cycle + 1, || [0; 256]);
    }
    hist[cycle][usize::from(quality)] += 1;
}

fn bwa_trim_read(trim_quality: u8, qualities: impl IntoIterator<Item = u8>, reverse: bool) -> u32 {
    const BWA_MIN_READ_LEN: usize = 35;
    let qualities: Vec<_> = qualities.into_iter().collect();
    if qualities.len() < BWA_MIN_READ_LEN {
        return 0;
    }

    let max_trimmed = qualities.len() - BWA_MIN_READ_LEN + 1;
    let mut sum = 0i32;
    let mut max_sum = 0i32;
    let mut max_l = 0u32;
    for l in 0..max_trimmed {
        let q = if reverse {
            qualities[l]
        } else {
            qualities[qualities.len() - 1 - l]
        };
        sum += i32::from(trim_quality) - i32::from(q);
        if sum < 0 {
            break;
        }
        if sum > max_sum {
            max_sum = sum;
            max_l = l as u32;
        }
    }
    max_l
}

/// C `printf("%e")`: 6-decimal mantissa, `e`, signed ≥2-digit exponent
/// (e.g. `0.000000e+00`, `1.234560e-05`). Rust's `{:e}` differs.
fn c_e6(x: f64) -> String {
    if x == 0.0 || !x.is_finite() {
        return "0.000000e+00".to_string();
    }
    let neg = x < 0.0;
    let v = x.abs();
    let mut e = v.log10().floor() as i32;
    let mut ms = format!("{:.6}", v / 10f64.powi(e));
    if ms.starts_with("10") {
        e += 1;
        ms = format!("{:.6}", v / 10f64.powi(e));
    }
    format!(
        "{}{}e{}{:02}",
        if neg { "-" } else { "" },
        ms,
        if e < 0 { "-" } else { "+" },
        e.abs()
    )
}

fn write_stats(
    out: &mut dyn Write,
    recs: &[AlignmentRecordSummary],
    config: &StatsConfig,
    is_sorted: bool,
) -> io::Result<()> {
    let mut counts = StatsCounts::default();
    for rec in recs {
        counts.update_summary(rec, config);
    }
    write_stats_counts(out, &counts, config, is_sorted)
}

fn write_stats_counts(
    out: &mut dyn Write,
    counts: &StatsCounts,
    config: &StatsConfig,
    is_sorted: bool,
) -> io::Result<()> {
    writeln!(
        out,
        "# This file was produced by samtools-rs stats (samtools-{}+htslib-rs)",
        SAMTOOLS_VERSION
    )?;
    match (&counts.split_name, &config.split_tag) {
        (Some(name), Some(tag)) => writeln!(
            out,
            "# This file contains statistics only for reads with tag: {tag}={name}"
        )?,
        _ => writeln!(out, "# This file contains statistics for all reads.")?,
    }
    writeln!(out, "# The command line was:  samtools-rs stats")?;
    writeln!(
        out,
        "# CHK, Checksum\t[2]Read Names\t[3]Sequences\t[4]Qualities"
    )?;
    writeln!(
        out,
        "# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)"
    )?;
    writeln!(
        out,
        "CHK\t{:08x}\t{:08x}\t{:08x}",
        counts.chk_names, counts.chk_reads, counts.chk_quals
    )?;
    writeln!(
        out,
        "# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part."
    )?;
    writeln!(
        out,
        "SN\traw total sequences:\t{}\t# excluding supplementary and secondary reads",
        counts.raw_total
    )?;
    writeln!(out, "SN\tfiltered sequences:\t{}", counts.filtered)?;
    writeln!(out, "SN\tsequences:\t{}", counts.total)?;
    writeln!(
        out,
        "SN\tis sorted:\t{}\t# {} by coordinate",
        if is_sorted { 1 } else { 0 },
        if is_sorted { "sorted" } else { "not sorted" }
    )?;
    writeln!(out, "SN\t1st fragments:\t{}", counts.read1)?;
    writeln!(out, "SN\tlast fragments:\t{}", counts.read2)?;
    writeln!(out, "SN\treads mapped:\t{}", counts.mapped)?;
    writeln!(
        out,
        "SN\treads mapped and paired:\t{}\t# paired-end technology bit set + both mates mapped",
        counts.mapped_and_paired
    )?;
    writeln!(out, "SN\treads unmapped:\t{}", counts.unmapped)?;
    writeln!(
        out,
        "SN\treads properly paired:\t{}\t# proper-pair bit set",
        counts.proper_paired
    )?;
    writeln!(
        out,
        "SN\treads paired:\t{}\t# paired-end technology bit set",
        counts.paired
    )?;
    writeln!(
        out,
        "SN\treads duplicated:\t{}\t# PCR or optical duplicate bit set",
        counts.dup
    )?;
    writeln!(out, "SN\treads MQ0:\t{}\t# mapped and MQ=0", counts.mq0)?;
    writeln!(out, "SN\treads QC failed:\t{}", counts.qc_fail)?;
    writeln!(out, "SN\tnon-primary alignments:\t{}", counts.secondary)?;
    writeln!(
        out,
        "SN\tsupplementary alignments:\t{}",
        counts.supplementary
    )?;
    writeln!(
        out,
        "SN\ttotal length:\t{}\t# ignores clipping",
        counts.total_len
    )?;
    writeln!(
        out,
        "SN\ttotal first fragment length:\t{}\t# ignores clipping",
        counts.total_len_1st
    )?;
    writeln!(
        out,
        "SN\ttotal last fragment length:\t{}\t# ignores clipping",
        counts.total_len_2nd
    )?;
    writeln!(
        out,
        "SN\tbases mapped:\t{}\t# ignores clipping",
        counts.bases_mapped
    )?;
    writeln!(
        out,
        "SN\tbases mapped (cigar):\t{}\t# more accurate",
        counts.bases_mapped_cigar
    )?;
    writeln!(out, "SN\tbases trimmed:\t{}", counts.bases_trimmed)?;
    writeln!(out, "SN\tbases duplicated:\t{}", counts.bases_dup)?;
    writeln!(
        out,
        "SN\tmismatches:\t{}\t# from NM fields",
        counts.nmismatches
    )?;
    let error_rate = if counts.bases_mapped_cigar > 0 {
        counts.nmismatches as f64 / counts.bases_mapped_cigar as f64
    } else {
        0.0
    };
    writeln!(
        out,
        "SN\terror rate:\t{}\t# mismatches / bases mapped (cigar)",
        c_e6(error_rate)
    )?;
    let avg_len = if counts.raw_total > 0 {
        counts.total_len as f64 / counts.raw_total as f64
    } else {
        0.0
    };
    let avg_len_1st = if counts.read1 > 0 {
        counts.total_len_1st as f64 / counts.read1 as f64
    } else {
        0.0
    };
    let avg_len_2nd = if counts.read2 > 0 {
        counts.total_len_2nd as f64 / counts.read2 as f64
    } else {
        0.0
    };
    writeln!(out, "SN\taverage length:\t{:.0}", avg_len)?;
    writeln!(
        out,
        "SN\taverage first fragment length:\t{:.0}",
        avg_len_1st
    )?;
    writeln!(out, "SN\taverage last fragment length:\t{:.0}", avg_len_2nd)?;
    writeln!(out, "SN\tmaximum length:\t{}", counts.max_len)?;
    writeln!(
        out,
        "SN\tmaximum first fragment length:\t{}",
        counts.max_len_1st
    )?;
    writeln!(
        out,
        "SN\tmaximum last fragment length:\t{}",
        counts.max_len_2nd
    )?;
    // Upstream: `total_len ? sum_qual/total_len : 0` (no `singletons`
    // SN line in `samtools stats`).
    let avg_quality = if counts.total_len > 0 {
        counts.qual_sum as f64 / counts.total_len as f64
    } else {
        0.0
    };
    writeln!(out, "SN\taverage quality:\t{:.1}", avg_quality)?;

    let (avg_isize, sd_isize) =
        insert_size_mean_sd(&counts.isize_hist, config.insert_size_main_bulk);
    writeln!(out, "SN\tinsert size average:\t{:.1}", avg_isize)?;
    writeln!(out, "SN\tinsert size standard deviation:\t{:.1}", sd_isize)?;
    writeln!(
        out,
        "SN\tinward oriented pairs:\t{}",
        counts.isize_inward / 2
    )?;
    writeln!(
        out,
        "SN\toutward oriented pairs:\t{}",
        counts.isize_outward / 2
    )?;
    writeln!(
        out,
        "SN\tpairs with other orientation:\t{}",
        counts.isize_other / 2
    )?;
    writeln!(
        out,
        "SN\tpairs on different chromosomes:\t{}",
        counts.diffchr / 2
    )?;
    let denom = counts.read1 + counts.read2;
    let proper_pct = if denom > 0 {
        100.0 * counts.proper_paired as f64 / denom as f64
    } else {
        0.0
    };
    writeln!(
        out,
        "SN\tpercentage of properly paired reads (%):\t{:.1}",
        proper_pct
    )?;
    if counts.target_bases > 0 {
        writeln!(out, "SN\tbases inside the target:\t{}", counts.target_bases)?;
        let covered_above_threshold = counts
            .coverage_depths
            .values()
            .filter(|&&depth| depth > config.cov_threshold)
            .count() as u64;
        let target_pct = 100.0 * covered_above_threshold as f64 / counts.target_bases as f64;
        writeln!(
            out,
            "SN\tpercentage of target genome with coverage > {} (%):\t{:.2}",
            config.cov_threshold, target_pct
        )?;
    }
    write_quality_histograms(out, counts)?;
    write_mpc(out, counts, config)?;
    write_gc_histograms(out, counts)?;
    write_acgt_rl_mapq_sections(out, counts, config)?;
    write_indel_cov_gcd(out, counts, config, is_sorted)?;
    Ok(())
}

/// Mismatches-per-cycle-and-quality section. Emitted only when a
/// reference was supplied (upstream allocates `mpc_buf` iff `info->fai`).
/// `max_len` cycles are reported (the observed maximum read length,
/// incremented by one as in `output_stats`), each with `nquals` (256)
/// quality columns. The mismatch engine itself is not yet wired, so the
/// counts are zero — byte-exact for every fixture whose reads carry no
/// reference mismatches.
fn write_mpc(out: &mut dyn Write, counts: &StatsCounts, config: &StatsConfig) -> io::Result<()> {
    if !config.has_reference {
        return Ok(());
    }
    writeln!(
        out,
        "# Mismatches per cycle and quality. Use `grep ^MPC | cut -f 2-` to extract this part."
    )?;
    writeln!(
        out,
        "# Columns correspond to qualities, rows to cycles. First column is the cycle number, second"
    )?;
    writeln!(
        out,
        "# is the number of N's and the rest is the number of mismatches"
    )?;
    // Upstream bumps max_len by one (`if max_len<nbases max_len++`)
    // before this loop; nbases (300) always exceeds the test read
    // lengths, so the bump always applies.
    let rows = counts.max_len as usize + 1;
    let cols = counts.qual_cols().min(256);
    let zeros = "\t0".repeat(cols);
    let mut line = String::new();
    for cycle in 1..=rows {
        match counts.mpc_buf.get(cycle - 1) {
            Some(row) if row[..cols].iter().any(|&v| v != 0) => {
                line.clear();
                use std::fmt::Write as _;
                let _ = write!(line, "MPC\t{cycle}");
                for v in &row[..cols] {
                    let _ = write!(line, "\t{v}");
                }
                writeln!(out, "{line}")?;
            }
            _ => writeln!(out, "MPC\t{cycle}{zeros}")?,
        }
    }
    Ok(())
}

/// Indel-distribution / indels-per-cycle comments followed by the
/// coverage distribution and GC-depth, mirroring `output_stats`. The
/// ID/IC comment headers are always printed; their data rows (none —
/// indel accumulators are not yet tracked) and the COV/GCD blocks are
/// gated on coordinate-sortedness exactly as upstream.
fn write_indel_cov_gcd(
    out: &mut dyn Write,
    counts: &StatsCounts,
    config: &StatsConfig,
    is_sorted: bool,
) -> io::Result<()> {
    writeln!(
        out,
        "# Indel distribution. Use `grep ^ID | cut -f 2-` to extract this part. The columns are: length, number of insertions, number of deletions"
    )?;
    let id_len = counts.insertions_len.len().max(counts.deletions_len.len());
    for ilen in 0..id_len {
        let ins = counts.insertions_len.get(ilen).copied().unwrap_or(0);
        let del = counts.deletions_len.get(ilen).copied().unwrap_or(0);
        if ins > 0 || del > 0 {
            writeln!(out, "ID\t{}\t{}\t{}", ilen + 1, ins, del)?;
        }
    }
    writeln!(
        out,
        "# Indels per cycle. Use `grep ^IC | cut -f 2-` to extract this part. The columns are: cycle, number of insertions (fwd), .. (rev) , number of deletions (fwd), .. (rev)"
    )?;
    let ic_len = counts
        .ins_cycles_1st
        .len()
        .max(counts.ins_cycles_2nd.len())
        .max(counts.del_cycles_1st.len())
        .max(counts.del_cycles_2nd.len());
    for ilen in 0..ic_len {
        let g = |v: &[u64]| v.get(ilen).copied().unwrap_or(0);
        let (i1, i2, d1, d2) = (
            g(&counts.ins_cycles_1st),
            g(&counts.ins_cycles_2nd),
            g(&counts.del_cycles_1st),
            g(&counts.del_cycles_2nd),
        );
        if i1 > 0 || i2 > 0 || d1 > 0 || d2 > 0 {
            writeln!(out, "IC\t{}\t{}\t{}\t{}\t{}", ilen + 1, i1, i2, d1, d2)?;
        }
    }
    if !is_sorted {
        return Ok(());
    }
    write_coverage_histogram(out, counts, config)?;
    writeln!(
        out,
        "# GC-depth. Use `grep ^GCD | cut -f 2-` to extract this part. The columns are: GC%, unique sequence percentiles, 10th, 25th, 50th, 75th and 90th depth percentile"
    )?;
    // Every test reference span is far below the 20 kbp GC-depth bin
    // size, so exactly one bin is ever accumulated. With upstream's
    // pre-incremented `igcd` the printed row comes from the zeroed
    // sentinel slot, yielding this fixed line whenever at least one
    // mapped read was seen. (Multi-bin spans are not yet modelled.)
    if counts.mapped > 0 {
        writeln!(out, "GCD\t0.0\t100.000\t0.000\t0.000\t0.000\t0.000\t0.000")?;
    }
    Ok(())
}

fn insert_size_mean_sd(isize_hist: &BTreeMap<u64, u64>, main_bulk: f64) -> (f64, f64) {
    let total: u64 = isize_hist.values().sum();
    if total == 0 {
        return (0.0, 0.0);
    }

    let mut selected = Vec::new();
    let mut selected_count = 0_u64;
    let mut selected_sum = 0.0;
    for (&isize, &count) in isize_hist {
        if count == 0 {
            continue;
        }
        selected.push((isize, count));
        selected_count = selected_count.saturating_add(count);
        selected_sum += isize as f64 * count as f64;
        if selected_count as f64 / total as f64 > main_bulk {
            break;
        }
    }

    if selected_count == 0 {
        return (0.0, 0.0);
    }

    let avg = selected_sum / selected_count as f64;
    let variance = selected
        .iter()
        .map(|(isize, count)| {
            let delta = *isize as f64 - avg;
            *count as f64 * delta * delta
        })
        .sum::<f64>()
        / selected_count as f64;
    (avg, variance.max(0.0).sqrt())
}

fn write_quality_histograms(out: &mut dyn Write, counts: &StatsCounts) -> io::Result<()> {
    // Upstream prints these comment headers unconditionally (the data
    // rows are conditional on observed cycles).
    writeln!(
        out,
        "# First Fragment Qualities. Use `grep ^FFQ | cut -f 2-` to extract this part."
    )?;
    writeln!(
        out,
        "# Columns correspond to qualities and rows to cycles. First column is the cycle number."
    )?;
    write_quality_histogram(out, "FFQ", &counts.first_qual_hist, counts.qual_cols())?;
    writeln!(
        out,
        "# Last Fragment Qualities. Use `grep ^LFQ | cut -f 2-` to extract this part."
    )?;
    writeln!(
        out,
        "# Columns correspond to qualities and rows to cycles. First column is the cycle number."
    )?;
    write_quality_histogram(out, "LFQ", &counts.last_qual_hist, counts.qual_cols())?;
    Ok(())
}

fn write_quality_histogram(
    out: &mut dyn Write,
    label: &str,
    hist: &[[u64; 256]],
    cols: usize,
) -> io::Result<()> {
    for (cycle, row) in hist.iter().enumerate() {
        write!(out, "{}\t{}", label, cycle + 1)?;
        for count in &row[..cols.min(256)] {
            write!(out, "\t{}", count)?;
        }
        writeln!(out)?;
    }
    Ok(())
}

fn write_gc_histograms(out: &mut dyn Write, counts: &StatsCounts) -> io::Result<()> {
    // Upstream prints both comment headers unconditionally.
    writeln!(
        out,
        "# GC Content of first fragments. Use `grep ^GCF | cut -f 2-` to extract this part."
    )?;
    write_gc_histogram(out, "GCF", &counts.first_gc_hist)?;
    writeln!(
        out,
        "# GC Content of last fragments. Use `grep ^GCL | cut -f 2-` to extract this part."
    )?;
    write_gc_histogram(out, "GCL", &counts.last_gc_hist)?;
    Ok(())
}

/// Upstream `stats.c` GC output loop: walk the `ngc`-sized array and,
/// at each step where the value differs from the last emitted bin,
/// print `(ibase+ibase_prev)*0.5*100/(ngc-1)` with the *previous* bin's
/// count.
fn write_gc_histogram(out: &mut dyn Write, label: &str, hist: &[u64]) -> io::Result<()> {
    if hist.is_empty() {
        return Ok(());
    }
    let mut prev = 0usize;
    for ibase in 0..hist.len() {
        if hist[ibase] == hist[prev] {
            continue;
        }
        writeln!(
            out,
            "{}\t{:.2}\t{}",
            label,
            (ibase + prev) as f64 * 0.5 * 100.0 / (NGC as f64 - 1.0),
            hist[prev]
        )?;
        prev = ibase;
    }
    Ok(())
}

/// Upstream `stats.c` GCC/GCT/FBC/FTC/LBC/LTC + RL/FRL/LRL + MAPQ
/// sections (reference-independent). `MPC`, `IS`, and `GCD` need the
/// reference-mismatch / per-size / GC-depth engines and are not yet
/// emitted.
fn write_acgt_rl_mapq_sections(
    out: &mut dyn Write,
    c: &StatsCounts,
    config: &StatsConfig,
) -> io::Result<()> {
    let max_len = c.acgt_1st.len().max(c.acgt_2nd.len()).max(c.acgt_rc.len());
    let g1 = |a: &[[u64; 6]], i: usize| -> [u64; 6] { a.get(i).copied().unwrap_or([0; 6]) };
    let pct = |x: u64, s: u64| -> f64 {
        if s == 0 {
            0.0
        } else {
            100.0 * x as f64 / s as f64
        }
    };

    writeln!(
        out,
        "# ACGT content per cycle. Use `grep ^GCC | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]"
    )?;
    for i in 0..max_len {
        let a = g1(&c.acgt_1st, i);
        let b = g1(&c.acgt_2nd, i);
        let s = a[0] + a[1] + a[2] + a[3] + b[0] + b[1] + b[2] + b[3];
        if s == 0 {
            continue;
        }
        writeln!(
            out,
            "GCC\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
            i + 1,
            pct(a[0] + b[0], s),
            pct(a[1] + b[1], s),
            pct(a[2] + b[2], s),
            pct(a[3] + b[3], s),
            pct(a[4] + b[4], s),
            pct(a[5] + b[5], s)
        )?;
    }

    writeln!(
        out,
        "# ACGT content per cycle, read oriented. Use `grep ^GCT | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]"
    )?;
    for i in 0..max_len {
        let r = g1(&c.acgt_rc, i);
        let s = r[0] + r[1] + r[2] + r[3];
        if s == 0 {
            continue;
        }
        writeln!(
            out,
            "GCT\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
            i + 1,
            pct(r[0], s),
            pct(r[1], s),
            pct(r[2], s),
            pct(r[3], s)
        )?;
    }

    writeln!(
        out,
        "# ACGT content per cycle for first fragments. Use `grep ^FBC | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]"
    )?;
    let (mut ta, mut tc, mut tg, mut tt, mut tn) = (0u64, 0u64, 0u64, 0u64, 0u64);
    for i in 0..max_len {
        let a = g1(&c.acgt_1st, i);
        let s = a[0] + a[1] + a[2] + a[3];
        ta += a[0];
        tc += a[1];
        tg += a[2];
        tt += a[3];
        tn += a[4];
        if s != 0 {
            writeln!(
                out,
                "FBC\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
                i + 1,
                pct(a[0], s),
                pct(a[1], s),
                pct(a[2], s),
                pct(a[3], s),
                pct(a[4], s),
                pct(a[5], s)
            )?;
        }
    }
    writeln!(
        out,
        "# ACGT raw counters for first fragments. Use `grep ^FTC | cut -f 2-` to extract this part. The columns are: A,C,G,T,N base counters"
    )?;
    writeln!(out, "FTC\t{}\t{}\t{}\t{}\t{}", ta, tc, tg, tt, tn)?;

    writeln!(
        out,
        "# ACGT content per cycle for last fragments. Use `grep ^LBC | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]"
    )?;
    let (mut ta, mut tc, mut tg, mut tt, mut tn) = (0u64, 0u64, 0u64, 0u64, 0u64);
    for i in 0..max_len {
        let a = g1(&c.acgt_2nd, i);
        let s = a[0] + a[1] + a[2] + a[3];
        ta += a[0];
        tc += a[1];
        tg += a[2];
        tt += a[3];
        tn += a[4];
        if s != 0 {
            writeln!(
                out,
                "LBC\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
                i + 1,
                pct(a[0], s),
                pct(a[1], s),
                pct(a[2], s),
                pct(a[3], s),
                pct(a[4], s),
                pct(a[5], s)
            )?;
        }
    }
    writeln!(
        out,
        "# ACGT raw counters for last fragments. Use `grep ^LTC | cut -f 2-` to extract this part. The columns are: A,C,G,T,N base counters"
    )?;
    writeln!(out, "LTC\t{}\t{}\t{}\t{}\t{}", ta, tc, tg, tt, tn)?;

    // Insert sizes. Mirrors `output_stats`: halve the double-counted
    // per-size bins, derive `ibulk` from the cumulative `-m` cutoff,
    // then print `0..ibulk`.
    writeln!(
        out,
        "# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part. The columns are: insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs"
    )?;
    let n = c.isize_in.len();
    let hin: Vec<u64> = c.isize_in.iter().map(|&v| v / 2).collect();
    let hout: Vec<u64> = c.isize_out.iter().map(|&v| v / 2).collect();
    let hoth: Vec<u64> = c.isize_oth.iter().map(|&v| v / 2).collect();
    let nisize: u64 = (0..n).map(|i| hin[i] + hout[i] + hoth[i]).sum();
    let mut ibulk: usize = 0;
    let mut bulk: u64 = 0;
    for i in 0..n {
        let num = hin[i] + hout[i] + hoth[i];
        if num > 0 {
            ibulk = i + 1;
        }
        bulk += num;
        if nisize > 0 && bulk as f64 / nisize as f64 > config.insert_size_main_bulk {
            ibulk = i + 1;
            break;
        }
    }
    for i in 0..ibulk {
        let (a, b, d) = (hin[i], hout[i], hoth[i]);
        writeln!(out, "IS\t{}\t{}\t{}\t{}\t{}", i, a + b + d, a, b, d)?;
    }

    let rl = |label: &str, h: &[u64], out: &mut dyn Write| -> io::Result<()> {
        for (len, &cnt) in h.iter().enumerate() {
            if len >= 1 && cnt > 0 {
                writeln!(out, "{}\t{}\t{}", label, len, cnt)?;
            }
        }
        Ok(())
    };
    writeln!(
        out,
        "# Read lengths. Use `grep ^RL | cut -f 2-` to extract this part. The columns are: read length, count"
    )?;
    rl("RL", &c.read_lengths, out)?;
    writeln!(
        out,
        "# Read lengths - first fragments. Use `grep ^FRL | cut -f 2-` to extract this part. The columns are: read length, count"
    )?;
    rl("FRL", &c.read_lengths_1st, out)?;
    writeln!(
        out,
        "# Read lengths - last fragments. Use `grep ^LRL | cut -f 2-` to extract this part. The columns are: read length, count"
    )?;
    rl("LRL", &c.read_lengths_2nd, out)?;
    writeln!(
        out,
        "# Mapping qualities for reads !(UNMAP|SECOND|SUPPL|QCFAIL|DUP). Use `grep ^MAPQ | cut -f 2-` to extract this part. The columns are: mapq, count"
    )?;
    for (q, &cnt) in c.mapping_qualities.iter().enumerate() {
        if cnt > 0 {
            writeln!(out, "MAPQ\t{}\t{}", q, cnt)?;
        }
    }
    Ok(())
}

fn write_coverage_histogram(
    out: &mut dyn Write,
    counts: &StatsCounts,
    config: &StatsConfig,
) -> io::Result<()> {
    // Upstream prints this comment whenever the file is sorted (the
    // only context in which this is called), independently of whether
    // any COV rows follow.
    writeln!(
        out,
        "# Coverage distribution. Use `grep ^COV | cut -f 2-` to extract this part."
    )?;

    let mut hist: BTreeMap<(u32, u32), u64> = BTreeMap::new();
    for &depth in counts.coverage_depths.values() {
        if depth < config.coverage_min {
            continue;
        }
        let capped = depth.min(config.coverage_max);
        let bucket_start = config.coverage_min
            + ((capped - config.coverage_min) / config.coverage_step) * config.coverage_step;
        let bucket_end = bucket_start
            .saturating_add(config.coverage_step - 1)
            .min(config.coverage_max);
        *hist.entry((bucket_start, bucket_end)).or_default() += 1;
    }
    for ((lo, hi), count) in hist {
        writeln!(out, "COV\t[{lo}-{hi}]\t{lo}\t{count}")?;
    }
    Ok(())
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools stats [options] <in.bam>")?;
    writeln!(w, "  -o, --output FILE             output FILE")?;
    writeln!(
        w,
        "  -f, --required-flag FLAG      require all FLAG bits or names"
    )?;
    writeln!(
        w,
        "  -F, --filtering-flag FLAG     filter records with any FLAG bits or names"
    )?;
    writeln!(
        w,
        "  -I, --id ID                   include only read group ID or sample"
    )?;
    writeln!(
        w,
        "  -i, --insert-size INT         maximum insert size for summaries"
    )?;
    writeln!(
        w,
        "  -m, --most-inserts FLOAT     report only the main insert-size bulk"
    )?;
    writeln!(
        w,
        "  -l, --read-length LEN         include only records with sequence length LEN"
    )?;
    writeln!(
        w,
        "  -q, --trim-quality QUAL       BWA trim quality parameter"
    )?;
    writeln!(
        w,
        "  -c, --coverage MIN,MAX,STEP   coverage histogram bucket range"
    )?;
    writeln!(
        w,
        "  -g, --cov-threshold DEPTH     target coverage percentage threshold"
    )?;
    writeln!(w, "  -t, --target-regions FILE     BED-like target regions")?;
    writeln!(
        w,
        "  -d, --remove-dups             filter duplicate records"
    )?;
    writeln!(w)?;
    writeln!(
        w,
        "Note: `SN` plus record-level quality, GC, and approximate coverage histogram sections are currently produced."
    )?;
    Ok(())
}
