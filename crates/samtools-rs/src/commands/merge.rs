//! `samtools merge` — merge multiple sorted BAM files.
//!
//! Mirrors `bam_merge` in `bam_sort.c`. This initial Rust port loads all
//! records from BAM/SAM inputs into memory and sorts by coordinate (or name
//! with `-n`) before writing the merged output. K-way streaming merge
//! is TODO. CRAM output is supported with a top-level `--reference`. `-R` and
//! `-L` restrict indexed BAM inputs by region/BED. Coordinate-sorted BAM
//! outputs can also write a BAI via `--write-index`.

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf};

use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools merge`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut name_sort = false;
    let mut output: Option<PathBuf> = None;
    let mut output_fmt = OutFmt::Bam;
    let mut positional: Vec<PathBuf> = Vec::new();
    let mut force = false;
    let mut local_write_index = false;
    let mut no_pg = false;
    let mut region: Option<String> = None;
    let mut bed: Option<PathBuf> = None;
    let mut tag_sort: Option<[u8; 2]> = None;
    let mut template_coordinate = false;
    let mut input_lists: Vec<PathBuf> = Vec::new();
    // `-s SEED` seeds the @RG/@PG ID-collision PRNG (default: HTSlib uses
    // time(); 0 keeps merges deterministic when unspecified).
    let mut seed: i64 = 0;
    // `-c` combine identical @RG IDs (don't suffix); `-p` same for @PG;
    // `-r` attach a filename-derived @RG to every record.
    let mut combine_rg = false;
    let mut combine_pg = false;
    let mut attach_rg = false;
    let mut customized_index = false;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        if let Some(v) = s.strip_prefix("--output-fmt=") {
            output_fmt = match parse_output_format(v) {
                Ok(fmt) => fmt,
                Err(e) => {
                    print_error("merge", e);
                    return ExitCode::from(1);
                }
            };
            continue;
        }
        match s {
            "-n" => name_sort = true,
            "-f" => force = true,
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "--output-fmt" | "-O" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("merge", format!("missing value for {}", s));
                    return ExitCode::from(1);
                };
                output_fmt = match parse_output_format(v) {
                    Ok(fmt) => fmt,
                    Err(e) => {
                        print_error("merge", e);
                        return ExitCode::from(1);
                    }
                };
            }
            "-t" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("merge", "missing value for -t");
                    return ExitCode::from(1);
                };
                tag_sort = match parse_tag(v) {
                    Ok(tag) => Some(tag),
                    Err(e) => {
                        print_error("merge", e);
                        return ExitCode::from(1);
                    }
                };
            }
            "-R" => {
                region = iter.next().and_then(|a| a.to_str().map(str::to_owned));
            }
            "-L" => {
                bed = iter.next().map(PathBuf::from);
            }
            "-X" | "--customized-index" => {
                customized_index = true;
            }
            "--template-coordinate" => {
                template_coordinate = true;
            }
            "-b" => {
                let Some(path) = iter.next() else {
                    print_error("merge", "missing value for -b");
                    return ExitCode::from(1);
                };
                input_lists.push(PathBuf::from(path));
            }
            "-@" | "--threads" | "-l" | "--compression-level" => {
                let _ = iter.next();
            }
            "--write-index" => local_write_index = true,
            "--no-PG" => {
                no_pg = true;
            }
            "-s" => {
                if let Some(v) = iter.next().and_then(|a| a.to_str())
                    && let Ok(n) = v.parse::<i64>()
                {
                    seed = n;
                }
            }
            "-u" => {}
            // `-c`/`-p`/`-r` (also grouped, e.g. `-cp`, `-rp`).
            _ if s.len() >= 2
                && s.starts_with('-')
                && !s.starts_with("--")
                && s[1..]
                    .bytes()
                    .all(|b| matches!(b, b'c' | b'p' | b'r' | b'u')) =>
            {
                for b in s[1..].bytes() {
                    match b {
                        b'c' => combine_rg = true,
                        b'p' => combine_pg = true,
                        b'r' => attach_rg = true,
                        _ => {}
                    }
                }
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error(
                    "merge",
                    format!("option `{}` is not yet supported in samtools-rs merge", s),
                );
                return ExitCode::from(1);
            }
            _ => positional.push(PathBuf::from(arg)),
        }
    }

    // Upstream synopsis: `samtools merge [options] <out.bam> <in1.bam> [<in2.bam>...]`.
    // If `-o` is given, all positionals are inputs; otherwise the first
    // positional is the output path.
    let (out_path, mut inputs): (Option<PathBuf>, Vec<PathBuf>) = if output.is_some() {
        (output, positional)
    } else if positional.is_empty() {
        let _ = print_usage();
        return ExitCode::from(1);
    } else {
        let mut iter = positional.into_iter();
        let out = iter.next();
        let inputs: Vec<_> = iter.collect();
        (out, inputs)
    };

    let mut listed_inputs = Vec::new();
    for list in &input_lists {
        match read_input_list(list) {
            Ok(list_inputs) => listed_inputs.extend(list_inputs),
            Err(e) => {
                print_error("merge", e.to_string());
                return ExitCode::from(1);
            }
        }
    }
    if !listed_inputs.is_empty() {
        listed_inputs.extend(inputs);
        inputs = listed_inputs;
    }

    let custom_index_paths = if customized_index {
        if inputs.len() % 2 != 0 {
            print_error(
                "merge",
                "-X requires one custom index path for each input alignment",
            );
            return ExitCode::from(1);
        }
        let index_paths = inputs.split_off(inputs.len() / 2);
        Some(index_paths)
    } else {
        None
    };

    if inputs.is_empty() {
        let _ = print_usage();
        return ExitCode::from(1);
    }

    let out_path = out_path.filter(|p| p.as_os_str() != "-");

    if let Some(p) = out_path.as_ref()
        && p.exists()
        && !force
    {
        print_error(
            "merge",
            format!(
                "output file \"{}\" exists. Use -f to overwrite.",
                p.display()
            ),
        );
        return ExitCode::from(1);
    }

    for path in &inputs {
        if path.as_os_str() != "-" && !path.exists() {
            print_hts_open_missing(path);
            print_error(
                "merge",
                format!(
                    "fail to open \"{}\": No such file or directory",
                    path.display()
                ),
            );
            return ExitCode::from(1);
        }
        let format = match sam_io::sam_open_format(path) {
            Ok(f) => f,
            Err(e) => {
                print_error("merge", e.to_string());
                return ExitCode::from(1);
            }
        };
        if !matches!(format.exact, Exact::Sam | Exact::Bam) {
            print_error(
                "merge",
                format!(
                    "only SAM and BAM input are currently supported (got {:?} for \"{}\")",
                    format.exact,
                    path.display()
                ),
            );
            return ExitCode::from(1);
        }
    }

    let write_index = local_write_index || current_global_args().write_index;
    let order = match (template_coordinate, tag_sort) {
        (true, _) => MergeOrder::TemplateCoordinate,
        (false, Some(tag)) => MergeOrder::Tag {
            tag,
            name_secondary: name_sort,
        },
        (false, None) if name_sort => MergeOrder::Name,
        (false, None) => MergeOrder::Coordinate,
    };

    if write_index {
        if out_path.is_none() {
            print_error("merge", "--write-index requires output file");
            return ExitCode::from(1);
        }
        if !matches!(order, MergeOrder::Coordinate) {
            print_error("merge", "--write-index requires coordinate sort output");
            return ExitCode::from(1);
        }
        if !matches!(output_fmt, OutFmt::Bam) {
            print_error("merge", "--write-index is only supported for BAM output");
            return ExitCode::from(1);
        }
    }

    if region.is_some() && bed.is_some() {
        print_error(
            "merge",
            "-R and -L are mutually exclusive in samtools-rs merge",
        );
        return ExitCode::from(1);
    }

    match run_merge(
        &inputs,
        out_path.as_deref(),
        order,
        output_fmt,
        MergeIdMode {
            combine_rg,
            combine_pg,
            attach_rg,
        },
        write_index,
        if no_pg { None } else { Some(args) },
        seed,
        custom_index_paths.as_deref(),
        match (region.as_deref(), bed.as_deref()) {
            (Some(r), None) => MergeRestriction::Region(r),
            (None, Some(path)) => MergeRestriction::Bed(path),
            _ => MergeRestriction::None,
        },
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("merge", "merge failed", &e);
            ExitCode::from(1)
        }
    }
}

fn parse_output_format(raw: &str) -> Result<OutFmt, String> {
    match raw.to_ascii_lowercase().as_str() {
        "sam" => Ok(OutFmt::Sam),
        "bam" => Ok(OutFmt::Bam),
        "cram" => Ok(OutFmt::Cram),
        _ => Err(format!("unsupported output format \"{}\"", raw)),
    }
}

fn parse_tag(raw: &str) -> Result<[u8; 2], String> {
    let bytes = raw.as_bytes();
    if bytes.len() == 2 {
        Ok([bytes[0], bytes[1]])
    } else {
        Err(format!(
            "merge tag must be exactly two bytes, got {:?}",
            raw
        ))
    }
}

fn read_input_list(path: &Path) -> io::Result<Vec<PathBuf>> {
    let file = File::open(path)?;
    let mut inputs = Vec::new();

    for line in BufReader::new(file).lines() {
        let line = line?;
        let line = line.trim();
        if !line.is_empty() {
            inputs.push(PathBuf::from(line));
        }
    }

    Ok(inputs)
}

#[derive(Clone, Copy)]
pub(crate) enum OutFmt {
    Sam,
    Bam,
    Cram,
}

pub(crate) enum MergeRestriction<'a> {
    None,
    Region(&'a str),
    Bed(&'a Path),
}

#[derive(Clone, Copy)]
pub(crate) enum MergeOrder {
    Coordinate,
    Name,
    Tag { tag: [u8; 2], name_secondary: bool },
    TemplateCoordinate,
}

/// `-c` (combine identical @RG IDs), `-p` (same for @PG), `-r` (attach a
/// filename-derived @RG to every record / one @RG per input).
#[derive(Clone, Copy, Default)]
pub(crate) struct MergeIdMode {
    pub combine_rg: bool,
    pub combine_pg: bool,
    pub attach_rg: bool,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_merge(
    inputs: &[PathBuf],
    output: Option<&Path>,
    order: MergeOrder,
    fmt: OutFmt,
    id_mode: MergeIdMode,
    write_index: bool,
    pg_argv: Option<&[OsString]>,
    seed: i64,
    custom_index_paths: Option<&[PathBuf]>,
    restriction: MergeRestriction<'_>,
) -> io::Result<()> {
    if let Some(paths) = custom_index_paths
        && paths.len() != inputs.len()
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "-X requires one custom index path for each input alignment",
        ));
    }

    let filter = match restriction {
        MergeRestriction::Region(r) => Some(RegionFilter::Regions(vec![
            r.parse::<htslib_rs::core::Region>().map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid -R region \"{r}\": {e}"),
                )
            })?,
        ])),
        MergeRestriction::Bed(path) => {
            let bed = crate::bedidx::load_bed_index(path)?;
            Some(RegionFilter::Regions(bed.to_htslib_regions()?))
        }
        MergeRestriction::None => None,
    };

    // Byte-faithful merged header text + per-file @RG/@PG ID translation
    // (upstream `samtools merge` PRNG-suffixes colliding IDs).
    let (merged_header_text, trans, force_rg) = reconcile_merge_headers(inputs, seed, id_mode)?;

    let (mut header, mut records) = read_records(
        &inputs[0],
        filter.as_ref(),
        custom_index_paths
            .and_then(|paths| paths.first())
            .map(PathBuf::as_path),
    )?;
    let first_reference_id_map: Vec<_> = (0..header.reference_sequences().len()).collect();
    remap_records(&mut records, &first_reference_id_map)?;
    let mut warned_missing_tags = HashSet::new();
    for rec in &mut records {
        remap_record_rg_pg(
            rec,
            &trans.rg[0],
            &trans.pg[0],
            force_rg[0].as_deref(),
            &mut warned_missing_tags,
        );
    }

    for (idx, path) in inputs[1..].iter().enumerate() {
        let (input_header, mut input_records) = read_records(
            path,
            filter.as_ref(),
            custom_index_paths
                .and_then(|paths| paths.get(idx + 1))
                .map(PathBuf::as_path),
        )?;
        let reference_id_map = merge_reference_sequences(&mut header, &input_header)?;
        merge_header_metadata(&mut header, &input_header)?;
        merge_read_groups(&mut header, &input_header)?;
        merge_programs(&mut header, &input_header)?;
        merge_comments(&mut header, &input_header);
        remap_records(&mut input_records, &reference_id_map)?;
        let mut warned_missing_tags = HashSet::new();
        for rec in &mut input_records {
            remap_record_rg_pg(
                rec,
                &trans.rg[idx + 1],
                &trans.pg[idx + 1],
                force_rg[idx + 1].as_deref(),
                &mut warned_missing_tags,
            );
        }
        records.append(&mut input_records);
    }

    match order {
        MergeOrder::Tag {
            tag,
            name_secondary,
        } => records.sort_by(|a, b| compare_by_tag(a, b, tag, name_secondary)),
        MergeOrder::TemplateCoordinate => {
            let lib_lookup = super::sort::library_lookup(&header);
            let mut keyed = Vec::with_capacity(records.len());
            for record in records.drain(..) {
                let key = super::sort::template_coordinate_key(&record, &lib_lookup)
                    .map_err(io::Error::other)?;
                keyed.push((key, record));
            }
            keyed.sort_by(|a, b| super::sort::template_coordinate_cmp(&a.0, &b.0));
            records = keyed.into_iter().map(|(_, record)| record).collect();
        }
        MergeOrder::Name => records.sort_by(natural_name_cmp),
        MergeOrder::Coordinate => records.sort_by_key(coordinate_key),
    }

    match order {
        MergeOrder::Tag {
            tag,
            name_secondary,
        } => set_sort_order(
            &mut header,
            "unsorted",
            None,
            Some(&format!(
                "unsorted:{}{}:{}",
                tag[0] as char,
                tag[1] as char,
                if name_secondary {
                    "queryname:natural"
                } else {
                    "coordinate"
                }
            )),
        ),
        MergeOrder::TemplateCoordinate => {
            set_sort_order(
                &mut header,
                "unsorted",
                Some("query"),
                Some("unsorted:template-coordinate"),
            );
        }
        MergeOrder::Name => {
            set_sort_order(&mut header, "queryname", None, Some("queryname:natural"))
        }
        MergeOrder::Coordinate => set_sort_order(&mut header, "coordinate", None, None),
    }

    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

    if matches!(fmt, OutFmt::Sam) {
        // SAM: emit the byte-faithful reconciled header text (preserving
        // @RG/@SQ field order and the input @HD verbatim) + the samtools
        // @PG, then records via the shared float-correct renderer.
        // Coordinate merge keeps input[0]'s @HD verbatim (no SO — matches
        // the upstream merge/* fixtures); name/tag merge sets SO/SS.
        let mut header_text = match order {
            MergeOrder::Coordinate => merged_header_text,
            MergeOrder::Name => apply_hd_so(
                &merged_header_text,
                "queryname",
                None,
                Some("queryname:natural"),
            ),
            MergeOrder::TemplateCoordinate => apply_hd_so(
                &merged_header_text,
                "unsorted",
                Some("query"),
                Some("unsorted:template-coordinate"),
            ),
            MergeOrder::Tag {
                tag,
                name_secondary,
            } => apply_hd_so(
                &merged_header_text,
                "unsorted",
                None,
                Some(&format!(
                    "unsorted:{}{}:{}",
                    tag[0] as char,
                    tag[1] as char,
                    if name_secondary {
                        "queryname:natural"
                    } else {
                        "coordinate"
                    }
                )),
            ),
        };
        if let Some(argv) = pg_argv {
            header_text =
                crate::pg::add_samtools_pg(&header_text, argv).map_err(io::Error::other)?;
        }
        let mut out: Box<dyn Write> = match output {
            Some(p) => Box::new(io::BufWriter::new(File::create(p)?)),
            None => Box::new(io::BufWriter::new(io::stdout().lock())),
        };
        out.write_all(header_text.as_bytes())?;
        for rec in &records {
            crate::sam_render::write_record(&mut out, &header, rec)?;
        }
        out.flush()?;
    } else {
        let mut writer = open_output(output, fmt, &header)?;
        for rec in &records {
            writer.write_record(&header, rec)?;
        }
        writer.finish()?;
    }

    if write_index {
        let Some(path) = output else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--write-index requires output file",
            ));
        };
        write_bam_index(path)?;
    }
    Ok(())
}

enum RegionFilter {
    Regions(Vec<htslib_rs::core::Region>),
}

fn read_records(
    input: &Path,
    filter: Option<&RegionFilter>,
    custom_index_path: Option<&Path>,
) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let format = sam_io::sam_open_format(input)?;
    match (format.exact, filter) {
        (Exact::Sam, None) => read_sam_records(input),
        (Exact::Sam, Some(_)) => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "samtools merge region filters require indexed BAM input (SAM is not supported)",
        )),
        (Exact::Bam, None) => read_bam_records(input),
        (Exact::Bam, Some(RegionFilter::Regions(regions))) => {
            read_bam_records_in_regions(input, custom_index_path, regions)
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "only SAM and BAM input are currently supported (CRAM TODO)",
        )),
    }
}

fn read_bam_records(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;
    let mut records = Vec::new();
    loop {
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        records.push(record);
    }
    Ok((header, records))
}

/// Reads BAM records overlapping `regions` using the input's BAI index.
/// Returns the records as `RecordBuf` for downstream sorting and writing.
fn read_bam_records_in_regions(
    input: &Path,
    custom_index_path: Option<&Path>,
    regions: &[htslib_rs::core::Region],
) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let header = htslib_rs::alignment_compat::read_bam_header_from_path(input)?;
    let query_input = match custom_index_path {
        Some(index_path) => explicit_indexed_path(input, index_path),
        None => input.to_path_buf(),
    };
    let regions = canonicalize_regions_for_header(&header, regions)?;
    let mut records = Vec::new();
    let mut seen = HashSet::new();
    for region in &regions {
        let bam_records =
            htslib_rs::alignment_compat::query_bam_records_from_path(&query_input, region)?;
        for bam_record in bam_records {
            let buf = sam::alignment::RecordBuf::try_from_alignment_record(&header, &bam_record)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            if seen.insert(record_key(&buf)) {
                records.push(buf);
            }
        }
    }
    Ok((header, records))
}

fn explicit_indexed_path(input: &Path, index: &Path) -> PathBuf {
    PathBuf::from(format!(
        "{}{}{}",
        input.to_string_lossy(),
        htslib_rs::index_compat::INDEX_DELIMITER,
        index.to_string_lossy()
    ))
}

fn canonicalize_regions_for_header(
    header: &sam::Header,
    regions: &[htslib_rs::core::Region],
) -> io::Result<Vec<htslib_rs::core::Region>> {
    regions
        .iter()
        .map(|region| {
            let name = canonical_reference_name(header, region.name()).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "region reference sequence does not exist: {}",
                        region.name()
                    ),
                )
            })?;

            Ok(htslib_rs::core::Region::new(name, region.interval()))
        })
        .collect()
}

fn canonical_reference_name(header: &sam::Header, name: &[u8]) -> Option<Vec<u8>> {
    let reference_sequences = header.reference_sequences();

    if let Some(i) = reference_sequences.get_index_of(name) {
        let (canonical_name, _) = reference_sequences.get_index(i)?;
        return Some(canonical_name.to_vec());
    }

    reference_sequences
        .iter()
        .find_map(|(canonical_name, reference_sequence)| {
            reference_sequence
                .other_fields()
                .get(&sam::header::record::value::map::reference_sequence::tag::ALTERNATIVE_NAMES)
                .filter(|alternative_names| {
                    alternative_names
                        .as_slice()
                        .split(|&b| b == b',')
                        .any(|alternative_name| alternative_name == name)
                })
                .map(|_| canonical_name.to_vec())
        })
}

fn record_key(record: &RecordBuf) -> (Vec<u8>, u16, Option<usize>, Option<usize>, String) {
    (
        record.name().map(|n| n.to_vec()).unwrap_or_default(),
        record.flags().bits(),
        record.reference_sequence_id(),
        record.alignment_start().map(usize::from),
        format!("{:?}", record.cigar().as_ref()),
    )
}

fn merge_reference_sequences(
    output_header: &mut sam::Header,
    input_header: &sam::Header,
) -> io::Result<Vec<usize>> {
    let mut reference_id_map = Vec::with_capacity(input_header.reference_sequences().len());

    for (name, input_reference_sequence) in input_header.reference_sequences() {
        let existing_index = output_header
            .reference_sequences()
            .iter()
            .position(|(existing_name, _)| existing_name == name);

        let new_id = if let Some(index) = existing_index {
            let (_, output_reference_sequence) = output_header
                .reference_sequences_mut()
                .get_index_mut(index)
                .expect("reference index from position must exist");

            if output_reference_sequence.length() != input_reference_sequence.length() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "conflicting @SQ length for {}: {} != {}",
                        String::from_utf8_lossy(name),
                        usize::from(output_reference_sequence.length()),
                        usize::from(input_reference_sequence.length())
                    ),
                ));
            }

            merge_reference_sequence_metadata(
                name,
                output_reference_sequence,
                input_reference_sequence,
            )?;

            index
        } else {
            let index = output_header.reference_sequences().len();
            output_header
                .reference_sequences_mut()
                .insert(name.clone(), input_reference_sequence.clone());
            index
        };

        reference_id_map.push(new_id);
    }

    Ok(reference_id_map)
}

fn merge_reference_sequence_metadata(
    name: &[u8],
    output_reference_sequence: &mut sam::header::record::value::Map<
        sam::header::record::value::map::ReferenceSequence,
    >,
    input_reference_sequence: &sam::header::record::value::Map<
        sam::header::record::value::map::ReferenceSequence,
    >,
) -> io::Result<()> {
    for (tag, input_value) in input_reference_sequence.other_fields() {
        if let Some(output_value) = output_reference_sequence.other_fields().get(tag) {
            if output_value != input_value {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "conflicting @SQ {} field for {}",
                        tag,
                        String::from_utf8_lossy(name)
                    ),
                ));
            }
        } else {
            output_reference_sequence
                .other_fields_mut()
                .insert(*tag, input_value.clone());
        }
    }

    Ok(())
}

fn merge_header_metadata(
    output_header: &mut sam::Header,
    input_header: &sam::Header,
) -> io::Result<()> {
    use sam::header::record::value::map::header::tag::{SORT_ORDER, SUBSORT_ORDER};

    let Some(input_hd) = input_header.header() else {
        return Ok(());
    };

    let Some(output_hd) = output_header.header_mut() else {
        *output_header.header_mut() = Some(input_hd.clone());
        return Ok(());
    };

    for (tag, input_value) in input_hd.other_fields() {
        if *tag == SORT_ORDER || *tag == SUBSORT_ORDER {
            continue;
        }

        if let Some(output_value) = output_hd.other_fields().get(tag) {
            if output_value != input_value {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("conflicting @HD {} field", tag),
                ));
            }
        } else {
            output_hd
                .other_fields_mut()
                .insert(*tag, input_value.clone());
        }
    }

    Ok(())
}

fn remap_records(records: &mut [RecordBuf], reference_id_map: &[usize]) -> io::Result<()> {
    for record in records {
        if let Some(tid) = record.reference_sequence_id() {
            let Some(&new_tid) = reference_id_map.get(tid) else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("record reference id {tid} is not present in its input header"),
                ));
            };
            *record.reference_sequence_id_mut() = Some(new_tid);
        }

        if let Some(tid) = record.mate_reference_sequence_id() {
            let Some(&new_tid) = reference_id_map.get(tid) else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("record mate reference id {tid} is not present in its input header"),
                ));
            };
            *record.mate_reference_sequence_id_mut() = Some(new_tid);
        }
    }

    Ok(())
}

fn merge_read_groups(
    output_header: &mut sam::Header,
    input_header: &sam::Header,
) -> io::Result<()> {
    // Reconciliation (PRNG ID-suffixing) is done on the *raw* header
    // text for SAM output; the noodles header is only used for record
    // ref-id remap / BAM output, so just keep the first definition of a
    // colliding @RG ID instead of erroring.
    for (id, input_read_group) in input_header.read_groups() {
        if output_header.read_groups().get(id).is_none() {
            output_header
                .read_groups_mut()
                .insert(id.clone(), input_read_group.clone());
        }
    }

    Ok(())
}

fn merge_programs(output_header: &mut sam::Header, input_header: &sam::Header) -> io::Result<()> {
    for (id, input_program) in input_header.programs().as_ref() {
        if output_header.programs().as_ref().get(id).is_none() {
            output_header
                .programs_mut()
                .as_mut()
                .insert(id.clone(), input_program.clone());
        }
    }

    Ok(())
}

/// Sets the `@HD` `SO:`/`GO:`/`SS:` fields in raw header text (preserving
/// the other fields/lines verbatim); used for name/tag/template merges.
fn apply_hd_so(raw: &str, so: &str, go: Option<&str>, ss: Option<&str>) -> String {
    let mut lines: Vec<String> = Vec::new();
    let mut had_hd = false;
    for line in raw.lines() {
        if line.starts_with("@HD") {
            had_hd = true;
            let kept: Vec<&str> = line
                .split('\t')
                .filter(|f| !f.starts_with("SO:") && !f.starts_with("GO:") && !f.starts_with("SS:"))
                .collect();
            let mut nl = kept.join("\t");
            nl.push_str(&format!("\tSO:{so}"));
            if let Some(go) = go {
                nl.push_str(&format!("\tGO:{go}"));
            }
            if let Some(ss) = ss {
                nl.push_str(&format!("\tSS:{ss}"));
            }
            lines.push(nl);
        } else {
            lines.push(line.to_string());
        }
    }
    if !had_hd {
        let go = go.map(|go| format!("\tGO:{go}")).unwrap_or_default();
        let hd = match ss {
            Some(ss) => format!("@HD\tVN:1.6\tSO:{so}{go}\tSS:{ss}"),
            None => format!("@HD\tVN:1.6\tSO:{so}{go}"),
        };
        lines.insert(0, hd);
    }
    let mut s = lines.join("\n");
    s.push('\n');
    s
}

/// Value of a `TAG:value` field within a tab-split header line.
fn header_field<'a>(line: &'a str, tag: &str) -> Option<&'a str> {
    line.split('\t').skip(1).find_map(|f| f.strip_prefix(tag))
}

/// Rewrites the `tag` field's value in a header line via `map` (only if
/// the field exists and a mapping is present), preserving field order.
fn rewrite_field(line: &str, tag: &str, map: &std::collections::HashMap<String, String>) -> String {
    line.split('\t')
        .map(|f| {
            if let Some(v) = f.strip_prefix(tag)
                && let Some(n) = map.get(v)
            {
                format!("{tag}{n}")
            } else {
                f.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\t")
}

/// Per-file `@RG`/`@PG` ID translation maps.
struct MergeTrans {
    rg: Vec<std::collections::HashMap<String, String>>,
    pg: Vec<std::collections::HashMap<String, String>>,
}

/// Builds the byte-faithful merged header text plus per-file `@RG`/`@PG`
/// ID translation maps, mirroring upstream `samtools merge` (`bam_sort.c`):
/// `@HD` from input[0] verbatim, `@SQ` unioned by `SN`, colliding
/// `@RG`/`@PG` IDs suffixed via the seeded `gen_unique_id` PRNG (in
/// header-line order, per file), `@RG.PG`/`@PG.PP` fields remapped, `@CO`
/// appended. `@PG` lines are emitted but the test harness strips them.
#[allow(clippy::type_complexity)]
fn reconcile_merge_headers(
    inputs: &[PathBuf],
    seed: i64,
    id_mode: MergeIdMode,
) -> io::Result<(String, MergeTrans, Vec<Option<String>>)> {
    use std::collections::{HashMap, HashSet};

    let mut rng = crate::rand48::Rand48::new(seed);
    let mut seen_rg: HashSet<String> = HashSet::new();
    let mut seen_pg: HashSet<String> = HashSet::new();
    // Final IDs already emitted, so `-c`/`-p` combine (don't duplicate).
    let mut emitted_rg: HashSet<String> = HashSet::new();
    let mut emitted_pg: HashSet<String> = HashSet::new();

    let mut hd: Option<String> = None;
    let mut seen_sn: HashSet<String> = HashSet::new();
    let mut sq: Vec<String> = Vec::new();
    let mut rg_lines: Vec<String> = Vec::new();
    let mut pg_lines: Vec<String> = Vec::new();
    let mut co: Vec<String> = Vec::new();

    let mut trans = MergeTrans {
        rg: Vec::with_capacity(inputs.len()),
        pg: Vec::with_capacity(inputs.len()),
    };
    let mut force_rg: Vec<Option<String>> = Vec::with_capacity(inputs.len());

    // `id, seen, combine` → final id (combine: keep existing, no PRNG draw).
    let assign = |id: &str,
                  seen: &mut HashSet<String>,
                  rng: &mut crate::rand48::Rand48,
                  combine: bool|
     -> String {
        if combine {
            seen.insert(id.to_string());
            id.to_string()
        } else {
            crate::rand48::gen_unique_id(id, seen, rng)
        }
    };

    for path in inputs {
        let exact = sam_io::sam_open_format(path)?.exact;
        let raw = crate::header_text::read_raw_header_text_with_format(path, exact)?;
        let mut rg_map: HashMap<String, String> = HashMap::new();
        let mut pg_map: HashMap<String, String> = HashMap::new();

        // `-r`: one filename-stem @RG per file (from its first @RG line),
        // all records forced to it.
        let stem = if id_mode.attach_rg {
            path.file_stem().map(|s| s.to_string_lossy().into_owned())
        } else {
            None
        };

        // Pass 1: build the @RG/@PG id maps in header-line order so the
        // gen_unique_id PRNG draw sequence matches upstream exactly.
        for line in raw.lines() {
            if !id_mode.attach_rg
                && line.starts_with("@RG\t")
                && let Some(id) = header_field(line, "ID:")
            {
                let new = assign(id, &mut seen_rg, &mut rng, id_mode.combine_rg);
                rg_map.insert(id.to_string(), new);
            } else if line.starts_with("@PG\t")
                && let Some(id) = header_field(line, "ID:")
            {
                let new = assign(id, &mut seen_pg, &mut rng, id_mode.combine_pg);
                pg_map.insert(id.to_string(), new);
            }
        }

        let mut emitted_file_rg = false;
        for line in raw.lines() {
            if let Some(rest) = line.strip_prefix("@HD") {
                if hd.is_none() {
                    hd = Some(format!("@HD{rest}"));
                }
            } else if line.starts_with("@SQ\t")
                && let Some(sn) = header_field(line, "SN:")
                && seen_sn.insert(sn.to_string())
            {
                sq.push(line.to_string());
            } else if line.starts_with("@RG\t") {
                if let Some(stem) = &stem {
                    // `-r`: emit the file's *first* @RG with ID→stem.
                    if !emitted_file_rg {
                        let mut m = HashMap::new();
                        if let Some(id) = header_field(line, "ID:") {
                            m.insert(id.to_string(), stem.clone());
                        }
                        rg_lines.push(rewrite_field(line, "ID:", &m));
                        emitted_file_rg = true;
                    }
                    continue;
                }
                let l = rewrite_field(line, "ID:", &rg_map);
                let final_id = header_field(line, "ID:")
                    .and_then(|i| rg_map.get(i).cloned())
                    .unwrap_or_default();
                if emitted_rg.insert(final_id) {
                    rg_lines.push(rewrite_field(&l, "PG:", &pg_map));
                }
            } else if line.starts_with("@PG\t") {
                let l = rewrite_field(line, "ID:", &pg_map);
                let final_id = header_field(line, "ID:")
                    .and_then(|i| pg_map.get(i).cloned())
                    .unwrap_or_default();
                if emitted_pg.insert(final_id) {
                    pg_lines.push(rewrite_field(&l, "PP:", &pg_map));
                }
            } else if line.starts_with("@CO\t") {
                co.push(line.to_string());
            }
        }

        // `-r` with no @RG in the file still gets one synthesized line.
        if let Some(stem) = &stem
            && !emitted_file_rg
        {
            rg_lines.push(format!("@RG\tID:{stem}"));
        }

        trans.rg.push(rg_map);
        trans.pg.push(pg_map);
        force_rg.push(stem);
    }

    let mut text = String::new();
    text.push_str(hd.as_deref().unwrap_or("@HD\tVN:1.6"));
    text.push('\n');
    for l in sq.iter().chain(&rg_lines).chain(&pg_lines).chain(&co) {
        text.push_str(l);
        text.push('\n');
    }
    Ok((text, trans, force_rg))
}

/// Remaps a record's `RG:Z:` / `PG:Z:` aux values through the trans maps.
fn remap_record_rg_pg(
    record: &mut RecordBuf,
    rg: &HashMap<String, String>,
    pg: &HashMap<String, String>,
    force_rg: Option<&str>,
    warned_missing_tags: &mut HashSet<([u8; 2], String)>,
) {
    use htslib_rs::sam::alignment::record::data::field::Tag;
    use htslib_rs::sam::alignment::record_buf::data::field::Value;
    // `-r`: force RG:Z: to the filename-derived id on every record.
    if let Some(stem) = force_rg {
        // `bam_translate` remaps PG before the forced RG is appended, so
        // records with both tags end with `PG:Z:<mapped> RG:Z:<stem>`.
        remap_aux_string_append_tail(record, [b'P', b'G'], pg, warned_missing_tags);
        let t = Tag::from([b'R', b'G']);
        // Upstream `-r` deletes any existing RG then appends the new one
        // at the tail (order-preserving, like bam_aux_del+bam_aux_append).
        let mut fields: Vec<_> = record
            .data()
            .iter()
            .filter(|(tg, _)| *tg != t)
            .map(|(tg, v)| (tg, v.clone()))
            .collect();
        fields.push((t, Value::String(bstr::BString::from(stem))));
        *record.data_mut() = fields.into_iter().collect();
        return;
    }
    for (tag, map) in [([b'R', b'G'], rg), ([b'P', b'G'], pg)] {
        remap_aux_string_append_tail(record, tag, map, warned_missing_tags);
    }
}

fn remap_aux_string_append_tail(
    record: &mut RecordBuf,
    tag: [u8; 2],
    map: &HashMap<String, String>,
    warned_missing_tags: &mut HashSet<([u8; 2], String)>,
) {
    use htslib_rs::sam::alignment::record::data::field::Tag;
    use htslib_rs::sam::alignment::record_buf::data::field::Value;

    let t = Tag::from(tag);
    let Some(Value::String(v)) = record.data().get(&t) else {
        return;
    };
    let cur = String::from_utf8_lossy(v).into_owned();

    let mut fields: Vec<_> = record
        .data()
        .iter()
        .filter(|(tg, _)| *tg != t)
        .map(|(tg, v)| (tg, v.clone()))
        .collect();

    if let Some(n) = map.get(&cur) {
        fields.push((t, Value::String(bstr::BString::from(n.clone()))));
    } else if warned_missing_tags.insert((tag, cur.clone())) {
        let tag_name = std::str::from_utf8(&tag).unwrap_or("??");
        let read_name = record
            .name()
            .map(|name| String::from_utf8_lossy(name).into_owned())
            .unwrap_or_else(|| "*".to_string());
        eprintln!(
            "[bam_translate] {tag_name} tag \"{cur}\" on read \"{read_name}\" encountered with no corresponding entry in header, tag lost. Unknown tags are only reported once per input file for each tag ID."
        );
    }
    // Otherwise upstream `bam_translate` drops a tag with no corresponding
    // header entry ("tag lost").
    *record.data_mut() = fields.into_iter().collect();
}

fn merge_comments(output_header: &mut sam::Header, input_header: &sam::Header) {
    output_header
        .comments_mut()
        .extend(input_header.comments().iter().cloned());
}

fn coordinate_key(r: &RecordBuf) -> (i32, i64) {
    let tid = r
        .reference_sequence_id()
        .map(|t| t as i32)
        .unwrap_or(i32::MAX);
    let pos = r.alignment_start().map(usize::from).unwrap_or(0) as i64;
    (tid, pos)
}

fn name_key(r: &RecordBuf) -> Vec<u8> {
    r.name().map(|s| s.to_vec()).unwrap_or_default()
}

fn compare_by_tag(a: &RecordBuf, b: &RecordBuf, tag: [u8; 2], name_sort: bool) -> Ordering {
    tag_sort_value(a, tag)
        .cmp(&tag_sort_value(b, tag))
        .then_with(|| {
            if name_sort {
                natural_name_cmp(a, b)
            } else {
                coordinate_key(a).cmp(&coordinate_key(b))
            }
        })
}

fn natural_name_cmp(a: &RecordBuf, b: &RecordBuf) -> Ordering {
    strnum_cmp(&name_key(a), &name_key(b))
}

fn strnum_cmp(a: &[u8], b: &[u8]) -> Ordering {
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        if a[i].is_ascii_digit() && b[j].is_ascii_digit() {
            let (si, sj) = (i, j);
            while i < a.len() && a[i].is_ascii_digit() {
                i += 1;
            }
            while j < b.len() && b[j].is_ascii_digit() {
                j += 1;
            }

            let na = trim_leading_zeroes(&a[si..i]);
            let nb = trim_leading_zeroes(&b[sj..j]);
            match na.len().cmp(&nb.len()).then_with(|| na.cmp(nb)) {
                Ordering::Equal => {}
                o => return o,
            }
        } else {
            match a[i].cmp(&b[j]) {
                Ordering::Equal => {
                    i += 1;
                    j += 1;
                }
                o => return o,
            }
        }
    }

    a.len().cmp(&b.len())
}

fn trim_leading_zeroes(s: &[u8]) -> &[u8] {
    let trimmed = s.iter().position(|&b| b != b'0').unwrap_or(s.len());
    &s[trimmed..]
}

#[derive(Clone, Debug)]
enum TagSortValue {
    Missing,
    Character(u8),
    Array(String),
    Text(Vec<u8>),
    Int(i64),
    Float(f32),
}

impl TagSortValue {
    fn rank(&self) -> u8 {
        match self {
            Self::Missing => 0,
            Self::Character(_) => b'A',
            Self::Array(_) => b'B',
            Self::Text(_) => b'H',
            Self::Int(_) => b'c',
            Self::Float(_) => b'f',
        }
    }
}

impl PartialEq for TagSortValue {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for TagSortValue {}

impl PartialOrd for TagSortValue {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TagSortValue {
    fn cmp(&self, other: &Self) -> Ordering {
        use TagSortValue::*;

        match (self, other) {
            (Missing, Missing) => Ordering::Equal,
            (Missing, _) => Ordering::Less,
            (_, Missing) => Ordering::Greater,
            (Int(a), Int(b)) => a.cmp(b),
            (Float(a), Float(b)) => a.total_cmp(b),
            (Int(a), Float(b)) => (*a as f32).total_cmp(b),
            (Float(a), Int(b)) => a.total_cmp(&(*b as f32)),
            (Character(a), Character(b)) => a.cmp(b),
            (Text(a), Text(b)) => a.cmp(b),
            (Array(a), Array(b)) => a.cmp(b),
            _ => self.rank().cmp(&other.rank()),
        }
    }
}

fn tag_sort_value(record: &RecordBuf, tag: [u8; 2]) -> TagSortValue {
    use sam::alignment::record_buf::data::field::Value;

    match record.data().get(&tag) {
        None => TagSortValue::Missing,
        Some(Value::Character(c)) => TagSortValue::Character(*c),
        Some(Value::Int8(n)) => TagSortValue::Int(i64::from(*n)),
        Some(Value::UInt8(n)) => TagSortValue::Int(i64::from(*n)),
        Some(Value::Int16(n)) => TagSortValue::Int(i64::from(*n)),
        Some(Value::UInt16(n)) => TagSortValue::Int(i64::from(*n)),
        Some(Value::Int32(n)) => TagSortValue::Int(i64::from(*n)),
        Some(Value::UInt32(n)) => TagSortValue::Int(i64::from(*n)),
        Some(Value::Float(n)) => TagSortValue::Float(*n),
        Some(Value::String(s)) | Some(Value::Hex(s)) => TagSortValue::Text(s.to_vec()),
        Some(Value::Array(array)) => TagSortValue::Array(format!("{:?}", array)),
    }
}

fn read_sam_records(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(input)?));
    let header = reader.read_header()?;
    let mut records = Vec::new();
    loop {
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        records.push(record);
    }
    Ok((header, records))
}

fn set_sort_order(header: &mut sam::Header, so: &str, go: Option<&str>, ss: Option<&str>) {
    use bstr::BString;
    use sam::header::record::value::map::header::tag::GROUP_ORDER;
    use sam::header::record::value::map::{self, Map};
    if let Some(hd) = header.header_mut() {
        hd.other_fields_mut()
            .insert(map::header::tag::SORT_ORDER, BString::from(so));
        match go {
            Some(go) => {
                hd.other_fields_mut().insert(GROUP_ORDER, BString::from(go));
            }
            None => {
                hd.other_fields_mut().shift_remove(&GROUP_ORDER);
            }
        }
        match ss {
            Some(ss) => {
                hd.other_fields_mut()
                    .insert(map::header::tag::SUBSORT_ORDER, BString::from(ss));
            }
            None => {
                hd.other_fields_mut()
                    .shift_remove(&map::header::tag::SUBSORT_ORDER);
            }
        }
    } else {
        let mut hd: Map<map::Header> = Map::default();
        hd.other_fields_mut()
            .insert(map::header::tag::SORT_ORDER, BString::from(so));
        if let Some(go) = go {
            hd.other_fields_mut().insert(GROUP_ORDER, BString::from(go));
        }
        if let Some(ss) = ss {
            hd.other_fields_mut()
                .insert(map::header::tag::SUBSORT_ORDER, BString::from(ss));
        }
        *header.header_mut() = Some(hd);
    }
}

trait MergeSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
    fn finish(self: Box<Self>) -> io::Result<()> {
        Ok(())
    }
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(File);
struct SamStdout(io::Stdout);
struct CramOut<W: Write> {
    writer: bam::io::Writer<bgzf::io::Writer<File>>,
    tmp_bam_path: crate::tmp_file::TempPath,
    reference: PathBuf,
    out: W,
}

impl MergeSink for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl MergeSink for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl MergeSink for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        // Shared renderer: htslib `%g` float aux spelling.
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}
impl MergeSink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}
impl<W: Write> MergeSink for CramOut<W> {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.writer.write_alignment_record(header, record)
    }

    fn finish(self: Box<Self>) -> io::Result<()> {
        let this = *self;
        let CramOut {
            writer,
            tmp_bam_path,
            reference,
            out,
        } = this;
        drop(writer);

        let result = htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference(
            tmp_bam_path.path(),
            &reference,
            out,
        )
        .map(|_| ());

        tmp_bam_path.close().ok();
        result
    }
}

fn open_output(
    out: Option<&Path>,
    fmt: OutFmt,
    header: &sam::Header,
) -> io::Result<Box<dyn MergeSink>> {
    match (out, fmt) {
        (Some(p), OutFmt::Sam) => {
            let mut file = File::create(p)?;
            crate::sam_render::write_header(&mut file, header)?;
            Ok(Box::new(SamFile(file)))
        }
        (Some(p), OutFmt::Bam) => {
            let file = File::create(p)?;
            let mut writer = bam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(BamFile(writer)))
        }
        (Some(p), OutFmt::Cram) => {
            let out = File::create(p)?;
            open_cram_output(header, out)
        }
        (None, OutFmt::Sam) => {
            let mut stdout = io::stdout();
            crate::sam_render::write_header(&mut stdout, header)?;
            Ok(Box::new(SamStdout(stdout)))
        }
        (None, OutFmt::Bam) => {
            let mut writer = bam::io::Writer::new(io::stdout());
            writer.write_header(header)?;
            Ok(Box::new(BamStdout(writer)))
        }
        (None, OutFmt::Cram) => open_cram_output(header, io::stdout()),
    }
}

fn open_cram_output<W>(header: &sam::Header, out: W) -> io::Result<Box<dyn MergeSink>>
where
    W: Write + 'static,
{
    let reference = current_global_args().reference.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM output requires top-level --reference FILE",
        )
    })?;
    crate::reference::ensure_fai_index(&reference, None)?;

    let (tmp_bam_file, tmp_bam_path) = crate::tmp_file::create_temp_file("merge", Some("bam"))?;
    let mut writer = bam::io::Writer::new(tmp_bam_file);
    writer.write_header(header)?;

    Ok(Box::new(CramOut {
        writer,
        tmp_bam_path,
        reference,
        out,
    }))
}

fn write_bam_index(path: &Path) -> io::Result<()> {
    let index = htslib_rs::index_compat::build_bai(path)?;
    htslib_rs::index_compat::write_bai(append_extension(path, "bai"), &index)
}

fn append_extension(path: &Path, ext: &str) -> PathBuf {
    let mut s = path.as_os_str().to_os_string();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(
        w,
        "Usage: samtools merge [options] <out.bam> <in1.bam|in1.sam> [<in2.bam|in2.sam> ...]"
    )?;
    writeln!(w, "Options:")?;
    writeln!(w, "  -n              name sort")?;
    writeln!(w, "  -f              force overwrite output")?;
    writeln!(w, "  -R REGION       restrict indexed BAM inputs to REGION")?;
    writeln!(
        w,
        "  -L BED          restrict indexed BAM inputs to BED intervals"
    )?;
    writeln!(w, "  -o FILE         output to FILE")?;
    writeln!(w, "  --output-fmt sam|bam|cram")?;
    writeln!(w, "  --write-index   write BAI index for BAM file output")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::BString;
    use htslib_rs::sam::alignment::record::data::field::Tag;
    use htslib_rs::sam::alignment::record_buf::data::field::Value;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn tmp_dir(name: &str) -> PathBuf {
        static NEXT_ID: AtomicUsize = AtomicUsize::new(0);

        let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!(
            "samtools-rs-merge-{}-{}-{}",
            name,
            std::process::id(),
            id
        ));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_sam_header(dir: &Path, name: &str, text: &str) -> PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, text).unwrap();
        path
    }

    fn record_with_string_tags(name: &str, tags: &[([u8; 2], &str)]) -> RecordBuf {
        let mut record = RecordBuf::default();
        *record.name_mut() = Some(BString::from(name));
        for (tag, value) in tags {
            record.data_mut().insert(
                Tag::from(*tag),
                Value::String(BString::from(value.as_bytes())),
            );
        }
        record
    }

    fn aux_string(record: &RecordBuf, tag: [u8; 2]) -> Option<String> {
        match record.data().get(&Tag::from(tag)) {
            Some(Value::String(s)) => Some(String::from_utf8_lossy(s).into_owned()),
            _ => None,
        }
    }

    #[test]
    fn reconcile_merge_headers_unions_sq_and_preserves_input_header_text() {
        let dir = tmp_dir("header-union");
        let init_text = "@HD\tVN:1.4\tSO:unknown\n@SQ\tSN:fish\tLN:133\tSP:frog\n";
        let translate_text = "@HD\tVN:1.4\tSO:unknown\n@SQ\tSN:donkey\tLN:133\n@SQ\tSN:fish\tLN:133\n@RG\tID:fish\tPU:trans\n";
        let first = write_sam_header(&dir, "first.sam", init_text);
        let second = write_sam_header(&dir, "second.sam", translate_text);

        let (merged, trans, force_rg) = reconcile_merge_headers(
            &[first.clone(), second.clone()],
            0x1234330e,
            MergeIdMode::default(),
        )
        .unwrap();

        assert_eq!(
            merged,
            concat!(
                "@HD\tVN:1.4\tSO:unknown\n",
                "@SQ\tSN:fish\tLN:133\tSP:frog\n",
                "@SQ\tSN:donkey\tLN:133\n",
                "@RG\tID:fish\tPU:trans\n",
            )
        );
        assert!(trans.rg[0].is_empty());
        assert_eq!(trans.rg[1].get("fish").map(String::as_str), Some("fish"));
        assert!(trans.pg.iter().all(HashMap::is_empty));
        assert_eq!(force_rg, vec![None, None]);
        assert_eq!(std::fs::read_to_string(&first).unwrap(), init_text);
        assert_eq!(std::fs::read_to_string(&second).unwrap(), translate_text);
    }

    #[test]
    fn reconcile_merge_headers_remaps_colliding_rg_pg_header_fields() {
        let dir = tmp_dir("rg-pg-collision");
        let first_text = concat!(
            "@HD\tVN:1.4\tSO:unknown\n",
            "@SQ\tSN:fish\tLN:133\tSP:frog\n",
            "@RG\tID:fish\tPU:out\n",
            "@PG\tXX:dummyx\tID:fish\tDS:out\n",
            "@PG\tPP:fish\tID:hook\tDS:out\n",
        );
        let second_text = concat!(
            "@HD\tVN:1.4\tSO:unknown\n",
            "@SQ\tSN:donkey\tLN:133\n",
            "@SQ\tSN:fish\tLN:133\n",
            "@RG\tID:fish\tPU:trans\n",
            "@PG\tXX:dummy\tID:fish\tDS:trans\n",
            "@PG\tPP:fish\tID:hook\tDS:trans\n",
        );
        let first = write_sam_header(&dir, "first.sam", first_text);
        let second = write_sam_header(&dir, "second.sam", second_text);

        let (merged, trans, _) =
            reconcile_merge_headers(&[first, second], 0x1234330e, MergeIdMode::default()).unwrap();

        let second_rg = trans.rg[1].get("fish").unwrap();
        let second_pg_fish = trans.pg[1].get("fish").unwrap();
        let second_pg_hook = trans.pg[1].get("hook").unwrap();

        assert_eq!(trans.rg[0].get("fish").map(String::as_str), Some("fish"));
        assert!(second_rg.starts_with("fish-"));
        assert_eq!(trans.pg[0].get("fish").map(String::as_str), Some("fish"));
        assert_eq!(trans.pg[0].get("hook").map(String::as_str), Some("hook"));
        assert!(second_pg_fish.starts_with("fish-"));
        assert!(second_pg_hook.starts_with("hook-"));

        assert!(merged.contains("@RG\tID:fish\tPU:out\n"));
        assert!(merged.contains(&format!("@RG\tID:{second_rg}\tPU:trans\n")));
        assert!(merged.contains("@PG\tXX:dummyx\tID:fish\tDS:out\n"));
        assert!(merged.contains("@PG\tPP:fish\tID:hook\tDS:out\n"));
        assert!(merged.contains(&format!("@PG\tXX:dummy\tID:{second_pg_fish}\tDS:trans\n")));
        assert!(merged.contains(&format!(
            "@PG\tPP:{second_pg_fish}\tID:{second_pg_hook}\tDS:trans\n"
        )));
    }

    #[test]
    fn reconcile_merge_headers_attach_rg_uses_filename_stem() {
        let dir = tmp_dir("attach-rg");
        let input_text = concat!(
            "@HD\tVN:1.4\tSO:unknown\n",
            "@SQ\tSN:fish\tLN:133\n",
            "@RG\tID:old\tSM:sample\n",
        );
        let input = write_sam_header(&dir, "filename.sam", input_text);

        let (merged, trans, force_rg) = reconcile_merge_headers(
            &[input],
            0x1234330e,
            MergeIdMode {
                attach_rg: true,
                ..MergeIdMode::default()
            },
        )
        .unwrap();

        assert!(trans.rg[0].is_empty());
        assert_eq!(force_rg, vec![Some("filename".to_string())]);
        assert!(merged.contains("@RG\tID:filename\tSM:sample\n"));
    }

    #[test]
    fn remap_record_rg_pg_translates_known_rg_pg_tags() {
        let mut record = record_with_string_tags(
            "123456789",
            &[([b'R', b'G'], "hello"), ([b'P', b'G'], "quail")],
        );
        let rg = HashMap::from([("hello".to_string(), "goodbye".to_string())]);
        let pg = HashMap::from([("quail".to_string(), "bird".to_string())]);
        let mut warned = HashSet::new();

        remap_record_rg_pg(&mut record, &rg, &pg, None, &mut warned);

        assert_eq!(
            aux_string(&record, [b'R', b'G']).as_deref(),
            Some("goodbye")
        );
        assert_eq!(aux_string(&record, [b'P', b'G']).as_deref(), Some("bird"));
        assert!(warned.is_empty());
    }

    #[test]
    fn remap_record_rg_pg_drops_unknown_tags_once() {
        let rg = HashMap::new();
        let pg = HashMap::new();
        let mut warned = HashSet::new();
        let mut record = record_with_string_tags(
            "123456789",
            &[([b'R', b'G'], "rg4hello"), ([b'P', b'G'], "pg5hello")],
        );

        remap_record_rg_pg(&mut record, &rg, &pg, None, &mut warned);

        assert_eq!(aux_string(&record, [b'R', b'G']), None);
        assert_eq!(aux_string(&record, [b'P', b'G']), None);
        assert!(warned.contains(&([b'R', b'G'], "rg4hello".to_string())));
        assert!(warned.contains(&([b'P', b'G'], "pg5hello".to_string())));

        let mut second = record_with_string_tags(
            "987654321",
            &[([b'R', b'G'], "rg4hello"), ([b'P', b'G'], "pg5hello")],
        );
        remap_record_rg_pg(&mut second, &rg, &pg, None, &mut warned);

        assert_eq!(aux_string(&second, [b'R', b'G']), None);
        assert_eq!(aux_string(&second, [b'P', b'G']), None);
        assert_eq!(warned.len(), 2);
    }

    #[test]
    fn remap_record_rg_pg_force_rg_appends_filename_rg() {
        let mut record = record_with_string_tags(
            "123456789",
            &[([b'R', b'G'], "old"), ([b'P', b'G'], "quail")],
        );
        let rg = HashMap::new();
        let pg = HashMap::from([("quail".to_string(), "bird".to_string())]);
        let mut warned = HashSet::new();

        remap_record_rg_pg(&mut record, &rg, &pg, Some("file"), &mut warned);

        assert_eq!(aux_string(&record, [b'P', b'G']).as_deref(), Some("bird"));
        assert_eq!(aux_string(&record, [b'R', b'G']).as_deref(), Some("file"));
        assert!(warned.is_empty());
    }
}
