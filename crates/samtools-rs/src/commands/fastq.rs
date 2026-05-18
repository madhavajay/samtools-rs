//! `samtools fastq` / `samtools fasta` / `samtools bam2fq` — convert
//! SAM/BAM records to FASTQ or FASTA text.
//!
//! Mirrors `main_bam2fq` in `bam_fastq.c`. Supports single-output mode (all
//! reads written to stdout, `-o FILE`, or `-0 FILE`), paired-output split
//! (`-1`/`-2`/`-s`/`-0`) with upstream-style name-grouped routing where
//! adjacent records sharing a qname pick the best per readpart and flush as
//! a unit (paired R1+R2 to `-1`/`-2`, R1-only or R2-only singletons to `-s`
//! when set or falling back to `-1`/`-2`, and READ_OTHER to `-0`).
//!
//! Also supports `--i1`/`--i2` index FASTQ extraction with
//! `--index-format` (default `i*i*`), `--quality-tag`, and
//! `--barcode-tag`, emitting one index record per adjacent qname-group
//! (upstream `flush_rec` → `output_index`) with the htslib-exact CASAVA
//! barcode normalization and, under `-i`, the CASAVA comment.
//!
//! The group's barcode is propagated across mates so an R2 (or other)
//! record lacking its own `BC` inherits the R1 mate's barcode in its
//! CASAVA comment (upstream `bam_fastq.c:952`).
//!
use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use htslib_rs::{bam, format::Exact, sam};

use crate::aux_list::parse_aux_list;
use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
use crate::io as sam_io;
use crate::sam_render::format_aux_float;

#[derive(Clone, Debug, Eq, PartialEq)]
enum AuxSelection {
    None,
    Tags(Vec<[u8; 2]>),
    All,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct TagFilter {
    tag: [u8; 2],
    values: Option<HashSet<String>>,
}

#[derive(Clone, Copy)]
struct FastqRenderOptions<'a> {
    append_read_number: bool,
    use_original_quality: bool,
    default_quality: Option<u8>,
    umi_tags: Option<&'a [[u8; 2]]>,
    casava: bool,
    barcode_tag: [u8; 2],
    aux_selection: &'a AuxSelection,
    strip_soft_clips: bool,
    soft_clip_backup_tag: Option<[u8; 2]>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct FlagFilters {
    require: u16,
    include_any: u16,
    exclude: u16,
    exclude_all: u16,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct FastqStats {
    processed_reads: usize,
    discarded_singletons: usize,
}

impl FlagFilters {
    fn is_enabled(self) -> bool {
        self.require != 0 || self.include_any != 0 || self.exclude != 0 || self.exclude_all != 0
    }
}

impl AuxSelection {
    fn is_enabled(&self) -> bool {
        !matches!(self, Self::None)
    }
}

/// Merge a new `-d` / `-D` filter into the existing list. Upstream rejects
/// repeated invocations with different tag names; same-tag invocations
/// union their value sets so `-d NM:13 -d NM:14` keeps records with
/// `NM:i:{13,14}`. A presence-only filter (`-d TAG`) "wins" if combined
/// with any value filter on the same tag (records pass on presence alone).
fn merge_tag_filter(filters: &mut Vec<TagFilter>, new_filter: TagFilter) -> Result<(), String> {
    if let Some(existing) = filters.iter_mut().find(|f| f.tag == new_filter.tag) {
        match (&mut existing.values, new_filter.values) {
            (Some(existing_set), Some(new_set)) => existing_set.extend(new_set),
            (Some(_), None) => existing.values = None,
            (None, _) => {}
        }
        return Ok(());
    }
    if let Some(other) = filters.iter().next()
        && other.tag != new_filter.tag
    {
        return Err(format!(
            "different tag \"{}{}\" specified after \"{}{}\"",
            char::from(new_filter.tag[0]),
            char::from(new_filter.tag[1]),
            char::from(other.tag[0]),
            char::from(other.tag[1]),
        ));
    }
    filters.push(new_filter);
    Ok(())
}

/// Union-merges `extra` into `selection`. `None` becomes `Tags(extra)`;
/// `All` is unchanged; `Tags(existing)` extends with non-duplicate tags
/// from `extra`. Matches upstream's accumulating `-t` / `-T` behavior.
fn merge_aux_selection(selection: &mut AuxSelection, extra: &[[u8; 2]]) {
    match selection {
        AuxSelection::None => {
            *selection = AuxSelection::Tags(extra.to_vec());
        }
        AuxSelection::All => {}
        AuxSelection::Tags(existing) => {
            for tag in extra {
                if !existing.iter().any(|t| t == tag) {
                    existing.push(*tag);
                }
            }
        }
    }
}

fn apply_aux_selection_arg(raw: &str, selection: &mut AuxSelection) -> Result<(), String> {
    match raw {
        "" | "*" => {
            *selection = AuxSelection::All;
            Ok(())
        }
        _ => {
            let tags = parse_aux_list(raw)
                .map_err(|e| format!("invalid -T value \"{raw}\": {e}"))?
                .into_iter()
                .collect::<Vec<_>>();
            merge_aux_selection(selection, &tags);
            Ok(())
        }
    }
}

/// Entry point for `samtools fastq` / `samtools fasta` / `samtools bam2fq`.
pub fn main(args: &[OsString]) -> ExitCode {
    let sub_name = args.first().and_then(|a| a.to_str()).unwrap_or("fastq");
    let fasta_mode = sub_name == "fasta";

    let mut output: Option<PathBuf> = None;
    let mut other_output: Option<PathBuf> = None;
    let mut read1_output: Option<PathBuf> = None;
    let mut read2_output: Option<PathBuf> = None;
    let mut singleton_output: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let mut require_flags = 0u16;
    let mut include_any_flags = 0u16;
    let mut exclude_flags = 0x900u16;
    let mut exclude_all_flags = 0u16;
    let mut aux_selection = AuxSelection::None;
    let mut tag_filters = Vec::new();
    let mut append_read_number_override: Option<bool> = None;
    let mut use_original_quality = false;
    let mut default_quality: Option<u8> = None;
    let mut umi_enabled = false;
    let mut umi_tags = vec![*b"OX", *b"RX"];
    let mut casava = false;
    let mut barcode_tag = *b"BC";
    let mut quality_tag = *b"QT";
    let mut index_file_1: Option<PathBuf> = None;
    let mut index_file_2: Option<PathBuf> = None;
    let mut index_format_arg: Option<String> = None;
    let mut strip_soft_clips = false;
    let mut soft_clip_backup_tag = Some(*b"s0");
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-o" => {
                output = iter.next().map(PathBuf::from);
            }
            "-0" => {
                other_output = iter.next().map(PathBuf::from);
            }
            "-1" => {
                read1_output = iter.next().map(PathBuf::from);
            }
            "-2" => {
                read2_output = iter.next().map(PathBuf::from);
            }
            "-s" => {
                singleton_output = iter.next().map(PathBuf::from);
            }
            "-f" | "--require-flags" => {
                require_flags = match parse_flag_arg(iter.next(), "-f", sub_name) {
                    Ok(flag) => flag,
                    Err(code) => return code,
                };
            }
            "--rf" | "--incl-flags" | "--include-flags" => {
                include_any_flags |= match parse_flag_arg(iter.next(), "--rf", sub_name) {
                    Ok(flag) => flag,
                    Err(code) => return code,
                };
            }
            "-F" | "--excl-flags" | "--exclude-flags" => {
                exclude_flags = match parse_flag_arg(iter.next(), "-F", sub_name) {
                    Ok(flag) => flag,
                    Err(code) => return code,
                };
            }
            "-G" => {
                exclude_all_flags = match parse_flag_arg(iter.next(), "-G", sub_name) {
                    Ok(flag) => flag,
                    Err(code) => return code,
                };
            }
            "-T" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error(sub_name, "missing value for -T");
                    return ExitCode::from(1);
                };
                if let Err(e) = apply_aux_selection_arg(raw, &mut aux_selection) {
                    print_error(sub_name, e);
                    return ExitCode::from(1);
                }
            }
            _ if s.starts_with("-T") && s.len() > 2 => {
                let raw = &s[2..];
                if let Err(e) = apply_aux_selection_arg(raw, &mut aux_selection) {
                    print_error(sub_name, e);
                    return ExitCode::from(1);
                }
            }
            "-t" => {
                merge_aux_selection(&mut aux_selection, &[*b"RG", *b"BC", *b"QT"]);
            }
            "-d" | "--tag" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error(sub_name, "missing value for -d");
                    return ExitCode::from(1);
                };
                match parse_tag_filter(raw) {
                    Ok(filter) => {
                        if let Err(e) = merge_tag_filter(&mut tag_filters, filter) {
                            print_error(sub_name, e);
                            return ExitCode::from(1);
                        }
                    }
                    Err(e) => {
                        print_error(sub_name, e);
                        return ExitCode::from(1);
                    }
                }
            }
            "-D" | "--tag-file" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error(sub_name, "missing value for -D");
                    return ExitCode::from(1);
                };
                match parse_tag_filter_file(raw) {
                    Ok(filter) => {
                        if let Err(e) = merge_tag_filter(&mut tag_filters, filter) {
                            print_error(sub_name, e);
                            return ExitCode::from(1);
                        }
                    }
                    Err(e) => {
                        print_error(sub_name, e);
                        return ExitCode::from(1);
                    }
                }
            }
            "-n" => {
                append_read_number_override = Some(false);
            }
            "-N" => {
                append_read_number_override = Some(true);
            }
            "-O" => {
                use_original_quality = true;
            }
            "-i" => {
                casava = true;
            }
            "-U" | "--UMI" | "--umi" => {
                umi_enabled = true;
            }
            "--UMI-tag" | "--umi-tag" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error(sub_name, "missing value for --UMI-tag");
                    return ExitCode::from(1);
                };
                umi_tags = match parse_ordered_tag_list(raw) {
                    Ok(tags) => tags,
                    Err(e) => {
                        print_error(sub_name, format!("invalid --UMI-tag value \"{raw}\": {e}"));
                        return ExitCode::from(1);
                    }
                };
            }
            "--barcode-tag" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error(sub_name, "missing value for --barcode-tag");
                    return ExitCode::from(1);
                };
                barcode_tag = match parse_filter_tag(raw) {
                    Ok(tag) => tag,
                    Err(e) => {
                        print_error(sub_name, e);
                        return ExitCode::from(1);
                    }
                };
            }
            "--quality-tag" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error(sub_name, "missing value for --quality-tag");
                    return ExitCode::from(1);
                };
                quality_tag = match parse_filter_tag(raw) {
                    Ok(tag) => tag,
                    Err(e) => {
                        print_error(sub_name, e);
                        return ExitCode::from(1);
                    }
                };
            }
            "--i1" | "--I1" => {
                index_file_1 = iter.next().map(PathBuf::from);
            }
            "--i2" | "--I2" => {
                index_file_2 = iter.next().map(PathBuf::from);
            }
            "--index-format" | "--if" | "--IF" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error(sub_name, "missing value for --index-format");
                    return ExitCode::from(1);
                };
                index_format_arg = Some(raw.to_string());
            }
            "--no-sc" => {
                strip_soft_clips = true;
            }
            "--no-sc-bkp" => {
                soft_clip_backup_tag = None;
            }
            "--sc-aux" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error(sub_name, "missing value for --sc-aux");
                    return ExitCode::from(1);
                };
                soft_clip_backup_tag = match parse_filter_tag(raw) {
                    Ok(tag) => Some(tag),
                    Err(e) => {
                        print_error(sub_name, e);
                        return ExitCode::from(1);
                    }
                };
            }
            "-v" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error(sub_name, "missing value for -v");
                    return ExitCode::from(1);
                };
                default_quality = match raw.parse::<u8>() {
                    Ok(q) if q <= 93 => Some(q),
                    _ => {
                        print_error(sub_name, format!("invalid -v value \"{}\"", raw));
                        return ExitCode::from(1);
                    }
                };
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "--help" => {
                let _ = print_usage(sub_name);
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error(
                    sub_name,
                    format!(
                        "option `{}` is not yet supported in samtools-rs {}",
                        s, sub_name
                    ),
                );
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                }
            }
        }
    }
    let umi_tags = umi_enabled.then_some(umi_tags);

    if index_file_2.is_some() && index_file_1.is_none() {
        print_error(sub_name, "Index one specified, but index two not given");
        return ExitCode::from(1);
    }

    let parsed_index_format = match index_format_arg.as_deref() {
        Some(spec) => match parse_index_format(spec) {
            Ok(items) => Some(items),
            Err(e) => {
                print_error(sub_name, e);
                return ExitCode::from(1);
            }
        },
        None => None,
    };

    let n_index_segments = parsed_index_format
        .as_ref()
        .map(|items| items.iter().filter(|item| item.is_index).count())
        .unwrap_or(0);
    if n_index_segments > 2 {
        print_error(sub_name, "Invalid index format: more than 2 indexes");
        return ExitCode::from(1);
    }
    if index_file_1.is_some() && index_format_arg.is_some() && n_index_segments == 0 {
        print_error(sub_name, "index_format not specified, but index file given");
        return ExitCode::from(1);
    }

    let effective_index_format = match parsed_index_format {
        Some(items) => items,
        None => parse_index_format("i*i*").expect("default index format parses"),
    };

    let stdin_input = input.as_ref().is_none_or(|path| path.as_os_str() == "-");

    let format = if stdin_input {
        None
    } else {
        let input = input.as_ref().expect("non-stdin input exists");
        if !input.exists() {
            print_hts_open_missing(input);
            print_error(
                "bam2fq",
                format!(
                    "Cannot read file \"{}\": No such file or directory",
                    input.display()
                ),
            );
            return ExitCode::from(1);
        }
        match sam_io::sam_open_format(input) {
            Ok(f) => Some(f),
            Err(e) => {
                print_error(sub_name, e.to_string());
                return ExitCode::from(1);
            }
        }
    };
    let input_exact = effective_fastq_input_exact(input.as_deref(), format.as_ref());

    let flag_filters = FlagFilters {
        require: require_flags,
        include_any: include_any_flags,
        exclude: exclude_flags,
        exclude_all: exclude_all_flags,
    };
    let filtering = flag_filters.is_enabled();
    let split_mode = read1_output.is_some()
        || read2_output.is_some()
        || singleton_output.is_some()
        || other_output.is_some();
    let singleton_only = singleton_output.is_some()
        && read1_output.is_none()
        && read2_output.is_none()
        && other_output.is_none();
    let other_only = other_output.is_some()
        && read1_output.is_none()
        && read2_output.is_none()
        && singleton_output.is_none();
    let append_read_number =
        append_read_number_override.unwrap_or(!split_mode || singleton_only || other_only);
    let append_index_read_number =
        append_read_number_override.unwrap_or(!split_mode || singleton_only);
    let render_options = FastqRenderOptions {
        append_read_number,
        use_original_quality,
        default_quality,
        umi_tags: umi_tags.as_deref(),
        casava,
        barcode_tag,
        aux_selection: &aux_selection,
        strip_soft_clips,
        soft_clip_backup_tag: strip_soft_clips.then_some(soft_clip_backup_tag).flatten(),
    };

    if split_mode {
        if fasta_mode && !tag_filters.is_empty() {
            print_error(
                sub_name,
                "-d/-D tag filtering is currently supported for FASTQ single-output mode only",
            );
            return ExitCode::from(1);
        }

        let singleton_set = singleton_output.is_some();
        let split = if stdin_input {
            let stdin = io::stdin().lock();
            let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(stdin));
            if fasta_mode {
                view_sam_reader_as_fasta_split(
                    &mut reader,
                    flag_filters,
                    append_read_number,
                    umi_tags.as_deref(),
                    casava,
                    barcode_tag,
                    singleton_set,
                )
            } else {
                view_sam_reader_as_fastq_split_with_aux(
                    &mut reader,
                    flag_filters,
                    render_options,
                    &tag_filters,
                    singleton_set,
                )
            }
        } else {
            let input = input.as_ref().expect("non-stdin input exists");
            match (input_exact.expect("non-stdin format exists"), fasta_mode) {
                (Exact::Sam, false) => view_sam_path_as_fastq_split(
                    input,
                    flag_filters,
                    render_options,
                    &tag_filters,
                    singleton_set,
                ),
                (Exact::Sam, true) => view_sam_path_as_fasta_split(
                    input,
                    flag_filters,
                    append_read_number,
                    umi_tags.as_deref(),
                    casava,
                    barcode_tag,
                    singleton_set,
                ),
                (Exact::Bam, false) => view_bam_path_as_fastq_split_with_aux(
                    input,
                    flag_filters,
                    render_options,
                    &tag_filters,
                    singleton_set,
                ),
                (Exact::Bam, true) => view_bam_path_as_fasta_split(
                    input,
                    flag_filters,
                    append_read_number,
                    umi_tags.as_deref(),
                    casava,
                    barcode_tag,
                    singleton_set,
                ),
                (Exact::Cram, false) => view_cram_path_as_fastq_split_with_aux(
                    input,
                    flag_filters,
                    render_options,
                    &tag_filters,
                    singleton_set,
                ),
                (Exact::Cram, true) => view_cram_path_as_fasta_split(
                    input,
                    flag_filters,
                    append_read_number,
                    umi_tags.as_deref(),
                    casava,
                    barcode_tag,
                    singleton_set,
                ),
                _ => {
                    print_error(
                        sub_name,
                        "only SAM, BAM, and CRAM input are currently supported",
                    );
                    return ExitCode::from(1);
                }
            }
        };

        let split = match split {
            Ok(split) => split,
            Err(e) => {
                print_error_errno(
                    sub_name,
                    format!(
                        "conversion failed for \"{}\"",
                        input
                            .as_ref()
                            .map(|path| path.display().to_string())
                            .unwrap_or_else(|| "-".to_string())
                    ),
                    &e,
                );
                return ExitCode::from(1);
            }
        };

        let paired_share_path = match (read1_output.as_ref(), read2_output.as_ref()) {
            (Some(a), Some(b)) => a == b,
            _ => false,
        };
        if paired_share_path && let Some(path) = read1_output.as_ref() {
            let interleaved = interleave_paired_records(&split.read1, &split.read2);
            if let Err(e) = write_text_file(path, &interleaved) {
                print_error_errno(sub_name, format!("open/write {}", path.display()), &e);
                return ExitCode::from(1);
            }
        } else {
            if let Some(path) = read1_output.as_ref()
                && let Err(e) = write_text_file(path, &split.read1)
            {
                print_error_errno(sub_name, format!("open/write {}", path.display()), &e);
                return ExitCode::from(1);
            }
            if let Some(path) = read2_output.as_ref()
                && let Err(e) = write_text_file(path, &split.read2)
            {
                print_error_errno(sub_name, format!("open/write {}", path.display()), &e);
                return ExitCode::from(1);
            }
        }
        if read1_output.is_none() && read2_output.is_none() && !singleton_only {
            let interleaved = interleave_paired_records(&split.read1, &split.read2);
            let mut out = match sam_io::open_text_output(None) {
                Ok(out) => out,
                Err(e) => {
                    print_error_errno(sub_name, "open stdout", &e);
                    return ExitCode::from(1);
                }
            };
            if let Err(e) = out.write_all(&interleaved)
                && e.kind() != io::ErrorKind::BrokenPipe
            {
                print_error_errno(sub_name, "write output", &e);
                return ExitCode::from(1);
            }
            if let Err(e) = sam_io::check_sam_close(&mut out)
                && e.kind() != io::ErrorKind::BrokenPipe
            {
                print_error_errno(sub_name, "close output", &e);
                return ExitCode::from(1);
            }
        }
        if let Some(path) = singleton_output.as_ref() {
            let payload = if singleton_only {
                let mut all = Vec::new();
                all.extend_from_slice(&split.read1);
                all.extend_from_slice(&split.read2);
                all.extend_from_slice(&split.singleton);
                all.extend_from_slice(&split.other);
                std::borrow::Cow::Owned(all)
            } else if other_output.is_some() {
                std::borrow::Cow::Borrowed(split.singleton.as_slice())
            } else {
                let mut merged = Vec::with_capacity(split.singleton.len() + split.other.len());
                merged.extend_from_slice(&split.singleton);
                merged.extend_from_slice(&split.other);
                std::borrow::Cow::Owned(merged)
            };
            if let Err(e) = write_text_file(path, payload.as_ref()) {
                print_error_errno(sub_name, format!("open/write {}", path.display()), &e);
                return ExitCode::from(1);
            }
        }
        if let Some(path) = other_output.as_ref() {
            let payload = if singleton_output.is_some() {
                std::borrow::Cow::Borrowed(split.other.as_slice())
            } else {
                let mut merged = Vec::with_capacity(split.singleton.len() + split.other.len());
                merged.extend_from_slice(&split.singleton);
                merged.extend_from_slice(&split.other);
                std::borrow::Cow::Owned(merged)
            };
            if let Err(e) = write_text_file(path, payload.as_ref()) {
                print_error_errno(sub_name, format!("open/write {}", path.display()), &e);
                return ExitCode::from(1);
            }
        }

        if (index_file_1.is_some() || index_file_2.is_some())
            && let Err(e) = emit_index_files(
                input.as_deref(),
                input_exact,
                stdin_input,
                flag_filters,
                IndexEmitOptions {
                    append_read_number: append_index_read_number,
                    use_original_quality,
                    default_quality,
                    umi_tags: umi_tags.as_deref(),
                    casava,
                    barcode_tag,
                    quality_tag,
                    index_format: &effective_index_format,
                    index_file_1: index_file_1.as_deref(),
                    index_file_2: index_file_2.as_deref(),
                    fasta_mode,
                },
            )
        {
            print_error_errno(sub_name, "index FASTQ output", &e);
            return ExitCode::from(1);
        }

        print_bam2fq_summary(FastqStats {
            processed_reads: count_split_records(&split, fasta_mode),
            discarded_singletons: count_fastx_records(&split.singleton, fasta_mode),
        });
        return ExitCode::SUCCESS;
    }

    let text = if stdin_input {
        if !tag_filters.is_empty() && fasta_mode {
            print_error(
                sub_name,
                "-d/-D tag filtering is currently supported for FASTQ single-output mode only",
            );
            return ExitCode::from(1);
        }
        let stdin = io::stdin().lock();
        let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(stdin));
        if !fasta_mode {
            view_sam_reader_as_fastq_text_with_aux(
                &mut reader,
                flag_filters,
                render_options,
                &tag_filters,
            )
        } else if fasta_mode {
            view_sam_reader_as_fasta_text(
                &mut reader,
                flag_filters,
                append_read_number,
                umi_tags.as_deref(),
                casava,
                barcode_tag,
            )
        } else {
            htslib_rs::alignment_compat::view_sam_as_fastq_text_from_reader_with_flag_filter_and_suffix(
                &mut reader,
                require_flags,
                exclude_flags,
                exclude_all_flags,
                append_read_number,
            )
        }
    } else {
        let input = input.as_ref().expect("non-stdin input exists");
        match (
            input_exact.expect("non-stdin format exists"),
            fasta_mode,
            filtering,
            use_original_quality
                || default_quality.is_some()
                || umi_tags.is_some()
                || casava
                || aux_selection.is_enabled()
                || !tag_filters.is_empty()
                || (!fasta_mode && flag_filters.include_any != 0),
        ) {
            (Exact::Sam, false, _, _) => view_sam_path_as_fastq_text_with_aux(
                input,
                flag_filters,
                render_options,
                &tag_filters,
            ),
            (Exact::Sam, true, _, true)
                if (umi_tags.is_some() || casava)
                    && !use_original_quality
                    && default_quality.is_none()
                    && !aux_selection.is_enabled()
                    && tag_filters.is_empty() =>
            {
                view_sam_path_as_fasta_text(
                    input,
                    flag_filters,
                    append_read_number,
                    umi_tags.as_deref(),
                    casava,
                    barcode_tag,
                )
            }
            (Exact::Bam, true, _, true)
                if (umi_tags.is_some() || casava)
                    && !use_original_quality
                    && default_quality.is_none()
                    && !aux_selection.is_enabled()
                    && tag_filters.is_empty() =>
            {
                view_bam_path_as_fasta_text(
                    input,
                    flag_filters,
                    append_read_number,
                    umi_tags.as_deref(),
                    casava,
                    barcode_tag,
                )
            }
            (Exact::Cram, true, _, true)
                if (umi_tags.is_some() || casava)
                    && !use_original_quality
                    && default_quality.is_none()
                    && !aux_selection.is_enabled()
                    && tag_filters.is_empty() =>
            {
                view_cram_path_as_fasta_text(
                    input,
                    flag_filters,
                    append_read_number,
                    umi_tags.as_deref(),
                    casava,
                    barcode_tag,
                )
            }
            (Exact::Sam, true, _, true)
            | (Exact::Bam, true, _, true)
            | (Exact::Cram, true, _, true) => {
                print_error(
                    sub_name,
                    "-d/-D tag filtering is currently supported for FASTQ single-output mode only",
                );
                return ExitCode::from(1);
            }
            (Exact::Sam, true, _, _) => view_sam_path_as_fasta_text(
                input,
                flag_filters,
                append_read_number,
                umi_tags.as_deref(),
                casava,
                barcode_tag,
            ),
            (Exact::Bam, false, _, _) => view_bam_path_as_fastq_text_with_aux(
                input,
                flag_filters,
                render_options,
                &tag_filters,
            ),
            (Exact::Bam, true, _, _) => view_bam_path_as_fasta_text(
                input,
                flag_filters,
                append_read_number,
                umi_tags.as_deref(),
                casava,
                barcode_tag,
            ),
            (Exact::Cram, false, _, _) => view_cram_path_as_fastq_text_with_aux(
                input,
                flag_filters,
                render_options,
                &tag_filters,
            ),
            (Exact::Cram, true, _, _) => view_cram_path_as_fasta_text(
                input,
                flag_filters,
                append_read_number,
                umi_tags.as_deref(),
                casava,
                barcode_tag,
            ),
            _ => {
                print_error(
                    sub_name,
                    "only SAM, BAM, and CRAM input are currently supported",
                );
                return ExitCode::from(1);
            }
        }
    };

    let text = match text {
        Ok(t) => t,
        Err(e) => {
            print_error_errno(
                sub_name,
                format!(
                    "conversion failed for \"{}\"",
                    input
                        .as_ref()
                        .map(|path| path.display().to_string())
                        .unwrap_or_else(|| "-".to_string())
                ),
                &e,
            );
            return ExitCode::from(1);
        }
    };

    let output = output.as_ref().or(other_output.as_ref());
    let mut out = match sam_io::open_text_output(output.map(PathBuf::as_path)) {
        Ok(out) => out,
        Err(e) => {
            print_error_errno(sub_name, "open -o output", &e);
            return ExitCode::from(1);
        }
    };
    if let Err(e) = out.write_all(text.as_bytes())
        && e.kind() != io::ErrorKind::BrokenPipe
    {
        print_error_errno(sub_name, "write output", &e);
        return ExitCode::from(1);
    }
    match sam_io::check_sam_close(&mut out) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => {}
        Err(e) => {
            print_error_errno(sub_name, "close output", &e);
            return ExitCode::from(1);
        }
    }

    if (index_file_1.is_some() || index_file_2.is_some())
        && let Err(e) = emit_index_files(
            input.as_deref(),
            input_exact,
            stdin_input,
            flag_filters,
            IndexEmitOptions {
                append_read_number,
                use_original_quality,
                default_quality,
                umi_tags: umi_tags.as_deref(),
                casava,
                barcode_tag,
                quality_tag,
                index_format: &effective_index_format,
                index_file_1: index_file_1.as_deref(),
                index_file_2: index_file_2.as_deref(),
                fasta_mode,
            },
        )
    {
        print_error_errno(sub_name, "index FASTQ output", &e);
        return ExitCode::from(1);
    }

    print_bam2fq_summary(FastqStats {
        processed_reads: count_fastx_records(text.as_bytes(), fasta_mode),
        discarded_singletons: 0,
    });
    ExitCode::SUCCESS
}

fn effective_fastq_input_exact(
    input: Option<&std::path::Path>,
    format: Option<&htslib_rs::format::Format>,
) -> Option<Exact> {
    let exact = format.map(|f| f.exact)?;
    if exact == Exact::Unknown
        && input
            .and_then(|path| path.extension())
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case("sam"))
    {
        return Some(Exact::Sam);
    }
    Some(exact)
}

fn write_text_file(path: &std::path::Path, text: &[u8]) -> io::Result<()> {
    let mut file = File::create(path)?;
    file.write_all(text)
}

fn print_bam2fq_summary(stats: FastqStats) {
    eprintln!(
        "[M::bam2fq_mainloop] discarded {} singletons",
        stats.discarded_singletons
    );
    eprintln!(
        "[M::bam2fq_mainloop] processed {} reads",
        stats.processed_reads
    );
}

fn count_split_records(split: &FastqSplitBuffers, fasta_mode: bool) -> usize {
    count_fastx_records(&split.read1, fasta_mode)
        + count_fastx_records(&split.read2, fasta_mode)
        + count_fastx_records(&split.singleton, fasta_mode)
        + count_fastx_records(&split.other, fasta_mode)
}

fn count_fastx_records(text: &[u8], fasta_mode: bool) -> usize {
    if text.is_empty() {
        return 0;
    }

    if fasta_mode {
        return text.iter().filter(|&&b| b == b'>').count();
    }

    text.split(|&b| b == b'\n')
        .enumerate()
        .filter(|(i, line)| i % 4 == 0 && line.first() == Some(&b'@'))
        .count()
}

fn view_sam_path_as_fastq_split(
    input: &std::path::Path,
    flag_filters: FlagFilters,
    options: FastqRenderOptions<'_>,
    tag_filters: &[TagFilter],
    singleton_set: bool,
) -> io::Result<FastqSplitBuffers> {
    let file = File::open(input)?;
    let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(file));
    view_sam_reader_as_fastq_split_with_aux(
        &mut reader,
        flag_filters,
        options,
        tag_filters,
        singleton_set,
    )
}

fn view_sam_path_as_fastq_text_with_aux(
    input: &std::path::Path,
    flag_filters: FlagFilters,
    options: FastqRenderOptions<'_>,
    tag_filters: &[TagFilter],
) -> io::Result<String> {
    let file = File::open(input)?;
    let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(file));
    view_sam_reader_as_fastq_text_with_aux(&mut reader, flag_filters, options, tag_filters)
}

fn view_sam_path_as_fasta_split(
    input: &std::path::Path,
    flag_filters: FlagFilters,
    append_read_number: bool,
    umi_tags: Option<&[[u8; 2]]>,
    casava: bool,
    barcode_tag: [u8; 2],
    singleton_set: bool,
) -> io::Result<FastqSplitBuffers> {
    let file = File::open(input)?;
    let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(file));
    view_sam_reader_as_fasta_split(
        &mut reader,
        flag_filters,
        append_read_number,
        umi_tags,
        casava,
        barcode_tag,
        singleton_set,
    )
}

fn view_sam_path_as_fasta_text(
    input: &std::path::Path,
    flag_filters: FlagFilters,
    append_read_number: bool,
    umi_tags: Option<&[[u8; 2]]>,
    casava: bool,
    barcode_tag: [u8; 2],
) -> io::Result<String> {
    let file = File::open(input)?;
    let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(file));
    view_sam_reader_as_fasta_text(
        &mut reader,
        flag_filters,
        append_read_number,
        umi_tags,
        casava,
        barcode_tag,
    )
}

fn view_bam_path_as_fastq_text_with_aux(
    input: &std::path::Path,
    flag_filters: FlagFilters,
    options: FastqRenderOptions<'_>,
    tag_filters: &[TagFilter],
) -> io::Result<String> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;
    let mut writer = Vec::new();
    let mut record = htslib_rs::sam::alignment::RecordBuf::default();

    loop {
        let n = reader.read_record_buf(&header, &mut record)?;
        if n == 0 {
            break;
        }

        if record_passes_flag_filter(&record, flag_filters)?
            && record_passes_tag_filters(&record, tag_filters)?
        {
            write_fastq_record_with_aux(&mut writer, &record, options)?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn view_bam_path_as_fasta_text(
    input: &std::path::Path,
    flag_filters: FlagFilters,
    append_read_number: bool,
    umi_tags: Option<&[[u8; 2]]>,
    casava: bool,
    barcode_tag: [u8; 2],
) -> io::Result<String> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;
    let mut writer = Vec::new();
    let mut record = htslib_rs::sam::alignment::RecordBuf::default();

    loop {
        let n = reader.read_record_buf(&header, &mut record)?;
        if n == 0 {
            break;
        }

        if record_passes_flag_filter(&record, flag_filters)? {
            write_fasta_record(
                &mut writer,
                &record,
                append_read_number,
                umi_tags,
                casava,
                barcode_tag,
            )?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn view_bam_path_as_fastq_split_with_aux(
    input: &std::path::Path,
    flag_filters: FlagFilters,
    options: FastqRenderOptions<'_>,
    tag_filters: &[TagFilter],
    singleton_set: bool,
) -> io::Result<FastqSplitBuffers> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;
    let mut split = FastqSplitBuffers::default();
    let mut grouper = GroupedSplitWriter::new(singleton_set, options.casava, options.barcode_tag);
    let mut record = htslib_rs::sam::alignment::RecordBuf::default();

    loop {
        let n = reader.read_record_buf(&header, &mut record)?;
        if n == 0 {
            break;
        }

        if record_passes_flag_filter(&record, flag_filters)?
            && record_passes_tag_filters(&record, tag_filters)?
        {
            let text = render_fastq_record_with_aux_to_vec(&record, options)?;
            grouper.add_text(&record, text, &mut split)?;
        }
    }
    grouper.flush(&mut split);

    Ok(split)
}

fn view_bam_path_as_fasta_split(
    input: &std::path::Path,
    flag_filters: FlagFilters,
    append_read_number: bool,
    umi_tags: Option<&[[u8; 2]]>,
    casava: bool,
    barcode_tag: [u8; 2],
    singleton_set: bool,
) -> io::Result<FastqSplitBuffers> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;
    let mut split = FastqSplitBuffers::default();
    let mut grouper = GroupedSplitWriter::new(singleton_set, casava, barcode_tag);
    let mut record = htslib_rs::sam::alignment::RecordBuf::default();

    loop {
        let n = reader.read_record_buf(&header, &mut record)?;
        if n == 0 {
            break;
        }

        if record_passes_flag_filter(&record, flag_filters)? {
            let text = render_fasta_record_to_vec(
                &record,
                append_read_number,
                umi_tags,
                casava,
                barcode_tag,
            )?;
            grouper.add_text(&record, text, &mut split)?;
        }
    }
    grouper.flush(&mut split);

    Ok(split)
}

fn view_cram_path_as_fastq_text_with_aux(
    input: &Path,
    flag_filters: FlagFilters,
    options: FastqRenderOptions<'_>,
    tag_filters: &[TagFilter],
) -> io::Result<String> {
    let (_header, records) = read_cram_records_for_fastq(input)?;
    let mut writer = Vec::new();

    for record in records {
        if record_passes_flag_filter(&record, flag_filters)?
            && record_passes_tag_filters(&record, tag_filters)?
        {
            write_fastq_record_with_aux(&mut writer, &record, options)?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn view_cram_path_as_fasta_text(
    input: &Path,
    flag_filters: FlagFilters,
    append_read_number: bool,
    umi_tags: Option<&[[u8; 2]]>,
    casava: bool,
    barcode_tag: [u8; 2],
) -> io::Result<String> {
    let (_header, records) = read_cram_records_for_fastq(input)?;
    let mut writer = Vec::new();

    for record in records {
        if record_passes_flag_filter(&record, flag_filters)? {
            write_fasta_record(
                &mut writer,
                &record,
                append_read_number,
                umi_tags,
                casava,
                barcode_tag,
            )?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn view_cram_path_as_fastq_split_with_aux(
    input: &Path,
    flag_filters: FlagFilters,
    options: FastqRenderOptions<'_>,
    tag_filters: &[TagFilter],
    singleton_set: bool,
) -> io::Result<FastqSplitBuffers> {
    let (_header, records) = read_cram_records_for_fastq(input)?;
    let mut split = FastqSplitBuffers::default();
    let mut grouper = GroupedSplitWriter::new(singleton_set, options.casava, options.barcode_tag);

    for record in records {
        if record_passes_flag_filter(&record, flag_filters)?
            && record_passes_tag_filters(&record, tag_filters)?
        {
            let text = render_fastq_record_with_aux_to_vec(&record, options)?;
            grouper.add_text(&record, text, &mut split)?;
        }
    }
    grouper.flush(&mut split);

    Ok(split)
}

fn view_cram_path_as_fasta_split(
    input: &Path,
    flag_filters: FlagFilters,
    append_read_number: bool,
    umi_tags: Option<&[[u8; 2]]>,
    casava: bool,
    barcode_tag: [u8; 2],
    singleton_set: bool,
) -> io::Result<FastqSplitBuffers> {
    let (_header, records) = read_cram_records_for_fastq(input)?;
    let mut split = FastqSplitBuffers::default();
    let mut grouper = GroupedSplitWriter::new(singleton_set, casava, barcode_tag);

    for record in records {
        if record_passes_flag_filter(&record, flag_filters)? {
            let text = render_fasta_record_to_vec(
                &record,
                append_read_number,
                umi_tags,
                casava,
                barcode_tag,
            )?;
            grouper.add_text(&record, text, &mut split)?;
        }
    }
    grouper.flush(&mut split);

    Ok(split)
}

fn read_cram_records_for_fastq(
    input: &Path,
) -> io::Result<(sam::Header, Vec<sam::alignment::RecordBuf>)> {
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(input)?;

    if let Some(reference) = crate::sam_global::current_global_args().reference {
        let records = htslib_rs::alignment_compat::query_cram_records_all_from_path_with_reference(
            input, &reference,
        )?;
        return Ok((header, records));
    }

    if let Some(reference) = fastq_cram_reference_from_header(&header)? {
        let records = htslib_rs::alignment_compat::query_cram_records_all_from_path_with_reference(
            input,
            reference.path(),
        )?;
        return Ok((header, records));
    }

    let records = htslib_rs::alignment_compat::query_cram_records_all_from_path(input)?;
    Ok((header, records))
}

struct FastqReferenceGuard {
    path: PathBuf,
    cleanup: Vec<PathBuf>,
}

impl FastqReferenceGuard {
    fn new(path: PathBuf) -> Self {
        Self {
            path,
            cleanup: Vec::new(),
        }
    }

    fn temporary(path: PathBuf, cleanup: Vec<PathBuf>) -> Self {
        Self { path, cleanup }
    }

    fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for FastqReferenceGuard {
    fn drop(&mut self) {
        for path in self.cleanup.iter().rev() {
            let _ = fs::remove_file(path);
        }
    }
}

fn fastq_cram_reference_from_header(
    header: &sam::Header,
) -> io::Result<Option<FastqReferenceGuard>> {
    let mut header_text = Vec::new();
    crate::sam_render::write_header(&mut header_text, header)?;
    let header_text = String::from_utf8(header_text)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    if let Some(reference) = reference_from_header_uri_text(&header_text) {
        return Ok(Some(FastqReferenceGuard::new(reference)));
    }

    let Some(ref_path) = std::env::var_os("REF_PATH") else {
        return Ok(None);
    };
    reference_from_ref_path_header(&header_text, &ref_path.to_string_lossy())
}

fn reference_from_header_uri_text(header_text: &str) -> Option<PathBuf> {
    for line in header_text.lines().filter(|line| line.starts_with("@SQ\t")) {
        for field in line.split('\t').skip(1) {
            let Some(uri) = field.strip_prefix("UR:") else {
                continue;
            };
            let path = PathBuf::from(uri.strip_prefix("file://").unwrap_or(uri));
            if path.is_file() {
                return Some(path);
            }
        }
    }
    None
}

fn reference_from_ref_path_header(
    header_text: &str,
    ref_path: &str,
) -> io::Result<Option<FastqReferenceGuard>> {
    let mut sequences = Vec::new();

    for line in header_text.lines().filter(|line| line.starts_with("@SQ\t")) {
        let mut name = None;
        let mut md5 = None;
        let mut len = None;
        for field in line.split('\t').skip(1) {
            if let Some(value) = field.strip_prefix("SN:") {
                name = Some(value);
            } else if let Some(value) = field.strip_prefix("LN:") {
                len = value.parse::<usize>().ok();
            } else if let Some(value) = field.strip_prefix("M5:") {
                md5 = Some(value);
            }
        }

        let (Some(name), Some(md5)) = (name, md5) else {
            continue;
        };
        if let Some(sequence) = read_ref_path_md5_sequence(ref_path, md5)? {
            sequences.push((name.to_string(), sequence));
        } else if let Some(len) = len {
            sequences.push((name.to_string(), "N".repeat(len)));
        }
    }

    if sequences.is_empty() {
        return Ok(None);
    }

    let fasta = temporary_reference_path("fastq-ref-path", "fa");
    {
        let mut out = File::create(&fasta)?;
        for (name, sequence) in &sequences {
            writeln!(out, ">{name}")?;
            out.write_all(sequence.as_bytes())?;
            out.write_all(b"\n")?;
        }
    }
    let fai = crate::reference::ensure_fai_index(&fasta, None)?;
    Ok(Some(FastqReferenceGuard::temporary(
        fasta.clone(),
        vec![fai, fasta],
    )))
}

fn read_ref_path_md5_sequence(ref_path: &str, md5: &str) -> io::Result<Option<String>> {
    for template in ref_path.split(':').filter(|part| !part.is_empty()) {
        let candidate = if template.contains("%s") {
            PathBuf::from(template.replace("%s", md5))
        } else {
            Path::new(template).join(md5)
        };
        match fs::read_to_string(&candidate) {
            Ok(text) => {
                let sequence: String = text.chars().filter(|c| !c.is_whitespace()).collect();
                if !sequence.is_empty() {
                    return Ok(Some(sequence));
                }
            }
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => return Err(e),
        }
    }
    Ok(None)
}

fn temporary_reference_path(stem: &str, ext: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        "samtools-rs-{stem}-{}-{nanos}.{ext}",
        std::process::id()
    ))
}

fn view_sam_reader_as_fastq_text_with_aux<R>(
    reader: &mut htslib_rs::sam::io::Reader<R>,
    flag_filters: FlagFilters,
    options: FastqRenderOptions<'_>,
    tag_filters: &[TagFilter],
) -> io::Result<String>
where
    R: io::BufRead,
{
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, flag_filters)?
            && record_passes_tag_filters(&record, tag_filters)?
        {
            write_fastq_record_with_aux(&mut writer, &record, options)?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn view_sam_reader_as_fasta_text<R>(
    reader: &mut htslib_rs::sam::io::Reader<R>,
    flag_filters: FlagFilters,
    append_read_number: bool,
    umi_tags: Option<&[[u8; 2]]>,
    casava: bool,
    barcode_tag: [u8; 2],
) -> io::Result<String>
where
    R: io::BufRead,
{
    let _header = reader.read_header()?;
    let mut writer = Vec::new();

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, flag_filters)? {
            write_fasta_record(
                &mut writer,
                &record,
                append_read_number,
                umi_tags,
                casava,
                barcode_tag,
            )?;
        }
    }

    String::from_utf8(writer).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn view_sam_reader_as_fastq_split_with_aux<R>(
    reader: &mut htslib_rs::sam::io::Reader<R>,
    flag_filters: FlagFilters,
    options: FastqRenderOptions<'_>,
    tag_filters: &[TagFilter],
    singleton_set: bool,
) -> io::Result<FastqSplitBuffers>
where
    R: io::BufRead,
{
    let _header = reader.read_header()?;
    let mut split = FastqSplitBuffers::default();
    let mut grouper = GroupedSplitWriter::new(singleton_set, options.casava, options.barcode_tag);

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, flag_filters)?
            && record_passes_tag_filters(&record, tag_filters)?
        {
            let text = render_fastq_record_with_aux_to_vec(&record, options)?;
            grouper.add_text(&record, text, &mut split)?;
        }
    }
    grouper.flush(&mut split);

    Ok(split)
}

fn view_sam_reader_as_fasta_split<R>(
    reader: &mut htslib_rs::sam::io::Reader<R>,
    flag_filters: FlagFilters,
    append_read_number: bool,
    umi_tags: Option<&[[u8; 2]]>,
    casava: bool,
    barcode_tag: [u8; 2],
    singleton_set: bool,
) -> io::Result<FastqSplitBuffers>
where
    R: io::BufRead,
{
    let _header = reader.read_header()?;
    let mut split = FastqSplitBuffers::default();
    let mut grouper = GroupedSplitWriter::new(singleton_set, casava, barcode_tag);

    for result in reader.records() {
        let record = result?;
        if record_passes_flag_filter(&record, flag_filters)? {
            let text = render_fasta_record_to_vec(
                &record,
                append_read_number,
                umi_tags,
                casava,
                barcode_tag,
            )?;
            grouper.add_text(&record, text, &mut split)?;
        }
    }
    grouper.flush(&mut split);

    Ok(split)
}

#[derive(Default)]
struct FastqSplitBuffers {
    read1: Vec<u8>,
    read2: Vec<u8>,
    singleton: Vec<u8>,
    other: Vec<u8>,
}

#[derive(Clone, Copy)]
enum ReadPart {
    Other = 0,
    Read1 = 1,
    Read2 = 2,
}

fn record_read_part<R>(record: &R) -> io::Result<ReadPart>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let flags = record.flags()?;
    Ok(if flags.is_first_segment() && !flags.is_last_segment() {
        ReadPart::Read1
    } else if flags.is_last_segment() && !flags.is_first_segment() {
        ReadPart::Read2
    } else {
        ReadPart::Other
    })
}

fn record_score<R>(record: &R) -> u8
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let qual = record.quality_scores();
    let has_quality = qual.iter().next().is_some_and(|res| match res {
        Ok(score) => score != 0xff,
        Err(_) => false,
    });
    if has_quality { 2 } else { 1 }
}

/// Buffers FASTQ/FASTA text for the current qname group, flushing into a
/// `FastqSplitBuffers` when the qname changes, applying upstream
/// `bam_fastq.c::flush_rec` routing.
///
/// When `singleton_set` is true (a `-s` file is configured), R1-only and
/// R2-only singletons go into `singleton`; otherwise they fall back to the
/// `read1` / `read2` buffers respectively, matching upstream's `fpse` /
/// `fpr[1]` / `fpr[2]` fallback.
#[derive(Default)]
struct GroupedSplitWriter {
    singleton_set: bool,
    /// `-i` (CASAVA): when set, the group's barcode is propagated to
    /// every mate's CASAVA comment that lacked its own (upstream copies
    /// the `BC` aux across mates before formatting — `bam_fastq.c:952`).
    casava: bool,
    barcode_tag: [u8; 2],
    current_qname: Option<Vec<u8>>,
    best_score: [u8; 3],
    pending_text: [Option<Vec<u8>>; 3],
    /// First non-empty raw barcode seen in the current qname group.
    group_barcode: Option<String>,
}

impl GroupedSplitWriter {
    fn new(singleton_set: bool, casava: bool, barcode_tag: [u8; 2]) -> Self {
        Self {
            singleton_set,
            casava,
            barcode_tag,
            ..Self::default()
        }
    }

    fn add_text<R>(
        &mut self,
        record: &R,
        text: Vec<u8>,
        split: &mut FastqSplitBuffers,
    ) -> io::Result<()>
    where
        R: htslib_rs::sam::alignment::Record + ?Sized,
    {
        let qname = record
            .name()
            .map(|name| name.to_vec())
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing record name"))?;
        if self.current_qname.as_deref() != Some(&qname) {
            self.flush(split);
            self.current_qname = Some(qname);
            self.group_barcode = None;
        }

        if self.casava
            && self.group_barcode.is_none()
            && let Some(bc) = fastq_string_tag(record, self.barcode_tag)?
            && !bc.is_empty()
        {
            self.group_barcode = Some(bc);
        }

        let part = record_read_part(record)? as usize;
        let score = record_score(record);
        if score > self.best_score[part] {
            self.best_score[part] = score;
            self.pending_text[part] = Some(text);
        }
        Ok(())
    }

    fn flush(&mut self, split: &mut FastqSplitBuffers) {
        let [s0, s1, s2] = self.best_score;
        let [mut t0, mut t1, mut t2] = std::mem::take(&mut self.pending_text);

        // Propagate the group's barcode into any mate whose CASAVA
        // comment had no barcode of its own (upstream copies the `BC`
        // aux across mates before formatting the comment).
        if self.casava
            && let Some(bc) = self.group_barcode.as_deref()
        {
            let bc_field = casava_barcode_field(Some(bc));
            for t in [&mut t0, &mut t1, &mut t2] {
                if let Some(text) = t.as_mut() {
                    fill_casava_barcode(text, &bc_field);
                }
            }
        }

        if s1 > 0 && s2 > 0 {
            if let Some(t) = t1 {
                split.read1.extend_from_slice(&t);
            }
            if let Some(t) = t2 {
                split.read2.extend_from_slice(&t);
            }
        } else if s1 > 0
            && let Some(t) = t1
        {
            if self.singleton_set {
                split.singleton.extend_from_slice(&t);
            } else {
                split.read1.extend_from_slice(&t);
            }
        } else if s2 > 0
            && let Some(t) = t2
        {
            if self.singleton_set {
                split.singleton.extend_from_slice(&t);
            } else {
                split.read2.extend_from_slice(&t);
            }
        }

        if s0 > 0
            && let Some(t) = t0
        {
            split.other.extend_from_slice(&t);
        }

        self.best_score = [0; 3];
    }
}

/// Interleaves matched paired records: pair `k` from `read1` is written
/// before pair `k` from `read2`. Both buffers must contain the same number
/// of records (one record per name-grouped paired flush). Records are
/// delimited by 4 lines each: header, sequence, `+` separator, quality
/// (FASTQ) or 2 lines (FASTA: header, sequence).
fn interleave_paired_records(read1: &[u8], read2: &[u8]) -> Vec<u8> {
    let r1 = split_fastx_records(read1);
    let r2 = split_fastx_records(read2);
    let n = r1.len().min(r2.len());
    let mut out =
        Vec::with_capacity(read1.len() + read2.len() + 2 * (r1.len() + r2.len() - 2 * n) * 4);
    for i in 0..n {
        out.extend_from_slice(r1[i]);
        out.extend_from_slice(r2[i]);
    }
    for record in &r1[n..] {
        out.extend_from_slice(record);
    }
    for record in &r2[n..] {
        out.extend_from_slice(record);
    }
    out
}

fn split_fastx_records(text: &[u8]) -> Vec<&[u8]> {
    let mut records = Vec::new();
    let mut start = 0;
    for (i, &b) in text.iter().enumerate() {
        if b == b'\n' && i + 1 < text.len() && matches!(text[i + 1], b'@' | b'>') {
            records.push(&text[start..=i]);
            start = i + 1;
        }
    }
    if start < text.len() {
        records.push(&text[start..]);
    }
    records
}

fn render_fastq_record_with_aux_to_vec<R>(
    record: &R,
    options: FastqRenderOptions<'_>,
) -> io::Result<Vec<u8>>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let mut buf = Vec::new();
    write_fastq_record_with_aux(&mut buf, record, options)?;
    Ok(buf)
}

fn render_fasta_record_to_vec<R>(
    record: &R,
    append_read_number: bool,
    umi_tags: Option<&[[u8; 2]]>,
    casava: bool,
    barcode_tag: [u8; 2],
) -> io::Result<Vec<u8>>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let mut buf = Vec::new();
    write_fasta_record(
        &mut buf,
        record,
        append_read_number,
        umi_tags,
        casava,
        barcode_tag,
    )?;
    Ok(buf)
}

fn record_passes_tag_filters<R>(record: &R, filters: &[TagFilter]) -> io::Result<bool>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    if filters.is_empty() {
        return Ok(true);
    }

    let mut found = vec![false; filters.len()];
    for result in record.data().iter() {
        let (tag, value) = result?;
        let tag_bytes = <[u8; 2]>::from(tag);
        let value_payload = if filters
            .iter()
            .any(|filter| filter.tag == tag_bytes && filter.values.is_some())
        {
            Some(aux_value_payload(value)?)
        } else {
            None
        };

        for (i, filter) in filters.iter().enumerate() {
            if found[i] || filter.tag != tag_bytes {
                continue;
            }

            found[i] = match &filter.values {
                Some(values) => value_payload
                    .as_ref()
                    .is_some_and(|payload| values.contains(payload)),
                None => true,
            };
        }
    }

    Ok(found.into_iter().all(|matched| matched))
}

fn record_passes_flag_filter<R>(record: &R, filters: FlagFilters) -> io::Result<bool>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let flag = record.flags()?.bits();
    Ok(
        (filters.require == 0 || (flag & filters.require) == filters.require)
            && (filters.include_any == 0 || (flag & filters.include_any) != 0)
            && (filters.exclude == 0 || (flag & filters.exclude) == 0)
            && (filters.exclude_all == 0 || (flag & filters.exclude_all) != filters.exclude_all),
    )
}

fn write_fastq_record_with_aux<W, R>(
    writer: &mut W,
    record: &R,
    options: FastqRenderOptions<'_>,
) -> io::Result<()>
where
    W: Write,
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let name = fastq_record_name(record)?;
    let name = append_fastq_umi(name, record, options.umi_tags)?;
    let name = append_fastq_read_number(name, record, options.append_read_number)?;
    let payload = fastq_payload(record, options)?;
    let sequence = payload.sequence;
    let quality = payload.quality;
    if options.strip_soft_clips && sequence.is_empty() {
        return Ok(());
    }

    if sequence.len() != quality.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("FASTQ quality length differs from sequence length for {name}"),
        ));
    }

    write!(writer, "@{name}")?;
    if options.casava {
        write!(
            writer,
            " {}",
            fastq_casava_comment(record, options.barcode_tag)?
        )?;
    } else {
        for field in fastq_aux_fields_with_soft_clip(
            record,
            options.aux_selection,
            payload.soft_clip_aux.as_ref(),
        )? {
            write!(writer, "\t{field}")?;
        }
    }
    writeln!(writer)?;
    writeln!(writer, "{sequence}")?;
    writeln!(writer, "+")?;
    writeln!(writer, "{quality}")?;

    Ok(())
}

fn write_fasta_record<W, R>(
    writer: &mut W,
    record: &R,
    append_read_number: bool,
    umi_tags: Option<&[[u8; 2]]>,
    casava: bool,
    barcode_tag: [u8; 2],
) -> io::Result<()>
where
    W: Write,
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let name = fastq_record_name(record)?;
    let name = append_fastq_umi(name, record, umi_tags)?;
    let name = append_fastq_read_number(name, record, append_read_number)?;
    let sequence = fastq_sequence_string(record);

    write!(writer, ">{name}")?;
    if casava {
        write!(writer, " {}", fastq_casava_comment(record, barcode_tag)?)?;
    }
    writeln!(writer)?;
    writeln!(writer, "{sequence}")?;

    Ok(())
}

struct FastqPayload {
    sequence: String,
    quality: String,
    soft_clip_aux: Option<SoftClipAux>,
}

struct SoftClipAux {
    tag: [u8; 2],
    field: String,
}

fn fastq_payload<R>(record: &R, options: FastqRenderOptions<'_>) -> io::Result<FastqPayload>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    if !options.strip_soft_clips {
        return Ok(FastqPayload {
            sequence: fastq_sequence_string(record),
            quality: fastq_quality_scores_string(
                record,
                options.use_original_quality,
                options.default_quality,
            )?,
            soft_clip_aux: None,
        });
    }

    let seq = record.sequence().iter().collect::<Vec<_>>();
    let qual = fastq_quality_scores_storage_string(
        record,
        options.use_original_quality,
        options.default_quality,
    )?
    .into_bytes();
    if seq.len() != qual.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTQ quality length differs from sequence length",
        ));
    }

    let mut kept_seq = Vec::with_capacity(seq.len());
    let mut kept_qual = Vec::with_capacity(qual.len());
    let mut clipped_seq = Vec::new();
    let mut clipped_qual = Vec::new();
    let mut cigar_ops = Vec::new();

    if record.cigar().as_ref().is_empty() {
        kept_seq.extend_from_slice(&seq);
        kept_qual.extend_from_slice(&qual);
    } else {
        use htslib_rs::sam::alignment::record::cigar::op::Kind;

        let mut cursor = 0usize;
        for op in record.cigar().as_ref() {
            let op = op?;
            let kind = op.kind();
            cigar_ops.push((op.len(), kind));
            if !kind.consumes_read() {
                continue;
            }

            let len = op.len();
            let end = cursor.checked_add(len).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CIGAR read length overflow")
            })?;
            if end > seq.len() || end > qual.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CIGAR consumes more read bases than available",
                ));
            }

            if kind == Kind::SoftClip {
                clipped_seq.extend_from_slice(&seq[cursor..end]);
                clipped_qual.extend_from_slice(&qual[cursor..end]);
            } else {
                kept_seq.extend_from_slice(&seq[cursor..end]);
                kept_qual.extend_from_slice(&qual[cursor..end]);
            }
            cursor = end;
        }
    }

    let is_reverse = record.flags()?.is_reverse_complemented();
    if is_reverse {
        kept_seq = kept_seq.into_iter().rev().map(complement_base).collect();
        clipped_seq = clipped_seq.into_iter().rev().map(complement_base).collect();
        kept_qual.reverse();
        clipped_qual.reverse();
        cigar_ops.reverse();
    }

    let soft_clip_aux = options
        .soft_clip_backup_tag
        .filter(|_| !clipped_seq.is_empty())
        .map(|tag| {
            let cigar = cigar_ops
                .iter()
                .map(|(len, kind)| format!("{len}{}", cigar_kind_char(*kind)))
                .collect::<String>();
            let clipped_seq = String::from_utf8_lossy(&clipped_seq);
            let clipped_qual = String::from_utf8_lossy(&clipped_qual);
            let field = format!(
                "{}{}:Z:{cigar}:{clipped_seq}:{clipped_qual}",
                char::from(tag[0]),
                char::from(tag[1])
            );
            SoftClipAux { tag, field }
        });

    Ok(FastqPayload {
        sequence: String::from_utf8_lossy(&kept_seq).into_owned(),
        quality: String::from_utf8(kept_qual)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        soft_clip_aux,
    })
}

fn cigar_kind_char(kind: htslib_rs::sam::alignment::record::cigar::op::Kind) -> char {
    use htslib_rs::sam::alignment::record::cigar::op::Kind;

    match kind {
        Kind::Match => 'M',
        Kind::Insertion => 'I',
        Kind::Deletion => 'D',
        Kind::Skip => 'N',
        Kind::SoftClip => 'S',
        Kind::HardClip => 'H',
        Kind::Pad => 'P',
        Kind::SequenceMatch => '=',
        Kind::SequenceMismatch => 'X',
    }
}

fn fastq_record_name<R>(record: &R) -> io::Result<String>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    record
        .name()
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing FASTQ record name"))
}

fn append_fastq_read_number<R>(
    mut name: String,
    record: &R,
    append_read_number: bool,
) -> io::Result<String>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    if append_read_number {
        let flags = record.flags()?;
        if flags.is_first_segment() {
            name.push_str("/1");
        } else if flags.is_last_segment() {
            name.push_str("/2");
        }
    }

    Ok(name)
}

fn append_fastq_umi<R>(
    mut name: String,
    record: &R,
    umi_tags: Option<&[[u8; 2]]>,
) -> io::Result<String>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let Some(umi_tags) = umi_tags else {
        return Ok(name);
    };
    let Some(umi) = fastq_umi_string(record, umi_tags)? else {
        return Ok(name);
    };

    let umi = umi
        .chars()
        .map(|c| if c.is_ascii_alphabetic() { c } else { '+' })
        .collect::<String>();
    if let Some(hash) = name.rfind('#') {
        name.insert_str(hash, &format!(":{umi}"));
    } else {
        name.push(':');
        name.push_str(&umi);
    }

    Ok(name)
}

fn fastq_umi_string<R>(record: &R, umi_tags: &[[u8; 2]]) -> io::Result<Option<String>>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    use htslib_rs::sam::alignment::record::data::field::{Tag, Value};

    let data = record.data();
    for tag in umi_tags {
        let tag = Tag::from(*tag);
        let Some(value) = data.get(&tag).transpose()? else {
            continue;
        };
        if let Value::String(s) = value {
            return Ok(Some(String::from_utf8_lossy(s).into_owned()));
        }
    }

    Ok(None)
}

fn fastq_casava_comment<R>(record: &R, barcode_tag: [u8; 2]) -> io::Result<String>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let flags = record.flags()?;
    let read_number = if flags.is_last_segment() { 2 } else { 1 };
    let filter = if flags.is_qc_fail() { "Y" } else { "N" };
    let barcode = casava_barcode_field(fastq_string_tag(record, barcode_tag)?.as_deref());

    Ok(format!("{read_number}:{filter}:0:{barcode}"))
}

/// If the FASTQ/FASTA header line of `text` ends with a CASAVA comment
/// whose barcode is the placeholder `0` (the record had no `BC` of its
/// own), replaces that `0` with `bc_field`. A record that already
/// carries its own barcode is left untouched (mirrors upstream copying
/// `BC` only into mates that lack it).
fn fill_casava_barcode(text: &mut Vec<u8>, bc_field: &str) {
    let hdr_end = text.iter().position(|&b| b == b'\n').unwrap_or(text.len());
    let Some(sp) = text[..hdr_end].iter().rposition(|&b| b == b' ') else {
        return;
    };
    let token_start = sp + 1;
    let token = &text[token_start..hdr_end];
    // Shape: `<rnum>:<filt>:0:<bc>` — digit, ':', Y/N, ':', '0', ':'.
    if token.len() < 6
        || !token[0].is_ascii_digit()
        || token[1] != b':'
        || token[3] != b':'
        || token[4] != b'0'
        || token[5] != b':'
    {
        return;
    }
    let bc_start = token_start + 6;
    if &text[bc_start..hdr_end] != b"0" {
        return; // record had its own barcode; leave it.
    }
    let mut new = Vec::with_capacity(text.len() - 1 + bc_field.len());
    new.extend_from_slice(&text[..bc_start]);
    new.extend_from_slice(bc_field.as_bytes());
    new.extend_from_slice(&text[hdr_end..]);
    *text = new;
}

/// Renders the barcode portion of a CASAVA comment exactly as htslib's
/// `fastq_format1` does: `0` when absent or when the first character is
/// not a sequence base; otherwise every non-alphabetic byte becomes `+`
/// and lowercase is upper-cased (so `ac-gt` → `AC+GT`).
fn casava_barcode_field(bc: Option<&str>) -> String {
    let Some(bc) = bc.filter(|b| !b.is_empty()) else {
        return "0".to_string();
    };
    let first = bc.as_bytes()[0];
    if !first.is_ascii_alphabetic() {
        return "0".to_string();
    }
    bc.chars()
        .map(|c| {
            if c.is_ascii_alphabetic() {
                c.to_ascii_uppercase()
            } else {
                '+'
            }
        })
        .collect()
}

fn fastq_string_tag<R>(record: &R, tag: [u8; 2]) -> io::Result<Option<String>>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    use htslib_rs::sam::alignment::record::data::field::{Tag, Value};

    let data = record.data();
    let Some(value) = data.get(&Tag::from(tag)).transpose()? else {
        return Ok(None);
    };

    match value {
        Value::String(s) => Ok(Some(String::from_utf8_lossy(s).into_owned())),
        _ => Ok(None),
    }
}

fn fastq_sequence_string<R>(record: &R) -> String
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let sequence = record.sequence();
    let mut bases = sequence.iter().collect::<Vec<_>>();
    let is_reverse = record
        .flags()
        .map(|flags| flags.is_reverse_complemented())
        .unwrap_or(false);

    if is_reverse {
        bases = bases.into_iter().rev().map(complement_base).collect();
    }

    String::from_utf8_lossy(&bases).into_owned()
}

fn fastq_quality_scores_string<R>(
    record: &R,
    use_original_quality: bool,
    default_quality: Option<u8>,
) -> io::Result<String>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let mut quality =
        fastq_quality_scores_storage_string(record, use_original_quality, default_quality)?;
    if record.flags()?.is_reverse_complemented() {
        quality = quality.chars().rev().collect();
    }
    Ok(quality)
}

fn fastq_quality_scores_storage_string<R>(
    record: &R,
    use_original_quality: bool,
    default_quality: Option<u8>,
) -> io::Result<String>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    if use_original_quality && let Some(oq) = original_quality_string(record)? {
        return Ok(oq);
    }

    let scores = record
        .quality_scores()
        .iter()
        .collect::<io::Result<Vec<_>>>()?;
    if scores.is_empty()
        && let Some(default_quality) = default_quality
    {
        let len = record.sequence().iter().count();
        return Ok(std::iter::repeat_n(char::from(default_quality + b'!'), len).collect());
    }
    let bytes = scores
        .into_iter()
        .map(|score| {
            score.checked_add(b'!').ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "FASTQ quality score overflow")
            })
        })
        .collect::<io::Result<Vec<_>>>()?;

    String::from_utf8(bytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn original_quality_string<R>(record: &R) -> io::Result<Option<String>>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    use htslib_rs::sam::alignment::record::data::field::{Tag, Value};

    let tag = Tag::from([b'O', b'Q']);
    let data = record.data();
    let Some(value) = data.get(&tag).transpose()? else {
        return Ok(None);
    };

    match value {
        Value::String(s) => Ok(Some(String::from_utf8_lossy(s).into_owned())),
        _ => Ok(None),
    }
}

fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'M' => b'K',
        b'R' => b'Y',
        b'W' => b'W',
        b'S' => b'S',
        b'Y' => b'R',
        b'K' => b'M',
        b'V' => b'B',
        b'H' => b'D',
        b'D' => b'H',
        b'B' => b'V',
        b'N' => b'N',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'm' => b'k',
        b'r' => b'y',
        b'w' => b'w',
        b's' => b's',
        b'y' => b'r',
        b'k' => b'm',
        b'v' => b'b',
        b'h' => b'd',
        b'd' => b'h',
        b'b' => b'v',
        b'n' => b'n',
        _ => base,
    }
}

fn fastq_aux_fields<R>(record: &R, aux_selection: &AuxSelection) -> io::Result<Vec<String>>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    if matches!(aux_selection, AuxSelection::None) {
        return Ok(Vec::new());
    }

    record
        .data()
        .iter()
        .filter_map(|result| match result {
            Ok((tag, value)) => {
                let tag_bytes = <[u8; 2]>::from(tag);
                match aux_selection {
                    AuxSelection::All => Some(Ok((tag_bytes, value))),
                    AuxSelection::Tags(tags) if tags.iter().any(|wanted| wanted == &tag_bytes) => {
                        Some(Ok((tag_bytes, value)))
                    }
                    AuxSelection::Tags(_) | AuxSelection::None => None,
                }
            }
            Err(e) => Some(Err(e)),
        })
        .filter_map(|result| match result {
            Ok((tag, value)) => match format_fastq_aux_field(&tag, value) {
                Ok(Some(field)) => Some(Ok(field)),
                Ok(None) => None,
                Err(e) => Some(Err(e)),
            },
            Err(e) => Some(Err(e)),
        })
        .collect()
}

fn fastq_aux_fields_with_soft_clip<R>(
    record: &R,
    aux_selection: &AuxSelection,
    soft_clip_aux: Option<&SoftClipAux>,
) -> io::Result<Vec<String>>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let Some(soft_clip_aux) = soft_clip_aux else {
        return fastq_aux_fields(record, aux_selection);
    };
    if !aux_selection_includes(aux_selection, soft_clip_aux.tag) {
        return fastq_aux_fields(record, aux_selection);
    }

    use htslib_rs::sam::alignment::record::data::field::Tag;

    let soft_clip_tag = Tag::from(soft_clip_aux.tag);
    let mut fields = record
        .data()
        .iter()
        .filter_map(|result| match result {
            Ok((tag, value)) => {
                if tag == soft_clip_tag {
                    return None;
                }
                let tag_bytes = <[u8; 2]>::from(tag);
                if !aux_selection_includes(aux_selection, tag_bytes) {
                    return None;
                }
                match format_fastq_aux_field(&tag_bytes, value) {
                    Ok(Some(field)) => Some(Ok(field)),
                    Ok(None) => None,
                    Err(e) => Some(Err(e)),
                }
            }
            Err(e) => Some(Err(e)),
        })
        .collect::<io::Result<Vec<_>>>()?;
    fields.push(soft_clip_aux.field.clone());
    Ok(fields)
}

fn aux_selection_includes(selection: &AuxSelection, tag: [u8; 2]) -> bool {
    match selection {
        AuxSelection::All => true,
        AuxSelection::Tags(tags) => tags.iter().any(|wanted| wanted == &tag),
        AuxSelection::None => false,
    }
}

fn format_fastq_aux_field(
    tag: &[u8; 2],
    value: htslib_rs::sam::alignment::record::data::field::Value<'_>,
) -> io::Result<Option<String>> {
    use htslib_rs::sam::alignment::record::data::field::{Value, value::Array};

    let tag = match std::str::from_utf8(tag) {
        Ok(tag) => tag,
        Err(_) => return Ok(None),
    };

    let field = match value {
        Value::Character(n) => format!("{tag}:A:{}", char::from(n)),
        Value::Int8(n) => format!("{tag}:i:{n}"),
        Value::UInt8(n) => format!("{tag}:i:{n}"),
        Value::Int16(n) => format!("{tag}:i:{n}"),
        Value::UInt16(n) => format!("{tag}:i:{n}"),
        Value::Int32(n) => format!("{tag}:i:{n}"),
        Value::UInt32(n) => format!("{tag}:i:{n}"),
        Value::Float(n) => format!("{tag}:f:{}", format_aux_float(n)),
        Value::String(s) => format!("{tag}:Z:{}", String::from_utf8_lossy(s)),
        Value::Hex(s) => format!("{tag}:H:{}", String::from_utf8_lossy(s)),
        Value::Array(Array::Int8(values)) => format!("{tag}:B:c,{}", join_array(values.iter())?),
        Value::Array(Array::UInt8(values)) => format!("{tag}:B:C,{}", join_array(values.iter())?),
        Value::Array(Array::Int16(values)) => format!("{tag}:B:s,{}", join_array(values.iter())?),
        Value::Array(Array::UInt16(values)) => format!("{tag}:B:S,{}", join_array(values.iter())?),
        Value::Array(Array::Int32(values)) => format!("{tag}:B:i,{}", join_array(values.iter())?),
        Value::Array(Array::UInt32(values)) => format!("{tag}:B:I,{}", join_array(values.iter())?),
        Value::Array(Array::Float(values)) => {
            format!("{tag}:B:f,{}", join_float_array(values.iter())?)
        }
    };

    Ok(Some(field))
}

fn aux_value_payload(
    value: htslib_rs::sam::alignment::record::data::field::Value<'_>,
) -> io::Result<String> {
    use htslib_rs::sam::alignment::record::data::field::{Value, value::Array};

    let payload = match value {
        Value::Character(n) => char::from(n).to_string(),
        Value::Int8(n) => n.to_string(),
        Value::UInt8(n) => n.to_string(),
        Value::Int16(n) => n.to_string(),
        Value::UInt16(n) => n.to_string(),
        Value::Int32(n) => n.to_string(),
        Value::UInt32(n) => n.to_string(),
        Value::Float(n) => format_aux_float(n),
        Value::String(s) | Value::Hex(s) => String::from_utf8_lossy(s).into_owned(),
        Value::Array(Array::Int8(values)) => format!("c,{}", join_array(values.iter())?),
        Value::Array(Array::UInt8(values)) => format!("C,{}", join_array(values.iter())?),
        Value::Array(Array::Int16(values)) => format!("s,{}", join_array(values.iter())?),
        Value::Array(Array::UInt16(values)) => format!("S,{}", join_array(values.iter())?),
        Value::Array(Array::Int32(values)) => format!("i,{}", join_array(values.iter())?),
        Value::Array(Array::UInt32(values)) => format!("I,{}", join_array(values.iter())?),
        Value::Array(Array::Float(values)) => format!("f,{}", join_float_array(values.iter())?),
    };

    Ok(payload)
}

fn join_array<N>(iter: Box<dyn Iterator<Item = io::Result<N>> + '_>) -> io::Result<String>
where
    N: std::fmt::Display,
{
    let mut values = Vec::new();
    for result in iter {
        values.push(result?.to_string());
    }
    Ok(values.join(","))
}

fn join_float_array(iter: Box<dyn Iterator<Item = io::Result<f32>> + '_>) -> io::Result<String> {
    let mut values = Vec::new();
    for result in iter {
        values.push(format_aux_float(result?));
    }
    Ok(values.join(","))
}

fn parse_flag_arg(arg: Option<&OsString>, opt: &str, sub_name: &str) -> Result<u16, ExitCode> {
    let Some(raw) = arg.and_then(|a| a.to_str()) else {
        print_error(sub_name, format!("missing value for {}", opt));
        return Err(ExitCode::from(1));
    };
    match crate::bam_flag::str_to_flag(raw) {
        Some(flag) => u16::try_from(flag).map_err(|_| {
            print_error(sub_name, format!("Could not parse \"{}\"", raw));
            ExitCode::from(1)
        }),
        None => {
            print_error(sub_name, format!("Could not parse \"{}\"", raw));
            Err(ExitCode::from(1))
        }
    }
}

fn parse_tag_filter(raw: &str) -> Result<TagFilter, String> {
    let (tag, value) = match raw.split_once(':') {
        Some((tag, value)) => (parse_filter_tag(tag)?, Some(value.to_string())),
        None => (parse_filter_tag(raw)?, None),
    };

    Ok(TagFilter {
        tag,
        values: value.map(|value| HashSet::from([value])),
    })
}

fn parse_tag_filter_file(raw: &str) -> Result<TagFilter, String> {
    let Some((tag, path)) = raw.split_once(':') else {
        return Err(format!("invalid -D value \"{raw}\": expected TAG:FILE"));
    };
    let tag = parse_filter_tag(tag)?;
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("failed to read -D value file \"{path}\": {e}"))?;
    let values = text
        .lines()
        .filter(|line| !line.is_empty())
        .map(ToOwned::to_owned)
        .collect();

    Ok(TagFilter {
        tag,
        values: Some(values),
    })
}

fn parse_filter_tag(raw: &str) -> Result<[u8; 2], String> {
    let bytes = raw.as_bytes();
    if bytes.len() != 2 {
        return Err(format!(
            "invalid tag filter \"{raw}\": expected two-character tag"
        ));
    }

    Ok([bytes[0], bytes[1]])
}

fn parse_ordered_tag_list(raw: &str) -> Result<Vec<[u8; 2]>, String> {
    let mut tags = Vec::new();
    for tag in raw.split(',') {
        let bytes = tag.as_bytes();
        if bytes.len() != 2 {
            return Err("auxiliary tags should be exactly two characters long".to_string());
        }
        if !tags.iter().any(|seen| seen == bytes) {
            tags.push([bytes[0], bytes[1]]);
        }
    }

    Ok(tags)
}

fn print_usage(sub: &str) -> io::Result<()> {
    let mut w = io::stderr().lock();
    let suffix = if sub == "fasta" { "fasta" } else { "fastq" };
    writeln!(
        w,
        "Usage: samtools {} [options] <in.bam>  > out.{}",
        sub, suffix
    )?;
    writeln!(w, "  -o FILE      write output to FILE (default stdout)")?;
    writeln!(
        w,
        "  -0 FILE      write all reads to FILE in single-output mode"
    )?;
    writeln!(w, "  -n           do not append /1 or /2 to read names")?;
    writeln!(w, "  -N           append /1 or /2 to read names")?;
    writeln!(w, "  -O           use OQ tag qualities when present")?;
    writeln!(w, "  -i           add Illumina CASAVA 1.8 fields")?;
    writeln!(w, "  --barcode-tag TAG")?;
    writeln!(w, "               aux tag to use for CASAVA barcodes [BC]")?;
    writeln!(
        w,
        "  -v INT       default quality score for missing qualities"
    )?;
    writeln!(
        w,
        "  -U, --UMI    append UMI aux tag sequence to read names"
    )?;
    writeln!(w, "  --UMI-tag TAGLIST")?;
    writeln!(
        w,
        "               aux tags to search for UMI sequence [OX,RX]"
    )?;
    writeln!(w, "  -T TAGLIST   copy aux tags to FASTQ comments")?;
    writeln!(w, "  -d, --tag TAG[:VAL] filter by aux tag presence/value")?;
    writeln!(
        w,
        "  -D, --tag-file TAG:FILE filter aux tag values from FILE"
    )?;
    writeln!(
        w,
        "  -f, --require-flags FLAG include reads with all FLAG bits set"
    )?;
    writeln!(
        w,
        "      --rf, --include-flags FLAG include reads with any FLAG bits set"
    )?;
    writeln!(
        w,
        "  -F, --exclude-flags FLAG exclude reads with any FLAG bits set"
    )?;
    writeln!(w, "  -G FLAG      exclude reads with all FLAG bits set")?;
    writeln!(w, "  --i1 FILE    write first index reads to FILE")?;
    writeln!(w, "  --i2 FILE    write second index reads to FILE")?;
    writeln!(w, "  --quality-tag TAG")?;
    writeln!(w, "               aux tag holding barcode qualities [QT]")?;
    writeln!(w, "  --index-format STR")?;
    writeln!(w, "               how to parse barcode/quality tags [i*i*]")?;
    Ok(())
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct IndexFormatItem {
    is_index: bool,
    len: Option<usize>,
}

fn parse_index_format(spec: &str) -> Result<Vec<IndexFormatItem>, String> {
    let mut items = Vec::new();
    let bytes = spec.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let kind = match bytes[i] {
            b'i' | b'I' => true,
            b'n' | b'N' => false,
            other => {
                return Err(format!(
                    "Unknown index-format code '{}' in \"{}\"",
                    char::from(other),
                    spec
                ));
            }
        };
        i += 1;

        let len = if i < bytes.len() && bytes[i].is_ascii_digit() {
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
            let digits = std::str::from_utf8(&bytes[start..i]).unwrap_or("0");
            let n = digits
                .parse::<usize>()
                .map_err(|_| format!("invalid length in index-format \"{}\"", spec))?;
            Some(n)
        } else if i < bytes.len() && bytes[i] == b'*' {
            i += 1;
            None
        } else if i >= bytes.len() {
            return Err(format!(
                "incomplete index-format \"{}\": expected length or '*'",
                spec
            ));
        } else {
            return Err(format!(
                "unexpected character '{}' in index-format \"{}\"",
                char::from(bytes[i]),
                spec
            ));
        };

        items.push(IndexFormatItem {
            is_index: kind,
            len,
        });
    }

    if items.iter().filter(|item| item.is_index).count() > 2 {
        return Err(format!(
            "Invalid index format: more than 2 indexes in \"{}\"",
            spec
        ));
    }

    Ok(items)
}

#[derive(Clone, Copy)]
struct IndexEmitOptions<'a> {
    append_read_number: bool,
    use_original_quality: bool,
    default_quality: Option<u8>,
    umi_tags: Option<&'a [[u8; 2]]>,
    casava: bool,
    barcode_tag: [u8; 2],
    quality_tag: [u8; 2],
    index_format: &'a [IndexFormatItem],
    index_file_1: Option<&'a std::path::Path>,
    index_file_2: Option<&'a std::path::Path>,
    fasta_mode: bool,
}

fn emit_index_files(
    input: Option<&std::path::Path>,
    exact: Option<Exact>,
    stdin_input: bool,
    flag_filters: FlagFilters,
    options: IndexEmitOptions<'_>,
) -> io::Result<()> {
    let mut i1_writer = match options.index_file_1 {
        Some(path) => Some(File::create(path)?),
        None => None,
    };
    let mut i2_writer = match options.index_file_2 {
        Some(path) => Some(File::create(path)?),
        None => None,
    };

    let render = |out_i1: &mut Option<File>,
                  out_i2: &mut Option<File>,
                  record: &dyn IndexRecord|
     -> io::Result<bool> {
        emit_index_for_record(out_i1.as_mut(), out_i2.as_mut(), record, options)
    };

    // Upstream `bam_fastq.c::flush_rec` emits at most ONE index record
    // per qname group (template), not per alignment record. Records
    // arrive name-grouped; emit for the first barcode-bearing record in
    // each group and skip the rest of that group.
    let mut group_qname: Option<String> = None;
    let mut group_emitted = false;

    if stdin_input {
        // For stdin we cannot re-iterate; skip index emission.
        return Ok(());
    }
    let input =
        input.ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing input path"))?;
    let exact =
        exact.ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing input format"))?;

    match exact {
        Exact::Sam => {
            let file = File::open(input)?;
            let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(file));
            let _header = reader.read_header()?;
            for result in reader.records() {
                let record = result?;
                if record_passes_flag_filter(&record, flag_filters)?
                    && record_index_eligible(&record)?
                {
                    let wrapped = SamIndexRecord(&record);
                    let qname = wrapped.name()?;
                    if group_qname.as_deref() != Some(qname.as_str()) {
                        group_qname = Some(qname);
                        group_emitted = false;
                    }
                    if !group_emitted && render(&mut i1_writer, &mut i2_writer, &wrapped)? {
                        group_emitted = true;
                    }
                }
            }
        }
        Exact::Bam => {
            let mut reader = bam::io::Reader::new(File::open(input)?);
            let header = reader.read_header()?;
            let mut record = htslib_rs::sam::alignment::RecordBuf::default();
            loop {
                let n = reader.read_record_buf(&header, &mut record)?;
                if n == 0 {
                    break;
                }
                if record_passes_flag_filter(&record, flag_filters)?
                    && record_index_eligible(&record)?
                {
                    let wrapped = RecordBufIndexRecord(&record);
                    let qname = wrapped.name()?;
                    if group_qname.as_deref() != Some(qname.as_str()) {
                        group_qname = Some(qname);
                        group_emitted = false;
                    }
                    if !group_emitted && render(&mut i1_writer, &mut i2_writer, &wrapped)? {
                        group_emitted = true;
                    }
                }
            }
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "index FASTQ output for CRAM input is not yet supported",
            ));
        }
    }

    Ok(())
}

fn record_index_eligible<R>(record: &R) -> io::Result<bool>
where
    R: htslib_rs::sam::alignment::Record + ?Sized,
{
    let flags = record.flags()?;
    Ok(!flags.is_last_segment())
}

trait IndexRecord {
    fn name(&self) -> io::Result<String>;
    fn flags(&self) -> io::Result<htslib_rs::sam::alignment::record::Flags>;
    fn barcode(&self, tag: [u8; 2]) -> io::Result<Option<String>>;
    fn quality_tag(&self, tag: [u8; 2]) -> io::Result<Option<String>>;
    fn umi(&self, tags: &[[u8; 2]]) -> io::Result<Option<String>>;
}

struct SamIndexRecord<'a>(&'a htslib_rs::sam::record::Record);
struct RecordBufIndexRecord<'a>(&'a htslib_rs::sam::alignment::RecordBuf);

impl IndexRecord for SamIndexRecord<'_> {
    fn name(&self) -> io::Result<String> {
        fastq_record_name(self.0)
    }
    fn flags(&self) -> io::Result<htslib_rs::sam::alignment::record::Flags> {
        self.0.flags()
    }
    fn barcode(&self, tag: [u8; 2]) -> io::Result<Option<String>> {
        fastq_string_tag(self.0, tag)
    }
    fn quality_tag(&self, tag: [u8; 2]) -> io::Result<Option<String>> {
        fastq_string_tag(self.0, tag)
    }
    fn umi(&self, tags: &[[u8; 2]]) -> io::Result<Option<String>> {
        fastq_umi_string(self.0, tags)
    }
}

impl IndexRecord for RecordBufIndexRecord<'_> {
    fn name(&self) -> io::Result<String> {
        fastq_record_name(self.0)
    }
    fn flags(&self) -> io::Result<htslib_rs::sam::alignment::record::Flags> {
        htslib_rs::sam::alignment::Record::flags(self.0)
    }
    fn barcode(&self, tag: [u8; 2]) -> io::Result<Option<String>> {
        fastq_string_tag(self.0, tag)
    }
    fn quality_tag(&self, tag: [u8; 2]) -> io::Result<Option<String>> {
        fastq_string_tag(self.0, tag)
    }
    fn umi(&self, tags: &[[u8; 2]]) -> io::Result<Option<String>> {
        fastq_umi_string(self.0, tags)
    }
}

fn emit_index_for_record(
    i1_writer: Option<&mut File>,
    i2_writer: Option<&mut File>,
    record: &dyn IndexRecord,
    options: IndexEmitOptions<'_>,
) -> io::Result<bool> {
    let Some(bc) = record.barcode(options.barcode_tag)? else {
        return Ok(false);
    };

    let qt_raw = record.quality_tag(options.quality_tag)?;
    let qt = match qt_raw.as_deref() {
        Some(q) if q.len() == bc.len() => Some(q.to_string()),
        _ => None,
    };

    let bc_bytes = bc.as_bytes();
    let qt_bytes = qt.as_deref().map(str::as_bytes);

    let name = record.name()?;
    let name = append_fastq_umi_with_str(name, record, options.umi_tags)?;
    let name = append_fastq_read_number_with_record(name, record, options.append_read_number)?;
    // With `-i` (CASAVA), the index record carries the same
    // ` <rnum>:<filt>:0:<barcode>` comment as the main reads, using the
    // representative record's normalized barcode (upstream `flush_rec`).
    let name = if options.casava {
        let flags = record.flags()?;
        let read_number = if flags.is_last_segment() { 2 } else { 1 };
        let filter = if flags.is_qc_fail() { "Y" } else { "N" };
        let bc_field = casava_barcode_field(record.barcode(options.barcode_tag)?.as_deref());
        format!("{name} {read_number}:{filter}:0:{bc_field}")
    } else {
        name
    };

    let writers = [i1_writer, i2_writer];
    let mut writer_idx = 0usize;
    let mut writers_arr: [Option<&mut File>; 2] = writers;
    let mut bc_cursor = 0usize;

    for item in options.index_format {
        if bc_cursor >= bc_bytes.len() {
            break;
        }

        let segment_end = match item.len {
            Some(n) => (bc_cursor + n).min(bc_bytes.len()),
            None => {
                let mut j = bc_cursor;
                while j < bc_bytes.len() && bc_bytes[j].is_ascii_alphabetic() {
                    j += 1;
                }
                j
            }
        };
        let advance_past_sep = item.len.is_none();

        if item.is_index {
            if writer_idx >= 2 {
                break;
            }
            if let Some(out) = writers_arr[writer_idx].as_deref_mut() {
                let seq = &bc_bytes[bc_cursor..segment_end];
                let qual = qt_bytes.map(|q| &q[bc_cursor..segment_end]);
                write_index_fastq_record(
                    out,
                    &name,
                    seq,
                    qual,
                    options.fasta_mode,
                    options.default_quality,
                    options.use_original_quality,
                )?;
            }
            writer_idx += 1;
        }

        bc_cursor = segment_end + if advance_past_sep { 1 } else { 0 };
    }

    Ok(true)
}

fn write_index_fastq_record(
    writer: &mut File,
    name: &str,
    seq: &[u8],
    qual: Option<&[u8]>,
    fasta_mode: bool,
    default_quality: Option<u8>,
    _use_original_quality: bool,
) -> io::Result<()> {
    if seq.is_empty() {
        return Ok(());
    }
    if fasta_mode {
        writeln!(writer, ">{name}")?;
        writer.write_all(seq)?;
        writer.write_all(b"\n")?;
        return Ok(());
    }
    writeln!(writer, "@{name}")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    match qual {
        Some(q) if q.len() == seq.len() => writer.write_all(q)?,
        _ => {
            let fill = default_quality.unwrap_or(1) + b'!';
            for _ in 0..seq.len() {
                writer.write_all(&[fill])?;
            }
        }
    }
    writer.write_all(b"\n")?;
    Ok(())
}

fn append_fastq_umi_with_str(
    name: String,
    record: &dyn IndexRecord,
    umi_tags: Option<&[[u8; 2]]>,
) -> io::Result<String> {
    let Some(umi_tags) = umi_tags else {
        return Ok(name);
    };
    let Some(umi) = record.umi(umi_tags)? else {
        return Ok(name);
    };
    let mut name = name;
    let umi = umi
        .chars()
        .map(|c| if c.is_ascii_alphabetic() { c } else { '+' })
        .collect::<String>();
    if let Some(hash) = name.rfind('#') {
        name.insert_str(hash, &format!(":{umi}"));
    } else {
        name.push(':');
        name.push_str(&umi);
    }
    Ok(name)
}

fn append_fastq_read_number_with_record(
    mut name: String,
    record: &dyn IndexRecord,
    append: bool,
) -> io::Result<String> {
    if append {
        let flags = record.flags()?;
        if flags.is_first_segment() {
            name.push_str("/1");
        } else if flags.is_last_segment() {
            name.push_str("/2");
        }
    }
    Ok(name)
}
