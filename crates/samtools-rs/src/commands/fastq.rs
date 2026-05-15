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
//! **Not yet supported:** CRAM input, exact `score`-vs-`b_score` propagation
//! through CASAVA barcode copying between paired ends.

use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use htslib_rs::{bam, format::Exact};

use crate::aux_list::parse_aux_list;
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

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
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct FlagFilters {
    require: u16,
    include_any: u16,
    exclude: u16,
    exclude_all: u16,
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
                match raw {
                    "" | "*" => aux_selection = AuxSelection::All,
                    _ => match parse_aux_list(raw) {
                        Ok(tags) => {
                            let tags: Vec<[u8; 2]> = tags.into_iter().collect();
                            merge_aux_selection(&mut aux_selection, &tags);
                        }
                        Err(e) => {
                            print_error(sub_name, format!("invalid -T value \"{raw}\": {e}"));
                            return ExitCode::from(1);
                        }
                    },
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
                    Ok(filter) => tag_filters.push(filter),
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
                    Ok(filter) => tag_filters.push(filter),
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

    let stdin_input = input.as_ref().is_none_or(|path| path.as_os_str() == "-");

    let format = if stdin_input {
        None
    } else {
        let input = input.as_ref().expect("non-stdin input exists");
        match sam_io::sam_open_format(input) {
            Ok(f) => Some(f),
            Err(e) => {
                print_error(sub_name, e.to_string());
                return ExitCode::from(1);
            }
        }
    };

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
    let append_read_number = append_read_number_override.unwrap_or(!split_mode || singleton_only);
    let render_options = FastqRenderOptions {
        append_read_number,
        use_original_quality,
        default_quality,
        umi_tags: umi_tags.as_deref(),
        casava,
        barcode_tag,
        aux_selection: &aux_selection,
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
            if fasta_mode && (flag_filters.include_any != 0 || umi_tags.is_some() || casava) {
                view_sam_reader_as_fasta_split(
                    &mut reader,
                    flag_filters,
                    append_read_number,
                    umi_tags.as_deref(),
                    casava,
                    barcode_tag,
                    singleton_set,
                )
            } else if fasta_mode {
                htslib_rs::alignment_compat::view_sam_as_fasta_split_text_from_reader_with_flag_filter_and_suffix(
                    &mut reader,
                    require_flags,
                    exclude_flags,
                    exclude_all_flags,
                    append_read_number,
                )
                .map(FastqSplitBuffers::from_fast_path)
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
            match (format.expect("non-stdin format exists").exact, fasta_mode) {
                (Exact::Sam, false) => view_sam_path_as_fastq_split(
                    input,
                    flag_filters,
                    render_options,
                    &tag_filters,
                    singleton_set,
                ),
                (Exact::Sam, true) => {
                    if flag_filters.include_any != 0 || umi_tags.is_some() || casava {
                        view_sam_path_as_fasta_split(
                            input,
                            flag_filters,
                            append_read_number,
                            umi_tags.as_deref(),
                            casava,
                            barcode_tag,
                            singleton_set,
                        )
                    } else {
                        htslib_rs::alignment_compat::view_sam_as_fasta_split_text_from_path_with_flag_filter_and_suffix(
                            input,
                            require_flags,
                            exclude_flags,
                            exclude_all_flags,
                            append_read_number,
                        )
                        .map(FastqSplitBuffers::from_fast_path)
                    }
                }
                (Exact::Bam, false) => view_bam_path_as_fastq_split_with_aux(
                    input,
                    flag_filters,
                    render_options,
                    &tag_filters,
                    singleton_set,
                ),
                (Exact::Bam, true) => {
                    if flag_filters.include_any != 0 || umi_tags.is_some() || casava {
                        view_bam_path_as_fasta_split(
                            input,
                            flag_filters,
                            append_read_number,
                            umi_tags.as_deref(),
                            casava,
                            barcode_tag,
                            singleton_set,
                        )
                    } else {
                        htslib_rs::alignment_compat::view_bam_as_fasta_split_text_from_path_with_flag_filter_and_suffix(
                            input,
                            require_flags,
                            exclude_flags,
                            exclude_all_flags,
                            append_read_number,
                        )
                        .map(FastqSplitBuffers::from_fast_path)
                    }
                }
                _ => {
                    print_error(
                        sub_name,
                        "only SAM and BAM input are currently supported (CRAM TODO)",
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
        if !fasta_mode
            && (use_original_quality
                || default_quality.is_some()
                || umi_tags.is_some()
                || casava
                || aux_selection.is_enabled()
                || !tag_filters.is_empty())
        {
            view_sam_reader_as_fastq_text_with_aux(
                &mut reader,
                flag_filters,
                render_options,
                &tag_filters,
            )
        } else if fasta_mode && (flag_filters.include_any != 0 || umi_tags.is_some() || casava) {
            view_sam_reader_as_fasta_text(
                &mut reader,
                flag_filters,
                append_read_number,
                umi_tags.as_deref(),
                casava,
                barcode_tag,
            )
        } else if fasta_mode {
            htslib_rs::alignment_compat::view_sam_as_fasta_text_from_reader_with_flag_filter_and_suffix(
                &mut reader,
                require_flags,
                exclude_flags,
                exclude_all_flags,
                append_read_number,
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
            format.expect("non-stdin format exists").exact,
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
        (Exact::Sam, false, _, true) => view_sam_path_as_fastq_text_with_aux(
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
        (Exact::Sam, true, _, true) | (Exact::Bam, true, _, true) => {
            print_error(
                sub_name,
                "-d/-D tag filtering is currently supported for FASTQ single-output mode only",
            );
            return ExitCode::from(1);
        }
        (Exact::Sam, false, false, false) => {
            htslib_rs::alignment_compat::view_sam_as_fastq_text_from_path_with_limit_and_suffix(
                input,
                None,
                append_read_number,
            )
        }
        (Exact::Sam, false, true, false) => {
            htslib_rs::alignment_compat::view_sam_as_fastq_text_from_path_with_flag_filter_and_suffix(
                input,
                require_flags,
                exclude_flags,
                exclude_all_flags,
                append_read_number,
            )
        }
        (Exact::Sam, true, false, _) => {
            htslib_rs::alignment_compat::view_sam_as_fasta_text_from_path_with_limit_and_suffix(
                input,
                None,
                append_read_number,
            )
        }
        (Exact::Sam, true, true, _) => {
            if flag_filters.include_any != 0 {
                view_sam_path_as_fasta_text(
                    input,
                    flag_filters,
                    append_read_number,
                    None,
                    false,
                    *b"BC",
                )
            } else {
                htslib_rs::alignment_compat::view_sam_as_fasta_text_from_path_with_flag_filter_and_suffix(
                input,
                require_flags,
                exclude_flags,
                exclude_all_flags,
                append_read_number,
            )
            }
        }
        (Exact::Bam, false, false, false) => {
            htslib_rs::alignment_compat::view_bam_as_fastq_text_from_path_with_limit_and_suffix(
                input,
                None,
                append_read_number,
            )
        }
        (Exact::Bam, false, true, false) => {
            htslib_rs::alignment_compat::view_bam_as_fastq_text_from_path_with_flag_filter_and_suffix(
                input,
                require_flags,
                exclude_flags,
                exclude_all_flags,
                append_read_number,
            )
        }
        (Exact::Bam, false, _, true) => view_bam_path_as_fastq_text_with_aux(
            input,
            flag_filters,
            render_options,
            &tag_filters,
        ),
        (Exact::Bam, true, false, _) => {
            htslib_rs::alignment_compat::view_bam_as_fasta_text_from_path_with_limit_and_suffix(
                input,
                None,
                append_read_number,
            )
        }
        (Exact::Bam, true, true, _) => {
            if flag_filters.include_any != 0 {
                view_bam_path_as_fasta_text(
                    input,
                    flag_filters,
                    append_read_number,
                    None,
                    false,
                    *b"BC",
                )
            } else {
                htslib_rs::alignment_compat::view_bam_as_fasta_text_from_path_with_flag_filter_and_suffix(
                input,
                require_flags,
                exclude_flags,
                exclude_all_flags,
                append_read_number,
            )
            }
        }
        _ => {
            print_error(
                sub_name,
                "only SAM and BAM input are currently supported (CRAM TODO)",
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
        Ok(()) => ExitCode::SUCCESS,
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno(sub_name, "close output", &e);
            ExitCode::from(1)
        }
    }
}

fn write_text_file(path: &std::path::Path, text: &[u8]) -> io::Result<()> {
    let mut file = File::create(path)?;
    file.write_all(text)
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
    let mut grouper = GroupedSplitWriter::new(singleton_set);
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
    let mut grouper = GroupedSplitWriter::new(singleton_set);
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
    let mut grouper = GroupedSplitWriter::new(singleton_set);

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
    let mut grouper = GroupedSplitWriter::new(singleton_set);

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

impl FastqSplitBuffers {
    fn from_fast_path(text: htslib_rs::alignment_compat::FastxSplitText) -> Self {
        Self {
            read1: text.read1.into_bytes(),
            read2: text.read2.into_bytes(),
            singleton: text.singleton.into_bytes(),
            other: Vec::new(),
        }
    }
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
    current_qname: Option<Vec<u8>>,
    best_score: [u8; 3],
    pending_text: [Option<Vec<u8>>; 3],
}

impl GroupedSplitWriter {
    fn new(singleton_set: bool) -> Self {
        Self {
            singleton_set,
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
        let [t0, t1, t2] = std::mem::take(&mut self.pending_text);

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
    let sequence = fastq_sequence_string(record);
    let quality = fastq_quality_scores_string(
        record,
        options.use_original_quality,
        options.default_quality,
    )?;

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
        for field in fastq_aux_fields(record, options.aux_selection)? {
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
    let barcode = fastq_string_tag(record, barcode_tag)?.unwrap_or_else(|| "0".to_string());

    Ok(format!("{read_number}:{filter}:0:{barcode}"))
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
    if use_original_quality && let Some(mut oq) = original_quality_string(record)? {
        if record.flags()?.is_reverse_complemented() {
            oq = oq.chars().rev().collect();
        }
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
    let mut bytes = scores
        .into_iter()
        .map(|score| {
            score.checked_add(b'!').ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "FASTQ quality score overflow")
            })
        })
        .collect::<io::Result<Vec<_>>>()?;

    if record.flags()?.is_reverse_complemented() {
        bytes.reverse();
    }

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

fn format_aux_float(n: f32) -> String {
    let abs = n.abs();
    if n != 0.0 && !(1e-4..1e6).contains(&abs) {
        format_htslib_exponent(n)
    } else {
        format!("{n}")
    }
}

fn format_htslib_exponent(n: f32) -> String {
    let raw = format!("{n:e}");
    let Some((mantissa, exponent)) = raw.split_once('e') else {
        return raw;
    };
    let value = exponent.parse::<i32>().unwrap_or(0);
    format!("{mantissa}e{value:+03}")
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
    Ok(())
}
