//! `samtools import` — convert FASTA / FASTQ records to SAM/BAM/CRAM.
//!
//! Mirrors `main_import` in `bam_import.c`. The initial Rust port supports
//! single FASTA/FASTQ and paired FASTQ → SAM/BAM (stdout or `-o FILE`). All emitted
//! SAM records are unmapped, with paired FASTQ records using the standard
//! unmapped read1/read2 flags.
//!
//! **Not yet supported:** paired singleton/other grouping (`-0` with paired
//! inputs), full read group validation, CRAM output.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use flate2::read::MultiGzDecoder;

use crate::diagnostics::{print_error, print_error_errno};

/// Entry point for `samtools import`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut output: Option<PathBuf> = None;
    let mut output_format = OutputFormat::Sam;
    let mut single_input: Option<PathBuf> = None;
    let mut interleaved_input: Option<PathBuf> = None;
    let mut read1_input: Option<PathBuf> = None;
    let mut read2_input: Option<PathBuf> = None;
    let mut index1_input: Option<PathBuf> = None;
    let mut index2_input: Option<PathBuf> = None;
    let mut index_on_both_reads = false;
    let mut positional_inputs: Vec<PathBuf> = Vec::new();
    let mut fastx_options = htslib_rs::fastq_compat::FastxToSamOptions::default();
    let mut read_group_header: Option<String> = None;
    let mut read_group_id_arg: Option<String> = None;
    let mut read_group_line: Option<String> = None;
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-o" => {
                output = iter.next().map(PathBuf::from);
            }
            "-O" | "--output-fmt" => {
                let Some(value) = iter.next() else {
                    print_error("import", "missing value for -O");
                    return ExitCode::from(1);
                };
                match OutputFormat::parse(&value.to_string_lossy()) {
                    Some(format) => output_format = format,
                    None => {
                        print_error(
                            "import",
                            format!("unsupported --output-fmt \"{}\"", value.to_string_lossy()),
                        );
                        return ExitCode::from(1);
                    }
                }
            }
            "--sam" => {
                output_format = OutputFormat::Sam;
            }
            "--bam" => {
                output_format = OutputFormat::Bam;
            }
            "--cram" => {
                output_format = OutputFormat::Cram;
            }
            "-1" => {
                read1_input = iter.next().map(PathBuf::from);
            }
            "--r1" => {
                read1_input = iter.next().map(PathBuf::from);
            }
            "-2" => {
                read2_input = iter.next().map(PathBuf::from);
            }
            "--r2" => {
                read2_input = iter.next().map(PathBuf::from);
            }
            "--i1" => {
                index1_input = iter.next().map(PathBuf::from);
            }
            "--i2" => {
                index2_input = iter.next().map(PathBuf::from);
            }
            "-s" => {
                interleaved_input = iter.next().map(PathBuf::from);
            }
            "-0" => {
                single_input = iter.next().map(PathBuf::from);
            }
            "-b" => {
                index_on_both_reads = true;
            }
            "-i" => {
                fastx_options.casava = true;
            }
            "-N" | "--name2" => {
                fastx_options.name2 = true;
            }
            "-U" | "--umi" | "--UMI" => {
                fastx_options.umi_tag = Some(String::from("RX"));
            }
            "--UMI-tag" | "--umi-tag" => {
                fastx_options.umi_tag = iter
                    .next()
                    .map(|value| value.to_string_lossy().into_owned());
            }
            "--barcode-tag" => {
                fastx_options.barcode_tag = iter
                    .next()
                    .map(|value| value.to_string_lossy().into_owned());
            }
            "--quality-tag" => {
                fastx_options.barcode_quality_tag = iter
                    .next()
                    .map(|value| value.to_string_lossy().into_owned());
            }
            "-T" => {
                if let Some(tags) = iter.next() {
                    configure_aux_tags(&mut fastx_options, &tags.to_string_lossy());
                }
            }
            "-R" => {
                if let Some(id) = iter.next() {
                    read_group_id_arg = Some(id.to_string_lossy().into_owned());
                }
            }
            "-r" => {
                if let Some(spec) = iter.next() {
                    append_read_group_spec(&mut read_group_line, &spec.to_string_lossy());
                }
            }
            "--no-PG" => {}
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error(
                    "import",
                    format!("option `{}` is not yet supported in samtools-rs import", s),
                );
                return ExitCode::from(1);
            }
            _ => {
                positional_inputs.push(PathBuf::from(arg));
            }
        }
    }

    match finalize_read_group(read_group_line, read_group_id_arg) {
        Ok(Some((header, id))) => {
            fastx_options.read_group_id = Some(id);
            read_group_header = Some(header);
        }
        Ok(None) => {}
        Err(e) => {
            print_error("import", e);
            return ExitCode::from(1);
        }
    }

    if output_format == OutputFormat::Cram {
        print_error(
            "import",
            "CRAM output is not yet supported in samtools-rs import",
        );
        return ExitCode::from(1);
    }

    if single_input.is_none()
        && interleaved_input.is_none()
        && read1_input.is_none()
        && read2_input.is_none()
    {
        match positional_inputs.as_slice() {
            [input] => match positional_fastq_looks_interleaved(input) {
                Ok(true) => interleaved_input = Some(input.clone()),
                Ok(false) => single_input = Some(input.clone()),
                Err(e) => {
                    print_error_errno("import", format!("inspect \"{}\"", input.display()), &e);
                    return ExitCode::from(1);
                }
            },
            [read1, read2] => {
                read1_input = Some(read1.clone());
                read2_input = Some(read2.clone());
            }
            [] => {}
            _ => {
                print_error(
                    "import",
                    "expected one single-end input or two paired inputs",
                );
                return ExitCode::from(1);
            }
        }
    } else if !positional_inputs.is_empty() {
        print_error(
            "import",
            "positional inputs cannot be combined with explicit input options",
        );
        return ExitCode::from(1);
    }

    let explicit_input_count = usize::from(single_input.is_some())
        + usize::from(interleaved_input.is_some())
        + usize::from(read1_input.is_some() || read2_input.is_some());
    if explicit_input_count > 1 {
        print_error(
            "import",
            "single-end, interleaved, and paired -1/-2 inputs are mutually exclusive",
        );
        return ExitCode::from(1);
    }

    let paired_inputs = match (read1_input, read2_input) {
        (Some(read1), Some(read2)) => Some((read1, read2)),
        (Some(_), None) | (None, Some(_)) => {
            print_error("import", "paired FASTQ import requires both -1 and -2");
            return ExitCode::from(1);
        }
        (None, None) => None,
    };

    let mut out: Box<dyn Write> = match output.as_ref() {
        Some(p) => match File::create(p) {
            Ok(f) => Box::new(f),
            Err(e) => {
                print_error_errno("import", "open -o output", &e);
                return ExitCode::from(1);
            }
        },
        None => Box::new(io::stdout().lock()),
    };

    if let Some(interleaved) = interleaved_input {
        let inputs = SingleFastxInputs {
            input: &interleaved,
            index1: index1_input.as_deref(),
            index2: index2_input.as_deref(),
        };
        if let Err(e) = stream_fastx_as_sam(
            inputs,
            &mut out,
            output_format,
            &fastx_options,
            read_group_header.as_deref(),
            Some(reverse_comment_for_single(
                " -n -o paired.fastq",
                index1_input.is_some(),
                index2_input.is_some(),
                &fastx_options,
            )),
        ) {
            print_error_errno(
                "import",
                format!("import \"{}\"", interleaved.display()),
                &e,
            );
            return ExitCode::from(1);
        }
    } else if let Some((read1, read2)) = paired_inputs {
        let inputs = PairedFastqInputs {
            read1: &read1,
            read2: &read2,
            index1: index1_input.as_deref(),
            index2: index2_input.as_deref(),
            index_on_both_reads,
        };
        if let Err(e) = stream_paired_fastq_as_sam(
            inputs,
            &mut out,
            output_format,
            &fastx_options,
            read_group_header.as_deref(),
        ) {
            print_error_errno(
                "import",
                format!("import \"{}\" \"{}\"", read1.display(), read2.display()),
                &e,
            );
            return ExitCode::from(1);
        }
    } else {
        let Some(input) = single_input else {
            let _ = print_usage();
            return ExitCode::from(1);
        };

        let inputs = SingleFastxInputs {
            input: &input,
            index1: index1_input.as_deref(),
            index2: index2_input.as_deref(),
        };
        let reverse_comment = if index1_input.is_some() || index2_input.is_some() {
            Some(reverse_comment_for_single(
                " -0 unpaired.fastq",
                index1_input.is_some(),
                index2_input.is_some(),
                &fastx_options,
            ))
        } else if fastx_options.umi_tag.is_some() {
            Some(reverse_comment_for_single(
                " -n -o paired.fastq",
                false,
                false,
                &fastx_options,
            ))
        } else {
            read_group_header
                .as_ref()
                .map(|_| "Reverse with: samtools fastq -0 single.fastq".to_string())
        };

        if let Err(e) = stream_fastx_as_sam(
            inputs,
            &mut out,
            output_format,
            &fastx_options,
            read_group_header.as_deref(),
            reverse_comment,
        ) {
            print_error_errno("import", format!("import \"{}\"", input.display()), &e);
            return ExitCode::from(1);
        }
    }
    ExitCode::SUCCESS
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum OutputFormat {
    Sam,
    Bam,
    Cram,
}

impl OutputFormat {
    fn parse(value: &str) -> Option<Self> {
        match value
            .split(',')
            .next()
            .unwrap_or("")
            .to_ascii_lowercase()
            .as_str()
        {
            "sam" => Some(Self::Sam),
            "bam" => Some(Self::Bam),
            "cram" => Some(Self::Cram),
            _ => None,
        }
    }
}

fn append_read_group_spec(line: &mut Option<String>, spec: &str) {
    let spec = spec.replace("\\t", "\t");
    match line {
        Some(line) => {
            line.push('\t');
            line.push_str(&spec);
        }
        None if spec.starts_with("@RG") => *line = Some(spec),
        None => *line = Some(format!("@RG\t{spec}")),
    }
}

fn finalize_read_group(
    line: Option<String>,
    id_arg: Option<String>,
) -> Result<Option<(String, String)>, &'static str> {
    if let Some(header) = line {
        return read_group_from_header(header).map(Some);
    }

    Ok(id_arg.map(|id| (format!("@RG\tID:{id}"), id)))
}

fn read_group_from_header(header: String) -> Result<(String, String), &'static str> {
    let id = header
        .split('\t')
        .find_map(|field| field.strip_prefix("ID:"))
        .map(String::from);

    match id {
        Some(id) => Ok((header, id)),
        None => Err("\"-r RG-LINE\" option contained no ID field"),
    }
}

fn configure_aux_tags(options: &mut htslib_rs::fastq_compat::FastxToSamOptions, tag_list: &str) {
    options.include_aux = true;
    if tag_list.is_empty() || tag_list == "*" {
        options.aux_tags = None;
    } else {
        options.aux_tags = Some(
            tag_list
                .split(',')
                .filter(|tag| tag.len() == 2)
                .map(String::from)
                .collect(),
        );
    }
}

struct SingleFastxInputs<'a> {
    input: &'a std::path::Path,
    index1: Option<&'a std::path::Path>,
    index2: Option<&'a std::path::Path>,
}

fn stream_fastx_as_sam(
    inputs: SingleFastxInputs<'_>,
    out: &mut dyn Write,
    output_format: OutputFormat,
    options: &htslib_rs::fastq_compat::FastxToSamOptions,
    read_group_header: Option<&str>,
    reverse_comment: Option<String>,
) -> io::Result<()> {
    let mut reader = open_text_reader(inputs.input)?;
    let first = first_non_whitespace(&mut *reader)?;
    let mut buf = Vec::new();

    match first {
        Some(b'>') => htslib_rs::fastq_compat::write_sam_from_fasta(reader, &mut buf, options)?,
        Some(b'@') if inputs.index1.is_some() || inputs.index2.is_some() => {
            let index1_reader = inputs.index1.map(open_text_reader).transpose()?;
            let index2_reader = inputs.index2.map(open_text_reader).transpose()?;
            htslib_rs::fastq_compat::write_sam_from_fastq_with_indexes(
                reader,
                index1_reader,
                index2_reader,
                &mut buf,
                options,
            )?;
        }
        Some(b'@') => htslib_rs::fastq_compat::write_sam_from_fastq(reader, &mut buf, options)?,
        Some(c) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected FASTA/FASTQ definition, got {:?}", char::from(c)),
            ));
        }
        None => {}
    }

    write_import_output(
        out,
        output_format,
        &buf,
        read_group_header,
        reverse_comment.as_deref(),
    )
}

struct PairedFastqInputs<'a> {
    read1: &'a std::path::Path,
    read2: &'a std::path::Path,
    index1: Option<&'a std::path::Path>,
    index2: Option<&'a std::path::Path>,
    index_on_both_reads: bool,
}

fn stream_paired_fastq_as_sam(
    inputs: PairedFastqInputs<'_>,
    out: &mut dyn Write,
    output_format: OutputFormat,
    options: &htslib_rs::fastq_compat::FastxToSamOptions,
    read_group_header: Option<&str>,
) -> io::Result<()> {
    let read1_reader = open_text_reader(inputs.read1)?;
    let read2_reader = open_text_reader(inputs.read2)?;
    let index1_reader = inputs.index1.map(open_text_reader).transpose()?;
    let index2_reader = inputs.index2.map(open_text_reader).transpose()?;
    let mut buf = Vec::new();

    if index1_reader.is_some() || index2_reader.is_some() {
        htslib_rs::fastq_compat::write_sam_from_paired_fastq_with_indexes(
            read1_reader,
            read2_reader,
            index1_reader,
            index2_reader,
            &mut buf,
            options,
            inputs.index_on_both_reads,
        )?;
    } else {
        htslib_rs::fastq_compat::write_sam_from_paired_fastq(
            read1_reader,
            read2_reader,
            &mut buf,
            options,
        )?;
    }

    let reverse_comment = if inputs.index1.is_some() || inputs.index2.is_some() {
        let mut comment = String::from("Reverse with: samtools fastq");
        let mut index_format = String::new();
        if inputs.index1.is_some() {
            comment.push_str(" --i1 I1.fastq");
            index_format.push_str("i*");
        }
        if inputs.index2.is_some() {
            comment.push_str(" --i2 I2.fastq");
            index_format.push_str("i*");
        }
        comment.push_str(" -1 R1.fastq -2 R2.fastq");
        if !index_format.is_empty() {
            comment.push_str(" --index-format=\"");
            comment.push_str(&index_format);
            comment.push('"');
        }
        Some(comment)
    } else {
        read_group_header
            .map(|_| String::from("Reverse with: samtools fastq -1 R1.fastq -2 R2.fastq"))
    };

    write_import_output(
        out,
        output_format,
        &buf,
        read_group_header,
        reverse_comment.as_deref(),
    )
}

fn write_import_output(
    out: &mut dyn Write,
    output_format: OutputFormat,
    body: &[u8],
    read_group_header: Option<&str>,
    reverse_comment: Option<&str>,
) -> io::Result<()> {
    let mut sam = Vec::new();
    if read_group_header.is_some() || reverse_comment.is_some() {
        writeln!(sam, "@HD\tVN:1.6\tSO:unsorted\tGO:query")?;
        if let Some(reverse_comment) = reverse_comment {
            writeln!(sam, "@CO\t{reverse_comment}")?;
        }
        if let Some(read_group_header) = read_group_header {
            writeln!(sam, "{read_group_header}")?;
        }
    }
    sam.extend_from_slice(body);

    match output_format {
        OutputFormat::Sam => out.write_all(&sam),
        OutputFormat::Bam => {
            htslib_rs::alignment_compat::write_bam_from_sam_reader(Cursor::new(sam), out)?;
            Ok(())
        }
        OutputFormat::Cram => Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "CRAM output is not yet supported",
        )),
    }
}

fn reverse_comment_for_single(
    read_part: &str,
    has_index1: bool,
    has_index2: bool,
    options: &htslib_rs::fastq_compat::FastxToSamOptions,
) -> String {
    let mut comment = String::from("Reverse with: samtools fastq");
    let mut index_format = String::new();
    if has_index1 {
        comment.push_str(" --i1 I1.fastq");
        index_format.push_str("i*");
    }
    if has_index2 {
        comment.push_str(" --i2 I2.fastq");
        index_format.push_str("i*");
    }
    comment.push_str(read_part);
    if options.casava {
        comment.push_str(" -i");
        index_format.push_str("i*i*");
    }
    if let Some(tag) = options.umi_tag.as_deref() {
        comment.push_str(" -U --UMI-tag ");
        comment.push_str(tag);
    }
    if !index_format.is_empty() {
        if options.casava && !has_index1 && !has_index2 {
            comment.push_str(" --index-format '");
            comment.push_str(&index_format);
            comment.push('\'');
        } else {
            comment.push_str(" --index-format=\"");
            comment.push_str(&index_format);
            comment.push('"');
        }
    }
    comment
}

fn open_text_reader(input: &std::path::Path) -> io::Result<Box<dyn BufRead>> {
    let file = File::open(input)?;
    // Detect gzip magic and route accordingly.
    let mut probe = [0u8; 2];
    let mut probe_file = File::open(input)?;
    let n = std::io::Read::read(&mut probe_file, &mut probe)?;
    drop(probe_file);
    if n >= 2 && probe[0] == 0x1f && probe[1] == 0x8b {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn first_non_whitespace(reader: &mut dyn BufRead) -> io::Result<Option<u8>> {
    loop {
        let available = reader.fill_buf()?;
        if available.is_empty() {
            return Ok(None);
        }
        if let Some(byte) = available.iter().copied().find(|b| !b.is_ascii_whitespace()) {
            return Ok(Some(byte));
        }
        let len = available.len();
        reader.consume(len);
    }
}

fn positional_fastq_looks_interleaved(input: &std::path::Path) -> io::Result<bool> {
    let mut reader = open_text_reader(input)?;
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Ok(false);
        }

        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }

        let Some(definition) = line.strip_prefix('@') else {
            return Ok(false);
        };
        let name = definition.split_whitespace().next().unwrap_or("");
        return Ok(name.ends_with("/1") || name.ends_with("/2"));
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools import [options] <in.fa>|<in.fq>")?;
    writeln!(
        w,
        "  -o FILE      write SAM output to FILE (default stdout)"
    )?;
    writeln!(w, "  -O sam|bam   output format [sam]")?;
    writeln!(
        w,
        "  -s FILE      single-end input (alias for positional arg)"
    )?;
    writeln!(
        w,
        "  -0 FILE      single-end input (accepted for upstream compatibility)"
    )?;
    writeln!(w, "  -1 FILE      paired FASTQ read 1 input")?;
    writeln!(w, "  -2 FILE      paired FASTQ read 2 input")?;
    writeln!(w, "  --i1 FILE    index FASTQ read 1 input")?;
    writeln!(w, "  --i2 FILE    index FASTQ read 2 input")?;
    writeln!(w, "  -i           parse CASAVA identifier")?;
    writeln!(w, "  -U, --umi    parse UMI from read name")?;
    writeln!(w, "  --UMI-tag TAG")?;
    writeln!(w, "               tag to use for UMI sequences [RX]")?;
    writeln!(w, "  --barcode-tag TAG")?;
    writeln!(
        w,
        "               tag to use with CASAVA barcode sequences [BC]"
    )?;
    writeln!(w, "  --quality-tag TAG")?;
    writeln!(w, "               tag to use with barcode qualities [QT]")?;
    writeln!(w, "  -N, --name2  use 2nd field as read name")?;
    writeln!(w, "  -T TAGLIST   preserve FASTQ definition aux tags")?;
    writeln!(w, "  -R STRING    add @RG ID:STRING and RG tags")?;
    writeln!(w, "  -r STRING    add complete or partial @RG line")?;
    Ok(())
}
