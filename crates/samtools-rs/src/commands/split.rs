//! `samtools split` — split an alignment file by `@RG` ID or aux tag value.
//!
//! Mirrors `main_split` in `bam_split.c`. The initial Rust port supports:
//!  - Input BAM, SAM, and whole-file CRAM.
//!  - Output template via `-f` with `%*` (input basename), `%!` (RG ID or
//!    tag value), `%#` (output index), and `%.` (extension).
//!  - `-u <path>` — write records with unknown/missing RG/tag to this file.
//!  - `-h <path>` — use an alternate SAM header for the unaccounted output.
//!  - `-d TAG` — split by a string or integer aux tag instead of header `@RG`.
//!  - `-M N` / `--max-split N` — cap dynamically-created `-d` outputs.
//!  - `--output-fmt sam|bam|cram` — CRAM output requires `--reference`.
//!  - `--no-PG` — suppress the default `@PG` line.
//!  - `--write-index` — build BAI indexes for BAM outputs.
//!  - `-p N` — width for `%#` padding.
//!  - `-` — read SAM or BAM from stdin.

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::cram;
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf};

use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools split`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut template = String::from("%*_%!.%.");
    let mut unaccounted: Option<PathBuf> = None;
    let mut unaccounted_header: Option<PathBuf> = None;
    let mut output_fmt: Option<OutFmt> = None;
    let mut pad_width: usize = 0;
    let mut split_tag: Option<[u8; 2]> = None;
    let mut max_split: usize = 100;
    let mut add_pg = true;
    let mut local_write_index = false;
    let mut input: Option<PathBuf> = None;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-f" => {
                template = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .map(|s| s.to_string())
                    .unwrap_or(template);
            }
            "-u" => {
                unaccounted = iter.next().map(PathBuf::from);
            }
            "-h" => {
                unaccounted_header = iter.next().map(PathBuf::from);
            }
            "-d" => {
                let Some(tag) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("split", "missing argument for -d");
                    return ExitCode::from(1);
                };
                let bytes = tag.as_bytes();
                if bytes.len() != 2 {
                    print_error("split", "TAG for -d must be exactly two characters");
                    return ExitCode::from(1);
                }
                split_tag = Some([bytes[0], bytes[1]]);
            }
            "-M" | "--max-split" => {
                let Some(value) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("split", format!("missing argument for {}", s));
                    return ExitCode::from(1);
                };
                max_split = match value.parse::<isize>() {
                    Ok(n) if n < 0 => usize::MAX,
                    Ok(n) if n > 0 => n as usize,
                    _ => {
                        print_error("split", format!("Invalid -M argument: \"{}\"", value));
                        return ExitCode::from(1);
                    }
                };
            }
            "-p" => {
                pad_width = iter
                    .next()
                    .and_then(|a| a.to_str())
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
            }
            "--output-fmt" | "-O" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("split", format!("missing argument for {}", s));
                    return ExitCode::from(1);
                };
                output_fmt = Some(match parse_output_format(v) {
                    Ok(fmt) => fmt,
                    Err(e) => {
                        print_error("split", e);
                        return ExitCode::from(1);
                    }
                });
            }
            "--no-PG" => add_pg = false,
            "--write-index" => local_write_index = true,
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error(
                    "split",
                    format!("option `{}` is not yet supported in samtools-rs split", s),
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

    let Some(input) = input else {
        let _ = print_usage();
        return ExitCode::from(1);
    };

    if input.as_os_str() != "-" {
        if !input.exists() {
            print_hts_open_missing(&input);
            print_error(
                "split",
                format!(
                    "Could not open \"{}\": No such file or directory",
                    input.display()
                ),
            );
            return ExitCode::from(1);
        }
        let format = match sam_io::sam_open_format(&input) {
            Ok(f) => f,
            Err(e) => {
                print_error("split", e.to_string());
                return ExitCode::from(1);
            }
        };
        if !matches!(format.exact, Exact::Sam | Exact::Bam | Exact::Cram) {
            print_error(
                "split",
                "only SAM, BAM, and CRAM input are currently supported",
            );
            return ExitCode::from(1);
        }
    }

    let globals = current_global_args();
    let write_index = local_write_index || globals.write_index;
    if write_index && matches!(output_fmt, Some(fmt) if !matches!(fmt, OutFmt::Bam)) {
        print_error("split", "--write-index is only supported for BAM output");
        return ExitCode::from(1);
    }
    if matches!(output_fmt, Some(OutFmt::Cram)) && globals.reference.is_none() {
        print_error(
            "split",
            "CRAM output requires a reference (use top-level --reference FILE)",
        );
        return ExitCode::from(1);
    }

    let options = SplitOptions {
        template: &template,
        unaccounted: unaccounted.as_deref(),
        unaccounted_header: unaccounted_header.as_deref(),
        fmt: output_fmt,
        reference: globals.reference.as_deref(),
        pad_width,
        split_tag,
        max_split,
        add_pg,
        write_index,
        argv: args,
    };

    match run_split(&input, options) {
        Ok(()) => ExitCode::SUCCESS,
        Err(SplitError::AlreadyReported) => ExitCode::from(1),
        Err(SplitError::Io(e)) => {
            print_error_errno("split", "split failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Copy)]
enum OutFmt {
    Sam,
    Bam,
    Cram,
}

fn parse_output_format(raw: &str) -> Result<OutFmt, String> {
    let head = raw.split(',').next().unwrap_or("").to_ascii_lowercase();
    match head.as_str() {
        "sam" => Ok(OutFmt::Sam),
        "bam" => Ok(OutFmt::Bam),
        "cram" => Ok(OutFmt::Cram),
        _ => Err(format!("unsupported output format \"{}\"", raw)),
    }
}

struct SplitOptions<'a> {
    template: &'a str,
    unaccounted: Option<&'a Path>,
    unaccounted_header: Option<&'a Path>,
    fmt: Option<OutFmt>,
    reference: Option<&'a Path>,
    pad_width: usize,
    split_tag: Option<[u8; 2]>,
    max_split: usize,
    add_pg: bool,
    write_index: bool,
    argv: &'a [OsString],
}

enum SplitError {
    Io(io::Error),
    AlreadyReported,
}

impl From<io::Error> for SplitError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

type SplitResult<T> = Result<T, SplitError>;

impl SplitOptions<'_> {
    fn render_fmt(&self) -> OutFmt {
        self.fmt.unwrap_or(OutFmt::Bam)
    }

    fn fmt_for_path(&self, path: &Path) -> io::Result<OutFmt> {
        let fmt = self.fmt.unwrap_or_else(|| infer_output_format(path));
        if self.write_index && !matches!(fmt, OutFmt::Bam) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--write-index is only supported for BAM output",
            ));
        }
        if matches!(fmt, OutFmt::Cram) && self.reference.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "CRAM output requires a reference (use top-level --reference FILE)",
            ));
        }
        Ok(fmt)
    }
}

fn infer_output_format(path: &Path) -> OutFmt {
    match path
        .extension()
        .and_then(|ext| ext.to_str())
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        Some("sam") => OutFmt::Sam,
        Some("cram") => OutFmt::Cram,
        _ => OutFmt::Bam,
    }
}

fn run_split(input: &Path, options: SplitOptions<'_>) -> SplitResult<()> {
    let (header, records) = if input.as_os_str() == "-" {
        read_stdin_records()?
    } else {
        let format = sam_io::sam_open_format(input)?;
        match format.exact {
            Exact::Sam => read_sam_records(input)?,
            Exact::Bam => read_bam_records(input)?,
            Exact::Cram => read_cram_records(input)?,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "only SAM, BAM, and CRAM input are currently supported",
                )
                .into());
            }
        }
    };

    let input_base = input
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("split");
    let initial_tags: Vec<String> = if options.split_tag.is_none() {
        header
            .read_groups()
            .iter()
            .map(|(id, _)| String::from_utf8_lossy(id).into_owned())
            .collect()
    } else {
        Vec::new()
    };

    let mut outputs: Vec<SplitOutput> = Vec::with_capacity(initial_tags.len());
    let mut id_to_index: HashMap<String, usize> = HashMap::new();
    let mut index_paths = Vec::new();
    for (idx, tag_value) in initial_tags.iter().enumerate() {
        let output_header = prepare_output_header(
            &header,
            options.split_tag,
            tag_value,
            options.add_pg,
            options.argv,
        )?;
        let path = render_template(
            options.template,
            input_base,
            tag_value,
            idx,
            options.render_fmt(),
            options.pad_width,
        )?;
        let path = PathBuf::from(path);
        let fmt = options.fmt_for_path(&path)?;
        track_index_path(&mut index_paths, options.write_index, fmt, &path);
        outputs.push(SplitOutput {
            header: output_header.clone(),
            sink: open_output(&path, fmt, &output_header, options.reference)?,
        });
        id_to_index.insert(tag_value.clone(), idx);
    }

    let override_header = match (options.unaccounted, options.unaccounted_header) {
        (Some(_), Some(path)) => Some(read_sam_header(path)?),
        _ => None,
    };
    let unaccounted_output_header = prepare_unaccounted_header(
        override_header.as_ref().unwrap_or(&header),
        options.add_pg,
        options.argv,
    )?;
    let mut unaccounted_out: Option<Box<dyn SplitSink>> = match options.unaccounted {
        Some(p) => {
            let fmt = options.fmt_for_path(p)?;
            track_index_path(&mut index_paths, options.write_index, fmt, p);
            Some(open_output(
                p,
                fmt,
                &unaccounted_output_header,
                options.reference,
            )?)
        }
        None => None,
    };

    for record in &records {
        let tag_value = record_split_tag_value(record, options.split_tag, options.pad_width);
        let target = match tag_value.as_ref() {
            Some(value) => match id_to_index.get(value).copied() {
                Some(i) => Some(SplitTarget::Output(i)),
                None if options.split_tag.is_some() && outputs.len() < options.max_split => {
                    let idx = outputs.len();
                    let output_header = prepare_output_header(
                        &header,
                        options.split_tag,
                        value,
                        options.add_pg,
                        options.argv,
                    )?;
                    let path = render_template(
                        options.template,
                        input_base,
                        value,
                        idx,
                        options.render_fmt(),
                        options.pad_width,
                    )?;
                    let path = PathBuf::from(path);
                    let fmt = options.fmt_for_path(&path)?;
                    track_index_path(&mut index_paths, options.write_index, fmt, &path);
                    outputs.push(SplitOutput {
                        header: output_header.clone(),
                        sink: open_output(&path, fmt, &output_header, options.reference)?,
                    });
                    id_to_index.insert(value.clone(), idx);
                    Some(SplitTarget::Output(idx))
                }
                None => unaccounted_out.as_ref().map(|_| SplitTarget::Unaccounted),
            },
            None => unaccounted_out.as_ref().map(|_| SplitTarget::Unaccounted),
        };
        match target {
            Some(SplitTarget::Output(i)) => {
                let output = &mut outputs[i];
                output.sink.write_record(&output.header, record)?;
            }
            Some(SplitTarget::Unaccounted) => {
                if let Some(sink) = unaccounted_out.as_deref_mut() {
                    sink.write_record(&unaccounted_output_header, record)?;
                }
            }
            None => {
                report_unaccounted_record(record, options.split_tag, tag_value.as_deref());
                return Err(SplitError::AlreadyReported);
            }
        }
    }

    for output in &mut outputs {
        output.sink.finish(&output.header)?;
    }
    if let Some(sink) = unaccounted_out.as_deref_mut() {
        sink.finish(&unaccounted_output_header)?;
    }
    drop(outputs);
    drop(unaccounted_out);
    for path in index_paths {
        write_bam_index(&path)?;
    }
    Ok(())
}

fn report_unaccounted_record(record: &RecordBuf, split_tag: Option<[u8; 2]>, value: Option<&str>) {
    let read_name = record
        .name()
        .map(|name| String::from_utf8_lossy(name.as_ref()))
        .map(std::borrow::Cow::into_owned)
        .unwrap_or_else(|| "*".to_string());

    if let Some(value) = value {
        eprintln!("Read \"{read_name}\" with unaccounted for tag \"{value}\".");
    } else {
        let tag = split_tag.unwrap_or(*b"RG");
        let tag = String::from_utf8_lossy(&tag);
        eprintln!("Read \"{read_name}\" has no {tag} tag.");
    }
}

struct SplitOutput {
    header: sam::Header,
    sink: Box<dyn SplitSink>,
}

enum SplitTarget {
    Output(usize),
    Unaccounted,
}

fn prepare_output_header(
    header: &sam::Header,
    split_tag: Option<[u8; 2]>,
    tag_value: &str,
    add_pg: bool,
    argv: &[OsString],
) -> io::Result<sam::Header> {
    let header = header_for_split_output(header, split_tag, tag_value);
    prepare_unaccounted_header(&header, add_pg, argv)
}

fn prepare_unaccounted_header(
    header: &sam::Header,
    add_pg: bool,
    argv: &[OsString],
) -> io::Result<sam::Header> {
    if add_pg {
        crate::pg::add_samtools_pg_to_header(header, argv)
    } else {
        Ok(header.clone())
    }
}

fn header_for_split_output(
    header: &sam::Header,
    split_tag: Option<[u8; 2]>,
    tag_value: &str,
) -> sam::Header {
    let tag = split_tag.unwrap_or(*b"RG");
    if tag != *b"RG" {
        return header.clone();
    }

    let mut output_header = header.clone();
    let tag_value_bytes = tag_value.as_bytes();
    let had_read_group = output_header.read_groups().contains_key(tag_value_bytes);
    output_header.read_groups_mut().retain(|id, _| {
        let id: &[u8] = id.as_ref();
        id == tag_value_bytes
    });

    if !had_read_group && split_tag.is_some() {
        output_header.read_groups_mut().insert(
            tag_value.to_string().into(),
            sam::header::record::value::Map::<sam::header::record::value::map::ReadGroup>::default(
            ),
        );
    }

    output_header
}

fn track_index_path(index_paths: &mut Vec<PathBuf>, write_index: bool, fmt: OutFmt, path: &Path) {
    if write_index && matches!(fmt, OutFmt::Bam) {
        index_paths.push(path.to_path_buf());
    }
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

fn read_bam_records(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    read_bam_records_from_reader(&mut bam::io::Reader::new(File::open(input)?))
}

fn read_bam_records_from_reader<R>(
    reader: &mut bam::io::Reader<R>,
) -> io::Result<(sam::Header, Vec<RecordBuf>)>
where
    R: Read,
{
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

fn read_sam_records(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    read_sam_records_from_reader(&mut sam::io::Reader::new(BufReader::new(File::open(
        input,
    )?)))
}

fn read_sam_records_from_reader<R>(
    reader: &mut sam::io::Reader<R>,
) -> io::Result<(sam::Header, Vec<RecordBuf>)>
where
    R: io::BufRead,
{
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

fn read_cram_records(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(input)?;
    let records = if let Some(reference) = current_global_args().reference {
        htslib_rs::alignment_compat::query_cram_records_all_from_path_with_reference(
            input, &reference,
        )?
    } else {
        htslib_rs::alignment_compat::query_cram_records_all_from_path(input)?
    };
    Ok((header, records))
}

fn read_stdin_records() -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let mut bytes = Vec::new();
    io::stdin().read_to_end(&mut bytes)?;
    if bytes.first() == Some(&b'@') {
        let text =
            String::from_utf8(bytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut reader = sam::io::Reader::new(BufReader::new(Cursor::new(text)));
        return read_sam_records_from_reader(&mut reader);
    }

    let mut reader = bam::io::Reader::new(Cursor::new(bytes));
    read_bam_records_from_reader(&mut reader)
}

fn read_sam_header(input: &Path) -> io::Result<sam::Header> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(input)?));
    reader.read_header()
}

fn record_split_tag_value(record: &RecordBuf, tag: Option<[u8; 2]>, pad: usize) -> Option<String> {
    let tag = tag.unwrap_or(*b"RG");
    record.data().get(&tag).and_then(|value| match value {
        sam::alignment::record_buf::data::field::Value::String(s) => {
            Some(String::from_utf8_lossy(s.as_ref()).into_owned())
        }
        sam::alignment::record_buf::data::field::Value::Hex(s) => {
            Some(String::from_utf8_lossy(s.as_ref()).into_owned())
        }
        sam::alignment::record_buf::data::field::Value::Int8(n) => Some(format_tag_int(*n, pad)),
        sam::alignment::record_buf::data::field::Value::UInt8(n) => Some(format_tag_int(*n, pad)),
        sam::alignment::record_buf::data::field::Value::Int16(n) => Some(format_tag_int(*n, pad)),
        sam::alignment::record_buf::data::field::Value::UInt16(n) => Some(format_tag_int(*n, pad)),
        sam::alignment::record_buf::data::field::Value::Int32(n) => Some(format_tag_int(*n, pad)),
        sam::alignment::record_buf::data::field::Value::UInt32(n) => Some(format_tag_int(*n, pad)),
        _ => None,
    })
}

fn format_tag_int<N>(value: N, pad: usize) -> String
where
    N: Into<i64>,
{
    let value = value.into();
    if pad == 0 {
        value.to_string()
    } else if value < 0 {
        format!("{:0width$}", value, width = pad + 1)
    } else {
        format!("{:0width$}", value, width = pad)
    }
}

fn render_template(
    template: &str,
    input_base: &str,
    tag_value: &str,
    idx: usize,
    fmt: OutFmt,
    pad: usize,
) -> io::Result<String> {
    let ext = match fmt {
        OutFmt::Sam => "sam",
        OutFmt::Bam => "bam",
        OutFmt::Cram => "cram",
    };
    let idx_str = if pad > 0 {
        format!("{:0width$}", idx, width = pad)
    } else {
        idx.to_string()
    };

    let mut out = String::with_capacity(template.len() + 16);
    let mut chars = template.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '%' {
            match chars.next() {
                Some('!') => out.push_str(tag_value),
                Some('#') => out.push_str(&idx_str),
                Some('.') => out.push_str(ext),
                Some('*') => out.push_str(input_base),
                Some('%') => out.push('%'),
                Some(other) => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid split output format escape `%{}`", other),
                    ));
                }
                None => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "invalid split output format: trailing `%`",
                    ));
                }
            }
        } else {
            out.push(c);
        }
    }
    Ok(out)
}

trait SplitSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;

    fn finish(&mut self, _header: &sam::Header) -> io::Result<()> {
        Ok(())
    }
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
impl SplitSink for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}

struct SamFile(File);
impl SplitSink for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        // Use the shared renderer so `f32` aux values get htslib's `%g`
        // spelling rather than noodles' plain decimals.
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}

struct CramFile(cram::io::Writer<File>);
impl SplitSink for CramFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }

    fn finish(&mut self, header: &sam::Header) -> io::Result<()> {
        self.0.try_finish(header)
    }
}

fn open_output(
    path: &Path,
    fmt: OutFmt,
    header: &sam::Header,
    reference: Option<&Path>,
) -> io::Result<Box<dyn SplitSink>> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)?;
    }
    let file = File::create(path)?;
    match fmt {
        OutFmt::Sam => {
            let mut file = file;
            crate::sam_render::write_header(&mut file, header)?;
            Ok(Box::new(SamFile(file)))
        }
        OutFmt::Bam => {
            let mut writer = bam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(BamFile(writer)))
        }
        OutFmt::Cram => {
            let reference = reference.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "CRAM output requires --reference",
                )
            })?;
            crate::reference::ensure_fai_index(reference, None)?;
            let repository =
                htslib_rs::alignment_compat::cram_reference_repository_from_fasta_path(reference)?;
            let mut writer = cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repository)
                .build_from_writer(file);
            writer.write_header(header)?;
            Ok(Box::new(CramFile(writer)))
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools split [options] <in.bam>")?;
    writeln!(w, "Options:")?;
    writeln!(
        w,
        "  -f STR     output filename template (default: %*_%!.%.)"
    )?;
    writeln!(w, "  -u FILE    write records with unknown @RG/tag to FILE")?;
    writeln!(w, "  -h FILE    use alternate SAM header for -u output")?;
    writeln!(w, "  -d TAG     split by aux TAG instead of @RG")?;
    writeln!(w, "  -M N       maximum number of outputs created by -d")?;
    writeln!(w, "  -p N       pad the %# index to N digits")?;
    writeln!(w, "  --output-fmt sam|bam|cram")?;
    writeln!(w, "  --no-PG    do not add a @PG line")?;
    writeln!(w, "  --write-index")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::BString;
    use htslib_rs::sam::alignment::record::data::field::Tag;
    use htslib_rs::sam::alignment::record_buf::data::field::Value;
    use std::io::Cursor;

    fn argv(args: &[&str]) -> Vec<OsString> {
        args.iter().map(OsString::from).collect()
    }

    fn exit_to_u8(code: ExitCode) -> u8 {
        format!("{:?}", code)
            .chars()
            .filter(|c| c.is_ascii_digit())
            .collect::<String>()
            .parse()
            .unwrap_or(255)
    }

    fn parse_header(text: &str) -> sam::Header {
        let mut reader = sam::io::Reader::new(BufReader::new(Cursor::new(text.as_bytes())));
        reader.read_header().unwrap()
    }

    fn string_tag_record(tag: [u8; 2], value: &str) -> RecordBuf {
        let mut record = RecordBuf::default();
        record
            .data_mut()
            .insert(Tag::from(tag), Value::String(BString::from(value)));
        record
    }

    #[test]
    fn render_template_matches_upstream_split_format_cases() {
        assert_eq!(
            render_template("%*_%#.%.", "basename", "1#2.3", 4, OutFmt::Bam, 0).unwrap(),
            "basename_4.bam"
        );
        assert_eq!(
            render_template("%*_%!.%.", "basename", "1#2.3", 4, OutFmt::Bam, 0).unwrap(),
            "basename_1#2.3.bam"
        );
        assert_eq!(
            render_template("%*_%#.%.", "basename", "1#2.3", 4, OutFmt::Sam, 0).unwrap(),
            "basename_4.sam"
        );
        assert_eq!(
            render_template("%*_%#.%.", "basename", "1#2.3", 4, OutFmt::Cram, 0).unwrap(),
            "basename_4.cram"
        );
        assert_eq!(
            render_template("%*_%#.%.", "basename", "1#2.3", 4, OutFmt::Bam, 5).unwrap(),
            "basename_00004.bam"
        );
        assert_eq!(
            render_template("%%%*_%#.%.%%", "basename", "1#2.3", 4, OutFmt::Bam, 0).unwrap(),
            "%basename_4.bam%"
        );
    }

    #[test]
    fn render_template_rejects_bad_percent_escapes() {
        assert!(render_template("%%%*_%#.%.%", "basename", "1#2.3", 4, OutFmt::Bam, 0).is_err());
        assert!(render_template("%s_%#.%.", "basename", "1#2.3", 4, OutFmt::Bam, 0).is_err());
    }

    #[test]
    fn header_for_split_output_filters_read_groups() {
        let header = parse_header("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:8\n@RG\tID:grp1\n@RG\tID:grp2\n");

        let filtered = header_for_split_output(&header, None, "grp1");
        assert!(filtered.read_groups().contains_key("grp1".as_bytes()));
        assert!(!filtered.read_groups().contains_key("grp2".as_bytes()));

        let synthesized = header_for_split_output(&header, Some(*b"RG"), "grp3");
        assert!(synthesized.read_groups().contains_key("grp3".as_bytes()));
        assert!(!synthesized.read_groups().contains_key("grp1".as_bytes()));

        let non_rg = header_for_split_output(&header, Some(*b"an"), "aardvark");
        assert_eq!(non_rg.read_groups().len(), 2);
    }

    #[test]
    fn record_split_tag_value_formats_string_and_integer_tags() {
        let string_record = string_tag_record(*b"an", "aardvark");
        assert_eq!(
            record_split_tag_value(&string_record, Some(*b"an"), 0).as_deref(),
            Some("aardvark")
        );

        let mut int_record = RecordBuf::default();
        int_record
            .data_mut()
            .insert(Tag::from(*b"nn"), Value::Int8(-2));
        assert_eq!(
            record_split_tag_value(&int_record, Some(*b"nn"), 2).as_deref(),
            Some("-02")
        );
    }

    #[test]
    fn split_parse_args_handles_help_and_early_errors() {
        assert_eq!(exit_to_u8(main(&argv(&["split", "--help"]))), 0);
        assert_eq!(exit_to_u8(main(&argv(&["split"]))), 1);
        assert_eq!(exit_to_u8(main(&argv(&["split", "-d", "R", "in.sam"]))), 1);
        assert_eq!(exit_to_u8(main(&argv(&["split", "-M", "0", "in.sam"]))), 1);
        assert_eq!(
            exit_to_u8(main(&argv(&["split", "--output-fmt", "vcf", "in.sam"]))),
            1
        );
        assert_eq!(
            exit_to_u8(main(&argv(&["split", "--unknown", "in.sam"]))),
            1
        );
    }
}
