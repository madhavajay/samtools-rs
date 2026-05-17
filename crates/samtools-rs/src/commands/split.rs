//! `samtools split` — split an alignment file by `@RG` ID or aux tag value.
//!
//! Mirrors `main_split` in `bam_split.c`. The initial Rust port supports:
//!  - Input BAM and SAM (CRAM TODO).
//!  - Output template via `-f` with `%*` (input basename), `%!` (RG ID or
//!    tag value), `%#` (output index), and `%.` (extension).
//!  - `-u <path>` — write records with unknown/missing RG/tag to this file.
//!  - `-h <path>` — use an alternate SAM header for the unaccounted output.
//!  - `-d TAG` — split by a string or integer aux tag instead of header `@RG`.
//!  - `-M N` / `--max-split N` — cap dynamically-created `-d` outputs.
//!  - `--output-fmt sam|bam` — only `bam` (default) and `sam` are honored.
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
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf};

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools split`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut template = String::from("%*_%!.%.");
    let mut unaccounted: Option<PathBuf> = None;
    let mut unaccounted_header: Option<PathBuf> = None;
    let mut output_fmt = OutFmt::Bam;
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
            "--output-fmt" => {
                let v = iter.next().and_then(|a| a.to_str()).unwrap_or("bam");
                output_fmt = match v.to_lowercase().as_str() {
                    "sam" => OutFmt::Sam,
                    "bam" => OutFmt::Bam,
                    _ => OutFmt::Bam,
                };
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
        let format = match sam_io::sam_open_format(&input) {
            Ok(f) => f,
            Err(e) => {
                print_error("split", e.to_string());
                return ExitCode::from(1);
            }
        };
        if !matches!(format.exact, Exact::Sam | Exact::Bam) {
            print_error(
                "split",
                "only SAM and BAM input are currently supported (CRAM TODO)",
            );
            return ExitCode::from(1);
        }
    }

    let options = SplitOptions {
        template: &template,
        unaccounted: unaccounted.as_deref(),
        unaccounted_header: unaccounted_header.as_deref(),
        fmt: output_fmt,
        pad_width,
        split_tag,
        max_split,
        add_pg,
        write_index: local_write_index || current_global_args().write_index,
        argv: args,
    };

    match run_split(&input, options) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("split", "split failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Copy)]
enum OutFmt {
    Sam,
    Bam,
}

struct SplitOptions<'a> {
    template: &'a str,
    unaccounted: Option<&'a Path>,
    unaccounted_header: Option<&'a Path>,
    fmt: OutFmt,
    pad_width: usize,
    split_tag: Option<[u8; 2]>,
    max_split: usize,
    add_pg: bool,
    write_index: bool,
    argv: &'a [OsString],
}

fn run_split(input: &Path, options: SplitOptions<'_>) -> io::Result<()> {
    let (header, records) = if input.as_os_str() == "-" {
        read_stdin_records()?
    } else {
        let format = sam_io::sam_open_format(input)?;
        match format.exact {
            Exact::Sam => read_sam_records(input)?,
            Exact::Bam => read_bam_records(input)?,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "only SAM and BAM input are currently supported (CRAM TODO)",
                ));
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
            options.fmt,
            options.pad_width,
        );
        let path = PathBuf::from(path);
        track_index_path(&mut index_paths, options.write_index, options.fmt, &path);
        outputs.push(SplitOutput {
            header: output_header.clone(),
            sink: open_output(&path, options.fmt, &output_header)?,
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
            track_index_path(&mut index_paths, options.write_index, options.fmt, p);
            Some(open_output(p, options.fmt, &unaccounted_output_header)?)
        }
        None => None,
    };

    for record in &records {
        let tag_value = record_split_tag_value(record, options.split_tag, options.pad_width);
        let target = match tag_value {
            Some(value) => match id_to_index.get(&value).copied() {
                Some(i) => Some(SplitTarget::Output(i)),
                None if options.split_tag.is_some() && outputs.len() < options.max_split => {
                    let idx = outputs.len();
                    let output_header = prepare_output_header(
                        &header,
                        options.split_tag,
                        &value,
                        options.add_pg,
                        options.argv,
                    )?;
                    let path = render_template(
                        options.template,
                        input_base,
                        &value,
                        idx,
                        options.fmt,
                        options.pad_width,
                    );
                    let path = PathBuf::from(path);
                    track_index_path(&mut index_paths, options.write_index, options.fmt, &path);
                    outputs.push(SplitOutput {
                        header: output_header.clone(),
                        sink: open_output(&path, options.fmt, &output_header)?,
                    });
                    id_to_index.insert(value, idx);
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
            None => {}
        }
    }

    drop(outputs);
    drop(unaccounted_out);
    for path in index_paths {
        write_bam_index(&path)?;
    }
    Ok(())
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
) -> String {
    let ext = match fmt {
        OutFmt::Sam => "sam",
        OutFmt::Bam => "bam",
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
                Some(other) => {
                    out.push('%');
                    out.push(other);
                }
                None => out.push('%'),
            }
        } else {
            out.push(c);
        }
    }
    out
}

trait SplitSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
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

fn open_output(path: &Path, fmt: OutFmt, header: &sam::Header) -> io::Result<Box<dyn SplitSink>> {
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
    writeln!(w, "  --output-fmt sam|bam")?;
    writeln!(w, "  --no-PG    do not add a @PG line")?;
    writeln!(w, "  --write-index")?;
    Ok(())
}
