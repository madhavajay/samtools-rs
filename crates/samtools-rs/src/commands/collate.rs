//! `samtools collate` / `samtools bamshuf` — group reads with the same name
//! together, with a (pseudo-)randomised order otherwise.
//!
//! Mirrors `main_bamshuf` in `bamshuf.c`. The upstream implementation uses
//! name-hash bucketing with on-disk temp files for memory bounding. This
//! initial Rust port performs an in-memory name sort for BAM/SAM/reference-backed
//! CRAM inputs, which gives the same per-name grouping result but does not scale
//! to inputs larger than memory.

use std::collections::{HashMap, VecDeque};
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use bstr::BString;
use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf};

use crate::bam_flag::{BAM_FREAD1, BAM_FREAD2, BAM_FSECONDARY, BAM_FSUPPLEMENTARY};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools collate` (alias `bamshuf`).
pub fn main(args: &[OsString]) -> ExitCode {
    let mut output_prefix: Option<String> = None;
    let mut to_stdout = false;
    let mut output_fmt = OutFmt::Bam;
    let mut input: Option<PathBuf> = None;
    let mut no_pg = false;
    let mut fast = false;
    let mut reads_store = 10_000usize;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        if let Some(v) = s.strip_prefix("--output-fmt=") {
            output_fmt = match parse_output_format(v) {
                Ok(fmt) => fmt,
                Err(e) => {
                    print_error("collate", e);
                    return ExitCode::from(1);
                }
            };
            continue;
        }
        match s {
            "-o" => {
                output_prefix = iter.next().and_then(|a| a.to_str()).map(|s| s.to_string());
            }
            "-O" => {
                to_stdout = true;
            }
            "--output-fmt" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("collate", "missing value for --output-fmt");
                    return ExitCode::from(1);
                };
                output_fmt = match parse_output_format(v) {
                    Ok(fmt) => fmt,
                    Err(e) => {
                        print_error("collate", e);
                        return ExitCode::from(1);
                    }
                };
            }
            "--no-PG" => {
                no_pg = true;
            }
            "-f" => {
                fast = true;
            }
            "-r" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("collate", "missing value for -r");
                    return ExitCode::from(1);
                };
                reads_store = match v.parse::<usize>() {
                    Ok(n) => n.max(2),
                    Err(_) => {
                        print_error("collate", format!("invalid -r value \"{}\"", v));
                        return ExitCode::from(1);
                    }
                };
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error(
                    "collate",
                    format!("option `{}` is not yet supported in samtools-rs collate", s),
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

    let format = match sam_io::sam_open_format(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error("collate", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam | Exact::Cram) {
        print_error(
            "collate",
            "only SAM, BAM, and reference-backed CRAM input are currently supported",
        );
        return ExitCode::from(1);
    }

    let output = if to_stdout {
        OutputTarget::Stdout
    } else {
        let prefix = output_prefix.unwrap_or_else(|| "collated".to_string());
        let ext = match output_fmt {
            OutFmt::Sam => "sam",
            OutFmt::Bam => "bam",
        };
        let path = if prefix.ends_with(&format!(".{}", ext)) {
            PathBuf::from(prefix)
        } else {
            PathBuf::from(format!("{}.{}", prefix, ext))
        };
        OutputTarget::File(path)
    };

    match run_collate(
        &input,
        output,
        output_fmt,
        fast,
        reads_store,
        if no_pg { None } else { Some(args) },
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("collate", "collate failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Copy)]
enum OutFmt {
    Sam,
    Bam,
}

enum OutputTarget {
    Stdout,
    File(PathBuf),
}

fn parse_output_format(raw: &str) -> Result<OutFmt, String> {
    match raw.to_ascii_lowercase().as_str() {
        "sam" => Ok(OutFmt::Sam),
        "bam" => Ok(OutFmt::Bam),
        _ => Err(format!("unsupported output format \"{}\"", raw)),
    }
}

fn run_collate(
    input: &Path,
    output: OutputTarget,
    fmt: OutFmt,
    fast: bool,
    reads_store: usize,
    pg_argv: Option<&[OsString]>,
) -> io::Result<()> {
    let format = sam_io::sam_open_format(input)?;
    let (mut header, records) = match format.exact {
        Exact::Sam => read_sam_records(input)?,
        Exact::Bam => read_bam_records(input)?,
        Exact::Cram => read_cram_records(input)?,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "only SAM, BAM, and reference-backed CRAM input are currently supported",
            ));
        }
    };

    set_collate_header(&mut header);

    let records = if fast {
        fast_collate_records(records, reads_store)
    } else {
        sort_records_by_name(records)
    };

    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

    let mut writer = open_output(&output, fmt, &header)?;
    for rec in &records {
        writer.write_record(&header, rec)?;
    }
    Ok(())
}

fn sort_records_by_name(mut records: Vec<RecordBuf>) -> Vec<RecordBuf> {
    records.sort_by_key(name_key);
    records
}

fn fast_collate_records(records: Vec<RecordBuf>, reads_store: usize) -> Vec<RecordBuf> {
    let reads_store = reads_store.max(2);
    let mut paired = Vec::new();
    let mut deferred = Vec::new();
    let mut stored: HashMap<Vec<u8>, RecordBuf> = HashMap::new();
    let mut order = VecDeque::new();

    for record in records {
        if !is_fast_collate_candidate(&record) {
            continue;
        }

        let name = name_key(&record);
        if let Some(mate) = stored.remove(&name) {
            let (r1, r2) = order_pair(record, mate);
            paired.push(r1);
            paired.push(r2);
        } else {
            if stored.len() >= reads_store {
                flush_oldest_stored_record(&mut stored, &mut order, &mut deferred);
            }
            order.push_back(name.clone());
            stored.insert(name, record);
        }
    }

    while !stored.is_empty() {
        flush_oldest_stored_record(&mut stored, &mut order, &mut deferred);
    }

    paired.extend(sort_records_by_name(deferred));
    paired
}

fn flush_oldest_stored_record(
    stored: &mut HashMap<Vec<u8>, RecordBuf>,
    order: &mut VecDeque<Vec<u8>>,
    deferred: &mut Vec<RecordBuf>,
) {
    while let Some(name) = order.pop_front() {
        if let Some(record) = stored.remove(&name) {
            deferred.push(record);
            break;
        }
    }
}

fn is_fast_collate_candidate(record: &RecordBuf) -> bool {
    let flags = record.flags().bits() as u32;
    let read_flag = flags & (BAM_FREAD1 | BAM_FREAD2);
    flags & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY) == 0
        && matches!(read_flag, BAM_FREAD1 | BAM_FREAD2)
}

fn order_pair(a: RecordBuf, b: RecordBuf) -> (RecordBuf, RecordBuf) {
    if a.flags().is_first_segment() {
        (a, b)
    } else {
        (b, a)
    }
}

fn name_key(record: &RecordBuf) -> Vec<u8> {
    record.name().map(|s| s.to_vec()).unwrap_or_default()
}

fn set_collate_header(header: &mut sam::Header) {
    use sam::header::record::value::map::{self, Map};

    if let Some(hd) = header.header_mut() {
        hd.other_fields_mut()
            .insert(map::header::tag::SORT_ORDER, BString::from("unsorted"));
        hd.other_fields_mut()
            .insert(map::header::tag::GROUP_ORDER, BString::from("query"));
    } else {
        let mut hd: Map<map::Header> = Map::default();
        hd.other_fields_mut()
            .insert(map::header::tag::SORT_ORDER, BString::from("unsorted"));
        hd.other_fields_mut()
            .insert(map::header::tag::GROUP_ORDER, BString::from("query"));
        *header.header_mut() = Some(hd);
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

fn read_sam_records(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(input)?));
    read_sam_records_from_reader(&mut reader)
}

fn read_cram_records(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let reference = current_global_args().reference.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM input requires top-level --reference FILE",
        )
    })?;
    let text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            input, reference, None,
        )?;
    let mut reader = sam::io::Reader::new(BufReader::new(Cursor::new(text)));
    read_sam_records_from_reader(&mut reader)
}

fn read_sam_records_from_reader<R>(
    reader: &mut sam::io::Reader<R>,
) -> io::Result<(sam::Header, Vec<RecordBuf>)>
where
    R: BufRead,
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

trait CollateSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(sam::io::Writer<File>);
struct SamStdout(sam::io::Writer<io::Stdout>);

impl CollateSink for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl CollateSink for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl CollateSink for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl CollateSink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}

fn open_output(
    out: &OutputTarget,
    fmt: OutFmt,
    header: &sam::Header,
) -> io::Result<Box<dyn CollateSink>> {
    match (out, fmt) {
        (OutputTarget::Stdout, OutFmt::Sam) => {
            let mut writer = sam::io::Writer::new(io::stdout());
            writer.write_header(header)?;
            Ok(Box::new(SamStdout(writer)))
        }
        (OutputTarget::Stdout, OutFmt::Bam) => {
            let mut writer = bam::io::Writer::new(io::stdout());
            writer.write_header(header)?;
            Ok(Box::new(BamStdout(writer)))
        }
        (OutputTarget::File(p), OutFmt::Sam) => {
            let file = File::create(p)?;
            let mut writer = sam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(SamFile(writer)))
        }
        (OutputTarget::File(p), OutFmt::Bam) => {
            let file = File::create(p)?;
            let mut writer = bam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(BamFile(writer)))
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(
        w,
        "Usage: samtools collate [options] <in.bam|in.sam|in.cram> [<out.prefix>]"
    )?;
    writeln!(w, "  -o PREFIX   output prefix or path")?;
    writeln!(w, "  -O          write to stdout")?;
    writeln!(
        w,
        "  -f          fast mode: output primary read pairs early"
    )?;
    writeln!(w, "  -r INT      working reads stored with -f")?;
    writeln!(w, "  --output-fmt sam|bam")?;
    Ok(())
}
