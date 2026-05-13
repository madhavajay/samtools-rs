//! `samtools merge` — merge multiple sorted BAM files.
//!
//! Mirrors `bam_merge` in `bam_sort.c`. This initial Rust port loads all
//! records from the inputs into memory and sorts by coordinate (or name
//! with `-n`) before writing the merged output. K-way streaming merge
//! and CRAM are TODO.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::{Exact, detect_path};
use htslib_rs::sam;

use crate::diagnostics::{print_error, print_error_errno};

/// Entry point for `samtools merge`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut name_sort = false;
    let mut output: Option<PathBuf> = None;
    let mut output_fmt = OutFmt::Bam;
    let mut positional: Vec<PathBuf> = Vec::new();
    let mut force = false;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-n" => name_sort = true,
            "-f" => force = true,
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "--output-fmt" | "-O" => {
                let v = iter.next().and_then(|a| a.to_str()).unwrap_or("bam");
                output_fmt = match v.to_lowercase().as_str() {
                    "sam" => OutFmt::Sam,
                    "bam" => OutFmt::Bam,
                    _ => OutFmt::Bam,
                };
            }
            "-@" | "--threads" | "-l" | "--compression-level" | "-R" => {
                let _ = iter.next();
            }
            "--no-PG" | "--write-index" | "-c" | "-p" | "-u" => {}
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
    let (out_path, inputs): (Option<PathBuf>, Vec<PathBuf>) = if output.is_some() {
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

    if inputs.is_empty() {
        let _ = print_usage();
        return ExitCode::from(1);
    }

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
        let format = match detect_path(path) {
            Ok(f) => f,
            Err(e) => {
                print_error(
                    "merge",
                    format!("failed to detect format of \"{}\": {}", path.display(), e),
                );
                return ExitCode::from(1);
            }
        };
        if format.exact != Exact::Bam {
            print_error(
                "merge",
                format!(
                    "only BAM input is currently supported (got {:?} for \"{}\")",
                    format.exact,
                    path.display()
                ),
            );
            return ExitCode::from(1);
        }
    }

    match run_merge(&inputs, out_path.as_deref(), name_sort, output_fmt) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("merge", "merge failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum OutFmt {
    Sam,
    Bam,
}

pub(crate) fn run_merge(
    inputs: &[PathBuf],
    output: Option<&Path>,
    name_sort: bool,
    fmt: OutFmt,
) -> io::Result<()> {
    let mut first_reader = bam::io::Reader::new(File::open(&inputs[0])?);
    let mut header = first_reader.read_header()?;

    let mut records: Vec<bam::Record> = Vec::new();
    let mut record = bam::Record::default();
    loop {
        let n = first_reader.read_record(&mut record)?;
        if n == 0 {
            break;
        }
        records.push(record.clone());
    }
    drop(first_reader);

    for path in &inputs[1..] {
        let mut reader = bam::io::Reader::new(File::open(path)?);
        let _h = reader.read_header()?;
        loop {
            let n = reader.read_record(&mut record)?;
            if n == 0 {
                break;
            }
            records.push(record.clone());
        }
    }

    if name_sort {
        records.sort_by(|a, b| {
            let an = a.name().map(|s| s.to_vec()).unwrap_or_default();
            let bn = b.name().map(|s| s.to_vec()).unwrap_or_default();
            an.cmp(&bn)
        });
    } else {
        records.sort_by(|a, b| {
            let key = |r: &bam::Record| -> (i32, i64) {
                let tid = r
                    .reference_sequence_id()
                    .and_then(|res| res.ok())
                    .map(|t| t as i32)
                    .unwrap_or(i32::MAX);
                let pos = r
                    .alignment_start()
                    .and_then(|res| res.ok())
                    .map(|p| usize::from(p) as i64)
                    .unwrap_or(0);
                (tid, pos)
            };
            key(a).cmp(&key(b))
        });
    }

    set_sort_order(
        &mut header,
        if name_sort { "queryname" } else { "coordinate" },
    );

    let mut writer = open_output(output, fmt, &header)?;
    for rec in &records {
        writer.write_record(&header, rec)?;
    }
    Ok(())
}

fn set_sort_order(header: &mut sam::Header, so: &str) {
    use bstr::BString;
    use sam::header::record::value::map::{self, Map};
    if let Some(hd) = header.header_mut() {
        hd.other_fields_mut()
            .insert(map::header::tag::SORT_ORDER, BString::from(so));
    } else {
        let mut hd: Map<map::Header> = Map::default();
        hd.other_fields_mut()
            .insert(map::header::tag::SORT_ORDER, BString::from(so));
        *header.header_mut() = Some(hd);
    }
}

trait BamLike {
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(sam::io::Writer<File>);
struct SamStdout(sam::io::Writer<io::Stdout>);

impl BamLike for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        self.0.write_record(header, record)
    }
}
impl BamLike for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        self.0.write_record(header, record)
    }
}
impl BamLike for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl BamLike for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}

fn open_output(
    out: Option<&Path>,
    fmt: OutFmt,
    header: &sam::Header,
) -> io::Result<Box<dyn BamLike>> {
    match (out, fmt) {
        (Some(p), OutFmt::Sam) => {
            let file = File::create(p)?;
            let mut writer = sam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(SamFile(writer)))
        }
        (Some(p), OutFmt::Bam) => {
            let file = File::create(p)?;
            let mut writer = bam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(BamFile(writer)))
        }
        (None, OutFmt::Sam) => {
            let mut writer = sam::io::Writer::new(io::stdout());
            writer.write_header(header)?;
            Ok(Box::new(SamStdout(writer)))
        }
        (None, OutFmt::Bam) => {
            let mut writer = bam::io::Writer::new(io::stdout());
            writer.write_header(header)?;
            Ok(Box::new(BamStdout(writer)))
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(
        w,
        "Usage: samtools merge [options] <out.bam> <in1.bam> [<in2.bam> ...]"
    )?;
    writeln!(w, "Options:")?;
    writeln!(w, "  -n              name sort")?;
    writeln!(w, "  -f              force overwrite output")?;
    writeln!(w, "  -o FILE         output to FILE")?;
    writeln!(w, "  --output-fmt sam|bam")?;
    Ok(())
}
