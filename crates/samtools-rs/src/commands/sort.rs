//! `samtools sort` — sort alignment records.
//!
//! Mirrors `main_sort` in `bam_sort.c`. The upstream implementation is the
//! largest single file in samtools (138k LOC) and supports external k-way
//! merge with temp files, name/coordinate/tag/template-coordinate sort,
//! and many auxiliary flags.
//!
//! This initial Rust port supports **in-memory coordinate sort or name sort
//! for BAM/SAM**, which is sufficient for small/medium inputs. Records are
//! sorted by `(reference_sequence_id, alignment_start)` for coordinate mode
//! or by `qname` for name mode, then written to the output.
//!
//! Supported flags:
//!  - `-n` — name sort (default is coordinate sort).
//!  - `-o FILE` — output file (default stdout).
//!  - `-O sam|bam`, `--output-fmt sam|bam` — output format (default: bam).
//!  - `-@`/`--threads`, `-m`/`--max-mem`, `-T`/`--temp` — accepted but ignored.
//!  - `--no-PG` — accepted, silently ignored.
//!  - `--write-index` — write a BAI next to coordinate-sorted BAM output.
//!
//! Not yet supported: external merge (large inputs spill to disk), tag
//! sort (`-t TAG`), template-coordinate sort (`-M`), minimiser sort (`-N`),
//! CRAM I/O.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf};

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools sort`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut name_sort = false;
    let mut output: Option<PathBuf> = None;
    let mut output_fmt = OutFmt::Bam;
    let mut input: Option<PathBuf> = None;
    let mut local_write_index = false;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-n" | "--name" => {
                name_sort = true;
            }
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "-O" | "--output-fmt" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("sort", format!("missing value for {}", s));
                    return ExitCode::from(1);
                };
                output_fmt = match parse_output_format(v) {
                    Ok(fmt) => fmt,
                    Err(e) => {
                        print_error("sort", e);
                        return ExitCode::from(1);
                    }
                };
            }
            "-@"
            | "--threads"
            | "-m"
            | "--max-mem"
            | "-T"
            | "--temp"
            | "-l"
            | "--compression-level"
            | "-K" => {
                let _ = iter.next();
            }
            "--write-index" => {
                local_write_index = true;
            }
            "--no-PG" | "-u" => {
                // Accepted but ignored for the in-memory port.
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error(
                    "sort",
                    format!("option `{}` is not yet supported in samtools-rs sort", s),
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
            print_error("sort", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam) {
        print_error(
            "sort",
            "only SAM and BAM input are currently supported (CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    let write_index = local_write_index || current_global_args().write_index;
    if write_index {
        if output.is_none() {
            print_error("sort", "--write-index requires -o FILE");
            return ExitCode::from(1);
        }
        if name_sort {
            print_error("sort", "--write-index requires coordinate sort output");
            return ExitCode::from(1);
        }
        if !matches!(output_fmt, OutFmt::Bam) {
            print_error("sort", "--write-index is only supported for BAM output");
            return ExitCode::from(1);
        }
    }

    match run_sort(
        &input,
        output.as_deref(),
        name_sort,
        output_fmt,
        write_index,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("sort", "sort failed", &e);
            ExitCode::from(1)
        }
    }
}

fn parse_output_format(raw: &str) -> Result<OutFmt, String> {
    match raw.to_ascii_lowercase().as_str() {
        "sam" => Ok(OutFmt::Sam),
        "bam" => Ok(OutFmt::Bam),
        _ => Err(format!("unsupported output format \"{}\"", raw)),
    }
}

#[derive(Clone, Copy)]
pub(crate) enum OutFmt {
    Sam,
    Bam,
}

pub(crate) fn run_sort(
    input: &Path,
    output: Option<&Path>,
    name_sort: bool,
    fmt: OutFmt,
    write_index: bool,
) -> io::Result<()> {
    let format = sam_io::sam_open_format(input)?;
    let (mut header, mut records) = match format.exact {
        Exact::Sam => read_sam_records(input)?,
        Exact::Bam => read_bam_records(input)?,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "only SAM and BAM input are currently supported (CRAM TODO)",
            ));
        }
    };

    if name_sort {
        records.sort_by(|a, b| {
            let an = a.name().map(|s| s.to_vec()).unwrap_or_default();
            let bn = b.name().map(|s| s.to_vec()).unwrap_or_default();
            an.cmp(&bn)
        });
    } else {
        // Coordinate sort: by (reference_sequence_id, alignment_start).
        records.sort_by(|a, b| {
            // Records with no reference (unmapped) sort to the end.
            let key = |r: &RecordBuf| -> (i32, i64) {
                let tid = r
                    .reference_sequence_id()
                    .map(|t| t as i32)
                    .unwrap_or(i32::MAX);
                let pos = r.alignment_start().map(usize::from).unwrap_or(0) as i64;
                (tid, pos)
            };
            key(a).cmp(&key(b))
        });
    }

    // Update @HD SO to reflect new sort order so downstream consumers can
    // tell. (Header is otherwise preserved verbatim.)
    set_sort_order(
        &mut header,
        if name_sort { "queryname" } else { "coordinate" },
    );

    {
        let mut writer = open_output(output, fmt, &header)?;
        for rec in &records {
            writer.write_record(&header, rec)?;
        }
    }

    if write_index {
        let Some(path) = output else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--write-index requires -o FILE",
            ));
        };
        write_bam_index(path)?;
    }
    Ok(())
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

trait SortSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(sam::io::Writer<File>);
struct SamStdout(sam::io::Writer<io::Stdout>);

impl SortSink for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl SortSink for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl SortSink for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl SortSink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}

fn open_output(
    out: Option<&Path>,
    fmt: OutFmt,
    header: &sam::Header,
) -> io::Result<Box<dyn SortSink>> {
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
    writeln!(w, "Usage: samtools sort [options] <in.bam|in.sam>")?;
    writeln!(
        w,
        "  -n              sort by read name (default: coordinate)"
    )?;
    writeln!(w, "  -o FILE         write output to FILE (default stdout)")?;
    writeln!(w, "  --output-fmt sam|bam")?;
    writeln!(
        w,
        "  -@/-m/-T/-K     accepted but currently ignored (in-memory sort only)"
    )?;
    writeln!(w, "  --write-index   write BAI index for BAM file output")?;
    Ok(())
}
