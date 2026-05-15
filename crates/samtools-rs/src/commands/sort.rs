//! `samtools sort` — sort alignment records.
//!
//! Mirrors `main_sort` in `bam_sort.c`. The upstream implementation is the
//! largest single file in samtools (138k LOC) and supports external k-way
//! merge with temp files, name/coordinate/tag/template-coordinate sort,
//! and many auxiliary flags.
//!
//! This initial Rust port supports **in-memory coordinate, name, or tag sort
//! for BAM/SAM/reference-backed CRAM**, which is sufficient for small/medium inputs. Records are
//! sorted by `(reference_sequence_id, alignment_start)` for coordinate mode,
//! by `qname` for name mode, or by `TAG` with coordinate/name secondary keys
//! for tag mode, then written to the output.
//!
//! Supported flags:
//!  - `-n` — name sort (default is coordinate sort).
//!  - `-t TAG` — sort by auxiliary tag, using coordinate/name as secondary key.
//!  - `-o FILE` — output file (default stdout).
//!  - `-O sam|bam`, `--output-fmt sam|bam` — output format (default: bam).
//!  - `-@`/`--threads`, `-m`/`--max-mem`, `-T`/`--temp` — accepted but ignored.
//!  - `--no-PG` — accepted, silently ignored.
//!  - `--write-index` — write a BAI next to coordinate-sorted BAM output.
//!
//! Not yet supported: external merge (large inputs spill to disk),
//! template-coordinate sort (`-M`), minimiser sort (`-N`), CRAM output.

use std::cmp::Ordering;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Write};
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
    let mut no_pg = false;
    let mut tag_sort: Option<[u8; 2]> = None;

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-n" | "--name" => {
                name_sort = true;
            }
            "-o" | "--output" => {
                output = match iter.next().and_then(|a| a.to_str()) {
                    // `-o -` means stdout (output stays None).
                    Some("-") | None => None,
                    Some(p) => Some(PathBuf::from(p)),
                };
            }
            "-t" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("sort", "missing value for -t");
                    return ExitCode::from(1);
                };
                tag_sort = match parse_tag(v) {
                    Ok(tag) => Some(tag),
                    Err(e) => {
                        print_error("sort", e);
                        return ExitCode::from(1);
                    }
                };
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
            "--no-PG" => {
                no_pg = true;
            }
            "-u" => {
                // Accepted but currently ignored (controls uncompressed output).
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            // Attached-value forms of the accepted-but-ignored options
            // (`-@4`, `-m768M`, `-l6`, `-Kprefix`, `-Tprefix`, `--threads=4`).
            _ if (s.starts_with("-@")
                || s.starts_with("-m")
                || s.starts_with("-l")
                || s.starts_with("-K")
                || s.starts_with("-T"))
                && s.len() > 2
                && !s.starts_with("--") =>
            {
                // value is in the same token; nothing to consume.
            }
            _ if s.starts_with("--threads=")
                || s.starts_with("--max-mem=")
                || s.starts_with("--compression-level=")
                || s.starts_with("--temp=") =>
            {
                // value embedded; ignored.
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
    if !matches!(format.exact, Exact::Sam | Exact::Bam | Exact::Cram) {
        print_error(
            "sort",
            "only SAM, BAM, and reference-backed CRAM input are currently supported",
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
        tag_sort,
        output_fmt,
        write_index,
        if no_pg { None } else { Some(args) },
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

fn parse_tag(raw: &str) -> Result<[u8; 2], String> {
    let bytes = raw.as_bytes();
    if bytes.len() == 2 {
        Ok([bytes[0], bytes[1]])
    } else {
        Err(format!("sort tag must be exactly two bytes, got {:?}", raw))
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
    tag_sort: Option<[u8; 2]>,
    fmt: OutFmt,
    write_index: bool,
    pg_argv: Option<&[OsString]>,
) -> io::Result<()> {
    let format = sam_io::sam_open_format(input)?;
    let (mut header, mut records) = match format.exact {
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

    if let Some(tag) = tag_sort {
        records.sort_by(|a, b| compare_by_tag(a, b, tag, name_sort));
    } else if name_sort {
        records.sort_by(name_cmp);
    } else {
        // Coordinate sort: by (reference_sequence_id, alignment_start).
        records.sort_by(|a, b| {
            // Records with no reference (unmapped) sort to the end.
            coordinate_key(a).cmp(&coordinate_key(b))
        });
    }

    // Sort-order tags for @HD.
    let (so, ss): (String, Option<String>) = if let Some(tag) = tag_sort {
        (
            "unsorted".to_string(),
            Some(format!(
                "unsorted:{}{}:{}",
                tag[0] as char,
                tag[1] as char,
                if name_sort {
                    "queryname:natural"
                } else {
                    "coordinate"
                }
            )),
        )
    } else if name_sort {
        (
            "queryname".to_string(),
            Some("queryname:natural".to_string()),
        )
    } else {
        ("coordinate".to_string(), None)
    };
    set_sort_order(&mut header, &so, ss.as_deref());

    // Emit the *raw* input header (preserving @SQ/@RG field order, @CO,
    // etc. — noodles' canonical writer reorders @RG fields) with the @HD
    // SO/SS applied and the samtools @PG appended.
    let mut header_text = apply_hd_sort_order(
        &crate::header_text::read_raw_header_text_with_format(input, format.exact)?,
        &so,
        ss.as_deref(),
    );
    if let Some(argv) = pg_argv {
        header_text = crate::pg::add_samtools_pg(&header_text, argv).map_err(io::Error::other)?;
    }

    {
        let mut writer = open_output(output, fmt, &header, &header_text)?;
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

/// Port of `bam_sort.c`'s `strnum_cmp` natural-order comparison: runs of
/// digits compare numerically (leading zeros skipped, then by length,
/// then by first differing digit); everything else byte-wise.
fn strnum_cmp(a: &[u8], b: &[u8]) -> Ordering {
    let (mut ia, mut ib) = (0usize, 0usize);
    let is_digit = |c: u8| c.is_ascii_digit();
    while ia < a.len() && ib < b.len() {
        let (ca, cb) = (a[ia], b[ib]);
        if !is_digit(ca) || !is_digit(cb) {
            if ca != cb {
                return ca.cmp(&cb);
            }
            ia += 1;
            ib += 1;
        } else {
            while ia < a.len() && a[ia] == b'0' {
                ia += 1;
            }
            while ib < b.len() && b[ib] == b'0' {
                ib += 1;
            }
            while ia < a.len() && ib < b.len() && is_digit(a[ia]) && a.get(ia) == b.get(ib) {
                ia += 1;
                ib += 1;
            }
            let diff =
                a.get(ia).copied().unwrap_or(0) as i32 - b.get(ib).copied().unwrap_or(0) as i32;
            while ia < a.len() && ib < b.len() && is_digit(a[ia]) && is_digit(b[ib]) {
                ia += 1;
                ib += 1;
            }
            if ia < a.len() && is_digit(a[ia]) {
                return Ordering::Greater;
            } else if ib < b.len() && is_digit(b[ib]) {
                return Ordering::Less;
            } else if diff != 0 {
                return diff.cmp(&0);
            }
        }
    }
    let ra = ia < a.len();
    let rb = ib < b.len();
    if ra {
        Ordering::Greater
    } else if rb {
        Ordering::Less
    } else {
        Ordering::Equal
    }
}

/// `bam_sort.c` QueryName secondary key:
/// `((f&0xc0)<<8)|((f&0x100)<<3)|((f&0x800)>>3)` — READ1, READ2,
/// (primary), SUPPLEMENTARY, SECONDARY.
fn qname_flag_key(r: &RecordBuf) -> u32 {
    let f = u32::from(u16::from(r.flags()));
    ((f & 0xc0) << 8) | ((f & 0x100) << 3) | ((f & 0x800) >> 3)
}

/// Full `bam_sort.c` QueryName comparator.
fn name_cmp(a: &RecordBuf, b: &RecordBuf) -> Ordering {
    strnum_cmp(&name_key(a), &name_key(b)).then_with(|| qname_flag_key(a).cmp(&qname_flag_key(b)))
}

fn compare_by_tag(a: &RecordBuf, b: &RecordBuf, tag: [u8; 2], name_sort: bool) -> Ordering {
    tag_sort_value(a, tag)
        .cmp(&tag_sort_value(b, tag))
        .then_with(|| {
            if name_sort {
                name_cmp(a, b)
            } else {
                coordinate_key(a).cmp(&coordinate_key(b))
            }
        })
        .then_with(|| name_key(a).cmp(&name_key(b)))
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

fn set_sort_order(header: &mut sam::Header, so: &str, ss: Option<&str>) {
    use bstr::BString;
    use sam::header::record::value::map::{self, Map};
    if let Some(hd) = header.header_mut() {
        hd.other_fields_mut()
            .insert(map::header::tag::SORT_ORDER, BString::from(so));
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
        if let Some(ss) = ss {
            hd.other_fields_mut()
                .insert(map::header::tag::SUBSORT_ORDER, BString::from(ss));
        }
        *header.header_mut() = Some(hd);
    }
}

trait SortSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(File);
struct SamStdout(io::Stdout);

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
        // Shared renderer: htslib `%g` float aux spelling.
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}
impl SortSink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}

/// Replaces/sets the `@HD` line's `SO:`/`SS:` fields in raw header text,
/// preserving every other line and field verbatim (so `@RG`/`@SQ`/`@CO`
/// keep their original byte form). Inserts an `@HD` if absent.
fn apply_hd_sort_order(raw: &str, so: &str, ss: Option<&str>) -> String {
    let mut lines: Vec<String> = Vec::new();
    let mut had_hd = false;
    for line in raw.lines() {
        if line.starts_with("@HD") {
            had_hd = true;
            let mut fields: Vec<&str> = line
                .split('\t')
                .filter(|f| !f.starts_with("SO:") && !f.starts_with("SS:"))
                .collect();
            let mut nl = fields.join("\t");
            if fields.is_empty() {
                nl.push_str("@HD");
            }
            nl.push_str(&format!("\tSO:{so}"));
            if let Some(ss) = ss {
                nl.push_str(&format!("\tSS:{ss}"));
            }
            lines.push(nl);
            let _ = &mut fields;
        } else {
            lines.push(line.to_string());
        }
    }
    if !had_hd {
        let hd = match ss {
            Some(ss) => format!("@HD\tVN:1.6\tSO:{so}\tSS:{ss}"),
            None => format!("@HD\tVN:1.6\tSO:{so}"),
        };
        lines.insert(0, hd);
    }
    let mut s = lines.join("\n");
    s.push('\n');
    s
}

fn open_output(
    out: Option<&Path>,
    fmt: OutFmt,
    header: &sam::Header,
    header_text: &str,
) -> io::Result<Box<dyn SortSink>> {
    match (out, fmt) {
        (Some(p), OutFmt::Sam) => {
            let mut file = File::create(p)?;
            file.write_all(header_text.as_bytes())?;
            Ok(Box::new(SamFile(file)))
        }
        (Some(p), OutFmt::Bam) => {
            let file = File::create(p)?;
            let mut writer = bam::io::Writer::new(file);
            writer.write_header(header)?;
            Ok(Box::new(BamFile(writer)))
        }
        (None, OutFmt::Sam) => {
            let mut stdout = io::stdout();
            stdout.write_all(header_text.as_bytes())?;
            Ok(Box::new(SamStdout(stdout)))
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
    writeln!(w, "Usage: samtools sort [options] <in.bam|in.sam|in.cram>")?;
    writeln!(
        w,
        "  -n              sort by read name (default: coordinate)"
    )?;
    writeln!(
        w,
        "  -t TAG          sort by auxiliary tag, then coordinate/name"
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
