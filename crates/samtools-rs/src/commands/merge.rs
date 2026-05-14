//! `samtools merge` — merge multiple sorted BAM files.
//!
//! Mirrors `bam_merge` in `bam_sort.c`. This initial Rust port loads all
//! records from BAM/SAM inputs into memory and sorts by coordinate (or name
//! with `-n`) before writing the merged output. K-way streaming merge
//! and CRAM are TODO. `-R` and `-L` restrict indexed BAM inputs by region/BED.
//! Coordinate-sorted BAM outputs can also write a BAI via `--write-index`.

use std::cmp::Ordering;
use std::collections::HashSet;
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
            "-@" | "--threads" | "-l" | "--compression-level" => {
                let _ = iter.next();
            }
            "--write-index" => local_write_index = true,
            "--no-PG" => {
                no_pg = true;
            }
            "-s" => {
                let _ = iter.next();
            }
            "-c" | "-p" | "-u" => {}
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
    let order = match tag_sort {
        Some(tag) => MergeOrder::Tag {
            tag,
            name_secondary: name_sort,
        },
        None if name_sort => MergeOrder::Name,
        None => MergeOrder::Coordinate,
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
        write_index,
        if no_pg { None } else { Some(args) },
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

#[derive(Clone, Copy)]
pub(crate) enum OutFmt {
    Sam,
    Bam,
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
}

pub(crate) fn run_merge(
    inputs: &[PathBuf],
    output: Option<&Path>,
    order: MergeOrder,
    fmt: OutFmt,
    write_index: bool,
    pg_argv: Option<&[OsString]>,
    restriction: MergeRestriction<'_>,
) -> io::Result<()> {
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

    let (mut header, mut records) = read_records(&inputs[0], filter.as_ref())?;
    let first_reference_id_map: Vec<_> = (0..header.reference_sequences().len()).collect();
    remap_records(&mut records, &first_reference_id_map)?;

    for path in &inputs[1..] {
        let (input_header, mut input_records) = read_records(path, filter.as_ref())?;
        let reference_id_map = merge_reference_sequences(&mut header, &input_header)?;
        merge_read_groups(&mut header, &input_header)?;
        merge_programs(&mut header, &input_header)?;
        merge_comments(&mut header, &input_header);
        remap_records(&mut input_records, &reference_id_map)?;
        records.append(&mut input_records);
    }

    match order {
        MergeOrder::Tag {
            tag,
            name_secondary,
        } => records.sort_by(|a, b| compare_by_tag(a, b, tag, name_secondary)),
        MergeOrder::Name => records.sort_by_key(name_key),
        MergeOrder::Coordinate => records.sort_by_key(coordinate_key),
    }

    if let MergeOrder::Tag {
        tag,
        name_secondary,
    } = order
    {
        set_sort_order(
            &mut header,
            "unsorted",
            Some(&format!(
                "unsorted:{}{}:{}",
                tag[0] as char,
                tag[1] as char,
                if name_secondary {
                    "queryname:lexicographical"
                } else {
                    "coordinate"
                }
            )),
        );
    } else {
        set_sort_order(
            &mut header,
            if matches!(order, MergeOrder::Name) {
                "queryname"
            } else {
                "coordinate"
            },
            None,
        );
    }

    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

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
            read_bam_records_in_regions(input, regions)
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
    regions: &[htslib_rs::core::Region],
) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let header = htslib_rs::alignment_compat::read_bam_header_from_path(input)?;
    let mut records = Vec::new();
    let mut seen = HashSet::new();
    for region in regions {
        let bam_records = htslib_rs::alignment_compat::query_bam_records_from_path(input, region)?;
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
                .reference_sequences()
                .get_index(index)
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
    for (id, input_read_group) in input_header.read_groups() {
        if let Some(output_read_group) = output_header.read_groups().get(id) {
            if output_read_group != input_read_group {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "conflicting @RG definition for {}",
                        String::from_utf8_lossy(id)
                    ),
                ));
            }
        } else {
            output_header
                .read_groups_mut()
                .insert(id.clone(), input_read_group.clone());
        }
    }

    Ok(())
}

fn merge_programs(output_header: &mut sam::Header, input_header: &sam::Header) -> io::Result<()> {
    for (id, input_program) in input_header.programs().as_ref() {
        if let Some(output_program) = output_header.programs().as_ref().get(id) {
            if output_program != input_program {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "conflicting @PG definition for {}",
                        String::from_utf8_lossy(id)
                    ),
                ));
            }
        } else {
            output_header
                .programs_mut()
                .as_mut()
                .insert(id.clone(), input_program.clone());
        }
    }

    Ok(())
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
                name_key(a).cmp(&name_key(b))
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

trait MergeSink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(sam::io::Writer<File>);
struct SamStdout(sam::io::Writer<io::Stdout>);

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
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl MergeSink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}

fn open_output(
    out: Option<&Path>,
    fmt: OutFmt,
    header: &sam::Header,
) -> io::Result<Box<dyn MergeSink>> {
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
    writeln!(w, "  --output-fmt sam|bam")?;
    writeln!(w, "  --write-index   write BAI index for BAM file output")?;
    Ok(())
}
