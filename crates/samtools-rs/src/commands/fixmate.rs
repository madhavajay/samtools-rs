//! `samtools fixmate` — fix mate-related flags and positions on paired records.
//!
//! Mirrors `bam_mate.c` in upstream samtools. Initial Rust port handles
//! **name-sorted BAM/SAM input**: adjacent records with the same `qname` are
//! paired up and their `FMUNMAP`/`FMREVERSE` flags + `mate_reference_sequence_id`
//! + `mate_alignment_start` are made consistent.
//!
//! **Not yet supported:** CRAM input/output.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use bstr::BString;
use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{
    self,
    alignment::{
        RecordBuf,
        record::{Flags, cigar::op::Kind, data::field::Tag},
        record_buf::data::field::Value,
    },
};

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sanitize::{SanitizeFlags, parse_sanitize_options, sanitize_record};

/// Entry point for `samtools fixmate`.
pub fn main(args: &[OsString]) -> ExitCode {
    let opts = match parse_args(args) {
        Ok(opts) => opts,
        Err(ParseError::Usage) => {
            let _ = print_usage();
            return ExitCode::SUCCESS;
        }
        Err(ParseError::Err(e)) => {
            print_error("fixmate", e);
            return ExitCode::from(1);
        }
    };

    let Some(input) = opts.input else {
        let _ = print_usage();
        return ExitCode::from(1);
    };

    let format = match sam_io::sam_open_format(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error("fixmate", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam) {
        print_error(
            "fixmate",
            "only SAM and BAM input are currently supported (CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    let pg_argv = if opts.no_pg { None } else { Some(args) };
    let settings = FixmateSettings {
        remove_reads: opts.remove_reads,
        mate_score: opts.mate_score,
        add_template_cigar: opts.add_template_cigar,
        sanitize_flags: opts.sanitize_flags.unwrap_or(SanitizeFlags::ALL),
    };
    match run_fixmate(
        &input,
        opts.output.as_deref(),
        opts.output_fmt,
        pg_argv,
        settings,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("fixmate", "fixmate failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct Opts {
    output: Option<PathBuf>,
    input: Option<PathBuf>,
    output_fmt: OutFmt,
    sanitize_flags: Option<SanitizeFlags>,
    no_pg: bool,
    /// `-r`: remove secondary and unmapped reads from the output, and
    /// clear `PROPER_PAIR` / `MATE_REVERSE` on a pair where one mate is
    /// unmapped (mirrors upstream's `remove_reads`).
    remove_reads: bool,
    /// `-m`: add `ms:i` mate score tags for duplicate marking.
    mate_score: bool,
    /// `-c`: add lowercase `ct:Z` template CIGAR tag.
    add_template_cigar: bool,
}

impl Default for Opts {
    fn default() -> Self {
        Self {
            output: None,
            input: None,
            output_fmt: OutFmt::Bam,
            sanitize_flags: None,
            no_pg: false,
            remove_reads: false,
            mate_score: false,
            add_template_cigar: false,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum ParseError {
    Usage,
    Err(String),
}

fn parse_args(args: &[OsString]) -> Result<Opts, ParseError> {
    let mut opts = Opts::default();
    let mut iter = args.iter().skip(1);

    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-O" | "--output-fmt" => {
                let v = next_value(&mut iter, s)?;
                opts.output_fmt = parse_output_fmt(&v);
            }
            "-@" | "--threads" | "-l" => {
                let _ = next_value(&mut iter, s)?;
            }
            "-z" | "--sanitize" => {
                let v = next_value(&mut iter, s)?;
                opts.sanitize_flags = Some(parse_sanitize_options(&v).map_err(ParseError::Err)?);
            }
            "--no-PG" => {
                opts.no_pg = true;
            }
            "-r" => {
                opts.remove_reads = true;
            }
            "-m" => {
                opts.mate_score = true;
            }
            "-c" => {
                opts.add_template_cigar = true;
            }
            _ if s.starts_with("-cO") => {
                opts.add_template_cigar = true;
                let v = if s.len() > 3 {
                    s[3..].to_owned()
                } else {
                    next_value(&mut iter, "-O")?
                };
                opts.output_fmt = parse_output_fmt(&v);
            }
            "-p" => {
                // Accepted but not yet implemented.
            }
            "--help" => return Err(ParseError::Usage),
            _ if s.starts_with('-') && s != "-" => {
                return Err(ParseError::Err(format!("unknown option {}", s)));
            }
            _ => {
                if opts.input.is_none() {
                    opts.input = Some(PathBuf::from(arg));
                } else if opts.output.is_none() && s != "-" {
                    // A `-` output operand means stdout (output stays None).
                    opts.output = Some(PathBuf::from(arg));
                }
            }
        }
    }

    Ok(opts)
}

fn parse_output_fmt(raw: &str) -> OutFmt {
    match raw.to_lowercase().as_str() {
        "sam" => OutFmt::Sam,
        "bam" => OutFmt::Bam,
        _ => OutFmt::Bam,
    }
}

fn next_value<'a, I>(iter: &mut I, option: &str) -> Result<String, ParseError>
where
    I: Iterator<Item = &'a OsString>,
{
    iter.next()
        .and_then(|a| a.to_str().map(str::to_owned))
        .ok_or_else(|| ParseError::Err(format!("missing value for {option}")))
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum OutFmt {
    Sam,
    Bam,
}

#[derive(Clone, Copy)]
struct FixmateSettings {
    remove_reads: bool,
    mate_score: bool,
    add_template_cigar: bool,
    sanitize_flags: SanitizeFlags,
}

fn run_fixmate(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    settings: FixmateSettings,
) -> io::Result<()> {
    let format = sam_io::sam_open_format(input)?;
    match format.exact {
        Exact::Sam => run_fixmate_sam(input, output, fmt, pg_argv, settings),
        Exact::Bam => run_fixmate_bam(input, output, fmt, pg_argv, settings),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "only SAM and BAM input are currently supported (CRAM TODO)",
        )),
    }
}

fn run_fixmate_bam(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    settings: FixmateSettings,
) -> io::Result<()> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let mut header = reader.read_header()?;
    reject_coordinate_sorted(&header)?;
    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }
    let mut sink = open_output(output, fmt, &header)?;
    let mut pending: Option<RecordBuf> = None;
    let mut next = RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut next)?;
        if n == 0 {
            break;
        }
        sanitize_record(&header, &mut next, settings.sanitize_flags);
        write_fixed_record(
            &header,
            sink.as_mut(),
            &mut pending,
            next.clone(),
            settings.remove_reads,
            settings.mate_score,
            settings.add_template_cigar,
        )?;
    }
    if let Some(rec) = pending
        && !skip_for_remove_reads(&rec, settings.remove_reads)
    {
        sink.write_record(&header, &rec)?;
    }
    Ok(())
}

fn run_fixmate_sam(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    settings: FixmateSettings,
) -> io::Result<()> {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(input)?));
    let mut header = reader.read_header()?;
    reject_coordinate_sorted(&header)?;
    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }
    let mut sink = open_output(output, fmt, &header)?;
    let mut pending: Option<RecordBuf> = None;
    let mut next = RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut next)?;
        if n == 0 {
            break;
        }
        sanitize_record(&header, &mut next, settings.sanitize_flags);
        write_fixed_record(
            &header,
            sink.as_mut(),
            &mut pending,
            next.clone(),
            settings.remove_reads,
            settings.mate_score,
            settings.add_template_cigar,
        )?;
    }
    if let Some(rec) = pending
        && !skip_for_remove_reads(&rec, settings.remove_reads)
    {
        sink.write_record(&header, &rec)?;
    }
    Ok(())
}

fn reject_coordinate_sorted(header: &sam::Header) -> io::Result<()> {
    let is_coordinate_sorted = header
        .header()
        .as_ref()
        .and_then(|hd| {
            hd.other_fields()
                .get(&sam::header::record::value::map::header::tag::SORT_ORDER)
        })
        .is_some_and(|value| {
            let bytes: &[u8] = value.as_ref();
            bytes.eq_ignore_ascii_case(b"coordinate")
        });

    if is_coordinate_sorted {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Coordinate sorted, require grouped/sorted by queryname.",
        ))
    } else {
        Ok(())
    }
}

/// `-r` skips secondary and unmapped alignments from the output.
fn skip_for_remove_reads(rec: &RecordBuf, remove_reads: bool) -> bool {
    if !remove_reads {
        return false;
    }
    let flags = rec.flags();
    flags.is_secondary() || flags.is_unmapped()
}

fn write_fixed_record(
    header: &sam::Header,
    sink: &mut dyn Sink,
    pending: &mut Option<RecordBuf>,
    next: RecordBuf,
    remove_reads: bool,
    mate_score: bool,
    add_template_cigar: bool,
) -> io::Result<()> {
    match pending.take() {
        None => *pending = Some(next),
        Some(prev) => {
            let prev_name = prev.name().map(|n| n.to_vec());
            let next_name = next.name().map(|n| n.to_vec());
            if prev_name == next_name && next_name.is_some() {
                let (a, b) = pair_fixmate(prev, next, remove_reads, mate_score, add_template_cigar);
                if !skip_for_remove_reads(&a, remove_reads) {
                    sink.write_record(header, &a)?;
                }
                if !skip_for_remove_reads(&b, remove_reads) {
                    sink.write_record(header, &b)?;
                }
                *pending = None;
            } else {
                if !skip_for_remove_reads(&prev, remove_reads) {
                    sink.write_record(header, &prev)?;
                }
                *pending = Some(next);
            }
        }
    }
    Ok(())
}

fn pair_fixmate(
    mut a: RecordBuf,
    mut b: RecordBuf,
    remove_reads: bool,
    mate_score: bool,
    add_template_cigar: bool,
) -> (RecordBuf, RecordBuf) {
    let a_tid = a.reference_sequence_id();
    let b_tid = b.reference_sequence_id();
    let a_pos = a.alignment_start();
    let b_pos = b.alignment_start();
    apply_mate_flags(&mut a, &b);
    apply_mate_flags(&mut b, &a);
    update_mate_aux_tags(&mut a, &b);
    update_mate_aux_tags(&mut b, &a);
    update_template_lengths(&mut a, &mut b);
    if add_template_cigar {
        update_template_cigar_tag(&mut a, &mut b);
    }
    if mate_score {
        update_mate_score_tag(&mut a, &b);
        update_mate_score_tag(&mut b, &a);
    }
    *a.mate_reference_sequence_id_mut() = b_tid;
    *b.mate_reference_sequence_id_mut() = a_tid;
    *a.mate_alignment_start_mut() = b_pos;
    *b.mate_alignment_start_mut() = a_pos;

    if remove_reads {
        // Mirror upstream: when one mate is unmapped, clear MATE_REVERSE
        // and PROPER_PAIR on the surviving mate to avoid leaving stale
        // orphan-flag combinations after the unmapped mate is dropped.
        if a.flags().is_unmapped() {
            let mut f = b.flags();
            f.remove(Flags::MATE_REVERSE_COMPLEMENTED);
            f.remove(Flags::PROPERLY_SEGMENTED);
            *b.flags_mut() = f;
        }
        if b.flags().is_unmapped() {
            let mut f = a.flags();
            f.remove(Flags::MATE_REVERSE_COMPLEMENTED);
            f.remove(Flags::PROPERLY_SEGMENTED);
            *a.flags_mut() = f;
        }
    }
    (a, b)
}

fn update_template_lengths(a: &mut RecordBuf, b: &mut RecordBuf) {
    let tlen = if a.flags().is_unmapped()
        || b.flags().is_unmapped()
        || a.reference_sequence_id().is_none()
        || a.reference_sequence_id() != b.reference_sequence_id()
    {
        None
    } else {
        match (five_prime_position(a), five_prime_position(b)) {
            (Some(a5), Some(b5)) => Some((b5 - a5, a5 - b5)),
            _ => None,
        }
    };

    match tlen {
        Some((a_tlen, b_tlen)) => {
            *a.template_length_mut() = i32::try_from(a_tlen).unwrap_or(0);
            *b.template_length_mut() = i32::try_from(b_tlen).unwrap_or(0);
        }
        None => {
            *a.template_length_mut() = 0;
            *b.template_length_mut() = 0;
        }
    }
}

fn five_prime_position(record: &RecordBuf) -> Option<i64> {
    if record.flags().is_reverse_complemented() {
        record.alignment_end().map(|pos| pos.get() as i64)
    } else {
        record.alignment_start().map(|pos| pos.get() as i64 - 1)
    }
}

fn update_template_cigar_tag(a: &mut RecordBuf, b: &mut RecordBuf) {
    let ct_tag = Tag::from([b'c', b't']);
    a.data_mut().remove(&ct_tag);
    b.data_mut().remove(&ct_tag);

    if a.flags().is_unmapped()
        || b.flags().is_unmapped()
        || a.reference_sequence_id().is_none()
        || a.reference_sequence_id() != b.reference_sequence_id()
    {
        return;
    }

    let (Some(a_start), Some(b_start)) = (a.alignment_start(), b.alignment_start()) else {
        return;
    };
    let (Some(a_end), Some(b_end)) = (a.alignment_end(), b.alignment_end()) else {
        return;
    };

    let (left, left_start, left_end, right, right_start) = if a_start > b_start {
        (b, b_start, b_end, a, a_start)
    } else {
        (a, a_start, a_end, b, b_start)
    };

    let gap = right_start.get() as isize - left_end.get() as isize - 1;
    let ct = format!(
        "{}{}{}{}T{}{}{}",
        segment_index(left),
        strand(left),
        format_cigar(left.cigar()),
        gap,
        segment_index(right),
        strand(right),
        format_cigar(right.cigar())
    );
    debug_assert!(left_start <= right_start);
    left.data_mut()
        .insert(ct_tag, Value::String(BString::from(ct)));
}

fn segment_index(record: &RecordBuf) -> char {
    if record.flags().is_first_segment() {
        '1'
    } else {
        '2'
    }
}

fn strand(record: &RecordBuf) -> char {
    if record.flags().is_reverse_complemented() {
        'R'
    } else {
        'F'
    }
}

fn apply_mate_flags(target: &mut RecordBuf, mate: &RecordBuf) {
    let mut flags = target.flags();
    flags.insert(Flags::SEGMENTED);
    if mate.flags().is_unmapped() {
        flags.insert(Flags::MATE_UNMAPPED);
    } else {
        flags.remove(Flags::MATE_UNMAPPED);
    }
    if mate.flags().is_reverse_complemented() {
        flags.insert(Flags::MATE_REVERSE_COMPLEMENTED);
    } else {
        flags.remove(Flags::MATE_REVERSE_COMPLEMENTED);
    }
    *target.flags_mut() = flags;
}

fn update_mate_score_tag(target: &mut RecordBuf, mate: &RecordBuf) {
    target
        .data_mut()
        .insert(Tag::from([b'm', b's']), Value::from(mate_score(mate)));
}

fn mate_score(record: &RecordBuf) -> u32 {
    record
        .quality_scores()
        .iter()
        .filter(|&quality| quality >= 15)
        .map(u32::from)
        .sum()
}

fn update_mate_aux_tags(target: &mut RecordBuf, mate: &RecordBuf) {
    let mc_tag = Tag::from([b'M', b'C']);
    let mq_tag = Tag::from([b'M', b'Q']);

    if mate.flags().is_unmapped() {
        target.data_mut().remove(&mc_tag);
        target.data_mut().remove(&mq_tag);
        return;
    }

    // Upstream `bam_mate.c` adds MQ before MC.
    if let Some(mapping_quality) = mate.mapping_quality() {
        target
            .data_mut()
            .insert(mq_tag, Value::from(mapping_quality.get()));
    } else {
        target.data_mut().remove(&mq_tag);
    }

    target.data_mut().insert(
        mc_tag,
        Value::String(BString::from(format_cigar(mate.cigar()))),
    );
}

fn format_cigar(cigar: &sam::alignment::record_buf::Cigar) -> String {
    if cigar.as_ref().is_empty() {
        return String::from("*");
    }

    let mut s = String::new();
    for op in cigar.as_ref() {
        s.push_str(&op.len().to_string());
        s.push(match op.kind() {
            Kind::Match => 'M',
            Kind::Insertion => 'I',
            Kind::Deletion => 'D',
            Kind::Skip => 'N',
            Kind::SoftClip => 'S',
            Kind::HardClip => 'H',
            Kind::Pad => 'P',
            Kind::SequenceMatch => '=',
            Kind::SequenceMismatch => 'X',
        });
    }
    s
}

trait Sink {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(File);
struct SamStdout(io::Stdout);

impl Sink for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl Sink for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.0.write_alignment_record(header, record)
    }
}
impl Sink for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        // Shared renderer: htslib `%g` float aux spelling.
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}
impl Sink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}

fn open_output(out: Option<&Path>, fmt: OutFmt, header: &sam::Header) -> io::Result<Box<dyn Sink>> {
    match (out, fmt) {
        (Some(p), OutFmt::Sam) => {
            let mut w = File::create(p)?;
            crate::sam_render::write_header(&mut w, header)?;
            Ok(Box::new(SamFile(w)))
        }
        (Some(p), OutFmt::Bam) => {
            let mut w = bam::io::Writer::new(File::create(p)?);
            w.write_header(header)?;
            Ok(Box::new(BamFile(w)))
        }
        (None, OutFmt::Sam) => {
            let mut w = io::stdout();
            crate::sam_render::write_header(&mut w, header)?;
            Ok(Box::new(SamStdout(w)))
        }
        (None, OutFmt::Bam) => {
            let mut w = bam::io::Writer::new(io::stdout());
            w.write_header(header)?;
            Ok(Box::new(BamStdout(w)))
        }
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools fixmate [options] <in.bam> [<out.bam>]")?;
    writeln!(w, "  -O sam|bam   output format (default: bam)")?;
    writeln!(w, "  -m           add mate score tags")?;
    writeln!(w, "  -c           add template CIGAR ct tag")?;
    writeln!(w, "  -z, --sanitize FLAG[,FLAG]")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(rest: &[&str]) -> Vec<OsString> {
        std::iter::once(OsString::from("fixmate"))
            .chain(rest.iter().map(OsString::from))
            .collect()
    }

    #[test]
    fn parses_sanitize_option_without_treating_value_as_input() {
        let opts = parse_args(&argv(&["-z", "on", "in.bam", "out.bam"])).unwrap();

        assert_eq!(opts.input.as_deref(), Some(Path::new("in.bam")));
        assert_eq!(opts.output.as_deref(), Some(Path::new("out.bam")));
        assert!(opts.sanitize_flags.unwrap().contains(SanitizeFlags::CIGAR));
    }

    #[test]
    fn rejects_missing_sanitize_value() {
        assert_eq!(
            parse_args(&argv(&["--sanitize"])).unwrap_err(),
            ParseError::Err(String::from("missing value for --sanitize"))
        );
    }

    #[test]
    fn rejects_invalid_sanitize_value() {
        assert_eq!(
            parse_args(&argv(&["--sanitize", "nope"])).unwrap_err(),
            ParseError::Err(String::from("unrecognised sanitize keyword \"nope\""))
        );
    }
}
