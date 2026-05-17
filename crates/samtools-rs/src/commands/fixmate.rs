//! `samtools fixmate` — fix mate-related flags and positions on paired records.
//!
//! Mirrors `bam_mate.c` in upstream samtools. Initial Rust port handles
//! **name-sorted BAM/SAM input**: adjacent records with the same `qname` are
//! paired up and their `FMUNMAP`/`FMREVERSE` flags + `mate_reference_sequence_id`
//! + `mate_alignment_start` are made consistent.
//!
//! CRAM output is supported when a top-level `--reference` is provided.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read, Write};
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
        record_buf::data::field::{Value, value::Array},
    },
};

use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
use crate::io as sam_io;
use crate::sam_global::current_global_args;
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

    if !is_stdin_path(&input) {
        if !input.exists() {
            print_hts_open_missing(&input);
            print_error(
                "fixmate",
                "cannot open input file: No such file or directory",
            );
            return ExitCode::from(1);
        }
        let format = match sam_io::sam_open_format(&input) {
            Ok(f) => f,
            Err(e) => {
                print_error("fixmate", e.to_string());
                return ExitCode::from(1);
            }
        };
        let exact = fixmate_input_exact(&input, format.exact);
        if !matches!(exact, Exact::Sam | Exact::Bam | Exact::Cram) {
            print_error(
                "fixmate",
                "only SAM, BAM, and reference-backed CRAM input are currently supported",
            );
            return ExitCode::from(1);
        }
    }

    let pg_argv = if opts.no_pg { None } else { Some(args) };
    let settings = FixmateSettings {
        remove_reads: opts.remove_reads,
        mate_score: opts.mate_score,
        add_template_cigar: opts.add_template_cigar,
        base_mods: opts.base_mods,
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
    /// `-M`: fix MM/ML/MN base-modification tags on secondary/supplementary
    /// alignments.
    base_mods: bool,
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
            base_mods: false,
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
            _ if s.starts_with("--output-fmt=") => {
                opts.output_fmt = parse_output_fmt(s.trim_start_matches("--output-fmt="));
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
            "-M" => {
                opts.base_mods = true;
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
        "cram" => OutFmt::Cram,
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
    Cram,
}

#[derive(Clone, Copy)]
struct FixmateSettings {
    remove_reads: bool,
    mate_score: bool,
    add_template_cigar: bool,
    base_mods: bool,
    sanitize_flags: SanitizeFlags,
}

fn run_fixmate(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    settings: FixmateSettings,
) -> io::Result<()> {
    if is_stdin_path(input) {
        return run_fixmate_stdin(output, fmt, pg_argv, settings);
    }

    let format = sam_io::sam_open_format(input)?;
    match fixmate_input_exact(input, format.exact) {
        Exact::Sam => run_fixmate_sam(input, output, fmt, pg_argv, settings),
        Exact::Bam => run_fixmate_bam(input, output, fmt, pg_argv, settings),
        Exact::Cram => run_fixmate_cram(input, output, fmt, pg_argv, settings),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "only SAM, BAM, and reference-backed CRAM input are currently supported",
        )),
    }
}

fn is_stdin_path(path: &Path) -> bool {
    path.as_os_str() == "-"
}

fn fixmate_input_exact(path: &Path, detected: Exact) -> Exact {
    if detected != Exact::Unknown {
        return detected;
    }

    match path
        .extension()
        .and_then(|ext| ext.to_str())
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        Some("sam") => Exact::Sam,
        Some("bam") => Exact::Bam,
        Some("cram") => Exact::Cram,
        _ => detected,
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
    let mut sink = open_output(output, fmt, &header, None)?;
    let mut template = Vec::new();
    let mut next = RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut next)?;
        if n == 0 {
            break;
        }
        sanitize_record(&header, &mut next, settings.sanitize_flags);
        write_or_buffer_template(
            &header,
            sink.as_mut(),
            &mut template,
            next.clone(),
            settings,
        )?;
    }
    flush_template(&header, sink.as_mut(), &mut template, settings)?;
    sink.finish()
}

fn run_fixmate_cram(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    settings: FixmateSettings,
) -> io::Result<()> {
    let reference = current_global_args().reference.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM input requires top-level --reference FILE",
        )
    })?;
    crate::reference::ensure_fai_index(&reference, None)?;
    let sam =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            input, &reference, None,
        )?;
    run_fixmate_sam_bytes(sam.into_bytes(), output, fmt, pg_argv, settings)
}

fn run_fixmate_sam(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    settings: FixmateSettings,
) -> io::Result<()> {
    let mut reader = crate::sam_compat::open_sam_reader_tolerant(input)?;
    let mut header = reader.read_header()?;
    reject_coordinate_sorted(&header)?;
    let mut raw_header =
        crate::header_text::read_raw_header_text_with_format(input, Exact::Sam).ok();
    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
        if let Some(raw) = raw_header.as_mut() {
            *raw = crate::pg::add_samtools_pg(raw, argv).map_err(io::Error::other)?;
        }
    }
    let mut sink = open_output(output, fmt, &header, raw_header.as_deref())?;
    let mut template = Vec::new();
    let mut next = RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut next)?;
        if n == 0 {
            break;
        }
        sanitize_record(&header, &mut next, settings.sanitize_flags);
        write_or_buffer_template(
            &header,
            sink.as_mut(),
            &mut template,
            next.clone(),
            settings,
        )?;
    }
    flush_template(&header, sink.as_mut(), &mut template, settings)?;
    sink.finish()
}

fn run_fixmate_stdin(
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    settings: FixmateSettings,
) -> io::Result<()> {
    let mut bytes = Vec::new();
    io::stdin().read_to_end(&mut bytes)?;

    if bytes.first() == Some(&b'@') {
        return run_fixmate_sam_bytes(bytes, output, fmt, pg_argv, settings);
    }

    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        "stdin fixmate input must be SAM text",
    ))
}

fn run_fixmate_sam_bytes(
    bytes: Vec<u8>,
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    settings: FixmateSettings,
) -> io::Result<()> {
    let raw_header = raw_sam_header_text_from_bytes(&bytes)?;
    let normalized = crate::sam_compat::normalize_sam_aux_int_types(&bytes);
    let reader = sam::io::Reader::new(BufReader::new(Cursor::new(normalized)));
    run_fixmate_sam_reader(reader, Some(raw_header), output, fmt, pg_argv, settings)
}

fn run_fixmate_sam_reader<R: BufRead>(
    mut reader: sam::io::Reader<R>,
    mut raw_header: Option<String>,
    output: Option<&Path>,
    fmt: OutFmt,
    pg_argv: Option<&[OsString]>,
    settings: FixmateSettings,
) -> io::Result<()> {
    let mut header = reader.read_header()?;
    reject_coordinate_sorted(&header)?;
    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
        if let Some(raw) = raw_header.as_mut() {
            *raw = crate::pg::add_samtools_pg(raw, argv).map_err(io::Error::other)?;
        }
    }
    let mut sink = open_output(output, fmt, &header, raw_header.as_deref())?;
    let mut template = Vec::new();
    let mut next = RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut next)?;
        if n == 0 {
            break;
        }
        sanitize_record(&header, &mut next, settings.sanitize_flags);
        write_or_buffer_template(
            &header,
            sink.as_mut(),
            &mut template,
            next.clone(),
            settings,
        )?;
    }
    flush_template(&header, sink.as_mut(), &mut template, settings)?;
    sink.finish()
}

fn raw_sam_header_text_from_bytes(bytes: &[u8]) -> io::Result<String> {
    let text =
        std::str::from_utf8(bytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut raw = String::new();
    for line in text.split_inclusive('\n') {
        if line.starts_with('@') {
            raw.push_str(line);
        } else {
            break;
        }
    }
    Ok(raw)
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

fn write_or_buffer_template(
    header: &sam::Header,
    sink: &mut dyn Sink,
    template: &mut Vec<RecordBuf>,
    next: RecordBuf,
    settings: FixmateSettings,
) -> io::Result<()> {
    if template.is_empty() || same_template(template.first().unwrap(), &next) {
        template.push(next);
    } else {
        flush_template(header, sink, template, settings)?;
        template.push(next);
    }
    Ok(())
}

fn same_template(a: &RecordBuf, b: &RecordBuf) -> bool {
    let a_name = a.name().map(|n| n.to_vec());
    let b_name = b.name().map(|n| n.to_vec());
    a_name == b_name && b_name.is_some()
}

fn flush_template(
    header: &sam::Header,
    sink: &mut dyn Sink,
    template: &mut Vec<RecordBuf>,
    settings: FixmateSettings,
) -> io::Result<()> {
    if template.is_empty() {
        return Ok(());
    }

    fix_template(template, settings);
    let template_length_overrides = template_length_overrides(template);

    for (rec, template_length) in template.drain(..).zip(template_length_overrides) {
        if !skip_for_remove_reads(&rec, settings.remove_reads) {
            if let Some(template_length) = template_length {
                sink.write_record_with_template_length(header, &rec, template_length)?;
            } else {
                sink.write_record(header, &rec)?;
            }
        }
    }

    Ok(())
}

fn fix_template(template: &mut [RecordBuf], settings: FixmateSettings) {
    if settings.base_mods {
        fix_template_base_mods(template);
    }

    let primary: Vec<_> = template
        .iter()
        .enumerate()
        .filter_map(|(i, rec)| {
            let flags = rec.flags();
            (!flags.is_secondary() && !flags.is_supplementary()).then_some(i)
        })
        .collect();

    if primary.len() >= 2 {
        let i = primary[0];
        let j = primary[1];
        let (a, b) = pair_fixmate(
            template[i].clone(),
            template[j].clone(),
            settings.remove_reads,
            settings.mate_score,
            settings.add_template_cigar,
        );
        template[i] = a;
        template[j] = b;
    }
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

#[derive(Clone, Debug)]
struct BaseModPrimary {
    seq: Vec<u8>,
    len: usize,
    reverse: bool,
}

fn fix_template_base_mods(template: &mut [RecordBuf]) {
    let mut primary: [Option<BaseModPrimary>; 2] = [None, None];

    for rec in template.iter_mut() {
        if rec.flags().is_secondary() || rec.flags().is_supplementary() {
            continue;
        }

        let slot = read_segment_slot(rec);
        fix_base_mod_record(None, rec);
        primary[slot] = Some(BaseModPrimary {
            seq: rec.sequence().as_ref().to_vec(),
            len: record_query_len(rec),
            reverse: rec.flags().is_reverse_complemented(),
        });
    }

    for rec in template.iter_mut() {
        if !(rec.flags().is_secondary() || rec.flags().is_supplementary()) {
            continue;
        }

        let slot = read_segment_slot(rec);
        fix_base_mod_record(primary[slot].as_ref(), rec);
    }
}

fn read_segment_slot(rec: &RecordBuf) -> usize {
    usize::from(rec.flags().is_last_segment())
}

fn fix_base_mod_record(primary: Option<&BaseModPrimary>, rec: &mut RecordBuf) {
    normalize_base_mod_tag_names(rec);

    if primary.is_none() && (rec.flags().is_secondary() || rec.flags().is_supplementary()) {
        if has_tag(rec, [b'M', b'M']) {
            if rec.sequence().is_empty() {
                validate_or_delete_mod_tags(rec);
            } else {
                delete_mod_tags(rec);
            }
        } else {
            delete_mod_tags(rec);
        }
        return;
    }

    let Some(mm) = string_tag(rec, [b'M', b'M']) else {
        delete_mod_tags(rec);
        return;
    };

    let seq_len = record_query_len(rec);
    let mn = int_tag(rec, [b'M', b'N']);
    let (end5, end3) = hard_clip_ends(rec);

    match primary {
        None => {
            if end5 == 0 && end3 == 0 && mn.is_none_or(|n| n < 0) {
                if seq_len > 0 {
                    set_mn_append(rec, seq_len);
                }
            } else if (end5 != 0 || end3 != 0) && mn != Some(seq_len as i64) {
                delete_mod_tags(rec);
                return;
            }
        }
        Some(primary) => {
            if mn == Some(seq_len as i64) {
                if rec.sequence().is_empty() && (end5 != 0 || end3 != 0) {
                    delete_mod_tags(rec);
                    return;
                }
                validate_or_delete_mod_tags(rec);
                return;
            }

            if primary.len != seq_len + end5 + end3 {
                delete_mod_tags(rec);
                return;
            }

            if (end5 != 0 || end3 != 0) && mn.is_none_or(|n| n < 0 || n as usize == primary.len) {
                if primary.seq.is_empty() {
                    delete_mod_tags(rec);
                    return;
                }
                let ml = ml_tag(rec);
                match trim_base_mods(
                    &primary.seq,
                    primary.reverse,
                    &mm,
                    ml.as_deref(),
                    end5,
                    end3,
                ) {
                    Some((new_mm, new_ml)) => {
                        set_string_preserve_order(rec, [b'M', b'M'], new_mm);
                        if let Some(new_ml) = new_ml {
                            set_ml_preserve_order(rec, new_ml);
                        }
                    }
                    None => {
                        delete_mod_tags(rec);
                        return;
                    }
                }
            }

            if seq_len > 0 {
                set_mn_append(rec, seq_len);
            }
        }
    }

    validate_or_delete_mod_tags(rec);
}

fn normalize_base_mod_tag_names(rec: &mut RecordBuf) {
    let mm = Tag::from([b'M', b'M']);
    let ml = Tag::from([b'M', b'L']);
    let draft_mm = Tag::from([b'M', b'm']);
    let draft_ml = Tag::from([b'M', b'l']);
    let fields: Vec<_> = rec
        .data()
        .iter()
        .map(|(tag, value)| {
            let tag = if tag == draft_mm {
                mm
            } else if tag == draft_ml {
                ml
            } else {
                tag
            };
            (tag, value.clone())
        })
        .collect();
    *rec.data_mut() = fields.into_iter().collect();
}

fn hard_clip_ends(rec: &RecordBuf) -> (usize, usize) {
    let ops = rec.cigar().as_ref();
    let left = ops
        .first()
        .filter(|op| op.kind() == Kind::HardClip)
        .map_or(0, |op| op.len());
    let right = ops
        .last()
        .filter(|op| op.kind() == Kind::HardClip)
        .map_or(0, |op| op.len());

    if rec.flags().is_reverse_complemented() {
        (right, left)
    } else {
        (left, right)
    }
}

fn record_query_len(rec: &RecordBuf) -> usize {
    let seq_len = rec.sequence().len();
    if seq_len > 0 {
        return seq_len;
    }

    rec.cigar()
        .as_ref()
        .iter()
        .filter(|op| {
            matches!(
                op.kind(),
                Kind::Match
                    | Kind::Insertion
                    | Kind::SoftClip
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
            )
        })
        .map(|op| op.len())
        .sum()
}

fn string_tag(rec: &RecordBuf, tag: [u8; 2]) -> Option<String> {
    match rec.data().get(&Tag::from(tag)) {
        Some(Value::String(s)) => Some(String::from_utf8_lossy(s).into_owned()),
        _ => None,
    }
}

fn int_tag(rec: &RecordBuf, tag: [u8; 2]) -> Option<i64> {
    rec.data().get(&Tag::from(tag)).and_then(Value::as_int)
}

fn ml_tag(rec: &RecordBuf) -> Option<Vec<u8>> {
    match rec.data().get(&Tag::from([b'M', b'L'])) {
        Some(Value::Array(Array::UInt8(values))) => Some(values.clone()),
        _ => None,
    }
}

fn has_tag(rec: &RecordBuf, tag: [u8; 2]) -> bool {
    rec.data().get(&Tag::from(tag)).is_some()
}

fn set_string_preserve_order(rec: &mut RecordBuf, tag: [u8; 2], value: String) {
    rec.data_mut()
        .insert(Tag::from(tag), Value::String(BString::from(value)));
}

fn set_ml_preserve_order(rec: &mut RecordBuf, values: Vec<u8>) {
    rec.data_mut()
        .insert(Tag::from([b'M', b'L']), Value::Array(Array::UInt8(values)));
}

fn set_mn_append(rec: &mut RecordBuf, len: usize) {
    aux_set_append(rec, Tag::from([b'M', b'N']), Value::from(len as u32));
}

fn delete_mod_tags(rec: &mut RecordBuf) {
    let mm = Tag::from([b'M', b'M']);
    let draft_mm = Tag::from([b'M', b'm']);
    let ml = Tag::from([b'M', b'L']);
    let draft_ml = Tag::from([b'M', b'l']);
    let mn = Tag::from([b'M', b'N']);
    let kept: Vec<_> = rec
        .data()
        .iter()
        .filter(|(tag, _)| {
            *tag != mm && *tag != draft_mm && *tag != ml && *tag != draft_ml && *tag != mn
        })
        .map(|(tag, value)| (tag, value.clone()))
        .collect();
    *rec.data_mut() = kept.into_iter().collect();
}

fn validate_or_delete_mod_tags(rec: &mut RecordBuf) {
    let Some(mm) = string_tag(rec, [b'M', b'M']) else {
        delete_mod_tags(rec);
        return;
    };

    match validate_mm(
        rec.sequence().as_ref(),
        rec.flags().is_reverse_complemented(),
        &mm,
        ml_tag(rec).as_deref(),
    ) {
        true => {}
        false => delete_mod_tags(rec),
    }
}

fn validate_mm(seq: &[u8], reverse: bool, mm: &str, ml: Option<&[u8]>) -> bool {
    let mut total_mods = 0usize;
    for group in mm.split_terminator(';') {
        let Some((fundamental, deltas)) = parse_mm_group(group) else {
            return false;
        };
        let fundamental = oriented_base(fundamental, reverse);
        let base_count = count_base(seq, fundamental);
        let mut n = 0usize;
        for delta in deltas {
            if !seq.is_empty() && base_count.saturating_sub(n) <= delta {
                return false;
            }
            n += delta + 1;
            total_mods += 1;
        }
    }

    ml.is_none_or(|values| values.len() == total_mods)
}

fn trim_base_mods(
    primary_seq: &[u8],
    primary_reverse: bool,
    mm: &str,
    ml: Option<&[u8]>,
    end5: usize,
    end3: usize,
) -> Option<(String, Option<Vec<u8>>)> {
    let mut new_mm = String::new();
    let mut new_ml = Vec::new();
    let mut ml_index = 0usize;

    for group in mm.split_terminator(';') {
        let (fundamental, deltas) = parse_mm_group(group)?;
        let fundamental = oriented_base(fundamental, primary_reverse);
        let counts5 = count_base(&primary_seq[..end5.min(primary_seq.len())], fundamental);
        let counts3_end = primary_seq.len().saturating_sub(end3);
        let counts3 = count_base(&primary_seq[..counts3_end], fundamental);
        let header = group
            .split_once(',')
            .map(|(header, _)| header)
            .unwrap_or(group);

        let mut old_pos = 0usize;
        let mut kept_positions = Vec::new();
        let mut kept_ml = Vec::new();
        for delta in deltas {
            old_pos += delta;
            let keep = old_pos >= counts5 && old_pos < counts3;
            if keep {
                kept_positions.push(old_pos - counts5);
                if let Some(ml) = ml {
                    kept_ml.push(*ml.get(ml_index)?);
                }
            }
            old_pos += 1;
            ml_index += 1;
        }

        new_mm.push_str(header);
        if kept_positions.is_empty() {
            new_mm.push(';');
        } else {
            new_mm.push(',');
            let mut prev = None;
            for (i, pos) in kept_positions.into_iter().enumerate() {
                if i > 0 {
                    new_mm.push(',');
                }
                let delta = match prev {
                    Some(prev) => pos - prev - 1,
                    None => pos,
                };
                new_mm.push_str(&delta.to_string());
                prev = Some(pos);
            }
            new_mm.push(';');
        }
        new_ml.extend(kept_ml);
    }

    if let Some(ml) = ml
        && ml_index != ml.len()
    {
        return None;
    }

    Some((new_mm, ml.map(|_| new_ml)))
}

fn parse_mm_group(group: &str) -> Option<(u8, Vec<usize>)> {
    let fundamental = group.as_bytes().first().copied()?.to_ascii_uppercase();
    let Some((_, rest)) = group.split_once(',') else {
        return Some((fundamental, Vec::new()));
    };
    let mut deltas = Vec::new();
    for raw in rest.split(',') {
        deltas.push(raw.parse().ok()?);
    }
    Some((fundamental, deltas))
}

fn oriented_base(base: u8, reverse: bool) -> u8 {
    if reverse { complement_base(base) } else { base }
}

fn complement_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b => b,
    }
}

fn count_base(seq: &[u8], base: u8) -> usize {
    seq.iter()
        .filter(|b| b.to_ascii_uppercase() == base)
        .count()
}

fn update_template_lengths(a: &mut RecordBuf, b: &mut RecordBuf) {
    match template_lengths(a, b) {
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

fn template_length_overrides(template: &[RecordBuf]) -> Vec<Option<i64>> {
    let mut overrides = vec![None; template.len()];
    let primary: Vec<_> = template
        .iter()
        .enumerate()
        .filter_map(|(i, rec)| {
            let flags = rec.flags();
            (!flags.is_secondary() && !flags.is_supplementary()).then_some(i)
        })
        .collect();

    if primary.len() >= 2 {
        let i = primary[0];
        let j = primary[1];
        if let Some((a_tlen, b_tlen)) = template_lengths(&template[i], &template[j]) {
            if i32::try_from(a_tlen).is_err() {
                overrides[i] = Some(a_tlen);
            }
            if i32::try_from(b_tlen).is_err() {
                overrides[j] = Some(b_tlen);
            }
        }
    }

    overrides
}

fn template_lengths(a: &RecordBuf, b: &RecordBuf) -> Option<(i64, i64)> {
    if a.flags().is_unmapped()
        || b.flags().is_unmapped()
        || a.reference_sequence_id().is_none()
        || a.reference_sequence_id() != b.reference_sequence_id()
    {
        return None;
    }

    match (five_prime_position(a), five_prime_position(b)) {
        (Some(a5), Some(b5)) => Some((b5 - a5, a5 - b5)),
        _ => None,
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
    // Order-preserving (noodles' remove is a swap_remove that would
    // reorder the surrounding MC/MQ tags).
    aux_del(a, ct_tag);
    aux_del(b, ct_tag);

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

/// Order-preserving `bam_aux_del`: drop `tag`, keep the rest in order
/// (noodles' `Data::remove` is a `swap_remove` that reorders).
fn aux_del(target: &mut RecordBuf, tag: Tag) {
    let kept: Vec<_> = target
        .data()
        .iter()
        .filter(|(t, _)| *t != tag)
        .map(|(t, v)| (t, v.clone()))
        .collect();
    *target.data_mut() = kept.into_iter().collect();
}

/// Mirrors HTSlib's `bam_aux_del` + `bam_aux_append`: remove any existing
/// `tag` (preserving the order of the others) then append the new value at
/// the end — so an updated tag moves to the tail, like upstream.
fn aux_set_append(target: &mut RecordBuf, tag: Tag, value: Value) {
    let mut fields: Vec<_> = target
        .data()
        .iter()
        .filter(|(t, _)| *t != tag)
        .map(|(t, v)| (t, v.clone()))
        .collect();
    fields.push((tag, value));
    *target.data_mut() = fields.into_iter().collect();
}

fn update_mate_aux_tags(target: &mut RecordBuf, mate: &RecordBuf) {
    let mc_tag = Tag::from([b'M', b'C']);
    let mq_tag = Tag::from([b'M', b'Q']);

    if mate.flags().is_unmapped() {
        aux_del(target, mq_tag);
        if target.flags().is_unmapped() {
            aux_del(target, mc_tag);
        } else {
            // bam_mate.c:197 — MC is added when *either* read is mapped;
            // an empty mate CIGAR formats as `*` (→ `MC:Z:*`).
            aux_set_append(
                target,
                mc_tag,
                Value::String(BString::from(format_cigar(mate.cigar()))),
            );
        }
        return;
    }

    // bam_mate.c:188-207 — del-then-append (moves the tag to the tail),
    // MQ before MC.
    if let Some(mapping_quality) = mate.mapping_quality() {
        aux_set_append(target, mq_tag, Value::from(mapping_quality.get()));
    } else {
        aux_del(target, mq_tag);
    }

    aux_set_append(
        target,
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

    fn write_record_with_template_length(
        &mut self,
        header: &sam::Header,
        record: &RecordBuf,
        template_length: i64,
    ) -> io::Result<()> {
        let _ = template_length;
        self.write_record(header, record)
    }

    fn finish(self: Box<Self>) -> io::Result<()> {
        Ok(())
    }
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);
struct SamFile(File);
struct SamStdout(io::Stdout);
struct CramOut<W: Write> {
    writer: bam::io::Writer<bgzf::io::Writer<File>>,
    tmp_bam_path: crate::tmp_file::TempPath,
    reference: PathBuf,
    out: W,
}

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

    fn write_record_with_template_length(
        &mut self,
        header: &sam::Header,
        record: &RecordBuf,
        template_length: i64,
    ) -> io::Result<()> {
        crate::sam_render::write_record_with_template_length(
            &mut self.0,
            header,
            record,
            template_length,
        )
    }
}
impl Sink for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        crate::sam_render::write_record(&mut self.0, header, record)
    }

    fn write_record_with_template_length(
        &mut self,
        header: &sam::Header,
        record: &RecordBuf,
        template_length: i64,
    ) -> io::Result<()> {
        crate::sam_render::write_record_with_template_length(
            &mut self.0,
            header,
            record,
            template_length,
        )
    }
}
impl<W: Write> Sink for CramOut<W> {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        use sam::alignment::io::Write as _;
        self.writer.write_alignment_record(header, record)
    }

    fn finish(self: Box<Self>) -> io::Result<()> {
        let this = *self;
        let CramOut {
            writer,
            tmp_bam_path,
            reference,
            out,
        } = this;
        drop(writer);

        let result = htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference(
            tmp_bam_path.path(),
            &reference,
            out,
        )
        .map(|_| ());

        tmp_bam_path.close().ok();
        result
    }
}

fn open_output(
    out: Option<&Path>,
    fmt: OutFmt,
    header: &sam::Header,
    raw_header: Option<&str>,
) -> io::Result<Box<dyn Sink>> {
    match (out, fmt) {
        (Some(p), OutFmt::Sam) => {
            let mut w = File::create(p)?;
            if let Some(raw) = raw_header {
                w.write_all(raw.as_bytes())?;
            } else {
                crate::sam_render::write_header(&mut w, header)?;
            }
            Ok(Box::new(SamFile(w)))
        }
        (Some(p), OutFmt::Bam) => {
            let mut w = bam::io::Writer::new(File::create(p)?);
            w.write_header(header)?;
            Ok(Box::new(BamFile(w)))
        }
        (Some(p), OutFmt::Cram) => {
            let out = File::create(p)?;
            open_cram_output(header, out)
        }
        (None, OutFmt::Sam) => {
            let mut w = io::stdout();
            if let Some(raw) = raw_header {
                w.write_all(raw.as_bytes())?;
            } else {
                crate::sam_render::write_header(&mut w, header)?;
            }
            Ok(Box::new(SamStdout(w)))
        }
        (None, OutFmt::Bam) => {
            let mut w = bam::io::Writer::new(io::stdout());
            w.write_header(header)?;
            Ok(Box::new(BamStdout(w)))
        }
        (None, OutFmt::Cram) => open_cram_output(header, io::stdout()),
    }
}

fn open_cram_output<W>(header: &sam::Header, out: W) -> io::Result<Box<dyn Sink>>
where
    W: Write + 'static,
{
    let reference = current_global_args().reference.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM output requires top-level --reference FILE",
        )
    })?;
    crate::reference::ensure_fai_index(&reference, None)?;

    let (tmp_bam_file, tmp_bam_path) = crate::tmp_file::create_temp_file("fixmate", Some("bam"))?;
    let mut writer = bam::io::Writer::new(tmp_bam_file);
    writer.write_header(header)?;

    Ok(Box::new(CramOut {
        writer,
        tmp_bam_path,
        reference,
        out,
    }))
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools fixmate [options] <in.bam> [<out.bam>]")?;
    writeln!(w, "  -O sam|bam|cram   output format (default: bam)")?;
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
