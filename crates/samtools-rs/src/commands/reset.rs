//! `samtools reset` — revert aligner changes in reads.
//!
//! Mirrors `main_reset` in `reset.c`. Initial Rust port operates on BAM
//! input and writes BAM output via `RecordBuf` mutation:
//!  - clears `reference_sequence_id`, `alignment_start`, `cigar`,
//!    `mapping_quality`, `mate_reference_sequence_id`,
//!    `mate_alignment_start`, `template_length`
//!  - clears flag bits that depend on alignment (FUNMAP set to 1,
//!    FSECONDARY/FSUPPLEMENTARY/FPROPER_PAIR/FMUNMAP/FREVERSE/FMREVERSE
//!    cleared)
//!  - drops common aligner-added aux tags (NM, MD, AS, XS, SA, MC, MQ,
//!    NH, HI) by default
//!
//! `--dupflag` preserves duplicate flags, SAM stdin input works, reverse-strand
//! sequence/quality re-reversal is implemented, legacy SAM `@HD VN:1` input is
//! accepted, `--no-RG` removes read-group headers and tags, and `--reject-PG` /
//! `--no-PG` remove program header chains. `-x`/`--keep-tag` aux-tag filtering
//! is honored.

use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{
    self,
    alignment::{
        RecordBuf,
        record::{Flags, MappingQuality},
    },
};

use crate::aux_list::{AuxTag, parse_aux_list};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

const DEFAULT_DROP_TAGS: &[&[u8; 2]] = &[
    b"NM", b"MD", b"AS", b"XS", b"SA", b"MC", b"MQ", b"NH", b"HI", b"ms",
];

/// Entry point for `samtools reset`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut output: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let mut output_fmt = OutFmt::Bam;
    let mut fmt_explicit = false;
    let mut extra_drop: Vec<AuxTag> = Vec::new();
    let mut keep_only: Option<HashSet<AuxTag>> = None;
    let mut preserve_duplicate = false;
    let mut remove_read_groups = false;
    let mut no_pg = false;
    let mut reject_programs = Vec::new();

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-O" | "--output-fmt" => {
                let v = iter.next().and_then(|a| a.to_str()).unwrap_or("bam");
                output_fmt = match v.to_lowercase().as_str() {
                    "sam" => OutFmt::Sam,
                    "bam" => OutFmt::Bam,
                    _ => OutFmt::Bam,
                };
                fmt_explicit = true;
            }
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "-x" | "--remove-tag" | "--remove-tags" => {
                if let Some(v) = iter.next().and_then(|a| a.to_str()) {
                    if let Some(rest) = v.strip_prefix('^') {
                        match parse_aux_list(rest) {
                            Ok(tags) => merge_keep_tags(&mut keep_only, tags),
                            Err(e) => {
                                print_error("reset", format!("invalid -x value \"{rest}\": {e}"));
                                return ExitCode::from(1);
                            }
                        }
                    } else {
                        match parse_aux_list(v) {
                            Ok(tags) => extra_drop.extend(tags),
                            Err(e) => {
                                print_error("reset", format!("invalid -x value \"{v}\": {e}"));
                                return ExitCode::from(1);
                            }
                        }
                    }
                }
            }
            "--keep-tag" | "--keep-tags" => {
                if let Some(v) = iter.next().and_then(|a| a.to_str()) {
                    match parse_aux_list(v) {
                        Ok(tags) => merge_keep_tags(&mut keep_only, tags),
                        Err(e) => {
                            print_error("reset", format!("invalid --keep-tag value \"{v}\": {e}"));
                            return ExitCode::from(1);
                        }
                    }
                }
            }
            "--dupflag" => {
                preserve_duplicate = true;
            }
            "--no-RG" => {
                remove_read_groups = true;
            }
            "--reject-PG" => {
                if let Some(v) = iter.next().and_then(|a| a.to_str()) {
                    reject_programs.push(v.to_string());
                }
            }
            "--no-PG" => {
                no_pg = true;
            }
            "-T" => {
                if matches!(s, "-T") {
                    let _ = iter.next();
                }
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("reset", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                }
            }
        }
    }

    let settings = ResetSettings {
        extra_drop: &extra_drop,
        keep_only: keep_only.as_ref(),
        preserve_duplicate,
        remove_read_groups,
        no_pg,
        reject_programs: &reject_programs,
        pg_argv: Some(args),
    };

    // Upstream `sam_open_mode`: infer the format from the `-o` filename
    // extension when `--output-fmt` was not given.
    if !fmt_explicit && let Some(p) = output.as_deref().and_then(|p| p.to_str()) {
        // `sam_open_mode`: BAM/CRAM only for those extensions; otherwise
        // (incl. `.sam` and no/unknown extension) plain SAM text.
        output_fmt = if p.ends_with(".bam") {
            OutFmt::Bam
        } else {
            OutFmt::Sam
        };
    }

    let result = match input {
        Some(input) if input != Path::new("-") => {
            let format = match sam_io::sam_open_format(&input) {
                Ok(f) => f,
                Err(e) => {
                    print_error("reset", e.to_string());
                    return ExitCode::from(1);
                }
            };
            if !matches!(format.exact, Exact::Sam | Exact::Bam) {
                print_error(
                    "reset",
                    "only SAM and BAM input are currently supported (CRAM TODO)",
                );
                return ExitCode::from(1);
            }

            run_reset(
                &input,
                format.exact,
                output.as_deref(),
                output_fmt,
                &settings,
            )
        }
        _ => run_reset_stdin(output.as_deref(), output_fmt, &settings),
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("reset", "reset failed", &e);
            ExitCode::from(1)
        }
    }
}

fn merge_keep_tags(keep_only: &mut Option<HashSet<AuxTag>>, tags: HashSet<AuxTag>) {
    keep_only.get_or_insert_with(HashSet::new).extend(tags);
}

#[derive(Clone, Copy)]
enum OutFmt {
    Sam,
    Bam,
}

struct ResetSettings<'a> {
    extra_drop: &'a [[u8; 2]],
    keep_only: Option<&'a HashSet<[u8; 2]>>,
    preserve_duplicate: bool,
    remove_read_groups: bool,
    no_pg: bool,
    reject_programs: &'a [String],
    pg_argv: Option<&'a [OsString]>,
}

fn run_reset(
    input: &Path,
    input_format: Exact,
    output: Option<&Path>,
    fmt: OutFmt,
    settings: &ResetSettings<'_>,
) -> io::Result<()> {
    match input_format {
        Exact::Sam => run_reset_sam(input, output, fmt, settings),
        Exact::Bam => run_reset_bam(input, output, fmt, settings),
        _ => unreachable!("input format checked by caller"),
    }
}

/// Upstream `reset` output header (`reset.c:307-324`): a **fresh**
/// header — `@HD\tVN:1.6` (`sam_hdr_init` + `SAM_FORMAT_VERSION`),
/// then the input's `@RG` lines verbatim unless `--no-RG`
/// (`getRGlines`), then the input's `@PG` lines verbatim with the
/// `--reject-PG` cut (`getPGlines`) plus the samtools `@PG`. `@SQ`,
/// `@CO`, the original `@HD`, and any other lines are dropped (all
/// reads become unmapped). Returns `None` if the raw header is
/// unavailable.
fn reset_raw_header(path: &Path, settings: &ResetSettings<'_>) -> Option<String> {
    let raw = crate::header_text::read_raw_header_text(path).ok()?;

    let field = |line: &str, tag: &str| -> Option<String> {
        line.split('\t')
            .skip(1)
            .find_map(|f| f.strip_prefix(tag).map(|v| v.to_string()))
    };
    let rejected: std::collections::HashSet<String> =
        settings.reject_programs.iter().cloned().collect();

    // `--reject-PG ID` (`reset.c:223`): iterate @PG lines in header
    // order; keep them until the first whose `ID` matches, then drop
    // that one and **all subsequent @PG** ("from this PG onwards").
    let mut pg_cut = false;
    let mut lines: Vec<String> = vec!["@HD\tVN:1.6".to_string()];
    for line in raw.lines() {
        if line.is_empty() {
            continue;
        }
        if line.starts_with("@RG") {
            if !settings.remove_read_groups {
                lines.push(line.to_string());
            }
        } else if line.starts_with("@PG") {
            if pg_cut {
                continue;
            }
            if let Some(id) = field(line, "ID:")
                && rejected.contains(&id)
            {
                pg_cut = true;
                continue;
            }
            lines.push(line.to_string());
        }
        // @HD / @SQ / @CO / anything else: dropped (fresh header).
    }
    let mut text = lines.join("\n");
    text.push('\n');
    if !settings.no_pg
        && let Some(argv) = settings.pg_argv
    {
        text = crate::pg::add_samtools_pg(&text, argv).ok()?;
    }
    Some(text)
}

fn run_reset_bam(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    settings: &ResetSettings<'_>,
) -> io::Result<()> {
    // Upstream rebuilds a fresh @HD/@RG/@PG header for *all* output
    // formats (no @SQ); use it for BAM too so binary output drops @SQ.
    let raw = reset_raw_header(input, settings);
    let mut reader = bam::io::Reader::new(File::open(input)?);
    run_reset_bam_reader(&mut reader, output, fmt, settings, raw.as_deref())
}

fn run_reset_bam_reader<R>(
    reader: &mut bam::io::Reader<R>,
    output: Option<&Path>,
    fmt: OutFmt,
    settings: &ResetSettings<'_>,
    raw_header: Option<&str>,
) -> io::Result<()>
where
    R: Read,
{
    let header = reader.read_header()?;
    let mut output_header = header.clone();
    reset_header(&mut output_header, settings)?;
    let mut sink = open_output(output, fmt, &output_header, raw_header)?;

    let mut record = RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut record)?;
        if n == 0 {
            break;
        }
        reset_record(&mut record, settings);
        sink.write_record(&output_header, &record)?;
    }
    Ok(())
}

fn run_reset_sam(
    input: &Path,
    output: Option<&Path>,
    fmt: OutFmt,
    settings: &ResetSettings<'_>,
) -> io::Result<()> {
    // Upstream rebuilds a fresh @HD/@RG/@PG header for *all* output
    // formats (no @SQ); use it for BAM too so binary output drops @SQ.
    let raw = reset_raw_header(input, settings);
    let bytes = normalize_legacy_sam_header_version(std::fs::read(input)?);
    let bytes = crate::sam_compat::normalize_sam_aux_int_types(&bytes);
    let mut reader = sam::io::Reader::new(BufReader::new(Cursor::new(bytes)));
    run_reset_sam_reader(&mut reader, output, fmt, settings, raw.as_deref())
}

fn run_reset_stdin(
    output: Option<&Path>,
    fmt: OutFmt,
    settings: &ResetSettings<'_>,
) -> io::Result<()> {
    let stdin = io::stdin();
    let mut input = Vec::new();
    stdin.lock().read_to_end(&mut input)?;

    if !looks_like_sam(&input) {
        let mut reader = bam::io::Reader::new(Cursor::new(input));
        return run_reset_bam_reader(&mut reader, output, fmt, settings, None);
    }

    let input = normalize_legacy_sam_header_version(input);
    let input = crate::sam_compat::normalize_sam_aux_int_types(&input);
    let mut reader = sam::io::Reader::new(BufReader::new(Cursor::new(input)));
    run_reset_sam_reader(&mut reader, output, fmt, settings, None)
}

fn looks_like_sam(input: &[u8]) -> bool {
    input
        .iter()
        .copied()
        .find(|b| !matches!(b, b' ' | b'\t' | b'\r' | b'\n'))
        == Some(b'@')
}

fn run_reset_sam_reader<R>(
    reader: &mut sam::io::Reader<R>,
    output: Option<&Path>,
    fmt: OutFmt,
    settings: &ResetSettings<'_>,
    raw_header: Option<&str>,
) -> io::Result<()>
where
    R: BufRead,
{
    let header = reader.read_header()?;
    let mut output_header = header.clone();
    reset_header(&mut output_header, settings)?;
    let mut sink = open_output(output, fmt, &output_header, raw_header)?;

    loop {
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        reset_record(&mut record, settings);
        sink.write_record(&output_header, &record)?;
    }
    Ok(())
}

fn reset_header(header: &mut sam::Header, settings: &ResetSettings<'_>) -> io::Result<()> {
    *header.header_mut() = Some(Default::default());
    header.reference_sequences_mut().clear();
    header.comments_mut().clear();

    if settings.remove_read_groups {
        header.read_groups_mut().clear();
    }

    if !settings.reject_programs.is_empty() {
        reject_header_programs(header, settings.reject_programs);
    }

    if !settings.no_pg
        && let Some(argv) = settings.pg_argv
    {
        *header = crate::pg::add_samtools_pg_to_header(header, argv)?;
    }
    Ok(())
}

fn reject_header_programs(header: &mut sam::Header, rejected_ids: &[String]) {
    use htslib_rs::sam::header::record::value::map::program::tag;

    let mut rejected: HashSet<Vec<u8>> = rejected_ids
        .iter()
        .map(|id| id.as_bytes().to_vec())
        .collect();

    loop {
        let mut changed = false;
        for (id, program) in header.programs().as_ref() {
            let id_bytes: &[u8] = id.as_ref();
            if rejected.contains(id_bytes) {
                continue;
            }

            if let Some(previous_id) = program.other_fields().get(&tag::PREVIOUS_PROGRAM_ID) {
                let previous_id: &[u8] = previous_id.as_ref();
                if rejected.contains(previous_id) {
                    rejected.insert(id_bytes.to_vec());
                    changed = true;
                }
            }
        }

        if !changed {
            break;
        }
    }

    header
        .programs_mut()
        .as_mut()
        .retain(|id, _| !rejected.contains(id.as_ref() as &[u8]));
}

fn normalize_legacy_sam_header_version(mut input: Vec<u8>) -> Vec<u8> {
    let line_end = input
        .iter()
        .position(|&b| b == b'\n')
        .map(|i| i + 1)
        .unwrap_or(input.len());

    if !input[..line_end].starts_with(b"@HD\t") {
        return input;
    }

    for i in 0..line_end.saturating_sub(5) {
        if !input[i..].starts_with(b"\tVN:1") {
            continue;
        }

        if matches!(input.get(i + 5), None | Some(b'\t' | b'\r' | b'\n')) {
            input.splice(i + 5..i + 5, b".0".iter().copied());
            break;
        }
    }

    input
}

fn reset_record(record: &mut RecordBuf, settings: &ResetSettings<'_>) {
    let was_reverse = record.flags().is_reverse_complemented();

    // Reset alignment fields.
    *record.reference_sequence_id_mut() = None;
    *record.alignment_start_mut() = None;
    *record.cigar_mut() = sam::alignment::record_buf::Cigar::default();
    *record.mapping_quality_mut() = Some(MappingQuality::MIN);
    *record.mate_reference_sequence_id_mut() = None;
    *record.mate_alignment_start_mut() = None;
    *record.template_length_mut() = 0;

    // Reset flag bits.
    let mut flags = record.flags();
    flags.remove(Flags::PROPERLY_SEGMENTED);
    flags.remove(Flags::SECONDARY);
    flags.remove(Flags::SUPPLEMENTARY);
    if !settings.preserve_duplicate {
        flags.remove(Flags::DUPLICATE);
    }
    flags.remove(Flags::REVERSE_COMPLEMENTED);
    flags.remove(Flags::MATE_REVERSE_COMPLEMENTED);
    flags.insert(Flags::UNMAPPED);
    if flags.is_segmented() {
        flags.insert(Flags::MATE_UNMAPPED);
    } else {
        flags.remove(Flags::MATE_UNMAPPED);
    }
    *record.flags_mut() = flags;

    if was_reverse {
        reverse_complement_record_sequence(record);
        record.quality_scores_mut().as_mut().reverse();
    }

    // Drop aligner-added aux tags.
    let data = record.data_mut();
    let mut to_drop: HashSet<[u8; 2]> = HashSet::new();
    for tag in DEFAULT_DROP_TAGS {
        to_drop.insert(**tag);
    }
    for tag in settings.extra_drop {
        to_drop.insert(*tag);
    }
    if settings.remove_read_groups {
        to_drop.insert(*b"RG");
    }
    // noodles `Data::remove` is a swap-remove (order-breaking); rebuild
    // the field list preserving the original aux order (matches HTSlib).
    let kept: Vec<_> = data
        .iter()
        .filter(|(t, _)| {
            let bytes: [u8; 2] = (*t).into();
            let drop = (settings.remove_read_groups && bytes == *b"RG")
                || match settings.keep_only {
                    Some(keep) => !keep.contains(&bytes),
                    None => to_drop.contains(&bytes),
                };
            !drop
        })
        .map(|(t, v)| (t, v.clone()))
        .collect();
    *data = kept.into_iter().collect();
}

fn reverse_complement_record_sequence(record: &mut RecordBuf) {
    let sequence = record.sequence_mut().as_mut();
    sequence.reverse();
    for base in sequence {
        *base = complement_base(*base);
    }
}

fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'M' => b'K',
        b'R' => b'Y',
        b'W' => b'W',
        b'S' => b'S',
        b'Y' => b'R',
        b'K' => b'M',
        b'V' => b'B',
        b'H' => b'D',
        b'D' => b'H',
        b'B' => b'V',
        b'N' => b'N',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'm' => b'k',
        b'r' => b'y',
        b'w' => b'w',
        b's' => b's',
        b'y' => b'r',
        b'k' => b'm',
        b'v' => b'b',
        b'h' => b'd',
        b'd' => b'h',
        b'b' => b'v',
        b'n' => b'n',
        _ => base,
    }
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
            let hdr = bam_header(header, raw_header)?;
            let mut w = bam::io::Writer::new(File::create(p)?);
            w.write_header(&hdr)?;
            Ok(Box::new(BamFile(w)))
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
            let hdr = bam_header(header, raw_header)?;
            let mut w = bam::io::Writer::new(io::stdout());
            w.write_header(&hdr)?;
            Ok(Box::new(BamStdout(w)))
        }
    }
}

/// Header for BAM output: the rebuilt `@HD/@RG/@PG` text parsed back
/// into a structured header (so binary output also drops `@SQ`), or
/// the in-memory header when no raw header is available (stdin path).
fn bam_header(fallback: &sam::Header, raw_header: Option<&str>) -> io::Result<sam::Header> {
    match raw_header {
        Some(raw) => {
            let mut reader = sam::io::Reader::new(BufReader::new(Cursor::new(raw.as_bytes())));
            reader.read_header()
        }
        None => Ok(fallback.clone()),
    }
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools reset [options] [in.bam|in.sam|-]")?;
    writeln!(w, "  -o FILE                 output FILE")?;
    writeln!(w, "  -O sam|bam              output format")?;
    writeln!(
        w,
        "  -x/--remove-tag TAG     drop the listed aux tags (comma-separated, ^ for keep)"
    )?;
    writeln!(w, "  --keep-tag TAG          only keep the listed aux tags")?;
    writeln!(
        w,
        "  --reject-PG ID          remove program header chain from ID"
    )?;
    writeln!(
        w,
        "  --no-PG                 remove all program header records"
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::{BufReader, Cursor};

    use super::*;

    #[test]
    fn reset_sam_reader_supports_stdin_style_input() {
        let input = concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t99\tchr1\t2\t60\t4M\t=\t6\t8\tACGT\t!!!!\tNM:i:1\tMD:Z:3A\tRG:Z:g1\n",
        );
        let tmp = std::env::temp_dir().join(format!(
            "samtools-rs-reset-stdin-{}-{}.sam",
            std::process::id(),
            std::thread::current().name().unwrap_or("unnamed")
        ));
        let extra_drop = Vec::new();
        let settings = ResetSettings {
            extra_drop: &extra_drop,
            keep_only: None,
            preserve_duplicate: false,
            remove_read_groups: false,
            no_pg: true,
            reject_programs: &[],
            pg_argv: None,
        };
        let mut reader = sam::io::Reader::new(BufReader::new(Cursor::new(input.as_bytes())));

        run_reset_sam_reader(&mut reader, Some(&tmp), OutFmt::Sam, &settings, None).unwrap();

        let text = std::fs::read_to_string(&tmp).unwrap();
        let _ = std::fs::remove_file(&tmp);
        assert!(text.starts_with("@HD\tVN:1.6\n"));
        assert!(!text.contains("@SQ\t"));
        let record = text.lines().find(|line| !line.starts_with('@')).unwrap();
        let fields: Vec<_> = record.split('\t').collect();
        assert_eq!(fields[1], "77");
        assert_eq!(fields[2], "*");
        assert_eq!(fields[3], "0");
        assert_eq!(fields[4], "0");
        assert_eq!(fields[5], "*");
        assert!(!record.contains("\tNM:i:"));
        assert!(!record.contains("\tMD:Z:"));
        assert!(record.contains("\tRG:Z:g1"));
    }

    #[test]
    fn reset_bam_reader_supports_stdin_style_input() {
        let input = concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t99\tchr1\t2\t60\t4M\t=\t6\t8\tACGT\t!!!!\tNM:i:1\tMD:Z:3A\n",
        );
        let mut sam_reader = sam::io::Reader::new(BufReader::new(Cursor::new(input.as_bytes())));
        let header = sam_reader.read_header().unwrap();
        let mut bam_bytes = Vec::new();
        {
            let mut writer = bam::io::Writer::new(&mut bam_bytes);
            writer.write_header(&header).unwrap();
            for result in sam_reader.records() {
                let record = result.unwrap();
                use sam::alignment::io::Write as _;
                writer.write_alignment_record(&header, &record).unwrap();
            }
        }

        let tmp = std::env::temp_dir().join(format!(
            "samtools-rs-reset-bam-stdin-{}-{}.sam",
            std::process::id(),
            std::thread::current().name().unwrap_or("unnamed")
        ));
        let extra_drop = Vec::new();
        let settings = ResetSettings {
            extra_drop: &extra_drop,
            keep_only: None,
            preserve_duplicate: false,
            remove_read_groups: false,
            no_pg: true,
            reject_programs: &[],
            pg_argv: None,
        };
        let mut reader = bam::io::Reader::new(Cursor::new(bam_bytes));

        run_reset_bam_reader(&mut reader, Some(&tmp), OutFmt::Sam, &settings, None).unwrap();

        let text = std::fs::read_to_string(&tmp).unwrap();
        let _ = std::fs::remove_file(&tmp);
        let record = text.lines().find(|line| !line.starts_with('@')).unwrap();
        let fields: Vec<_> = record.split('\t').collect();
        assert_eq!(fields[1], "77");
        assert_eq!(fields[2], "*");
        assert_eq!(fields[3], "0");
        assert_eq!(fields[4], "0");
        assert_eq!(fields[5], "*");
        assert!(!record.contains("\tNM:i:"));
        assert!(!record.contains("\tMD:Z:"));
    }
}
