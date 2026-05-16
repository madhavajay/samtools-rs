//! `samtools reference` — reconstruct references from alignments and MD tags.
//!
//! This is a partial port of `reference.c`'s MD:Z path. It supports SAM and
//! BAM input, optional region output, and FASTA output. Embedded CRAM reference
//! extraction is intentionally deferred because it needs CRAM container/block
//! internals that are not exposed in this samtools-rs-only pass.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::core::Region;
use htslib_rs::format::Exact;
use htslib_rs::sam;

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

#[derive(Debug, Default)]
struct ReferenceOptions {
    output: Option<PathBuf>,
    quiet: bool,
    embedded: bool,
    region: Option<Region>,
}

#[derive(Debug)]
enum ParseOutcome {
    Help,
    Error,
}

#[derive(Clone, Debug)]
struct RefBuf {
    name: String,
    seq: Vec<u8>,
    touched: bool,
}

#[derive(Clone, Debug)]
struct RegionTarget {
    tid: usize,
    start: usize,
    end: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CigarStep {
    Match { query_pos: usize, ref_pos: usize },
    Del { ref_pos: usize },
}

/// Entry point for `samtools reference`.
pub fn main(args: &[OsString]) -> ExitCode {
    let (opts, input) = match parse_args(args) {
        Ok(v) => v,
        Err(ParseOutcome::Help) => return ExitCode::SUCCESS,
        Err(ParseOutcome::Error) => return ExitCode::from(1),
    };

    if opts.embedded {
        print_error(
            "reference",
            "-e/--embedded requires CRAM container internals not exposed in samtools-rs",
        );
        return ExitCode::from(1);
    }

    let input = input.unwrap_or_else(|| PathBuf::from("-"));
    let mut writer = match sam_io::open_text_output(opts.output.as_deref()) {
        Ok(writer) => writer,
        Err(e) => {
            print_error_errno("reference", "open -o output", &e);
            return ExitCode::from(1);
        }
    };

    if let Err(e) = reference_path(&input, &opts, &mut writer) {
        print_error_errno(
            "reference",
            format!("failed to process \"{}\"", input.display()),
            &e,
        );
        return ExitCode::from(1);
    }

    if let Err(e) = sam_io::check_sam_close(&mut writer) {
        print_error_errno("reference", "close output", &e);
        return ExitCode::from(1);
    }

    ExitCode::SUCCESS
}

fn parse_args(args: &[OsString]) -> Result<(ReferenceOptions, Option<PathBuf>), ParseOutcome> {
    let mut opts = ReferenceOptions::default();
    let mut input = None;
    let mut iter = args.iter().skip(1);

    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-o" | "--output" => {
                let Some(value) = iter.next() else {
                    print_error("reference", "missing value for -o");
                    return Err(ParseOutcome::Error);
                };
                opts.output = Some(PathBuf::from(value));
            }
            "-q" | "--quiet" => opts.quiet = true,
            "-e" | "--embedded" => opts.embedded = true,
            "-r" | "--region" => {
                let Some(raw) = iter.next().and_then(|arg| arg.to_str()) else {
                    print_error("reference", "missing value for -r");
                    return Err(ParseOutcome::Error);
                };
                opts.region = Some(parse_region(raw)?);
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "--help" | "-h" => {
                let _ = print_usage();
                return Err(ParseOutcome::Help);
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("reference", format!("unknown option {s}"));
                return Err(ParseOutcome::Error);
            }
            _ => {
                if input.is_some() {
                    print_error("reference", "multiple input files are not supported");
                    return Err(ParseOutcome::Error);
                }
                input = Some(PathBuf::from(arg));
            }
        }
    }

    Ok((opts, input))
}

fn parse_region(raw: &str) -> Result<Region, ParseOutcome> {
    raw.parse::<Region>().map_err(|e| {
        print_error("reference", format!("invalid region \"{raw}\": {e}"));
        ParseOutcome::Error
    })
}

fn reference_path(input: &Path, opts: &ReferenceOptions, writer: &mut dyn Write) -> io::Result<()> {
    if input.as_os_str() == "-" {
        let stdin = io::stdin().lock();
        let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(stdin));
        let header = reader.read_header()?;
        let mut refs = init_refs(&header);
        let target = opts
            .region
            .as_ref()
            .map(|region| region_target(&header, region))
            .transpose()?;
        for result in reader.records() {
            let record = result?;
            update_refs(&header, &mut refs, &record, target.as_ref())?;
        }
        return dump_refs(writer, &refs, target.as_ref(), opts.quiet);
    }

    match sam_io::sam_open_format(input)?.exact {
        Exact::Sam => {
            let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(File::open(input)?));
            let header = reader.read_header()?;
            let mut refs = init_refs(&header);
            let target = opts
                .region
                .as_ref()
                .map(|region| region_target(&header, region))
                .transpose()?;
            for result in reader.records() {
                let record = result?;
                update_refs(&header, &mut refs, &record, target.as_ref())?;
            }
            dump_refs(writer, &refs, target.as_ref(), opts.quiet)
        }
        Exact::Bam => reference_bam_path(input, opts, writer),
        Exact::Cram => reference_cram_path(input, opts, writer),
        _ => Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "only SAM and BAM input are currently supported",
        )),
    }
}

/// `reference` MD path for CRAM input. The fixture CRAMs are built
/// with `embed_ref=1`, so the embedded reference travels in the
/// container and noodles decodes full SEQ with no external reference;
/// `update_refs` then reconstructs the reference from MD:Z + CIGAR +
/// SEQ exactly as the SAM/BAM paths do. The `-e` embedded-extraction
/// mode still needs CRAM container internals (TODO-NEXT #3).
fn reference_cram_path(
    input: &Path,
    opts: &ReferenceOptions,
    writer: &mut dyn Write,
) -> io::Result<()> {
    if opts.embedded {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "CRAM embedded-reference extraction (-e) needs CRAM container internals",
        ));
    }
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(input)?;
    let mut refs = init_refs(&header);
    let target = opts
        .region
        .as_ref()
        .map(|region| region_target(&header, region))
        .transpose()?;
    // The vendored noodles-cram now decodes an embedded reference
    // (`embed_ref`) directly, so the MD path needs no external
    // reference for embed_ref CRAM; fall back to `--reference` only
    // when one is supplied (reference-compressed CRAM).
    let records = match crate::sam_global::current_global_args().reference {
        Some(reference) => {
            htslib_rs::alignment_compat::query_cram_records_all_from_path_with_reference(
                input, reference,
            )?
        }
        None => htslib_rs::alignment_compat::query_cram_records_all_from_path(input)?,
    };
    for record in records {
        update_refs(&header, &mut refs, &record, target.as_ref())?;
    }
    dump_refs(writer, &refs, target.as_ref(), opts.quiet)
}

fn reference_bam_path(
    input: &Path,
    opts: &ReferenceOptions,
    writer: &mut dyn Write,
) -> io::Result<()> {
    let header = htslib_rs::alignment_compat::read_bam_header_from_path(input)?;
    let mut refs = init_refs(&header);
    let target = opts
        .region
        .as_ref()
        .map(|region| region_target(&header, region))
        .transpose()?;

    if let Some(region) = opts.region.as_ref() {
        match htslib_rs::alignment_compat::query_bam_records_from_path(input, region) {
            Ok(records) => {
                for record in records {
                    update_refs(&header, &mut refs, &record, target.as_ref())?;
                }
                return dump_refs(writer, &refs, target.as_ref(), opts.quiet);
            }
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => return Err(e),
        }
    }

    let mut reader = htslib_rs::bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;
    let mut refs = init_refs(&header);
    let target = opts
        .region
        .as_ref()
        .map(|region| region_target(&header, region))
        .transpose()?;
    let mut record = sam::alignment::RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut record)?;
        if n == 0 {
            break;
        }
        update_refs(&header, &mut refs, &record, target.as_ref())?;
    }
    dump_refs(writer, &refs, target.as_ref(), opts.quiet)
}

fn init_refs(header: &sam::Header) -> Vec<RefBuf> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, reference_sequence)| RefBuf {
            name: String::from_utf8_lossy(name.as_ref()).into_owned(),
            seq: vec![b'N'; usize::from(reference_sequence.length())],
            touched: false,
        })
        .collect()
}

fn region_target(header: &sam::Header, region: &Region) -> io::Result<RegionTarget> {
    let tid = header
        .reference_sequences()
        .get_index_of(region.name())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "region reference sequence does not exist: {}",
                    String::from_utf8_lossy(region.name())
                ),
            )
        })?;
    let (_, def) = header
        .reference_sequences()
        .get_index(tid)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid reference ID"))?;
    let ref_len = usize::from(def.length());
    let interval = region.interval();
    let start = interval.start().map(usize::from).unwrap_or(1);
    let end = interval.end().map(usize::from).unwrap_or(ref_len);
    if start == 0 || end < start {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid region interval: {region}"),
        ));
    }
    Ok(RegionTarget { tid, start, end })
}

fn update_refs<R>(
    header: &sam::Header,
    refs: &mut [RefBuf],
    record: &R,
    target: Option<&RegionTarget>,
) -> io::Result<()>
where
    R: sam::alignment::Record + ?Sized,
{
    let tid = match record.reference_sequence_id(header).transpose()? {
        Some(tid) => tid,
        None => return Ok(()),
    };
    if let Some(target) = target
        && (target.tid != tid || !record_overlaps(record, target)?)
    {
        return Ok(());
    }
    let Some(ref_buf) = refs.get_mut(tid) else {
        return Ok(());
    };
    if build_ref(record, &mut ref_buf.seq)? {
        ref_buf.touched = true;
    }
    Ok(())
}

fn record_overlaps<R>(record: &R, target: &RegionTarget) -> io::Result<bool>
where
    R: sam::alignment::Record + ?Sized,
{
    let Some(start) = record.alignment_start().transpose()? else {
        return Ok(false);
    };
    let Some(end) = record.alignment_end().transpose()? else {
        return Ok(false);
    };
    Ok(usize::from(start) <= target.end && target.start <= usize::from(end))
}

fn build_ref<R>(record: &R, ref_seq: &mut [u8]) -> io::Result<bool>
where
    R: sam::alignment::Record + ?Sized,
{
    let Some(md) = md_tag(record)? else {
        return Ok(false);
    };
    let Some(start) = record.alignment_start().transpose()? else {
        return Ok(false);
    };
    let ref_start = usize::from(start) - 1;
    let sequence = record.sequence().iter().collect::<Vec<_>>();
    let steps = cigar_steps(record, ref_start)?;
    let mut step_index = 0usize;
    let mut md_index = 0usize;
    let md = md.as_bytes();

    while md_index < md.len() {
        if md[md_index].is_ascii_digit() {
            let mut len = 0usize;
            while md_index < md.len() && md[md_index].is_ascii_digit() {
                len = len
                    .checked_mul(10)
                    .and_then(|n| n.checked_add((md[md_index] - b'0') as usize))
                    .ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "MD length overflow")
                    })?;
                md_index += 1;
            }
            for _ in 0..len {
                match next_step(&steps, &mut step_index)? {
                    CigarStep::Match { query_pos, ref_pos } => {
                        if ref_pos < ref_seq.len()
                            && let Some(base) = sequence.get(query_pos).copied()
                        {
                            ref_seq[ref_pos] = base.to_ascii_uppercase();
                        }
                    }
                    CigarStep::Del { .. } => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "MD:Z and CIGAR are incompatible",
                        ));
                    }
                }
            }
        } else if md[md_index] == b'^' {
            md_index += 1;
            while md_index < md.len() && md[md_index].is_ascii_alphabetic() {
                match next_step(&steps, &mut step_index)? {
                    CigarStep::Del { ref_pos } => {
                        if ref_pos < ref_seq.len() {
                            ref_seq[ref_pos] = md[md_index].to_ascii_uppercase();
                        }
                    }
                    CigarStep::Match { .. } => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "MD:Z and CIGAR are incompatible",
                        ));
                    }
                }
                md_index += 1;
            }
        } else if md[md_index].is_ascii_alphabetic() {
            match next_step(&steps, &mut step_index)? {
                CigarStep::Match { ref_pos, .. } => {
                    if ref_pos < ref_seq.len() {
                        ref_seq[ref_pos] = md[md_index].to_ascii_uppercase();
                    }
                }
                CigarStep::Del { .. } => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "MD:Z and CIGAR are incompatible",
                    ));
                }
            }
            md_index += 1;
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid MD:Z character {:?}", md[md_index] as char),
            ));
        }
    }

    Ok(true)
}

fn md_tag<R>(record: &R) -> io::Result<Option<String>>
where
    R: sam::alignment::Record + ?Sized,
{
    use sam::alignment::record::data::field::{Tag, Value};
    let tag = Tag::from([b'M', b'D']);
    let data = record.data();
    let Some(value) = data.get(&tag).transpose()? else {
        return Ok(None);
    };
    match value {
        Value::String(s) => Ok(Some(s.to_string())),
        _ => Ok(None),
    }
}

fn cigar_steps<R>(record: &R, ref_start: usize) -> io::Result<Vec<CigarStep>>
where
    R: sam::alignment::Record + ?Sized,
{
    use sam::alignment::record::cigar::op::Kind;
    let mut query_pos = 0usize;
    let mut ref_pos = ref_start;
    let mut steps = Vec::new();

    for result in record.cigar().iter() {
        let op = result?;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..op.len() {
                    steps.push(CigarStep::Match { query_pos, ref_pos });
                    query_pos += 1;
                    ref_pos += 1;
                }
            }
            Kind::Deletion => {
                for _ in 0..op.len() {
                    steps.push(CigarStep::Del { ref_pos });
                    ref_pos += 1;
                }
            }
            Kind::Insertion | Kind::SoftClip => query_pos += op.len(),
            Kind::Skip => ref_pos += op.len(),
            Kind::HardClip | Kind::Pad => {}
        }
    }

    Ok(steps)
}

fn next_step(steps: &[CigarStep], index: &mut usize) -> io::Result<CigarStep> {
    let Some(step) = steps.get(*index).copied() else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "MD:Z and CIGAR are incompatible",
        ));
    };
    *index += 1;
    Ok(step)
}

fn dump_refs(
    writer: &mut dyn Write,
    refs: &[RefBuf],
    target: Option<&RegionTarget>,
    quiet: bool,
) -> io::Result<()> {
    if let Some(target) = target {
        let Some(ref_buf) = refs.get(target.tid) else {
            return Ok(());
        };
        let start = target.start.saturating_sub(1);
        let end = target.end.min(ref_buf.seq.len());
        let seq = if start < end {
            &ref_buf.seq[start..end]
        } else {
            &[][..]
        };
        writeln!(
            writer,
            ">{}:{}-{} length: {}",
            ref_buf.name,
            target.start,
            target.end,
            seq.len()
        )?;
        write_fasta_sequence(writer, seq)?;
        if !quiet {
            eprintln!(
                "Dump ref {} len {}, coverage {:.2}%",
                target.tid,
                seq.len(),
                percent_non_n(seq)
            );
        }
        return Ok(());
    }

    for (tid, ref_buf) in refs.iter().enumerate() {
        if !ref_buf.touched {
            continue;
        }
        writeln!(writer, ">{}", ref_buf.name)?;
        write_fasta_sequence(writer, &ref_buf.seq)?;
        if !quiet {
            eprintln!(
                "Dump ref {} len {}, coverage {:.2}%",
                tid,
                ref_buf.seq.len(),
                percent_non_n(&ref_buf.seq)
            );
        }
    }
    Ok(())
}

fn write_fasta_sequence(writer: &mut dyn Write, seq: &[u8]) -> io::Result<()> {
    for chunk in seq.chunks(60) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn percent_non_n(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let n = seq
        .iter()
        .filter(|base| base.eq_ignore_ascii_case(&b'N'))
        .count();
    100.0 - n as f64 * 100.0 / seq.len() as f64
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stdout();
    writeln!(
        w,
        "Usage: samtools reference [-@ N] [-r region] [-e] [-q] [-o out.fa] [in.cram]"
    )
}
