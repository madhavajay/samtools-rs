//! `samtools checksum` — order-agnostic sequence-content checksums.
//!
//! This is a partial port of `bam_checksum.c`. It supports SAM/BAM input for
//! the default checksum columns and common filters. CRAM input needs an
//! all-record CRAM iterator, so that piece remains deferred.

use std::collections::BTreeMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use bstr::ByteSlice;
use htslib_rs::format::Exact;
use htslib_rs::sam;

use crate::bam_flag::{
    BAM_FPAIRED, BAM_FQCFAIL, BAM_FREAD1, BAM_FREAD2, BAM_FSECONDARY, BAM_FSUPPLEMENTARY,
    flag_to_str, str_to_flag,
};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sanitize::{SanitizeFlags, parse_sanitize_options, sanitize_record};

const PRIME: u64 = (1u64 << 31) - 1;

#[derive(Clone, Debug)]
struct ChecksumOptions {
    require_flags: u16,
    exclude_flags: u16,
    flag_mask: u16,
    rev_comp: bool,
    nrec: u64,
    verbose: bool,
    show_qc: bool,
    output: Option<PathBuf>,
    tags: String,
    merge: bool,
    tabs: bool,
    in_order: u8,
    check_pos: bool,
    check_cigar: bool,
    check_mate: bool,
    compat: bool,
    sanitize_flags: SanitizeFlags,
}

impl Default for ChecksumOptions {
    fn default() -> Self {
        Self {
            require_flags: 0,
            exclude_flags: (BAM_FSECONDARY | BAM_FSUPPLEMENTARY) as u16,
            flag_mask: (BAM_FPAIRED | BAM_FREAD1 | BAM_FREAD2) as u16,
            rev_comp: true,
            nrec: 0,
            verbose: false,
            show_qc: false,
            output: None,
            tags: "BC,FI,QT,RT,TC".to_string(),
            merge: false,
            tabs: false,
            in_order: 0,
            check_pos: false,
            check_cigar: false,
            check_mate: false,
            compat: false,
            sanitize_flags: SanitizeFlags::empty(),
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Sums {
    seq: [u64; 3],
    name: [u64; 3],
    qual: [u64; 3],
    aux: [u64; 3],
    pos: [u64; 3],
    cigar: [u64; 3],
    mate: [u64; 3],
    count: [u64; 3],
}

impl Default for Sums {
    fn default() -> Self {
        Self {
            seq: [1; 3],
            name: [1; 3],
            qual: [1; 3],
            aux: [1; 3],
            pos: [1; 3],
            cigar: [1; 3],
            mate: [1; 3],
            count: [0; 3],
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Crcs {
    seq: u32,
    name: u32,
    qual: u32,
    aux: u32,
    pos: u32,
    cigar: u32,
    mate: u32,
}

/// Entry point for `samtools checksum`.
pub fn main(args: &[OsString]) -> ExitCode {
    let (opts, inputs) = match parse_args(args) {
        Ok(v) => v,
        Err(ParseOutcome::Help) => return ExitCode::SUCCESS,
        Err(ParseOutcome::Error) => return ExitCode::from(1),
    };

    let mut writer = match sam_io::open_text_output(opts.output.as_deref()) {
        Ok(writer) => writer,
        Err(e) => {
            print_error_errno("checksum", "open -o output", &e);
            return ExitCode::from(1);
        }
    };

    if opts.merge {
        if inputs.is_empty() {
            print_error("checksum", "-m requires at least one checksum input file");
            return ExitCode::from(1);
        }
        if let Err(e) = combine_paths(&inputs, &opts, &mut writer) {
            print_error_errno("checksum", "failed to merge checksum files", &e);
            return ExitCode::from(1);
        }
        if let Err(e) = sam_io::check_sam_close(&mut writer) {
            print_error_errno("checksum", "close output", &e);
            return ExitCode::from(1);
        }
        return ExitCode::SUCCESS;
    }

    let inputs = if inputs.is_empty() {
        vec![PathBuf::from("-")]
    } else {
        inputs
    };

    let mut ret = ExitCode::SUCCESS;
    for input in inputs {
        if let Err(e) = checksum_path(&input, &opts, &mut writer) {
            print_error_errno(
                "checksum",
                format!("error reading from \"{}\"", input.display()),
                &e,
            );
            ret = ExitCode::from(1);
        }
    }

    if let Err(e) = sam_io::check_sam_close(&mut writer) {
        print_error_errno("checksum", "close output", &e);
        return ExitCode::from(1);
    }

    ret
}

enum ParseOutcome {
    Help,
    Error,
}

fn parse_args(args: &[OsString]) -> Result<(ChecksumOptions, Vec<PathBuf>), ParseOutcome> {
    let mut opts = ChecksumOptions::default();
    let mut inputs = Vec::new();
    let mut iter = args.iter().skip(1).peekable();

    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        if let Some(shorts) = s.strip_prefix('-')
            && shorts.len() > 1
            && !shorts.starts_with('-')
            && shorts
                .bytes()
                .all(|b| matches!(b, b'c' | b'q' | b'v' | b'T' | b'O' | b'P' | b'C' | b'M'))
        {
            for b in shorts.bytes() {
                match b {
                    b'c' => opts.rev_comp = false,
                    b'q' => opts.show_qc = true,
                    b'v' => opts.verbose = true,
                    b'T' => opts.tabs = true,
                    b'O' => opts.in_order = opts.in_order.saturating_add(1),
                    b'P' => opts.check_pos = true,
                    b'C' => opts.check_cigar = true,
                    b'M' => opts.check_mate = true,
                    _ => unreachable!(),
                }
            }
            continue;
        }
        match s {
            "-F" | "--exclude-flags" => {
                opts.exclude_flags = parse_flag_value(iter.next(), "-F")?;
            }
            "-f" | "--require-flags" => {
                opts.require_flags = parse_flag_value(iter.next(), "-f")?;
            }
            "-b" | "--flag-mask" => {
                opts.flag_mask = parse_flag_value(iter.next(), "-b")?;
            }
            "-c" | "--no-rev-comp" => opts.rev_comp = false,
            "-N" | "--count" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("checksum", "missing value for -N");
                    return Err(ParseOutcome::Error);
                };
                opts.nrec = match raw.parse() {
                    Ok(n) => n,
                    Err(_) => {
                        print_error("checksum", format!("invalid -N value \"{raw}\""));
                        return Err(ParseOutcome::Error);
                    }
                };
            }
            "-o" | "--output" => opts.output = iter.next().map(PathBuf::from),
            "-q" | "--show-qc" => opts.show_qc = true,
            "-v" | "--verbose" => opts.verbose = true,
            "-t" | "--tags" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("checksum", "missing value for -t");
                    return Err(ParseOutcome::Error);
                };
                opts.tags = raw.to_string();
                if !valid_tags(raw) {
                    print_error("checksum", "Bad tag string. Should be XX,YY,... syntax");
                    return Err(ParseOutcome::Error);
                }
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "-m" | "--merge" => opts.merge = true,
            "-T" | "--tabs" => opts.tabs = true,
            "-O" | "--in-order" => opts.in_order = opts.in_order.saturating_add(1),
            "-P" | "--check-pos" => opts.check_pos = true,
            "-C" | "--check-cigar" => opts.check_cigar = true,
            "-M" | "--check-mate" => opts.check_mate = true,
            "-B" | "--bamseqchksum" => {
                opts.compat = true;
                opts.show_qc = true;
            }
            "-a" | "--all" => apply_all_options(&mut opts),
            "-z" | "--sanitize" => {
                let Some(raw) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("checksum", "missing value for -z");
                    return Err(ParseOutcome::Error);
                };
                opts.sanitize_flags = match parse_sanitize_options(raw) {
                    Ok(flags) => flags,
                    Err(e) => {
                        print_error("checksum", e);
                        return Err(ParseOutcome::Error);
                    }
                };
            }
            "--help" | "-h" => {
                let _ = print_usage();
                return Err(ParseOutcome::Help);
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("checksum", format!("unknown option {s}"));
                return Err(ParseOutcome::Error);
            }
            _ => inputs.push(PathBuf::from(arg)),
        }
    }

    Ok((opts, inputs))
}

fn apply_all_options(opts: &mut ChecksumOptions) {
    opts.require_flags = 0;
    opts.exclude_flags = 0;
    opts.flag_mask = 0x0fff;
    opts.rev_comp = false;
    opts.in_order = 1;
    opts.check_pos = true;
    opts.check_cigar = true;
    opts.check_mate = true;
    opts.sanitize_flags = SanitizeFlags::ALL_WITH_CIGARX;
    opts.tags = "*,cF,MD,NM".to_string();
}

fn parse_flag_value(value: Option<&OsString>, option: &str) -> Result<u16, ParseOutcome> {
    let Some(raw) = value.and_then(|a| a.to_str()) else {
        print_error("checksum", format!("missing value for {option}"));
        return Err(ParseOutcome::Error);
    };
    let Some(flag) = str_to_flag(raw) else {
        print_error("checksum", format!("could not parse flag {raw}"));
        return Err(ParseOutcome::Error);
    };
    Ok(flag as u16)
}

fn valid_tags(raw: &str) -> bool {
    let tags = raw.split(',').collect::<Vec<_>>();
    tags.iter().enumerate().all(|(i, tag)| {
        tag.len() == 2 || (*tag == "*" && i == 0 && (tags.len() == 1 || raw.starts_with("*,")))
    })
}

fn checksum_path(input: &Path, opts: &ChecksumOptions, writer: &mut dyn Write) -> io::Result<()> {
    let mut all = Sums::default();
    let mut no_rg = Sums::default();
    let mut groups = BTreeMap::<String, Sums>::new();

    let stdin_input = input.as_os_str() == "-";
    if stdin_input {
        let stdin = io::stdin().lock();
        let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(stdin));
        let header = reader.read_header()?;
        checksum_sam_reader(
            &mut reader,
            &header,
            opts,
            &mut all,
            &mut no_rg,
            &mut groups,
        )?;
    } else {
        let format = sam_io::sam_open_format(input)?;
        match format.exact {
            Exact::Sam => {
                let file = File::open(input)?;
                let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(file));
                let header = reader.read_header()?;
                checksum_sam_reader(
                    &mut reader,
                    &header,
                    opts,
                    &mut all,
                    &mut no_rg,
                    &mut groups,
                )?;
            }
            Exact::Bam => {
                let mut reader = htslib_rs::bam::io::Reader::new(File::open(input)?);
                let header = reader.read_header()?;
                let mut record = sam::alignment::RecordBuf::default();
                let mut seen = 0u64;
                loop {
                    let n = reader.read_record_buf(&header, &mut record)?;
                    if n == 0 {
                        break;
                    }
                    if update_record_buf(
                        &header,
                        &mut record,
                        opts,
                        &mut all,
                        &mut no_rg,
                        &mut groups,
                    )? {
                        seen += 1;
                        if opts.nrec != 0 && seen == opts.nrec {
                            break;
                        }
                    }
                }
            }
            Exact::Cram => {
                // TODO-NEXT #2: whole-CRAM via the htslib-rs all-record
                // iterator, decoded against the global --reference.
                let Some(reference) = crate::sam_global::current_global_args().reference else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "CRAM checksum requires top-level --reference FILE",
                    ));
                };
                let header = htslib_rs::alignment_compat::read_cram_header_from_path(input)?;
                let mut seen = 0u64;
                for mut record in
                    htslib_rs::alignment_compat::query_cram_records_all_from_path_with_reference(
                        input, &reference,
                    )?
                {
                    if update_record_buf(
                        &header,
                        &mut record,
                        opts,
                        &mut all,
                        &mut no_rg,
                        &mut groups,
                    )? {
                        seen += 1;
                        if opts.nrec != 0 && seen == opts.nrec {
                            break;
                        }
                    }
                }
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "only SAM and BAM input are currently supported",
                ));
            }
        }
    }

    write_report(writer, input, opts, &all, &no_rg, &groups)
}

fn checksum_sam_reader<R>(
    reader: &mut htslib_rs::sam::io::Reader<R>,
    header: &sam::Header,
    opts: &ChecksumOptions,
    all: &mut Sums,
    no_rg: &mut Sums,
    groups: &mut BTreeMap<String, Sums>,
) -> io::Result<()>
where
    R: io::BufRead,
{
    let mut seen = 0u64;
    let mut record = sam::alignment::RecordBuf::default();
    loop {
        let n = reader.read_record_buf(header, &mut record)?;
        if n == 0 {
            break;
        }
        if update_record_buf(header, &mut record, opts, all, no_rg, groups)? {
            seen += 1;
            if opts.nrec != 0 && seen == opts.nrec {
                break;
            }
        }
    }
    Ok(())
}

fn update_record_buf(
    header: &sam::Header,
    record: &mut sam::alignment::RecordBuf,
    opts: &ChecksumOptions,
    all: &mut Sums,
    no_rg: &mut Sums,
    groups: &mut BTreeMap<String, Sums>,
) -> io::Result<bool> {
    let original_flag = record.flags().bits();
    if original_flag & opts.exclude_flags != 0 {
        return Ok(false);
    }
    if (original_flag & opts.require_flags) != opts.require_flags {
        return Ok(false);
    }

    sanitize_record(header, record, opts.sanitize_flags);
    let flag = record.flags().bits();
    let crcs = record_crcs(header, record, opts, flag)?;
    let qcfail = flag & BAM_FQCFAIL as u16 != 0;
    let group = read_group(record)?;

    let group_count = if let Some(group) = group {
        let sums = groups.entry(group).or_default();
        let count = sums.count[0];
        sums_update(qcfail, sums, &crcs, opts, count);
        count
    } else {
        let count = no_rg.count[0];
        sums_update(qcfail, no_rg, &crcs, opts, count);
        count
    };

    sums_update(qcfail, all, &crcs, opts, group_count);
    Ok(true)
}

fn record_crcs<R>(
    header: &sam::Header,
    record: &R,
    opts: &ChecksumOptions,
    flag: u16,
) -> io::Result<Crcs>
where
    R: sam::alignment::Record + ?Sized,
{
    let masked_flag = (flag & opts.flag_mask) as u8;
    let seq = sequence_bytes(record, opts.rev_comp)?;
    let qual = quality_bytes(record, opts.rev_comp)?;
    let name = record
        .name()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing record name"))?;

    let crc0 = crc32(0, &[]);
    let seq_crc = crc32(crc32(crc0, &[masked_flag]), &seq);

    let mut name_with_nul = Vec::with_capacity(name.len() + 1);
    name_with_nul.extend_from_slice(name);
    name_with_nul.push(0);
    let name_crc = crc32(crc32(crc32(crc0, &name_with_nul), &[masked_flag]), &seq);
    let qual_crc = crc32(seq_crc, &qual);

    let aux = aux_bytes(record, &opts.tags)?;
    let aux_crc = crc32(seq_crc, &aux);
    let pos_crc = if opts.check_pos {
        crc32(seq_crc, &position_bytes(header, record)?)
    } else {
        1
    };
    let cigar_crc = if opts.check_cigar {
        let mapq = record
            .mapping_quality()
            .transpose()?
            .map(|mapping_quality| mapping_quality.get())
            .unwrap_or(255);
        let mapq_crc = crc32(seq_crc, &(u32::from(mapq)).to_le_bytes());
        crc32(mapq_crc, &cigar_bytes(record)?)
    } else {
        1
    };
    let mate_crc = if opts.check_mate {
        crc32(seq_crc, &mate_bytes(header, record)?)
    } else {
        1
    };

    Ok(Crcs {
        seq: seq_crc,
        name: name_crc,
        qual: qual_crc,
        aux: aux_crc,
        pos: pos_crc,
        cigar: cigar_crc,
        mate: mate_crc,
    })
}

fn position_bytes<R>(header: &sam::Header, record: &R) -> io::Result<[u8; 12]>
where
    R: sam::alignment::Record + ?Sized,
{
    let tid = record
        .reference_sequence_id(header)
        .transpose()?
        .map(|id| id as i32)
        .unwrap_or(-1);
    let pos = record
        .alignment_start()
        .transpose()?
        .map(|start| usize::from(start) as i64 - 1)
        .unwrap_or(-1);

    let mut bytes = [0; 12];
    bytes[..4].copy_from_slice(&(tid as u32).to_le_bytes());
    bytes[4..].copy_from_slice(&(pos as u64).to_le_bytes());
    Ok(bytes)
}

fn mate_bytes<R>(header: &sam::Header, record: &R) -> io::Result<[u8; 12]>
where
    R: sam::alignment::Record + ?Sized,
{
    let mate_tid = record
        .mate_reference_sequence_id(header)
        .transpose()?
        .map(|id| id as i32)
        .unwrap_or(-1);
    let mate_pos = record
        .mate_alignment_start()
        .transpose()?
        .map(|start| usize::from(start) as i64 - 1)
        .unwrap_or(-1);

    let mut bytes = [0; 12];
    bytes[..4].copy_from_slice(&(mate_tid as u32).to_le_bytes());
    bytes[4..].copy_from_slice(&(mate_pos as u64).to_le_bytes());
    Ok(bytes)
}

fn cigar_bytes<R>(record: &R) -> io::Result<Vec<u8>>
where
    R: sam::alignment::Record + ?Sized,
{
    use sam::alignment::record::cigar::op::Kind;

    let mut bytes = Vec::with_capacity(record.cigar().len() * 4);
    for result in record.cigar().iter() {
        let op = result?;
        let op_code = match op.kind() {
            Kind::Match => 0,
            Kind::Insertion => 1,
            Kind::Deletion => 2,
            Kind::Skip => 3,
            Kind::SoftClip => 4,
            Kind::HardClip => 5,
            Kind::Pad => 6,
            Kind::SequenceMatch => 7,
            Kind::SequenceMismatch => 8,
        };
        let len =
            u32::try_from(op.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let word = (len << 4) | op_code;
        bytes.extend_from_slice(&word.to_le_bytes());
    }
    Ok(bytes)
}

fn sequence_bytes<R>(record: &R, rev_comp: bool) -> io::Result<Vec<u8>>
where
    R: sam::alignment::Record + ?Sized,
{
    let mut bases = record.sequence().iter().collect::<Vec<_>>();
    if rev_comp && record.flags()?.is_reverse_complemented() {
        bases = bases.into_iter().rev().map(complement_base).collect();
    }
    Ok(bases)
}

fn quality_bytes<R>(record: &R, rev_comp: bool) -> io::Result<Vec<u8>>
where
    R: sam::alignment::Record + ?Sized,
{
    let mut scores = record
        .quality_scores()
        .iter()
        .map(|result| {
            result.and_then(|score| {
                score.checked_add(b'!').ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "quality score overflow")
                })
            })
        })
        .collect::<io::Result<Vec<_>>>()?;
    if rev_comp && record.flags()?.is_reverse_complemented() {
        scores.reverse();
    }
    Ok(scores)
}

fn aux_bytes<R>(record: &R, tags: &str) -> io::Result<Vec<u8>>
where
    R: sam::alignment::Record + ?Sized,
{
    let tag_list = tags
        .split(',')
        .filter_map(|tag| tag.as_bytes().try_into().ok())
        .collect::<Vec<[u8; 2]>>();
    let include_all = tags == "*" || tags.starts_with("*,");
    let excluded = if include_all {
        tag_list.as_slice()
    } else {
        &[][..]
    };
    let mut fields = BTreeMap::<[u8; 2], Vec<u8>>::new();

    for result in record.data().iter() {
        let (tag, value) = result?;
        let tag_bytes = <[u8; 2]>::from(tag);
        let keep = if include_all {
            (b'0'..=b'z').contains(&tag_bytes[0])
                && (b'0'..=b'z').contains(&tag_bytes[1])
                && !excluded.iter().any(|excluded| excluded == &tag_bytes)
        } else {
            tag_list.iter().any(|wanted| wanted == &tag_bytes)
        };
        if keep && let Some(encoded) = encode_aux_value(&tag_bytes, value)? {
            fields.insert(tag_bytes, encoded);
        }
    }

    let mut out = Vec::new();
    if include_all {
        for encoded in fields.values() {
            out.extend_from_slice(encoded);
        }
    } else {
        for tag in tag_list {
            if let Some(encoded) = fields.get(&tag) {
                out.extend_from_slice(encoded);
            }
        }
    }
    Ok(out)
}

fn encode_aux_value(
    tag: &[u8; 2],
    value: sam::alignment::record::data::field::Value<'_>,
) -> io::Result<Option<Vec<u8>>> {
    use sam::alignment::record::data::field::Value;

    let mut out = Vec::new();
    out.extend_from_slice(tag);
    match value {
        Value::Character(c) => {
            out.push(b'A');
            out.push(c);
        }
        Value::Int8(n) => encode_aux_int(&mut out, i64::from(n)),
        Value::UInt8(n) => encode_aux_int(&mut out, i64::from(n)),
        Value::Int16(n) => encode_aux_int(&mut out, i64::from(n)),
        Value::UInt16(n) => encode_aux_int(&mut out, i64::from(n)),
        Value::Int32(n) => encode_aux_int(&mut out, i64::from(n)),
        Value::UInt32(n) => encode_aux_int(&mut out, i64::from(n)),
        Value::Float(n) => {
            out.push(b'f');
            out.extend_from_slice(&n.to_le_bytes());
        }
        Value::String(s) => {
            out.push(b'Z');
            out.extend_from_slice(s.as_bytes());
            out.push(0);
        }
        Value::Hex(s) => {
            out.push(b'H');
            out.extend_from_slice(s.as_bytes());
            out.push(0);
        }
        Value::Array(array) => encode_aux_array(&mut out, array)?,
    }
    Ok(Some(out))
}

fn encode_aux_int(out: &mut Vec<u8>, n: i64) {
    if n >= 0 {
        if n <= u8::MAX.into() {
            out.push(b'C');
            out.push(n as u8);
        } else if n <= u16::MAX.into() {
            out.push(b'S');
            out.extend_from_slice(&(n as u16).to_le_bytes());
        } else {
            out.push(b'I');
            out.extend_from_slice(&(n as u32).to_le_bytes());
        }
    } else if n >= i8::MIN.into() && n <= i8::MAX.into() {
        out.push(b'c');
        out.push(n as i8 as u8);
    } else if n >= i16::MIN.into() && n <= i16::MAX.into() {
        out.push(b's');
        out.extend_from_slice(&(n as i16).to_le_bytes());
    } else {
        out.push(b'i');
        out.extend_from_slice(&(n as i32).to_le_bytes());
    }
}

fn encode_aux_array(
    out: &mut Vec<u8>,
    array: sam::alignment::record::data::field::value::Array<'_>,
) -> io::Result<()> {
    use sam::alignment::record::data::field::value::Array;

    out.push(b'B');
    match array {
        Array::Int8(values) => {
            out.push(b'c');
            out.extend_from_slice(&(values.len() as u32).to_le_bytes());
            for value in values.iter() {
                out.push(value? as u8);
            }
        }
        Array::UInt8(values) => {
            out.push(b'C');
            out.extend_from_slice(&(values.len() as u32).to_le_bytes());
            for value in values.iter() {
                out.push(value?);
            }
        }
        Array::Int16(values) => {
            out.push(b's');
            out.extend_from_slice(&(values.len() as u32).to_le_bytes());
            for value in values.iter() {
                out.extend_from_slice(&value?.to_le_bytes());
            }
        }
        Array::UInt16(values) => {
            out.push(b'S');
            out.extend_from_slice(&(values.len() as u32).to_le_bytes());
            for value in values.iter() {
                out.extend_from_slice(&value?.to_le_bytes());
            }
        }
        Array::Int32(values) => {
            out.push(b'i');
            out.extend_from_slice(&(values.len() as u32).to_le_bytes());
            for value in values.iter() {
                out.extend_from_slice(&value?.to_le_bytes());
            }
        }
        Array::UInt32(values) => {
            out.push(b'I');
            out.extend_from_slice(&(values.len() as u32).to_le_bytes());
            for value in values.iter() {
                out.extend_from_slice(&value?.to_le_bytes());
            }
        }
        Array::Float(values) => {
            out.push(b'f');
            out.extend_from_slice(&(values.len() as u32).to_le_bytes());
            for value in values.iter() {
                out.extend_from_slice(&value?.to_le_bytes());
            }
        }
    }
    Ok(())
}

fn read_group<R>(record: &R) -> io::Result<Option<String>>
where
    R: sam::alignment::Record + ?Sized,
{
    use sam::alignment::record::data::field::{Tag, Value};
    let tag = Tag::from([b'R', b'G']);
    let data = record.data();
    let Some(value) = data.get(&tag).transpose()? else {
        return Ok(None);
    };
    match value {
        Value::String(s) => Ok(Some(s.to_string())),
        _ => Ok(None),
    }
}

fn sums_update(qcfail: bool, sums: &mut Sums, crcs: &Crcs, opts: &ChecksumOptions, count: u64) {
    let count_crc = if opts.in_order == 0 {
        0
    } else {
        let order_count = if opts.in_order == 1 {
            count
        } else {
            sums.count[0]
        };
        crc32(0, &order_count.to_le_bytes())
    };

    sums_update_row(0, sums, crcs, count_crc);
    if opts.show_qc && !qcfail {
        sums_update_row(1, sums, crcs, count_crc);
    }
    if opts.show_qc && qcfail {
        sums_update_row(2, sums, crcs, count_crc);
    }
}

fn sums_update_row(row: usize, sums: &mut Sums, crcs: &Crcs, count_crc: u32) {
    sums.seq[row] = update_hash(sums.seq[row], count_crc ^ crcs.seq);
    sums.name[row] = update_hash(sums.name[row], count_crc ^ crcs.name);
    sums.qual[row] = update_hash(sums.qual[row], count_crc ^ crcs.qual);
    sums.aux[row] = update_hash(sums.aux[row], count_crc ^ crcs.aux);
    sums.pos[row] = update_hash(sums.pos[row], count_crc ^ crcs.pos);
    sums.cigar[row] = update_hash(sums.cigar[row], count_crc ^ crcs.cigar);
    sums.mate[row] = update_hash(sums.mate[row], count_crc ^ crcs.mate);
    sums.count[row] += 1;
}

fn sums_update_row_n(row: usize, sums: &mut Sums, crcs: &Crcs, n: u64) {
    sums.seq[row] = update_hash(sums.seq[row], crcs.seq);
    sums.name[row] = update_hash(sums.name[row], crcs.name);
    sums.qual[row] = update_hash(sums.qual[row], crcs.qual);
    sums.aux[row] = update_hash(sums.aux[row], crcs.aux);
    sums.pos[row] = update_hash(sums.pos[row], crcs.pos);
    sums.cigar[row] = update_hash(sums.cigar[row], crcs.cigar);
    sums.mate[row] = update_hash(sums.mate[row], crcs.mate);
    sums.count[row] += n;
}

fn combine_paths(
    inputs: &[PathBuf],
    opts: &ChecksumOptions,
    writer: &mut dyn Write,
) -> io::Result<()> {
    let mut all = Sums::default();
    let mut no_rg = Sums::default();
    let mut groups = BTreeMap::<String, Sums>::new();
    let mut merged_opts = opts.clone();

    for input in inputs {
        parse_checksum_file(input, &mut merged_opts, &mut all, &mut no_rg, &mut groups)?;
    }

    write_report(
        writer,
        Path::new("merge"),
        &merged_opts,
        &all,
        &no_rg,
        &groups,
    )
}

fn parse_checksum_file(
    input: &Path,
    opts: &mut ChecksumOptions,
    all: &mut Sums,
    no_rg: &mut Sums,
    groups: &mut BTreeMap<String, Sums>,
) -> io::Result<()> {
    use std::io::BufRead;

    let file = File::open(input)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        if let Some(tags) = line.strip_prefix("# Aux tags:") {
            opts.tags = tags.trim().to_string();
            continue;
        }
        if let Some(flags) = line.strip_prefix("# BAM flags:") {
            opts.flag_mask = str_to_flag(flags.trim()).unwrap_or(0) as u16;
            continue;
        }
        if line.starts_with("# Group") {
            opts.check_pos = line.split_whitespace().any(|token| token == "+chr/pos");
            opts.check_cigar = line.split_whitespace().any(|token| token == "+cigar");
            opts.check_mate = line.split_whitespace().any(|token| token == "+mate");
            continue;
        }
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let cols = line.split_whitespace().collect::<Vec<_>>();
        let expected_cols = 8
            + usize::from(opts.check_pos)
            + usize::from(opts.check_cigar)
            + usize::from(opts.check_mate);
        if cols.len() != expected_cols {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("incorrect number of columns in {}", input.display()),
            ));
        }
        if cols[0] == "all" {
            continue;
        }

        let row = match cols[1] {
            "all" => 0,
            "pass" => {
                opts.show_qc = true;
                1
            }
            "fail" => {
                opts.show_qc = true;
                2
            }
            other => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid QC column {other:?}"),
                ));
            }
        };
        let count = cols[2]
            .parse::<u64>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let crcs = Crcs {
            seq: parse_hex(cols[3])?,
            name: parse_hex(cols[4])?,
            qual: parse_hex(cols[5])?,
            aux: parse_hex(cols[6])?,
            pos: if opts.check_pos {
                parse_hex(cols[7])?
            } else {
                1
            },
            cigar: if opts.check_cigar {
                parse_hex(cols[7 + usize::from(opts.check_pos)])?
            } else {
                1
            },
            mate: if opts.check_mate {
                parse_hex(cols[7 + usize::from(opts.check_pos) + usize::from(opts.check_cigar)])?
            } else {
                1
            },
        };

        if cols[0] == "-" {
            sums_update_row_n(row, no_rg, &crcs, count);
        } else {
            sums_update_row_n(
                row,
                groups.entry(cols[0].to_string()).or_default(),
                &crcs,
                count,
            );
        }
        sums_update_row_n(row, all, &crcs, count);
    }

    Ok(())
}

fn parse_hex(s: &str) -> io::Result<u32> {
    u32::from_str_radix(s, 16).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn update_hash(hash: u64, mut crc: u32) -> u64 {
    crc &= PRIME as u32;
    if crc == 0 || crc == PRIME as u32 {
        crc = 1;
    }
    (hash * u64::from(crc)) % PRIME
}

fn write_report(
    writer: &mut dyn Write,
    input: &Path,
    opts: &ChecksumOptions,
    all: &Sums,
    no_rg: &Sums,
    groups: &BTreeMap<String, Sums>,
) -> io::Result<()> {
    if opts.compat {
        return write_bamseqchksum_report(writer, all, no_rg, groups);
    }

    let sep = if opts.tabs { "\t" } else { " " };
    let aux_sep = if opts.tabs { "\t" } else { "          " };
    let flag_sep = if opts.tabs { "\t" } else { "         " };
    writeln!(writer, "# Checksum 1.0 for file:{sep}{}", input.display())?;
    writeln!(writer, "# Aux tags:{aux_sep}{}", opts.tags)?;
    writeln!(
        writer,
        "# BAM flags:{flag_sep}{}",
        flag_to_str(opts.flag_mask as u32)
    )?;
    writeln!(writer)?;
    if opts.tabs {
        write!(writer, "# Group\tQC\tcount\tflag+seq\t+name\t+qual\t+aux")?;
        if opts.check_pos {
            write!(writer, "\t+chr/pos")?;
        }
        if opts.check_cigar {
            write!(writer, "\t+cigar")?;
        }
        if opts.check_mate {
            write!(writer, "\t+mate")?;
        }
        writeln!(writer, "\tcombined")?;
    } else {
        write!(
            writer,
            "# Group    QC          count  flag+seq  +name     +qual     +aux    "
        )?;
        if opts.check_pos {
            write!(writer, "  +chr/pos")?;
        }
        if opts.check_cigar {
            write!(writer, "  +cigar  ")?;
        }
        if opts.check_mate {
            write!(writer, "  +mate   ")?;
        }
        writeln!(writer, "  combined")?;
    }

    write_sums(writer, "all", all, opts)?;
    if opts.verbose || no_rg.count[0] + no_rg.count[1] != 0 {
        write_sums(writer, "-", no_rg, opts)?;
    }
    for (group, sums) in groups {
        write_sums(writer, group, sums, opts)?;
    }
    Ok(())
}

fn write_bamseqchksum_report(
    writer: &mut dyn Write,
    all: &Sums,
    no_rg: &Sums,
    groups: &BTreeMap<String, Sums>,
) -> io::Result<()> {
    writeln!(
        writer,
        "###\tset\tcount\t\tb_seq\tname_b_seq\tb_seq_qual\tb_seq_tags(BC,FI,QT,RT,TC)"
    )?;
    write_bamseqchksum_sums(writer, "all", all)?;
    write_bamseqchksum_sums(writer, "", no_rg)?;
    for (group, sums) in groups {
        write_bamseqchksum_sums(writer, group, sums)?;
    }
    Ok(())
}

fn write_bamseqchksum_sums(writer: &mut dyn Write, group: &str, sums: &Sums) -> io::Result<()> {
    for (row, qc) in ["all", "pass"].iter().enumerate() {
        writeln!(
            writer,
            "{}\t{}\t{}\t{:x}\t{:x}\t{:x}\t{:x}",
            group,
            qc,
            sums.count[row],
            sums.seq[row],
            sums.name[row],
            sums.qual[row],
            sums.aux[row]
        )?;
    }
    Ok(())
}

fn write_sums(
    writer: &mut dyn Write,
    group: &str,
    sums: &Sums,
    opts: &ChecksumOptions,
) -> io::Result<()> {
    for (row, qc) in ["all", "pass", "fail"].iter().enumerate() {
        if row > 0 && !opts.show_qc {
            continue;
        }
        if !opts.verbose && sums.count[row] == 0 {
            continue;
        }
        let combined = combined_checksum(sums, row, opts);
        if opts.tabs {
            write!(
                writer,
                "{}\t{}\t{}\t{:x}\t{:x}\t{:x}\t{:x}",
                group,
                qc,
                sums.count[row],
                sums.seq[row],
                sums.name[row],
                sums.qual[row],
                sums.aux[row]
            )?;
            if opts.check_pos {
                write!(writer, "\t{:x}", sums.pos[row])?;
            }
            if opts.check_cigar {
                write!(writer, "\t{:x}", sums.cigar[row])?;
            }
            if opts.check_mate {
                write!(writer, "\t{:x}", sums.mate[row])?;
            }
            writeln!(writer, "\t{:x}", combined)?;
        } else {
            write!(
                writer,
                "{:<10} {:<4} {:12}  {:08x}  {:08x}  {:08x}  {:08x}",
                group,
                qc,
                sums.count[row],
                sums.seq[row],
                sums.name[row],
                sums.qual[row],
                sums.aux[row]
            )?;
            if opts.check_pos {
                write!(writer, "  {:08x}", sums.pos[row])?;
            }
            if opts.check_cigar {
                write!(writer, "  {:08x}", sums.cigar[row])?;
            }
            if opts.check_mate {
                write!(writer, "  {:08x}", sums.mate[row])?;
            }
            writeln!(writer, "  {:08x}", combined)?;
        }
    }
    Ok(())
}

fn combined_checksum(sums: &Sums, row: usize, opts: &ChecksumOptions) -> u64 {
    let mut h = 1;
    h = update_hash(h, (sums.count[row] >> 32) as u32);
    h = update_hash(h, sums.count[row] as u32);
    h = update_hash(h, sums.seq[row] as u32);
    h = update_hash(h, sums.name[row] as u32);
    h = update_hash(h, sums.seq[row] as u32);
    h = update_hash(h, sums.aux[row] as u32);
    if opts.check_pos {
        h = update_hash(h, sums.pos[row] as u32);
    }
    if opts.check_cigar {
        h = update_hash(h, sums.cigar[row] as u32);
    }
    if opts.check_mate {
        h = update_hash(h, sums.mate[row] as u32);
    }
    h
}

fn crc32(initial: u32, bytes: &[u8]) -> u32 {
    let mut crc = libdeflater::Crc::with_initial(initial);
    crc.update(bytes);
    crc.sum()
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

fn print_usage() -> io::Result<()> {
    let mut w = io::stdout();
    writeln!(w, "Usage: samtools checksum [options] [file.bam ...]")?;
    writeln!(w)?;
    writeln!(w, "Options:")?;
    writeln!(
        w,
        "  -F, --exclude-flags FLAG    Filter if any FLAGs are present [0x900]"
    )?;
    writeln!(
        w,
        "  -f, --require-flags FLAG    Filter unless all FLAGs are present [0]"
    )?;
    writeln!(
        w,
        "  -b, --flag-mask FLAG        BAM FLAGs to use in checksums [0x0c1]"
    )?;
    writeln!(
        w,
        "  -c, --no-rev-comp           Do not reverse-complement sequences"
    )?;
    writeln!(
        w,
        "  -t, --tags STR[,STR]        Select tags to checksum [BC,FI,QT,RT,TC]"
    )?;
    writeln!(
        w,
        "  -N, --count INT             Stop after INT records [0]"
    )?;
    writeln!(
        w,
        "  -o, --output FILE           Write report to FILE [stdout]"
    )?;
    writeln!(
        w,
        "  -q, --show-qc               Also show QC pass/fail lines"
    )?;
    writeln!(
        w,
        "  -T, --tabs                  Format output as tab delimited text"
    )?;
    writeln!(
        w,
        "  -P, --check-pos             Also checksum CHR / POS [off]"
    )?;
    writeln!(
        w,
        "  -C, --check-cigar           Also checksum MAPQ / CIGAR [off]"
    )?;
    writeln!(
        w,
        "  -M, --check-mate            Also checksum PNEXT / RNEXT / TLEN [off]"
    )?;
    writeln!(
        w,
        "  -z, --sanitize FLAGS        Perform sanity checks and fix records [off]"
    )?;
    writeln!(
        w,
        "  -B, --bamseqchksum          Report in bamseqchksum format"
    )?;
    writeln!(w, "  -v, --verbose               Show zero-count rows")?;
    Ok(())
}
