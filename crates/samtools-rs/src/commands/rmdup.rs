//! `samtools rmdup` — remove PCR duplicates (deprecated; `markdup` is preferred).
//!
//! Mirrors `bam_rmdup` / `bam_rmdupse` in upstream samtools. Upstream's
//! implementation is paired-aware and works on coordinate-sorted BAMs.
//!
//! This Rust port implements single-end and adjacent paired-end duplicate
//! removal for SAM/BAM/reference-backed CRAM input. SE records are keyed by
//! `(reference_sequence_id, alignment_start, reverse-flag)`, while PE records
//! are paired by qname and keyed by the canonical pair of end coordinates. The
//! record or pair with the highest mapping quality score is retained.
//! `-O sam|bam|cram` / `--output-fmt[=]FMT` selects output format, and CRAM
//! input/output requires `-T` / `--reference` or top-level `--reference`.

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Cursor, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam::{self, alignment::RecordBuf, alignment::io::Write as _};

use crate::bam_flag::{BAM_FMUNMAP, BAM_FPAIRED, BAM_FREVERSE, BAM_FUNMAP};
use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools rmdup`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut output: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let mut no_pg = false;
    let mut single_end = false;
    let mut output_fmt: Option<OutFmt> = None;
    let mut reference: Option<PathBuf> = None;
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-S" | "-s" => {
                single_end = true;
            }
            "-O" | "--output-fmt" => {
                let value = iter.next().and_then(|a| a.to_str()).unwrap_or("bam");
                output_fmt = match parse_output_format(value) {
                    Ok(fmt) => Some(fmt),
                    Err(e) => {
                        print_error("rmdup", e);
                        return ExitCode::from(1);
                    }
                };
            }
            _ if s.starts_with("--output-fmt=") => {
                let value = &s["--output-fmt=".len()..];
                output_fmt = match parse_output_format(value) {
                    Ok(fmt) => Some(fmt),
                    Err(e) => {
                        print_error("rmdup", e);
                        return ExitCode::from(1);
                    }
                };
            }
            "-T" | "--reference" => {
                reference = iter.next().map(PathBuf::from);
            }
            _ if s.starts_with("--reference=") => {
                reference = Some(PathBuf::from(&s["--reference=".len()..]));
            }
            "--no-PG" => {
                no_pg = true;
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("rmdup", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                } else if output.is_none() && s != "-" {
                    // A `-` output operand means stdout (output stays None).
                    output = Some(PathBuf::from(arg));
                }
            }
        }
    }

    let Some(input) = input else {
        let _ = print_usage();
        return ExitCode::from(1);
    };

    if input.as_os_str() != "-" && !input.exists() {
        print_hts_open_missing(&input);
        print_error(
            "rmdup",
            format!(
                "failed to open \"{}\" for input: No such file or directory",
                input.display()
            ),
        );
        return ExitCode::from(1);
    }

    let format = match sam_io::sam_open_format(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error("rmdup", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam | Exact::Cram) {
        print_error(
            "rmdup",
            "only SAM, BAM, and reference-backed CRAM input are currently supported",
        );
        return ExitCode::from(1);
    }

    let output_fmt =
        output_fmt.unwrap_or_else(|| infer_output_format(format.exact, output.as_deref()));
    let pg_argv = if no_pg { None } else { Some(args) };
    let result = match format.exact {
        Exact::Sam => run_sam_rmdup(
            &input,
            output.as_deref(),
            output_fmt,
            reference.as_deref(),
            pg_argv,
            single_end,
        ),
        Exact::Bam => run_bam_rmdup(
            &input,
            output.as_deref(),
            output_fmt,
            reference.as_deref(),
            pg_argv,
            single_end,
        ),
        Exact::Cram => run_cram_rmdup(
            &input,
            output.as_deref(),
            output_fmt,
            reference.as_deref(),
            pg_argv,
            single_end,
        ),
        _ => unreachable!("format checked above"),
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("rmdup", "rmdup failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum OutFmt {
    Sam,
    Bam,
    Cram,
}

fn parse_output_format(raw: &str) -> Result<OutFmt, String> {
    let head = raw.split(',').next().unwrap_or("").to_ascii_lowercase();
    match head.as_str() {
        "sam" => Ok(OutFmt::Sam),
        "bam" => Ok(OutFmt::Bam),
        "cram" => Ok(OutFmt::Cram),
        _ => Err(format!("unsupported output format \"{}\"", raw)),
    }
}

fn infer_output_format(input_fmt: Exact, output: Option<&Path>) -> OutFmt {
    if output
        .and_then(|p| p.extension())
        .and_then(|s| s.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("cram"))
    {
        return OutFmt::Cram;
    }

    match input_fmt {
        Exact::Sam => OutFmt::Sam,
        _ => OutFmt::Bam,
    }
}

fn rmdup_reference(local: Option<&Path>) -> Option<PathBuf> {
    local
        .map(Path::to_path_buf)
        .or_else(|| current_global_args().reference)
}

fn run_bam_rmdup(
    input: &Path,
    output: Option<&Path>,
    output_fmt: OutFmt,
    reference: Option<&Path>,
    pg_argv: Option<&[OsString]>,
    single_end: bool,
) -> io::Result<()> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let mut header = reader.read_header()?;
    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

    let mut records: Vec<RecordBuf> = Vec::new();
    let mut record = RecordBuf::default();
    loop {
        let n = reader.read_record_buf(&header, &mut record)?;
        if n == 0 {
            break;
        }
        records.push(record.clone());
    }

    let result = duplicate_keep_mask_for_records(&records, single_end);
    emit_rmdup_diagnostics(&header, &result.stats, single_end);
    write_kept_records(
        &header,
        &records,
        &result.keep,
        output,
        output_fmt,
        reference,
    )
}

fn run_cram_rmdup(
    input: &Path,
    output: Option<&Path>,
    output_fmt: OutFmt,
    reference: Option<&Path>,
    pg_argv: Option<&[OsString]>,
    single_end: bool,
) -> io::Result<()> {
    let reference = rmdup_reference(reference).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "CRAM input requires -T/--reference FILE or top-level --reference FILE",
        )
    })?;
    crate::reference::ensure_fai_index(&reference, None)?;
    let sam =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            input, &reference, None,
        )?;
    let mut reader = sam::io::Reader::new(BufReader::new(Cursor::new(sam.into_bytes())));
    let mut header = reader.read_header()?;
    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

    let mut records: Vec<RecordBuf> = Vec::new();
    loop {
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        records.push(record);
    }

    let result = duplicate_keep_mask_for_records(&records, single_end);
    emit_rmdup_diagnostics(&header, &result.stats, single_end);
    write_kept_records(
        &header,
        &records,
        &result.keep,
        output,
        output_fmt,
        Some(&reference),
    )
}

fn run_sam_rmdup(
    input: &Path,
    output: Option<&Path>,
    output_fmt: OutFmt,
    reference: Option<&Path>,
    pg_argv: Option<&[OsString]>,
    single_end: bool,
) -> io::Result<()> {
    let mut reader = crate::sam_compat::open_sam_reader_tolerant(input)?;
    let mut header = reader.read_header()?;
    if let Some(argv) = pg_argv {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

    let mut records: Vec<RecordBuf> = Vec::new();
    loop {
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        records.push(record);
    }

    let result = duplicate_keep_mask_for_records(&records, single_end);
    emit_rmdup_diagnostics(&header, &result.stats, single_end);
    write_kept_records(
        &header,
        &records,
        &result.keep,
        output,
        output_fmt,
        reference,
    )
}

struct RmdupResult {
    keep: Vec<bool>,
    stats: RmdupStats,
}

#[derive(Default)]
struct RmdupStats {
    total_units: u64,
    removed_units: u64,
    reference_ids: Vec<usize>,
}

fn duplicate_keep_mask_for_records(records: &[RecordBuf], single_end: bool) -> RmdupResult {
    type PosKey = (i32, i64, bool);
    type PairKey = (PosKey, PosKey);
    type PairIdx = (usize, usize);
    let mut se_best: HashMap<PosKey, usize> = HashMap::new();
    let mut pair_pending: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut pair_best: HashMap<PairKey, PairIdx> = HashMap::new();
    let mut keep = vec![false; records.len()];
    let mut stats = RmdupStats::default();
    for (i, rec) in records.iter().enumerate() {
        let flag = rec.flags().bits() as u32;
        if flag & BAM_FUNMAP != 0 {
            keep[i] = true;
            continue;
        }
        let tid = rec.reference_sequence_id().map(|t| t as i32).unwrap_or(-1);
        if tid >= 0 && !stats.reference_ids.contains(&(tid as usize)) {
            stats.reference_ids.push(tid as usize);
        }
        let pos = rec.alignment_start().map(usize::from).unwrap_or(0) as i64;
        let rev = flag & BAM_FREVERSE != 0;
        let mapq = rec.mapping_quality().map(u8::from).unwrap_or(0);
        let me = (tid, pos, rev);

        let paired_both_mapped = !single_end && flag & BAM_FPAIRED != 0 && flag & BAM_FMUNMAP == 0;
        if paired_both_mapped {
            let name = rec.name().map(|n| n.to_vec()).unwrap_or_default();
            match pair_pending.remove(&name) {
                None => {
                    pair_pending.insert(name, i);
                }
                Some(first_idx) => {
                    stats.total_units += 1;
                    let first = pos_key(&records[first_idx]);
                    let key = if first <= me {
                        (first, me)
                    } else {
                        (me, first)
                    };
                    let score = pair_score(records, first_idx, i);
                    match pair_best.get(&key).copied() {
                        Some((prev_first, prev_second)) => {
                            let prev_score = pair_score(records, prev_first, prev_second);
                            if score > prev_score {
                                keep[prev_first] = false;
                                keep[prev_second] = false;
                                keep[first_idx] = true;
                                keep[i] = true;
                                pair_best.insert(key, (first_idx, i));
                            }
                            stats.removed_units += 1;
                        }
                        None => {
                            keep[first_idx] = true;
                            keep[i] = true;
                            pair_best.insert(key, (first_idx, i));
                        }
                    }
                }
            }
            continue;
        }

        stats.total_units += 1;
        match se_best.get(&me) {
            Some(&idx) => {
                let prev_mapq = records[idx].mapping_quality().map(u8::from).unwrap_or(0);
                if mapq > prev_mapq {
                    keep[idx] = false;
                    keep[i] = true;
                    se_best.insert(me, i);
                }
                stats.removed_units += 1;
            }
            None => {
                keep[i] = true;
                se_best.insert(me, i);
            }
        }
    }

    for idx in pair_pending.into_values() {
        keep[idx] = true;
        if !single_end {
            stats.total_units += 1;
        }
    }
    RmdupResult { keep, stats }
}

fn emit_rmdup_diagnostics(header: &sam::Header, stats: &RmdupStats, single_end: bool) {
    if !single_end {
        for &tid in &stats.reference_ids {
            if let Some((name, _)) = header.reference_sequences().get_index(tid) {
                eprintln!(
                    "[bam_rmdup_core] processing reference {}...",
                    String::from_utf8_lossy(name)
                );
            }
        }
    }

    let fraction = if stats.total_units > 0 {
        stats.removed_units as f64 / stats.total_units as f64
    } else {
        0.0
    };
    let prefix = if single_end {
        "bam_rmdupse_core"
    } else {
        "bam_rmdup_core"
    };
    eprintln!(
        "[{prefix}] {} / {} = {:.4} in library '\t'",
        stats.removed_units, stats.total_units, fraction
    );
}

trait BamLike {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct BamFile(bam::io::Writer<bgzf::io::Writer<File>>);
struct BamStdout(bam::io::Writer<bgzf::io::Writer<io::Stdout>>);

impl BamLike for BamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        self.0.write_alignment_record(header, record)
    }
}
impl BamLike for BamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        self.0.write_alignment_record(header, record)
    }
}

fn open_bam_output(out: Option<&Path>, header: &sam::Header) -> io::Result<Box<dyn BamLike>> {
    match out {
        Some(p) => {
            let mut writer = bam::io::Writer::new(File::create(p)?);
            writer.write_header(header)?;
            Ok(Box::new(BamFile(writer)))
        }
        None => {
            let mut writer = bam::io::Writer::new(io::stdout());
            writer.write_header(header)?;
            Ok(Box::new(BamStdout(writer)))
        }
    }
}

trait SamLike {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()>;
}

struct SamFile(File);
struct SamStdout(io::Stdout);

impl SamLike for SamFile {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        // Shared renderer: htslib `%g` float aux spelling.
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}
impl SamLike for SamStdout {
    fn write_record(&mut self, header: &sam::Header, record: &RecordBuf) -> io::Result<()> {
        crate::sam_render::write_record(&mut self.0, header, record)
    }
}

fn open_sam_output(out: Option<&Path>, header: &sam::Header) -> io::Result<Box<dyn SamLike>> {
    match out {
        Some(p) => {
            let mut writer = File::create(p)?;
            crate::sam_render::write_header(&mut writer, header)?;
            Ok(Box::new(SamFile(writer)))
        }
        None => {
            let mut writer = io::stdout();
            crate::sam_render::write_header(&mut writer, header)?;
            Ok(Box::new(SamStdout(writer)))
        }
    }
}

fn write_kept_records(
    header: &sam::Header,
    records: &[RecordBuf],
    keep: &[bool],
    output: Option<&Path>,
    fmt: OutFmt,
    reference: Option<&Path>,
) -> io::Result<()> {
    match fmt {
        OutFmt::Sam => {
            let mut writer = open_sam_output(output, header)?;
            for (i, rec) in records.iter().enumerate() {
                if keep[i] {
                    writer.write_record(header, rec)?;
                }
            }
            Ok(())
        }
        OutFmt::Bam => {
            let mut writer = open_bam_output(output, header)?;
            for (i, rec) in records.iter().enumerate() {
                if keep[i] {
                    writer.write_record(header, rec)?;
                }
            }
            Ok(())
        }
        OutFmt::Cram => {
            let reference = rmdup_reference(reference).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "CRAM output requires -T/--reference FILE or top-level --reference FILE",
                )
            })?;
            write_cram_output(header, records, keep, output, &reference)
        }
    }
}

fn write_cram_output(
    header: &sam::Header,
    records: &[RecordBuf],
    keep: &[bool],
    output: Option<&Path>,
    reference: &Path,
) -> io::Result<()> {
    crate::reference::ensure_fai_index(reference, None)?;
    let (tmp_file, tmp_path) = crate::tmp_file::create_temp_file("rmdup", Some("bam"))?;
    {
        use sam::alignment::io::Write as _;
        let mut writer = bam::io::Writer::new(tmp_file);
        writer.write_header(header)?;
        for (i, record) in records.iter().enumerate() {
            if keep[i] {
                writer.write_alignment_record(header, record)?;
            }
        }
    }

    let result = match output {
        Some(path) => {
            let out = File::create(path)?;
            htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference(
                tmp_path.path(),
                reference,
                out,
            )
            .map(|_| ())
        }
        None => {
            let out = io::stdout().lock();
            htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference(
                tmp_path.path(),
                reference,
                out,
            )
            .map(|_| ())
        }
    };
    tmp_path.close().ok();
    result
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(
        w,
        "Usage: samtools rmdup [-sS] [-O sam|bam|cram] [-T ref.fa] <in.bam|in.sam|in.cram> [<out>]"
    )?;
    writeln!(
        w,
        "  -s    treat reads as single-end (this port: single-end only)"
    )?;
    writeln!(w, "  -S    treat paired-end as single-end (alias of -s)")?;
    writeln!(w)?;
    writeln!(w, "NOTE: rmdup is deprecated; prefer `samtools markdup`.")?;
    Ok(())
}

fn pos_key(record: &RecordBuf) -> (i32, i64, bool) {
    let flag = record.flags().bits() as u32;
    (
        record
            .reference_sequence_id()
            .map(|t| t as i32)
            .unwrap_or(-1),
        record.alignment_start().map(usize::from).unwrap_or(0) as i64,
        flag & BAM_FREVERSE != 0,
    )
}

fn pair_score(records: &[RecordBuf], first: usize, second: usize) -> u32 {
    u32::from(records[first].mapping_quality().map(u8::from).unwrap_or(0))
        + u32::from(records[second].mapping_quality().map(u8::from).unwrap_or(0))
}
