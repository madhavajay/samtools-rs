//! `samtools faidx` — FASTA index build and region extraction.
//!
//! Mirrors `faidx_main` in `faidx.c`. This initial port covers index builds
//! and local uncompressed region extraction. BGZI support and the long tail of
//! formatting options are TODO.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools faidx`.
pub fn main(args: &[OsString]) -> ExitCode {
    match run_faidx(args, false) {
        Ok(code) => code,
        Err(e) => {
            print_error("faidx", e);
            ExitCode::from(1)
        }
    }
}

/// Entry point for `samtools fqidx` (FASTQ index build).
pub fn fqidx_main(args: &[OsString]) -> ExitCode {
    match run_faidx(args, true) {
        Ok(code) => code,
        Err(e) => {
            print_error("fqidx", e);
            ExitCode::from(1)
        }
    }
}

fn run_faidx(args: &[OsString], is_fastq: bool) -> Result<ExitCode, String> {
    let mut opts = Opts {
        is_fastq,
        line_len: 60,
        ..Opts::default()
    };
    let mut i = 1;
    while i < args.len() {
        let Some(s) = args[i].to_str() else {
            opts.positional.push(PathBuf::from(&args[i]));
            i += 1;
            continue;
        };
        match s {
            "--fai-idx" => {
                i += 1;
                let v = args
                    .get(i)
                    .map(PathBuf::from)
                    .ok_or_else(|| "missing value for --fai-idx".to_string())?;
                opts.fai_path = Some(v);
                i += 1;
            }
            "--gzi-idx" => {
                i += 1;
                let v = args
                    .get(i)
                    .map(PathBuf::from)
                    .ok_or_else(|| "missing value for --gzi-idx".to_string())?;
                opts.gzi_path = Some(v);
                i += 1;
            }
            "-f" | "--fastq" => {
                opts.is_fastq = true;
                i += 1;
            }
            "-r" | "--region-file" => {
                i += 1;
                let v = args
                    .get(i)
                    .map(PathBuf::from)
                    .ok_or_else(|| "missing value for -r".to_string())?;
                opts.region_file = Some(v);
                i += 1;
            }
            "-o" | "--output" => {
                i += 1;
                let v = args
                    .get(i)
                    .map(PathBuf::from)
                    .ok_or_else(|| "missing value for -o".to_string())?;
                opts.output = Some(v);
                i += 1;
            }
            "--length" | "-n" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| "missing value for --length".to_string())?;
                opts.line_len = v
                    .parse::<usize>()
                    .map_err(|_| "invalid value for --length".to_string())?;
                if opts.line_len == 0 {
                    return Err("invalid value for --length".to_string());
                }
                i += 1;
            }
            "--write-index" => {
                opts.write_index = true;
                i += 1;
            }
            "-i" | "--reverse-complement" => {
                opts.reverse_complement = true;
                i += 1;
            }
            "--mark-strand" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| "missing value for --mark-strand".to_string())?;
                opts.mark_strand = parse_mark_strand(v)?;
                i += 1;
            }
            "-@" | "--threads" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| format!("missing value for {s}"))?;
                opts.threads = Some(parse_thread_count(v)?);
                i += 1;
            }
            _ if s.starts_with("-@") && s.len() > 2 => {
                opts.threads = Some(parse_thread_count(&s[2..])?);
                i += 1;
            }
            _ if s.starts_with("--threads=") => {
                let v = s
                    .split_once('=')
                    .map(|(_, value)| value)
                    .unwrap_or_default();
                opts.threads = Some(parse_thread_count(v)?);
                i += 1;
            }
            "--output-fmt-opt" => {
                i += 1;
                args.get(i)
                    .ok_or_else(|| format!("missing value for {s}"))?;
                i += 1;
            }
            _ if s.starts_with("--output-fmt-opt=") => {
                i += 1;
            }
            "--continue" => {
                opts.continue_on_missing = true;
                i += 1;
            }
            "-c" => {
                opts.continue_on_missing = true;
                i += 1;
            }
            "--help" | "-h" => {
                let _ = print_usage(opts.is_fastq);
                return Ok(ExitCode::SUCCESS);
            }
            _ if s.starts_with('-') && s != "-" => {
                return Err(format!(
                    "option `{}` is not yet supported in samtools-rs faidx",
                    s
                ));
            }
            _ => {
                opts.positional.push(PathBuf::from(&args[i]));
                i += 1;
            }
        }
    }

    if opts.positional.is_empty() {
        let _ = print_usage(opts.is_fastq);
        return Ok(ExitCode::SUCCESS);
    }

    // Build-only path: a single positional arg and no region file.
    if opts.positional.len() == 1 && opts.region_file.is_none() {
        let input = &opts.positional[0];
        if !input.exists() {
            eprintln!(
                "[E::fai_build3_core] Failed to open the file {} : No such file or directory",
                input.display()
            );
            eprintln!("[faidx] Could not build fai index {}.fai", input.display());
            return Ok(ExitCode::from(1));
        }
        match build_index(
            input,
            opts.fai_path.as_deref(),
            opts.gzi_path.as_deref(),
            opts.is_fastq,
            opts.worker_count(),
        ) {
            Ok(()) => Ok(ExitCode::SUCCESS),
            Err(e) => {
                print_error_errno(
                    if opts.is_fastq { "fqidx" } else { "faidx" },
                    format!("Could not build fai index {}.fai", input.display()),
                    &e,
                );
                Ok(ExitCode::from(1))
            }
        }
    } else {
        let input = opts.positional.remove(0);
        opts.regions.extend(
            opts.positional
                .iter()
                .map(|p| p.to_string_lossy().into_owned()),
        );
        run_retrieval(&input, &opts).map(|()| ExitCode::SUCCESS)
    }
}

#[derive(Default)]
struct Opts {
    fai_path: Option<PathBuf>,
    gzi_path: Option<PathBuf>,
    region_file: Option<PathBuf>,
    output: Option<PathBuf>,
    regions: Vec<String>,
    positional: Vec<PathBuf>,
    line_len: usize,
    is_fastq: bool,
    write_index: bool,
    continue_on_missing: bool,
    reverse_complement: bool,
    mark_strand: MarkStrand,
    threads: Option<usize>,
}

impl Opts {
    fn worker_count(&self) -> Option<NonZero<usize>> {
        self.threads
            .or_else(|| current_global_args().threads)
            .and_then(NonZero::new)
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
enum MarkStrand {
    #[default]
    Rc,
    Sign,
    No,
    Custom {
        forward: String,
        reverse: String,
    },
}

fn build_index(
    input: &Path,
    fai_path: Option<&Path>,
    gzi_path: Option<&Path>,
    is_fastq: bool,
    worker_count: Option<NonZero<usize>>,
) -> std::io::Result<()> {
    let src = std::fs::read(input)?;
    let fai_out = fai_path
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| append_extension(input, "fai"));
    let out_file = File::create(&fai_out)?;

    let reader: Box<dyn BufRead> = match htslib_rs::bgzf_compat::detect_compression_kind(&src) {
        htslib_rs::bgzf_compat::CompressionKind::Uncompressed => {
            Box::new(BufReader::new(Cursor::new(src)))
        }
        htslib_rs::bgzf_compat::CompressionKind::Gzip => {
            let data = htslib_rs::bgzf_compat::read_auto_with_worker_count(&src, worker_count)?;
            Box::new(BufReader::new(Cursor::new(data)))
        }
        htslib_rs::bgzf_compat::CompressionKind::Bgzf => {
            let mut cursor = Cursor::new(&src);
            let gzi = htslib_rs::bgzf_compat::build_gzi(&mut cursor)?;
            let gzi_out = gzi_path
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| append_extension(input, "gzi"));
            htslib_rs::bgzf_compat::write_gzi_to_path(gzi_out, &gzi)?;

            let data = htslib_rs::bgzf_compat::read_auto_with_worker_count(&src, worker_count)?;
            Box::new(BufReader::new(Cursor::new(data)))
        }
    };

    if is_fastq {
        let index = htslib_rs::faidx_compat::build_fastq_index(reader)?;
        htslib_rs::faidx_compat::write_fastq_index(out_file, &index)?;
    } else {
        let index = htslib_rs::faidx_compat::build_index(reader)?;
        htslib_rs::faidx_compat::write_index(out_file, &index)?;
    }
    Ok(())
}

fn run_retrieval(input: &Path, opts: &Opts) -> Result<(), String> {
    let mut regions = opts.regions.clone();
    if let Some(path) = opts.region_file.as_deref() {
        regions.extend(read_region_file(path)?);
    }

    if regions.is_empty() {
        return Err("no regions specified".into());
    }

    let mut out = Vec::new();
    let worker_count = opts.worker_count();

    if opts.is_fastq {
        let index = read_or_build_fastq_index(input, opts.fai_path.as_deref(), worker_count)?;
        let mut reader = read_indexed_source(input, worker_count).map_err(|e| e.to_string())?;
        for region in regions {
            if fastq_region_is_missing(&index, &region) && opts.continue_on_missing {
                warn_missing_reference(&region);
                warn_missing_reference(&region);
                warn_failed_fetch(&region);
                warn_missing_reference(&region);
                warn_failed_fetch(&region);
                let name = format_region_name(&region, opts.reverse_complement, &opts.mark_strand);
                write_retrieval_record(&mut out, b'@', &name, &[], Some(&[]), opts.line_len, false)
                    .map_err(|e| e.to_string())?;
                continue;
            }
            if let Some(action) = fastq_retrieval_bounds_action(&index, &region) {
                match action {
                    RetrievalBoundsAction::Zero => {
                        warn_zero_length_region(&region);
                        warn_zero_length_region(&region);
                        let name =
                            format_region_name(&region, opts.reverse_complement, &opts.mark_strand);
                        write_retrieval_record(
                            &mut out,
                            b'@',
                            &name,
                            &[],
                            Some(&[]),
                            opts.line_len,
                            false,
                        )
                        .map_err(|e| e.to_string())?;
                        continue;
                    }
                    RetrievalBoundsAction::Truncated => {
                        warn_truncated_region(&region);
                        warn_truncated_region(&region);
                    }
                }
            }

            let sequence =
                htslib_rs::faidx_compat::fetch_fastq_region_sequence(&mut reader, &index, &region);
            let quality =
                htslib_rs::faidx_compat::fetch_fastq_region_quality(&mut reader, &index, &region);
            match (sequence, quality) {
                (Ok(mut sequence), Ok(mut quality)) => {
                    let name =
                        format_region_name(&region, opts.reverse_complement, &opts.mark_strand);
                    if opts.reverse_complement {
                        reverse_complement_in_place(&mut sequence);
                        quality.reverse();
                    }
                    write_retrieval_record(
                        &mut out,
                        b'@',
                        &name,
                        &sequence,
                        Some(&quality),
                        opts.line_len,
                        false,
                    )
                    .map_err(|e| e.to_string())?;
                }
                (Err(e), _) | (_, Err(e)) if opts.continue_on_missing => {
                    print_error("fqidx", format!("failed to retrieve \"{region}\": {e}"));
                }
                (Err(e), _) | (_, Err(e)) => return Err(e.to_string()),
            }
        }
    } else {
        let index = read_or_build_index(input, opts.fai_path.as_deref(), worker_count)?;
        let mut reader = read_indexed_source(input, worker_count).map_err(|e| e.to_string())?;
        for region in regions {
            if region_is_missing(&index, &region) && opts.continue_on_missing {
                warn_missing_reference(&region);
                warn_missing_reference(&region);
                warn_failed_fetch(&region);
                let name = format_region_name(&region, opts.reverse_complement, &opts.mark_strand);
                write_retrieval_record(&mut out, b'>', &name, &[], None, opts.line_len, false)
                    .map_err(|e| e.to_string())?;
                continue;
            }
            if let Some(action) = retrieval_bounds_action(&index, &region) {
                match action {
                    RetrievalBoundsAction::Zero => {
                        warn_zero_length_region(&region);
                        let name =
                            format_region_name(&region, opts.reverse_complement, &opts.mark_strand);
                        write_retrieval_record(
                            &mut out,
                            b'>',
                            &name,
                            &[],
                            None,
                            opts.line_len,
                            false,
                        )
                        .map_err(|e| e.to_string())?;
                        continue;
                    }
                    RetrievalBoundsAction::Truncated => warn_truncated_region(&region),
                }
            }

            match htslib_rs::faidx_compat::fetch_region_sequence(&mut reader, &index, &region) {
                Ok(mut sequence) => {
                    let name =
                        format_region_name(&region, opts.reverse_complement, &opts.mark_strand);
                    if opts.reverse_complement {
                        reverse_complement_in_place(&mut sequence);
                    }
                    write_retrieval_record(
                        &mut out,
                        b'>',
                        &name,
                        &sequence,
                        None,
                        opts.line_len,
                        false,
                    )
                    .map_err(|e| e.to_string())?;
                }
                Err(e) if opts.continue_on_missing => {
                    print_error("faidx", format!("failed to retrieve \"{region}\": {e}"));
                }
                Err(e) => return Err(e.to_string()),
            }
        }
    }

    write_output(opts.output.as_deref(), &out, worker_count).map_err(|e| e.to_string())?;

    if opts.write_index
        && let Some(output) = opts.output.as_deref()
    {
        build_index(output, None, None, opts.is_fastq, worker_count).map_err(|e| e.to_string())?;
    }

    Ok(())
}

fn read_or_build_index(
    input: &Path,
    fai_path: Option<&Path>,
    worker_count: Option<NonZero<usize>>,
) -> Result<htslib_rs::faidx_compat::Index, String> {
    let path = fai_path
        .map(Path::to_path_buf)
        .unwrap_or_else(|| append_extension(input, "fai"));
    if path.exists() {
        let file = File::open(&path).map_err(|e| e.to_string())?;
        return htslib_rs::faidx_compat::read_index(BufReader::new(file))
            .map_err(|e| e.to_string());
    }

    build_index(input, Some(&path), None, false, worker_count).map_err(|e| e.to_string())?;
    let file = File::open(&path).map_err(|e| e.to_string())?;
    htslib_rs::faidx_compat::read_index(BufReader::new(file)).map_err(|e| e.to_string())
}

fn read_or_build_fastq_index(
    input: &Path,
    fai_path: Option<&Path>,
    worker_count: Option<NonZero<usize>>,
) -> Result<htslib_rs::faidx_compat::FastqIndex, String> {
    let path = fai_path
        .map(Path::to_path_buf)
        .unwrap_or_else(|| append_extension(input, "fai"));
    if path.exists() {
        let file = File::open(&path).map_err(|e| e.to_string())?;
        return htslib_rs::faidx_compat::read_fastq_index(BufReader::new(file))
            .map_err(|e| e.to_string());
    }

    build_index(input, Some(&path), None, true, worker_count).map_err(|e| e.to_string())?;
    let file = File::open(&path).map_err(|e| e.to_string())?;
    htslib_rs::faidx_compat::read_fastq_index(BufReader::new(file)).map_err(|e| e.to_string())
}

fn read_region_file(path: &Path) -> Result<Vec<String>, String> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let reader = BufReader::new(file);
    let mut regions = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(region) = trimmed.split_ascii_whitespace().next() {
            regions.push(region.to_string());
        }
    }

    Ok(regions)
}

fn parse_thread_count(raw: &str) -> Result<usize, String> {
    raw.parse::<usize>()
        .map_err(|_| format!("invalid thread count \"{raw}\""))
}

fn read_indexed_source(
    input: &Path,
    worker_count: Option<NonZero<usize>>,
) -> io::Result<Cursor<Vec<u8>>> {
    let src = std::fs::read(input)?;
    let data = htslib_rs::bgzf_compat::read_auto_with_worker_count(&src, worker_count)?;

    Ok(Cursor::new(data))
}

fn write_output(
    output: Option<&Path>,
    data: &[u8],
    worker_count: Option<NonZero<usize>>,
) -> io::Result<()> {
    if let Some(path) = output {
        if is_bgzf_output_path(path) {
            let encoded = htslib_rs::bgzf_compat::write_all_with_kind_and_worker_count(
                data,
                htslib_rs::bgzf_compat::CompressionKind::Bgzf,
                worker_count,
            )?;
            std::fs::write(path, encoded)
        } else {
            let mut out = sam_io::open_text_output(Some(path))?;
            out.write_all(data)?;
            sam_io::check_sam_close(&mut out)
        }
    } else {
        let mut out = sam_io::open_text_output(None)?;
        out.write_all(data)?;
        sam_io::check_sam_close(&mut out)
    }
}

enum RetrievalBoundsAction {
    Zero,
    Truncated,
}

fn retrieval_bounds_action(
    index: &htslib_rs::faidx_compat::Index,
    region: &str,
) -> Option<RetrievalBoundsAction> {
    let parsed = parse_region_bounds(region)?;
    let len = htslib_rs::faidx_compat::sequence_len(index, parsed.name)?;

    bounds_action(parsed, len)
}

fn fastq_retrieval_bounds_action(
    index: &htslib_rs::faidx_compat::FastqIndex,
    region: &str,
) -> Option<RetrievalBoundsAction> {
    let parsed = parse_region_bounds(region)?;
    let len = htslib_rs::faidx_compat::fastq_sequence_len(index, parsed.name)?;

    bounds_action(parsed, len)
}

fn region_is_missing(index: &htslib_rs::faidx_compat::Index, region: &str) -> bool {
    htslib_rs::faidx_compat::sequence_len(index, region_reference_name(region)).is_none()
}

fn fastq_region_is_missing(index: &htslib_rs::faidx_compat::FastqIndex, region: &str) -> bool {
    htslib_rs::faidx_compat::fastq_sequence_len(index, region_reference_name(region)).is_none()
}

fn region_reference_name(region: &str) -> &str {
    region
        .rsplit_once(':')
        .map(|(name, _)| name)
        .unwrap_or(region)
}

fn bounds_action(parsed: ParsedRegionBounds<'_>, len: u64) -> Option<RetrievalBoundsAction> {
    let start = parsed.start.unwrap_or(1);
    let end = parsed.end.unwrap_or(len);

    if start == 0 || start > len || start > end {
        Some(RetrievalBoundsAction::Zero)
    } else if end > len {
        Some(RetrievalBoundsAction::Truncated)
    } else {
        None
    }
}

struct ParsedRegionBounds<'a> {
    name: &'a str,
    start: Option<u64>,
    end: Option<u64>,
}

fn parse_region_bounds(region: &str) -> Option<ParsedRegionBounds<'_>> {
    let (name, interval) = region.rsplit_once(':')?;
    if name.is_empty() {
        return None;
    }

    let (start, end) = interval
        .split_once('-')
        .map_or((interval, ""), |(start, end)| (start, end));

    Some(ParsedRegionBounds {
        name,
        start: parse_region_position(start),
        end: parse_region_position(end),
    })
}

fn parse_region_position(value: &str) -> Option<u64> {
    if value.is_empty() {
        return None;
    }

    value
        .bytes()
        .filter(|b| *b != b',')
        .try_fold(0_u64, |n, b| {
            if b.is_ascii_digit() {
                Some(n.saturating_mul(10).saturating_add(u64::from(b - b'0')))
            } else {
                None
            }
        })
}

fn warn_zero_length_region(region: &str) {
    eprintln!("[faidx] Zero length sequence: {region}");
}

fn warn_truncated_region(region: &str) {
    eprintln!("[faidx] Truncated sequence: {region}");
}

fn warn_missing_reference(region: &str) {
    eprintln!(
        "[W::fai_get_val] Reference {region} not found in FASTA file, returning empty sequence"
    );
}

fn warn_failed_fetch(region: &str) {
    eprintln!("[faidx] Failed to fetch sequence in {region}");
}

fn is_bgzf_output_path(path: &Path) -> bool {
    matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("gz" | "bgz" | "bgzf")
    )
}

fn write_retrieval_record<W>(
    writer: &mut W,
    marker: u8,
    region: &str,
    sequence: &[u8],
    quality: Option<&[u8]>,
    line_len: usize,
    include_length: bool,
) -> io::Result<()>
where
    W: Write + ?Sized,
{
    if include_length {
        writeln!(
            writer,
            "{}{} length: {}",
            char::from(marker),
            region,
            sequence.len()
        )?;
    } else {
        writeln!(writer, "{}{}", char::from(marker), region)?;
    }
    write_wrapped(writer, sequence, line_len)?;

    if let Some(quality) = quality {
        writeln!(writer, "+")?;
        write_wrapped(writer, quality, line_len)?;
    }

    Ok(())
}

fn parse_mark_strand(value: &str) -> Result<MarkStrand, String> {
    match value {
        "rc" => Ok(MarkStrand::Rc),
        "sign" => Ok(MarkStrand::Sign),
        "no" => Ok(MarkStrand::No),
        _ => {
            let mut fields = value.split(',').map(str::trim);
            match (fields.next(), fields.next(), fields.next(), fields.next()) {
                (Some("custom"), Some(forward), Some(reverse), None)
                    if !forward.is_empty() && !reverse.is_empty() =>
                {
                    Ok(MarkStrand::Custom {
                        forward: forward.to_string(),
                        reverse: reverse.to_string(),
                    })
                }
                _ => Err(format!("invalid --mark-strand value `{value}`")),
            }
        }
    }
}

fn format_region_name(region: &str, reverse_complement: bool, mark_strand: &MarkStrand) -> String {
    if !reverse_complement {
        if let MarkStrand::Custom { forward, .. } = mark_strand {
            let _ = forward;
        }
        return region.to_string();
    }

    match mark_strand {
        MarkStrand::Rc => format!("{region}/rc"),
        MarkStrand::Sign => format!("{region}(-)"),
        MarkStrand::No => region.to_string(),
        MarkStrand::Custom { forward, reverse } => {
            let _ = forward;
            format!("{region} {reverse}")
        }
    }
}

fn reverse_complement_in_place(sequence: &mut [u8]) {
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

fn write_wrapped<W>(writer: &mut W, data: &[u8], line_len: usize) -> io::Result<()>
where
    W: Write + ?Sized,
{
    for chunk in data.chunks(line_len) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn append_extension(path: &Path, ext: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

fn print_usage(is_fastq: bool) -> std::io::Result<()> {
    use std::io::Write as _;
    let mut w = std::io::stdout().lock();
    let cmd = if is_fastq { "fqidx" } else { "faidx" };
    writeln!(
        w,
        "Usage: samtools {} <file.fa|file.fa.gz> [<region> [...]]",
        cmd
    )?;
    writeln!(
        w,
        "  Builds <file.fa>.fai (and emits sequences for regions, TODO)."
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::ffi::OsString;
    use std::fs;
    use std::path::PathBuf;

    use crate::sam_global::{SamGlobalArgs, set_current_global_args};

    use super::run_faidx;

    fn argv(args: &[impl AsRef<str>]) -> Vec<OsString> {
        args.iter().map(|s| OsString::from(s.as_ref())).collect()
    }

    fn tmp_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn faidx_local_threads_use_bgzf_worker_output_path() {
        let dir = tmp_dir("samtools-rs-faidx-local-threads");
        let input = dir.join("ref.fa");
        let output = dir.join("out.fa.bgz");
        fs::write(&input, b">sq\nACGTACGT\n").unwrap();

        run_faidx(&argv(&["faidx", input.to_str().unwrap()]), false).unwrap();
        run_faidx(
            &argv(&[
                "faidx",
                "-@2",
                "-o",
                output.to_str().unwrap(),
                input.to_str().unwrap(),
                "sq:1-4",
            ]),
            false,
        )
        .unwrap();

        let encoded = fs::read(&output).unwrap();
        assert_eq!(
            htslib_rs::bgzf_compat::detect_compression_kind(&encoded),
            htslib_rs::bgzf_compat::CompressionKind::Bgzf
        );
        assert_eq!(
            htslib_rs::bgzf_compat::read_auto_with_worker_count(
                &encoded,
                std::num::NonZero::new(2)
            )
            .unwrap(),
            b">sq:1-4\nACGT\n"
        );

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn faidx_global_threads_feed_bgzf_output_path() {
        let dir = tmp_dir("samtools-rs-faidx-global-threads");
        let input = dir.join("ref.fa");
        let output = dir.join("out.fa.bgz");
        fs::write(&input, b">sq\nACGTACGT\n").unwrap();
        run_faidx(&argv(&["faidx", input.to_str().unwrap()]), false).unwrap();

        set_current_global_args(SamGlobalArgs {
            threads: Some(2),
            ..SamGlobalArgs::default()
        });
        let result = run_faidx(
            &argv(&[
                "faidx",
                "-o",
                output.to_str().unwrap(),
                input.to_str().unwrap(),
                "sq:5-8",
            ]),
            false,
        );
        set_current_global_args(SamGlobalArgs::default());

        result.unwrap();
        let encoded = fs::read(&output).unwrap();
        assert_eq!(
            htslib_rs::bgzf_compat::read_auto_with_worker_count(
                &encoded,
                std::num::NonZero::new(2)
            )
            .unwrap(),
            b">sq:5-8\nACGT\n"
        );

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn fqidx_threads_use_bgzf_worker_output_path() {
        let dir = tmp_dir("samtools-rs-fqidx-threads");
        let input = dir.join("reads.fq");
        let local_output = dir.join("local.fq.bgz");
        let global_output = dir.join("global.fq.bgz");
        fs::write(&input, b"@r1\nACGTACGT\n+\nabcdefgh\n").unwrap();

        run_faidx(&argv(&["fqidx", input.to_str().unwrap()]), true).unwrap();
        run_faidx(
            &argv(&[
                "fqidx",
                "-@2",
                "-o",
                local_output.to_str().unwrap(),
                input.to_str().unwrap(),
                "r1:1-4",
            ]),
            true,
        )
        .unwrap();

        set_current_global_args(SamGlobalArgs {
            threads: Some(2),
            ..SamGlobalArgs::default()
        });
        let result = run_faidx(
            &argv(&[
                "fqidx",
                "-o",
                global_output.to_str().unwrap(),
                input.to_str().unwrap(),
                "r1:5-8",
            ]),
            true,
        );
        set_current_global_args(SamGlobalArgs::default());
        result.unwrap();

        let local_encoded = fs::read(&local_output).unwrap();
        let global_encoded = fs::read(&global_output).unwrap();
        for encoded in [&local_encoded, &global_encoded] {
            assert_eq!(
                htslib_rs::bgzf_compat::detect_compression_kind(encoded),
                htslib_rs::bgzf_compat::CompressionKind::Bgzf
            );
        }

        assert_eq!(
            htslib_rs::bgzf_compat::read_auto_with_worker_count(
                &local_encoded,
                std::num::NonZero::new(2)
            )
            .unwrap(),
            b"@r1:1-4\nACGT\n+\nabcd\n"
        );
        assert_eq!(
            htslib_rs::bgzf_compat::read_auto_with_worker_count(
                &global_encoded,
                std::num::NonZero::new(2)
            )
            .unwrap(),
            b"@r1:5-8\nACGT\n+\nefgh\n"
        );

        fs::remove_dir_all(dir).unwrap();
    }
}
