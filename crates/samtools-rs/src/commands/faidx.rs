//! `samtools faidx` — FASTA index build and region extraction.
//!
//! Mirrors `faidx_main` in `faidx.c`. This initial port covers index builds
//! and local uncompressed region extraction. BGZI support and the long tail of
//! formatting options are TODO.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use crate::diagnostics::{print_error, print_error_errno};

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
        line_len: 50,
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
            "-@" | "--threads" | "--output-fmt-opt" => {
                i += 1;
                args.get(i)
                    .ok_or_else(|| format!("missing value for {s}"))?;
                i += 1;
            }
            "--continue" => {
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
        match build_index(input, opts.fai_path.as_deref(), opts.is_fastq) {
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

fn build_index(input: &Path, fai_path: Option<&Path>, is_fastq: bool) -> std::io::Result<()> {
    let file = File::open(input)?;
    let reader = BufReader::new(file);
    let fai_out = fai_path
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| append_extension(input, "fai"));
    let out_file = File::create(&fai_out)?;
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

    let mut out: Box<dyn Write> = match opts.output.as_deref() {
        Some(path) => Box::new(File::create(path).map_err(|e| e.to_string())?),
        None => Box::new(io::stdout().lock()),
    };

    if opts.is_fastq {
        let index = read_or_build_fastq_index(input, opts.fai_path.as_deref())?;
        let mut reader = File::open(input).map_err(|e| e.to_string())?;
        for region in regions {
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
                        should_include_length(&region),
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
        let index = read_or_build_index(input, opts.fai_path.as_deref())?;
        let mut reader = File::open(input).map_err(|e| e.to_string())?;
        for region in regions {
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
                        should_include_length(&region),
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

    Ok(())
}

fn read_or_build_index(
    input: &Path,
    fai_path: Option<&Path>,
) -> Result<htslib_rs::faidx_compat::Index, String> {
    let path = fai_path
        .map(Path::to_path_buf)
        .unwrap_or_else(|| append_extension(input, "fai"));
    if path.exists() {
        let file = File::open(&path).map_err(|e| e.to_string())?;
        return htslib_rs::faidx_compat::read_index(BufReader::new(file))
            .map_err(|e| e.to_string());
    }

    build_index(input, Some(&path), false).map_err(|e| e.to_string())?;
    let file = File::open(&path).map_err(|e| e.to_string())?;
    htslib_rs::faidx_compat::read_index(BufReader::new(file)).map_err(|e| e.to_string())
}

fn read_or_build_fastq_index(
    input: &Path,
    fai_path: Option<&Path>,
) -> Result<htslib_rs::faidx_compat::FastqIndex, String> {
    let path = fai_path
        .map(Path::to_path_buf)
        .unwrap_or_else(|| append_extension(input, "fai"));
    if path.exists() {
        let file = File::open(&path).map_err(|e| e.to_string())?;
        return htslib_rs::faidx_compat::read_fastq_index(BufReader::new(file))
            .map_err(|e| e.to_string());
    }

    build_index(input, Some(&path), true).map_err(|e| e.to_string())?;
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

fn should_include_length(region: &str) -> bool {
    region.contains(':')
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
