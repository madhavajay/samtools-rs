//! `samtools index` — build a BAI/CSI/CRAI index for a SAM/BAM/CRAM file.
//!
//! Mirrors `bam_index.c`. Supports:
//!  - `-b` / `--bai`   — BAI (default for BAM)
//!  - `-c` / `--csi`   — CSI
//!  - `-m` / `--min-shift INT` — CSI min-shift (also implies `-c`)
//!  - `-M` — interpret all arguments as files to index (vs legacy `<in> <out.idx>`)
//!  - `-o` / `--output FILE` — write index to FILE
//!  - `-@` / `--threads INT` — accepted, not yet wired through

use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::format::{Category, Exact, detect_path};

use crate::diagnostics::{print_error, print_error_errno};

/// Entry point for `samtools index`.
pub fn main(args: &[OsString]) -> ExitCode {
    let opts = match parse_args(args) {
        Ok(o) => o,
        Err(ParseError::Usage(code)) => {
            print_usage(if code == 0 {
                Stream::Stdout
            } else {
                Stream::Stderr
            });
            return ExitCode::from(code);
        }
        Err(ParseError::Err(msg)) => {
            print_error("index", msg);
            print_usage(Stream::Stderr);
            return ExitCode::from(1);
        }
    };

    let mut inputs = opts.inputs.clone();
    let mut explicit_idx = opts.output.clone();

    // Legacy synopsis: `samtools index <in> <out.idx>`.
    if inputs.len() == 2 && explicit_idx.is_none() && nonexistent_or_index(&inputs[1]) {
        explicit_idx = Some(inputs.pop().unwrap());
    }

    if inputs.len() > 1 && !opts.multiple {
        print_error(
            "index",
            "use -M to enable indexing more than one alignment file",
        );
        return ExitCode::from(1);
    }
    if explicit_idx.is_some() && inputs.len() > 1 {
        print_error("index", "can't use -o with multiple input alignment files");
        return ExitCode::from(1);
    }

    for src in &inputs {
        if let Err(e) = build_index(src, explicit_idx.as_deref(), opts.csi, opts.min_shift) {
            print_error_errno(
                "index",
                format!("failed to create index for \"{}\"", src.display()),
                &e,
            );
            return ExitCode::from(1);
        }
    }
    ExitCode::SUCCESS
}

#[derive(Default)]
struct Opts {
    csi: bool,
    min_shift: Option<u8>,
    multiple: bool,
    output: Option<PathBuf>,
    inputs: Vec<PathBuf>,
}

enum ParseError {
    Usage(u8),
    Err(String),
}

fn parse_args(args: &[OsString]) -> Result<Opts, ParseError> {
    let mut opts = Opts::default();
    let mut i = 1;
    while i < args.len() {
        let Some(s) = args[i].to_str() else {
            opts.inputs.push(PathBuf::from(&args[i]));
            i += 1;
            continue;
        };
        match s {
            "-b" | "--bai" => {
                opts.csi = false;
                i += 1;
            }
            "-c" | "--csi" => {
                opts.csi = true;
                i += 1;
            }
            "-M" => {
                opts.multiple = true;
                i += 1;
            }
            "-m" | "--min-shift" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("invalid min_shift".into()))?;
                let parsed: u8 = v
                    .parse()
                    .map_err(|_| ParseError::Err("invalid min_shift".into()))?;
                opts.csi = true;
                opts.min_shift = Some(parsed);
                i += 1;
            }
            "-o" | "--output" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| ParseError::Err("missing value for -o".into()))?;
                opts.output = Some(PathBuf::from(v));
                i += 1;
            }
            "-@" | "--threads" => {
                i += 1;
                args.get(i)
                    .and_then(|a| a.to_str())
                    .and_then(|s| s.parse::<u32>().ok())
                    .ok_or_else(|| ParseError::Err("invalid thread count".into()))?;
                i += 1;
            }
            "--help" => return Err(ParseError::Usage(0)),
            _ if s.starts_with('-') && s != "-" => {
                return Err(ParseError::Err(format!("unknown option {}", s)));
            }
            _ => {
                opts.inputs.push(PathBuf::from(&args[i]));
                i += 1;
            }
        }
    }
    if opts.inputs.is_empty() {
        return Err(ParseError::Usage(0));
    }
    Ok(opts)
}

fn nonexistent_or_index(path: &Path) -> bool {
    if !path.exists() {
        return true;
    }
    matches!(
        detect_path(path).map(|f| f.category),
        Ok(Category::IndexFile)
    )
}

fn build_index(
    src: &Path,
    fn_idx: Option<&Path>,
    csi: bool,
    min_shift: Option<u8>,
) -> std::io::Result<()> {
    let format = detect_path(src)
        .map_err(|e| std::io::Error::other(format!("failed to detect format: {e}")))?;
    match format.exact {
        Exact::Bam => {
            if csi {
                let index = match min_shift {
                    Some(m) => htslib_rs::index_compat::build_bam_csi_with_min_shift(src, m)?,
                    None => htslib_rs::index_compat::build_bam_csi(src)?,
                };
                let out = fn_idx
                    .map(|p| p.to_path_buf())
                    .unwrap_or_else(|| append_extension(src, "csi"));
                htslib_rs::index_compat::write_csi(out, &index)?;
            } else {
                let index = htslib_rs::index_compat::build_bai(src)?;
                let out = fn_idx
                    .map(|p| p.to_path_buf())
                    .unwrap_or_else(|| append_extension(src, "bai"));
                htslib_rs::index_compat::write_bai(out, &index)?;
            }
            Ok(())
        }
        Exact::Sam => {
            if csi {
                let index = match min_shift {
                    Some(m) => htslib_rs::index_compat::build_sam_csi_with_min_shift(src, m)?,
                    None => htslib_rs::index_compat::build_sam_csi(src)?,
                };
                let out = fn_idx
                    .map(|p| p.to_path_buf())
                    .unwrap_or_else(|| append_extension(src, "csi"));
                htslib_rs::index_compat::write_csi(out, &index)?;
            } else {
                let index = htslib_rs::index_compat::build_sam_bai(src)?;
                let out = fn_idx
                    .map(|p| p.to_path_buf())
                    .unwrap_or_else(|| append_extension(src, "bai"));
                htslib_rs::index_compat::write_bai(out, &index)?;
            }
            Ok(())
        }
        Exact::Cram => {
            let index = htslib_rs::index_compat::build_cram_crai(src)?;
            let out = fn_idx
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| append_extension(src, "crai"));
            htslib_rs::index_compat::write_cram_crai(out, &index)?;
            Ok(())
        }
        _ => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "{} is in a format that cannot be usefully indexed",
                src.display()
            ),
        )),
    }
}

fn append_extension(path: &Path, ext: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

enum Stream {
    Stdout,
    Stderr,
}

fn print_usage(stream: Stream) {
    use std::io::Write as _;
    let mut buf: Box<dyn std::io::Write> = match stream {
        Stream::Stdout => Box::new(std::io::stdout().lock()),
        Stream::Stderr => Box::new(std::io::stderr().lock()),
    };
    let _ = writeln!(
        buf,
        "Usage: samtools index -M [-bc] [-m INT] <in1.bam> <in2.bam>..."
    );
    let _ = writeln!(
        buf,
        "   or: samtools index [-bc] [-m INT] <in.bam> [out.index]"
    );
    let _ = writeln!(buf, "Options:");
    let _ = writeln!(
        buf,
        "  -b, --bai            Generate BAI-format index for BAM files [default]"
    );
    let _ = writeln!(
        buf,
        "  -c, --csi            Generate CSI-format index for BAM files"
    );
    let _ = writeln!(
        buf,
        "  -m, --min-shift INT  Set minimum interval size for CSI indices to 2^INT [14]"
    );
    let _ = writeln!(
        buf,
        "  -M                   Interpret all filename arguments as files to be indexed"
    );
    let _ = writeln!(
        buf,
        "  -o, --output FILE    Write index to FILE [alternative to <out.index> in args]"
    );
    let _ = writeln!(
        buf,
        "  -@, --threads INT    Sets the number of additional threads [0]"
    );
}
