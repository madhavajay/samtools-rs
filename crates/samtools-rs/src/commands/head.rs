//! `samtools head` — header (and optionally first N records) of a SAM/BAM/CRAM.
//!
//! Mirrors `main_head` in `sam_view.c`. For byte-for-byte parity with
//! upstream, the header is emitted in the original file order (not the
//! noodles canonical order) by extracting the raw header bytes directly.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use flate2::read::MultiGzDecoder;
use htslib_rs::format::Exact;

use crate::diagnostics::{print_error, print_error_errno};
use crate::header_text::read_raw_header_text;
use crate::io as sam_io;
use crate::sam_global::current_global_args;

/// Entry point for `samtools head`.
pub fn main(args: &[OsString]) -> ExitCode {
    let opts = match parse_args(args) {
        Ok(o) => o,
        Err(ParseError::Usage) => {
            let _ = write_usage(&mut io::stdout());
            return ExitCode::SUCCESS;
        }
        Err(ParseError::Err(msg)) => {
            print_error("head", msg);
            let _ = write_usage(&mut io::stderr());
            return ExitCode::from(1);
        }
    };

    let Some(path) = opts.path.as_ref() else {
        print_error("head", "reading from standard input is not yet supported");
        return ExitCode::from(1);
    };

    let header_text = match read_raw_header_text(path) {
        Ok(t) => t,
        Err(e) => {
            print_error_errno(
                "head",
                format!("failed to read the header from \"{}\"", path.display()),
                &e,
            );
            return ExitCode::from(1);
        }
    };

    let mut stdout = io::stdout().lock();
    if opts.all_headers {
        let _ = stdout.write_all(header_text.as_bytes());
    } else {
        let mut count = 0u64;
        let mut end = 0usize;
        let bytes = header_text.as_bytes();
        while count < opts.nheaders {
            match memchr::memchr(b'\n', &bytes[end..]) {
                Some(i) => {
                    end += i + 1;
                    count += 1;
                }
                None => break,
            }
        }
        if end > 0 {
            let _ = stdout.write_all(&bytes[..end]);
        }
    }

    if opts.nrecords > 0
        && let Err(e) = write_first_records(&mut stdout, path, opts.nrecords)
    {
        print_error_errno("head", "couldn't format record", &e);
        return ExitCode::from(1);
    }

    ExitCode::SUCCESS
}

struct Opts {
    all_headers: bool,
    nheaders: u64,
    nrecords: u64,
    path: Option<PathBuf>,
}

enum ParseError {
    Usage,
    Err(String),
}

fn parse_args(args: &[OsString]) -> Result<Opts, ParseError> {
    let mut all_headers = true;
    let mut nheaders: u64 = 0;
    let mut nrecords: u64 = 0;
    let mut path: Option<PathBuf> = None;

    let mut i = 1;
    while i < args.len() {
        let arg = &args[i];
        let Some(s) = arg.to_str() else {
            if path.is_none() {
                path = Some(PathBuf::from(arg));
                i += 1;
                continue;
            }
            return Err(ParseError::Err("too many positional arguments".to_string()));
        };

        if s == "-h" || s == "--headers" {
            i += 1;
            let v = args
                .get(i)
                .and_then(|a| a.to_str())
                .ok_or_else(|| ParseError::Err("missing value for -h".into()))?;
            all_headers = false;
            nheaders = parse_u64(v)?;
            i += 1;
        } else if let Some(rest) = s.strip_prefix("--headers=") {
            all_headers = false;
            nheaders = parse_u64(rest)?;
            i += 1;
        } else if s == "-n" || s == "--records" {
            i += 1;
            let v = args
                .get(i)
                .and_then(|a| a.to_str())
                .ok_or_else(|| ParseError::Err("missing value for -n".into()))?;
            nrecords = parse_u64(v)?;
            i += 1;
        } else if let Some(rest) = s.strip_prefix("--records=") {
            nrecords = parse_u64(rest)?;
            i += 1;
        } else if s == "--help" {
            return Err(ParseError::Usage);
        } else if s.starts_with('-') && s != "-" {
            return Err(ParseError::Err(format!("unknown option {}", s)));
        } else {
            if path.is_some() {
                return Err(ParseError::Err("too many positional arguments".into()));
            }
            path = Some(PathBuf::from(arg));
            i += 1;
        }
    }

    Ok(Opts {
        all_headers,
        nheaders,
        nrecords,
        path,
    })
}

fn parse_u64(s: &str) -> Result<u64, ParseError> {
    s.parse::<u64>()
        .map_err(|_| ParseError::Err(format!("expected integer, got \"{}\"", s)))
}

fn is_bgzf_path(path: &Path) -> io::Result<bool> {
    let mut file = File::open(path)?;
    let mut hdr = [0u8; 2];
    let n = file.read(&mut hdr)?;
    Ok(n >= 2 && hdr[0] == 0x1f && hdr[1] == 0x8b)
}

/// Write the first `n` alignment records from `path` as SAM text. The header
/// has already been emitted by the caller.
fn write_first_records<W: Write>(out: &mut W, path: &Path, n: u64) -> io::Result<()> {
    let limit = usize::try_from(n).unwrap_or(usize::MAX);
    let format = sam_io::sam_open_format(path)?;
    match format.exact {
        Exact::Sam => stream_sam_records(out, path, limit),
        Exact::Bam => {
            // htslib-rs renders BAM-as-SAM with the noodles canonical header
            // followed by record lines. Strip the header before emitting.
            let text = htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(
                path,
                Some(limit),
            )?;
            let tail = strip_header_lines(text.as_bytes());
            out.write_all(tail)
        }
        Exact::Cram => {
            let reference = current_global_args().reference.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "CRAM record extraction requires --reference",
                )
            })?;
            let text =
                htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                    path,
                    reference,
                    Some(limit),
                )?;
            let tail = strip_header_lines(text.as_bytes());
            out.write_all(tail)
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported format",
        )),
    }
}

fn stream_sam_records<W: Write>(out: &mut W, path: &Path, n: usize) -> io::Result<()> {
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if is_bgzf_path(path)? {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut count = 0usize;
    let mut line = Vec::with_capacity(1024);
    let mut reader = reader;
    while count < n {
        line.clear();
        // Reuse a Vec<u8> rather than String to preserve binary fields.
        let read = read_until_newline(&mut reader, &mut line)?;
        if read == 0 {
            break;
        }
        if line.starts_with(b"@") {
            // Skip remaining header lines (shouldn't happen if header was
            // already consumed, but defensive).
            continue;
        }
        out.write_all(&line)?;
        count += 1;
    }
    Ok(())
}

fn read_until_newline<R: BufRead>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize> {
    let start = buf.len();
    let n = reader.read_until(b'\n', buf)?;
    Ok(buf.len() - start + n - n)
}

fn strip_header_lines(bytes: &[u8]) -> &[u8] {
    let mut tail = bytes;
    while let Some(pos) = memchr::memchr(b'\n', tail) {
        if tail.starts_with(b"@") {
            tail = &tail[pos + 1..];
        } else {
            break;
        }
    }
    tail
}

fn write_usage<W: Write>(w: &mut W) -> io::Result<()> {
    writeln!(w, "Usage: samtools head [OPTION]... [FILE]")?;
    writeln!(w, "Options:")?;
    writeln!(w, "  -h, --headers INT   Display INT header lines [all]")?;
    writeln!(
        w,
        "  -n, --records INT   Display INT alignment record lines [none]"
    )?;
    Ok(())
}
