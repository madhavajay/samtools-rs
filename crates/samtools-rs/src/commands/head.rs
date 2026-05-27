//! `samtools head` — header (and optionally first N records) of a SAM/BAM/CRAM.
//!
//! Mirrors `main_head` in `sam_view.c`. For byte-for-byte parity with
//! upstream, the header is emitted in the original file order (not the
//! noodles canonical order) by extracting the raw header bytes directly.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use flate2::read::MultiGzDecoder;
use htslib_rs::format::Exact;

use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
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

    if opts
        .path
        .as_ref()
        .is_none_or(|path| path.as_os_str() == "-")
    {
        let mut stdin = io::stdin().lock();
        let mut input = Vec::new();
        if let Err(e) = stdin.read_to_end(&mut input) {
            print_error_errno("head", "couldn't read from standard input", &e);
            return ExitCode::from(1);
        }

        let mut stdout = io::stdout().lock();
        if let Err(e) = write_stdin_alignment_head(
            &mut stdout,
            &input,
            opts.all_headers,
            opts.nheaders,
            opts.nrecords,
            current_global_args().reference.as_deref(),
        ) {
            print_error_errno("head", "couldn't read from standard input", &e);
            return ExitCode::from(1);
        }
        return ExitCode::SUCCESS;
    }

    let path = opts.path.as_ref().expect("stdin path handled above");
    if !path.exists() {
        print_hts_open_missing(path);
        print_error(
            "head",
            format!(
                "failed to open \"{}\" for reading: No such file or directory",
                path.display()
            ),
        );
        return ExitCode::from(1);
    }

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

fn write_sam_head_from_reader<W, R>(
    out: &mut W,
    mut reader: R,
    all_headers: bool,
    nheaders: u64,
    nrecords: u64,
) -> io::Result<()>
where
    W: Write,
    R: BufRead,
{
    let mut header_count = 0u64;
    let mut record_count = 0u64;
    let mut in_header = true;
    let mut line = Vec::with_capacity(1024);

    loop {
        line.clear();
        let read = read_until_newline(&mut reader, &mut line)?;
        if read == 0 {
            break;
        }

        if in_header && line.starts_with(b"@") {
            if all_headers || header_count < nheaders {
                out.write_all(&line)?;
            }
            header_count += 1;
            continue;
        }

        in_header = false;
        if record_count >= nrecords {
            break;
        }

        out.write_all(&line)?;
        record_count += 1;
    }

    Ok(())
}

fn write_stdin_alignment_head<W: Write>(
    out: &mut W,
    input: &[u8],
    all_headers: bool,
    nheaders: u64,
    nrecords: u64,
    reference: Option<&Path>,
) -> io::Result<()> {
    let decoded = decode_stdin_alignment_bytes(input)?;

    if decoded.starts_with(b"BAM\x01") {
        let header_text = read_bam_header_text_from_uncompressed_bytes(&decoded)?;
        write_header_text(out, &header_text, all_headers, nheaders)?;

        if nrecords > 0 {
            let limit = usize::try_from(nrecords).unwrap_or(usize::MAX);
            let text =
                htslib_rs::alignment_compat::view_bam_as_sam_text(Cursor::new(input), Some(limit))?;
            out.write_all(strip_header_lines(text.as_bytes()))?;
        }

        return Ok(());
    }

    if input.starts_with(b"CRAM") {
        let header_text = read_cram_header_text_from_bytes(input)?;
        write_header_text(out, &header_text, all_headers, nheaders)?;

        if nrecords > 0 {
            let reference = reference.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "CRAM record extraction requires --reference",
                )
            })?;
            let limit = usize::try_from(nrecords).unwrap_or(usize::MAX);
            let text = htslib_rs::alignment_compat::view_cram_as_sam_text_with_reference(
                Cursor::new(input),
                reference,
                Some(limit),
            )?;
            out.write_all(strip_header_lines(text.as_bytes()))?;
        }

        return Ok(());
    }

    write_sam_head_from_reader(
        out,
        BufReader::new(Cursor::new(decoded)),
        all_headers,
        nheaders,
        nrecords,
    )
}

fn decode_stdin_alignment_bytes(input: &[u8]) -> io::Result<Vec<u8>> {
    if input.starts_with(&[0x1f, 0x8b]) {
        let mut decoder = MultiGzDecoder::new(Cursor::new(input));
        let mut decoded = Vec::new();
        decoder.read_to_end(&mut decoded)?;
        Ok(decoded)
    } else {
        Ok(input.to_vec())
    }
}

fn read_bam_header_text_from_uncompressed_bytes(input: &[u8]) -> io::Result<String> {
    let mut reader = Cursor::new(input);
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if &magic != b"BAM\x01" {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "not BAM"));
    }

    let mut len_bytes = [0u8; 4];
    reader.read_exact(&mut len_bytes)?;
    let l_text = u32::from_le_bytes(len_bytes) as usize;
    let mut text = vec![0u8; l_text];
    reader.read_exact(&mut text)?;
    while text.last() == Some(&0) {
        text.pop();
    }

    String::from_utf8(text).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_cram_header_text_from_bytes(input: &[u8]) -> io::Result<String> {
    let mut reader = htslib_rs::cram::io::Reader::new(Cursor::new(input));
    let mut header_reader = reader.header_reader();

    header_reader.read_magic_number()?;
    header_reader.read_format_version()?;
    header_reader.read_file_id()?;

    let mut container_reader = header_reader.container_reader()?;
    let mut raw_sam_header_reader = container_reader.raw_sam_header_reader()?;
    let mut raw_header = String::new();
    raw_sam_header_reader.read_to_string(&mut raw_header)?;
    raw_sam_header_reader.discard_to_end()?;

    Ok(raw_header)
}

fn write_header_text<W: Write>(
    out: &mut W,
    header_text: &str,
    all_headers: bool,
    nheaders: u64,
) -> io::Result<()> {
    if all_headers {
        return out.write_all(header_text.as_bytes());
    }

    let mut count = 0u64;
    let mut end = 0usize;
    let bytes = header_text.as_bytes();
    while count < nheaders {
        match memchr::memchr(b'\n', &bytes[end..]) {
            Some(i) => {
                end += i + 1;
                count += 1;
            }
            None => break,
        }
    }

    if end > 0 {
        out.write_all(&bytes[..end])?;
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

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::{
        decode_stdin_alignment_bytes, write_sam_head_from_reader, write_stdin_alignment_head,
    };

    fn fixtures_dir() -> PathBuf {
        let manifest = env!("CARGO_MANIFEST_DIR");
        PathBuf::from(manifest)
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("repos")
            .join("samtools")
            .join("test")
            .join("dat")
    }

    fn htslib_fixtures_dir() -> PathBuf {
        let manifest = env!("CARGO_MANIFEST_DIR");
        PathBuf::from(manifest)
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("repos")
            .join("htslib-rs")
            .join("htslib")
            .join("test")
    }

    #[test]
    fn stdin_sam_head_emits_all_headers_only_by_default() {
        let input = b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:8\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\n";
        let mut out = Vec::new();

        write_sam_head_from_reader(&mut out, &input[..], true, 0, 0).unwrap();

        assert_eq!(out, b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:8\n");
    }

    #[test]
    fn stdin_sam_head_limits_headers_and_records() {
        let input = b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:8\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\nr2\t0\tref\t2\t20\t2M\t*\t0\t0\tTG\t##\n";
        let mut out = Vec::new();

        write_sam_head_from_reader(&mut out, &input[..], false, 1, 1).unwrap();

        assert_eq!(
            out,
            b"@HD\tVN:1.6\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\n"
        );
    }

    #[test]
    fn stdin_bam_head_emits_raw_header_and_limited_records() {
        let input = std::fs::read(fixtures_dir().join("test_input_1_a.bam")).unwrap();
        let mut out = Vec::new();

        write_stdin_alignment_head(&mut out, &input, false, 2, 1, None).unwrap();

        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.iter().filter(|line| line.starts_with('@')).count(), 2);
        assert_eq!(
            lines.iter().filter(|line| !line.starts_with('@')).count(),
            1
        );
    }

    #[test]
    fn stdin_cram_head_emits_raw_header_and_limited_records_with_reference() {
        let fixtures = htslib_fixtures_dir();
        let reference = fixtures.join("ce.fa");
        let input = std::fs::read(fixtures.join("range.cram")).unwrap();
        let mut out = Vec::new();

        write_stdin_alignment_head(&mut out, &input, false, 2, 1, Some(&reference)).unwrap();

        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.iter().filter(|line| line.starts_with('@')).count(), 2);
        assert_eq!(
            lines.iter().filter(|line| !line.starts_with('@')).count(),
            1
        );
    }

    #[test]
    fn stdin_cram_records_require_reference() {
        let input = std::fs::read(htslib_fixtures_dir().join("range.cram")).unwrap();
        let mut out = Vec::new();

        let err = write_stdin_alignment_head(&mut out, &input, true, 0, 1, None).unwrap_err();

        assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
        assert!(err.to_string().contains("--reference"));
    }

    #[test]
    fn gzip_sam_stdin_is_decoded_before_head_processing() {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        use std::io::Write as _;

        let input = b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:8\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\n";
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(input).unwrap();
        let compressed = encoder.finish().unwrap();

        let decoded = decode_stdin_alignment_bytes(&compressed).unwrap();

        assert_eq!(decoded, input);
    }
}
