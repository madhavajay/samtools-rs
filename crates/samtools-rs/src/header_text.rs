//! Helpers for extracting the raw header text from SAM/BAM/CRAM inputs.
//!
//! samtools subcommands (`head`, `view`, `reheader`, ...) emit the header
//! exactly as it was stored in the input file for byte-for-byte parity with
//! upstream. The noodles canonical serializer reorders header lines, so we
//! pull the bytes directly: SAM lines stream verbatim, BAM is read as
//! length-prefixed text after the magic, and CRAM is read from the embedded
//! SAM header container.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{Category, Exact, detect_path};

/// Read the raw header text from a SAM/BAM/CRAM file as it was stored on
/// disk, preserving the original line order.
pub fn read_raw_header_text(path: &Path) -> io::Result<String> {
    let format =
        detect_path(path).map_err(|e| io::Error::other(format!("failed to detect format: {e}")))?;
    if format.category != Category::SequenceData {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{} is not sequence data", path.display()),
        ));
    }
    read_raw_header_text_with_format(path, format.exact)
}

/// Variant of [`read_raw_header_text`] for callers that have already
/// classified the input format and want to skip the redundant detection.
pub fn read_raw_header_text_with_format(path: &Path, exact: Exact) -> io::Result<String> {
    match exact {
        Exact::Sam => read_sam_header_text(path),
        Exact::Bam => read_bam_header_text(path),
        Exact::Cram => htslib_rs::alignment_compat::read_raw_cram_header_text(path),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported format",
        )),
    }
}

fn read_sam_header_text(path: &Path) -> io::Result<String> {
    let file = File::open(path)?;
    let mut reader: Box<dyn BufRead> = if is_bgzf_path(path)? {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut out = String::new();
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        if line.starts_with('@') {
            out.push_str(&line);
        } else {
            break;
        }
    }
    Ok(out)
}

fn read_bam_header_text(path: &Path) -> io::Result<String> {
    let file = File::open(path)?;
    let mut reader = MultiGzDecoder::new(file);
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

fn is_bgzf_path(path: &Path) -> io::Result<bool> {
    let mut file = File::open(path)?;
    let mut hdr = [0u8; 2];
    let n = file.read(&mut hdr)?;
    Ok(n >= 2 && hdr[0] == 0x1f && hdr[1] == 0x8b)
}
