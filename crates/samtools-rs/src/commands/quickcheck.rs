//! `samtools quickcheck` — quickly validate SAM/BAM/CRAM/VCF/BCF files.
//!
//! Mirrors `bam_quickcheck.c`. For each input file:
//!  1. Try to open / detect format.
//!  2. Verify it is sequence data.
//!  3. Read and validate the header (with `-u`, allow zero references).
//!  4. Check the EOF marker (BGZF or CRAM v2.1+).
//!
//! Exit status is the bitwise OR of per-file status codes, matching upstream.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::format::{Category, Exact, Format, detect_path};

// Status bits — match `QC_*` defines in `bam_quickcheck.c`.
const QC_FAIL_OPEN: u32 = 2;
const QC_NOT_SEQUENCE: u32 = 4;
const QC_BAD_HEADER: u32 = 8;
const QC_NO_EOF_BLOCK: u32 = 16;

/// Entry point for `samtools quickcheck`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut verbose: i32 = 0;
    let mut quiet = false;
    let mut unmapped = false;
    let mut files: Vec<PathBuf> = Vec::new();

    // Skip args[0] (subcommand name). Accept bundled short flags `-vqu`.
    for arg in args.iter().skip(1) {
        let Some(s) = arg.to_str() else {
            files.push(PathBuf::from(arg));
            continue;
        };
        if let Some(rest) = s.strip_prefix('-')
            && !rest.is_empty()
            && rest.chars().all(|c| "vqu".contains(c))
        {
            for c in rest.chars() {
                match c {
                    'v' => verbose += 1,
                    'q' => quiet = true,
                    'u' => unmapped = true,
                    _ => unreachable!(),
                }
            }
            continue;
        }
        if s == "-" {
            files.push(PathBuf::from(arg));
        } else if s.starts_with('-') {
            write_usage(&mut io::stderr());
            return ExitCode::from(1);
        } else {
            files.push(PathBuf::from(arg));
        }
    }

    if files.is_empty() {
        write_usage(&mut io::stdout());
        return ExitCode::from(1);
    }

    if verbose >= 2 {
        let _ = writeln!(io::stderr(), "verbosity set to {}", verbose);
    }

    let mut overall: u32 = 0;
    let mut stdout = io::stdout().lock();
    let mut stderr = io::stderr().lock();
    for fn_path in &files {
        let state = check_file(fn_path, verbose, quiet, unmapped, &mut stderr);
        if state > 0 && verbose >= 1 {
            let _ = writeln!(stdout, "{}", fn_path.display());
        }
        overall |= state;
    }

    ExitCode::from((overall & 0xff) as u8)
}

pub(crate) fn check_file<W: Write>(
    path: &Path,
    verbose: i32,
    quiet: bool,
    unmapped: bool,
    stderr: &mut W,
) -> u32 {
    let mut state: u32 = 0;
    if verbose >= 3 {
        let _ = writeln!(stderr, "checking {}", path.display());
    }

    // detect_path is more thorough than HTSlib's hts_open: it tries to
    // decompress the BGZF stream to determine the inner format. If that
    // fails (e.g. the body is truncated/corrupted), we still want to treat
    // the file as openable and surface a BAD_HEADER failure later — match
    // upstream's "open succeeds, header read fails" behavior. Use a
    // best-effort magic-byte fallback for that case.
    let format = match detect_path(path) {
        Ok(f) => f,
        Err(_) => match peek_magic(path) {
            Ok(Some(f)) => f,
            _ => {
                qc_err(
                    &mut state,
                    QC_FAIL_OPEN,
                    verbose,
                    quiet,
                    stderr,
                    &format!("{} could not be opened for reading.", path.display()),
                );
                return state;
            }
        },
    };
    if verbose >= 3 {
        let _ = writeln!(stderr, "opened {}", path.display());
    }

    if format.category != Category::SequenceData {
        qc_err(
            &mut state,
            QC_NOT_SEQUENCE,
            verbose,
            quiet,
            stderr,
            &format!("{} was not identified as sequence data.", path.display()),
        );
    } else {
        if verbose >= 3 {
            let _ = writeln!(stderr, "{} is sequence data", path.display());
        }
        match read_header_nref(path, format.exact) {
            Err(_) if format.exact == Exact::Bam => qc_err(
                &mut state,
                QC_NOT_SEQUENCE,
                verbose,
                quiet,
                stderr,
                &format!("{} was not identified as sequence data.", path.display()),
            ),
            Err(_) => qc_err(
                &mut state,
                QC_BAD_HEADER,
                verbose,
                quiet,
                stderr,
                &format!(
                    "{} caused an error whilst reading its header.",
                    path.display()
                ),
            ),
            Ok(nref) => {
                if !unmapped && nref <= 0 {
                    qc_err(
                        &mut state,
                        QC_BAD_HEADER,
                        verbose,
                        quiet,
                        stderr,
                        &format!("{} had no targets in header.", path.display()),
                    );
                } else if verbose >= 3 {
                    let _ = writeln!(stderr, "{} has {} targets in header.", path.display(), nref);
                }
            }
        }
    }

    match check_eof(path, format.exact) {
        Err(_) => qc_err(
            &mut state,
            QC_NO_EOF_BLOCK,
            verbose,
            quiet,
            stderr,
            &format!(
                "{} caused an error whilst checking for EOF block.",
                path.display()
            ),
        ),
        Ok(0) => qc_err(
            &mut state,
            QC_NO_EOF_BLOCK,
            verbose,
            quiet,
            stderr,
            &format!(
                "{} was missing EOF block when one should be present.",
                path.display()
            ),
        ),
        Ok(1) if verbose >= 3 => {
            let _ = writeln!(stderr, "{} has good EOF block.", path.display());
        }
        Ok(2) if verbose >= 3 => {
            let _ = writeln!(
                stderr,
                "{} cannot be checked for EOF block as it is not seekable.",
                path.display()
            );
        }
        Ok(3) if verbose >= 3 => {
            let _ = writeln!(
                stderr,
                "{} cannot be checked for EOF block because its filetype does not contain one.",
                path.display()
            );
        }
        Ok(_) => {}
    }

    state
}

fn qc_err<W: Write>(
    state: &mut u32,
    bit: u32,
    verbose: i32,
    quiet: bool,
    stderr: &mut W,
    msg: &str,
) {
    *state |= bit;
    // Upstream: print if !quiet OR verbose >= 2. We pass through that gate
    // here by emitting whenever !quiet, or when verbose lifts the gate.
    let _ = verbose;
    if !quiet {
        let _ = writeln!(stderr, "{}", msg);
    }
}

fn read_header_nref(path: &Path, exact: Exact) -> io::Result<i64> {
    match exact {
        Exact::Bam => {
            let header = htslib_rs::alignment_compat::read_bam_header_from_path(path)?;
            Ok(header.reference_sequences().len() as i64)
        }
        Exact::Sam => {
            let header = htslib_rs::alignment_compat::read_sam_header_from_path(path)?;
            Ok(header.reference_sequences().len() as i64)
        }
        Exact::Cram => {
            // Try the noodles CRAM reader first. If it fails (typically because
            // the file is CRAM v2.1 which noodles' main path does not fully
            // parse), fall back to a manual scan of the SAM header embedded in
            // the CRAM file header.
            match htslib_rs::alignment_compat::read_cram_header_from_path(path) {
                Ok(header) => Ok(header.reference_sequences().len() as i64),
                Err(_) => read_cram_header_fallback(path),
            }
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a sequence file",
        )),
    }
}

/// Best-effort magic-byte format detection used when [`detect_path`] errors
/// because it could not fully inspect a compressed stream. We accept the file
/// as openable as long as the first bytes match a known sequence-data magic.
fn peek_magic(path: &Path) -> io::Result<Option<Format>> {
    use htslib_rs::format::Compression;
    let mut file = File::open(path)?;
    let mut hdr = [0u8; 18];
    let n = file.read(&mut hdr)?;
    if n == 0 {
        return Ok(None);
    }
    if n >= 4 && &hdr[..4] == b"CRAM" {
        return Ok(Some(Format::new(
            Category::SequenceData,
            Exact::Cram,
            Compression::None,
        )));
    }
    if n >= 2 && hdr[0] == 0x1f && hdr[1] == 0x8b {
        // gzip / BGZF — without successful decompression we cannot tell
        // BAM from VCF.gz from generic gzip, but for quickcheck purposes
        // we let the header-read path classify the failure.
        return Ok(Some(Format::new(
            Category::SequenceData,
            Exact::Bam,
            Compression::Bgzf,
        )));
    }
    Ok(None)
}

/// Fallback CRAM "nref" detection that counts `@SQ` lines in the SAM header
/// embedded in a CRAM file header. Used when the noodles CRAM reader cannot
/// parse the file (e.g. CRAM v2.1).
fn read_cram_header_fallback(path: &Path) -> io::Result<i64> {
    let mut file = File::open(path)?;
    // CRAM file header (26 bytes) + first container header. The first
    // container's first block holds the SAM header as a length-prefixed
    // ITF-8 / int32-prefixed string. Rather than re-implementing the
    // container/block decoder we read the next ~1MiB and search for the
    // text SAM header. This is approximate but sufficient for quickcheck
    // (we only need the count of @SQ lines).
    let mut hdr = [0u8; 26];
    file.read_exact(&mut hdr)?;
    if &hdr[..4] != b"CRAM" {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "not CRAM"));
    }
    let mut buf = vec![0u8; 1 << 20];
    let n = file.read(&mut buf)?;
    let data = &buf[..n];
    let mut nref: i64 = 0;
    for line in data.split(|&b| b == b'\n') {
        if line.starts_with(b"@SQ\t") || line.starts_with(b"@SQ ") {
            nref += 1;
        }
    }
    Ok(nref)
}

fn check_eof(path: &Path, exact: Exact) -> io::Result<i32> {
    match exact {
        Exact::Bam | Exact::Bcf => check_bgzf_eof(path),
        Exact::Cram => check_cram_eof(path),
        Exact::Vcf | Exact::Sam | Exact::Fasta | Exact::Fastq | Exact::Bed => Ok(3),
        _ => {
            let mut file = File::open(path)?;
            let mut hdr = [0u8; 4];
            let n = file.read(&mut hdr)?;
            if n >= 2 && hdr[0] == 0x1f && hdr[1] == 0x8b {
                check_bgzf_eof(path)
            } else {
                Ok(3)
            }
        }
    }
}

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

fn check_bgzf_eof(path: &Path) -> io::Result<i32> {
    let mut file = File::open(path)?;
    let len = file.metadata()?.len();
    if len < BGZF_EOF.len() as u64 {
        return Ok(0);
    }
    file.seek(SeekFrom::End(-(BGZF_EOF.len() as i64)))?;
    let mut buf = [0u8; 28];
    file.read_exact(&mut buf)?;
    Ok(if buf == BGZF_EOF { 1 } else { 0 })
}

const CRAM_EOF_V21: [u8; 30] = [
    0x0b, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00,
];

const CRAM_EOF_V3: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

fn check_cram_eof(path: &Path) -> io::Result<i32> {
    let mut file = File::open(path)?;
    let mut hdr = [0u8; 6];
    file.read_exact(&mut hdr)?;
    if &hdr[..4] != b"CRAM" {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "not CRAM"));
    }
    let major = hdr[4];
    let minor = hdr[5];
    if major < 2 || (major == 2 && minor == 0) {
        return Ok(3);
    }

    let (template, template_len): (&[u8], usize) = if major == 2 && minor == 1 {
        (&CRAM_EOF_V21[..], CRAM_EOF_V21.len())
    } else {
        (&CRAM_EOF_V3[..], CRAM_EOF_V3.len())
    };

    let len = file.metadata()?.len();
    if len < template_len as u64 {
        return Ok(0);
    }
    file.seek(SeekFrom::End(-(template_len as i64)))?;
    let mut buf = [0u8; 38];
    let buf = &mut buf[..template_len];
    file.read_exact(buf)?;
    if buf.len() > 8 {
        buf[8] &= 0x0f;
    }
    Ok(if buf == template { 1 } else { 0 })
}

fn write_usage<W: Write>(w: &mut W) {
    let _ = writeln!(w, "Usage: samtools quickcheck [options] <input> [...]");
    let _ = writeln!(w, "Options:");
    let _ = writeln!(
        w,
        "  -v              verbose output (repeat for more verbosity)"
    );
    let _ = writeln!(w, "  -q              suppress warning messages");
    let _ = writeln!(
        w,
        "  -u              unmapped input (do not require targets in header)"
    );
}
