//! Shared I/O helpers for samtools-style open/close behavior.
//!
//! This is a partial Rust analogue of samtools' common open/close helpers. It
//! centralizes format detection and text-output finalization while fuller
//! HTSlib-style writer state, auto-indexing, and stdout autoflush semantics are
//! built out.

use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use htslib_rs::format::{Category, Compression, Exact, Format, detect_path};

/// Detects an input format and normalizes the error text for command callers.
pub fn sam_open_format(path: &Path) -> io::Result<Format> {
    detect_path(path).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("failed to detect format of \"{}\": {}", path.display(), e),
        )
    })
}

/// Resolves a samtools-style output mode from an explicit format, output path,
/// and default format.
///
/// This is a deliberately small `sam_open_mode` analogue: it resolves the
/// exact format portion of writer state but does not yet create format-specific
/// writers, attach options, or enable auto-indexing.
pub fn sam_open_mode(
    output_path: Option<&Path>,
    explicit: Option<Exact>,
    default: Exact,
) -> io::Result<Format> {
    let exact = explicit
        .or_else(|| output_path.and_then(exact_from_output_extension))
        .unwrap_or(default);
    if exact == Exact::Unknown {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "could not determine output format",
        ));
    }

    Ok(Format::new(
        category_for_exact(exact),
        exact,
        Compression::None,
    ))
}

fn exact_from_output_extension(path: &Path) -> Option<Exact> {
    let ext = path.extension()?.to_str()?.to_ascii_lowercase();
    match ext.as_str() {
        "sam" => Some(Exact::Sam),
        "bam" => Some(Exact::Bam),
        "cram" => Some(Exact::Cram),
        "vcf" => Some(Exact::Vcf),
        "bcf" => Some(Exact::Bcf),
        "bed" => Some(Exact::Bed),
        "fa" | "fasta" => Some(Exact::Fasta),
        "fq" | "fastq" => Some(Exact::Fastq),
        _ => None,
    }
}

fn category_for_exact(exact: Exact) -> Category {
    match exact {
        Exact::Sam | Exact::Bam | Exact::Cram => Category::SequenceData,
        Exact::Vcf | Exact::Bcf => Category::VariantData,
        Exact::Bai
        | Exact::Crai
        | Exact::Csi
        | Exact::Gzi
        | Exact::Tbi
        | Exact::Fai
        | Exact::Fqi => Category::IndexFile,
        Exact::Bed => Category::RegionList,
        _ => Category::Unknown,
    }
}

/// Opens a text output path, or stdout when no path is provided.
pub fn open_text_output(path: Option<&Path>) -> io::Result<Box<dyn Write>> {
    match path {
        Some(path) => Ok(Box::new(File::create(path)?)),
        None => Ok(Box::new(io::stdout().lock())),
    }
}

/// Flushes stdout before reporting an error path that may follow stdout writes.
pub fn autoflush_if_stdout() -> io::Result<()> {
    io::stdout().lock().flush()
}

/// Mirrors the close-check behavior used by samtools: flush buffered output and
/// surface any write error to the caller.
pub fn check_sam_close<W>(writer: &mut W) -> io::Result<()>
where
    W: Write + ?Sized,
{
    writer.flush()
}

/// Writes all bytes, then performs the shared close/flush check.
pub fn write_all_and_close<W>(writer: &mut W, bytes: &[u8]) -> io::Result<()>
where
    W: Write + ?Sized,
{
    writer.write_all(bytes)?;
    check_sam_close(writer)
}

#[cfg(test)]
mod tests {
    use std::io;
    use std::path::Path;

    use htslib_rs::format::Exact;

    use super::{
        check_sam_close, open_text_output, sam_open_format, sam_open_mode, write_all_and_close,
    };

    struct FailingFlush;

    impl io::Write for FailingFlush {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Err(io::Error::other("flush failed"))
        }
    }

    #[test]
    fn close_surfaces_flush_error() {
        let mut writer = FailingFlush;

        let err = check_sam_close(&mut writer).unwrap_err();

        assert_eq!(err.to_string(), "flush failed");
    }

    #[test]
    fn write_all_and_close_writes_to_vec() {
        let mut out = Vec::new();

        write_all_and_close(&mut out, b"abc").unwrap();

        assert_eq!(out, b"abc");
    }

    #[test]
    fn detects_missing_input_as_invalid_input() {
        let err = sam_open_format(Path::new("definitely-missing.sam")).unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert!(err.to_string().contains("failed to detect format"));
    }

    #[test]
    fn opens_text_output_file() {
        let path = std::env::temp_dir().join(format!(
            "samtools-rs-open-text-output-{}.txt",
            std::process::id()
        ));
        let _ = std::fs::remove_file(&path);

        {
            let mut out = open_text_output(Some(&path)).unwrap();
            write_all_and_close(&mut out, b"ok\n").unwrap();
        }

        assert_eq!(std::fs::read_to_string(&path).unwrap(), "ok\n");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn open_mode_uses_explicit_format_first() {
        let format = sam_open_mode(Some(Path::new("out.sam")), Some(Exact::Bam), Exact::Sam)
            .expect("format");

        assert_eq!(format.exact, Exact::Bam);
    }

    #[test]
    fn open_mode_infers_from_extension() {
        let format = sam_open_mode(Some(Path::new("out.CRAM")), None, Exact::Sam).expect("format");

        assert_eq!(format.exact, Exact::Cram);
    }

    #[test]
    fn open_mode_falls_back_to_default() {
        let format =
            sam_open_mode(Some(Path::new("out.unknown")), None, Exact::Sam).expect("format");

        assert_eq!(format.exact, Exact::Sam);
    }

    #[test]
    fn open_mode_rejects_unknown_default() {
        let err = sam_open_mode(None, None, Exact::Unknown).unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }
}
