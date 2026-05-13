//! `samtools idxstats` — print per-reference mapped/unmapped counts.
//!
//! Mirrors `bam_idxstats` in `bam_index.c`. Output format:
//!
//! ```text
//! <chr>\t<len>\t<mapped>\t<unmapped>\n
//! ...
//! *\t0\t0\t<unplaced-unmapped>\n
//! ```
//!
//! Indexed BAM inputs use the associated BAI/CSI metadata. SAM and unindexed
//! BAM inputs fall back to a streaming "slow" pass.

use std::ffi::OsString;
use std::io::{self, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use htslib_rs::csi::BinningIndex;
use htslib_rs::format::{Exact, detect_path};

use crate::bam_flag::BAM_FUNMAP;
use crate::diagnostics::{print_error, print_error_errno};
use crate::header_text::read_raw_header_text_with_format;

/// Entry point for `samtools idxstats`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut input: Option<PathBuf> = None;
    let mut explicit_index: Option<PathBuf> = None;
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-X" => {
                // legacy: include explicit index path as the next arg
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("idxstats", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                } else if explicit_index.is_none() {
                    explicit_index = Some(PathBuf::from(arg));
                }
            }
        }
    }

    let Some(input) = input else {
        let _ = print_usage();
        return ExitCode::from(1);
    };

    let format = match detect_path(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error("idxstats", format!("failed to detect format: {}", e));
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam) {
        print_error(
            "idxstats",
            "only SAM and BAM input is currently supported (CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    let header_text = match read_raw_header_text_with_format(&input, format.exact) {
        Ok(t) => t,
        Err(e) => {
            print_error_errno(
                "idxstats",
                format!("failed to read header for \"{}\"", input.display()),
                &e,
            );
            return ExitCode::from(1);
        }
    };

    // Collect (name, length) for each @SQ line, in original order.
    let refs: Vec<(String, i64)> = header_text
        .lines()
        .filter(|l| l.starts_with("@SQ\t"))
        .filter_map(parse_sq_line)
        .collect();

    let mut stdout = io::stdout().lock();
    if format.exact == Exact::Bam {
        match read_bam_index(&input, explicit_index.as_ref()) {
            Ok(index) => {
                let _ = write_index_stats(&mut stdout, &refs, index.as_ref());
                return ExitCode::SUCCESS;
            }
            Err(e) => {
                if explicit_index.is_some() {
                    print_error_errno(
                        "idxstats",
                        format!("failed to load index for \"{}\"", input.display()),
                        &e,
                    );
                    return ExitCode::from(1);
                }
            }
        }
    }

    let summaries = match format.exact {
        Exact::Sam => htslib_rs::alignment_compat::summarize_sam_records_from_path(&input),
        Exact::Bam => htslib_rs::alignment_compat::summarize_bam_records_from_path(&input),
        _ => unreachable!(),
    };
    let summaries = match summaries {
        Ok(v) => v,
        Err(e) => {
            print_error_errno(
                "idxstats",
                format!("error reading from \"{}\"", input.display()),
                &e,
            );
            return ExitCode::from(1);
        }
    };

    let (counts, no_coor) = slow_counts(
        refs.len(),
        summaries
            .iter()
            .map(|r| (r.flags_u16() as u32, r.reference_sequence_id())),
    );
    let _ = write_counts(&mut stdout, &refs, &counts, no_coor);

    ExitCode::SUCCESS
}

fn read_bam_index(
    input: &PathBuf,
    explicit_index: Option<&PathBuf>,
) -> io::Result<Box<dyn BinningIndex>> {
    if let Some(path) = explicit_index {
        let is_csi = path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case("csi"));
        if is_csi {
            htslib_rs::index_compat::read_csi(path).map(|index| Box::new(index) as _)
        } else {
            htslib_rs::index_compat::read_bai(path).map(|index| Box::new(index) as _)
        }
    } else {
        htslib_rs::index_compat::read_associated_bam_index(input)
    }
}

fn write_index_stats<W>(
    writer: &mut W,
    refs: &[(String, i64)],
    index: &dyn BinningIndex,
) -> io::Result<()>
where
    W: Write + ?Sized,
{
    let ref_metadata: Vec<(u64, u64)> = index
        .reference_sequences()
        .map(|rs| {
            rs.metadata()
                .map(|m| (m.mapped_record_count(), m.unmapped_record_count()))
                .unwrap_or((0, 0))
        })
        .collect();
    let no_coor = index.unplaced_unmapped_record_count().unwrap_or(0);
    write_counts(writer, refs, &ref_metadata, no_coor)
}

fn write_counts<W>(
    writer: &mut W,
    refs: &[(String, i64)],
    counts: &[(u64, u64)],
    no_coor: u64,
) -> io::Result<()>
where
    W: Write + ?Sized,
{
    for (i, (name, len)) in refs.iter().enumerate() {
        let (mapped, unmapped) = counts.get(i).copied().unwrap_or((0, 0));
        writeln!(writer, "{}\t{}\t{}\t{}", name, len, mapped, unmapped)?;
    }
    writeln!(writer, "*\t0\t0\t{}", no_coor)?;

    Ok(())
}

fn slow_counts<I>(reference_count: usize, records: I) -> (Vec<(u64, u64)>, u64)
where
    I: IntoIterator<Item = (u32, Option<usize>)>,
{
    let mut counts = vec![(0, 0); reference_count];
    let mut no_coor = 0;

    for (flags, tid) in records {
        let is_unmapped = flags & BAM_FUNMAP != 0;
        match tid.and_then(|i| counts.get_mut(i)) {
            Some((_, unmapped)) if is_unmapped => *unmapped += 1,
            Some((mapped, _)) => *mapped += 1,
            None if is_unmapped => no_coor += 1,
            None => {}
        }
    }

    (counts, no_coor)
}

fn parse_sq_line(line: &str) -> Option<(String, i64)> {
    let mut sn: Option<&str> = None;
    let mut ln: Option<i64> = None;
    for field in line.split('\t').skip(1) {
        if let Some(v) = field.strip_prefix("SN:") {
            sn = Some(v);
        } else if let Some(v) = field.strip_prefix("LN:") {
            ln = v.parse().ok();
        }
    }
    Some((sn?.to_string(), ln?))
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools idxstats [options] <in.bam>")?;
    writeln!(
        w,
        "  -X           Include customized index file (next positional)"
    )?;
    writeln!(w, "  -@ N         Number of additional threads (TODO)")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{slow_counts, write_counts};
    use crate::bam_flag::BAM_FUNMAP;

    #[test]
    fn slow_counts_separates_mapped_unmapped_and_no_coordinate() {
        let (counts, no_coor) = slow_counts(
            2,
            [
                (0, Some(0)),
                (BAM_FUNMAP, Some(0)),
                (0, Some(1)),
                (BAM_FUNMAP, None),
                (0, None),
            ],
        );

        assert_eq!(counts, vec![(1, 1), (1, 0)]);
        assert_eq!(no_coor, 1);
    }

    #[test]
    fn write_counts_matches_idxstats_text_shape() {
        let refs = vec![("chr1".to_string(), 8), ("chr2".to_string(), 4)];
        let counts = vec![(2, 1), (0, 3)];
        let mut out = Vec::new();

        write_counts(&mut out, &refs, &counts, 5).unwrap();

        assert_eq!(
            String::from_utf8(out).unwrap(),
            "chr1\t8\t2\t1\nchr2\t4\t0\t3\n*\t0\t0\t5\n"
        );
    }
}
