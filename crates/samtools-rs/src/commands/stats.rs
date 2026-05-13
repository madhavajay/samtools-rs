//! `samtools stats` — alignment statistics summary.
//!
//! Mirrors `stats.c` (123k LOC). Upstream produces a large blob with many
//! `SN`/`FFQ`/`LFQ`/`COV`/`GCF`/etc. sections. This initial Rust port emits
//! only the basic `SN` summary numbers that can be computed from
//! `AlignmentRecordSummary` without per-base access or pileup.
//!
//! **Pending:** quality histograms (FFQ/LFQ), insert size distributions
//! (IS), GC content (GCF), coverage histograms (COV), per-cycle stats,
//! BAQ adjustments, region restriction, reference-based mismatch counts.

use std::ffi::OsString;
use std::io::{self, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use htslib_rs::alignment_compat::AlignmentRecordSummary;
use htslib_rs::format::{Exact, detect_path};

use crate::bam_flag::{
    BAM_FDUP, BAM_FMUNMAP, BAM_FPAIRED, BAM_FPROPER_PAIR, BAM_FQCFAIL, BAM_FREAD1, BAM_FREAD2,
    BAM_FSECONDARY, BAM_FUNMAP,
};
use crate::diagnostics::{print_error, print_error_errno};
use crate::version::SAMTOOLS_VERSION;

/// Entry point for `samtools stats`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "-@" | "--threads" | "-r" | "--reference" | "-t" | "--target-regions" | "-l"
            | "--read-length" | "-I" | "--id" | "-S" | "--split" | "-P" | "--split-prefix"
            | "-g" | "--cov-threshold" | "-G" => {
                let _ = iter.next();
            }
            "-d" | "--remove-dups" | "-s" | "--sparse" | "-x" | "--sam" | "--remove-overlaps"
            | "--no-PG" => {
                // Accepted but not yet implemented.
            }
            "--help" | "-h" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("stats", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
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
            print_error(
                "stats",
                format!("failed to detect format of \"{}\": {}", input.display(), e),
            );
            return ExitCode::from(1);
        }
    };

    let summaries = match format.exact {
        Exact::Sam => htslib_rs::alignment_compat::summarize_sam_records_from_path(&input),
        Exact::Bam => htslib_rs::alignment_compat::summarize_bam_records_from_path(&input),
        _ => {
            print_error(
                "stats",
                "only SAM and BAM input are currently supported (CRAM TODO)",
            );
            return ExitCode::from(1);
        }
    };
    let summaries = match summaries {
        Ok(v) => v,
        Err(e) => {
            print_error_errno(
                "stats",
                format!("error reading from \"{}\"", input.display()),
                &e,
            );
            return ExitCode::from(1);
        }
    };

    let mut writer: Box<dyn Write> = match output.as_ref() {
        Some(p) => match std::fs::File::create(p) {
            Ok(f) => Box::new(f),
            Err(e) => {
                print_error_errno("stats", "open -o output", &e);
                return ExitCode::from(1);
            }
        },
        None => Box::new(io::stdout().lock()),
    };
    if let Err(e) = write_stats(&mut writer, &summaries) {
        print_error_errno("stats", "write failed", &e);
        return ExitCode::from(1);
    }
    ExitCode::SUCCESS
}

fn write_stats(out: &mut dyn Write, recs: &[AlignmentRecordSummary]) -> io::Result<()> {
    let mut total = 0u64;
    let mut mapped = 0u64;
    let mut unmapped = 0u64;
    let mut paired = 0u64;
    let mut proper_paired = 0u64;
    let mut read1 = 0u64;
    let mut read2 = 0u64;
    let mut mapped_and_paired = 0u64;
    let mut dup = 0u64;
    let mut qc_fail = 0u64;
    let mut secondary = 0u64;
    let mut mq0 = 0u64;
    let mut singletons = 0u64;
    let total_length: u64 = 0;
    let max_length: u64 = 0;
    let mut diffchr: u64 = 0;

    for rec in recs {
        let flag = rec.flags_u16() as u32;
        if flag & BAM_FSECONDARY != 0 {
            secondary += 1;
            continue;
        }
        total += 1;

        // We don't have sequence length on AlignmentRecordSummary directly;
        // skip per-record length contributions for now.
        let _ = (&total_length, &max_length);

        if flag & BAM_FUNMAP == 0 {
            mapped += 1;
            if let Some(q) = rec.mapping_quality()
                && q == 0
            {
                mq0 += 1;
            }
        } else {
            unmapped += 1;
        }
        if flag & BAM_FPAIRED != 0 {
            paired += 1;
            if flag & BAM_FREAD1 != 0 {
                read1 += 1;
            }
            if flag & BAM_FREAD2 != 0 {
                read2 += 1;
            }
            if flag & BAM_FPROPER_PAIR != 0 {
                proper_paired += 1;
            }
            if flag & BAM_FUNMAP == 0 && flag & BAM_FMUNMAP == 0 {
                mapped_and_paired += 1;
                if rec.reference_sequence_id() != rec.mate_reference_sequence_id()
                    && rec.reference_sequence_id().is_some()
                    && rec.mate_reference_sequence_id().is_some()
                {
                    diffchr += 1;
                }
            }
            if flag & BAM_FMUNMAP != 0 && flag & BAM_FUNMAP == 0 {
                singletons += 1;
            }
        }
        if flag & BAM_FDUP != 0 {
            dup += 1;
        }
        if flag & BAM_FQCFAIL != 0 {
            qc_fail += 1;
        }
    }

    writeln!(
        out,
        "# This file was produced by samtools-rs stats (samtools-{}+htslib-rs)",
        SAMTOOLS_VERSION
    )?;
    writeln!(out, "# This file contains statistics for all reads.")?;
    writeln!(out, "SN\traw total sequences:\t{}", total)?;
    writeln!(out, "SN\tfiltered sequences:\t0")?;
    writeln!(out, "SN\tsequences:\t{}", total)?;
    writeln!(out, "SN\t1st fragments:\t{}", read1)?;
    writeln!(out, "SN\tlast fragments:\t{}", read2)?;
    writeln!(out, "SN\treads mapped:\t{}", mapped)?;
    writeln!(out, "SN\treads mapped and paired:\t{}", mapped_and_paired)?;
    writeln!(out, "SN\treads unmapped:\t{}", unmapped)?;
    writeln!(out, "SN\treads properly paired:\t{}", proper_paired)?;
    writeln!(out, "SN\treads paired:\t{}", paired)?;
    writeln!(out, "SN\treads duplicated:\t{}", dup)?;
    writeln!(out, "SN\treads MQ0:\t{}", mq0)?;
    writeln!(out, "SN\treads QC failed:\t{}", qc_fail)?;
    writeln!(out, "SN\tnon-primary alignments:\t{}", secondary)?;
    writeln!(out, "SN\tsupplementary alignments:\t0")?;
    writeln!(out, "SN\tsingletons:\t{}", singletons)?;
    writeln!(out, "SN\tpairs on different chromosomes:\t{}", diffchr / 2)?;
    Ok(())
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools stats [options] <in.bam>")?;
    writeln!(w, "  -o FILE      output FILE")?;
    writeln!(w)?;
    writeln!(
        w,
        "Note: only the `SN` summary lines are currently produced."
    )?;
    Ok(())
}
