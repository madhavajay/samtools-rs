//! `samtools flagstat` — flag statistics for a SAM/BAM/CRAM file.
//!
//! Mirrors `bam_flagstat` in `bam_stat.c`. Output format is the upstream
//! "default" format unless `-O json|tsv` is given.

use std::ffi::OsString;
use std::io::{self, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use htslib_rs::format::Exact;

use crate::bam_flag::{
    BAM_FDUP, BAM_FMUNMAP, BAM_FPAIRED, BAM_FPROPER_PAIR, BAM_FQCFAIL, BAM_FREAD1, BAM_FREAD2,
    BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP,
};
use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;
use crate::sam_global::current_global_args;

#[derive(Default, Clone, Copy, PartialEq, Eq)]
enum OutFmt {
    #[default]
    Default,
    Json,
    Tsv,
}

#[derive(Default)]
struct Counts {
    // [QC-passed, QC-failed]
    n_reads: [u64; 2],
    n_mapped: [u64; 2],
    n_pair_all: [u64; 2],
    n_pair_map: [u64; 2],
    n_pair_good: [u64; 2],
    n_sgltn: [u64; 2],
    n_read1: [u64; 2],
    n_read2: [u64; 2],
    n_dup: [u64; 2],
    n_diffchr: [u64; 2],
    n_diffhigh: [u64; 2],
    n_secondary: [u64; 2],
    n_supp: [u64; 2],
    n_primary: [u64; 2],
    n_pmapped: [u64; 2],
    n_pdup: [u64; 2],
}

/// Entry point for `samtools flagstat`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut out_fmt = OutFmt::Default;
    let mut input: Option<PathBuf> = None;
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-O" | "--output-fmt" => {
                let v = iter.next().and_then(|a| a.to_str()).unwrap_or("");
                out_fmt = match v.to_lowercase().as_str() {
                    "json" => OutFmt::Json,
                    "tsv" => OutFmt::Tsv,
                    _ => OutFmt::Default,
                };
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("flagstat", format!("unknown option {}", s));
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
        let _ = writeln!(io::stderr(), "Usage: samtools flagstat [options] <in.bam>");
        return ExitCode::from(1);
    };

    let format = match sam_io::sam_open_format(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error("flagstat", format!("Cannot open input file: {}", e));
            return ExitCode::from(1);
        }
    };

    let summaries = match format.exact {
        Exact::Sam => htslib_rs::alignment_compat::summarize_sam_records_from_path(&input),
        Exact::Bam => htslib_rs::alignment_compat::summarize_bam_records_from_path(&input),
        Exact::Cram => {
            // flagstat only inspects flags (reference-independent in
            // CRAM), so a reference is optional — fall back to the
            // synthesizing path, matching `samtools flagstat foo.cram`.
            match current_global_args().reference {
                Some(reference) => {
                    htslib_rs::alignment_compat::summarize_cram_records_from_path_with_reference(
                        &input, reference,
                    )
                }
                None => {
                    htslib_rs::alignment_compat::summarize_cram_records_from_path_synthesizing_reference(
                        &input,
                    )
                }
            }
        }
        _ => {
            print_error(
                "flagstat",
                "only SAM, BAM, and CRAM input is currently supported",
            );
            return ExitCode::from(1);
        }
    };
    let summaries = match summaries {
        Ok(v) => v,
        Err(e) => {
            print_error_errno(
                "flagstat",
                format!("error reading from \"{}\"", input.display()),
                &e,
            );
            return ExitCode::from(1);
        }
    };

    let mut counts = Counts::default();
    for rec in &summaries {
        let flag = rec.flags_u16() as u32;
        let qual = rec.mapping_quality().unwrap_or(0);
        let tid = rec.reference_sequence_id();
        let mtid = rec.mate_reference_sequence_id();
        accumulate(&mut counts, flag, qual, tid, mtid);
    }

    let mut stdout = io::stdout().lock();
    match out_fmt {
        OutFmt::Default => write_default(&mut stdout, &counts),
        OutFmt::Json => write_json(&mut stdout, &counts),
        OutFmt::Tsv => write_tsv(&mut stdout, &counts),
    };
    ExitCode::SUCCESS
}

fn accumulate(s: &mut Counts, flag: u32, qual: u8, tid: Option<usize>, mtid: Option<usize>) {
    let w = if flag & BAM_FQCFAIL != 0 { 1 } else { 0 };
    s.n_reads[w] += 1;
    if flag & BAM_FSECONDARY != 0 {
        s.n_secondary[w] += 1;
    } else if flag & BAM_FSUPPLEMENTARY != 0 {
        s.n_supp[w] += 1;
    } else {
        s.n_primary[w] += 1;
        if flag & BAM_FPAIRED != 0 {
            s.n_pair_all[w] += 1;
            if flag & BAM_FPROPER_PAIR != 0 && flag & BAM_FUNMAP == 0 {
                s.n_pair_good[w] += 1;
            }
            if flag & BAM_FREAD1 != 0 {
                s.n_read1[w] += 1;
            }
            if flag & BAM_FREAD2 != 0 {
                s.n_read2[w] += 1;
            }
            if flag & BAM_FMUNMAP != 0 && flag & BAM_FUNMAP == 0 {
                s.n_sgltn[w] += 1;
            }
            if flag & BAM_FUNMAP == 0 && flag & BAM_FMUNMAP == 0 {
                s.n_pair_map[w] += 1;
                if tid != mtid && tid.is_some() && mtid.is_some() {
                    s.n_diffchr[w] += 1;
                    if qual >= 5 {
                        s.n_diffhigh[w] += 1;
                    }
                }
            }
        }
        if flag & BAM_FUNMAP == 0 {
            s.n_pmapped[w] += 1;
        }
        if flag & BAM_FDUP != 0 {
            s.n_pdup[w] += 1;
        }
    }
    if flag & BAM_FUNMAP == 0 {
        s.n_mapped[w] += 1;
    }
    if flag & BAM_FDUP != 0 {
        s.n_dup[w] += 1;
    }
}

fn percent_f32(n: u64, total: u64) -> String {
    // C uses `(float)n / total * 100.0` then `%.2f%%`. Promote via f32 to
    // match upstream precision.
    if total == 0 {
        "N/A".to_string()
    } else {
        let v = (n as f32) / (total as f32) * 100.0;
        format!("{:.2}%", v)
    }
}

fn write_default<W: Write>(w: &mut W, s: &Counts) {
    let _ = writeln!(
        w,
        "{} + {} in total (QC-passed reads + QC-failed reads)",
        s.n_reads[0], s.n_reads[1]
    );
    let _ = writeln!(w, "{} + {} primary", s.n_primary[0], s.n_primary[1]);
    let _ = writeln!(w, "{} + {} secondary", s.n_secondary[0], s.n_secondary[1]);
    let _ = writeln!(w, "{} + {} supplementary", s.n_supp[0], s.n_supp[1]);
    let _ = writeln!(w, "{} + {} duplicates", s.n_dup[0], s.n_dup[1]);
    let _ = writeln!(w, "{} + {} primary duplicates", s.n_pdup[0], s.n_pdup[1]);
    let _ = writeln!(
        w,
        "{} + {} mapped ({} : {})",
        s.n_mapped[0],
        s.n_mapped[1],
        percent_f32(s.n_mapped[0], s.n_reads[0]),
        percent_f32(s.n_mapped[1], s.n_reads[1])
    );
    let _ = writeln!(
        w,
        "{} + {} primary mapped ({} : {})",
        s.n_pmapped[0],
        s.n_pmapped[1],
        percent_f32(s.n_pmapped[0], s.n_primary[0]),
        percent_f32(s.n_pmapped[1], s.n_primary[1])
    );
    let _ = writeln!(
        w,
        "{} + {} paired in sequencing",
        s.n_pair_all[0], s.n_pair_all[1]
    );
    let _ = writeln!(w, "{} + {} read1", s.n_read1[0], s.n_read1[1]);
    let _ = writeln!(w, "{} + {} read2", s.n_read2[0], s.n_read2[1]);
    let _ = writeln!(
        w,
        "{} + {} properly paired ({} : {})",
        s.n_pair_good[0],
        s.n_pair_good[1],
        percent_f32(s.n_pair_good[0], s.n_pair_all[0]),
        percent_f32(s.n_pair_good[1], s.n_pair_all[1])
    );
    let _ = writeln!(
        w,
        "{} + {} with itself and mate mapped",
        s.n_pair_map[0], s.n_pair_map[1]
    );
    let _ = writeln!(
        w,
        "{} + {} singletons ({} : {})",
        s.n_sgltn[0],
        s.n_sgltn[1],
        percent_f32(s.n_sgltn[0], s.n_pair_all[0]),
        percent_f32(s.n_sgltn[1], s.n_pair_all[1])
    );
    let _ = writeln!(
        w,
        "{} + {} with mate mapped to a different chr",
        s.n_diffchr[0], s.n_diffchr[1]
    );
    let _ = writeln!(
        w,
        "{} + {} with mate mapped to a different chr (mapQ>=5)",
        s.n_diffhigh[0], s.n_diffhigh[1]
    );
}

fn write_tsv<W: Write>(w: &mut W, s: &Counts) {
    let _ = writeln!(
        w,
        "{}\t{}\ttotal (QC-passed reads + QC-failed reads)",
        s.n_reads[0], s.n_reads[1]
    );
    let _ = writeln!(w, "{}\t{}\tprimary", s.n_primary[0], s.n_primary[1]);
    let _ = writeln!(w, "{}\t{}\tsecondary", s.n_secondary[0], s.n_secondary[1]);
    let _ = writeln!(w, "{}\t{}\tsupplementary", s.n_supp[0], s.n_supp[1]);
    let _ = writeln!(w, "{}\t{}\tduplicates", s.n_dup[0], s.n_dup[1]);
    let _ = writeln!(w, "{}\t{}\tprimary duplicates", s.n_pdup[0], s.n_pdup[1]);
    let _ = writeln!(w, "{}\t{}\tmapped", s.n_mapped[0], s.n_mapped[1]);
    let _ = writeln!(
        w,
        "{}\t{}\tmapped %",
        percent_f32(s.n_mapped[0], s.n_reads[0]),
        percent_f32(s.n_mapped[1], s.n_reads[1])
    );
    let _ = writeln!(w, "{}\t{}\tprimary mapped", s.n_pmapped[0], s.n_pmapped[1]);
    let _ = writeln!(
        w,
        "{}\t{}\tprimary mapped %",
        percent_f32(s.n_pmapped[0], s.n_primary[0]),
        percent_f32(s.n_pmapped[1], s.n_primary[1])
    );
    let _ = writeln!(
        w,
        "{}\t{}\tpaired in sequencing",
        s.n_pair_all[0], s.n_pair_all[1]
    );
    let _ = writeln!(w, "{}\t{}\tread1", s.n_read1[0], s.n_read1[1]);
    let _ = writeln!(w, "{}\t{}\tread2", s.n_read2[0], s.n_read2[1]);
    let _ = writeln!(
        w,
        "{}\t{}\tproperly paired",
        s.n_pair_good[0], s.n_pair_good[1]
    );
    let _ = writeln!(
        w,
        "{}\t{}\tproperly paired %",
        percent_f32(s.n_pair_good[0], s.n_pair_all[0]),
        percent_f32(s.n_pair_good[1], s.n_pair_all[1])
    );
    let _ = writeln!(
        w,
        "{}\t{}\twith itself and mate mapped",
        s.n_pair_map[0], s.n_pair_map[1]
    );
    let _ = writeln!(w, "{}\t{}\tsingletons", s.n_sgltn[0], s.n_sgltn[1]);
    let _ = writeln!(
        w,
        "{}\t{}\tsingletons %",
        percent_f32(s.n_sgltn[0], s.n_pair_all[0]),
        percent_f32(s.n_sgltn[1], s.n_pair_all[1])
    );
    let _ = writeln!(
        w,
        "{}\t{}\twith mate mapped to a different chr",
        s.n_diffchr[0], s.n_diffchr[1]
    );
    let _ = writeln!(
        w,
        "{}\t{}\twith mate mapped to a different chr (mapQ>=5)",
        s.n_diffhigh[0], s.n_diffhigh[1]
    );
}

fn write_json<W: Write>(w: &mut W, s: &Counts) {
    // Match upstream's JSON shape character-for-character.
    let _ = write!(
        w,
        "{{\n \"QC-passed reads\": {{ \n  \"total\": {}, \n  \"primary\": {}, \n  \"secondary\": {}, \n  \"supplementary\": {}, \n  \"duplicates\": {}, \n  \"primary duplicates\": {}, \n  \"mapped\": {}, \n  \"mapped %\": {}, \n  \"primary mapped\": {}, \n  \"primary mapped %\": {}, \n  \"paired in sequencing\": {}, \n  \"read1\": {}, \n  \"read2\": {}, \n  \"properly paired\": {}, \n  \"properly paired %\": {}, \n  \"with itself and mate mapped\": {}, \n  \"singletons\": {}, \n  \"singletons %\": {}, \n  \"with mate mapped to a different chr\": {}, \n  \"with mate mapped to a different chr (mapQ >= 5)\": {} \n }},\n \"QC-failed reads\": {{ \n  \"total\": {}, \n  \"primary\": {}, \n  \"secondary\": {}, \n  \"supplementary\": {}, \n  \"duplicates\": {}, \n  \"primary duplicates\": {}, \n  \"mapped\": {}, \n  \"mapped %\": {}, \n  \"primary mapped\": {}, \n  \"primary mapped %\": {}, \n  \"paired in sequencing\": {}, \n  \"read1\": {}, \n  \"read2\": {}, \n  \"properly paired\": {}, \n  \"properly paired %\": {}, \n  \"with itself and mate mapped\": {}, \n  \"singletons\": {}, \n  \"singletons %\": {}, \n  \"with mate mapped to a different chr\": {}, \n  \"with mate mapped to a different chr (mapQ >= 5)\": {} \n }}\n}}\n",
        s.n_reads[0],
        s.n_primary[0],
        s.n_secondary[0],
        s.n_supp[0],
        s.n_dup[0],
        s.n_pdup[0],
        s.n_mapped[0],
        json_pct(s.n_mapped[0], s.n_reads[0]),
        s.n_pmapped[0],
        json_pct(s.n_pmapped[0], s.n_primary[0]),
        s.n_pair_all[0],
        s.n_read1[0],
        s.n_read2[0],
        s.n_pair_good[0],
        json_pct(s.n_pair_good[0], s.n_pair_all[0]),
        s.n_pair_map[0],
        s.n_sgltn[0],
        json_pct(s.n_sgltn[0], s.n_pair_all[0]),
        s.n_diffchr[0],
        s.n_diffhigh[0],
        s.n_reads[1],
        s.n_primary[1],
        s.n_secondary[1],
        s.n_supp[1],
        s.n_dup[1],
        s.n_pdup[1],
        s.n_mapped[1],
        json_pct(s.n_mapped[1], s.n_reads[1]),
        s.n_pmapped[1],
        json_pct(s.n_pmapped[1], s.n_primary[1]),
        s.n_pair_all[1],
        s.n_read1[1],
        s.n_read2[1],
        s.n_pair_good[1],
        json_pct(s.n_pair_good[1], s.n_pair_all[1]),
        s.n_pair_map[1],
        s.n_sgltn[1],
        json_pct(s.n_sgltn[1], s.n_pair_all[1]),
        s.n_diffchr[1],
        s.n_diffhigh[1],
    );
}

fn json_pct(n: u64, total: u64) -> String {
    if total == 0 {
        "null".to_string()
    } else {
        let v = (n as f32) / (total as f32) * 100.0;
        format!("{:.2}", v)
    }
}
