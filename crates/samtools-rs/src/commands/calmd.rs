//! `samtools calmd` / `samtools fillmd` — recalculate MD/NM tags and BAQ.
//!
//! Mirrors `bam_fillmd` in `bam_md.c`. The full upstream version recomputes
//! MD/NM tags against a reference and optionally applies BAQ (Base Alignment
//! Quality). The initial Rust port wraps the BAQ paths already available in
//! `htslib_rs::alignment_compat`:
//!  - default (`-r` not set): recomputes MD/NM tags when a reference is
//!    supplied, otherwise emits the input as SAM. SAM, BAM, and reference-
//!    backed CRAM input are supported for this record-text path.
//!  - `-r`: recalculate BAQ (`recalculate_baq_from_sam_path`).
//!  - `-E`: extended BAQ (`recalculate_extended_baq_from_sam_path`).
//!  - `-d`: drop all aux tags except `RG` from the SAM-text output.
//!  - `-e`: change matching bases to `=`.
//!  - `-q`: reduce base quality resolution.
//!
//!  - `-A` (with `-r`): apply the recalculated BAQ to the quality
//!    string (`recalculate_and_apply_baq_from_sam_path`).
//!  - `-b`/`-u`: BAM output (compressed / uncompressed), so
//!    `calmd -uAr in.sam ref.fa` emits a BGZF stream like upstream.
//!  - `-C cap`: cap MAPQ using the upstream `sam_cap_mapq` algorithm
//!    when `cap > 10`.
//!  - `-n max_nm`: when the recomputed NM is at least `max_nm`, matching
//!    bases are masked to `N` and their qualities are set to zero.
//!  - `-O cram` / `--output-fmt=cram`: CRAM output, encoded from the
//!    recalculated SAM stream with the supplied reference.
//!  - glued short-option clusters (`-uAr`) are split like `getopt`.
//!
//! **Pending:** full upstream MD/BAQ byte parity beyond the promoted harness
//! case.

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::format::Exact;

use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
use crate::io as sam_io;

/// Splits glued short-option clusters (`-uAr` → `-u -A -r`) the way
/// `getopt` does. Value-taking options (`-C`, `-n`, `-@`) consume the
/// rest of the cluster as their argument (or the next token), so
/// `-rC20` becomes `-r -C 20`. `--long`, a lone `-`, and non-option
/// operands pass through unchanged.
fn expand_short_clusters(args: &[OsString]) -> Vec<OsString> {
    const TAKES_VALUE: &[char] = &['C', 'n', '@', 'O'];
    let mut out: Vec<OsString> = Vec::with_capacity(args.len());
    for arg in args {
        let s = arg.to_str().unwrap_or("");
        let is_cluster = s.len() > 2
            && s.starts_with('-')
            && !s.starts_with("--")
            && s[1..].chars().all(|c| c.is_ascii_alphanumeric());
        if !is_cluster {
            out.push(arg.clone());
            continue;
        }
        let chars: Vec<char> = s[1..].chars().collect();
        let mut i = 0;
        while i < chars.len() {
            let c = chars[i];
            out.push(OsString::from(format!("-{c}")));
            if TAKES_VALUE.contains(&c) {
                let rest: String = chars[i + 1..].iter().collect();
                if !rest.is_empty() {
                    out.push(OsString::from(rest));
                }
                break;
            }
            i += 1;
        }
    }
    out
}

/// Entry point for `samtools calmd` / `samtools fillmd`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut realn = false;
    let mut extended = false;
    let mut use_equal = false;
    let mut always_apply = false;
    let mut bam_out = false;
    let mut uncompressed = false;
    let mut drop_aux_except_rg = false;
    let mut bin_quality = false;
    let mut update_md_nm = true;
    let mut reference: Option<PathBuf> = None;
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut no_pg = false;
    let mut quiet = false;
    let mut output_fmt: Option<OutFmt> = None;
    let mut max_nm = 0usize;
    let mut cap_q = 0i32;

    let args = expand_short_clusters(args);
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        if let Some(v) = s.strip_prefix("--output-fmt=") {
            output_fmt = match parse_output_format(v) {
                Ok(fmt) => Some(fmt),
                Err(e) => {
                    print_error("calmd", e);
                    return ExitCode::from(1);
                }
            };
            continue;
        }
        match s {
            "-r" => realn = true,
            "-E" => extended = true,
            "-e" => use_equal = true,
            "--no-PG" => no_pg = true,
            "-d" => {
                drop_aux_except_rg = true;
            }
            "-A" => always_apply = true,
            "-b" => {
                bam_out = true;
                output_fmt = Some(OutFmt::Bam);
            }
            "-u" => {
                bam_out = true;
                uncompressed = true;
                output_fmt = Some(OutFmt::Bam);
            }
            "-O" | "--output-fmt" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("calmd", "missing value for --output-fmt");
                    return ExitCode::from(1);
                };
                output_fmt = match parse_output_format(v) {
                    Ok(fmt) => Some(fmt),
                    Err(e) => {
                        print_error("calmd", e);
                        return ExitCode::from(1);
                    }
                };
                bam_out = matches!(output_fmt, Some(OutFmt::Bam));
            }
            "-n" => {
                if let Some(value) = iter.next().and_then(|a| a.to_str()) {
                    max_nm = value.parse::<usize>().unwrap_or(0);
                }
            }
            "-C" => {
                if let Some(value) = iter.next().and_then(|a| a.to_str()) {
                    cap_q = value.parse::<i32>().unwrap_or(0);
                }
            }
            "-Q" => quiet = true,
            "-q" => bin_quality = true,
            "-N" => update_md_nm = false,
            "-S" | "-h" => {
                // Accepted for CLI compatibility. The current SAM-text path
                // does not need their values beyond -C/-n, which are handled
                // above.
            }
            "-T" | "--reference" => {
                reference = iter.next().map(PathBuf::from);
            }
            "-o" | "--output" => {
                output = iter.next().map(PathBuf::from);
            }
            "-@" | "--threads" => {
                let _ = iter.next();
            }
            "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("calmd", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                } else if reference.is_none() {
                    reference = Some(PathBuf::from(arg));
                }
            }
        }
    }

    if output_fmt.is_none()
        && let Some(path) = output.as_deref()
        && path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case("cram"))
    {
        output_fmt = Some(OutFmt::Cram);
    }
    let output_fmt = output_fmt.unwrap_or(OutFmt::Sam);
    bam_out = bam_out || matches!(output_fmt, OutFmt::Bam);

    let Some(input) = input else {
        let _ = print_usage();
        return ExitCode::from(1);
    };

    if input.as_os_str() != "-" && !input.exists() {
        print_hts_open_missing(&input);
        print_error(
            "calmd",
            format!(
                "Failed to open input file '{}': No such file or directory",
                input.display()
            ),
        );
        return ExitCode::from(1);
    }

    let format = match sam_io::sam_open_format(&input) {
        Ok(f) => f,
        Err(e) => {
            print_error("calmd", e.to_string());
            return ExitCode::from(1);
        }
    };
    if !matches!(format.exact, Exact::Sam | Exact::Bam | Exact::Cram) {
        print_error("calmd", "only SAM, BAM, and CRAM input are supported");
        return ExitCode::from(1);
    }
    if format.exact == Exact::Cram && reference.is_none() {
        print_error("calmd", "CRAM input requires -T/--reference or ref.fa");
        return ExitCode::from(1);
    }
    let text = if realn {
        let Some(reference) = reference.as_ref() else {
            print_error("calmd", "-r/--reference required for BAQ recalculation");
            return ExitCode::from(1);
        };
        if format.exact == Exact::Sam {
            run_baq_from_sam_path(&input, reference, extended, always_apply)
        } else {
            input_as_sam_text(&input, format.exact, Some(reference.as_path())).and_then(|text| {
                run_baq_from_sam_text(&text, Some(reference), extended, always_apply)
            })
        }
    } else {
        input_as_sam_text(&input, format.exact, reference.as_deref())
    };

    let mut text = match text {
        Ok(t) => t,
        Err(e) => {
            print_error_errno(
                "calmd",
                format!("calmd failed for \"{}\"", input.display()),
                &e,
            );
            return ExitCode::from(1);
        }
    };

    if cap_q > 10 {
        let Some(reference) = reference.as_ref() else {
            print_error("calmd", "-C requires -T/--reference or ref.fa");
            return ExitCode::from(1);
        };
        match cap_mapping_qualities_sam_text(&text, reference, cap_q) {
            Ok(modified) => text = modified,
            Err(e) => {
                print_error_errno("calmd", "cap MAPQ", &e);
                return ExitCode::from(1);
            }
        }
    }

    if let Some(reference) = reference.as_ref() {
        match recalculate_md_nm_sam_text(&text, reference, max_nm, quiet, use_equal, update_md_nm) {
            Ok(modified) => text = modified,
            Err(e) => {
                print_error_errno("calmd", "recalculate MD/NM", &e);
                return ExitCode::from(1);
            }
        }
    }
    if drop_aux_except_rg {
        text = drop_aux_except_rg_sam_text(&text);
    }
    if bin_quality {
        text = bin_quality_sam_text(&text);
    }

    let text = if no_pg {
        text
    } else {
        match inject_pg_into_sam_text(&text, &args) {
            Ok(modified) => modified,
            Err(e) => {
                print_error_errno("calmd", "inject @PG line", &e);
                return ExitCode::from(1);
            }
        }
    };

    if bam_out {
        let mut out = match sam_io::open_text_output(output.as_deref()) {
            Ok(out) => out,
            Err(e) => {
                print_error_errno("calmd", "open -o output", &e);
                return ExitCode::from(1);
            }
        };
        use htslib_rs::bgzf::io::writer::CompressionLevel;
        let level = if uncompressed {
            CompressionLevel::try_from(0).unwrap_or_default()
        } else {
            CompressionLevel::default()
        };
        let mut reader = htslib_rs::sam::io::Reader::new(io::Cursor::new(text.into_bytes()));
        let result = htslib_rs::alignment_compat::write_bam_from_sam_reader_with_compression_level(
            &mut reader,
            &mut out,
            level,
        );
        match result {
            Ok(_) => {}
            Err(e) if e.kind() == io::ErrorKind::BrokenPipe => {}
            Err(e) => {
                print_error_errno("calmd", "write BAM output", &e);
                return ExitCode::from(1);
            }
        }
        if let Err(e) = out.flush()
            && e.kind() != io::ErrorKind::BrokenPipe
        {
            print_error_errno("calmd", "flush BAM output", &e);
            return ExitCode::from(1);
        }
        return ExitCode::SUCCESS;
    }

    if matches!(output_fmt, OutFmt::Cram) {
        let Some(reference) = reference.as_ref() else {
            print_error("calmd", "CRAM output requires -T/--reference or ref.fa");
            return ExitCode::from(1);
        };
        let out: Box<dyn Write> = match output.as_deref() {
            Some(path) => match fs::File::create(path) {
                Ok(file) => Box::new(file),
                Err(e) => {
                    print_error_errno("calmd", "open -o output", &e);
                    return ExitCode::from(1);
                }
            },
            None => Box::new(io::stdout().lock()),
        };
        if let Err(e) = write_cram_output_from_sam_text(&text, reference, out) {
            print_error_errno("calmd", "write CRAM output", &e);
            return ExitCode::from(1);
        }
        return ExitCode::SUCCESS;
    }

    let mut out = match sam_io::open_text_output(output.as_deref()) {
        Ok(out) => out,
        Err(e) => {
            print_error_errno("calmd", "open -o output", &e);
            return ExitCode::from(1);
        }
    };
    if let Err(e) = sam_io::write_all_and_close(&mut out, text.as_bytes())
        && e.kind() != io::ErrorKind::BrokenPipe
    {
        print_error_errno("calmd", "write output", &e);
        return ExitCode::from(1);
    }
    ExitCode::SUCCESS
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum OutFmt {
    Sam,
    Bam,
    Cram,
}

fn parse_output_format(raw: &str) -> Result<OutFmt, String> {
    match raw.to_ascii_lowercase().as_str() {
        "sam" => Ok(OutFmt::Sam),
        "bam" => Ok(OutFmt::Bam),
        "cram" => Ok(OutFmt::Cram),
        _ => Err(format!("unsupported output format \"{}\"", raw)),
    }
}

fn run_baq_from_sam_path(
    input: &Path,
    reference: &Path,
    extended: bool,
    apply: bool,
) -> io::Result<String> {
    if extended {
        htslib_rs::alignment_compat::recalculate_extended_baq_from_sam_path(input, reference)
    } else if apply {
        htslib_rs::alignment_compat::recalculate_and_apply_baq_from_sam_path(input, reference)
    } else {
        htslib_rs::alignment_compat::recalculate_baq_from_sam_path(input, reference)
    }
}

fn run_baq_from_sam_text(
    text: &str,
    reference: Option<&Path>,
    extended: bool,
    apply: bool,
) -> io::Result<String> {
    let (mut tmp_sam, tmp_sam_path) = crate::tmp_file::create_temp_file("calmd-baq", Some("sam"))?;
    tmp_sam.write_all(text.as_bytes())?;
    tmp_sam.flush()?;
    drop(tmp_sam);

    let result = if let Some(reference) = reference {
        run_baq_from_sam_path(tmp_sam_path.path(), reference, extended, apply)
    } else {
        htslib_rs::alignment_compat::apply_existing_baq_from_sam_path(tmp_sam_path.path())
    };
    tmp_sam_path.close().ok();
    result
}

fn write_cram_output_from_sam_text<W>(text: &str, reference: &Path, out: W) -> io::Result<()>
where
    W: Write,
{
    crate::reference::ensure_fai_index(reference, None)?;
    let (mut tmp_sam, tmp_sam_path) = crate::tmp_file::create_temp_file("calmd", Some("sam"))?;
    tmp_sam.write_all(text.as_bytes())?;
    tmp_sam.flush()?;
    drop(tmp_sam);

    let result = htslib_rs::alignment_compat::write_cram_from_sam_path_with_reference(
        tmp_sam_path.path(),
        reference,
        out,
    )
    .map(|_| ());

    tmp_sam_path.close().ok();
    result
}

fn cap_mapping_qualities_sam_text(
    text: &str,
    reference: &Path,
    threshold: i32,
) -> io::Result<String> {
    let (mut tmp_sam, tmp_sam_path) = crate::tmp_file::create_temp_file("calmd-cap", Some("sam"))?;
    tmp_sam.write_all(text.as_bytes())?;
    tmp_sam.flush()?;
    drop(tmp_sam);

    let caps = htslib_rs::alignment_compat::cap_mapping_qualities_from_sam_path(
        tmp_sam_path.path(),
        reference,
        threshold,
    )?;
    let mut caps = caps.into_iter();
    let mut out = String::with_capacity(text.len());

    for line in text.split_inclusive('\n') {
        let (line_body, newline) = line
            .strip_suffix('\n')
            .map(|body| (body, "\n"))
            .unwrap_or((line, ""));

        if line_body.starts_with('@') || line_body.is_empty() {
            out.push_str(line_body);
            out.push_str(newline);
            continue;
        }

        let Some(cap) = caps.next() else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "missing MAPQ cap for SAM record",
            ));
        };

        write_record_with_capped_mapq(&mut out, line_body, cap)?;
        out.push_str(newline);
    }

    if caps.next().is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "extra MAPQ caps for SAM text",
        ));
    }

    tmp_sam_path.close().ok();
    Ok(out)
}

fn recalculate_md_nm_sam_text(
    text: &str,
    reference: &Path,
    max_nm: usize,
    quiet: bool,
    use_equal: bool,
    update_md_nm: bool,
) -> io::Result<String> {
    let references = read_fasta(reference)?;
    let mut out = String::with_capacity(text.len());

    for line in text.split_inclusive('\n') {
        let (line_body, newline) = line
            .strip_suffix('\n')
            .map(|body| (body, "\n"))
            .unwrap_or((line, ""));

        if line_body.starts_with('@') || line_body.is_empty() {
            out.push_str(line_body);
            out.push_str(newline);
            continue;
        }

        let fields: Vec<&str> = line_body.split('\t').collect();
        if fields.len() < 11 || fields[2] == "*" || fields[5] == "*" {
            out.push_str(line_body);
            out.push_str(newline);
            continue;
        }
        let flag = fields[1].parse::<u16>().map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid FLAG {}: {e}", fields[1]),
            )
        })?;
        if flag & 0x4 != 0 {
            if max_nm > 0 || use_equal {
                let Some(reference_sequence) = references.get(fields[2]) else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("reference sequence {} not found", fields[2]),
                    ));
                };
                let start = fields[3].parse::<usize>().map_err(|e| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid POS {}: {e}", fields[3]),
                    )
                })?;
                if start > 0 {
                    let (_, nm, matching_read_positions) =
                        calculate_md_nm(fields[5], fields[9], reference_sequence, start - 1, true)?;
                    let (sequence_override, quality_override) = apply_sequence_quality_overrides(
                        fields[9],
                        fields[10],
                        max_nm,
                        nm,
                        &matching_read_positions,
                        use_equal,
                    )?;
                    write_record_with_sequence_quality(
                        &mut out,
                        &fields,
                        sequence_override.as_deref(),
                        quality_override.as_deref(),
                    );
                    out.push_str(newline);
                    continue;
                }
            }
            out.push_str(line_body);
            out.push_str(newline);
            continue;
        }

        let Some(reference_sequence) = references.get(fields[2]) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("reference sequence {} not found", fields[2]),
            ));
        };

        let start = fields[3].parse::<usize>().map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid POS {}: {e}", fields[3]),
            )
        })?;
        if start == 0 {
            out.push_str(line_body);
            out.push_str(newline);
            continue;
        }

        let (md, nm, matching_read_positions) = calculate_md_nm(
            fields[5],
            fields[9],
            reference_sequence,
            start - 1,
            max_nm > 0 || use_equal,
        )?;
        let (masked_sequence, masked_quality) = apply_sequence_quality_overrides(
            fields[9],
            fields[10],
            max_nm,
            nm,
            &matching_read_positions,
            use_equal,
        )?;
        write_record_with_md_nm(
            &mut out,
            &fields,
            &md,
            nm,
            masked_sequence.as_deref(),
            masked_quality.as_deref(),
            MdNmWriteOptions {
                quiet,
                update_md_nm,
            },
        );
        out.push_str(newline);
    }

    Ok(out)
}

fn input_as_sam_text(input: &Path, exact: Exact, reference: Option<&Path>) -> io::Result<String> {
    match exact {
        Exact::Sam => htslib_rs::alignment_compat::view_sam_text_from_path_with_limit(input, None),
        Exact::Bam => {
            htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(input, None)
        }
        Exact::Cram => {
            let reference = reference.ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "CRAM input requires reference")
            })?;
            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                input, reference, None,
            )
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "unsupported alignment format",
        )),
    }
}

fn drop_aux_except_rg_sam_text(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        let (line_body, newline) = line
            .strip_suffix('\n')
            .map(|body| (body, "\n"))
            .unwrap_or((line, ""));
        if line_body.starts_with('@') || line_body.is_empty() {
            out.push_str(line_body);
            out.push_str(newline);
            continue;
        }

        let mut first = true;
        for (i, field) in line_body.split('\t').enumerate() {
            if i >= 11 && !field.starts_with("RG:Z:") {
                continue;
            }
            if !first {
                out.push('\t');
            }
            first = false;
            out.push_str(field);
        }
        out.push_str(newline);
    }
    out
}

fn bin_quality_sam_text(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        let (line_body, newline) = line
            .strip_suffix('\n')
            .map(|body| (body, "\n"))
            .unwrap_or((line, ""));
        if line_body.starts_with('@') || line_body.is_empty() {
            out.push_str(line_body);
            out.push_str(newline);
            continue;
        }

        for (i, field) in line_body.split('\t').enumerate() {
            if i > 0 {
                out.push('\t');
            }
            if i == 10 && field != "*" {
                for quality in field.bytes() {
                    let score = quality.saturating_sub(b'!');
                    let binned = if score >= 3 {
                        score / 10 * 10 + 7
                    } else {
                        score
                    };
                    out.push(char::from(binned + b'!'));
                }
            } else {
                out.push_str(field);
            }
        }
        out.push_str(newline);
    }
    out
}

fn read_fasta(path: &Path) -> io::Result<HashMap<String, Vec<u8>>> {
    let text = fs::read_to_string(path)?;
    let mut references = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = current_name.replace(
                rest.split_whitespace()
                    .next()
                    .unwrap_or_default()
                    .to_string(),
            ) {
                references.insert(name, std::mem::take(&mut current_seq));
            }
        } else if !line.starts_with(';') {
            current_seq.extend(line.bytes().filter(|b| !b.is_ascii_whitespace()).map(|b| {
                if b.is_ascii_lowercase() {
                    b.to_ascii_uppercase()
                } else {
                    b
                }
            }));
        }
    }

    if let Some(name) = current_name {
        references.insert(name, current_seq);
    }

    Ok(references)
}

fn calculate_md_nm(
    cigar: &str,
    sequence: &str,
    reference: &[u8],
    start: usize,
    collect_matching_positions: bool,
) -> io::Result<(String, usize, Vec<usize>)> {
    let read = sequence.as_bytes();
    let mut read_i = 0usize;
    let mut ref_i = start;
    let mut nm = 0usize;
    let mut md = String::new();
    let mut match_count = 0usize;
    let mut matching_read_positions = Vec::new();
    let mut n = 0usize;

    for b in cigar.bytes() {
        if b.is_ascii_digit() {
            n = n * 10 + usize::from(b - b'0');
            continue;
        }
        if n == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid CIGAR operation in {cigar}"),
            ));
        }

        match b {
            b'M' | b'=' | b'X' => {
                for _ in 0..n {
                    let read_base = read.get(read_i).copied().ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "CIGAR consumes past read")
                    })?;
                    let ref_base = reference.get(ref_i).copied().ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "CIGAR consumes past reference")
                    })?;
                    if bases_match(read_base, ref_base) {
                        match_count += 1;
                        if collect_matching_positions {
                            matching_read_positions.push(read_i);
                        }
                    } else {
                        md.push_str(&match_count.to_string());
                        match_count = 0;
                        md.push(char::from(ref_base.to_ascii_uppercase()));
                        nm += 1;
                    }
                    read_i += 1;
                    ref_i += 1;
                }
            }
            b'I' => {
                read_i += n;
                nm += n;
                if read_i > read.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CIGAR consumes past read",
                    ));
                }
            }
            b'D' => {
                md.push_str(&match_count.to_string());
                match_count = 0;
                md.push('^');
                for _ in 0..n {
                    let ref_base = reference.get(ref_i).copied().ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "CIGAR consumes past reference")
                    })?;
                    md.push(char::from(ref_base.to_ascii_uppercase()));
                    ref_i += 1;
                }
                nm += n;
            }
            b'N' => {
                ref_i += n;
                if ref_i > reference.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CIGAR consumes past reference",
                    ));
                }
            }
            b'S' => {
                read_i += n;
                if read_i > read.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CIGAR consumes past read",
                    ));
                }
            }
            b'H' | b'P' => {}
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("unsupported CIGAR operation {}", char::from(b)),
                ));
            }
        }
        n = 0;
    }

    if n != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("trailing CIGAR length in {cigar}"),
        ));
    }

    md.push_str(&match_count.to_string());
    Ok((md, nm, matching_read_positions))
}

fn bases_match(read_base: u8, ref_base: u8) -> bool {
    read_base == b'=' || (read_base != b'N' && read_base.eq_ignore_ascii_case(&ref_base))
}

fn write_record_with_capped_mapq(out: &mut String, line: &str, cap: i32) -> io::Result<()> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 11 {
        out.push_str(line);
        return Ok(());
    }

    let current_mapq = fields[4].parse::<i32>().map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid MAPQ {}: {e}", fields[4]),
        )
    })?;
    let capped_mapq = if current_mapq > cap {
        if cap < 0 { 255 } else { cap.min(255) }
    } else {
        current_mapq
    };

    for (i, field) in fields.iter().enumerate() {
        if i > 0 {
            out.push('\t');
        }
        if i == 4 {
            out.push_str(&capped_mapq.to_string());
        } else {
            out.push_str(field);
        }
    }

    Ok(())
}

fn apply_sequence_quality_overrides(
    sequence: &str,
    quality: &str,
    max_nm: usize,
    nm: usize,
    matching_read_positions: &[usize],
    use_equal: bool,
) -> io::Result<(Option<String>, Option<String>)> {
    let mask_matches = max_nm > 0 && nm >= max_nm;
    if !mask_matches && !use_equal {
        return Ok((None, None));
    }

    let mut sequence = sequence.as_bytes().to_vec();
    for &position in matching_read_positions {
        if position >= sequence.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CIGAR consumes past read",
            ));
        }
        sequence[position] = if mask_matches { b'N' } else { b'=' };
    }

    let quality = if quality == "*" || !mask_matches {
        None
    } else {
        let mut quality = quality.as_bytes().to_vec();
        for &position in matching_read_positions {
            if position < quality.len() {
                quality[position] = b'!';
            }
        }
        Some(String::from_utf8(quality).map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid SAM quality string: {e}"),
            )
        })?)
    };

    let sequence = String::from_utf8(sequence).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid SAM sequence string: {e}"),
        )
    })?;

    Ok((Some(sequence), quality))
}

fn write_record_with_sequence_quality(
    out: &mut String,
    fields: &[&str],
    sequence_override: Option<&str>,
    quality_override: Option<&str>,
) {
    for (i, field) in fields.iter().enumerate() {
        if i > 0 {
            out.push('\t');
        }
        match i {
            9 => out.push_str(sequence_override.unwrap_or(field)),
            10 => out.push_str(quality_override.unwrap_or(field)),
            _ => out.push_str(field),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct MdNmWriteOptions {
    quiet: bool,
    update_md_nm: bool,
}

fn write_record_with_md_nm(
    out: &mut String,
    fields: &[&str],
    md: &str,
    nm: usize,
    sequence_override: Option<&str>,
    quality_override: Option<&str>,
    options: MdNmWriteOptions,
) {
    if !options.update_md_nm {
        write_record_with_sequence_quality(out, fields, sequence_override, quality_override);
        return;
    }

    let mut append_nm = true;
    let mut append_md = true;
    let mut rendered_fields = Vec::with_capacity(fields.len() + 2);

    for (i, field) in fields.iter().enumerate() {
        if i >= 11 {
            if let Some(old_md) = field.strip_prefix("MD:Z:") {
                if old_md.eq_ignore_ascii_case(md) {
                    append_md = false;
                    rendered_fields.push((*field).to_string());
                }
                continue;
            }
            if let Some(old_nm) = field.strip_prefix("NM:i:") {
                if old_nm.parse::<usize>().is_ok_and(|old_nm| old_nm == nm) {
                    append_nm = false;
                    rendered_fields.push((*field).to_string());
                } else if !options.quiet {
                    eprintln!(
                        "[bam_fillmd1] different NM for read '{}': {} -> {}",
                        fields[0], old_nm, nm
                    );
                }
                continue;
            }
        }
        let value = match i {
            9 => sequence_override.unwrap_or(field),
            10 => quality_override.unwrap_or(field),
            _ => field,
        };
        rendered_fields.push(value.to_string());
    }

    if append_nm {
        rendered_fields.push(format!("NM:i:{nm}"));
    }
    if append_md {
        if !options.quiet
            && let Some(old_md) = fields
                .iter()
                .skip(11)
                .find_map(|field| field.strip_prefix("MD:Z:"))
        {
            eprintln!(
                "[bam_fillmd1] different MD for read '{}': '{}' -> '{}'",
                fields[0], old_md, md
            );
        }
        rendered_fields.push(format!("MD:Z:{md}"));
    }

    out.push_str(&rendered_fields.join("\t"));
}

/// Inserts samtools' `@PG` chain entry into a SAM text blob. Splits the
/// header (lines starting with `@`) from the body, applies the shared
/// `pg::add_samtools_pg` helper to the header, then rejoins.
fn inject_pg_into_sam_text(text: &str, argv: &[OsString]) -> io::Result<String> {
    let mut header_end = 0;
    for line in text.split_inclusive('\n') {
        if line.starts_with('@') {
            header_end += line.len();
        } else {
            break;
        }
    }
    let header_slice = &text[..header_end];
    let body_slice = &text[header_end..];
    let new_header = crate::pg::add_samtools_pg(header_slice, argv)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut combined = String::with_capacity(new_header.len() + body_slice.len());
    combined.push_str(&new_header);
    combined.push_str(body_slice);
    Ok(combined)
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(w, "Usage: samtools calmd [options] <in.sam> [ref.fa]")?;
    writeln!(w, "  -r          recalculate BAQ (requires --reference)")?;
    writeln!(w, "  -E          extended BAQ (with -r)")?;
    writeln!(w, "  -e          change identical bases to '='")?;
    writeln!(w, "  -q          bin base qualities")?;
    writeln!(w, "  -T FILE     reference FASTA")?;
    writeln!(w, "  -o FILE     output FILE")?;
    Ok(())
}
