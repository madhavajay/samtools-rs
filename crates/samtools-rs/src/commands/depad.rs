//! `samtools depad` — convert padded SAM alignments to unpadded SAM.
//!
//! This is a first samtools-rs slice: SAM/BAM input with `-T` padded FASTA
//! reference and SAM/BAM output. CRAM input/output remains pending.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::format::Exact;

use crate::diagnostics::{print_error, print_error_errno};
use crate::io as sam_io;

/// Entry point for `samtools depad`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut reference: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut output_format = OutputFormat::Bam;
    let mut no_pg = false;
    let mut positional = Vec::new();

    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_string_lossy();
        match s.as_ref() {
            "-T" => {
                let Some(path) = iter.next() else {
                    print_error("depad", "missing value for -T");
                    let _ = print_usage();
                    return ExitCode::from(1);
                };
                reference = Some(PathBuf::from(path));
            }
            "-o" => {
                let Some(path) = iter.next() else {
                    print_error("depad", "missing value for -o");
                    let _ = print_usage();
                    return ExitCode::from(1);
                };
                output = Some(PathBuf::from(path));
            }
            "-s" => {
                output_format = OutputFormat::Sam;
            }
            "-u" | "-1" => {
                output_format = OutputFormat::Bam;
            }
            "-O" | "--output-fmt" => {
                let Some(format) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("depad", "missing value for -O/--output-fmt");
                    let _ = print_usage();
                    return ExitCode::from(1);
                };
                output_format = OutputFormat::parse(format).unwrap_or_else(|| {
                    print_error("depad", format!("unsupported output format {}", format));
                    OutputFormat::Invalid
                });
                if output_format == OutputFormat::Invalid {
                    return ExitCode::from(1);
                }
            }
            "--no-PG" => {
                no_pg = true;
            }
            "-h" | "--help" => {
                let _ = print_usage();
                return ExitCode::SUCCESS;
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("depad", format!("unknown option {}", s));
                let _ = print_usage();
                return ExitCode::from(1);
            }
            _ => positional.push(PathBuf::from(arg)),
        }
    }

    if positional.len() != 1 {
        let _ = print_usage();
        return ExitCode::from(1);
    }

    let Some(reference) = reference.as_deref() else {
        print_error("depad", "padded reference FASTA is required with -T");
        return ExitCode::from(1);
    };

    if matches!(
        output
            .as_deref()
            .and_then(|p| p.extension())
            .and_then(|e| e.to_str())
            .map(str::to_ascii_lowercase)
            .as_deref(),
        Some("sam")
    ) {
        output_format = OutputFormat::Sam;
    } else if matches!(
        output
            .as_deref()
            .and_then(|p| p.extension())
            .and_then(|e| e.to_str())
            .map(str::to_ascii_lowercase)
            .as_deref(),
        Some("bam")
    ) {
        output_format = OutputFormat::Bam;
    }

    let input = &positional[0];
    let format = match sam_io::sam_open_format(input) {
        Ok(f) => f,
        Err(e) => {
            print_error("depad", e.to_string());
            return ExitCode::from(1);
        }
    };

    if !matches!(format.exact, Exact::Sam | Exact::Bam) {
        print_error(
            "depad",
            "only SAM and BAM input are currently supported (CRAM TODO)",
        );
        return ExitCode::from(1);
    }

    let result = run_depad(
        input,
        format.exact,
        reference,
        output.as_deref(),
        output_format,
        !no_pg,
        args,
    );
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("depad", "depad failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OutputFormat {
    Sam,
    Bam,
    Invalid,
}

impl OutputFormat {
    fn parse(value: &str) -> Option<Self> {
        let head = value.split(',').next().unwrap_or(value);
        match head.to_ascii_lowercase().as_str() {
            "sam" => Some(Self::Sam),
            "bam" => Some(Self::Bam),
            _ => None,
        }
    }
}

fn run_depad(
    input: &Path,
    input_format: Exact,
    reference: &Path,
    output: Option<&Path>,
    output_format: OutputFormat,
    add_pg: bool,
    argv: &[OsString],
) -> io::Result<()> {
    let sam = depadded_sam_bytes(input, input_format, reference, add_pg, argv)?;
    match output_format {
        OutputFormat::Sam => {
            let mut writer = sam_io::open_text_output(output)?;
            writer.write_all(&sam)?;
            sam_io::check_sam_close(&mut writer)
        }
        OutputFormat::Bam => {
            let writer: Box<dyn Write> = match output {
                Some(path) => Box::new(File::create(path)?),
                None => Box::new(io::stdout()),
            };
            htslib_rs::alignment_compat::write_bam_from_sam_reader(
                BufReader::new(Cursor::new(sam)),
                writer,
            )?;
            Ok(())
        }
        OutputFormat::Invalid => unreachable!("validated during argument parsing"),
    }
}

fn depadded_sam_bytes(
    input: &Path,
    input_format: Exact,
    reference: &Path,
    add_pg: bool,
    argv: &[OsString],
) -> io::Result<Vec<u8>> {
    let references = read_padded_references(reference)?;
    let mut output = Vec::new();
    match input_format {
        Exact::Sam => {
            let mut reader = BufReader::new(File::open(input)?);
            depad_sam_reader(&mut reader, &references, add_pg, argv, &mut output)?;
        }
        Exact::Bam => {
            let sam = htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(
                input, None,
            )?;
            let mut reader = BufReader::new(Cursor::new(sam.into_bytes()));
            depad_sam_reader(&mut reader, &references, add_pg, argv, &mut output)?;
        }
        _ => unreachable!("validated by caller"),
    }
    Ok(output)
}

fn depad_sam_reader<R, W>(
    reader: &mut R,
    references: &[PaddedReference],
    add_pg: bool,
    argv: &[OsString],
    writer: &mut W,
) -> io::Result<()>
where
    R: BufRead,
    W: Write + ?Sized,
{
    let mut header = String::new();
    let mut records = Vec::new();
    let mut line = String::new();

    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }

        if line.starts_with('@') && records.is_empty() {
            header.push_str(&depad_header_line(&line, references)?);
        } else {
            records.push(depad_record_line(&line, references)?);
        }
    }

    if add_pg {
        header = crate::pg::add_samtools_pg(&header, argv)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    }

    writer.write_all(header.as_bytes())?;
    for record in records {
        writer.write_all(record.as_bytes())?;
    }

    Ok(())
}

fn depad_header_line(line: &str, references: &[PaddedReference]) -> io::Result<String> {
    if !line.starts_with("@SQ\t") {
        return Ok(line.to_string());
    }

    let trimmed = line.trim_end_matches(['\r', '\n']);
    let mut name = None;
    let mut rest = Vec::new();
    for field in trimmed.split('\t').skip(1) {
        if let Some(sn) = field.strip_prefix("SN:") {
            name = Some(sn);
        } else if field.starts_with("LN:") || field.starts_with("M5:") {
            continue;
        } else {
            rest.push(field);
        }
    }

    let Some(name) = name else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "@SQ line is missing SN",
        ));
    };
    let reference = find_reference(references, name)?;

    let mut out = format!("@SQ\tSN:{}\tLN:{}", name, reference.unpadded_len);
    for field in rest {
        out.push('\t');
        out.push_str(field);
    }
    out.push('\n');
    Ok(out)
}

fn depad_record_line(line: &str, references: &[PaddedReference]) -> io::Result<String> {
    if line.trim().is_empty() {
        return Ok(line.to_string());
    }

    let newline = if line.ends_with('\n') { "\n" } else { "" };
    let trimmed = line.trim_end_matches(['\r', '\n']);
    let mut fields = trimmed.split('\t').map(str::to_string).collect::<Vec<_>>();
    if fields.len() < 11 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "SAM record has fewer than 11 fields",
        ));
    }

    if fields[2] == "*" || fields[3] == "0" || fields[5] == "*" {
        return Ok(line.to_string());
    }

    let pos = fields[3].parse::<usize>().map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid SAM position: {}", e),
        )
    })?;
    let reference = find_reference(references, &fields[2])?;
    fields[3] = reference.unpadded_position(pos)?.to_string();
    fields[5] = depad_cigar(&fields[5], reference, pos)?;

    let mut out = fields.join("\t");
    out.push_str(newline);
    Ok(out)
}

fn depad_cigar(cigar: &str, reference: &PaddedReference, start: usize) -> io::Result<String> {
    let mut ref_pos = start;
    let mut out = Vec::new();
    let mut pending_pads = reference.pads_before_in_run(start)?;
    let mut emitted_pad_run_base = false;

    for (len, op) in parse_cigar(cigar)? {
        match op {
            'M' | '=' | 'X' => {
                for _ in 0..len {
                    if reference.is_pad(ref_pos)? {
                        if pending_pads > 0 {
                            push_cigar_op(&mut out, pending_pads, 'P');
                            pending_pads = 0;
                        }
                        push_cigar_op(&mut out, 1, 'I');
                        emitted_pad_run_base = true;
                    } else {
                        pending_pads = 0;
                        emitted_pad_run_base = false;
                        push_cigar_op(&mut out, 1, 'M');
                    }
                    ref_pos += 1;
                }
            }
            'D' | 'N' => {
                for _ in 0..len {
                    if reference.is_pad(ref_pos)? {
                        if emitted_pad_run_base {
                            push_cigar_op(&mut out, 1, 'P');
                        } else {
                            pending_pads += 1;
                        }
                    } else {
                        pending_pads = 0;
                        emitted_pad_run_base = false;
                        push_cigar_op(&mut out, 1, 'D');
                    }
                    ref_pos += 1;
                }
            }
            'I' | 'P' => {
                pending_pads = 0;
                emitted_pad_run_base = false;
                push_cigar_op(&mut out, len, op);
            }
            'S' | 'H' => {
                push_cigar_op(&mut out, len, op);
            }
            other => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("unsupported CIGAR op {}", other),
                ));
            }
        }
    }

    if out.is_empty() {
        return Ok(String::from("*"));
    }

    Ok(out
        .into_iter()
        .map(|(len, op)| format!("{}{}", len, op))
        .collect::<String>())
}

fn parse_cigar(cigar: &str) -> io::Result<Vec<(usize, char)>> {
    let mut ops = Vec::new();
    let mut len = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            len.push(ch);
            continue;
        }

        if len.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing CIGAR length before {}", ch),
            ));
        }
        let n = len.parse::<usize>().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidData, format!("invalid CIGAR: {}", e))
        })?;
        ops.push((n, ch));
        len.clear();
    }

    if !len.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "CIGAR ended with a length and no operation",
        ));
    }

    Ok(ops)
}

fn push_cigar_op(ops: &mut Vec<(usize, char)>, len: usize, op: char) {
    if len == 0 {
        return;
    }
    if let Some((last_len, last_op)) = ops.last_mut()
        && *last_op == op
    {
        *last_len += len;
        return;
    }
    ops.push((len, op));
}

#[derive(Debug)]
struct PaddedReference {
    name: String,
    is_pad: Vec<bool>,
    unpadded_before: Vec<usize>,
    unpadded_len: usize,
}

impl PaddedReference {
    fn new(name: String, sequence: String) -> Self {
        let is_pad = sequence.bytes().map(|b| b == b'*').collect::<Vec<_>>();
        let mut unpadded_before = Vec::with_capacity(is_pad.len() + 1);
        let mut count = 0;
        unpadded_before.push(0);
        for pad in &is_pad {
            if !pad {
                count += 1;
            }
            unpadded_before.push(count);
        }

        Self {
            name,
            is_pad,
            unpadded_before,
            unpadded_len: count,
        }
    }

    fn is_pad(&self, one_based_pos: usize) -> io::Result<bool> {
        self.is_pad
            .get(one_based_pos.saturating_sub(1))
            .copied()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "reference position {} is outside padded reference {}",
                        one_based_pos, self.name
                    ),
                )
            })
    }

    fn unpadded_position(&self, one_based_pos: usize) -> io::Result<usize> {
        if one_based_pos == 0 || one_based_pos > self.is_pad.len() + 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference position {} is outside padded reference {}",
                    one_based_pos, self.name
                ),
            ));
        }
        Ok(self.unpadded_before[one_based_pos - 1] + 1)
    }

    fn pads_before_in_run(&self, one_based_pos: usize) -> io::Result<usize> {
        if one_based_pos == 0 || one_based_pos > self.is_pad.len() + 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference position {} is outside padded reference {}",
                    one_based_pos, self.name
                ),
            ));
        }

        let mut count = 0;
        let mut idx = one_based_pos.saturating_sub(1);
        while idx > 0 && self.is_pad[idx - 1] {
            count += 1;
            idx -= 1;
        }
        Ok(count)
    }
}

fn find_reference<'a>(
    references: &'a [PaddedReference],
    name: &str,
) -> io::Result<&'a PaddedReference> {
    references.iter().find(|r| r.name == name).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("reference {} was not found in padded FASTA", name),
        )
    })
}

fn read_padded_references(path: &Path) -> io::Result<Vec<PaddedReference>> {
    let mut refs = Vec::new();
    let mut name: Option<String> = None;
    let mut sequence = String::new();

    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = name.take() {
                refs.push(PaddedReference::new(name, std::mem::take(&mut sequence)));
            }
            let Some(id) = rest.split_whitespace().next() else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "empty FASTA definition line",
                ));
            };
            name = Some(id.to_string());
        } else {
            sequence.push_str(line.trim());
        }
    }

    if let Some(name) = name {
        refs.push(PaddedReference::new(name, sequence));
    }

    if refs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "padded FASTA contains no references",
        ));
    }

    Ok(refs)
}

fn print_usage() -> io::Result<()> {
    let mut w = io::stderr().lock();
    writeln!(
        w,
        "Usage: samtools depad -T REF.fa [options] <in.sam|in.bam>"
    )?;
    writeln!(w, "  -T FILE       padded reference FASTA")?;
    writeln!(w, "  -s            write SAM output")?;
    writeln!(w, "  -u, -1        write BAM output")?;
    writeln!(w, "  -o FILE       write output to FILE")?;
    writeln!(w, "  --no-PG       do not add a new @PG header line")?;
    Ok(())
}
