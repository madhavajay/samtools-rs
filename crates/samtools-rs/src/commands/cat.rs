//! `samtools cat` — concatenate alignment files of the same format.
//!
//! Mirrors `main_cat` in `bam_cat.c`. The upstream implementation
//! concatenates BAM files at the BGZF block level for speed and supports
//! `-b <fofn>` (input file list), `-h <hdr>` (replace header), `-o <out>`
//! (output), and `-p N/M` for CRAM.
//!
//! This Rust port implements record-level concatenation (decompress +
//! re-encode) for BAM and the upstream-fixtured CRAM paths.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::bam;
use htslib_rs::bgzf;
use htslib_rs::format::Exact;
use htslib_rs::sam;

use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
use crate::io as sam_io;

/// Entry point for `samtools cat`.
pub fn main(args: &[OsString]) -> ExitCode {
    let mut header_file: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut no_pg = false;
    let mut region: Option<String> = None;
    let mut part: Option<Part> = None;
    let mut input_lists: Vec<PathBuf> = Vec::new();
    let mut inputs: Vec<PathBuf> = Vec::new();
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        match s {
            "-h" => {
                let v = iter.next().map(PathBuf::from);
                header_file = v;
            }
            "-o" => {
                let v = iter.next().map(PathBuf::from);
                output = v;
            }
            "-b" => {
                let Some(v) = iter.next() else {
                    print_error("cat", "missing value for -b");
                    return ExitCode::from(1);
                };
                input_lists.push(PathBuf::from(v));
            }
            "--no-PG" => {
                no_pg = true;
            }
            "-r" => {
                region = iter.next().and_then(|a| a.to_str().map(str::to_owned));
            }
            "-p" => {
                let Some(v) = iter.next().and_then(|a| a.to_str()) else {
                    print_error("cat", "missing value for -p");
                    return ExitCode::from(1);
                };
                match parse_part(v) {
                    Ok(parsed) => part = Some(parsed),
                    Err(e) => {
                        print_error("cat", e);
                        return ExitCode::from(1);
                    }
                }
            }
            "-q" | "-f" => {
                // Reserved upstream flags not yet supported.
                print_error(
                    "cat",
                    format!("option `{}` is not yet supported in samtools-rs cat", s),
                );
                return ExitCode::from(1);
            }
            _ if s.starts_with('-') && s != "-" => {
                print_error("cat", format!("unknown option {}", s));
                return ExitCode::from(1);
            }
            _ => {
                inputs.push(PathBuf::from(arg));
            }
        }
    }

    if !input_lists.is_empty() {
        let mut expanded = Vec::new();
        for list in &input_lists {
            match read_input_list(list) {
                Ok(mut list_inputs) => expanded.append(&mut list_inputs),
                Err(e) => {
                    print_error("cat", e.to_string());
                    return ExitCode::from(1);
                }
            }
        }
        expanded.extend(inputs);
        inputs = expanded;
    }

    if inputs.is_empty() {
        let _ = writeln!(
            io::stderr(),
            "Usage: samtools cat [-b list] [-h hdr] [-o out] in1.bam ..."
        );
        return ExitCode::from(1);
    }

    // Determine input format from the first file.
    if inputs[0].as_os_str() != "-" && !inputs[0].exists() {
        print_hts_open_missing(&inputs[0]);
        print_error(
            "cat",
            format!(
                "failed to open file '{}': No such file or directory",
                inputs[0].display()
            ),
        );
        return ExitCode::from(1);
    }
    let format = match sam_io::sam_open_format(&inputs[0]) {
        Ok(f) => f,
        Err(e) => {
            print_error("cat", e.to_string());
            return ExitCode::from(1);
        }
    };

    let result = match format.exact {
        Exact::Sam => {
            let _ = writeln!(io::stderr(), "[main_cat] ERROR: input is not BAM or CRAM");
            return ExitCode::from(1);
        }
        Exact::Bam => run_bam_cat(
            &inputs,
            header_file.as_deref(),
            output.as_deref(),
            !no_pg,
            args,
            region.as_deref(),
        ),
        Exact::Cram => run_cram_cat(
            &inputs,
            header_file.as_deref(),
            output.as_deref(),
            !no_pg,
            args,
            region.as_deref(),
            part,
        ),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "unsupported input format for cat",
        )),
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("cat", "concatenation failed", &e);
            ExitCode::from(1)
        }
    }
}

#[derive(Clone, Copy)]
struct Part {
    index: usize,
    total: usize,
}

fn parse_part(raw: &str) -> Result<Part, String> {
    let Some((index, total)) = raw.split_once('/') else {
        return Err(format!("malformed region {raw}. Should be e.g. '1/10'"));
    };
    let index = index
        .parse::<usize>()
        .map_err(|_| format!("malformed region {raw}. Should be e.g. '1/10'"))?;
    let total = total
        .parse::<usize>()
        .map_err(|_| format!("malformed region {raw}. Should be e.g. '1/10'"))?;
    if index == 0 || total == 0 || index > total {
        return Err(format!("malformed region {raw}. Should be e.g. '1/10'"));
    }
    Ok(Part { index, total })
}

fn read_input_list(path: &Path) -> io::Result<Vec<PathBuf>> {
    let file = File::open(path)?;
    let mut inputs = Vec::new();

    for line in BufReader::new(file).lines() {
        let line = line?;
        let line = line.trim();
        if !line.is_empty() {
            inputs.push(PathBuf::from(line));
        }
    }

    Ok(inputs)
}

fn run_bam_cat(
    inputs: &[PathBuf],
    header_file: Option<&Path>,
    output: Option<&Path>,
    add_pg: bool,
    argv: &[OsString],
    region: Option<&str>,
) -> io::Result<()> {
    let parsed_region = match region {
        Some(r) => Some(r.parse::<htslib_rs::core::Region>().map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("invalid -r region \"{r}\": {e}"),
            )
        })?),
        None => None,
    };

    // Pick the header. If -h <hdr.sam>, parse that; otherwise use the first
    // input file's header.
    let mut header: sam::Header = match header_file {
        Some(p) => {
            let mut reader = sam::io::Reader::new(BufReader::new(File::open(p)?));
            reader.read_header()?
        }
        None => {
            let mut reader = bam::io::Reader::new(File::open(&inputs[0])?);
            reader.read_header()?
        }
    };
    if add_pg {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

    // Open output writer.
    let mut writer: Box<dyn BamSink> = match output {
        Some(p) => Box::new(FileSink::new(p)?),
        None => Box::new(StdoutSink::new()?),
    };
    writer.write_header(&header)?;

    for input in inputs {
        let mut reader = bam::io::Reader::new(File::open(input)?);
        let input_header = reader.read_header()?;
        match parsed_region.as_ref() {
            Some(region) => {
                for record in
                    htslib_rs::alignment_compat::query_bam_records_from_path(input, region)?
                {
                    writer.write_record(&input_header, &record)?;
                }
            }
            None => {
                let mut record = bam::Record::default();
                loop {
                    let n = reader.read_record(&mut record)?;
                    if n == 0 {
                        break;
                    }
                    writer.write_record(&input_header, &record)?;
                }
            }
        }
    }
    writer.finish()
}

fn run_cram_cat(
    inputs: &[PathBuf],
    header_file: Option<&Path>,
    output: Option<&Path>,
    add_pg: bool,
    argv: &[OsString],
    region: Option<&str>,
    part: Option<Part>,
) -> io::Result<()> {
    for input in inputs {
        let format = sam_io::sam_open_format(input)?;
        if format.exact != Exact::Cram {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "all cat inputs must use the same format",
            ));
        }
    }

    let mut header: sam::Header = match header_file {
        Some(p) => {
            let mut reader = sam::io::Reader::new(BufReader::new(File::open(p)?));
            reader.read_header()?
        }
        None => htslib_rs::alignment_compat::read_cram_header_from_path(&inputs[0])?,
    };
    if add_pg {
        header = crate::pg::add_samtools_pg_to_header(&header, argv)?;
    }

    let selected_inputs = inputs_for_part(inputs, part);
    let parsed_region = match region {
        Some(r) => Some(r.parse::<htslib_rs::core::Region>().map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("invalid -r region \"{r}\": {e}"),
            )
        })?),
        None => None,
    };

    // Spool the concatenated stream as SAM into a temp file, then re-encode
    // it as a real CRAM. Writing SAM text directly to a `.cram`-named output
    // (the previous behaviour) produced a file that failed format detection,
    // so `cat`-ing the result again (e.g. `cat -p` parts) was rejected.
    let mut reference: Option<PathBuf> = None;
    let (tmp_file, tmp_path) = crate::tmp_file::create_temp_file("cat", Some("sam"))?;
    {
        let mut spool: Box<dyn Write> = Box::new(tmp_file);
        crate::sam_render::write_header(&mut spool, &header)?;

        for input in selected_inputs {
            let input_reference =
                reference_from_header_uri(input, Exact::Cram)?.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("CRAM input {} requires @SQ UR tags", input.display()),
                    )
                })?;
            let text = if let Some(region) = parsed_region.as_ref() {
                htslib_rs::alignment_compat::view_cram_regions_as_sam_text_from_path_with_reference(
                    input,
                    &input_reference,
                    std::slice::from_ref(region),
                    false,
                )?
            } else {
                htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                    input,
                    &input_reference,
                    None,
                )?
            };
            write_sam_record_lines(&mut spool, &crate::sam_render::fix_sam_text(&text))?;
            if reference.is_none() {
                reference = Some(input_reference);
            }
        }
        spool.flush()?;
    }

    let reference = reference.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "cat: no CRAM inputs selected for this part",
        )
    })?;
    crate::reference::ensure_fai_index(&reference, None)?;

    match output {
        Some(path) => {
            let out = File::create(path)?;
            htslib_rs::alignment_compat::write_cram_from_sam_path_with_reference(
                tmp_path.path(),
                &reference,
                out,
            ).map(|_| ())?;
        }
        None => {
            let out = io::stdout().lock();
            htslib_rs::alignment_compat::write_cram_from_sam_path_with_reference(
                tmp_path.path(),
                &reference,
                out,
            ).map(|_| ())?;
        }
    }

    Ok(())
}

fn inputs_for_part(inputs: &[PathBuf], part: Option<Part>) -> &[PathBuf] {
    let Some(part) = part else {
        return inputs;
    };
    let start = (part.index - 1) * inputs.len() / part.total;
    let end = part.index * inputs.len() / part.total;
    &inputs[start..end]
}

fn reference_from_header_uri(input: &Path, exact: Exact) -> io::Result<Option<PathBuf>> {
    let header_text = crate::header_text::read_raw_header_text_with_format(input, exact)?;
    for line in header_text.lines().filter(|line| line.starts_with("@SQ\t")) {
        for field in line.split('\t').skip(1) {
            let Some(uri) = field.strip_prefix("UR:") else {
                continue;
            };
            let path = uri.strip_prefix("file://").unwrap_or(uri);
            let path = PathBuf::from(path);
            if path.is_file() {
                return Ok(Some(path));
            }
        }
    }
    Ok(None)
}

fn write_sam_record_lines(out: &mut dyn Write, text: &str) -> io::Result<()> {
    for line in text.lines().filter(|line| !line.starts_with('@')) {
        out.write_all(line.as_bytes())?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

trait BamSink {
    fn write_header(&mut self, header: &sam::Header) -> io::Result<()>;
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()>;
    fn finish(self: Box<Self>) -> io::Result<()>;
}

struct FileSink {
    writer: bam::io::Writer<bgzf::io::Writer<File>>,
}

impl FileSink {
    fn new(path: &Path) -> io::Result<Self> {
        Ok(Self {
            writer: bam::io::Writer::new(File::create(path)?),
        })
    }
}

impl BamSink for FileSink {
    fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.writer.write_header(header)
    }
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        self.writer.write_record(header, record)
    }
    fn finish(self: Box<Self>) -> io::Result<()> {
        Ok(())
    }
}

struct StdoutSink {
    writer: bam::io::Writer<bgzf::io::Writer<io::Stdout>>,
}

impl StdoutSink {
    fn new() -> io::Result<Self> {
        Ok(Self {
            writer: bam::io::Writer::new(io::stdout()),
        })
    }
}

impl BamSink for StdoutSink {
    fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.writer.write_header(header)
    }
    fn write_record(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        self.writer.write_record(header, record)
    }
    fn finish(self: Box<Self>) -> io::Result<()> {
        Ok(())
    }
}
