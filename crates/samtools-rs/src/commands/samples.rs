//! `samtools samples` — list samples (`@RG SM:`) in SAM/BAM/CRAM files.
//!
//! Mirrors `main_samples` in `bam_samples.c`. The basic mode walks each
//! input file's header, dedups `@RG SM:` (or arbitrary `-T TAG`) values,
//! and prints `<sample>\t<filename>` per unique sample. Files with no
//! matching `@RG` line emit `.\t<filename>`.
//!
use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use crate::diagnostics::{print_error, print_error_errno, print_hts_open_missing};
use crate::header_text::read_raw_header_text;
use crate::io as sam_io;
use crate::reference::matching_reference;

/// Entry point for `samtools samples`.
pub fn main(args: &[OsString]) -> ExitCode {
    let opts = match parse_args(args) {
        Ok(o) => o,
        Err(ParseError::Usage) => {
            let _ = write_usage(&mut io::stdout());
            return ExitCode::SUCCESS;
        }
        Err(ParseError::Err(msg)) => {
            print_error("samples", msg);
            let _ = write_usage(&mut io::stderr());
            return ExitCode::from(1);
        }
    };

    let mut opts = opts;
    if let Err(e) = expand_fasta_lists(&mut opts) {
        print_error_errno("samples", "reading -F FASTA list", &e);
        return ExitCode::from(1);
    }

    let inputs = match collect_input_specs(&opts, io::stdin().lock()) {
        Ok(inputs) => inputs,
        Err(e) => {
            print_error("samples", e);
            return ExitCode::from(1);
        }
    };

    let mut out = match sam_io::open_text_output(opts.output.as_deref()) {
        Ok(out) => out,
        Err(e) => {
            print_error_errno("samples", "open output for writing", &e);
            return ExitCode::from(1);
        }
    };

    if opts.print_header {
        let _ = write!(out, "#{}\tPATH", opts.tag);
        if opts.test_index {
            let _ = write!(out, "\tINDEX");
        }
        if !opts.fa_paths.is_empty() {
            let _ = write!(out, "\tREFERENCE");
        }
        let _ = writeln!(out);
    }

    let mut overall = ExitCode::SUCCESS;
    for input in &inputs {
        if input.path.as_os_str() != "-" && !input.path.exists() {
            print_hts_open_missing(&input.path);
            print_error(
                "samples",
                format!(
                    "Failed to open \"{}\" for reading: No such file or directory",
                    input.path.display()
                ),
            );
            overall = ExitCode::from(1);
            continue;
        }
        if let Err(e) = print_samples_for(&mut out, input, &opts) {
            print_error_errno(
                "samples",
                format!("failed to process \"{}\"", input.path.display()),
                &e,
            );
            overall = ExitCode::from(1);
        }
    }
    match sam_io::check_sam_close(&mut out) {
        Ok(()) => overall,
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("samples", "close output", &e);
            ExitCode::from(1)
        }
    }
}

struct Opts {
    tag: String,
    print_header: bool,
    output: Option<PathBuf>,
    inputs: Vec<PathBuf>,
    /// `-i` — print an extra column showing whether the input has an
    /// associated index (`Y` or `N`).
    test_index: bool,
    /// `-f FILE` (repeatable) — FASTA paths to match `@SQ` dictionaries against.
    fa_paths: Vec<PathBuf>,
    /// `-F FILE` (repeatable) — files containing FASTA paths, one per line.
    fa_list_paths: Vec<PathBuf>,
    /// `-X` — inputs are paired with explicit index paths.
    custom_index: bool,
}

impl Default for Opts {
    fn default() -> Self {
        Self {
            tag: "SM".to_string(),
            print_header: false,
            output: None,
            inputs: Vec::new(),
            test_index: false,
            fa_paths: Vec::new(),
            fa_list_paths: Vec::new(),
            custom_index: false,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct InputSpec {
    path: PathBuf,
    index_path: Option<PathBuf>,
}

enum ParseError {
    Usage,
    Err(String),
}

fn parse_args(args: &[OsString]) -> Result<Opts, ParseError> {
    let mut opts = Opts::default();
    let mut i = 1;
    while i < args.len() {
        let Some(s) = args[i].to_str() else {
            opts.inputs.push(PathBuf::from(&args[i]));
            i += 1;
            continue;
        };
        match s {
            "-h" => {
                opts.print_header = true;
                i += 1;
            }
            "-?" | "--help" => return Err(ParseError::Usage),
            "-o" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| ParseError::Err("missing value for -o".into()))?;
                opts.output = Some(PathBuf::from(v));
                i += 1;
            }
            "-T" => {
                i += 1;
                let v = args
                    .get(i)
                    .and_then(|a| a.to_str())
                    .ok_or_else(|| ParseError::Err("missing value for -T".into()))?;
                if v.len() != 2 {
                    return Err(ParseError::Err(format!(
                        "Length of tag \"{}\" is not 2.",
                        v
                    )));
                }
                opts.tag = v.to_string();
                i += 1;
            }
            "-F" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| ParseError::Err("missing value for -F".into()))?;
                opts.fa_list_paths.push(PathBuf::from(v));
                i += 1;
            }
            "-X" => {
                opts.custom_index = true;
                i += 1;
            }
            "-i" => {
                opts.test_index = true;
                i += 1;
            }
            "-f" => {
                i += 1;
                let v = args
                    .get(i)
                    .map(PathBuf::from)
                    .ok_or_else(|| ParseError::Err("missing value for -f".into()))?;
                opts.fa_paths.push(v);
                i += 1;
            }
            _ if s.starts_with('-') && s != "-" => {
                return Err(ParseError::Err(format!("unknown option {}", s)));
            }
            _ => {
                opts.inputs.push(PathBuf::from(&args[i]));
                i += 1;
            }
        }
    }
    Ok(opts)
}

fn collect_input_specs<R: BufRead>(opts: &Opts, stdin: R) -> Result<Vec<InputSpec>, String> {
    if opts.inputs.is_empty() {
        return collect_input_specs_from_stdin(opts.custom_index, stdin);
    }

    if !opts.custom_index {
        return Ok(opts
            .inputs
            .iter()
            .map(|p| InputSpec {
                path: p.clone(),
                index_path: None,
            })
            .collect());
    }

    if !opts.inputs.len().is_multiple_of(2) {
        return Err(
            "Odd number of filenames detected! Each BAM file should have an index file".to_string(),
        );
    }

    let n = opts.inputs.len() / 2;
    Ok((0..n)
        .map(|i| InputSpec {
            path: opts.inputs[i].clone(),
            index_path: Some(opts.inputs[i + n].clone()),
        })
        .collect())
}

fn collect_input_specs_from_stdin<R: BufRead>(
    custom_index: bool,
    stdin: R,
) -> Result<Vec<InputSpec>, String> {
    let mut inputs = Vec::new();
    for line in stdin.lines() {
        let line = line.map_err(|e| format!("Cannot read from stdin: {}", e))?;
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        if custom_index {
            let Some((path, index_path)) = line.split_once('\t') else {
                return Err(format!(
                    "Expected path-to-bam(tab)path-to-index but got \"{}\"",
                    line
                ));
            };
            if index_path.is_empty() {
                return Err(format!(
                    "Expected path-to-bam(tab)path-to-index but got \"{}\"",
                    line
                ));
            }
            inputs.push(InputSpec {
                path: PathBuf::from(path),
                index_path: Some(PathBuf::from(index_path)),
            });
        } else {
            inputs.push(InputSpec {
                path: PathBuf::from(line),
                index_path: None,
            });
        }
    }
    Ok(inputs)
}

fn expand_fasta_lists(opts: &mut Opts) -> io::Result<()> {
    for p in &opts.fa_list_paths {
        let f = File::open(p)?;
        for line in BufReader::new(f).lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() {
                opts.fa_paths.push(PathBuf::from(trimmed));
            }
        }
    }
    Ok(())
}

fn print_samples_for<W: Write>(out: &mut W, input: &InputSpec, opts: &Opts) -> io::Result<()> {
    let fname = &input.path;
    let header_text = read_raw_header_text(fname)?;
    let mut seen: HashSet<String> = HashSet::new();
    let needle = format!("\t{}:", opts.tag);
    let mut sq_dict: Vec<(String, u64)> = Vec::new();
    for line in header_text.lines() {
        if line.starts_with("@SQ\t") {
            sq_dict.push(parse_sq(line));
        }
        if !line.starts_with("@RG") {
            continue;
        }
        if let Some(start) = line.find(&needle) {
            let after = &line[start + needle.len()..];
            let end = after.find('\t').unwrap_or(after.len());
            let value = &after[..end];
            seen.insert(value.to_string());
        }
    }

    let index_suffix = if opts.test_index {
        let has_index = index_present(fname, input.index_path.as_deref());
        if has_index { "\tY" } else { "\tN" }
    } else {
        ""
    };

    let ref_suffix = if opts.fa_paths.is_empty() {
        String::new()
    } else {
        let m = matching_reference(&opts.fa_paths, &sq_dict)
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| ".".to_string());
        format!("\t{}", m)
    };
    let suffix = format!("{}{}", index_suffix, ref_suffix);
    let index_suffix = suffix.as_str();

    if seen.is_empty() {
        writeln!(out, ".\t{}{}", fname.display(), index_suffix)?;
    } else {
        let mut samples: Vec<_> = seen.into_iter().collect();
        samples.sort();
        for s in samples {
            writeln!(out, "{}\t{}{}", s, fname.display(), index_suffix)?;
        }
    }
    Ok(())
}

fn parse_sq(line: &str) -> (String, u64) {
    let mut sn = String::new();
    let mut ln: u64 = 0;
    for field in line.split('\t').skip(1) {
        if let Some(v) = field.strip_prefix("SN:") {
            sn = v.to_string();
        } else if let Some(v) = field.strip_prefix("LN:")
            && let Ok(n) = v.parse()
        {
            ln = n;
        }
    }
    (sn, ln)
}

fn index_present(fname: &Path, custom_index: Option<&Path>) -> bool {
    use htslib_rs::index_compat::{IndexFormat, locate_associated_index};
    let candidates = [
        IndexFormat::Bai,
        IndexFormat::Csi,
        IndexFormat::Crai,
        IndexFormat::Tbi,
    ];
    let resolves = |base: &Path| {
        candidates
            .into_iter()
            .any(|fmt| locate_associated_index(base, fmt).is_some())
    };

    if let Some(index_path) = custom_index {
        // Upstream passes the custom path to `sam_index_load3`, which
        // accepts an exact index file, a directory holding the index, or
        // a prefix to which the standard suffix is appended. A bare
        // `.exists()` only caught the first; emulate the other two via
        // the shared resolver so non-default index locations register.
        if index_path.is_file() {
            return true;
        }
        let base = if index_path.is_dir() {
            match fname.file_name() {
                Some(name) => index_path.join(name),
                None => return false,
            }
        } else {
            index_path.to_path_buf()
        };
        return resolves(&base);
    }

    resolves(fname)
}

fn write_usage<W: Write>(w: &mut W) -> io::Result<()> {
    writeln!(
        w,
        "Usage: samtools samples [options] <in.sam>|<in.bam>|<in.cram> [...]"
    )?;
    writeln!(
        w,
        "       samtools samples [options] -X f1.bam f2.bam f1.bam.bai f2.bai"
    )?;
    writeln!(w, "Options:")?;
    writeln!(w, "  -h            print a header line")?;
    writeln!(w, "  -T TAG        @RG tag to extract (default: SM)")?;
    writeln!(w, "  -o FILE       write output to FILE")?;
    writeln!(
        w,
        "  -i            add a column showing whether the file is indexed"
    )?;
    writeln!(
        w,
        "  -f FILE.fa    add an indexed FASTA reference to match @SQ lines"
    )?;
    writeln!(
        w,
        "  -F FILE       read FASTA reference paths from FILE (one per line)"
    )?;
    writeln!(w, "  -X            use custom index file paths")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn custom_index_argv_pairs_first_half_with_second_half() {
        let opts = Opts {
            custom_index: true,
            inputs: vec![
                PathBuf::from("a.bam"),
                PathBuf::from("b.bam"),
                PathBuf::from("a.bai"),
                PathBuf::from("b.bai"),
            ],
            ..Opts::default()
        };

        let inputs = collect_input_specs(&opts, Cursor::new(Vec::<u8>::new())).unwrap();
        assert_eq!(
            inputs,
            vec![
                InputSpec {
                    path: PathBuf::from("a.bam"),
                    index_path: Some(PathBuf::from("a.bai")),
                },
                InputSpec {
                    path: PathBuf::from("b.bam"),
                    index_path: Some(PathBuf::from("b.bai")),
                },
            ]
        );
    }

    #[test]
    fn custom_index_stdin_requires_tab_pair() {
        let err = collect_input_specs_from_stdin(true, Cursor::new("a.bam\n")).unwrap_err();
        assert!(err.contains("Expected path-to-bam(tab)path-to-index"));
    }
}
