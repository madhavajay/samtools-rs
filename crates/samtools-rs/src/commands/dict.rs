//! `samtools dict` — create a sequence dictionary file from a FASTA.
//!
//! Mirrors `dict.c`. Reads a FASTA (optionally gzip/BGZF-compressed),
//! and for each sequence:
//!  - filters to printable ASCII characters,
//!  - uppercases them,
//!  - computes MD5,
//!  - emits `@SQ SN:<name>\tLN:<len>\tM5:<md5>` with optional `AH`, `AN`,
//!    `UR`, `AS`, `SP` tags.

use std::collections::HashSet;
use std::ffi::OsString;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use flate2::read::MultiGzDecoder;
use md5::{Digest, Md5};

use crate::diagnostics::{print_error, print_error_errno};

/// Entry point for `samtools dict`.
pub fn main(args: &[OsString]) -> ExitCode {
    let opts = match parse_args(args) {
        Ok(o) => o,
        Err(ParseError::Usage) => {
            let _ = write_usage(&mut io::stderr());
            return ExitCode::from(1);
        }
        Err(ParseError::Err(msg)) => {
            print_error("dict", msg);
            let _ = write_usage(&mut io::stderr());
            return ExitCode::from(1);
        }
    };

    let input_path = opts.input.clone().unwrap_or_else(|| PathBuf::from("-"));

    let is_alt = match opts.alt_path.as_ref() {
        Some(p) => match read_alt_names(p) {
            Ok(s) => Some(s),
            Err(e) => {
                print_error_errno("dict", format!("Cannot open {}", p.display()), &e);
                return ExitCode::from(1);
            }
        },
        None => None,
    };

    let mut out: Box<dyn Write> = match opts.output_path.as_ref() {
        Some(p) => match File::create(p) {
            Ok(f) => Box::new(f),
            Err(e) => {
                print_error_errno(
                    "dict",
                    format!("Cannot open {} for writing", p.display()),
                    &e,
                );
                return ExitCode::from(1);
            }
        },
        None => Box::new(io::stdout().lock()),
    };

    if opts.header {
        let _ = writeln!(out, "@HD\tVN:1.0\tSO:unsorted");
    }

    match write_dict(&mut out, &input_path, &opts, is_alt.as_ref()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            print_error_errno("dict", format!("Cannot open {}", input_path.display()), &e);
            ExitCode::from(1)
        }
    }
}

struct Opts {
    input: Option<PathBuf>,
    output_path: Option<PathBuf>,
    alt_path: Option<PathBuf>,
    assembly: Option<String>,
    species: Option<String>,
    uri: Option<String>,
    alias: bool,
    header: bool,
}

impl Default for Opts {
    fn default() -> Self {
        Self {
            input: None,
            output_path: None,
            alt_path: None,
            assembly: None,
            species: None,
            uri: None,
            alias: false,
            header: true,
        }
    }
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
            if opts.input.is_none() {
                opts.input = Some(PathBuf::from(&args[i]));
                i += 1;
                continue;
            }
            return Err(ParseError::Err("too many positional arguments".into()));
        };

        match s {
            "-h" | "--help" | "-?" => return Err(ParseError::Usage),
            "-H" | "--no-header" => {
                opts.header = false;
                i += 1;
            }
            "-A" | "--alias" | "--alternative-name" => {
                opts.alias = true;
                i += 1;
            }
            "-a" | "--assembly" => {
                i += 1;
                opts.assembly = Some(value(args, i, "-a")?);
                i += 1;
            }
            "-s" | "--species" => {
                i += 1;
                opts.species = Some(value(args, i, "-s")?);
                i += 1;
            }
            "-u" | "--uri" => {
                i += 1;
                opts.uri = Some(value(args, i, "-u")?);
                i += 1;
            }
            "-o" | "--output" => {
                i += 1;
                opts.output_path = Some(PathBuf::from(value(args, i, "-o")?));
                i += 1;
            }
            "-l" | "--alt" => {
                i += 1;
                opts.alt_path = Some(PathBuf::from(value(args, i, "-l")?));
                i += 1;
            }
            _ if s.starts_with("--") => {
                return Err(ParseError::Err(format!("unknown option {}", s)));
            }
            _ if s.starts_with('-')
                && s.len() > 2
                && !"alosuT".contains(s.chars().nth(1).unwrap_or(' ')) =>
            {
                // Bundled short flags like `-AH`. None of the option flags
                // here are value-taking when bundled, so split into singles.
                let mut new_args: Vec<OsString> = Vec::with_capacity(s.len() - 1);
                for c in s.chars().skip(1) {
                    new_args.push(OsString::from(format!("-{}", c)));
                }
                let mut rebuilt = args[..i].to_vec();
                rebuilt.extend(new_args);
                rebuilt.extend(args[i + 1..].iter().cloned());
                return parse_args(&rebuilt);
            }
            _ if s.starts_with('-') && s != "-" => {
                return Err(ParseError::Err(format!("unknown option {}", s)));
            }
            _ => {
                if opts.input.is_some() {
                    return Err(ParseError::Err("too many positional arguments".into()));
                }
                opts.input = Some(PathBuf::from(&args[i]));
                i += 1;
            }
        }
    }
    Ok(opts)
}

fn value(args: &[OsString], i: usize, name: &str) -> Result<String, ParseError> {
    args.get(i)
        .and_then(|a| a.to_str())
        .map(|s| s.to_string())
        .ok_or_else(|| ParseError::Err(format!("missing value for {}", name)))
}

fn read_alt_names(path: &Path) -> io::Result<HashSet<String>> {
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if is_bgzf_path(path)? {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut set = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('@') {
            continue;
        }
        let name = line.split('\t').next().unwrap_or("");
        if !name.is_empty() {
            set.insert(name.to_string());
        }
    }
    Ok(set)
}

fn is_bgzf_path(path: &Path) -> io::Result<bool> {
    let mut file = File::open(path)?;
    let mut hdr = [0u8; 2];
    let n = file.read(&mut hdr)?;
    Ok(n >= 2 && hdr[0] == 0x1f && hdr[1] == 0x8b)
}

fn write_dict<W: Write>(
    out: &mut W,
    path: &Path,
    opts: &Opts,
    is_alt: Option<&HashSet<String>>,
) -> io::Result<()> {
    let from_stdin = path.to_str() == Some("-");
    let reader: Box<dyn BufRead> = if from_stdin {
        Box::new(BufReader::new(io::stdin().lock()))
    } else {
        let file = File::open(path)?;
        if is_bgzf_path(path)? {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        }
    };

    let abs_path = if from_stdin {
        path.to_path_buf()
    } else {
        std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf())
    };

    let mut name: Option<String> = None;
    let mut seq_filtered: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(prev_name) = name.take() {
                write_sq(
                    out,
                    &prev_name,
                    &seq_filtered,
                    opts,
                    is_alt,
                    path,
                    &abs_path,
                )?;
                seq_filtered.clear();
            }
            let new_name = rest
                .split_ascii_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            name = Some(new_name);
        } else if name.is_some() {
            for &b in line.as_bytes() {
                if (b'!'..=b'~').contains(&b) {
                    seq_filtered.push(b.to_ascii_uppercase());
                }
            }
        }
    }
    if let Some(prev_name) = name.take() {
        write_sq(
            out,
            &prev_name,
            &seq_filtered,
            opts,
            is_alt,
            path,
            &abs_path,
        )?;
    }
    Ok(())
}

fn write_sq<W: Write>(
    out: &mut W,
    name: &str,
    seq: &[u8],
    opts: &Opts,
    is_alt: Option<&HashSet<String>>,
    path: &Path,
    abs_path: &Path,
) -> io::Result<()> {
    let mut hasher = Md5::new();
    hasher.update(seq);
    let digest = hasher.finalize();
    let mut hex = String::with_capacity(32);
    for b in digest.iter() {
        let _ = write!(hex, "{:02x}", b);
    }
    write!(out, "@SQ\tSN:{}\tLN:{}\tM5:{}", name, seq.len(), hex)?;

    if let Some(alt) = is_alt
        && alt.contains(name)
    {
        write!(out, "\tAH:*")?;
    }

    if opts.alias {
        // Upstream: strip the `chr` prefix if present, emit the unprefixed
        // name (or `chr<name>` if absent), and *then* append `,chrMT,MT` /
        // `,chrM,M` based on the canonicalised name (with `chr` stripped).
        let stripped = name.strip_prefix("chr").unwrap_or(name);
        if name.starts_with("chr") {
            write!(out, "\tAN:{}", stripped)?;
        } else {
            write!(out, "\tAN:chr{}", stripped)?;
        }
        if stripped == "M" {
            write!(out, ",chrMT,MT")?;
        } else if stripped == "MT" {
            write!(out, ",chrM,M")?;
        }
    }

    if let Some(uri) = opts.uri.as_ref() {
        write!(out, "\tUR:{}", uri)?;
    } else if path.to_str() != Some("-") {
        write!(out, "\tUR:file://{}", abs_path.display())?;
    }

    if let Some(a) = opts.assembly.as_ref() {
        write!(out, "\tAS:{}", a)?;
    }
    if let Some(sp) = opts.species.as_ref() {
        write!(out, "\tSP:{}", sp)?;
    }
    writeln!(out)?;
    Ok(())
}

fn write_usage<W: Write>(w: &mut W) -> io::Result<()> {
    writeln!(w)?;
    writeln!(
        w,
        "About:   Create a sequence dictionary file from a fasta file"
    )?;
    writeln!(w, "Usage:   samtools dict [options] <file.fa|file.fa.gz>")?;
    writeln!(w)?;
    writeln!(w, "Options: -a, --assembly STR    assembly")?;
    writeln!(w, "         -A, --alias, --alternative-name")?;
    writeln!(
        w,
        "                               add AN tag by adding/removing 'chr'"
    )?;
    writeln!(w, "         -H, --no-header       do not print @HD line")?;
    writeln!(
        w,
        "         -l, --alt FILE        add AH:* tag to alternate locus sequences"
    )?;
    writeln!(
        w,
        "         -o, --output FILE     file to write out dict file [stdout]"
    )?;
    writeln!(w, "         -s, --species STR     species")?;
    writeln!(
        w,
        "         -u, --uri STR         URI [file:///abs/path/to/file.fa]"
    )?;
    writeln!(w)?;
    Ok(())
}
