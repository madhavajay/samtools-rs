//! Shared samtools global option parsing.
//!
//! This mirrors the common long options listed in upstream `sam_opts.h`.
//! The parser is intentionally limited to top-level options that appear before
//! the subcommand so command-local parsers still own their own option surface.

use std::cell::RefCell;
use std::ffi::OsString;
use std::path::PathBuf;

/// Parsed values for upstream's standard `sam_global_args` option family.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SamGlobalArgs {
    pub input_fmt: Option<String>,
    pub input_fmt_options: Vec<String>,
    pub output_fmt: Option<String>,
    pub output_fmt_options: Vec<String>,
    pub reference: Option<PathBuf>,
    pub threads: Option<usize>,
    pub write_index: bool,
    pub verbosity: Option<i32>,
}

thread_local! {
    static CURRENT_GLOBAL_ARGS: RefCell<SamGlobalArgs> = RefCell::new(SamGlobalArgs::default());
}

/// Stores globals parsed by the top-level dispatcher for command I/O paths.
pub fn set_current_global_args(globals: SamGlobalArgs) {
    CURRENT_GLOBAL_ARGS.with(|cell| *cell.borrow_mut() = globals);
}

/// Returns the globals parsed by the top-level dispatcher.
pub fn current_global_args() -> SamGlobalArgs {
    CURRENT_GLOBAL_ARGS.with(|cell| cell.borrow().clone())
}

/// Applies top-level samtools global options and returns argv without them.
pub fn apply_top_level_global_args(
    args: Vec<OsString>,
) -> Result<(Vec<OsString>, SamGlobalArgs), String> {
    let (filtered, globals) = parse_top_level_global_args(args)?;

    if let Some(value) = globals.verbosity {
        htslib_rs::log_compat::set_hts_verbose(value);
    }

    Ok((filtered, globals))
}

/// Parses and strips recognized global options before the subcommand.
pub fn parse_top_level_global_args(
    args: Vec<OsString>,
) -> Result<(Vec<OsString>, SamGlobalArgs), String> {
    let mut out = Vec::with_capacity(args.len());
    let mut globals = SamGlobalArgs::default();
    let mut iter = args.into_iter().peekable();

    if let Some(program) = iter.next() {
        out.push(program);
    }

    let mut before_subcommand = true;
    while let Some(arg) = iter.next() {
        if !before_subcommand {
            out.push(arg);
            continue;
        }

        let Some(s) = arg.to_str() else {
            out.push(arg);
            before_subcommand = false;
            continue;
        };

        if s == "--" {
            out.push(arg);
            before_subcommand = false;
            continue;
        }

        if let Some((name, value)) = s.strip_prefix("--").and_then(split_long_option) {
            if parse_global_option(name, Some(value), &mut iter, &mut globals)? {
                continue;
            }
        } else if let Some(name) = s.strip_prefix("--")
            && parse_global_option(name, None, &mut iter, &mut globals)?
        {
            continue;
        }

        out.push(arg);
        before_subcommand = false;
    }

    Ok((out, globals))
}

fn split_long_option(s: &str) -> Option<(&str, &str)> {
    s.split_once('=')
}

fn parse_global_option<I>(
    name: &str,
    inline_value: Option<&str>,
    iter: &mut I,
    globals: &mut SamGlobalArgs,
) -> Result<bool, String>
where
    I: Iterator<Item = OsString>,
{
    match name {
        "input-fmt" => {
            globals.input_fmt = Some(required_value(name, inline_value, iter)?);
            Ok(true)
        }
        "input-fmt-option" => {
            globals
                .input_fmt_options
                .push(required_value(name, inline_value, iter)?);
            Ok(true)
        }
        "output-fmt" => {
            globals.output_fmt = Some(required_value(name, inline_value, iter)?);
            Ok(true)
        }
        "output-fmt-option" => {
            globals
                .output_fmt_options
                .push(required_value(name, inline_value, iter)?);
            Ok(true)
        }
        "reference" => {
            globals.reference = Some(PathBuf::from(required_value(name, inline_value, iter)?));
            Ok(true)
        }
        "threads" => {
            let raw = required_value(name, inline_value, iter)?;
            globals.threads = Some(parse_threads(&raw)?);
            Ok(true)
        }
        "write-index" => {
            if inline_value.is_some() {
                return Err(String::from("--write-index does not take a value"));
            }
            globals.write_index = true;
            Ok(true)
        }
        "verbosity" => {
            let raw = required_value(name, inline_value, iter)?;
            globals.verbosity = Some(parse_verbosity(&raw)?);
            Ok(true)
        }
        _ => Ok(false),
    }
}

fn required_value<I>(name: &str, inline_value: Option<&str>, iter: &mut I) -> Result<String, String>
where
    I: Iterator<Item = OsString>,
{
    if let Some(value) = inline_value {
        return Ok(value.to_owned());
    }

    iter.next()
        .and_then(|a| a.to_str().map(str::to_owned))
        .ok_or_else(|| format!("missing value for --{name}"))
}

fn parse_threads(raw: &str) -> Result<usize, String> {
    raw.parse::<usize>()
        .map_err(|_| format!("invalid --threads value \"{raw}\""))
}

pub fn parse_verbosity(raw: &str) -> Result<i32, String> {
    raw.parse::<i32>()
        .map_err(|_| format!("invalid --verbosity value \"{raw}\""))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    static LOG_LOCK: Mutex<()> = Mutex::new(());

    #[test]
    fn parses_and_strips_top_level_global_args() {
        let args = vec![
            OsString::from("samtools"),
            OsString::from("--threads"),
            OsString::from("4"),
            OsString::from("--reference=ref.fa"),
            OsString::from("--input-fmt"),
            OsString::from("bam"),
            OsString::from("--input-fmt-option"),
            OsString::from("decode_md=0"),
            OsString::from("--output-fmt=sam"),
            OsString::from("--output-fmt-option"),
            OsString::from("level=1"),
            OsString::from("--write-index"),
            OsString::from("view"),
            OsString::from("in.bam"),
        ];

        let (filtered, globals) = parse_top_level_global_args(args).unwrap();

        assert_eq!(
            filtered,
            vec![
                OsString::from("samtools"),
                OsString::from("view"),
                OsString::from("in.bam")
            ]
        );
        assert_eq!(globals.input_fmt.as_deref(), Some("bam"));
        assert_eq!(globals.input_fmt_options, vec!["decode_md=0"]);
        assert_eq!(globals.output_fmt.as_deref(), Some("sam"));
        assert_eq!(globals.output_fmt_options, vec!["level=1"]);
        assert_eq!(
            globals.reference.as_deref(),
            Some(std::path::Path::new("ref.fa"))
        );
        assert_eq!(globals.threads, Some(4));
        assert!(globals.write_index);
    }

    #[test]
    fn leaves_command_local_global_named_options_for_subcommands() {
        let args = vec![
            OsString::from("samtools"),
            OsString::from("view"),
            OsString::from("--output-fmt"),
            OsString::from("bam"),
            OsString::from("in.sam"),
        ];

        let (filtered, globals) = parse_top_level_global_args(args.clone()).unwrap();

        assert_eq!(filtered, args);
        assert_eq!(globals, SamGlobalArgs::default());
    }

    #[test]
    fn rejects_missing_required_values() {
        let args = vec![OsString::from("samtools"), OsString::from("--reference")];

        assert_eq!(
            parse_top_level_global_args(args).unwrap_err(),
            "missing value for --reference"
        );
    }

    #[test]
    fn stores_and_returns_current_global_args() {
        let _guard = LOG_LOCK.lock().unwrap();
        let globals = SamGlobalArgs {
            reference: Some(PathBuf::from("ref.fa")),
            threads: Some(2),
            ..SamGlobalArgs::default()
        };

        set_current_global_args(globals.clone());
        assert_eq!(current_global_args(), globals);
        set_current_global_args(SamGlobalArgs::default());
    }

    #[test]
    fn rejects_invalid_threads() {
        let args = vec![
            OsString::from("samtools"),
            OsString::from("--threads=loud"),
            OsString::from("view"),
        ];

        assert_eq!(
            parse_top_level_global_args(args).unwrap_err(),
            "invalid --threads value \"loud\""
        );
    }

    #[test]
    fn applies_top_level_verbosity() {
        let _guard = LOG_LOCK.lock().unwrap();
        let original = htslib_rs::log_compat::hts_verbose();
        let args = vec![
            OsString::from("samtools"),
            OsString::from("--verbosity"),
            OsString::from("6"),
            OsString::from("view"),
        ];

        let (filtered, globals) = apply_top_level_global_args(args).unwrap();

        assert_eq!(htslib_rs::log_compat::hts_verbose(), 6);
        assert_eq!(globals.verbosity, Some(6));
        assert_eq!(
            filtered,
            vec![OsString::from("samtools"), OsString::from("view")]
        );
        htslib_rs::log_compat::set_hts_verbose(original);
    }
}
