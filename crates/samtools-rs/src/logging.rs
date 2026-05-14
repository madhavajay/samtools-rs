//! Logging bridge to `htslib-rs`.

use std::ffi::OsString;

use crate::sam_global::parse_verbosity;

/// Applies top-level `--verbosity` options and returns argv without them.
pub fn apply_global_verbosity(args: Vec<OsString>) -> Result<Vec<OsString>, String> {
    let mut out = Vec::with_capacity(args.len());
    let mut iter = args.into_iter();

    if let Some(program) = iter.next() {
        out.push(program);
    }

    while let Some(arg) = iter.next() {
        let Some(s) = arg.to_str() else {
            out.push(arg);
            continue;
        };

        if s == "--verbosity" {
            let value = iter
                .next()
                .and_then(|a| a.to_str().map(str::to_owned))
                .ok_or_else(|| String::from("missing value for --verbosity"))?;
            set_verbosity(&value)?;
        } else if let Some(value) = s.strip_prefix("--verbosity=") {
            set_verbosity(value)?;
        } else {
            out.push(arg);
        }
    }

    Ok(out)
}

fn set_verbosity(raw: &str) -> Result<(), String> {
    let value = parse_verbosity(raw)?;
    htslib_rs::log_compat::set_hts_verbose(value);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    static LOG_LOCK: Mutex<()> = Mutex::new(());

    #[test]
    fn applies_separate_verbosity_arg() {
        let _guard = LOG_LOCK.lock().unwrap();
        let original = htslib_rs::log_compat::hts_verbose();
        let args = vec![
            OsString::from("samtools"),
            OsString::from("--verbosity"),
            OsString::from("5"),
            OsString::from("view"),
        ];

        let filtered = apply_global_verbosity(args).unwrap();

        assert_eq!(htslib_rs::log_compat::hts_verbose(), 5);
        assert_eq!(
            filtered,
            vec![OsString::from("samtools"), OsString::from("view")]
        );
        htslib_rs::log_compat::set_hts_verbose(original);
    }

    #[test]
    fn applies_equals_verbosity_arg() {
        let _guard = LOG_LOCK.lock().unwrap();
        let original = htslib_rs::log_compat::hts_verbose();
        let args = vec![
            OsString::from("samtools"),
            OsString::from("--verbosity=1"),
            OsString::from("quickcheck"),
        ];

        let filtered = apply_global_verbosity(args).unwrap();

        assert_eq!(htslib_rs::log_compat::hts_verbose(), 1);
        assert_eq!(
            filtered,
            vec![OsString::from("samtools"), OsString::from("quickcheck")]
        );
        htslib_rs::log_compat::set_hts_verbose(original);
    }

    #[test]
    fn rejects_invalid_verbosity_arg() {
        let _guard = LOG_LOCK.lock().unwrap();
        let args = vec![
            OsString::from("samtools"),
            OsString::from("--verbosity"),
            OsString::from("loud"),
            OsString::from("view"),
        ];

        assert_eq!(
            apply_global_verbosity(args).unwrap_err(),
            "invalid --verbosity value \"loud\""
        );
    }
}
