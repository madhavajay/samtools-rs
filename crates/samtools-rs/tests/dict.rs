//! Integration tests for `samtools dict`.
//!
//! Runs the library entry point against the upstream `samtools/test/dat/`
//! fixtures and asserts byte-for-byte parity with the checked-in expected
//! output files.

use std::ffi::OsString;
use std::path::PathBuf;

use samtools_rs::commands::dict;

fn fixtures_dir() -> PathBuf {
    let manifest = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("samtools")
        .join("test")
        .join("dat")
}

fn run_dict_to_file(args: &[&str]) -> Vec<u8> {
    let tmp = tempdir();
    let out = tmp.join("dict.out");
    let argv: Vec<OsString> = std::iter::once(OsString::from("dict"))
        .chain(args.iter().map(OsString::from))
        .chain(["-o".into(), out.to_string_lossy().into_owned().into()])
        .collect();
    let code = dict::main(&argv);
    let s = format!("{:?}", code);
    assert!(
        s.contains("ExitCode(unix_exit_status(0))")
            || s.contains("ExitCode(0)")
            || s == "ExitCode(unix_exit_status(0))",
        "dict failed: {:?}",
        code
    );
    std::fs::read(&out).expect("read dict output")
}

fn tempdir() -> PathBuf {
    let base = std::env::temp_dir().join(format!("samtools-rs-dict-{}", std::process::id()));
    std::fs::create_dir_all(&base).unwrap();
    base
}

#[test]
fn dict_matches_dict_out() {
    let dir = fixtures_dir();
    let input = dir.join("dict.fa");
    let expected = std::fs::read(dir.join("dict.out")).unwrap();
    let actual = run_dict_to_file(&[
        "-a",
        "hf37d5",
        "-s",
        "Homo floresiensis",
        "-u",
        "ftp://example.com/hf37d5.fa.gz",
        input.to_str().unwrap(),
    ]);
    assert_eq!(
        String::from_utf8_lossy(&actual),
        String::from_utf8_lossy(&expected),
        "dict output should match checked-in dict.out byte-for-byte"
    );
}
