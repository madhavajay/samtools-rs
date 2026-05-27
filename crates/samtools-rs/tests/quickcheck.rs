//! Integration tests for `samtools quickcheck`.
//!
//! Drives the library entry point directly with absolute paths to the
//! upstream `repos/samtools/test/quickcheck/` fixtures, then asserts the
//! per-file exit code. The CLI-level `quickcheck` test locks the
//! byte-for-byte `-v` output against `quickcheck/all.expected`.

use std::ffi::OsString;
use std::path::PathBuf;
use std::process::ExitCode;

use samtools_rs::commands::quickcheck;

fn fixtures_dir() -> PathBuf {
    let manifest = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("repos")
        .join("samtools")
        .join("test")
        .join("quickcheck")
}

fn run_cli(args: &[&str]) -> ExitCode {
    let argv: Vec<OsString> = std::iter::once(OsString::from("quickcheck"))
        .chain(args.iter().map(OsString::from))
        .collect();
    quickcheck::main(&argv)
}

fn exit_to_u8(code: ExitCode) -> u8 {
    // ExitCode doesn't expose the inner u8 directly; round-trip via Debug.
    let s = format!("{:?}", code);
    let n = s.chars().filter(|c| c.is_ascii_digit()).collect::<String>();
    n.parse().unwrap_or(255)
}

#[test]
fn ok_bam_passes() {
    let path = fixtures_dir().join("3.quickcheck.ok.bam");
    let code = run_cli(&[path.to_str().unwrap()]);
    assert_eq!(exit_to_u8(code), 0);
}

#[test]
fn missing_eof_fails() {
    let path = fixtures_dir().join("1.quickcheck.badeof.bam");
    let code = run_cli(&[path.to_str().unwrap()]);
    assert_eq!(exit_to_u8(code), 16);
}

#[test]
fn bad_header_fails() {
    let path = fixtures_dir().join("2.quickcheck.badheader.bam");
    let code = run_cli(&[path.to_str().unwrap()]);
    assert_eq!(exit_to_u8(code), 4);
}

#[test]
fn cram21_ok_passes() {
    let path = fixtures_dir().join("6.quickcheck.cram21.ok.cram");
    let code = run_cli(&[path.to_str().unwrap()]);
    assert_eq!(exit_to_u8(code), 0);
}

#[test]
fn cram30_truncated_fails() {
    let path = fixtures_dir().join("9.quickcheck.cram30.truncated.cram");
    let code = run_cli(&[path.to_str().unwrap()]);
    assert_eq!(exit_to_u8(code), 16);
}

#[test]
fn notargets_fails_without_u_passes_with_u() {
    let path = fixtures_dir().join("10.quickcheck.notargets.bam");
    let code = run_cli(&[path.to_str().unwrap()]);
    assert_eq!(exit_to_u8(code), 8);
    let code = run_cli(&["-u", path.to_str().unwrap()]);
    assert_eq!(exit_to_u8(code), 0);
}
