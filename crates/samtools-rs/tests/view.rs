//! Integration tests for `samtools view`.

use std::ffi::OsString;
use std::path::PathBuf;
use std::process::ExitCode;

use samtools_rs::commands::view;

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

fn exit_to_u8(code: ExitCode) -> u8 {
    format!("{:?}", code)
        .chars()
        .filter(|c| c.is_ascii_digit())
        .collect::<String>()
        .parse()
        .unwrap_or(255)
}

fn run(args: &[&str]) -> u8 {
    let argv: Vec<OsString> = std::iter::once(OsString::from("view"))
        .chain(args.iter().map(OsString::from))
        .collect();
    exit_to_u8(view::main(&argv))
}

#[test]
fn view_count_succeeds() {
    let p = fixtures_dir().join("view.001.sam");
    assert_eq!(run(&["-c", p.to_str().unwrap()]), 0);
}

#[test]
fn view_header_only_succeeds() {
    let p = fixtures_dir().join("view.001.sam");
    assert_eq!(run(&["-H", p.to_str().unwrap()]), 0);
}

#[test]
fn view_cram_header_only_succeeds() {
    let p = fixtures_dir().join("test_input_1_a.cram");
    assert_eq!(run(&["-H", p.to_str().unwrap()]), 0);
}

#[test]
fn view_filter_unmapped_succeeds() {
    let p = fixtures_dir().join("view.001.sam");
    // -f 4 (require unmapped) should succeed; output may be empty or contain
    // unmapped records.
    assert_eq!(run(&["-c", "-f", "4", p.to_str().unwrap()]), 0);
}
