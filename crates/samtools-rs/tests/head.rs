//! Integration tests for `samtools head`.

use std::ffi::OsString;
use std::path::PathBuf;
use std::process::ExitCode;

use samtools_rs::commands::head;

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
    let argv: Vec<OsString> = std::iter::once(OsString::from("head"))
        .chain(args.iter().map(OsString::from))
        .collect();
    exit_to_u8(head::main(&argv))
}

#[test]
fn head_all_headers_succeeds() {
    let p = fixtures_dir().join("view.001.sam");
    assert_eq!(run(&[p.to_str().unwrap()]), 0);
}

#[test]
fn head_h_5_succeeds() {
    let p = fixtures_dir().join("view.001.sam");
    assert_eq!(run(&["-h", "5", p.to_str().unwrap()]), 0);
}

#[test]
fn head_cram_headers_succeeds() {
    let p = fixtures_dir().join("test_input_1_a.cram");
    assert_eq!(run(&[p.to_str().unwrap()]), 0);
    assert_eq!(run(&["-h", "2", p.to_str().unwrap()]), 0);
}
