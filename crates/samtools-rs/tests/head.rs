//! Integration tests for `samtools head`.

use std::ffi::OsString;
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::Mutex;

use samtools_rs::commands::head;
use samtools_rs::run as samtools_run;

static GLOBAL_ARGS_LOCK: Mutex<()> = Mutex::new(());

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
        .join("dat")
}

fn htslib_fixtures_dir() -> PathBuf {
    let manifest = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("repos")
        .join("htslib-rs")
        .join("repos")
        .join("htslib")
        .join("test")
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

fn argv(name: &str, rest: &[&str]) -> Vec<OsString> {
    std::iter::once(OsString::from(name))
        .chain(rest.iter().map(OsString::from))
        .collect()
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

#[test]
fn head_cram_records_use_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "head",
                "-n",
                "1",
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );
}

#[test]
fn head_cram_records_without_reference_fail_cleanly() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let cram = htslib_fixtures_dir().join("range.cram");

    assert_ne!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &["head", "-n", "1", cram.to_str().unwrap()],
        ))),
        0
    );
}
