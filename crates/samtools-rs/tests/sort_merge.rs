//! Integration tests for `samtools sort`, `samtools merge`, and
//! `samtools collate`, all of which share the same in-memory sort backbone.

use std::ffi::OsString;
use std::path::PathBuf;
use std::process::ExitCode;

use samtools_rs::commands::{collate, merge, sort};

fn fixtures_dir() -> PathBuf {
    let manifest = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("samtools")
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

fn tmp_dir(name: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!("samtools-rs-{}-{}", name, std::process::id()));
    std::fs::create_dir_all(&p).unwrap();
    p
}

fn sample_bam() -> PathBuf {
    fixtures_dir().join("checksum").join("chk1.bam")
}

#[test]
fn sort_coordinate_succeeds() {
    let tmp = tmp_dir("sort1");
    let out = tmp.join("sorted.bam");
    let argv: Vec<OsString> = ["sort", "-o"]
        .iter()
        .map(OsString::from)
        .chain([
            out.to_string_lossy().into_owned().into(),
            sample_bam().to_string_lossy().into_owned().into(),
        ])
        .collect();
    assert_eq!(exit_to_u8(sort::main(&argv)), 0);
    assert!(out.exists());
    assert!(out.metadata().unwrap().len() > 0);
}

#[test]
fn sort_name_succeeds() {
    let tmp = tmp_dir("sort2");
    let out = tmp.join("named.bam");
    let argv: Vec<OsString> = ["sort", "-n", "-o"]
        .iter()
        .map(OsString::from)
        .chain([
            out.to_string_lossy().into_owned().into(),
            sample_bam().to_string_lossy().into_owned().into(),
        ])
        .collect();
    assert_eq!(exit_to_u8(sort::main(&argv)), 0);
}

#[test]
fn merge_two_succeeds() {
    let tmp = tmp_dir("merge1");
    let out = tmp.join("merged.bam");
    let bam = sample_bam();
    let argv: Vec<OsString> = ["merge", "-f", "-o"]
        .iter()
        .map(OsString::from)
        .chain([
            out.to_string_lossy().into_owned().into(),
            bam.to_string_lossy().into_owned().into(),
            bam.to_string_lossy().into_owned().into(),
        ])
        .collect();
    assert_eq!(exit_to_u8(merge::main(&argv)), 0);
}

#[test]
fn collate_stdout_succeeds() {
    let argv: Vec<OsString> = [
        "collate",
        "-O",
        "--output-fmt",
        "sam",
        sample_bam().to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(collate::main(&argv)), 0);
}
