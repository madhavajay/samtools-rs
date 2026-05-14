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
fn sort_write_index_builds_bai_for_coordinate_bam_output() {
    let tmp = tmp_dir("sort-write-index");
    let out = tmp.join("sorted.bam");
    let argv: Vec<OsString> = ["sort", "--write-index", "-o"]
        .iter()
        .map(OsString::from)
        .chain([
            out.to_string_lossy().into_owned().into(),
            sample_bam().to_string_lossy().into_owned().into(),
        ])
        .collect();
    assert_eq!(exit_to_u8(sort::main(&argv)), 0);
    assert!(out.exists());
    assert!(tmp.join("sorted.bam.bai").exists());
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
fn sort_sam_input_by_name_to_sam_output() {
    let tmp = tmp_dir("sort-sam-name");
    let sam = tmp.join("in.sam");
    let out = tmp.join("named.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "z\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "sort",
        "-n",
        "--output-fmt",
        "sam",
        "-o",
        out.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(sort::main(&argv)), 0);

    let text = std::fs::read_to_string(out).unwrap();
    let names: Vec<_> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(names, ["a", "z"]);
    assert!(text.contains("@HD\tVN:1.6\tSO:queryname\n"));
}

#[test]
fn sort_short_output_format_consumes_value() {
    let tmp = tmp_dir("sort-short-output-fmt");
    let sam = tmp.join("in.sam");
    let out = tmp.join("sorted.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "b\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "sort",
        "-O",
        "sam",
        "-o",
        out.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(sort::main(&argv)), 0);

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("@HD\tVN:1.6\tSO:coordinate\n"));
    assert!(text.contains("\na\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n"));
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
fn merge_write_index_builds_bai_for_coordinate_bam_output() {
    let tmp = tmp_dir("merge-write-index");
    let out = tmp.join("merged.bam");
    let bam = sample_bam();
    let argv: Vec<OsString> = ["merge", "-f", "--write-index", "-o"]
        .iter()
        .map(OsString::from)
        .chain([
            out.to_string_lossy().into_owned().into(),
            bam.to_string_lossy().into_owned().into(),
            bam.to_string_lossy().into_owned().into(),
        ])
        .collect();
    assert_eq!(exit_to_u8(merge::main(&argv)), 0);
    assert!(out.exists());
    assert!(tmp.join("merged.bam.bai").exists());
}

#[test]
fn merge_sam_inputs_to_sam_output() {
    let tmp = tmp_dir("merge-sam");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:8\n";
    std::fs::write(
        &sam_a,
        format!("{header}b\t0\tchr1\t4\t60\t4M\t*\t0\t0\tTGCA\t####\n"),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        format!("{header}a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n"),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "merge",
        "-f",
        "--output-fmt",
        "sam",
        "-o",
        out.to_str().unwrap(),
        sam_a.to_str().unwrap(),
        sam_b.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(merge::main(&argv)), 0);

    let text = std::fs::read_to_string(out).unwrap();
    let names: Vec<_> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(names, ["a", "b"]);
    assert!(text.contains("@HD\tVN:1.6\tSO:coordinate\n"));
}

#[test]
fn merge_short_output_format_consumes_value() {
    let tmp = tmp_dir("merge-short-output-fmt");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:8\n";
    std::fs::write(
        &sam_a,
        format!("{header}b\t0\tchr1\t4\t60\t4M\t*\t0\t0\tTGCA\t####\n"),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        format!("{header}a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n"),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "merge",
        "-f",
        "-O",
        "sam",
        "-o",
        out.to_str().unwrap(),
        sam_a.to_str().unwrap(),
        sam_b.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(merge::main(&argv)), 0);

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("@HD\tVN:1.6\tSO:coordinate\n"));
    assert!(text.contains("\na\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n"));
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

#[test]
fn collate_sam_input_groups_by_name_to_sam_output() {
    let tmp = tmp_dir("collate-sam-name");
    let sam = tmp.join("in.sam");
    let out_prefix = tmp.join("grouped");
    let out = tmp.join("grouped.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "z\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "z\t0\tchr1\t3\t60\t4M\t*\t0\t0\tCCCC\t$$$$\n",
        ),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "collate",
        "--output-fmt",
        "sam",
        "-o",
        out_prefix.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(collate::main(&argv)), 0);

    let text = std::fs::read_to_string(out).unwrap();
    let names: Vec<_> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(names, ["a", "z", "z"]);
}

#[test]
fn collate_rejects_invalid_output_format() {
    let argv: Vec<OsString> = [
        "collate",
        "--output-fmt",
        "cram",
        sample_bam().to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();

    assert_eq!(exit_to_u8(collate::main(&argv)), 1);
}
