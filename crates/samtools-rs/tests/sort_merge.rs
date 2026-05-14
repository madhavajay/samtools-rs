//! Integration tests for `samtools sort`, `samtools merge`, and
//! `samtools collate`, all of which share the same in-memory sort backbone.

use std::ffi::OsString;
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::Mutex;

use samtools_rs::commands::{collate, merge, sort};
use samtools_rs::run as samtools_run;

static GLOBAL_ARGS_LOCK: Mutex<()> = Mutex::new(());

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

fn htslib_fixtures_dir() -> PathBuf {
    let manifest = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("htslib-rs")
        .join("htslib")
        .join("test")
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
fn sort_cram_input_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("sort-cram");
    let out = tmp.join("sorted.sam");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    let argv: Vec<OsString> = [
        "samtools",
        "--reference",
        reference.to_str().unwrap(),
        "sort",
        "-n",
        "-O",
        "sam",
        "-o",
        out.to_str().unwrap(),
        cram.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(samtools_run(argv)), 0);

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("@HD\t"));
    assert!(text.contains("SO:queryname"));
    let names: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap())
        .collect();
    assert!(!names.is_empty());
    assert!(names.windows(2).all(|w| w[0] <= w[1]));
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
fn sort_tag_sorts_missing_first_then_numeric_tag_then_coordinate() {
    let tmp = tmp_dir("sort-tag-numeric");
    let sam = tmp.join("in.sam");
    let out = tmp.join("tag.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:20\n",
            "high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tAS:i:10\n",
            "low2\t0\tchr1\t3\t60\t4M\t*\t0\t0\tACGT\t!!!!\tAS:i:2\n",
            "missing\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "low1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tAS:i:2\n",
        ),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "sort",
        "-t",
        "AS",
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
    assert_eq!(names, ["missing", "low1", "low2", "high"]);
    assert!(text.starts_with("@HD\tVN:1.6\tSO:unsorted\tSS:unsorted:AS:coordinate\n"));
}

#[test]
fn sort_tag_with_name_sort_uses_name_secondary() {
    let tmp = tmp_dir("sort-tag-name");
    let sam = tmp.join("in.sam");
    let out = tmp.join("tag-name.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:20\n",
            "z\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:grp\n",
            "b\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "a\t0\tchr1\t3\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:grp\n",
        ),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "sort",
        "-n",
        "-t",
        "RG",
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
    assert_eq!(names, ["b", "a", "z"]);
    assert!(
        text.starts_with("@HD\tVN:1.6\tSO:unsorted\tSS:unsorted:RG:queryname:lexicographical\n")
    );
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
fn merge_unions_different_sq_headers_and_remaps_records() {
    let tmp = tmp_dir("merge-union-sq");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    std::fs::write(
        &sam_a,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "a\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr2\tLN:9\n",
            "b\t0\tchr2\t1\t60\t4M\t*\t0\t0\tCCCC\t####\n",
        ),
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
    assert!(text.contains("@SQ\tSN:chr1\tLN:8\n"));
    assert!(text.contains("@SQ\tSN:chr2\tLN:9\n"));
    assert!(text.contains("\na\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\n"));
    assert!(text.contains("\nb\t0\tchr2\t1\t60\t4M\t*\t0\t0\tCCCC\t####\n"));
}

#[test]
fn merge_rejects_conflicting_sq_lengths() {
    let tmp = tmp_dir("merge-conflicting-sq");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    std::fs::write(
        &sam_a,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "a\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:9\n",
            "b\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t####\n",
        ),
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
    assert_eq!(exit_to_u8(merge::main(&argv)), 1);
    assert!(!out.exists());
}

#[test]
fn merge_unions_read_group_headers() {
    let tmp = tmp_dir("merge-union-rg");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    std::fs::write(
        &sam_a,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:rg1\tSM:s1\n",
            "a\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\tRG:Z:rg1\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:rg2\tSM:s2\n",
            "b\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t####\tRG:Z:rg2\n",
        ),
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
    assert!(text.contains("@RG\tID:rg1\tSM:s1\n"));
    assert!(text.contains("@RG\tID:rg2\tSM:s2\n"));
    assert!(text.contains("\na\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\tRG:Z:rg1\n"));
    assert!(text.contains("\nb\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t####\tRG:Z:rg2\n"));
}

#[test]
fn merge_rejects_conflicting_read_group_headers() {
    let tmp = tmp_dir("merge-conflicting-rg");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    std::fs::write(
        &sam_a,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:rg1\tSM:s1\n",
            "a\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\tRG:Z:rg1\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:rg1\tSM:s2\n",
            "b\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t####\tRG:Z:rg1\n",
        ),
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
    assert_eq!(exit_to_u8(merge::main(&argv)), 1);
    assert!(!out.exists());
}

#[test]
fn merge_appends_comments_from_all_headers() {
    let tmp = tmp_dir("merge-comments");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    std::fs::write(
        &sam_a,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@CO\tfirst input comment\n",
            "a\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@CO\tsecond input comment\n",
            "b\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t####\n",
        ),
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
    assert!(text.contains("@CO\tfirst input comment\n"));
    assert!(text.contains("@CO\tsecond input comment\n"));
}

#[test]
fn merge_preserves_later_header_metadata_when_first_input_lacks_hd() {
    let tmp = tmp_dir("merge-later-hd");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    std::fs::write(
        &sam_a,
        concat!(
            "@SQ\tSN:chr1\tLN:8\n",
            "a\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        concat!(
            "@HD\tVN:1.5\tGO:query\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "b\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t####\n",
        ),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "merge",
        "-f",
        "--no-PG",
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
    let hd = text.lines().find(|line| line.starts_with("@HD")).unwrap();
    assert!(hd.contains("VN:1.5"));
    assert!(hd.contains("GO:query"));
    assert!(hd.contains("SO:coordinate"));
}

#[test]
fn merge_unions_compatible_header_metadata_fields() {
    let tmp = tmp_dir("merge-hd-fields");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    std::fs::write(
        &sam_a,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "a\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        concat!(
            "@HD\tVN:1.6\tGO:query\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "b\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t####\n",
        ),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "merge",
        "-f",
        "--no-PG",
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
    let hd = text.lines().find(|line| line.starts_with("@HD")).unwrap();
    assert!(hd.contains("VN:1.6"));
    assert!(hd.contains("GO:query"));
    assert!(hd.contains("SO:coordinate"));
}

#[test]
fn merge_rejects_conflicting_header_metadata_fields() {
    let tmp = tmp_dir("merge-hd-conflict");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    std::fs::write(
        &sam_a,
        concat!(
            "@HD\tVN:1.6\tGO:query\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "a\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        concat!(
            "@HD\tVN:1.6\tGO:reference\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "b\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t####\n",
        ),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "merge",
        "-f",
        "--no-PG",
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
    assert_eq!(exit_to_u8(merge::main(&argv)), 1);
    assert!(!out.exists());
}

#[test]
fn merge_unions_program_headers() {
    let tmp = tmp_dir("merge-union-pg");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    std::fs::write(
        &sam_a,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@PG\tID:pg1\tPN:tool1\n",
            "a\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@PG\tID:pg2\tPN:tool2\n",
            "b\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t####\n",
        ),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "merge",
        "-f",
        "--no-PG",
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
    assert!(text.contains("@PG\tID:pg1\tPN:tool1\n"));
    assert!(text.contains("@PG\tID:pg2\tPN:tool2\n"));
}

#[test]
fn merge_rejects_conflicting_program_headers() {
    let tmp = tmp_dir("merge-conflicting-pg");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    std::fs::write(
        &sam_a,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@PG\tID:pg1\tPN:tool1\n",
            "a\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@PG\tID:pg1\tPN:tool2\n",
            "b\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t####\n",
        ),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "merge",
        "-f",
        "--no-PG",
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
    assert_eq!(exit_to_u8(merge::main(&argv)), 1);
    assert!(!out.exists());
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
fn merge_tag_sort_orders_by_aux_tag_and_accepts_s_option() {
    let tmp = tmp_dir("merge-tag");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:8\n";
    std::fs::write(
        &sam_a,
        format!(
            "{header}missing\t0\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\t!!!!\n\
             high\t0\tchr1\t2\t60\t4M\t*\t0\t0\tCCCC\t####\tZZ:i:7\n"
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        format!("{header}low\t0\tchr1\t3\t60\t4M\t*\t0\t0\tGGGG\t$$$$\tZZ:i:3\n"),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "merge",
        "-f",
        "-s",
        "1",
        "-t",
        "ZZ",
        "--output-fmt=sam",
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
    assert!(text.starts_with("@HD\tVN:1.6\tSO:unsorted\tSS:unsorted:ZZ:coordinate\n"));
    let names: Vec<_> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(names, ["missing", "low", "high"]);
}

#[test]
fn merge_dash_output_writes_sam_to_stdout() {
    let tmp = tmp_dir("merge-dash-output");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
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
        "SAM",
        "-",
        sam_a.to_str().unwrap(),
        sam_b.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(merge::main(&argv)), 0);
    assert!(!tmp.join("-").exists());
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
fn collate_positional_prefix_writes_format_extension() {
    let tmp = tmp_dir("collate-prefix");
    let sam = tmp.join("in.sam");
    let out_prefix = tmp.join("legacy");
    let out = tmp.join("legacy.sam");
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
        "collate",
        "--output-fmt=sam",
        sam.to_str().unwrap(),
        out_prefix.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(collate::main(&argv)), 0);

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("@HD\tVN:1.6\tSO:unsorted\tGO:query\n"));
    let names: Vec<_> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(names, ["a", "b"]);
}

#[test]
fn collate_rejects_output_file_with_stdout() {
    let argv: Vec<OsString> = [
        "collate",
        "-O",
        "-o",
        "out.sam",
        sample_bam().to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();

    assert_eq!(exit_to_u8(collate::main(&argv)), 1);
}

#[test]
fn collate_fast_mode_pairs_primary_reads_and_drops_supplementary() {
    let tmp = tmp_dir("collate-fast");
    let sam = tmp.join("in.sam");
    let out = tmp.join("fast.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:20\n",
            "pair\t147\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\n",
            "supp\t2113\tchr1\t6\t60\t4M\t*\t0\t0\tCCCC\t$$$$\n",
            "solo\t65\tchr1\t9\t60\t4M\t*\t0\t0\tGGGG\t%%%%\n",
            "pair\t99\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "collate",
        "--output-fmt=sam",
        "-f",
        "-r",
        "2",
        "-o",
        out.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(collate::main(&argv)), 0);

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("@HD\tVN:1.6\tSO:unsorted\tGO:query\n"));
    let records: Vec<(&str, u16)> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| {
            let mut fields = line.split('\t');
            let name = fields.next().unwrap();
            let flag = fields.next().unwrap().parse::<u16>().unwrap();
            (name, flag)
        })
        .collect();
    assert_eq!(records, [("pair", 99), ("pair", 147), ("solo", 65)]);
}

#[test]
fn collate_cram_input_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("collate-cram");
    let out_prefix = tmp.join("collated");
    let out = tmp.join("collated.sam");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    let argv: Vec<OsString> = [
        "samtools",
        "--reference",
        reference.to_str().unwrap(),
        "collate",
        "--output-fmt",
        "sam",
        "-o",
        out_prefix.to_str().unwrap(),
        cram.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(samtools_run(argv)), 0);

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("@HD\t"));
    let names: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap())
        .collect();
    assert!(!names.is_empty());
    assert!(names.windows(2).all(|w| w[0] <= w[1]));
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

#[test]
fn sort_adds_pg_line_by_default_and_omits_with_no_pg() {
    let tmp = tmp_dir("sort-pg");
    let sam = tmp.join("in.sam");
    let default_out = tmp.join("default.sam");
    let no_pg_out = tmp.join("no_pg.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@PG\tID:upstream\tPN:upstream\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    let argv_default: Vec<OsString> = [
        "sort",
        "-O",
        "sam",
        "-o",
        default_out.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(sort::main(&argv_default)), 0);
    let default_text = std::fs::read_to_string(&default_out).unwrap();
    assert!(default_text.contains("PN:samtools"));
    assert!(default_text.contains("PP:upstream"));

    let argv_no_pg: Vec<OsString> = [
        "sort",
        "--no-PG",
        "-O",
        "sam",
        "-o",
        no_pg_out.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(sort::main(&argv_no_pg)), 0);
    let no_pg_text = std::fs::read_to_string(&no_pg_out).unwrap();
    assert!(!no_pg_text.contains("PN:samtools"));
    assert!(no_pg_text.contains("@PG\tID:upstream"));
}

#[test]
fn merge_adds_pg_line_by_default_and_omits_with_no_pg() {
    let tmp = tmp_dir("merge-pg");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let default_out = tmp.join("default.sam");
    let no_pg_out = tmp.join("no_pg.sam");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:8\n";
    std::fs::write(
        &sam_a,
        format!("{header}a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n"),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        format!("{header}b\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\n"),
    )
    .unwrap();

    let argv_default: Vec<OsString> = [
        "merge",
        "-f",
        "--output-fmt",
        "sam",
        "-o",
        default_out.to_str().unwrap(),
        sam_a.to_str().unwrap(),
        sam_b.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(merge::main(&argv_default)), 0);
    assert!(
        std::fs::read_to_string(&default_out)
            .unwrap()
            .contains("PN:samtools")
    );

    let argv_no_pg: Vec<OsString> = [
        "merge",
        "-f",
        "--no-PG",
        "--output-fmt",
        "sam",
        "-o",
        no_pg_out.to_str().unwrap(),
        sam_a.to_str().unwrap(),
        sam_b.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(merge::main(&argv_no_pg)), 0);
    assert!(
        !std::fs::read_to_string(&no_pg_out)
            .unwrap()
            .contains("PN:samtools")
    );
}

#[test]
fn collate_adds_pg_line_by_default_and_omits_with_no_pg() {
    let tmp = tmp_dir("collate-pg");
    let sam = tmp.join("in.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    let default_prefix = tmp.join("default");
    let default_out = tmp.join("default.sam");
    let argv_default: Vec<OsString> = [
        "collate",
        "--output-fmt",
        "sam",
        "-o",
        default_prefix.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(collate::main(&argv_default)), 0);
    assert!(
        std::fs::read_to_string(&default_out)
            .unwrap()
            .contains("PN:samtools")
    );

    let no_pg_prefix = tmp.join("no_pg");
    let no_pg_out = tmp.join("no_pg.sam");
    let argv_no_pg: Vec<OsString> = [
        "collate",
        "--no-PG",
        "--output-fmt",
        "sam",
        "-o",
        no_pg_prefix.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(collate::main(&argv_no_pg)), 0);
    assert!(
        !std::fs::read_to_string(&no_pg_out)
            .unwrap()
            .contains("PN:samtools")
    );
}

#[test]
fn merge_r_region_restricts_to_indexed_records() {
    let tmp = tmp_dir("merge-region");
    let bam = sample_bam();
    // Build a BAI alongside the input so query_bam_records_from_path can hit it.
    let indexed = tmp.join("indexed.bam");
    std::fs::copy(&bam, &indexed).unwrap();
    let idx_argv: Vec<OsString> = ["index", indexed.to_str().unwrap()]
        .iter()
        .map(OsString::from)
        .collect();
    assert_eq!(exit_to_u8(samtools_rs::commands::index::main(&idx_argv)), 0);

    let full_out = tmp.join("full.sam");
    let region_out = tmp.join("region.sam");
    let argv_full: Vec<OsString> = [
        "merge",
        "-f",
        "--output-fmt",
        "sam",
        "-o",
        full_out.to_str().unwrap(),
        indexed.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(merge::main(&argv_full)), 0);

    let argv_region: Vec<OsString> = [
        "merge",
        "-f",
        "-R",
        "17:1-2000",
        "--output-fmt",
        "sam",
        "-o",
        region_out.to_str().unwrap(),
        indexed.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(merge::main(&argv_region)), 0);

    let full = std::fs::read_to_string(&full_out).unwrap();
    let region = std::fs::read_to_string(&region_out).unwrap();
    let full_records = full.lines().filter(|l| !l.starts_with('@')).count();
    let region_records = region.lines().filter(|l| !l.starts_with('@')).count();
    assert!(full_records > region_records);
    assert!(region_records > 0);
}

#[test]
fn merge_l_bed_restricts_to_indexed_records_and_deduplicates_overlaps() {
    let tmp = tmp_dir("merge-bed");
    let bam = sample_bam();
    let indexed = tmp.join("indexed.bam");
    std::fs::copy(&bam, &indexed).unwrap();
    let idx_argv: Vec<OsString> = ["index", indexed.to_str().unwrap()]
        .iter()
        .map(OsString::from)
        .collect();
    assert_eq!(exit_to_u8(samtools_rs::commands::index::main(&idx_argv)), 0);

    let bed = tmp.join("regions.bed");
    std::fs::write(&bed, "17\t0\t2000\n17\t1000\t2500\n").unwrap();

    let full_out = tmp.join("full.sam");
    let bed_out = tmp.join("bed.sam");
    let argv_full: Vec<OsString> = [
        "merge",
        "-f",
        "--output-fmt",
        "sam",
        "-o",
        full_out.to_str().unwrap(),
        indexed.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(merge::main(&argv_full)), 0);

    let argv_bed: Vec<OsString> = [
        "merge",
        "-f",
        "-L",
        bed.to_str().unwrap(),
        "--output-fmt",
        "sam",
        "-o",
        bed_out.to_str().unwrap(),
        indexed.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(merge::main(&argv_bed)), 0);

    let full = std::fs::read_to_string(&full_out).unwrap();
    let bed_text = std::fs::read_to_string(&bed_out).unwrap();
    let full_records = full.lines().filter(|l| !l.starts_with('@')).count();
    let bed_records: Vec<_> = bed_text.lines().filter(|l| !l.starts_with('@')).collect();
    assert!(full_records > bed_records.len());
    assert!(!bed_records.is_empty());

    let unique: std::collections::HashSet<_> = bed_records.iter().copied().collect();
    assert_eq!(unique.len(), bed_records.len());
}
