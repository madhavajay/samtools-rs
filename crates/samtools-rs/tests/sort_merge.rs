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
    // Upstream `-n` uses natural (`strnum_cmp`) order, not lexical, so
    // assert the output equals its natural-sorted self.
    let natural = |x: &&str, y: &&str| -> std::cmp::Ordering {
        let (a, b) = (x.as_bytes(), y.as_bytes());
        let (mut i, mut j) = (0usize, 0usize);
        while i < a.len() && j < b.len() {
            if a[i].is_ascii_digit() && b[j].is_ascii_digit() {
                let (s, t) = (i, j);
                while i < a.len() && a[i].is_ascii_digit() {
                    i += 1;
                }
                while j < b.len() && b[j].is_ascii_digit() {
                    j += 1;
                }
                let na = std::str::from_utf8(&a[s..i])
                    .unwrap()
                    .trim_start_matches('0');
                let nb = std::str::from_utf8(&b[t..j])
                    .unwrap()
                    .trim_start_matches('0');
                match na.len().cmp(&nb.len()).then_with(|| na.cmp(nb)) {
                    std::cmp::Ordering::Equal => {}
                    o => return o,
                }
            } else {
                match a[i].cmp(&b[j]) {
                    std::cmp::Ordering::Equal => {
                        i += 1;
                        j += 1;
                    }
                    o => return o,
                }
            }
        }
        a.len().cmp(&b.len())
    };
    let mut sorted = names.clone();
    sorted.sort_by(natural);
    assert_eq!(names, sorted);
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
    // Upstream `-n` emits SS:queryname:natural (verified vs
    // sort/name.sort.expected.sam).
    assert!(text.contains("@HD\tVN:1.6\tSO:queryname\tSS:queryname:natural\n"));
}

#[test]
fn sort_sam_output_uses_htslib_float_aux_spelling() {
    // A large-exponent f32 aux value must be re-spelled the htslib `%g`
    // way (`6.022e+23`) rather than noodles' plain decimal, now that
    // sort's SAM sink routes through sam_render::write_record.
    let tmp = tmp_dir("sort-sam-float");
    let sam = tmp.join("in.sam");
    let out = tmp.join("sorted.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:unsorted\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tfa:f:6.022e+23\n",
        ),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "sort",
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
    assert!(
        text.contains("fa:f:6.022e+23"),
        "expected htslib float spelling, got:\n{text}"
    );
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
    // sort (unlike merge) sets @HD SO:coordinate for a coordinate sort.
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
    assert!(text.starts_with("@HD\tVN:1.6\tSO:unsorted\tSS:unsorted:RG:queryname:natural\n"));
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
    let mut names: Vec<_> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap().to_string())
        .collect();
    // collate intentionally shuffles by qname hash; assert the set.
    names.sort();
    assert_eq!(names, ["a", "b"]);
    // Upstream merge preserves input[0]'s @HD verbatim (no forced
    // SO:coordinate) — verified vs the merge/* fixtures.
    assert!(text.starts_with("@HD\tVN:1.6\n"));
}

#[test]
fn merge_reads_inputs_from_file_list() {
    let tmp = tmp_dir("merge-input-list");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let list = tmp.join("inputs.txt");
    let out = tmp.join("merged.sam");
    std::fs::write(
        &sam_a,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "b\t0\tchr1\t4\t60\t4M\t*\t0\t0\tTGCA\t####\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &list,
        format!("{}\n\n{}\n", sam_a.display(), sam_b.display()),
    )
    .unwrap();

    let argv: Vec<OsString> = [
        "merge",
        "-f",
        "--no-PG",
        "-b",
        list.to_str().unwrap(),
        "--output-fmt",
        "sam",
        "-o",
        out.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(merge::main(&argv)), 0);

    let text = std::fs::read_to_string(out).unwrap();
    let mut names: Vec<_> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap().to_string())
        .collect();
    // collate intentionally shuffles by qname hash; assert the set.
    names.sort();
    assert_eq!(names, ["a", "b"]);
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
fn merge_unions_compatible_sq_metadata_fields() {
    let tmp = tmp_dir("merge-sq-fields");
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
            "@SQ\tSN:chr1\tLN:8\tM5:0123456789abcdef0123456789abcdef\n",
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
    // Upstream keeps the *first* @SQ definition for an SN verbatim (it
    // does not graft M5/etc from later inputs).
    assert!(text.contains("@SQ\tSN:chr1\tLN:8\n"));
    assert_eq!(text.matches("@SQ\tSN:chr1").count(), 1);
}

#[test]
fn merge_rejects_conflicting_sq_metadata_fields() {
    let tmp = tmp_dir("merge-conflicting-sq-fields");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("merged.sam");
    std::fs::write(
        &sam_a,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\tM5:0123456789abcdef0123456789abcdef\n",
            "a\t0\tchr1\t2\t60\t4M\t*\t0\t0\tAAAA\t!!!!\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\tM5:fedcba9876543210fedcba9876543210\n",
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
fn merge_reconciles_conflicting_read_group_headers() {
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
    // Upstream `samtools merge` does NOT reject a conflicting @RG ID — it
    // reconciles by suffixing the colliding one (gen_unique_id PRNG).
    assert_eq!(exit_to_u8(merge::main(&argv)), 0);
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("@RG\tID:rg1\tSM:s1"));
    let suffixed = text
        .lines()
        .any(|l| l.starts_with("@RG\tID:rg1-") && l.contains("SM:s2"));
    assert!(suffixed, "conflicting @RG should be suffixed:\n{text}");
    // Record `b` (from input b) has its RG:Z: remapped to the suffixed id.
    let b = text.lines().find(|l| l.starts_with("b\t")).unwrap();
    assert!(b.contains("RG:Z:rg1-"), "record RG not remapped: {b}");
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
    // input[0] has no @HD → the first @HD found (input b) is used
    // verbatim; merge does not graft SO.
    assert!(hd.contains("VN:1.5"));
    assert!(hd.contains("GO:query"));
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
    // Upstream uses input[0]'s @HD verbatim; it does not union later
    // inputs' @HD metadata (no grafted GO:query) nor graft SO.
    assert_eq!(hd, "@HD\tVN:1.6");
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
fn merge_reconciles_conflicting_program_headers() {
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
    // Upstream reconciles conflicting @PG IDs (suffix) rather than erroring.
    assert_eq!(exit_to_u8(merge::main(&argv)), 0);
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("@PG\tID:pg1\tPN:tool1"));
    assert!(
        text.lines()
            .any(|l| l.starts_with("@PG\tID:pg1-") && l.contains("PN:tool2")),
        "conflicting @PG should be suffixed:\n{text}"
    );
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
    // Upstream merge preserves input[0]'s @HD verbatim (no forced SO).
    assert!(text.starts_with("@HD\tVN:1.6\n"));
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
    let mut names: Vec<_> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap().to_string())
        .collect();
    // collate intentionally shuffles by qname hash; assert the set.
    names.sort();
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
fn collate_matches_upstream_test_collate_fixtures() {
    // Byte-exact (modulo the harness-stripped @PG) vs the entire
    // upstream test_collate harness: -o (==collate.expected), positional
    // prefix, fast `-f`, and fast `-f -r 4` (ring-eviction deferral).
    let d = fixtures_dir();
    let tmp = tmp_dir("collate-fixtures");
    let np = |s: &str| -> String {
        s.lines()
            .filter(|l| !l.starts_with("@PG\t"))
            .map(|l| format!("{l}\n"))
            .collect()
    };
    let din = d.join("dat/test_input_1_d.sam");
    let fin = d.join("collate/fast_collate.sam");
    let cases: &[(&[&str], &str)] = &[
        (
            &["--output-fmt=sam", "-o", "@OUT@", "@DIN@"],
            "collate/collate.expected.sam",
        ),
        (
            &["--output-fmt=sam", "-f", "@FIN@", "-o", "@OUT@"],
            "collate/1_fast_collate.sam.expected",
        ),
        (
            &["--output-fmt=sam", "-f", "-r", "4", "@FIN@", "-o", "@OUT@"],
            "collate/2_fast_collate_with_tmp_used.sam.expected",
        ),
    ];
    for (i, (args, expected)) in cases.iter().enumerate() {
        let out = tmp.join(format!("c{i}.sam"));
        let v: Vec<OsString> = std::iter::once(OsString::from("collate"))
            .chain(args.iter().map(|a| {
                OsString::from(match *a {
                    "@OUT@" => out.to_str().unwrap().to_string(),
                    "@DIN@" => din.to_str().unwrap().to_string(),
                    "@FIN@" => fin.to_str().unwrap().to_string(),
                    other => other.to_string(),
                })
            }))
            .collect();
        assert_eq!(exit_to_u8(collate::main(&v)), 0, "{expected}");
        assert_eq!(
            np(&std::fs::read_to_string(&out).unwrap()),
            np(&std::fs::read_to_string(d.join(expected)).unwrap()),
            "collate {expected}"
        );
    }
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
fn collate_accepts_temporary_file_count_option() {
    let tmp = tmp_dir("collate-temp-count");
    let sam = tmp.join("in.sam");
    let out = tmp.join("collated.sam");
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
        "-n",
        "3",
        "-o",
        out.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]
    .iter()
    .map(OsString::from)
    .collect();
    assert_eq!(exit_to_u8(collate::main(&argv)), 0);

    let text = std::fs::read_to_string(out).unwrap();
    let mut names: Vec<_> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap().to_string())
        .collect();
    // collate intentionally shuffles by qname hash; assert the set.
    names.sort();
    assert_eq!(names, ["a", "b"]);
}

#[test]
fn collate_rejects_invalid_temporary_file_count() {
    let argv: Vec<OsString> = ["collate", "-n", "0", sample_bam().to_str().unwrap()]
        .iter()
        .map(OsString::from)
        .collect();

    assert_eq!(exit_to_u8(collate::main(&argv)), 1);
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
    // collate groups records by qname (hash-bucketed shuffle, not
    // sorted); assert each qname's records are contiguous.
    let mut seen = std::collections::HashSet::new();
    let mut prev = "";
    for n in &names {
        if *n != prev {
            assert!(seen.insert(*n), "qname {n} not contiguous after shuffle");
            prev = n;
        }
    }
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

#[test]
fn sort_matches_upstream_test_sort_fixtures() {
    // TODO.md §13: basic test_sort group byte-exact (modulo @PG, which
    // the harness strips via ignore_pg_header). Covers coordinate, `-n`
    // natural name, `-t TAG`, and `-n -t TAG` orders + raw @HD SO/SS.
    let d = fixtures_dir();
    let tmp = tmp_dir("sort-fixtures");
    let np = |s: &str| -> String {
        s.lines()
            .filter(|l| !l.starts_with("@PG\t"))
            .map(|l| format!("{l}\n"))
            .collect()
    };
    let cases: &[(&[&str], &str, &str)] = &[
        (
            &["-m", "10M"],
            "dat/test_input_1_a.bam",
            "sort/pos.sort.expected.sam",
        ),
        (
            &["-n", "-m", "10M"],
            "dat/test_input_1_a.bam",
            "sort/name.sort.expected.sam",
        ),
        (
            &["-n", "-m", "10M"],
            "dat/sort_name_input_1.sam",
            "sort/name3.sort.expected.sam",
        ),
        (
            &["-t", "RG", "-m", "10M"],
            "dat/test_input_1_a.bam",
            "sort/tag.rg.sort.expected.sam",
        ),
        (
            &["-n", "-t", "RG", "-m", "10M"],
            "dat/test_input_1_a.bam",
            "sort/tag.rg.n.sort.expected.sam",
        ),
        // Exercises SAM `AS:I:` (uint32) aux integer-synonym tolerance.
        (
            &["-t", "AS", "-m", "10M"],
            "dat/test_input_1_d.sam",
            "sort/tag.as.sort.expected.sam",
        ),
    ];
    for (args, input, expected) in cases {
        let out = tmp.join(expected.replace('/', "_"));
        let mut v: Vec<String> = vec!["sort".into()];
        v.extend(args.iter().map(|s| s.to_string()));
        v.push("-O".into());
        v.push("SAM".into());
        v.push("-o".into());
        v.push(out.to_str().unwrap().into());
        v.push(d.join(input).to_str().unwrap().into());
        let refs: Vec<&OsString> = Vec::new();
        let _ = refs;
        let argv: Vec<OsString> = v.iter().map(OsString::from).collect();
        assert_eq!(exit_to_u8(sort::main(&argv)), 0, "{expected}");
        assert_eq!(
            np(&std::fs::read_to_string(&out).unwrap()),
            np(&std::fs::read_to_string(d.join(expected)).unwrap()),
            "sort {expected} args={args:?}"
        );
    }
}

#[test]
fn merge_reconciles_rg_pg_byte_exact_vs_upstream() {
    // TODO-NEXT merge: `merge -s 1` @RG/@PG PRNG reconciliation +
    // raw-header preservation, byte-exact vs upstream test_merge
    // fixtures (modulo @PG, which the harness strips).
    let d = fixtures_dir();
    let tmp = tmp_dir("merge-fixtures");
    let np = |s: &str| -> String {
        s.lines()
            .filter(|l| !l.starts_with("@PG\t"))
            .map(|l| format!("{l}\n"))
            .collect()
    };
    // (extra_flags, inputs, expected). merge/5 exercises `-r`
    // (filename-derived @RG attached to every record); merge/6 `-cp`
    // (combine identical @RG and @PG IDs instead of suffixing).
    let cases: &[(&[&str], &[&str], &str)] = &[
        (
            &[],
            &[
                "dat/test_input_1_a.sam",
                "dat/test_input_1_b.sam",
                "dat/test_input_1_c.sam",
            ],
            "merge/2.merge.expected.sam",
        ),
        (
            &[],
            &["dat/test_input_1_b.bam"],
            "merge/4.merge.expected.sam",
        ),
        (
            &["-r"],
            &[
                "dat/test_input_1_a.sam",
                "dat/test_input_1_b.sam",
                "dat/test_input_1_c.sam",
            ],
            "merge/5.merge.expected.sam",
        ),
        (
            &["-cp"],
            &["dat/test_input_1_a.sam", "dat/test_input_1_b.sam"],
            "merge/6.merge.expected.sam",
        ),
        (
            &[],
            &[
                "dat/test_input_1_a_regex.sam",
                "dat/test_input_1_b_regex.sam",
            ],
            "merge/7.merge.expected.sam",
        ),
    ];
    for (flags, ins, expected) in cases {
        let out = tmp.join(expected.replace('/', "_"));
        let mut v: Vec<String> = vec!["merge".into()];
        for f in *flags {
            v.push((*f).into());
        }
        v.extend([
            "-s".into(),
            "1".into(),
            "-O".into(),
            "sam".into(),
            "-o".into(),
            out.to_str().unwrap().into(),
        ]);
        for i in *ins {
            v.push(d.join(i).to_str().unwrap().into());
        }
        let argv: Vec<OsString> = v.iter().map(OsString::from).collect();
        assert_eq!(exit_to_u8(merge::main(&argv)), 0, "{expected}");
        assert_eq!(
            np(&std::fs::read_to_string(&out).unwrap()),
            np(&std::fs::read_to_string(d.join(expected)).unwrap()),
            "merge {expected}"
        );
    }
}
