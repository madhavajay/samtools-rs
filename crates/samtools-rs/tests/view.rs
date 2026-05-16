//! Integration tests for `samtools view`.

use std::ffi::OsString;
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use samtools_rs::commands::view;
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
        .join("dat")
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

fn tmp_dir(name: &str) -> PathBuf {
    static NEXT_TMP_ID: AtomicUsize = AtomicUsize::new(0);

    let id = NEXT_TMP_ID.fetch_add(1, Ordering::Relaxed);
    let p = std::env::temp_dir().join(format!(
        "samtools-rs-view-{}-{}-{}",
        name,
        std::process::id(),
        id
    ));
    std::fs::create_dir_all(&p).unwrap();
    p
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

fn argv(name: &str, rest: &[&str]) -> Vec<OsString> {
    std::iter::once(OsString::from(name))
        .chain(rest.iter().map(OsString::from))
        .collect()
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
fn view_qname_file_filters_records_by_name() {
    let tmp = tmp_dir("qname-filter");
    let input = tmp.join("in.sam");
    let names_keep = tmp.join("keep.txt");
    let names_drop = tmp.join("drop.txt");
    let out_keep = tmp.join("keep.sam");
    let out_drop = tmp.join("drop.sam");

    std::fs::write(
        &input,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "alpha\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "beta\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "gamma\t0\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\t****\n",
        ),
    )
    .unwrap();
    std::fs::write(&names_keep, "alpha\ngamma\n").unwrap();
    std::fs::write(&names_drop, "beta\n").unwrap();

    assert_eq!(
        run(&[
            "-h",
            "-N",
            names_keep.to_str().unwrap(),
            "-o",
            out_keep.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    let kept = std::fs::read_to_string(&out_keep).unwrap();
    assert!(kept.contains("\nalpha\t"));
    assert!(!kept.contains("\nbeta\t"));
    assert!(kept.contains("\ngamma\t"));

    let neg_arg = format!("^{}", names_drop.display());
    assert_eq!(
        run(&[
            "-h",
            "-N",
            &neg_arg,
            "-o",
            out_drop.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    let kept_neg = std::fs::read_to_string(&out_drop).unwrap();
    assert!(kept_neg.contains("\nalpha\t"));
    assert!(!kept_neg.contains("\nbeta\t"));
    assert!(kept_neg.contains("\ngamma\t"));
}

#[test]
fn view_r_and_dash_cap_r_filter_by_read_group() {
    let tmp = tmp_dir("rg-filter");
    let input = tmp.join("in.sam");
    let rg_file = tmp.join("rgs.txt");
    let out_one = tmp.join("one.sam");
    let out_many = tmp.join("many.sam");
    let out_no_rg = tmp.join("no_rg.sam");

    std::fs::write(
        &input,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:grp1\tSM:s1\n",
            "@RG\tID:grp2\tSM:s2\n",
            "@RG\tID:grp3\tSM:s3\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:grp1\n",
            "b\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\tRG:Z:grp2\n",
            "c\t0\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\t****\tRG:Z:grp3\n",
            "d\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t&&&&\n",
        ),
    )
    .unwrap();
    std::fs::write(&rg_file, "grp1\ngrp3\n").unwrap();

    assert_eq!(
        run(&[
            "-r",
            "grp2",
            "-o",
            out_one.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    let one = std::fs::read_to_string(&out_one).unwrap();
    assert!(one.contains("\nb\t") || one.starts_with("b\t"));
    assert!(!one.contains("\na\t") && !one.starts_with("a\t"));
    assert!(!one.contains("\nc\t") && !one.starts_with("c\t"));
    assert!(!one.contains("\nd\t") && !one.starts_with("d\t"));

    assert_eq!(
        run(&[
            "-R",
            rg_file.to_str().unwrap(),
            "-o",
            out_many.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    let many = std::fs::read_to_string(&out_many).unwrap();
    assert!(many.contains("a\t"));
    assert!(!many.contains("\nb\t") && !many.starts_with("b\t"));
    assert!(many.contains("c\t"));
    assert!(!many.contains("\nd\t") && !many.starts_with("d\t"));

    assert_eq!(
        run(&[
            "-n",
            "-o",
            out_no_rg.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    let no_rg = std::fs::read_to_string(&out_no_rg).unwrap();
    assert!(no_rg.contains("a\t"));
    assert!(no_rg.contains("b\t"));
    assert!(no_rg.contains("c\t"));
    assert!(!no_rg.contains("\nd\t") && !no_rg.starts_with("d\t"));
}

#[test]
fn view_d_and_dash_cap_d_filter_by_aux_tag() {
    let tmp = tmp_dir("aux-tag-filter");
    let input = tmp.join("in.sam");
    let tag_file = tmp.join("nm.txt");
    let out_present = tmp.join("present.sam");
    let out_value = tmp.join("value.sam");
    let out_file = tmp.join("file.sam");

    std::fs::write(
        &input,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tNM:i:0\n",
            "b\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\tNM:i:1\n",
            "c\t0\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\t****\tNM:i:2\n",
            "d\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t&&&&\n",
        ),
    )
    .unwrap();
    std::fs::write(&tag_file, "0\n2\n").unwrap();

    assert_eq!(
        run(&[
            "-d",
            "NM",
            "-o",
            out_present.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    let present = std::fs::read_to_string(&out_present).unwrap();
    assert!(present.contains("a\t"));
    assert!(present.contains("b\t"));
    assert!(present.contains("c\t"));
    assert!(!present.contains("\nd\t") && !present.starts_with("d\t"));

    assert_eq!(
        run(&[
            "-d",
            "NM:1",
            "-o",
            out_value.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    let value_filter = std::fs::read_to_string(&out_value).unwrap();
    assert!(value_filter.contains("b\t"));
    assert!(!value_filter.contains("\na\t") && !value_filter.starts_with("a\t"));
    assert!(!value_filter.contains("\nc\t") && !value_filter.starts_with("c\t"));
    assert!(!value_filter.contains("\nd\t") && !value_filter.starts_with("d\t"));

    let tag_file_arg = format!("NM:{}", tag_file.display());
    assert_eq!(
        run(&[
            "-D",
            &tag_file_arg,
            "-o",
            out_file.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    let file_filter = std::fs::read_to_string(&out_file).unwrap();
    assert!(file_filter.contains("a\t"));
    assert!(!file_filter.contains("\nb\t") && !file_filter.starts_with("b\t"));
    assert!(file_filter.contains("c\t"));
    assert!(!file_filter.contains("\nd\t") && !file_filter.starts_with("d\t"));
}

#[test]
fn view_filter_unmapped_succeeds() {
    let p = fixtures_dir().join("view.001.sam");
    // -f 4 (require unmapped) should succeed; output may be empty or contain
    // unmapped records.
    assert_eq!(run(&["-c", "-f", "4", p.to_str().unwrap()]), 0);
}

#[test]
fn view_unselected_sam_output_splits_filter_results() {
    let tmp = tmp_dir("unselected-sam");
    let input = tmp.join("input.sam");
    let selected = tmp.join("selected.sam");
    let unselected = tmp.join("unselected.sam");

    std::fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\trg:Z:a\nr2\t0\tref\t2\t0\t2M\t*\t0\t0\tTG\t##\trg:Z:b\nr3\t4\t*\t0\t0\t*\t*\t0\t0\tNN\t!!\trg:Z:c\n",
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-h",
            "-q",
            "10",
            "-U",
            unselected.to_str().unwrap(),
            "-o",
            selected.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );

    let selected_text = std::fs::read_to_string(selected).unwrap();
    let unselected_text = std::fs::read_to_string(unselected).unwrap();
    assert!(selected_text.starts_with("@HD\t"));
    assert!(unselected_text.starts_with("@HD\t"));
    assert!(selected_text.contains("\nr1\t"));
    assert!(!selected_text.contains("\nr2\t"));
    assert!(!selected_text.contains("\nr3\t"));
    assert!(unselected_text.contains("\nr2\t"));
    assert!(unselected_text.contains("\nr3\t"));
    assert!(!unselected_text.contains("\nr1\t"));
}

#[test]
fn view_unselected_sam_output_splits_expr_results() {
    let tmp = tmp_dir("unselected-sam-expr");
    let input = tmp.join("input.sam");
    let selected = tmp.join("selected.sam");
    let unselected = tmp.join("unselected.sam");

    std::fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\trg:Z:a\nr2\t0\tref\t2\t0\t2M\t*\t0\t0\tTG\t##\trg:Z:b\n",
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-h",
            "-e",
            "mapq >= 10",
            "-U",
            unselected.to_str().unwrap(),
            "-o",
            selected.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );

    let selected_text = std::fs::read_to_string(selected).unwrap();
    let unselected_text = std::fs::read_to_string(unselected).unwrap();
    assert!(selected_text.starts_with("@HD\t"));
    assert!(unselected_text.starts_with("@HD\t"));
    assert!(selected_text.contains("\nr1\t"));
    assert!(!selected_text.contains("\nr2\t"));
    assert!(unselected_text.contains("\nr2\t"));
    assert!(!unselected_text.contains("\nr1\t"));
}

#[test]
fn view_unmap_unselected_sam_output_keeps_failed_records() {
    let tmp = tmp_dir("unmap-unselected-sam");
    let input = tmp.join("input.sam");
    let output = tmp.join("output.sam");

    std::fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\nr2\t0\tref\t2\t0\t2M\t*\t0\t99\tTG\t##\n",
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-q",
            "10",
            "-p",
            "-o",
            output.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );

    let text = std::fs::read_to_string(output).unwrap();
    let records: Vec<&str> = text.lines().collect();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].split('\t').next(), Some("r1"));
    assert_eq!(records[0].split('\t').nth(1), Some("0"));
    assert_eq!(records[1].split('\t').next(), Some("r2"));
    let failed_fields: Vec<&str> = records[1].split('\t').collect();
    assert_eq!(failed_fields[1], "4");
    assert_eq!(failed_fields[4], "0");
    assert_eq!(failed_fields[5], "*");
    assert_eq!(failed_fields[8], "0");
}

#[test]
fn view_unmap_unselected_sam_output_keeps_expr_failed_records() {
    let tmp = tmp_dir("unmap-unselected-sam-expr");
    let input = tmp.join("input.sam");
    let output = tmp.join("output.sam");

    std::fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\nr2\t0\tref\t2\t0\t2M\t*\t0\t99\tTG\t##\n",
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-e",
            "mapq >= 10",
            "-p",
            "-o",
            output.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );

    let text = std::fs::read_to_string(output).unwrap();
    let records: Vec<&str> = text.lines().collect();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].split('\t').next(), Some("r1"));
    assert_eq!(records[0].split('\t').nth(1), Some("0"));
    assert_eq!(records[1].split('\t').next(), Some("r2"));
    let failed_fields: Vec<&str> = records[1].split('\t').collect();
    assert_eq!(failed_fields[1], "4");
    assert_eq!(failed_fields[4], "0");
    assert_eq!(failed_fields[5], "*");
    assert_eq!(failed_fields[8], "0");
}

#[test]
fn view_bam_sam_unselected_output_splits_expr_results() {
    let tmp = tmp_dir("bam-unselected-sam-expr");
    let input = tmp.join("input.sam");
    let bam = tmp.join("input.bam");
    let selected = tmp.join("selected.sam");
    let unselected = tmp.join("unselected.sam");

    std::fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\nr2\t0\tref\t2\t0\t2M\t*\t0\t0\tTG\t##\n",
    )
    .unwrap();
    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &input,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-h",
            "-e",
            "mapq >= 10",
            "-U",
            unselected.to_str().unwrap(),
            "-o",
            selected.to_str().unwrap(),
            bam.to_str().unwrap(),
        ]),
        0
    );

    let selected_text = std::fs::read_to_string(selected).unwrap();
    let unselected_text = std::fs::read_to_string(unselected).unwrap();
    assert!(selected_text.contains("\nr1\t"));
    assert!(!selected_text.contains("\nr2\t"));
    assert!(unselected_text.contains("\nr2\t"));
    assert!(!unselected_text.contains("\nr1\t"));
}

#[test]
fn view_sam_to_bam_output_honors_mapq_filter() {
    let tmp = tmp_dir("sam-bam-line-filter");
    let out = tmp.join("out.bam");
    let sam = fixtures_dir().join("view.001.sam");

    assert_eq!(
        run(&[
            "-b",
            "-q",
            "10",
            "-o",
            out.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 10));
}

#[test]
fn view_sam_to_bam_output_strips_tags() {
    let tmp = tmp_dir("sam-bam-tag-filter");
    let input = tmp.join("input.sam");
    let out = tmp.join("out.bam");

    std::fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\trg:Z:a\taa:i:1\n",
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-b",
            "-x",
            "rg",
            "-o",
            out.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    assert!(text.contains("\tr1\t") || text.contains("\nr1\t"));
    assert!(!text.contains("\trg:Z:a"));
    assert!(text.contains("\taa:i:1"));
}

#[test]
fn view_sam_to_bam_output_strips_tags_with_expr() {
    let tmp = tmp_dir("sam-bam-tag-filter-expr");
    let input = tmp.join("input.sam");
    let out = tmp.join("out.bam");

    std::fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:ref1\tLN:56\nr1\t0\tref1\t1\t20\t2M\t*\t0\t0\tAC\t!!\trg:Z:a\taa:i:1\nr2\t0\tref1\t2\t0\t2M\t*\t0\t0\tTG\t##\trg:Z:b\taa:i:2\n",
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-b",
            "-e",
            "mapq >= 10",
            "-x",
            "rg",
            "-o",
            out.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    assert!(text.contains("\tr1\t") || text.contains("\nr1\t"));
    assert!(!text.contains("\nr2\t"));
    assert!(!text.contains("\trg:Z:a"));
    assert!(!text.contains("\trg:Z:b"));
    assert!(text.contains("\taa:i:1"));
}

#[test]
fn view_sanitize_mutates_sam_output_records() {
    let tmp = tmp_dir("sanitize-sam");
    let input = tmp.join("input.sam");
    let out = tmp.join("out.sam");

    std::fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t4\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\tNM:i:1\tMD:Z:1A\n",
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-z",
            "all",
            "-o",
            out.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    let record = text.lines().find(|line| !line.starts_with('@')).unwrap();
    let fields: Vec<&str> = record.split('\t').collect();
    assert_eq!(fields[0], "r1");
    assert_eq!(fields[4], "0");
    assert_eq!(fields[5], "*");
    assert!(!record.contains("\tNM:i:"));
    assert!(!record.contains("\tMD:Z:"));
}

#[test]
fn view_sam_to_cram_output_honors_mapq_filter() {
    let tmp = tmp_dir("sam-cram-line-filter");
    let out = tmp.join("out.cram");
    let sam = fixtures_dir().join("view.001.sam");
    let reference = fixtures_dir().join("view.001.fa");

    assert_eq!(
        run(&[
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "-q",
            "10",
            "-o",
            out.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            &out, &reference, None,
        )
        .unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 10));
}

#[test]
fn view_sam_to_cram_output_strips_tags_with_expr() {
    let tmp = tmp_dir("sam-cram-tag-filter-expr");
    let input = tmp.join("input.sam");
    let out = tmp.join("out.cram");
    let reference = fixtures_dir().join("view.001.fa");

    std::fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:ref1\tLN:56\nr1\t0\tref1\t1\t20\t2M\t*\t0\t0\tAC\t!!\trg:Z:a\taa:i:1\nr2\t0\tref1\t2\t0\t2M\t*\t0\t0\tTG\t##\trg:Z:b\taa:i:2\n",
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "-e",
            "mapq >= 10",
            "-x",
            "rg",
            "-o",
            out.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            &out, &reference, None,
        )
        .unwrap();
    assert!(text.contains("\tr1\t") || text.contains("\nr1\t"));
    assert!(!text.contains("\nr2\t"));
    assert!(!text.contains("\trg:Z:a"));
    assert!(!text.contains("\trg:Z:b"));
    assert!(text.contains("\taa:i:1"));
}

#[test]
fn view_bam_to_bam_output_honors_mapq_filter() {
    let tmp = tmp_dir("bam-bam-line-filter");
    let bam = tmp.join("input.bam");
    let out = tmp.join("out.bam");
    let sam = fixtures_dir().join("view.001.sam");

    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-b",
            "-q",
            "10",
            "-o",
            out.to_str().unwrap(),
            bam.to_str().unwrap(),
        ]),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 10));
    assert!(
        htslib_rs::alignment_compat::count_bam_records_from_path(&out).unwrap()
            < htslib_rs::alignment_compat::count_bam_records_from_path(&bam).unwrap()
    );
}

#[test]
fn view_bam_expr_count_succeeds() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-expr-count");
    let bam = tmp.join("input.bam");
    let out = tmp.join("count.txt");
    let sam = fixtures_dir().join("view.001.sam");

    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-c",
                "-e",
                "mapq >= 1",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
            ],
        ))),
        0
    );

    let count = std::fs::read_to_string(out)
        .unwrap()
        .trim()
        .parse::<usize>()
        .unwrap();
    let total = htslib_rs::alignment_compat::count_bam_records_from_path(&bam).unwrap();
    assert!(count > 0);
    assert!(count < total);
}

#[test]
fn view_bam_region_expr_count_succeeds() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-region-expr-count");
    let out = tmp.join("count.txt");
    let bam = htslib_fixtures_dir().join("range.bam");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-c",
                "-e",
                "mapq >= 0",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
                "CHROMOSOME_II:2980-2980",
            ],
        ))),
        0
    );

    let count = std::fs::read_to_string(out)
        .unwrap()
        .trim()
        .parse::<usize>()
        .unwrap();
    assert!(count > 0);
}

#[test]
fn view_dash_l_filters_by_read_group_library() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("view-l-library");
    let sam = tmp.join("in.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "@RG\tID:rg1\tLB:libA\tSM:s1\n",
            "@RG\tID:rg2\tLB:libB\tSM:s2\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:rg1\n",
            "r2\t0\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:rg2\n",
            "r3\t0\tchr1\t9\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:rg1\n",
            "r4\t0\tchr1\t13\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    let count = |args: &[&str]| -> usize {
        let out = tmp.join(format!("c{}.txt", args.join("_")));
        let mut full = vec!["-c", "-o", out.to_str().unwrap()];
        full.extend_from_slice(args);
        full.push(sam.to_str().unwrap());
        assert_eq!(run(&full), 0);
        std::fs::read_to_string(&out)
            .unwrap()
            .trim()
            .parse()
            .unwrap()
    };

    // libA → rg1 only (r1, r3); libB → rg2 (r2); unknown library → none.
    // The RG-less r4 is always excluded under -l.
    assert_eq!(count(&["-l", "libA"]), 2);
    assert_eq!(count(&["-l", "libB"]), 1);
    assert_eq!(count(&["-l", "libZ"]), 0);
    // --library long form behaves identically.
    assert_eq!(count(&["--library", "libA"]), 2);
}

#[test]
fn view_dash_cap_x_accepts_legacy_custom_index_synopsis() {
    // `view -X in.bam in.bam.bai region` — the second positional is the
    // explicit index path. We accept it as a no-op (our region query
    // builds/finds the index itself) and the region still applies, so
    // the count must match the non-`-X` invocation.
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("view-x-index");
    let bam = htslib_fixtures_dir().join("range.bam");
    let bai = htslib_fixtures_dir().join("range.bam.bai");
    let region = "CHROMOSOME_II:2980-2980";

    let plain = tmp.join("plain.txt");
    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-c",
                "-o",
                plain.to_str().unwrap(),
                bam.to_str().unwrap(),
                region,
            ],
        ))),
        0
    );

    let xed = tmp.join("xed.txt");
    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-X",
                "-c",
                "-o",
                xed.to_str().unwrap(),
                bam.to_str().unwrap(),
                bai.to_str().unwrap(),
                region,
            ],
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(&plain).unwrap(),
        std::fs::read_to_string(&xed).unwrap()
    );
}

#[test]
fn view_bam_region_expr_sam_output_succeeds() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-region-expr-sam");
    let out = tmp.join("view.sam");
    let bam = htslib_fixtures_dir().join("range.bam");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-h",
                "-e",
                "mapq >= 20",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
                "CHROMOSOME_II:2980-2980",
            ],
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("@HD\t"));
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 20));
}

#[test]
fn view_bam_expr_sam_output_succeeds() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-expr-sam");
    let bam = tmp.join("input.bam");
    let out = tmp.join("view.sam");
    let sam = fixtures_dir().join("view.001.sam");

    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-h",
                "-e",
                "mapq >= 10",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
            ],
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("@HD\t"));
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 10));
}

#[test]
fn view_bam_expr_bam_output_succeeds() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-expr-bam");
    let bam = tmp.join("input.bam");
    let out = tmp.join("view.bam");
    let sam = fixtures_dir().join("view.001.sam");

    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-b",
                "-e",
                "mapq >= 10",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
            ],
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 10));
    assert!(
        htslib_rs::alignment_compat::count_bam_records_from_path(&out).unwrap()
            < htslib_rs::alignment_compat::count_bam_records_from_path(&bam).unwrap()
    );
}

#[test]
fn view_bam_expr_cram_output_succeeds() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-expr-cram");
    let bam = tmp.join("input.bam");
    let out = tmp.join("view.cram");
    let sam = fixtures_dir().join("view.001.sam");
    let reference = fixtures_dir().join("view.001.fa");

    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-C",
                "-T",
                reference.to_str().unwrap(),
                "-e",
                "mapq >= 10",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
            ],
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            &out, &reference, None,
        )
        .unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 10));
}

#[test]
fn view_cram_expr_count_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-expr-count");
    let out = tmp.join("count.txt");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-c",
                "-e",
                "mapq >= 0",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    let count = std::fs::read_to_string(out)
        .unwrap()
        .trim()
        .parse::<usize>()
        .unwrap();
    assert!(count > 0);
}

#[test]
fn view_cram_region_expr_count_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-region-expr-count");
    let out = tmp.join("count.txt");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-c",
                "-e",
                "mapq >= 20",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
                "CHROMOSOME_II:2980-2980",
            ],
        ))),
        0
    );

    let count = std::fs::read_to_string(out)
        .unwrap()
        .trim()
        .parse::<usize>()
        .unwrap();
    assert!(count > 0);
}

#[test]
fn view_cram_region_expr_sam_output_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-region-expr-sam");
    let out = tmp.join("view.sam");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-h",
                "-e",
                "mapq >= 20",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
                "CHROMOSOME_II:2980-2980",
            ],
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("@HD\t"));
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 20));
}

#[test]
fn view_cram_expr_sam_output_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-expr-sam");
    let out = tmp.join("view.sam");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-h",
                "-e",
                "mapq >= 20",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("@HD\t"));
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 20));
}

#[test]
fn view_cram_sam_unselected_output_splits_expr_results() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-unselected-sam-expr");
    let selected = tmp.join("selected.sam");
    let unselected = tmp.join("unselected.sam");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-h",
                "-e",
                "mapq >= 20",
                "-U",
                unselected.to_str().unwrap(),
                "-o",
                selected.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    let selected_text = std::fs::read_to_string(selected).unwrap();
    let unselected_text = std::fs::read_to_string(unselected).unwrap();
    assert!(selected_text.starts_with("@HD\t"));
    assert!(unselected_text.starts_with("@HD\t"));

    let selected_records: Vec<&str> = selected_text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .collect();
    let unselected_records: Vec<&str> = unselected_text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .collect();
    assert!(!selected_records.is_empty());
    assert!(!unselected_records.is_empty());
    assert!(
        selected_records
            .iter()
            .all(|line| line.split('\t').nth(4).unwrap().parse::<u8>().unwrap() >= 20)
    );
    assert!(
        unselected_records.iter().any(|line| line
            .split('\t')
            .nth(4)
            .unwrap()
            .parse::<u8>()
            .unwrap()
            < 20)
    );
}

#[test]
fn view_cram_expr_bam_output_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-expr-bam");
    let out = tmp.join("view.bam");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-b",
                "-e",
                "mapq >= 20",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 20));
}

#[test]
fn view_bam_region_expr_bam_output_succeeds() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-region-expr-bam");
    let out = tmp.join("region.bam");
    let bam = htslib_fixtures_dir().join("range.bam");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-b",
                "-e",
                "mapq >= 20",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
                "CHROMOSOME_II:2980-2980",
            ],
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 20));
}

#[test]
fn view_bam_region_bam_output_honors_mapq_filter() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-region-bam-line-filter");
    let out = tmp.join("region.bam");
    let bam = htslib_fixtures_dir().join("range.bam");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-b",
                "-q",
                "20",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
                "CHROMOSOME_II:2980-2980",
            ],
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    let records: Vec<&str> = text.lines().filter(|line| !line.starts_with('@')).collect();
    assert!(!records.is_empty());
    for record in records {
        let fields: Vec<&str> = record.split('\t').collect();
        assert_eq!(fields[2], "CHROMOSOME_II");
        assert!(fields[4].parse::<u8>().unwrap() >= 20);
    }
}

#[test]
fn view_cram_region_expr_bam_output_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-region-expr-bam");
    let out = tmp.join("region.bam");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-b",
                "-e",
                "mapq >= 20",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
                "CHROMOSOME_II:2980-2980",
            ],
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 20));
}

#[test]
fn view_cram_expr_cram_output_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-expr-cram");
    let out = tmp.join("view.cram");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-C",
                "-e",
                "mapq >= 20",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            &out, &reference, None,
        )
        .unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 20));
}

#[test]
fn view_bam_region_expr_cram_output_succeeds() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-region-expr-cram");
    let out = tmp.join("region.cram");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let bam = fixtures.join("range.bam");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-C",
                "-T",
                reference.to_str().unwrap(),
                "-e",
                "mapq >= 20",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
                "CHROMOSOME_II:2980-2980",
            ],
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            &out, &reference, None,
        )
        .unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 20));
}

#[test]
fn view_cram_region_expr_cram_output_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-region-expr-cram");
    let out = tmp.join("region.cram");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-C",
                "-e",
                "mapq >= 20",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
                "CHROMOSOME_II:2980-2980",
            ],
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            &out, &reference, None,
        )
        .unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 20));
}

#[test]
fn view_cram_count_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-count");
    let out = tmp.join("count.txt");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-c",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    let count = std::fs::read_to_string(out)
        .unwrap()
        .trim()
        .parse::<usize>()
        .unwrap();
    assert!(count > 0);
}

#[test]
fn view_cram_filtered_count_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-filtered-count");
    let out = tmp.join("count.txt");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-c",
                "-q",
                "20",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    let count = std::fs::read_to_string(out)
        .unwrap()
        .trim()
        .parse::<usize>()
        .unwrap();
    assert!(count > 0);
}

#[test]
fn view_cram_sam_output_uses_local_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-sam");
    let out = tmp.join("view.sam");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-h",
                "-T",
                reference.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("@HD\t"));
    assert!(text.contains("HS18_09653"));
}

#[test]
fn view_cram_to_bam_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-bam");
    let out = tmp.join("view.bam");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-b",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    assert!(out.metadata().unwrap().len() > 0);
    assert!(htslib_rs::alignment_compat::count_bam_records_from_path(&out).unwrap() > 0);
}

#[test]
fn view_cram_to_bam_output_honors_mapq_filter() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-bam-line-filter");
    let out = tmp.join("view.bam");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-b",
                "-q",
                "20",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 20));
}

#[test]
fn view_cram_region_to_bam_uses_local_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-region-bam");
    let out = tmp.join("region.bam");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-b",
                "-T",
                reference.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
                "CHROMOSOME_II:2980-2980",
            ],
        ))),
        0
    );

    assert!(out.metadata().unwrap().len() > 0);
    assert!(htslib_rs::alignment_compat::count_bam_records_from_path(&out).unwrap() > 0);
}

#[test]
fn view_cram_to_cram_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-cram");
    let out = tmp.join("view.cram");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-C",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    assert!(out.metadata().unwrap().len() > 0);
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(&out).unwrap();
    assert!(!header.reference_sequences().is_empty());
}

#[test]
fn view_cram_to_cram_output_honors_mapq_filter() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-cram-line-filter");
    let out = tmp.join("view.cram");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-C",
                "-q",
                "20",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            &out, &reference, None,
        )
        .unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 20));
}

#[test]
fn view_bam_to_cram_uses_local_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-cram");
    let bam = tmp.join("input.bam");
    let out = tmp.join("view.cram");
    let sam = fixtures_dir().join("view.001.sam");
    let reference = fixtures_dir().join("view.001.fa");

    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-C",
                "-T",
                reference.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
            ],
        ))),
        0
    );

    assert!(out.metadata().unwrap().len() > 0);
    let header = htslib_rs::alignment_compat::read_cram_header_from_path(&out).unwrap();
    assert!(!header.reference_sequences().is_empty());
}

#[test]
fn view_bam_to_cram_output_honors_mapq_filter() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-cram-line-filter");
    let bam = tmp.join("input.bam");
    let out = tmp.join("view.cram");
    let sam = fixtures_dir().join("view.001.sam");
    let reference = fixtures_dir().join("view.001.fa");

    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-C",
                "-T",
                reference.to_str().unwrap(),
                "-q",
                "10",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
            ],
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            &out, &reference, None,
        )
        .unwrap();
    let mapqs: Vec<u8> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert!(!mapqs.is_empty());
    assert!(mapqs.iter().all(|mapq| *mapq >= 10));
    assert!(mapqs.len() < htslib_rs::alignment_compat::count_bam_records_from_path(&bam).unwrap());
}

#[test]
fn view_bam_region_to_cram_uses_local_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-region-cram");
    let bam = tmp.join("input.bam");
    let bai = tmp.join("input.bam.bai");
    let out = tmp.join("region.cram");
    let sam = fixtures_dir().join("view.001.sam");
    let reference = fixtures_dir().join("view.001.fa");

    htslib_rs::alignment_compat::write_bam_from_sam_path_with_bai(&sam, &bam, &bai).unwrap();

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-C",
                "-T",
                reference.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
                "ref1:1-5",
            ],
        ))),
        0
    );

    assert!(out.metadata().unwrap().len() > 0);
    assert!(
        htslib_rs::alignment_compat::benchmark_cram_view_from_path_with_reference(&out, &reference)
            .unwrap()
            > 0
    );
}

#[test]
fn view_cram_region_to_cram_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("cram-region-cram");
    let out = tmp.join("region.cram");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-C",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
                "CHROMOSOME_II:2980-2980",
            ],
        ))),
        0
    );

    assert!(out.metadata().unwrap().len() > 0);
    assert!(
        htslib_rs::alignment_compat::benchmark_cram_view_from_path_with_reference(&out, &reference)
            .unwrap()
            > 0
    );
}

#[test]
fn view_cram_records_without_reference_fail_cleanly() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let cram = htslib_fixtures_dir().join("range.cram");

    assert_ne!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &["view", "-c", cram.to_str().unwrap()],
        ))),
        0
    );
}

#[test]
fn view_header_only_appends_pg_line_by_default() {
    let tmp = tmp_dir("view-pg-header-only");
    let sam = tmp.join("input.sam");
    let out = tmp.join("hdr.sam");
    std::fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\n@PG\tID:upstream\tPN:upstream\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\n",
    )
    .unwrap();

    assert_eq!(
        run(&["-H", "-o", out.to_str().unwrap(), sam.to_str().unwrap()]),
        0
    );
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("@PG\tID:upstream"));
    assert!(text.contains("PN:samtools"));
    assert!(text.contains("PP:upstream"));
}

#[test]
fn view_no_pg_suppresses_pg_line_in_header_only() {
    let tmp = tmp_dir("view-no-pg-header-only");
    let sam = tmp.join("input.sam");
    let out = tmp.join("hdr.sam");
    std::fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\n@PG\tID:upstream\tPN:upstream\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\n",
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-H",
            "--no-PG",
            "-o",
            out.to_str().unwrap(),
            sam.to_str().unwrap()
        ]),
        0
    );
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("@PG\tID:upstream"));
    assert!(!text.contains("PN:samtools"));
}

#[test]
fn view_h_flag_sam_output_appends_pg_line() {
    let tmp = tmp_dir("view-pg-sam-output");
    let sam = tmp.join("input.sam");
    let out = tmp.join("full.sam");
    std::fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\n",
    )
    .unwrap();

    assert_eq!(
        run(&["-h", "-o", out.to_str().unwrap(), sam.to_str().unwrap()]),
        0
    );
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("PN:samtools"));
    assert!(text.contains("\nr1\t"));
}

#[test]
fn view_p_unmap_unselected_routes_into_bam_output_for_sam_input() {
    let tmp = tmp_dir("view-unmap-bam");
    let input = tmp.join("input.sam");
    let bam_out = tmp.join("out.bam");
    let sam_out = tmp.join("out.sam");

    std::fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\nr2\t0\tref\t2\t0\t2M\t*\t0\t0\tTG\t##\n",
    )
    .unwrap();

    // -p with binary output: records below MAPQ threshold should be marked
    // unmapped in the resulting BAM. Round-trip BAM -> SAM via view.
    assert_eq!(
        run(&[
            "-b",
            "-p",
            "-h",
            "-q",
            "10",
            "-o",
            bam_out.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    assert_eq!(
        run(&[
            "-h",
            "-o",
            sam_out.to_str().unwrap(),
            bam_out.to_str().unwrap()
        ]),
        0
    );
    let text = std::fs::read_to_string(&sam_out).unwrap();
    let bodies: Vec<&str> = text.lines().filter(|l| !l.starts_with('@')).collect();
    assert_eq!(bodies.len(), 2);
    // r2 (mapq=0) gets unmapped flag set (4) and CIGAR/MAPQ cleared.
    let r1 = bodies.iter().find(|l| l.starts_with("r1\t")).unwrap();
    let r2 = bodies.iter().find(|l| l.starts_with("r2\t")).unwrap();
    let r1_flags: u32 = r1.split('\t').nth(1).unwrap().parse().unwrap();
    let r2_flags: u32 = r2.split('\t').nth(1).unwrap().parse().unwrap();
    assert_eq!(r1_flags & 4, 0);
    assert_eq!(r2_flags & 4, 4);
}

#[test]
fn view_u_unselected_routes_into_bam_output_for_sam_input() {
    let tmp = tmp_dir("view-unselected-bam");
    let input = tmp.join("input.sam");
    let sel_bam = tmp.join("sel.bam");
    let unsel_bam = tmp.join("unsel.bam");
    let sel_sam = tmp.join("sel.sam");
    let unsel_sam = tmp.join("unsel.sam");

    std::fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\nr2\t0\tref\t2\t0\t2M\t*\t0\t0\tTG\t##\n",
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-b",
            "-h",
            "-q",
            "10",
            "-U",
            unsel_bam.to_str().unwrap(),
            "-o",
            sel_bam.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    assert_eq!(
        run(&[
            "-h",
            "-o",
            sel_sam.to_str().unwrap(),
            sel_bam.to_str().unwrap()
        ]),
        0
    );
    assert_eq!(
        run(&[
            "-h",
            "-o",
            unsel_sam.to_str().unwrap(),
            unsel_bam.to_str().unwrap()
        ]),
        0
    );

    let sel_text = std::fs::read_to_string(&sel_sam).unwrap();
    let unsel_text = std::fs::read_to_string(&unsel_sam).unwrap();
    assert!(sel_text.contains("\nr1\t"));
    assert!(!sel_text.contains("\nr2\t"));
    assert!(unsel_text.contains("\nr2\t"));
    assert!(!unsel_text.contains("\nr1\t"));
}

/// `view -b in.sam` must embed the samtools `@PG` in the **binary**
/// BAM header (completed library batch #4 / TODO.md). Round-trip via `view -h`:
/// without `--no-PG` the BAM carries a `PN:samtools` `@PG`; with
/// `--no-PG` it does not. Records are unaffected.
#[test]
fn view_b_embeds_pg_in_binary_bam_header() {
    let tmp = tmp_dir("view-b-pg");
    let sam = tmp.join("in.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:10\nr1\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();

    // -b without --no-PG: BAM header gains a samtools @PG.
    let bam = tmp.join("out.bam");
    assert_eq!(
        run(&["-b", sam.to_str().unwrap(), "-o", bam.to_str().unwrap()]),
        0
    );
    let hdr = tmp.join("hdr.sam");
    assert_eq!(
        run(&[
            "-H",
            "--no-PG",
            bam.to_str().unwrap(),
            "-o",
            hdr.to_str().unwrap(),
        ]),
        0
    );
    let hdr_text = std::fs::read_to_string(&hdr).unwrap();
    assert!(
        hdr_text
            .lines()
            .any(|l| l.starts_with("@PG") && l.contains("PN:samtools")),
        "binary BAM header must carry the samtools @PG, got:\n{hdr_text}"
    );

    // -b --no-PG: no @PG in the BAM header.
    let bam2 = tmp.join("out.nopg.bam");
    assert_eq!(
        run(&[
            "-b",
            "--no-PG",
            sam.to_str().unwrap(),
            "-o",
            bam2.to_str().unwrap(),
        ]),
        0
    );
    let hdr2 = tmp.join("hdr2.sam");
    assert_eq!(
        run(&[
            "-H",
            "--no-PG",
            bam2.to_str().unwrap(),
            "-o",
            hdr2.to_str().unwrap(),
        ]),
        0
    );
    assert!(
        !std::fs::read_to_string(&hdr2).unwrap().contains("@PG"),
        "view -b --no-PG must not add a @PG"
    );

    // Records survive the SAM->BAM @PG-injection round-trip.
    let recs = tmp.join("recs.sam");
    assert_eq!(
        run(&[
            "--no-PG",
            bam.to_str().unwrap(),
            "-o",
            recs.to_str().unwrap(),
        ]),
        0
    );
    let rt = std::fs::read_to_string(&recs).unwrap();
    assert!(
        rt.lines()
            .any(|l| l.starts_with("r1\t") && l.contains("ACGT"))
    );

    // BAM-input -> BAM-output also injects the @PG into the binary
    // header (records streamed unchanged), suppressed by --no-PG.
    let nopg_bam = tmp.join("base.nopg.bam");
    assert_eq!(
        run(&[
            "-b",
            "--no-PG",
            sam.to_str().unwrap(),
            "-o",
            nopg_bam.to_str().unwrap(),
        ]),
        0
    );
    let bb = tmp.join("bam2bam.bam");
    assert_eq!(
        run(&["-b", nopg_bam.to_str().unwrap(), "-o", bb.to_str().unwrap()]),
        0
    );
    let bbh = tmp.join("bb.sam");
    assert_eq!(
        run(&[
            "-H",
            "--no-PG",
            bb.to_str().unwrap(),
            "-o",
            bbh.to_str().unwrap(),
        ]),
        0
    );
    assert!(
        std::fs::read_to_string(&bbh)
            .unwrap()
            .lines()
            .any(|l| l.starts_with("@PG") && l.contains("PN:samtools")),
        "BAM->BAM must inject the samtools @PG into the binary header"
    );
    let bbrecs = tmp.join("bb.recs.sam");
    assert_eq!(
        run(&[
            "--no-PG",
            bb.to_str().unwrap(),
            "-o",
            bbrecs.to_str().unwrap(),
        ]),
        0
    );
    assert!(
        std::fs::read_to_string(&bbrecs)
            .unwrap()
            .lines()
            .any(|l| l.starts_with("r1\t") && l.contains("ACGT")),
        "BAM->BAM @PG injection must preserve records"
    );
}

/// The SAM-text-intermediate binary paths (here: `view -b -z`
/// sanitizer on BAM input) also inject the samtools `@PG` into the
/// binary BAM header (completed library batch #4).
#[test]
fn view_b_sanitizer_bam_path_embeds_pg() {
    let tmp = tmp_dir("view-b-sanitize-pg");
    let sam = tmp.join("in.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:10\nr1\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();
    let bam = tmp.join("in.bam");
    assert_eq!(
        run(&[
            "-b",
            "--no-PG",
            sam.to_str().unwrap(),
            "-o",
            bam.to_str().unwrap(),
        ]),
        0
    );
    // BAM input + sanitizer -> BAM output (SAM-text intermediate path).
    let out = tmp.join("out.bam");
    assert_eq!(
        run(&[
            "-b",
            "-z",
            "all",
            bam.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ]),
        0
    );
    let hdr = tmp.join("h.sam");
    assert_eq!(
        run(&[
            "-H",
            "--no-PG",
            out.to_str().unwrap(),
            "-o",
            hdr.to_str().unwrap(),
        ]),
        0
    );
    assert!(
        std::fs::read_to_string(&hdr)
            .unwrap()
            .lines()
            .any(|l| l.starts_with("@PG") && l.contains("PN:samtools")),
        "sanitizer BAM->BAM path must inject the samtools @PG"
    );
}

/// BAM-input filter and region binary copies inject the samtools
/// `@PG` (routed via the SAM-text path when `@PG` is wanted),
/// suppressed by `--no-PG` which keeps the fast binary copy
/// (completed library batch #4).
#[test]
fn view_b_bam_filter_and_region_paths_embed_pg() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("view-b-filter-region-pg");
    let sam = tmp.join("in.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:20\n\
         r1\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
         r2\t0\tc1\t5\t10\t4M\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();
    let bam = tmp.join("in.bam");
    assert_eq!(
        run(&[
            "-b",
            "--no-PG",
            sam.to_str().unwrap(),
            "-o",
            bam.to_str().unwrap(),
        ]),
        0
    );
    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &["index", bam.to_str().unwrap()]
        ))),
        0
    );

    let has_pg = |p: &std::path::Path| -> bool {
        let h = tmp.join("h.sam");
        assert_eq!(
            run(&[
                "-H",
                "--no-PG",
                p.to_str().unwrap(),
                "-o",
                h.to_str().unwrap()
            ]),
            0
        );
        std::fs::read_to_string(&h)
            .unwrap()
            .lines()
            .any(|l| l.starts_with("@PG") && l.contains("PN:samtools"))
    };

    // BAM filter -> BAM
    let f = tmp.join("f.bam");
    assert_eq!(
        run(&[
            "-b",
            "-e",
            "mapq>=20",
            bam.to_str().unwrap(),
            "-o",
            f.to_str().unwrap(),
        ]),
        0
    );
    assert!(has_pg(&f), "BAM filter -> BAM must inject @PG");

    // BAM region -> BAM
    let r = tmp.join("r.bam");
    assert_eq!(
        run(&[
            "-b",
            bam.to_str().unwrap(),
            "c1:1-3",
            "-o",
            r.to_str().unwrap(),
        ]),
        0
    );
    assert!(has_pg(&r), "BAM region -> BAM must inject @PG");

    // --no-PG keeps it absent.
    let fn_ = tmp.join("fn.bam");
    assert_eq!(
        run(&[
            "-b",
            "--no-PG",
            "-e",
            "mapq>=20",
            bam.to_str().unwrap(),
            "-o",
            fn_.to_str().unwrap(),
        ]),
        0
    );
    assert!(!has_pg(&fn_), "view -b --no-PG must not add a @PG");
}

/// `view -C` (SAM->CRAM) likewise embeds the samtools `@PG` in the
/// CRAM header unless `--no-PG` (completed library batch #4).
#[test]
fn view_c_embeds_pg_in_binary_cram_header() {
    // Uses the same known-good fixtures as the other SAM->CRAM tests.
    let sam = fixtures_dir().join("view.001.sam");
    let reference = fixtures_dir().join("view.001.fa");
    let r = reference.to_str().unwrap();
    let tmp = tmp_dir("view-c-pg");

    // Read the CRAM header back via the htslib-rs helper (the same
    // approach the other SAM->CRAM tests use), so the assertion is on
    // the bytes actually stored in the binary CRAM header.
    let cram_header = |p: &std::path::Path| -> String {
        let text =
            htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
                p,
                &reference,
                Some(0),
            )
            .unwrap();
        text.lines()
            .take_while(|l| l.starts_with('@'))
            .map(|l| format!("{l}\n"))
            .collect()
    };

    let cram = tmp.join("out.cram");
    assert_eq!(
        run(&[
            "-C",
            "-T",
            r,
            sam.to_str().unwrap(),
            "-o",
            cram.to_str().unwrap(),
        ]),
        0
    );
    assert!(
        cram_header(&cram)
            .lines()
            .any(|l| l.starts_with("@PG") && l.contains("PN:samtools")),
        "binary CRAM header must carry the samtools @PG"
    );

    let cram2 = tmp.join("out.nopg.cram");
    assert_eq!(
        run(&[
            "-C",
            "--no-PG",
            "-T",
            r,
            sam.to_str().unwrap(),
            "-o",
            cram2.to_str().unwrap(),
        ]),
        0
    );
    assert!(
        !cram_header(&cram2).contains("PN:samtools"),
        "view -C --no-PG must not add a samtools @PG"
    );
}
