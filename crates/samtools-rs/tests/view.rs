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
fn view_save_counts_reports_processed_accepted_and_rejected() {
    let tmp = tmp_dir("save-counts");
    let input = tmp.join("in.sam");
    let output = tmp.join("out.sam");
    let counts = tmp.join("counts.json");

    std::fs::write(
        &input,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "b\t128\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "c\t0\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\t****\n",
            "d\t129\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t&&&&\n",
        ),
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-f",
            "128",
            "--save-counts",
            counts.to_str().unwrap(),
            "--no-PG",
            "-o",
            output.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );

    assert_eq!(
        std::fs::read_to_string(&counts).unwrap(),
        concat!(
            "{\n",
            "    \"records_processed\" : 4,\n",
            "    \"records_filter_accepted\" : 2,\n",
            "    \"records_filter_rejected\" : 2\n",
            "}\n",
        )
    );

    let selected = std::fs::read_to_string(&output).unwrap();
    assert!(!selected.contains("\na\t") && !selected.starts_with("a\t"));
    assert!(selected.contains("\nb\t") || selected.starts_with("b\t"));
    assert!(!selected.contains("\nc\t") && !selected.starts_with("c\t"));
    assert!(selected.contains("\nd\t") || selected.starts_with("d\t"));
}

#[test]
fn view_count_save_counts_cram_without_reference_uses_summary_path() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("save-counts-cram-no-ref");
    let output = tmp.join("count.txt");
    let counts = tmp.join("counts.json");
    let cram = htslib_fixtures_dir().join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-c",
                "-f",
                "64",
                "--save-counts",
                counts.to_str().unwrap(),
                "-o",
                output.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    assert_eq!(std::fs::read_to_string(&output).unwrap(), "55\n");
    assert_eq!(
        std::fs::read_to_string(&counts).unwrap(),
        concat!(
            "{\n",
            "    \"records_processed\" : 112,\n",
            "    \"records_filter_accepted\" : 55,\n",
            "    \"records_filter_rejected\" : 57\n",
            "}\n"
        )
    );
}

#[test]
fn view_count_save_counts_no_reference_cram_supports_summary_expr_filters() {
    let tmp = tmp_dir("save-counts-cram-expr-no-ref");
    let cram = fixtures_dir().join("test_input_1_a.cram");

    for (expr, expected_count, accepted, rejected) in [
        ("mapq>=20", "14\n", 14, 1),
        ("flag.proper_pair", "4\n", 4, 11),
    ] {
        let label = expr.replace(|c: char| !c.is_ascii_alphanumeric(), "_");
        let output = tmp.join(format!("{label}.txt"));
        let counts = tmp.join(format!("{label}.json"));

        assert_eq!(
            run(&[
                "-c",
                "-e",
                expr,
                "--save-counts",
                counts.to_str().unwrap(),
                "-o",
                output.to_str().unwrap(),
                cram.to_str().unwrap(),
            ]),
            0,
            "expression {expr}"
        );

        assert_eq!(std::fs::read_to_string(&output).unwrap(), expected_count);
        assert_eq!(
            std::fs::read_to_string(&counts).unwrap(),
            format!(
                "{{\n    \"records_processed\" : 15,\n    \"records_filter_accepted\" : {accepted},\n    \"records_filter_rejected\" : {rejected}\n}}\n"
            )
        );
    }
}

#[test]
fn view_count_save_counts_no_reference_cram_supports_aux_rg_library_filters() {
    let tmp = tmp_dir("save-counts-cram-aux-no-ref");
    let reference = tmp.join("ref.fa");
    let sam = tmp.join("input.sam");
    let cram = tmp.join("input.cram");

    std::fs::write(&reference, ">chr1\nACGTACGTACGT\n").unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:12\n",
            "@RG\tID:rg1\tLB:lib1\n",
            "@RG\tID:rg2\tLB:lib2\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:rg1\tNM:i:0\tXX:Z:keep\n",
            "b\t0\tchr1\t2\t10\t4M\t*\t0\t0\tCGTA\t####\tRG:Z:rg2\tNM:i:1\tXX:Z:drop\n",
            "c\t4\t*\t0\t0\t*\t*\t0\t0\tNN\t!!\tYY:i:7\n",
        ),
    )
    .unwrap();

    assert_eq!(
        run(&[
            "--no-PG",
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "-o",
            cram.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );

    for (label, args, expected_count, accepted, rejected) in [
        (
            "read-group-with-no-rg-pass-through",
            vec!["-r", "rg1"],
            "2\n",
            2,
            1,
        ),
        (
            "read-group-exclude-no-rg",
            vec!["-r", "rg1", "-n"],
            "1\n",
            1,
            2,
        ),
        ("library", vec!["-l", "lib2"], "1\n", 1, 2),
        ("aux-string", vec!["-d", "XX:keep"], "1\n", 1, 2),
        ("aux-int", vec!["-d", "YY:7"], "1\n", 1, 2),
    ] {
        let output = tmp.join(format!("{label}.txt"));
        let counts = tmp.join(format!("{label}.json"));
        let mut view_args = vec![
            "-c",
            "--save-counts",
            counts.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ];
        view_args.extend(args);
        view_args.push(cram.to_str().unwrap());

        assert_eq!(run(&view_args), 0, "case {label}");
        assert_eq!(std::fs::read_to_string(&output).unwrap(), expected_count);
        assert_eq!(
            std::fs::read_to_string(&counts).unwrap(),
            format!(
                "{{\n    \"records_processed\" : 3,\n    \"records_filter_accepted\" : {accepted},\n    \"records_filter_rejected\" : {rejected}\n}}\n"
            ),
            "case {label}"
        );
    }
}

#[test]
fn view_min_query_length_filters_by_cigar_query_consuming_ops() {
    let tmp = tmp_dir("min-query-len");
    let input = tmp.join("in.sam");
    let output = tmp.join("out.sam");

    std::fs::write(
        &input,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "m5\t0\tchr1\t1\t60\t5M\t*\t0\t0\tAAAAA\t!!!!!\n",
            "ins\t0\tchr1\t1\t60\t3M2I\t*\t0\t0\tAAAAC\t#####\n",
            "del\t0\tchr1\t1\t60\t4M3D\t*\t0\t0\tAAAA\t****\n",
            "soft\t0\tchr1\t1\t60\t2S4M\t*\t0\t0\tTTAAAA\t&&&&&&\n",
            "hard\t0\tchr1\t1\t60\t2H4M\t*\t0\t0\tAAAA\t++++\n",
        ),
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-m",
            "5",
            "--no-PG",
            "-o",
            output.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );

    let selected = std::fs::read_to_string(&output).unwrap();
    assert!(selected.contains("\nm5\t") || selected.starts_with("m5\t"));
    assert!(selected.contains("\nins\t") || selected.starts_with("ins\t"));
    assert!(!selected.contains("\ndel\t") && !selected.starts_with("del\t"));
    assert!(selected.contains("\nsoft\t") || selected.starts_with("soft\t"));
    assert!(!selected.contains("\nhard\t") && !selected.starts_with("hard\t"));
}

#[test]
fn view_filtered_sam_input_normalizes_float_aux_arrays() {
    let tmp = tmp_dir("filtered-float-aux");
    let input = tmp.join("in.sam");
    let output = tmp.join("out.sam");

    std::fs::write(
        &input,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "drop\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "keep\t128\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\t",
            "bg:B:f,2.71828,0.0000000000000000000000000000000006626,2997900000\n",
        ),
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-f",
            "128",
            "--no-PG",
            "-o",
            output.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );

    let selected = std::fs::read_to_string(&output).unwrap();
    assert!(selected.contains("keep\t128\t"));
    assert!(selected.contains("bg:B:f,2.71828,6.626e-34,2.9979e+09"));
    assert!(!selected.contains("2997900000"));
    assert!(!selected.contains("\ndrop\t") && !selected.starts_with("drop\t"));
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
    let out_rg_and_no_rg = tmp.join("rg_and_no_rg.sam");
    let out_header = tmp.join("header.sam");

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
    assert!(one.contains("\nd\t") || one.starts_with("d\t"));

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
    assert!(many.contains("\nd\t") || many.starts_with("d\t"));

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

    assert_eq!(
        run(&[
            "-r",
            "grp2",
            "-n",
            "-o",
            out_rg_and_no_rg.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    let rg_and_no_rg = std::fs::read_to_string(&out_rg_and_no_rg).unwrap();
    assert!(rg_and_no_rg.contains("\nb\t") || rg_and_no_rg.starts_with("b\t"));
    assert!(!rg_and_no_rg.contains("\nd\t") && !rg_and_no_rg.starts_with("d\t"));

    assert_eq!(
        run(&[
            "-h",
            "-r",
            "grp2",
            "-o",
            out_header.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    let with_header = std::fs::read_to_string(&out_header).unwrap();
    assert!(!with_header.contains("@RG\tID:grp1\t"));
    assert!(with_header.contains("@RG\tID:grp2\t"));
    assert!(!with_header.contains("@RG\tID:grp3\t"));
    assert!(with_header.contains("\nb\t"));
    assert!(with_header.contains("\nd\t"));
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
fn view_expr_supports_htslib_flag_names_for_sam_and_bam() {
    let tmp = tmp_dir("expr-flag-names");
    let input = tmp.join("input.sam");
    let bam = tmp.join("input.bam");
    let sam_out = tmp.join("sam.out");
    let bam_out = tmp.join("bam.out");

    std::fs::write(
        &input,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:ref\tLN:100\n",
            "proper\t99\tref\t1\t30\t4M\t=\t20\t19\tACGT\t!!!!\n",
            "plain\t0\tref\t2\t30\t4M\t*\t0\t0\tTGCA\t####\n",
        ),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &input,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-e",
            "flag.proper_pair",
            "--no-PG",
            "-o",
            sam_out.to_str().unwrap(),
            input.to_str().unwrap(),
        ]),
        0
    );
    assert_eq!(
        run(&[
            "-e",
            "flag.proper_pair",
            "--no-PG",
            "-o",
            bam_out.to_str().unwrap(),
            bam.to_str().unwrap(),
        ]),
        0
    );

    let sam_text = std::fs::read_to_string(sam_out).unwrap();
    let bam_text = std::fs::read_to_string(bam_out).unwrap();
    assert!(sam_text.contains("proper\t"));
    assert!(!sam_text.contains("plain\t"));
    assert!(bam_text.contains("proper\t"));
    assert!(!bam_text.contains("plain\t"));
}

#[test]
fn view_expr_supports_cigar_derived_symbols_on_sam_line_path() {
    let tmp = tmp_dir("expr-cigar-symbols");
    let input = tmp.join("input.sam");
    let rlen_out = tmp.join("rlen.sam");
    let sclen_out = tmp.join("sclen.sam");
    let hclen_out = tmp.join("hclen.sam");
    let endpos_out = tmp.join("endpos.sam");

    std::fs::write(
        &input,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:ref\tLN:100\n",
            "plain\t0\tref\t2\t30\t4M\t*\t0\t0\tTGCA\t####\n",
            "soft\t0\tref\t5\t30\t2S4M\t*\t0\t0\tAATGCA\t!!!!!!\n",
            "hard\t0\tref\t10\t30\t3H4M\t*\t0\t0\tACGT\t!!!!\n",
            "long\t0\tref\t20\t30\t8M2D\t*\t0\t0\tAAAAAAAA\t!!!!!!!!\n",
        ),
    )
    .unwrap();

    for (expr, out) in [
        ("rlen>=10", &rlen_out),
        ("sclen>0", &sclen_out),
        ("hclen>0", &hclen_out),
        ("endpos>=29", &endpos_out),
    ] {
        assert_eq!(
            run(&[
                "-e",
                expr,
                "--no-PG",
                "-o",
                out.to_str().unwrap(),
                input.to_str().unwrap(),
            ]),
            0,
            "expression {expr}"
        );
    }

    assert!(
        std::fs::read_to_string(&rlen_out)
            .unwrap()
            .contains("long\t")
    );
    assert!(
        std::fs::read_to_string(&sclen_out)
            .unwrap()
            .contains("soft\t")
    );
    assert!(
        std::fs::read_to_string(&hclen_out)
            .unwrap()
            .contains("hard\t")
    );
    assert!(
        std::fs::read_to_string(&endpos_out)
            .unwrap()
            .contains("long\t")
    );
}

#[test]
fn view_expr_supports_mate_symbols_on_sam_line_path() {
    let tmp = tmp_dir("expr-mate-symbols");
    let input = tmp.join("input.sam");
    let mpos_out = tmp.join("mpos.sam");
    let pnext_out = tmp.join("pnext.sam");
    let tlen_out = tmp.join("tlen.sam");
    let rnext_out = tmp.join("rnext.sam");
    let mrname_out = tmp.join("mrname.sam");
    let ncigar_out = tmp.join("ncigar.sam");

    std::fs::write(
        &input,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:ref\tLN:100\n",
            "pair\t99\tref\t2\t30\t4M\t=\t20\t22\tTGCA\t####\n",
            "pair\t147\tref\t20\t30\t4M\t=\t2\t-22\tACGT\t!!!!\n",
            "single\t0\tref\t50\t30\t4M\t*\t0\t0\tNNNN\t!!!!\n",
        ),
    )
    .unwrap();

    for (expr, out) in [
        ("mpos>0", &mpos_out),
        ("pnext>0", &pnext_out),
        ("tlen!=0", &tlen_out),
        ("rnext==\"ref\"", &rnext_out),
        ("mrname==\"ref\"", &mrname_out),
        ("ncigar==1", &ncigar_out),
    ] {
        assert_eq!(
            run(&[
                "-e",
                expr,
                "--no-PG",
                "-o",
                out.to_str().unwrap(),
                input.to_str().unwrap(),
            ]),
            0,
            "expression {expr}"
        );
    }

    for out in [&mpos_out, &pnext_out, &tlen_out, &rnext_out, &mrname_out] {
        let text = std::fs::read_to_string(out).unwrap();
        assert!(text.contains("pair\t"));
        assert!(!text.contains("single\t"));
    }

    let text = std::fs::read_to_string(&ncigar_out).unwrap();
    assert!(text.contains("pair\t"));
    assert!(text.contains("single\t"));
}

#[test]
fn view_count_mode_supports_paired_and_mate_expressions() {
    // The SAM-output expression path already covers mate/paired
    // symbols; this locks the same symbols in `-c` count mode (and
    // BAM input), the combination the Task 3 follow-up called out.
    let tmp = tmp_dir("count-paired-mate-expr");
    let sam = tmp.join("input.sam");
    let bam = tmp.join("input.bam");

    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:200\n",
            // r1: a proper pair on chr1 (two records), TLEN +/-35.
            "r1\t99\tchr1\t10\t60\t5M\t=\t40\t35\tACGTA\t!!!!!\n",
            "r1\t147\tchr1\t40\t60\t5M\t=\t10\t-35\tTTTTT\t!!!!!\n",
            // r2: an unpaired single record.
            "r2\t0\tchr1\t80\t30\t5M\t*\t0\t0\tGGGGG\t!!!!!\n",
        ),
    )
    .unwrap();

    // Stage a BAM copy so the count path is exercised on binary input
    // too (not just the SAM line path).
    assert_eq!(
        run(&[
            "-b",
            "--no-PG",
            "-o",
            bam.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );

    for input in [&sam, &bam] {
        let path = input.to_str().unwrap();
        let count = |expr: &str| -> String {
            let out = tmp.join("c.txt");
            assert_eq!(
                run(&["-c", "-e", expr, "-o", out.to_str().unwrap(), path]),
                0,
                "count -e {expr} on {path}"
            );
            std::fs::read_to_string(&out).unwrap().trim().to_string()
        };

        assert_eq!(count("flag.proper_pair"), "2", "proper_pair count");
        assert_eq!(count("flag.paired"), "2", "paired count");
        assert_eq!(count("tlen>0"), "1", "tlen>0 count");
        assert_eq!(count("mpos==40"), "1", "mate-position count");
        assert_eq!(count("pnext==40"), "1", "pnext alias count");
        assert_eq!(count("mrname==\"chr1\""), "2", "mate-reference count");
    }
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

    let header_out = tmp.join("lib_b_header.sam");
    assert_eq!(
        run(&[
            "-h",
            "-l",
            "libB",
            "-o",
            header_out.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );
    let header_text = std::fs::read_to_string(header_out).unwrap();
    assert!(header_text.contains("@RG\tID:rg1\tLB:libA\tSM:s1\n"));
    assert!(header_text.contains("@RG\tID:rg2\tLB:libB\tSM:s2\n"));
    assert!(!header_text.contains("\nr1\t"));
    assert!(header_text.contains("\nr2\t"));
    assert!(!header_text.contains("\nr3\t"));
    assert!(!header_text.contains("\nr4\t"));
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
fn view_dash_cap_x_uses_index_at_nondefault_path() {
    // `view -X in.bam custom/dir/in.bam.bai region` — the explicit index
    // lives in a directory with no default-location index next to the BAM.
    // The region query must still succeed using the provided index, and
    // must not create an index beside the source BAM.
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("view-x-nondefault");
    let data_dir = tmp.join("data");
    let idx_dir = tmp.join("idx");
    std::fs::create_dir_all(&data_dir).unwrap();
    std::fs::create_dir_all(&idx_dir).unwrap();

    let bam = data_dir.join("range.bam");
    std::fs::copy(htslib_fixtures_dir().join("range.bam"), &bam).unwrap();
    let custom_bai = idx_dir.join("range.bam.bai");
    std::fs::copy(htslib_fixtures_dir().join("range.bam.bai"), &custom_bai).unwrap();
    // Ensure there is no default-location index beside the BAM.
    assert!(!data_dir.join("range.bam.bai").exists());

    let region = "CHROMOSOME_II:2980-2980";
    let expected = tmp.join("expected.txt");
    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "-c",
                "-o",
                expected.to_str().unwrap(),
                htslib_fixtures_dir().join("range.bam").to_str().unwrap(),
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
                custom_bai.to_str().unwrap(),
                region,
            ],
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(&expected).unwrap(),
        std::fs::read_to_string(&xed).unwrap()
    );
    // The source BAM directory must stay clean (no leaked default index).
    assert!(!data_dir.join("range.bam.bai").exists());
    assert!(!data_dir.join("range.bam.csi").exists());
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
fn view_bam_to_cram_uses_header_ur_reference_without_dash_t() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bam-cram-header-ur");
    let reference = tmp.join("ref.fa");
    let sam = tmp.join("input.sam");
    let bam = tmp.join("input.bam");
    let out = tmp.join("view.cram");

    std::fs::write(&reference, ">ref1\nACGTACGT\n").unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();
    std::fs::write(
        &sam,
        format!(
            "@HD\tVN:1.6\n@SQ\tSN:ref1\tLN:8\tUR:file://{}\n\
             r1\t0\tref1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            reference.display()
        ),
    )
    .unwrap();

    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();

    assert_eq!(
        run(&["-C", "-o", out.to_str().unwrap(), bam.to_str().unwrap()]),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            &out, &reference, None,
        )
        .unwrap();
    assert!(text.contains("\nr1\t"));

    assert_eq!(run(&["-c", out.to_str().unwrap()]), 0);

    let roundtrip_bam = tmp.join("roundtrip.bam");
    assert_eq!(
        run(&[
            "-b",
            "-o",
            roundtrip_bam.to_str().unwrap(),
            out.to_str().unwrap(),
        ]),
        0
    );
    let roundtrip_text = htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(
        &roundtrip_bam,
        None,
    )
    .unwrap();
    assert!(roundtrip_text.contains("\nr1\t"));
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
fn view_cram_count_without_reference_uses_summary_path() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let cram = htslib_fixtures_dir().join("range.cram");

    assert_eq!(
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
fn view_glued_ho_writes_header_to_output_path() {
    let tmp = tmp_dir("view-glued-ho");
    let sam = tmp.join("input.sam");
    let out = tmp.join("out.bam");
    std::fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\n",
    )
    .unwrap();

    assert_eq!(
        run(&["-ho", out.to_str().unwrap(), sam.to_str().unwrap()]),
        0
    );
    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    assert!(text.starts_with("@HD\t"));
    assert!(text.contains("\nr1\t"));
}

#[test]
fn view_bam_output_resolves_reference_aliases() {
    let tmp = tmp_dir("view-reference-alias");
    let sam = tmp.join("input.sam");
    let out = tmp.join("out.bam");
    std::fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:r3\tLN:50\tAN:ref3\nr1\t0\tref3\t1\t30\t1M\t*\t0\t0\tA\t!\n",
    )
    .unwrap();

    assert_eq!(
        run(&[
            "-ho",
            out.to_str().unwrap(),
            "--no-PG",
            sam.to_str().unwrap()
        ]),
        0
    );
    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    assert!(text.contains("\n@SQ\tSN:r3\tLN:50\tAN:ref3\n"));
    assert!(text.contains("\nr1\t0\tr3\t1\t30\t1M\t*\t0\t0\tA\t!"));
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

#[test]
fn view_bam_input_unselected_routes_into_bam_output() {
    let tmp = tmp_dir("view-bam-unselected-bam");
    let sam = tmp.join("input.sam");
    let input_bam = tmp.join("input.bam");
    let sel_bam = tmp.join("selected.bam");
    let unsel_bam = tmp.join("unselected.bam");

    std::fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\nr2\t0\tref\t2\t0\t2M\t*\t0\t0\tTG\t##\nr3\t4\t*\t0\t0\t*\t*\t0\t0\tNN\t!!\n",
    )
    .unwrap();
    assert_eq!(
        run(&[
            "--no-PG",
            "-b",
            "-o",
            input_bam.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );

    assert_eq!(
        run(&[
            "--no-PG",
            "-b",
            "-q",
            "10",
            "-U",
            unsel_bam.to_str().unwrap(),
            "-o",
            sel_bam.to_str().unwrap(),
            input_bam.to_str().unwrap(),
        ]),
        0
    );

    let selected = htslib_rs::alignment_compat::summarize_bam_records_from_path(&sel_bam).unwrap();
    let unselected =
        htslib_rs::alignment_compat::summarize_bam_records_from_path(&unsel_bam).unwrap();
    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].name_bytes(), Some(&b"r1"[..]));
    assert_eq!(unselected.len(), 2);
    assert_eq!(unselected[0].name_bytes(), Some(&b"r2"[..]));
    assert_eq!(unselected[1].name_bytes(), Some(&b"r3"[..]));
}

#[test]
fn view_bam_input_unmap_unselected_routes_into_bam_output() {
    let tmp = tmp_dir("view-bam-unmap-bam");
    let sam = tmp.join("input.sam");
    let input_bam = tmp.join("input.bam");
    let out_bam = tmp.join("out.bam");

    std::fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\nr2\t0\tref\t2\t0\t2M\t*\t0\t0\tTG\t##\n",
    )
    .unwrap();
    assert_eq!(
        run(&[
            "--no-PG",
            "-b",
            "-o",
            input_bam.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );

    assert_eq!(
        run(&[
            "--no-PG",
            "-b",
            "-p",
            "-q",
            "10",
            "-o",
            out_bam.to_str().unwrap(),
            input_bam.to_str().unwrap(),
        ]),
        0
    );

    let records = htslib_rs::alignment_compat::summarize_bam_records_from_path(&out_bam).unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].name_bytes(), Some(&b"r1"[..]));
    assert_eq!(records[0].flags_u16() & 0x4, 0);
    assert_eq!(records[1].name_bytes(), Some(&b"r2"[..]));
    assert_eq!(records[1].flags_u16() & 0x4, 0x4);
    assert_eq!(records[1].mapping_quality(), Some(0));
}

#[test]
fn view_cram_input_unselected_routes_into_bam_output_with_reference() {
    let tmp = tmp_dir("view-cram-unselected-bam");
    let reference = tmp.join("ref.fa");
    let sam = tmp.join("input.sam");
    let input_cram = tmp.join("input.cram");
    let sel_bam = tmp.join("selected.bam");
    let unsel_bam = tmp.join("unselected.bam");

    std::fs::write(&reference, ">ref\nACGTACGTACGT\n").unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();
    std::fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:12\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\nr2\t0\tref\t2\t0\t2M\t*\t0\t0\tCG\t##\n",
    )
    .unwrap();
    assert_eq!(
        run(&[
            "--no-PG",
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "-o",
            input_cram.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );

    assert_eq!(
        run(&[
            "--no-PG",
            "-b",
            "-T",
            reference.to_str().unwrap(),
            "-q",
            "10",
            "-U",
            unsel_bam.to_str().unwrap(),
            "-o",
            sel_bam.to_str().unwrap(),
            input_cram.to_str().unwrap(),
        ]),
        0
    );

    let selected = htslib_rs::alignment_compat::summarize_bam_records_from_path(&sel_bam).unwrap();
    let unselected =
        htslib_rs::alignment_compat::summarize_bam_records_from_path(&unsel_bam).unwrap();
    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].name_bytes(), Some(&b"r1"[..]));
    assert_eq!(unselected.len(), 1);
    assert_eq!(unselected[0].name_bytes(), Some(&b"r2"[..]));
}

#[test]
fn view_bam_input_unselected_routes_into_cram_output_with_reference() {
    let tmp = tmp_dir("view-bam-unselected-cram");
    let reference = tmp.join("ref.fa");
    let sam = tmp.join("input.sam");
    let input_bam = tmp.join("input.bam");
    let sel_cram = tmp.join("selected.cram");
    let unsel_cram = tmp.join("unselected.cram");

    std::fs::write(&reference, ">ref\nACGTACGTACGT\n").unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();
    std::fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:12\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\nr2\t0\tref\t2\t0\t2M\t*\t0\t0\tCG\t##\n",
    )
    .unwrap();
    assert_eq!(
        run(&[
            "--no-PG",
            "-b",
            "-o",
            input_bam.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );

    assert_eq!(
        run(&[
            "--no-PG",
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "-q",
            "10",
            "-U",
            unsel_cram.to_str().unwrap(),
            "-o",
            sel_cram.to_str().unwrap(),
            input_bam.to_str().unwrap(),
        ]),
        0
    );

    let selected = htslib_rs::alignment_compat::summarize_cram_records_from_path_with_reference(
        &sel_cram, &reference,
    )
    .unwrap();
    let unselected = htslib_rs::alignment_compat::summarize_cram_records_from_path_with_reference(
        &unsel_cram,
        &reference,
    )
    .unwrap();
    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].name_bytes(), Some(&b"r1"[..]));
    assert_eq!(unselected.len(), 1);
    assert_eq!(unselected[0].name_bytes(), Some(&b"r2"[..]));
}

#[test]
fn view_cram_input_unmap_unselected_routes_into_cram_output_with_reference() {
    let tmp = tmp_dir("view-cram-unmap-cram");
    let reference = tmp.join("ref.fa");
    let sam = tmp.join("input.sam");
    let input_cram = tmp.join("input.cram");
    let out_cram = tmp.join("out.cram");

    std::fs::write(&reference, ">ref\nACGTACGTACGT\n").unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();
    std::fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:12\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\nr2\t0\tref\t2\t0\t2M\t*\t0\t0\tCG\t##\n",
    )
    .unwrap();
    assert_eq!(
        run(&[
            "--no-PG",
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "-o",
            input_cram.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );

    assert_eq!(
        run(&[
            "--no-PG",
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "-p",
            "-q",
            "10",
            "-o",
            out_cram.to_str().unwrap(),
            input_cram.to_str().unwrap(),
        ]),
        0
    );

    let records = htslib_rs::alignment_compat::summarize_cram_records_from_path_with_reference(
        &out_cram, &reference,
    )
    .unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].name_bytes(), Some(&b"r1"[..]));
    assert_eq!(records[0].flags_u16() & 0x4, 0);
    assert_eq!(records[1].name_bytes(), Some(&b"r2"[..]));
    assert_eq!(records[1].flags_u16() & 0x4, 0x4);
    assert_eq!(records[1].mapping_quality(), Some(0));
}

#[test]
fn view_bam_input_binary_output_strips_aux_tags() {
    let tmp = tmp_dir("view-bam-strip-tags");
    let sam = tmp.join("input.sam");
    let input_bam = tmp.join("input.bam");
    let out_bam = tmp.join("out.bam");

    std::fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\tNM:i:0\tXX:Z:keep\n",
    )
    .unwrap();
    assert_eq!(
        run(&[
            "--no-PG",
            "-b",
            "-o",
            input_bam.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );

    assert_eq!(
        run(&[
            "--no-PG",
            "-b",
            "-x",
            "NM",
            "-o",
            out_bam.to_str().unwrap(),
            input_bam.to_str().unwrap(),
        ]),
        0
    );

    let records = htslib_rs::alignment_compat::summarize_bam_records_from_path(&out_bam).unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].aux_value(*b"NM"), None);
    assert_eq!(records[0].aux_value(*b"XX"), Some(&b"keep"[..]));
}

#[test]
fn view_cram_input_binary_output_keeps_only_requested_aux_tags() {
    let tmp = tmp_dir("view-cram-keep-tags");
    let reference = tmp.join("ref.fa");
    let sam = tmp.join("input.sam");
    let input_cram = tmp.join("input.cram");
    let out_cram = tmp.join("out.cram");
    let rendered = tmp.join("rendered.sam");

    std::fs::write(&reference, ">ref\nACGTACGTACGT\n").unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();
    std::fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:12\n@RG\tID:rg1\nr1\t0\tref\t1\t20\t2M\t*\t0\t0\tAC\t!!\tRG:Z:rg1\tNM:i:0\tXX:Z:drop\n",
    )
    .unwrap();
    assert_eq!(
        run(&[
            "--no-PG",
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "-o",
            input_cram.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );

    assert_eq!(
        run(&[
            "--no-PG",
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "--keep-tag",
            "RG",
            "-o",
            out_cram.to_str().unwrap(),
            input_cram.to_str().unwrap(),
        ]),
        0
    );

    let records = htslib_rs::alignment_compat::summarize_cram_records_from_path_with_reference(
        &out_cram, &reference,
    )
    .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].aux_value(*b"RG"), Some(&b"rg1"[..]));
    assert_eq!(records[0].aux_value(*b"XX"), None);

    assert_eq!(
        run(&[
            "--no-PG",
            "-T",
            reference.to_str().unwrap(),
            "-o",
            rendered.to_str().unwrap(),
            out_cram.to_str().unwrap(),
        ]),
        0
    );
    let rendered = std::fs::read_to_string(&rendered).unwrap();
    assert_eq!(rendered.matches("\tRG:Z:rg1").count(), 1);
    assert!(!rendered.contains("\tXX:Z:drop"));
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

fn cram_container_count(path: &std::path::Path) -> usize {
    let file = std::fs::File::open(path).unwrap();
    let mut reader = htslib_rs::cram::io::reader::Builder::default().build_from_reader(file);
    reader.read_header().unwrap();

    let mut containers = 0;
    let mut container = htslib_rs::cram::io::reader::Container::default();
    while reader.read_container(&mut container).unwrap() != 0 {
        containers += 1;
    }
    containers
}

#[test]
fn view_cram_seqs_per_slice_partitions_into_multiple_containers() {
    let tmp = tmp_dir("cram-seqs-per-slice");
    let sam = htslib_fixtures_dir().join("ce#1000.sam");
    let reference = htslib_fixtures_dir().join("ce.fa");

    // Default: the ~1000-record file collapses into one container.
    let default_out = tmp.join("default.cram");
    assert_eq!(
        run(&[
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "-o",
            default_out.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );
    assert_eq!(cram_container_count(&default_out), 1);

    // `-O cram,seqs_per_slice=100` must cut a new slice/container
    // every 100 records, yielding multiple containers.
    let chunked_out = tmp.join("chunked.cram");
    assert_eq!(
        run(&[
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "-O",
            "cram,seqs_per_slice=100",
            "-o",
            chunked_out.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );
    assert!(
        cram_container_count(&chunked_out) > 1,
        "seqs_per_slice=100 should produce more than one container"
    );

    // The same option via --output-fmt-option must behave identically.
    let opt_out = tmp.join("opt.cram");
    assert_eq!(
        run(&[
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "--output-fmt-option",
            "seqs_per_slice=100",
            "-o",
            opt_out.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );
    assert!(cram_container_count(&opt_out) > 1);
}

#[test]
fn view_cram_seqs_per_slice_applies_on_filtered_output() {
    // The filter CRAM-output path must also honor seqs_per_slice (the
    // Task 2 follow-up: filter/region writers, not just full-file).
    let tmp = tmp_dir("cram-seqs-per-slice-filtered");
    let sam = htslib_fixtures_dir().join("ce#1000.sam");
    let reference = htslib_fixtures_dir().join("ce.fa");

    // `-e mapq>=0` keeps every record; with seqs_per_slice=100 the
    // filtered CRAM output must still span multiple containers.
    let out = tmp.join("filtered.cram");
    assert_eq!(
        run(&[
            "-C",
            "-T",
            reference.to_str().unwrap(),
            "-e",
            "mapq>=0",
            "-O",
            "cram,seqs_per_slice=100",
            "-o",
            out.to_str().unwrap(),
            sam.to_str().unwrap(),
        ]),
        0
    );
    assert!(
        cram_container_count(&out) > 1,
        "seqs_per_slice=100 must apply on the filtered CRAM write path"
    );
}
