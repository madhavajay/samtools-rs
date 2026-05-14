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
