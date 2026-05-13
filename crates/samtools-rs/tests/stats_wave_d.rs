//! Integration smoke-tests for the Wave D commands (`depth`, `coverage`,
//! `bedcov`) and the indexed-region modes of `view`.

use std::ffi::OsString;
use std::io::Read;
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering};

use flate2::read::MultiGzDecoder;
use samtools_rs::commands::{bedcov, coverage, depth, index, view};
use samtools_rs::native;

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

fn exit_to_u8(code: ExitCode) -> u8 {
    format!("{:?}", code)
        .chars()
        .filter(|c| c.is_ascii_digit())
        .collect::<String>()
        .parse()
        .unwrap_or(255)
}

fn tmp_dir(name: &str) -> PathBuf {
    static NEXT_TMP_ID: AtomicUsize = AtomicUsize::new(0);

    let id = NEXT_TMP_ID.fetch_add(1, Ordering::Relaxed);
    let p = std::env::temp_dir().join(format!(
        "samtools-rs-waveD-{}-{}-{}",
        name,
        std::process::id(),
        id
    ));
    std::fs::create_dir_all(&p).unwrap();
    p
}

fn argv(name: &str, rest: &[&str]) -> Vec<OsString> {
    std::iter::once(OsString::from(name))
        .chain(rest.iter().map(OsString::from))
        .collect()
}

fn indexed_bam() -> PathBuf {
    let tmp = tmp_dir("idx");
    let bam = tmp.join("in.bam");
    std::fs::copy(fixtures_dir().join("checksum").join("chk1.bam"), &bam).unwrap();
    assert_eq!(
        exit_to_u8(index::main(&argv("index", &[bam.to_str().unwrap()]))),
        0
    );
    bam
}

#[test]
fn depth_outputs_some_positions() {
    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[indexed_bam().to_str().unwrap()]
        ))),
        0
    );
}

#[test]
fn depth_region_outputs_only_requested_interval() {
    let p = indexed_bam();
    let tmp = tmp_dir("depth-region");
    let out = tmp.join("depth.tsv");
    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-r",
                "17:1-10000",
                "-o",
                out.to_str().unwrap(),
                p.to_str().unwrap()
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(!text.is_empty());
    for line in text.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields[0], "17");
        let pos = fields[1].parse::<usize>().unwrap();
        assert!((1..=10000).contains(&pos));
    }
}

#[test]
fn depth_bed_outputs_only_bed_interval() {
    let p = indexed_bam();
    let tmp = tmp_dir("depth-bed");
    let bed = tmp.join("regions.bed");
    let out = tmp.join("depth.tsv");
    std::fs::write(&bed, "17\t0\t10000\n").unwrap();
    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-b",
                bed.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                p.to_str().unwrap()
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(!text.is_empty());
    for line in text.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields[0], "17");
        let pos = fields[1].parse::<usize>().unwrap();
        assert!((1..=10000).contains(&pos));
    }
}

#[test]
fn coverage_outputs_per_ref() {
    let p = indexed_bam();
    assert_eq!(
        exit_to_u8(coverage::main(&argv("coverage", &[p.to_str().unwrap()]))),
        0
    );
}

#[test]
fn coverage_region_outputs_requested_interval() {
    let p = indexed_bam();
    let tmp = tmp_dir("coverage-region");
    let out = tmp.join("coverage.tsv");
    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-r",
                "17:1-10000",
                "-o",
                out.to_str().unwrap(),
                p.to_str().unwrap()
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    let rows: Vec<_> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(rows.len(), 1);
    assert!(rows[0].starts_with("17\t1\t10000\t"));
}

#[test]
fn bedcov_with_one_region() {
    let bam = indexed_bam();
    let tmp = tmp_dir("bedcov");
    let bed = tmp.join("r.bed");
    std::fs::write(&bed, "17\t0\t10000\n").unwrap();
    assert_eq!(
        exit_to_u8(bedcov::main(&argv(
            "bedcov",
            &[bed.to_str().unwrap(), bam.to_str().unwrap()]
        ))),
        0
    );
}

#[test]
fn bedcov_count_column_succeeds() {
    let bam = indexed_bam();
    let tmp = tmp_dir("bedcov-count");
    let bed = tmp.join("r.bed");
    std::fs::write(&bed, "17\t0\t10000\n").unwrap();
    assert_eq!(
        exit_to_u8(bedcov::main(&argv(
            "bedcov",
            &["-c", bed.to_str().unwrap(), bam.to_str().unwrap()]
        ))),
        0
    );
}

#[test]
fn bedcov_depth_column_succeeds() {
    let bam = indexed_bam();
    let tmp = tmp_dir("bedcov-depth");
    let bed = tmp.join("r.bed");
    std::fs::write(&bed, "17\t0\t10000\n").unwrap();
    assert_eq!(
        exit_to_u8(bedcov::main(&argv(
            "bedcov",
            &["-d", "2", bed.to_str().unwrap(), bam.to_str().unwrap()]
        ))),
        0
    );
}

#[test]
fn view_region_query_succeeds() {
    let p = indexed_bam();
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &["-c", p.to_str().unwrap(), "17:1-1000000"]
        ))),
        0
    );
}

#[test]
fn view_region_bam_output_succeeds() {
    let p = indexed_bam();
    let tmp = tmp_dir("view-bam-region");
    let out = tmp.join("slice.bam");
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &[
                "-P",
                "-b",
                p.to_str().unwrap(),
                "17:1-10000",
                "-o",
                out.to_str().unwrap()
            ]
        ))),
        0
    );

    assert!(out.exists());
    assert!(htslib_rs::alignment_compat::count_bam_records_from_path(&out).unwrap() > 0);
}

#[test]
fn view_bed_bam_output_succeeds() {
    let p = indexed_bam();
    let tmp = tmp_dir("view-bam-bed");
    let bed = tmp.join("regions.bed");
    let out = tmp.join("slice.bam");
    std::fs::write(&bed, "17\t0\t10000\n").unwrap();

    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &[
                "-P",
                "-b",
                p.to_str().unwrap(),
                "-L",
                bed.to_str().unwrap(),
                "-o",
                out.to_str().unwrap()
            ]
        ))),
        0
    );

    assert!(out.exists());
    assert!(htslib_rs::alignment_compat::count_bam_records_from_path(&out).unwrap() > 0);
}

#[test]
fn native_view_region_writes_bam_output() {
    let p = indexed_bam();
    let tmp = tmp_dir("native-view-region");
    let out = tmp.join("slice.bam");

    native::view_region(&p, "17:1-10000", &out, Some(1), None).unwrap();

    assert!(out.exists());
    assert!(htslib_rs::alignment_compat::count_bam_records_from_path(&out).unwrap() > 0);
}

#[test]
fn native_index_writes_default_and_explicit_bai() {
    let tmp = tmp_dir("native-index");
    let bam = tmp.join("in.bam");
    std::fs::copy(fixtures_dir().join("checksum").join("chk1.bam"), &bam).unwrap();

    let default_bai = native::index(&bam, Option::<&PathBuf>::None, Some(1)).unwrap();
    assert_eq!(default_bai, tmp.join("in.bam.bai"));
    assert!(default_bai.exists());

    let explicit_bai = tmp.join("explicit.bai");
    let written = native::index(&bam, Some(&explicit_bai), Some(1)).unwrap();
    assert_eq!(written, explicit_bai);
    assert!(written.exists());
}

#[test]
fn native_bam_to_fastq_pair_writes_gzip_outputs() {
    let tmp = tmp_dir("native-fastq");
    let bam = tmp.join("in.bam");
    std::fs::copy(fixtures_dir().join("checksum").join("chk1.bam"), &bam).unwrap();

    let r1 = tmp.join("sample_R1.fastq.gz");
    let r2 = tmp.join("sample_R2.fastq.gz");
    let other = tmp.join("sample_other.fastq.gz");
    let singleton = tmp.join("sample_single.fastq.gz");

    native::bam_to_fastq_pair(
        &bam,
        &r1,
        &r2,
        Some(other.as_path()),
        Some(singleton.as_path()),
        true,
        Some(1),
    )
    .unwrap();

    let r1_text = read_gzip_text(&r1);
    let r2_text = read_gzip_text(&r2);
    let other_text = read_gzip_text(&other);
    let singleton_text = read_gzip_text(&singleton);

    assert!(r1_text.starts_with('@'));
    assert!(r2_text.starts_with('@'));
    assert!(other_text.is_empty());
    assert!(singleton_text.starts_with('@') || singleton_text.is_empty());
}

#[test]
fn native_fastq_routes_paired_singleton_and_other_outputs() {
    let tmp = tmp_dir("native-fastq-routing");
    let sam = tmp.join("in.sam");
    let bam = tmp.join("in.bam");
    let r1 = tmp.join("r1.fastq");
    let r2 = tmp.join("r2.fastq");
    let other = tmp.join("other.fastq");
    let singleton = tmp.join("singleton.fastq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "pair\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\n",
            "pair\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\n",
            "solo\t4\t*\t0\t0\t*\t*\t0\t0\tNNNN\t$$$$\n",
        ),
    )
    .unwrap();
    let bam_file = std::fs::File::create(&bam).unwrap();
    htslib_rs::alignment_compat::write_bam_from_sam_path(&sam, bam_file).unwrap();

    native::bam_to_fastq_pair(
        &bam,
        &r1,
        &r2,
        Some(other.as_path()),
        Some(singleton.as_path()),
        true,
        Some(1),
    )
    .unwrap();

    assert_eq!(
        std::fs::read_to_string(r1).unwrap(),
        "@pair\nACGT\n+\n!!!!\n"
    );
    assert_eq!(
        std::fs::read_to_string(r2).unwrap(),
        "@pair\nTGCA\n+\n####\n"
    );
    assert_eq!(std::fs::read_to_string(other).unwrap(), "");
    assert_eq!(
        std::fs::read_to_string(singleton).unwrap(),
        "@solo\nNNNN\n+\n$$$$\n"
    );
}

#[test]
fn native_depth_and_summary_return_structured_values() {
    let p = indexed_bam();

    let depths = native::depth(&p, "17:1-10000", true, Some(1)).unwrap();
    let summary = native::depth_summary(&p, "17:1-10000", true, Some(1)).unwrap();

    assert!(!depths.is_empty());
    assert_eq!(depths[0].reference_name, "17");
    assert_eq!(depths.len(), 10000);
    assert_eq!(depths.first().unwrap().position, 1);
    assert_eq!(depths.last().unwrap().position, 10000);
    assert!(summary.max >= summary.min);
    assert!(summary.mean >= 0.0);
    assert!(summary.median >= 0.0);
}

#[test]
fn native_sort_merge_and_quickcheck_wrappers_work() {
    let tmp = tmp_dir("native-p1");
    let bam_a = tmp.join("a.bam");
    let bam_b = tmp.join("b.bam");
    std::fs::copy(fixtures_dir().join("checksum").join("chk1.bam"), &bam_a).unwrap();
    std::fs::copy(fixtures_dir().join("checksum").join("chk1.bam"), &bam_b).unwrap();

    let sorted = tmp.join("sorted.bam");
    native::sort(&bam_a, &sorted, true, Some(1)).unwrap();
    assert!(htslib_rs::alignment_compat::count_bam_records_from_path(&sorted).unwrap() > 0);

    let merged = tmp.join("merged.bam");
    native::merge(&merged, &[&bam_a, &bam_b], true, Some(1)).unwrap();
    assert!(htslib_rs::alignment_compat::count_bam_records_from_path(&merged).unwrap() > 0);
    assert_eq!(
        htslib_rs::alignment_compat::read_bam_header_from_path(&bam_a)
            .unwrap()
            .reference_sequences()
            .len(),
        htslib_rs::alignment_compat::read_bam_header_from_path(&merged)
            .unwrap()
            .reference_sequences()
            .len()
    );

    let coord_sorted = tmp.join("coord-sorted.bam");
    native::sort(&bam_a, &coord_sorted, false, Some(1)).unwrap();
    let bai = native::index(&coord_sorted, Option::<&PathBuf>::None, Some(1)).unwrap();
    assert!(bai.exists());

    native::quickcheck(&bam_a, true).unwrap();
    native::quickcheck(
        fixtures_dir().join("quickcheck/6.quickcheck.cram21.ok.cram"),
        true,
    )
    .unwrap();
}

#[test]
fn native_cram_region_view_outputs_compatible_bam() {
    let tmp = tmp_dir("native-cram-view");
    let cram = htslib_fixtures_dir().join("range.cram");
    let reference = htslib_fixtures_dir().join("ce.fa");
    let out = tmp.join("slice.bam");
    let sorted = tmp.join("slice.sorted.bam");
    let merged = tmp.join("slice.merged.bam");
    let r1 = tmp.join("slice.r1.fastq");
    let r2 = tmp.join("slice.r2.fastq");

    assert!(native::view_region(&cram, "CHROMOSOME_II:2980-2980", &out, Some(1), None).is_err());
    native::view_region(
        &cram,
        "CHROMOSOME_II:2980-2980",
        &out,
        Some(1),
        Some(reference.as_path()),
    )
    .unwrap();

    assert!(htslib_rs::alignment_compat::count_bam_records_from_path(&out).unwrap() > 0);
    assert!(
        native::index(&out, Option::<&PathBuf>::None, Some(1))
            .unwrap()
            .exists()
    );
    native::sort(&out, &sorted, false, Some(1)).unwrap();
    native::merge(&merged, &[&out, &sorted], true, Some(1)).unwrap();
    native::bam_to_fastq_pair(&out, &r1, &r2, None, None, true, Some(1)).unwrap();
    assert!(r1.exists());
    assert!(r2.exists());
}

#[test]
fn native_extract_unmapped_pairs_supports_bam_and_cram_reference_path() {
    let tmp = tmp_dir("native-unmapped");
    let sam = tmp.join("in.sam");
    let bam = tmp.join("in.bam");
    let reference = tmp.join("ref.fa");
    let cram = tmp.join("in.cram");
    let bam_out = tmp.join("unmapped.bam");
    let cram_out = tmp.join("unmapped-from-cram.bam");
    std::fs::write(&reference, ">chr1\nACGT\n").unwrap();
    std::fs::write(reference.with_extension("fa.fai"), "chr1\t4\t6\t4\t5\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:4\n",
            "pair\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!\n",
            "pair\t141\t*\t0\t0\t*\t*\t0\t0\tTGCA\t####\n",
            "mapped\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_cram_from_sam_path_with_reference(
        &sam,
        &reference,
        std::fs::File::create(&cram).unwrap(),
    )
    .unwrap();

    native::extract_unmapped_pairs(&bam, &bam_out, 12, Some(1), None).unwrap();
    native::extract_unmapped_pairs(&cram, &cram_out, 12, Some(1), Some(reference.as_path()))
        .unwrap();

    assert_eq!(
        htslib_rs::alignment_compat::count_bam_records_from_path(&bam_out).unwrap(),
        2
    );
    assert_eq!(
        htslib_rs::alignment_compat::count_bam_records_from_path(&cram_out).unwrap(),
        2
    );
}

fn read_gzip_text(path: &PathBuf) -> String {
    let file = std::fs::File::open(path).unwrap();
    let mut decoder = MultiGzDecoder::new(file);
    let mut text = String::new();
    decoder.read_to_string(&mut text).unwrap();
    text
}

#[test]
fn view_filter_expr_for_sam() {
    let p = fixtures_dir().join("dat").join("mpileup.1.sam");
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &["-c", "-e", "pos<1000", p.to_str().unwrap()]
        ))),
        0
    );
}
