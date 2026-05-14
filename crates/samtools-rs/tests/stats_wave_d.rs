//! Integration smoke-tests for the Wave D commands (`depth`, `coverage`,
//! `bedcov`) and the indexed-region modes of `view`.

use std::ffi::OsString;
use std::io::Read;
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use flate2::read::MultiGzDecoder;
use samtools_rs::commands::{bedcov, coverage, depth, index, stats, view};
use samtools_rs::native;
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
fn depth_sam_input_supports_region_restriction() {
    let tmp = tmp_dir("depth-sam-region");
    let sam = tmp.join("in.sam");
    let out = tmp.join("depth.tsv");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "r2\t0\tchr1\t4\t60\t3M\t*\t0\t0\tTGC\t###\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-r",
                "chr1:3-5",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "chr1\t3\t1\nchr1\t4\t2\nchr1\t5\t2\n"
    );
}

#[test]
fn depth_flag_filters_match_default_exclusion_controls() {
    let tmp = tmp_dir("depth-flag-filters");
    let sam = tmp.join("in.sam");
    let default_out = tmp.join("default.tsv");
    let include_dup_out = tmp.join("include-dup.tsv");
    let only_dup_out = tmp.join("only-dup.tsv");
    let require_dup_out = tmp.join("require-dup.tsv");
    let exclude_reverse_out = tmp.join("exclude-reverse.tsv");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:6\n",
            "normal\t0\tchr1\t1\t60\t2M\t*\t0\t0\tAA\t!!\n",
            "dup\t1024\tchr1\t2\t60\t2M\t*\t0\t0\tCC\t!!\n",
            "rev\t16\tchr1\t4\t60\t2M\t*\t0\t0\tGG\t!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &["-o", default_out.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&default_out).unwrap(),
        "chr1\t1\t1\nchr1\t2\t1\nchr1\t4\t1\nchr1\t5\t1\n"
    );

    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-g",
                "DUP",
                "-o",
                include_dup_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&include_dup_out).unwrap(),
        "chr1\t1\t1\nchr1\t2\t2\nchr1\t3\t1\nchr1\t4\t1\nchr1\t5\t1\n"
    );

    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-g",
                "DUP",
                "--incl-flags",
                "DUP",
                "-o",
                only_dup_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&only_dup_out).unwrap(),
        "chr1\t2\t1\nchr1\t3\t1\n"
    );

    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-g",
                "DUP",
                "--require-flags",
                "DUP",
                "-o",
                require_dup_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&require_dup_out).unwrap(),
        "chr1\t2\t1\nchr1\t3\t1\n"
    );

    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-G",
                "REVERSE",
                "-o",
                exclude_reverse_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(exclude_reverse_out).unwrap(),
        "chr1\t1\t1\nchr1\t2\t1\n"
    );
}

#[test]
fn depth_min_read_len_filters_short_alignments() {
    let tmp = tmp_dir("depth-min-read-len");
    let sam = tmp.join("in.sam");
    let out = tmp.join("depth.tsv");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "short\t0\tchr1\t1\t60\t2M\t*\t0\t0\tAA\t!!\n",
            "long\t0\tchr1\t4\t60\t2M1I2M\t*\t0\t0\tCCGTT\t!!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-l",
                "5",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "chr1\t4\t1\nchr1\t5\t1\nchr1\t6\t1\nchr1\t7\t1\n"
    );
}

#[test]
fn depth_multi_input_outputs_one_column_per_input() {
    let p = indexed_bam();
    let tmp = tmp_dir("depth-multi");
    let out = tmp.join("depth.tsv");
    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-r",
                "17:1-10000",
                "-o",
                out.to_str().unwrap(),
                p.to_str().unwrap(),
                p.to_str().unwrap()
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(!text.is_empty());
    for line in text.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 4);
        assert_eq!(fields[0], "17");
        assert_eq!(fields[2], fields[3]);
    }
}

#[test]
fn depth_header_lists_input_columns() {
    let p = indexed_bam();
    let tmp = tmp_dir("depth-header");
    let out = tmp.join("depth.tsv");
    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-H",
                "-r",
                "17:1-10000",
                "-o",
                out.to_str().unwrap(),
                p.to_str().unwrap(),
                p.to_str().unwrap()
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    let mut lines = text.lines();
    assert_eq!(
        lines.next().unwrap(),
        format!("#CHROM\tPOS\t{}\t{}", p.display(), p.display())
    );
    assert!(lines.next().unwrap().starts_with("17\t"));
}

#[test]
fn depth_reads_input_paths_from_file_list() {
    let p = indexed_bam();
    let tmp = tmp_dir("depth-file-list");
    let list = tmp.join("inputs.txt");
    let out = tmp.join("depth.tsv");
    std::fs::write(&list, format!("{}\n{}\n", p.display(), p.display())).unwrap();

    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-H",
                "-f",
                list.to_str().unwrap(),
                "-r",
                "17:1-10000",
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    let mut lines = text.lines();
    assert_eq!(
        lines.next().unwrap(),
        format!("#CHROM\tPOS\t{}\t{}", p.display(), p.display())
    );
    for line in lines {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 4);
        assert_eq!(fields[0], "17");
        assert_eq!(fields[2], fields[3]);
    }
}

#[test]
fn depth_cram_region_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("depth-cram");
    let out = tmp.join("depth.tsv");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "depth",
                "-r",
                "CHROMOSOME_II:2980-2980",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "CHROMOSOME_II\t2980\t1\n"
    );
}

#[test]
fn depth_cram_without_reference_fails_cleanly() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let cram = htslib_fixtures_dir().join("range.cram");

    assert_ne!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "depth",
                "-r",
                "CHROMOSOME_II:2980-2980",
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );
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
fn coverage_multi_input_aggregates_per_reference() {
    let p = indexed_bam();
    let tmp = tmp_dir("coverage-multi");
    let single_out = tmp.join("single.tsv");
    let multi_out = tmp.join("multi.tsv");

    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-r",
                "17:1-10000",
                "-o",
                single_out.to_str().unwrap(),
                p.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-r",
                "17:1-10000",
                "-o",
                multi_out.to_str().unwrap(),
                p.to_str().unwrap(),
                p.to_str().unwrap()
            ]
        ))),
        0
    );

    let single_text = std::fs::read_to_string(single_out).unwrap();
    let multi_text = std::fs::read_to_string(multi_out).unwrap();
    let single_rows: Vec<_> = single_text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect();
    let multi_rows: Vec<_> = multi_text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect();
    assert_eq!(single_rows.len(), 1);
    assert_eq!(multi_rows.len(), 1);

    let single_fields: Vec<_> = single_rows[0].split('\t').collect();
    let multi_fields: Vec<_> = multi_rows[0].split('\t').collect();
    assert_eq!(multi_fields[0], "17");
    assert_eq!(
        multi_fields[3].parse::<u64>().unwrap(),
        single_fields[3].parse::<u64>().unwrap() * 2
    );
    assert!(multi_fields[6].parse::<f64>().unwrap() > single_fields[6].parse::<f64>().unwrap());
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
    let fields: Vec<_> = rows[0].split('\t').collect();
    assert!(fields[7].parse::<f64>().unwrap() > 0.0);
}

#[test]
fn coverage_sam_input_supports_region_restriction() {
    let tmp = tmp_dir("coverage-sam-region");
    let sam = tmp.join("in.sam");
    let out = tmp.join("coverage.tsv");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "r2\t0\tchr1\t3\t30\t4M\t*\t0\t0\tTGCA\tIIII\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-r",
                "chr1:3-5",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    let rows: Vec<_> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(
        rows,
        ["chr1\t3\t5\t2\t3\t100.000000\t1.666667\t24.000000\t45.000000"]
    );
}

#[test]
fn coverage_flag_filters_match_default_exclusion_controls() {
    let tmp = tmp_dir("coverage-flag-filters");
    let sam = tmp.join("in.sam");
    let default_out = tmp.join("default.tsv");
    let include_dup_out = tmp.join("include-dup.tsv");
    let only_dup_out = tmp.join("only-dup.tsv");
    let exclude_reverse_out = tmp.join("exclude-reverse.tsv");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "normal\t0\tchr1\t1\t60\t2M\t*\t0\t0\tAA\tII\n",
            "dup\t1024\tchr1\t2\t60\t2M\t*\t0\t0\tCC\tII\n",
            "rev\t16\tchr1\t4\t60\t2M\t*\t0\t0\tGG\tII\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-H",
                "-o",
                default_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&default_out).unwrap(),
        "chr1\t1\t8\t2\t4\t50.000000\t0.500000\t40.000000\t60.000000\n"
    );

    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-H",
                "--ff",
                "0",
                "-o",
                include_dup_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&include_dup_out).unwrap(),
        "chr1\t1\t8\t3\t5\t62.500000\t0.750000\t40.000000\t60.000000\n"
    );

    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-H",
                "--ff",
                "0",
                "--rf",
                "DUP",
                "-o",
                only_dup_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&only_dup_out).unwrap(),
        "chr1\t1\t8\t1\t2\t25.000000\t0.250000\t40.000000\t60.000000\n"
    );

    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-H",
                "--ff",
                "REVERSE",
                "-o",
                exclude_reverse_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(exclude_reverse_out).unwrap(),
        "chr1\t1\t8\t2\t3\t37.500000\t0.500000\t40.000000\t60.000000\n"
    );
}

#[test]
fn coverage_reads_input_paths_from_bam_list() {
    let tmp = tmp_dir("coverage-bam-list");
    let sam = tmp.join("in.sam");
    let list = tmp.join("inputs.txt");
    let out = tmp.join("coverage.tsv");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "normal\t0\tchr1\t1\t60\t2M\t*\t0\t0\tAA\tII\n",
            "rev\t16\tchr1\t4\t60\t2M\t*\t0\t0\tGG\tII\n",
        ),
    )
    .unwrap();
    std::fs::write(&list, format!("{}\n{}\n", sam.display(), sam.display())).unwrap();

    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-H",
                "-b",
                list.to_str().unwrap(),
                "-o",
                out.to_str().unwrap()
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "chr1\t1\t8\t4\t4\t50.000000\t1.000000\t40.000000\t60.000000\n"
    );
}

#[test]
fn coverage_min_read_len_filters_short_alignments() {
    let tmp = tmp_dir("coverage-min-read-len");
    let sam = tmp.join("in.sam");
    let out = tmp.join("coverage.tsv");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "short\t0\tchr1\t1\t60\t2M\t*\t0\t0\tAA\tII\n",
            "long\t0\tchr1\t4\t60\t2M1I2M\t*\t0\t0\tCCGTT\tIIIII\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-H",
                "-l",
                "5",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "chr1\t1\t8\t1\t4\t50.000000\t0.500000\t40.000000\t60.000000\n"
    );
}

#[test]
fn coverage_max_depth_caps_reported_depth_metrics() {
    let tmp = tmp_dir("coverage-max-depth");
    let sam = tmp.join("in.sam");
    let default_out = tmp.join("default.tsv");
    let capped_out = tmp.join("capped.tsv");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:4\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII\n",
            "r2\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\tIIII\n",
            "r3\t0\tchr1\t1\t60\t4M\t*\t0\t0\tGGGG\tIIII\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-H",
                "-o",
                default_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-H",
                "-d",
                "1",
                "-o",
                capped_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(default_out).unwrap(),
        "chr1\t1\t4\t3\t4\t100.000000\t3.000000\t40.000000\t60.000000\n"
    );
    assert_eq!(
        std::fs::read_to_string(capped_out).unwrap(),
        "chr1\t1\t4\t3\t4\t100.000000\t1.000000\t40.000000\t60.000000\n"
    );
}

#[test]
fn coverage_min_depth_and_base_quality_filter_covbases() {
    let p = indexed_bam();
    let tmp = tmp_dir("coverage-filters");
    let default_out = tmp.join("coverage-default.tsv");
    let filtered_out = tmp.join("coverage-filtered.tsv");

    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-r",
                "17:1-10000",
                "-o",
                default_out.to_str().unwrap(),
                p.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-r",
                "17:1-10000",
                "--min-depth",
                "2",
                "-Q",
                "30",
                "-o",
                filtered_out.to_str().unwrap(),
                p.to_str().unwrap()
            ]
        ))),
        0
    );

    let default_row = std::fs::read_to_string(default_out)
        .unwrap()
        .lines()
        .find(|line| !line.starts_with('#'))
        .unwrap()
        .to_owned();
    let filtered_row = std::fs::read_to_string(filtered_out)
        .unwrap()
        .lines()
        .find(|line| !line.starts_with('#'))
        .unwrap()
        .to_owned();
    let default_fields: Vec<_> = default_row.split('\t').collect();
    let filtered_fields: Vec<_> = filtered_row.split('\t').collect();
    let default_covbases = default_fields[4].parse::<usize>().unwrap();
    let filtered_covbases = filtered_fields[4].parse::<usize>().unwrap();

    assert!(filtered_covbases < default_covbases);
    assert!(filtered_fields[7].parse::<f64>().unwrap() >= 30.0);
}

#[test]
fn coverage_cram_region_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("coverage-cram");
    let out = tmp.join("coverage.tsv");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "coverage",
                "-r",
                "CHROMOSOME_II:2980-2980",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(
        text.contains(
            "CHROMOSOME_II\t2980\t2980\t1\t1\t100.000000\t1.000000\t35.000000\t60.000000\n"
        )
    );
}

#[test]
fn coverage_cram_without_reference_fails_cleanly() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let cram = htslib_fixtures_dir().join("range.cram");

    assert_ne!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "coverage",
                "-r",
                "CHROMOSOME_II:2980-2980",
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );
}

#[test]
fn stats_cram_uses_top_level_reference() {
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
                "stats",
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );
}

#[test]
fn stats_cram_without_reference_fails_cleanly() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let cram = htslib_fixtures_dir().join("range.cram");

    assert_ne!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &["stats", cram.to_str().unwrap()],
        ))),
        0
    );
}

#[test]
fn stats_cram_region_and_target_file_restrict_summary_counts() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");
    let tmp = tmp_dir("stats-cram-region");
    let targets = tmp.join("targets.txt");
    let full_out = tmp.join("full.stats");
    let region_out = tmp.join("region.stats");
    let target_out = tmp.join("target.stats");

    std::fs::write(&targets, "CHROMOSOME_II 2980 2980\n").unwrap();

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "stats",
                "-o",
                full_out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "stats",
                "-o",
                region_out.to_str().unwrap(),
                cram.to_str().unwrap(),
                "CHROMOSOME_II:2980-2980",
            ],
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "stats",
                "-t",
                targets.to_str().unwrap(),
                "-o",
                target_out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );

    let full_total = stats_sn_value(&std::fs::read_to_string(full_out).unwrap(), "sequences");
    let region_total = stats_sn_value(&std::fs::read_to_string(region_out).unwrap(), "sequences");
    let target_total = stats_sn_value(&std::fs::read_to_string(target_out).unwrap(), "sequences");

    assert!(region_total > 0);
    assert!(region_total < full_total);
    assert_eq!(target_total, region_total);
}

#[test]
fn stats_bam_positional_region_restricts_summary_counts() {
    let bam = indexed_bam();
    let tmp = tmp_dir("stats-region");
    let full_out = tmp.join("full.stats");
    let region_out = tmp.join("region.stats");

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &["-o", full_out.to_str().unwrap(), bam.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-o",
                region_out.to_str().unwrap(),
                bam.to_str().unwrap(),
                "17:2726-2730"
            ]
        ))),
        0
    );

    let full_total = stats_sn_value(&std::fs::read_to_string(full_out).unwrap(), "sequences");
    let region_total = stats_sn_value(&std::fs::read_to_string(region_out).unwrap(), "sequences");

    assert!(region_total > 0);
    assert!(region_total < full_total);
}

#[test]
fn stats_bam_target_file_restricts_summary_counts() {
    let bam = indexed_bam();
    let tmp = tmp_dir("stats-targets");
    let targets = tmp.join("targets.txt");
    let full_out = tmp.join("full.stats");
    let target_out = tmp.join("target.stats");

    std::fs::write(&targets, "# comments are ignored\n17 2726 2730\n").unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &["-o", full_out.to_str().unwrap(), bam.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-t",
                targets.to_str().unwrap(),
                "-o",
                target_out.to_str().unwrap(),
                bam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let full_total = stats_sn_value(&std::fs::read_to_string(full_out).unwrap(), "sequences");
    let target_total = stats_sn_value(&std::fs::read_to_string(target_out).unwrap(), "sequences");

    assert!(target_total > 0);
    assert!(target_total < full_total);
}

#[test]
fn stats_remove_dups_filters_duplicate_marked_reads() {
    let tmp = tmp_dir("stats-remove-dups");
    let sam = tmp.join("dups.sam");
    let all_out = tmp.join("all.stats");
    let dedup_out = tmp.join("dedup.stats");
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!
r2\t1024\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &["-o", all_out.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-d",
                "-o",
                dedup_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let all_text = std::fs::read_to_string(all_out).unwrap();
    let dedup_text = std::fs::read_to_string(dedup_out).unwrap();

    assert_eq!(stats_sn_value(&all_text, "raw total sequences"), 2);
    assert_eq!(stats_sn_value(&all_text, "filtered sequences"), 0);
    assert_eq!(stats_sn_value(&all_text, "sequences"), 2);
    assert_eq!(stats_sn_value(&all_text, "reads duplicated"), 1);
    assert_eq!(stats_sn_value(&dedup_text, "raw total sequences"), 2);
    assert_eq!(stats_sn_value(&dedup_text, "filtered sequences"), 1);
    assert_eq!(stats_sn_value(&dedup_text, "sequences"), 1);
    assert_eq!(stats_sn_value(&dedup_text, "reads duplicated"), 0);
}

#[test]
fn stats_target_regions_deduplicate_overlaps_for_sam_and_bam() {
    let fixtures = fixtures_dir().join("stat");
    let sam = fixtures.join("11_target.sam");
    let bam = fixtures.join("11_target.bam");
    let targets = fixtures.join("11.stats.targets");
    let tmp = tmp_dir("stats-target-dedup");
    let sam_out = tmp.join("sam.stats");
    let bam_out = tmp.join("bam.stats");

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-t",
                targets.to_str().unwrap(),
                "-o",
                sam_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-o",
                bam_out.to_str().unwrap(),
                bam.to_str().unwrap(),
                "ref1:10-24",
                "ref1:30-46",
                "ref1:39-56",
            ]
        ))),
        0
    );

    let sam_total = stats_sn_value(&std::fs::read_to_string(sam_out).unwrap(), "sequences");
    let bam_total = stats_sn_value(&std::fs::read_to_string(bam_out).unwrap(), "sequences");

    assert_eq!(sam_total, 26);
    assert_eq!(bam_total, 26);
}

fn stats_sn_value(text: &str, key: &str) -> u64 {
    text.lines()
        .find_map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            (fields.len() >= 3 && fields[0] == "SN" && fields[1] == format!("{key}:"))
                .then(|| fields[2].parse().unwrap())
        })
        .unwrap()
}

fn stats_sn_text<'a>(text: &'a str, key: &str) -> &'a str {
    text.lines()
        .find_map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            (fields.len() >= 3 && fields[0] == "SN" && fields[1] == format!("{key}:"))
                .then(|| fields[2])
        })
        .unwrap_or("")
}

fn quality_hist_value(text: &str, label: &str, cycle: usize, quality: usize) -> u64 {
    text.lines()
        .find_map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            (fields.len() >= quality + 3
                && fields[0] == label
                && fields[1].parse::<usize>().unwrap() == cycle)
                .then(|| fields[quality + 2].parse().unwrap())
        })
        .unwrap()
}

fn gc_hist_value(text: &str, label: &str, percent: &str) -> u64 {
    text.lines()
        .find_map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            (fields.len() >= 3 && fields[0] == label && fields[1] == percent)
                .then(|| fields[2].parse().unwrap())
        })
        .unwrap_or(0)
}

#[test]
fn stats_emits_insert_size_and_supplementary_sn_lines() {
    let tmp = tmp_dir("stats-insert-size");
    let sam = tmp.join("paired.sam");
    let out = tmp.join("paired.stats");
    // Three records:
    //  - r1 first/forward, r1 mate/reverse: classic FR pair with TLEN 100
    //  - supp: a supplementary alignment of r1 that must NOT contribute to
    //          the IS bins or to raw totals
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
r1\t99\tchr1\t1\t60\t10M\t=\t91\t100\tACGTACGTAC\t!!!!!!!!!!
r1\t147\tchr1\t91\t60\t10M\t=\t1\t-100\tACGTACGTAC\t!!!!!!!!!!
r1\t2147\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t!!!!!!!!!!
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &["-o", out.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert_eq!(stats_sn_value(&text, "raw total sequences"), 2);
    assert_eq!(stats_sn_value(&text, "supplementary alignments"), 1);
    assert_eq!(stats_sn_value(&text, "inward oriented pairs"), 1);
    assert_eq!(stats_sn_value(&text, "outward oriented pairs"), 0);
    assert_eq!(stats_sn_value(&text, "pairs with other orientation"), 0);
    assert_eq!(stats_sn_text(&text, "insert size average"), "100.0");
    assert_eq!(
        stats_sn_text(&text, "insert size standard deviation"),
        "0.0"
    );
    assert_eq!(
        stats_sn_text(&text, "percentage of properly paired reads (%)"),
        "100.0"
    );
}

#[test]
fn stats_classifies_outward_and_other_orientation() {
    let tmp = tmp_dir("stats-orientation");
    let sam = tmp.join("oriented.sam");
    let out = tmp.join("oriented.stats");
    // Two pairs on the same chromosome, both mapped:
    //   pair A: outward — read1 is reverse (5'-most reverse), read2 forward
    //   pair B: same-direction FF — both forward → "other"
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
a\t83\tchr1\t1\t60\t10M\t=\t91\t100\tACGTACGTAC\t!!!!!!!!!!
a\t163\tchr1\t91\t60\t10M\t=\t1\t-100\tACGTACGTAC\t!!!!!!!!!!
b\t65\tchr1\t1\t60\t10M\t=\t91\t100\tACGTACGTAC\t!!!!!!!!!!
b\t129\tchr1\t91\t60\t10M\t=\t1\t-100\tACGTACGTAC\t!!!!!!!!!!
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &["-o", out.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert_eq!(stats_sn_value(&text, "inward oriented pairs"), 0);
    assert_eq!(stats_sn_value(&text, "outward oriented pairs"), 1);
    assert_eq!(stats_sn_value(&text, "pairs with other orientation"), 1);
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
fn bedcov_sam_input_supports_depth_and_count_columns() {
    let tmp = tmp_dir("bedcov-sam");
    let bed = tmp.join("r.bed");
    let sam = tmp.join("in.sam");
    std::fs::write(&bed, "chr1\t2\t5\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "r2\t0\tchr1\t3\t30\t4M\t*\t0\t0\tTGCA\tIIII\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(bedcov::main(&argv(
            "bedcov",
            &[
                "-d",
                "2",
                "-c",
                bed.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
}

#[test]
fn bedcov_cram_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bedcov-cram");
    let bed = tmp.join("r.bed");
    let fixtures = htslib_fixtures_dir();
    let reference = fixtures.join("ce.fa");
    let cram = fixtures.join("range.cram");
    std::fs::write(&bed, "CHROMOSOME_II\t2979\t2980\n").unwrap();

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "bedcov",
                bed.to_str().unwrap(),
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );
}

#[test]
fn bedcov_cram_without_reference_fails_cleanly() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("bedcov-cram-no-ref");
    let bed = tmp.join("r.bed");
    let cram = htslib_fixtures_dir().join("range.cram");
    std::fs::write(&bed, "CHROMOSOME_II\t2979\t2980\n").unwrap();

    assert_ne!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &["bedcov", bed.to_str().unwrap(), cram.to_str().unwrap()],
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

#[test]
fn stats_is_sorted_reflects_header_sort_order() {
    let tmp = tmp_dir("stats-is-sorted");
    let sam_sorted = tmp.join("sorted.sam");
    let sam_unsorted = tmp.join("unsorted.sam");
    let out_sorted = tmp.join("sorted.stats");
    let out_unsorted = tmp.join("unsorted.stats");
    std::fs::write(
        &sam_sorted,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!
",
    )
    .unwrap();
    std::fs::write(
        &sam_unsorted,
        "\
@HD\tVN:1.6\tSO:queryname
@SQ\tSN:chr1\tLN:100
r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-o",
                out_sorted.to_str().unwrap(),
                sam_sorted.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-o",
                out_unsorted.to_str().unwrap(),
                sam_unsorted.to_str().unwrap()
            ]
        ))),
        0
    );

    let sorted_text = std::fs::read_to_string(&out_sorted).unwrap();
    let unsorted_text = std::fs::read_to_string(&out_unsorted).unwrap();
    assert_eq!(stats_sn_value(&sorted_text, "is sorted"), 1);
    assert_eq!(stats_sn_value(&unsorted_text, "is sorted"), 0);
}

#[test]
fn coverage_histogram_emits_ascii_plot() {
    let tmp = tmp_dir("coverage-hist");
    let sam = tmp.join("in.sam");
    let out = tmp.join("hist.txt");
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:80
r1\t0\tchr1\t1\t60\t40M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
r2\t0\tchr1\t1\t60\t40M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(coverage::main(&argv(
            "coverage",
            &[
                "-m",
                "-w",
                "20",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("chr1"));
    // Histogram rows begin with `>` and a percentage. The bottom row is the
    // 0.00% threshold; ensure it's present and contains at least one filled
    // glyph (`:` or `.`).
    let bottom = text
        .lines()
        .find(|l| l.contains("0.00%"))
        .expect("histogram has a 0% row");
    assert!(bottom.contains(':') || bottom.contains('.'));
}

#[test]
fn stats_emits_sequence_length_sn_lines() {
    let tmp = tmp_dir("stats-seq-len");
    let sam = tmp.join("len.sam");
    let out = tmp.join("len.stats");
    // Two paired reads with sequence length 10; average quality is the
    // mean of `!` (ASCII 33 → Phred 0) values.
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
r1\t99\tchr1\t1\t60\t10M\t=\t91\t100\tACGTACGTAC\tIIIIIIIIII
r1\t147\tchr1\t91\t60\t10M\t=\t1\t-100\tACGTACGTAC\tIIIIIIIIII
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &["-o", out.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert_eq!(stats_sn_value(&text, "total length"), 20);
    assert_eq!(stats_sn_value(&text, "total first fragment length"), 10);
    assert_eq!(stats_sn_value(&text, "total last fragment length"), 10);
    assert_eq!(stats_sn_value(&text, "maximum length"), 10);
    assert_eq!(stats_sn_text(&text, "average length"), "10");
    // Quality 'I' is ASCII 73 → Phred 73-33=40.
    assert_eq!(stats_sn_text(&text, "average quality"), "40.0");
}

#[test]
fn stats_emits_first_and_last_fragment_quality_histograms() {
    let tmp = tmp_dir("stats-quality-hist");
    let sam = tmp.join("qual.sam");
    let out = tmp.join("qual.stats");
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
r1\t99\tchr1\t1\t60\t10M\t=\t91\t100\tACGTACGTAC\tIIIIIIIIII
r1\t147\tchr1\t91\t60\t10M\t=\t1\t-100\tTGCATGCATG\t!!!!!!!!!!
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &["-o", out.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert_eq!(quality_hist_value(&text, "FFQ", 1, 40), 1);
    assert_eq!(quality_hist_value(&text, "FFQ", 10, 40), 1);
    assert_eq!(quality_hist_value(&text, "LFQ", 1, 0), 1);
    assert_eq!(quality_hist_value(&text, "LFQ", 10, 0), 1);
}

#[test]
fn stats_emits_first_and_last_fragment_gc_histograms() {
    let tmp = tmp_dir("stats-gc-hist");
    let sam = tmp.join("gc.sam");
    let all_out = tmp.join("all.stats");
    let dedup_out = tmp.join("dedup.stats");
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
r1\t99\tchr1\t1\t60\t10M\t=\t91\t100\tGGGGGAAAAA\tIIIIIIIIII
r1\t147\tchr1\t91\t60\t10M\t=\t1\t-100\tCCCCAAAAAA\tIIIIIIIIII
dup\t1123\tchr1\t201\t60\t10M\t=\t291\t100\tGGGGGGGGGG\tIIIIIIIIII
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &["-o", all_out.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-d",
                "-o",
                dedup_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    let all_text = std::fs::read_to_string(all_out).unwrap();
    assert_eq!(gc_hist_value(&all_text, "GCF", "50.00"), 1);
    assert_eq!(gc_hist_value(&all_text, "GCF", "100.00"), 1);
    assert_eq!(gc_hist_value(&all_text, "GCL", "40.00"), 1);

    let dedup_text = std::fs::read_to_string(dedup_out).unwrap();
    assert_eq!(gc_hist_value(&dedup_text, "GCF", "50.00"), 1);
    assert_eq!(gc_hist_value(&dedup_text, "GCF", "100.00"), 0);
    assert_eq!(gc_hist_value(&dedup_text, "GCL", "40.00"), 1);
}

#[test]
fn stats_emits_cigar_walk_coverage_histogram() {
    let tmp = tmp_dir("stats-cov");
    let sam = tmp.join("cov.sam");
    let out = tmp.join("cov.stats");
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:10
r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII
r2\t0\tchr1\t3\t60\t4M\t*\t0\t0\tCCCC\tIIII
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &["-o", out.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.contains("# Coverage distribution."));
    assert!(text.contains("COV\t[1-1]\t1\t4\n"));
    assert!(text.contains("COV\t[2-2]\t2\t2\n"));
}

#[test]
fn stats_coverage_option_groups_cov_bins() {
    let tmp = tmp_dir("stats-cov-bins");
    let sam = tmp.join("cov.sam");
    let out = tmp.join("cov.stats");
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:10
r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII
r2\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\tIIII
r3\t0\tchr1\t3\t60\t4M\t*\t0\t0\tGGGG\tIIII
r4\t0\tchr1\t5\t60\t2M\t*\t0\t0\tTT\tII
r5\t0\tchr1\t5\t60\t2M\t*\t0\t0\tAA\tII
r6\t0\tchr1\t5\t60\t2M\t*\t0\t0\tCC\tII
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-c",
                "2,4,2",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(!text.contains("COV\t[1-1]"));
    assert!(text.contains("COV\t[2-3]\t2\t4\n"));
    assert!(text.contains("COV\t[4-4]\t4\t2\n"));
}

#[test]
fn stats_cov_threshold_reports_target_percentage() {
    let tmp = tmp_dir("stats-cov-threshold");
    let sam = tmp.join("cov.sam");
    let targets = tmp.join("targets.bed");
    let out = tmp.join("cov.stats");
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:10
r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII
r2\t0\tchr1\t3\t60\t4M\t*\t0\t0\tCCCC\tIIII
",
    )
    .unwrap();
    std::fs::write(&targets, "chr1\t1\t6\n").unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-t",
                targets.to_str().unwrap(),
                "-g",
                "1",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.contains("SN\tbases inside the target:\t6\n"));
    assert!(text.contains("SN\tpercentage of target genome with coverage > 1 (%):\t33.33\n"));
}

#[test]
fn stats_filters_required_and_filtering_flags() {
    let tmp = tmp_dir("stats-flag-filters");
    let sam = tmp.join("flags.sam");
    let required_out = tmp.join("required.stats");
    let filtered_out = tmp.join("filtered.stats");
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
read1\t65\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII
read2\t0\tchr1\t5\t60\t4M\t*\t0\t0\tCCCC\tIIII
read3\t4\t*\t0\t0\t*\t*\t0\t0\tGGGG\tIIII
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-f",
                "READ1",
                "-o",
                required_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
    let required_text = std::fs::read_to_string(required_out).unwrap();
    assert!(required_text.contains("SN\traw total sequences:\t3\n"));
    assert!(required_text.contains("SN\tfiltered sequences:\t2\n"));
    assert!(required_text.contains("SN\tsequences:\t1\n"));
    assert!(required_text.contains("SN\t1st fragments:\t1\n"));

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-F",
                "UNMAP",
                "-o",
                filtered_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
    let filtered_text = std::fs::read_to_string(filtered_out).unwrap();
    assert!(filtered_text.contains("SN\traw total sequences:\t3\n"));
    assert!(filtered_text.contains("SN\tfiltered sequences:\t1\n"));
    assert!(filtered_text.contains("SN\tsequences:\t2\n"));
    assert!(filtered_text.contains("SN\treads mapped:\t2\n"));
}

#[test]
fn stats_filters_by_exact_read_length() {
    let tmp = tmp_dir("stats-read-length-filter");
    let sam = tmp.join("lengths.sam");
    let out = tmp.join("lengths.stats");
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
short\t0\tchr1\t1\t60\t3M\t*\t0\t0\tAAA\tIII
long\t0\tchr1\t5\t60\t5M\t*\t0\t0\tCCCCC\tIIIII
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-l",
                "5",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.contains("SN\traw total sequences:\t2\n"));
    assert!(text.contains("SN\tfiltered sequences:\t1\n"));
    assert!(text.contains("SN\tsequences:\t1\n"));
    assert!(text.contains("SN\ttotal length:\t5\n"));
    assert!(text.contains("SN\tmaximum length:\t5\n"));
}

#[test]
fn stats_emits_bases_mapped_and_error_rate_sn_lines() {
    let tmp = tmp_dir("stats-error-rate");
    let sam = tmp.join("nm.sam");
    let out = tmp.join("nm.stats");
    // Two mapped reads with NM:i:1 and NM:i:0 → 1 mismatch over 20
    // CIGAR-derived mapped bases gives error rate 0.05.
    std::fs::write(
        &sam,
        "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:1
r2\t0\tchr1\t11\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0
",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &["-o", out.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert_eq!(stats_sn_value(&text, "bases mapped"), 20);
    assert_eq!(stats_sn_value(&text, "bases mapped (cigar)"), 20);
    assert_eq!(stats_sn_value(&text, "mismatches"), 1);
    let error_rate = stats_sn_text(&text, "error rate");
    // 1/20 = 0.05 → upstream prints `5.000000e-02`.
    assert!(error_rate.starts_with("5.0"));
    assert!(error_rate.contains('e'));
}
