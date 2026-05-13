//! Smoke-tests for `cat`, `reheader`, `fastq`, `samples`, `idxstats`,
//! `flagstat`, `index`, `faidx`, `import`, `bedcov`, `rmdup`, `split`.

use std::ffi::OsString;
use std::path::PathBuf;
use std::process::ExitCode;

use samtools_rs::commands::{
    bedcov, cat, faidx, fastq, flagstat, fqidx, idxstats, import, index, reheader, rmdup, samples,
    split,
};

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
    let p = std::env::temp_dir().join(format!("samtools-rs-misc-{}-{}", name, std::process::id()));
    std::fs::create_dir_all(&p).unwrap();
    p
}

fn argv(name: &str, rest: &[&str]) -> Vec<OsString> {
    std::iter::once(OsString::from(name))
        .chain(rest.iter().map(OsString::from))
        .collect()
}

fn sample_bam() -> PathBuf {
    fixtures_dir().join("checksum").join("chk1.bam")
}

#[test]
fn flagstat_succeeds() {
    let p = sample_bam();
    assert_eq!(
        exit_to_u8(flagstat::main(&argv("flagstat", &[p.to_str().unwrap()]))),
        0
    );
}

#[test]
fn idxstats_succeeds() {
    let p = sample_bam();
    // Build index next to a copy so the fixture stays clean.
    let tmp = tmp_dir("idx");
    let bam = tmp.join("in.bam");
    std::fs::copy(&p, &bam).unwrap();
    assert_eq!(
        exit_to_u8(index::main(&argv("index", &[bam.to_str().unwrap()]))),
        0
    );
    assert_eq!(
        exit_to_u8(idxstats::main(&argv("idxstats", &[bam.to_str().unwrap()]))),
        0
    );
}

#[test]
fn idxstats_bam_without_index_uses_slow_path() {
    let tmp = tmp_dir("idx-slow-bam");
    let bam = tmp.join("in.bam");
    std::fs::copy(sample_bam(), &bam).unwrap();
    assert_eq!(
        exit_to_u8(idxstats::main(&argv("idxstats", &[bam.to_str().unwrap()]))),
        0
    );
}

#[test]
fn idxstats_sam_uses_slow_path() {
    let tmp = tmp_dir("idx-slow-sam");
    let sam = tmp.join("in.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@SQ\tSN:chr2\tLN:4\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "r2\t4\tchr1\t1\t0\t*\t*\t0\t0\t*\t*\n",
            "r3\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n",
        ),
    )
    .unwrap();
    assert_eq!(
        exit_to_u8(idxstats::main(&argv("idxstats", &[sam.to_str().unwrap()]))),
        0
    );
}

#[test]
fn samples_succeeds() {
    let p = sample_bam();
    assert_eq!(
        exit_to_u8(samples::main(&argv("samples", &[p.to_str().unwrap()]))),
        0
    );
}

#[test]
fn samples_reports_index_presence() {
    let tmp = tmp_dir("samples-index");
    let bam = tmp.join("in.bam");
    let out = tmp.join("samples.txt");
    std::fs::copy(sample_bam(), &bam).unwrap();

    assert_eq!(
        exit_to_u8(samples::main(&argv(
            "samples",
            &["-i", "-o", out.to_str().unwrap(), bam.to_str().unwrap()]
        ))),
        0
    );
    assert!(std::fs::read_to_string(&out).unwrap().ends_with("\tN\n"));

    assert_eq!(
        exit_to_u8(index::main(&argv("index", &[bam.to_str().unwrap()]))),
        0
    );
    assert_eq!(
        exit_to_u8(samples::main(&argv(
            "samples",
            &["-i", "-o", out.to_str().unwrap(), bam.to_str().unwrap()]
        ))),
        0
    );
    assert!(std::fs::read_to_string(out).unwrap().ends_with("\tY\n"));
}

#[test]
fn samples_custom_index_pair_reports_index_presence() {
    let tmp = tmp_dir("samples-custom-index");
    let bam = tmp.join("in.bam");
    let out = tmp.join("samples.txt");
    let custom_index = tmp.join("in.custom.bai");
    std::fs::copy(sample_bam(), &bam).unwrap();

    assert_eq!(
        exit_to_u8(index::main(&argv("index", &[bam.to_str().unwrap()]))),
        0
    );
    std::fs::rename(tmp.join("in.bam.bai"), &custom_index).unwrap();

    assert_eq!(
        exit_to_u8(samples::main(&argv(
            "samples",
            &[
                "-X",
                "-i",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
                custom_index.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert!(std::fs::read_to_string(out).unwrap().ends_with("\tY\n"));
}

#[test]
fn samples_cram_header_succeeds() {
    let tmp = tmp_dir("samples-cram");
    let out = tmp.join("samples.txt");
    let cram = fixtures_dir().join("dat/test_input_1_a.cram");

    assert_eq!(
        exit_to_u8(samples::main(&argv(
            "samples",
            &["-o", out.to_str().unwrap(), cram.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        format!(".\t{}\n", cram.display())
    );
}

#[test]
fn samples_matches_reference_from_fasta_and_list() {
    let tmp = tmp_dir("samples-ref");
    let sam = tmp.join("in.sam");
    let fa = tmp.join("ref.fa");
    let fa_list = tmp.join("refs.txt");
    let out = tmp.join("samples.txt");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:4\n",
            "@RG\tID:g1\tSM:s1\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();
    std::fs::write(&fa, ">chr1\nACGT\n").unwrap();
    std::fs::write(&fa_list, format!("{}\n", fa.display())).unwrap();

    assert_eq!(
        exit_to_u8(samples::main(&argv(
            "samples",
            &[
                "-f",
                fa.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        format!("s1\t{}\t{}\n", sam.display(), fa.display())
    );

    assert_eq!(
        exit_to_u8(samples::main(&argv(
            "samples",
            &[
                "-F",
                fa_list.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        format!("s1\t{}\t{}\n", sam.display(), fa.display())
    );
}

#[test]
fn cat_two_succeeds() {
    let tmp = tmp_dir("cat");
    let out = tmp.join("cat.bam");
    let p = sample_bam();
    assert_eq!(
        exit_to_u8(cat::main(&argv(
            "cat",
            &[
                p.to_str().unwrap(),
                p.to_str().unwrap(),
                "-o",
                out.to_str().unwrap()
            ]
        ))),
        0
    );
    assert!(out.exists());
}

#[test]
fn reheader_succeeds() {
    let tmp = tmp_dir("reh");
    let hdr = fixtures_dir().join("reheader").join("hdr.sam");
    let p = sample_bam();
    // reheader writes to stdout — exercise the path without capturing.
    assert_eq!(
        exit_to_u8(reheader::main(&argv(
            "reheader",
            &[hdr.to_str().unwrap(), p.to_str().unwrap()]
        ))),
        0
    );
    let _ = tmp; // keep around
}

#[test]
fn fastq_from_sam_succeeds() {
    let p = fixtures_dir().join("dat").join("view.001.sam");
    assert_eq!(
        exit_to_u8(fastq::main(&argv("fastq", &[p.to_str().unwrap()]))),
        0
    );
}

#[test]
fn fastq_filters_by_required_flags() {
    let tmp = tmp_dir("fastq-filter");
    let sam = tmp.join("in.sam");
    let out = tmp.join("reads.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "mapped\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tTGCA\t####\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-f",
                "4",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@unmapped\nTGCA\n+\n####\n"
    );
}

#[test]
fn fastq_zero_output_writes_single_stream_file() {
    let tmp = tmp_dir("fastq-zero");
    let sam = tmp.join("in.sam");
    let out = tmp.join("reads.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r0\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &["-0", out.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@r0\nACGT\n+\n!!!!\n"
    );
}

#[test]
fn fastq_appends_pair_suffixes_unless_suppressed() {
    let tmp = tmp_dir("fastq-suffix");
    let sam = tmp.join("in.sam");
    let suffixed = tmp.join("suffixed.fq");
    let unsuffixed = tmp.join("unsuffixed.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "pair\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\n",
            "pair\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &["-o", suffixed.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-n",
                "-o",
                unsuffixed.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(suffixed).unwrap(),
        "@pair/1\nACGT\n+\n!!!!\n@pair/2\nTGCA\n+\n####\n"
    );
    assert_eq!(
        std::fs::read_to_string(unsuffixed).unwrap(),
        "@pair\nACGT\n+\n!!!!\n@pair\nTGCA\n+\n####\n"
    );
}

#[test]
fn fastq_splits_read1_read2_and_singleton_outputs() {
    let tmp = tmp_dir("fastq-split");
    let sam = tmp.join("in.sam");
    let r1 = tmp.join("r1.fq");
    let r2 = tmp.join("r2.fq");
    let singleton = tmp.join("single.fq");
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

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-1",
                r1.to_str().unwrap(),
                "-2",
                r2.to_str().unwrap(),
                "-s",
                singleton.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(r1).unwrap(),
        "@pair\nACGT\n+\n!!!!\n"
    );
    assert_eq!(
        std::fs::read_to_string(r2).unwrap(),
        "@pair\nTGCA\n+\n####\n"
    );
    assert_eq!(
        std::fs::read_to_string(singleton).unwrap(),
        "@solo\nNNNN\n+\n$$$$\n"
    );
}

#[test]
fn fastq_zero_routes_unpaired_reads_in_split_mode() {
    let tmp = tmp_dir("fastq-split-zero");
    let sam = tmp.join("in.sam");
    let r1 = tmp.join("r1.fq");
    let r2 = tmp.join("r2.fq");
    let other = tmp.join("other.fq");
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

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-1",
                r1.to_str().unwrap(),
                "-2",
                r2.to_str().unwrap(),
                "-0",
                other.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(r1).unwrap(),
        "@pair\nACGT\n+\n!!!!\n"
    );
    assert_eq!(
        std::fs::read_to_string(r2).unwrap(),
        "@pair\nTGCA\n+\n####\n"
    );
    assert_eq!(
        std::fs::read_to_string(other).unwrap(),
        "@solo\nNNNN\n+\n$$$$\n"
    );
}

#[test]
fn faidx_builds_index() {
    let tmp = tmp_dir("fai");
    let src = fixtures_dir().join("dat").join("dict.fa");
    let copy = tmp.join("ref.fa");
    std::fs::copy(&src, &copy).unwrap();
    assert_eq!(
        exit_to_u8(faidx::main(&argv("faidx", &[copy.to_str().unwrap()]))),
        0
    );
    assert!(tmp.join("ref.fa.fai").exists());
}

#[test]
fn faidx_extracts_region_to_file() {
    let tmp = tmp_dir("fai-region");
    let fa = tmp.join("ref.fa");
    let out = tmp.join("out.fa");
    std::fs::write(&fa, ">chr1\nACGTACGTACGT\n>chr2\nTTTTCCCC\n").unwrap();
    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "--length",
                "4",
                "-o",
                out.to_str().unwrap(),
                fa.to_str().unwrap(),
                "chr1:3-10"
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        ">chr1:3-10 length: 8\nGTAC\nGTAC\n"
    );
}

#[test]
fn faidx_extracts_regions_from_file() {
    let tmp = tmp_dir("fai-region-file");
    let fa = tmp.join("ref.fa");
    let regions = tmp.join("regions.txt");
    let out = tmp.join("out.fa");
    std::fs::write(&fa, ">chr1\nACGTACGTACGT\n>chr2\nTTTTCCCC\n").unwrap();
    std::fs::write(&regions, "chr1:1-4\nchr2:5-8\n").unwrap();
    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "-r",
                regions.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                fa.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        ">chr1:1-4 length: 4\nACGT\n>chr2:5-8 length: 4\nCCCC\n"
    );
}

#[test]
fn fqidx_extracts_region_to_file() {
    let tmp = tmp_dir("fqi-region");
    let fq = tmp.join("reads.fq");
    let out = tmp.join("out.fq");
    std::fs::write(&fq, "@r1\nACGTACGT\n+\nabcdefgh\n").unwrap();
    assert_eq!(
        exit_to_u8(fqidx::main(&argv(
            "fqidx",
            &[
                "--length",
                "4",
                "-o",
                out.to_str().unwrap(),
                fq.to_str().unwrap(),
                "r1:2-7"
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@r1:2-7 length: 6\nCGTA\nCG\n+\nbcde\nfg\n"
    );
}

#[test]
fn faidx_fastq_mode_extracts_fastq_region() {
    let tmp = tmp_dir("fai-fastq-region");
    let fq = tmp.join("reads.fq");
    let regions = tmp.join("regions.txt");
    let out = tmp.join("out.fq");
    std::fs::write(&fq, "@r1\nACGTACGT\n+\nabcdefgh\n").unwrap();
    std::fs::write(&regions, "r1:1-4\n").unwrap();
    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "-f",
                fq.to_str().unwrap(),
                "-r",
                regions.to_str().unwrap(),
                "-o",
                out.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@r1:1-4 length: 4\nACGT\n+\nabcd\n"
    );
}

#[test]
fn faidx_reverse_complement_marks_default_rc() {
    let tmp = tmp_dir("fai-rc-default");
    let fa = tmp.join("ref.fa");
    let out = tmp.join("out.fa");
    std::fs::write(&fa, ">rc\nACGTMRWSYKVHDBN\n").unwrap();
    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "-i",
                "-o",
                out.to_str().unwrap(),
                fa.to_str().unwrap(),
                "rc"
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        ">rc/rc\nNVHDBMRSWYKACGT\n"
    );
}

#[test]
fn faidx_reverse_complement_marks_sign_and_no() {
    let tmp = tmp_dir("fai-rc-sign-no");
    let fa = tmp.join("ref.fa");
    let sign_out = tmp.join("sign.fa");
    let no_out = tmp.join("no.fa");
    std::fs::write(&fa, ">rc\nACGTMRWSYKVHDBN\n").unwrap();
    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "--mark-strand",
                "sign",
                "-i",
                "-o",
                sign_out.to_str().unwrap(),
                fa.to_str().unwrap(),
                "rc",
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "--mark-strand",
                "no",
                "-i",
                "-o",
                no_out.to_str().unwrap(),
                fa.to_str().unwrap(),
                "rc",
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(sign_out).unwrap(),
        ">rc(-)\nNVHDBMRSWYKACGT\n"
    );
    assert_eq!(
        std::fs::read_to_string(no_out).unwrap(),
        ">rc\nNVHDBMRSWYKACGT\n"
    );
}

#[test]
fn faidx_reverse_complement_marks_custom() {
    let tmp = tmp_dir("fai-rc-custom");
    let fa = tmp.join("ref.fa");
    let out = tmp.join("out.fa");
    std::fs::write(&fa, ">rc\nACGTMRWSYKVHDBN\n").unwrap();
    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "--mark-strand",
                "custom, forward, reverse",
                "-i",
                "-o",
                out.to_str().unwrap(),
                fa.to_str().unwrap(),
                "rc",
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        ">rc reverse\nNVHDBMRSWYKACGT\n"
    );
}

#[test]
fn fqidx_reverse_complement_reverses_quality() {
    let tmp = tmp_dir("fqi-rc-default");
    let fq = tmp.join("reads.fq");
    let out = tmp.join("out.fq");
    std::fs::write(&fq, "@rc\nACGTMRWSYKVHDBN\n+\nabcdefghijklmno\n").unwrap();
    assert_eq!(
        exit_to_u8(fqidx::main(&argv(
            "fqidx",
            &[
                "-i",
                "-o",
                out.to_str().unwrap(),
                fq.to_str().unwrap(),
                "rc"
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@rc/rc\nNVHDBMRSWYKACGT\n+\nonmlkjihgfedcba\n"
    );
}

#[test]
fn fqidx_reverse_complement_marks_no_and_custom() {
    let tmp = tmp_dir("fqi-rc-marks");
    let fq = tmp.join("reads.fq");
    let no_out = tmp.join("no.fq");
    let custom_out = tmp.join("custom.fq");
    std::fs::write(&fq, "@rc\nACGTMRWSYKVHDBN\n+\nabcdefghijklmno\n").unwrap();
    assert_eq!(
        exit_to_u8(fqidx::main(&argv(
            "fqidx",
            &[
                "--mark-strand",
                "no",
                "-i",
                "-o",
                no_out.to_str().unwrap(),
                fq.to_str().unwrap(),
                "rc",
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(fqidx::main(&argv(
            "fqidx",
            &[
                "--mark-strand",
                "custom, forward, reverse",
                "-i",
                "-o",
                custom_out.to_str().unwrap(),
                fq.to_str().unwrap(),
                "rc",
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(no_out).unwrap(),
        "@rc\nNVHDBMRSWYKACGT\n+\nonmlkjihgfedcba\n"
    );
    assert_eq!(
        std::fs::read_to_string(custom_out).unwrap(),
        "@rc reverse\nNVHDBMRSWYKACGT\n+\nonmlkjihgfedcba\n"
    );
}

#[test]
fn import_fastq_to_sam() {
    let tmp = tmp_dir("imp");
    let fq = tmp.join("in.fq");
    std::fs::write(&fq, "@r1\nACGT\n+\n!!!!\n@r2\nTTTT\n+\n####\n").unwrap();
    let out = tmp.join("out.sam");
    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &["-o", out.to_str().unwrap(), fq.to_str().unwrap()]
        ))),
        0
    );
    let body = std::fs::read_to_string(&out).unwrap();
    assert!(body.contains("r1\t4\t*\t0\t0\t*"));
}

#[test]
fn import_fastq_to_bam() {
    let tmp = tmp_dir("imp-bam");
    let fq = tmp.join("in.fq");
    std::fs::write(&fq, "@r1\nACGT\n+\n!!!!\n@r2\nTTTT\n+\n####\n").unwrap();
    let out = tmp.join("out.bam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "-O",
                "bam",
                "-o",
                out.to_str().unwrap(),
                fq.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        htslib_rs::alignment_compat::count_bam_records_from_path(&out).unwrap(),
        2
    );
}

#[test]
fn import_fasta_to_sam() {
    let tmp = tmp_dir("imp-fa");
    let fa = tmp.join("in.fa");
    std::fs::write(&fa, ">r1\nACGT\n>r2\nTT\nTT\n").unwrap();
    let out = tmp.join("out.sam");
    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &["-o", out.to_str().unwrap(), fa.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\nr2\t4\t*\t0\t0\t*\t*\t0\t0\tTTTT\t*\n"
    );
}

#[test]
fn import_accepts_zero_input_and_no_pg() {
    let tmp = tmp_dir("imp-zero");
    let fq = tmp.join("in.fq");
    std::fs::write(&fq, "@r0\nACGT\n+\n!!!!\n").unwrap();
    let out = tmp.join("out.sam");
    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                "-0",
                fq.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "r0\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!\n"
    );
}

#[test]
fn import_paired_fastq_to_sam() {
    let tmp = tmp_dir("imp-paired");
    let r1 = tmp.join("r1.fq");
    let r2 = tmp.join("r2.fq");
    std::fs::write(&r1, "@p\nAC\n+\n!!\n").unwrap();
    std::fs::write(&r2, "@p\nTG\n+\n##\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "-1",
                r1.to_str().unwrap(),
                "-2",
                r2.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "p\t77\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\np\t141\t*\t0\t0\t*\t*\t0\t0\tTG\t##\n"
    );
}

#[test]
fn import_parses_fastq_metadata_options() {
    let tmp = tmp_dir("imp-meta");
    let fq = tmp.join("in.fq");
    std::fs::write(
        &fq,
        "@instrument:run:flowcell:1:1101:100:200:ACGT#0 1:Y:0:barcode\nA\n+\n!\n",
    )
    .unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "-i",
                "-U",
                "--barcode-tag",
                "XB",
                "-o",
                out.to_str().unwrap(),
                fq.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n\
@CO\tReverse with: samtools fastq -n -o paired.fastq -i -U --UMI-tag RX --index-format 'i*i*'\n\
instrument:run:flowcell:1:1101:100:200#0\t589\t*\t0\t0\t*\t*\t0\t0\tA\t!\tXB:Z:barcode\tRX:Z:ACGT\n"
    );
}

#[test]
fn import_umi_writes_reverse_comment() {
    let tmp = tmp_dir("imp-umi-comment");
    let fq = tmp.join("in.fq");
    std::fs::write(&fq, "@r:ACGT\nAC\n+\n!!\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "-U",
                "--UMI-tag",
                "OX",
                "-o",
                out.to_str().unwrap(),
                fq.to_str().unwrap()
            ]
        ))),
        0
    );
    assert!(
        std::fs::read_to_string(&out)
            .unwrap()
            .starts_with("@HD\tVN:1.6\tSO:unsorted\tGO:query\n@CO\tReverse with: samtools fastq -n -o paired.fastq -U --UMI-tag OX\n")
    );
}

#[test]
fn import_preserves_selected_fastq_aux_tags() {
    let tmp = tmp_dir("imp-aux");
    let fq = tmp.join("in.fq");
    std::fs::write(
        &fq,
        "@r1\tXX:i:10\tXZ:i:20\tAA:Z:keep\tAB:Z:drop\nAC\n+\n!!\n",
    )
    .unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "-T",
                "XZ,AA",
                "-o",
                out.to_str().unwrap(),
                fq.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "r1\t4\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\tXZ:i:20\tAA:Z:keep\n"
    );
}

#[test]
fn import_preserves_float_aux_tag_with_upstream_exponent_format() {
    let tmp = tmp_dir("imp-aux-float");
    let fq = tmp.join("in.fq");
    std::fs::write(&fq, "@r1\tFF:f:-1e20\nAC\n+\n!!\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &["-T", "*", "-o", out.to_str().unwrap(), fq.to_str().unwrap()]
        ))),
        0
    );
    assert!(
        std::fs::read_to_string(&out)
            .unwrap()
            .contains("\tFF:f:-1e+20\n")
    );
}

#[test]
fn import_paired_positional_fastq_with_read_group() {
    let tmp = tmp_dir("imp-rg");
    let r1 = tmp.join("r1.fq");
    let r2 = tmp.join("r2.fq");
    std::fs::write(&r1, "@p\nAC\n+\n!!\n").unwrap();
    std::fs::write(&r2, "@p\nTG\n+\n##\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                r1.to_str().unwrap(),
                r2.to_str().unwrap(),
                "-R",
                "rgid",
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n\
@CO\tReverse with: samtools fastq -1 R1.fastq -2 R2.fastq\n\
@RG\tID:rgid\n\
p\t77\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\tRG:Z:rgid\n\
p\t141\t*\t0\t0\t*\t*\t0\t0\tTG\t##\tRG:Z:rgid\n"
    );
}

#[test]
fn import_read_group_line_extracts_id() {
    let tmp = tmp_dir("imp-rg-line");
    let fq = tmp.join("in.fq");
    std::fs::write(&fq, "@r\nAC\n+\n!!\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                "-0",
                fq.to_str().unwrap(),
                "-r",
                "ID:rgid\\tSM:sample",
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n\
@CO\tReverse with: samtools fastq -0 single.fastq\n\
@RG\tID:rgid\tSM:sample\n\
r\t4\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\tRG:Z:rgid\n"
    );
}

#[test]
fn import_repeated_read_group_line_options_accumulate() {
    let tmp = tmp_dir("imp-rg-line-repeat");
    let fq = tmp.join("in.fq");
    std::fs::write(&fq, "@r\nAC\n+\n!!\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                "-0",
                fq.to_str().unwrap(),
                "-r",
                "SM:sample",
                "-r",
                "ID:rgid",
                "-r",
                "LB:lib",
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n\
@CO\tReverse with: samtools fastq -0 single.fastq\n\
@RG\tSM:sample\tID:rgid\tLB:lib\n\
r\t4\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\tRG:Z:rgid\n"
    );
}

#[test]
fn import_read_group_line_takes_precedence_over_id_option() {
    let tmp = tmp_dir("imp-rg-precedence");
    let fq = tmp.join("in.fq");
    std::fs::write(&fq, "@r\nAC\n+\n!!\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                "-R",
                "ignored",
                "-0",
                fq.to_str().unwrap(),
                "-r",
                "ID:rgid",
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n\
@CO\tReverse with: samtools fastq -0 single.fastq\n\
@RG\tID:rgid\n\
r\t4\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\tRG:Z:rgid\n"
    );
}

#[test]
fn import_rejects_read_group_line_without_id() {
    let tmp = tmp_dir("imp-rg-line-invalid");
    let fq = tmp.join("in.fq");
    std::fs::write(&fq, "@r\nAC\n+\n!!\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                "-0",
                fq.to_str().unwrap(),
                "-r",
                "SM:sample",
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        1
    );
    assert!(!out.exists());
}

#[test]
fn import_interleaved_fastq_preserves_aux_tags() {
    let tmp = tmp_dir("imp-interleaved");
    let fq = tmp.join("interleaved.fq");
    std::fs::write(
        &fq,
        "@p/1\tBC:Z:AAA+CCC\nAC\n+\n!!\n@p/2\tBC:Z:AAA+CCC\nTG\n+\n##\n",
    )
    .unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                "-s",
                fq.to_str().unwrap(),
                "-T",
                "",
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n\
@CO\tReverse with: samtools fastq -n -o paired.fastq\n\
p\t77\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\tBC:Z:AAA+CCC\n\
p\t141\t*\t0\t0\t*\t*\t0\t0\tTG\t##\tBC:Z:AAA+CCC\n"
    );
}

#[test]
fn import_interleaved_casava_writes_reverse_comment() {
    let tmp = tmp_dir("imp-interleaved-casava");
    let fq = tmp.join("interleaved.fq");
    std::fs::write(
        &fq,
        "@p/1 1:N:0:AAA+CCC\nAC\n+\n!!\n@p/2 2:N:0:AAA+CCC\nTG\n+\n##\n",
    )
    .unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                fq.to_str().unwrap(),
                "-i",
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert!(std::fs::read_to_string(&out).unwrap().contains(
        "@CO\tReverse with: samtools fastq -n -o paired.fastq -i --index-format 'i*i*'\n"
    ));
}

#[test]
fn import_positional_interleaved_fastq_preserves_aux_tags() {
    let tmp = tmp_dir("imp-pos-interleaved");
    let fq = tmp.join("interleaved.fq");
    std::fs::write(
        &fq,
        "@p/1\tBC:Z:AAA+CCC\nAC\n+\n!!\n@p/2\tBC:Z:AAA+CCC\nTG\n+\n##\n",
    )
    .unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                fq.to_str().unwrap(),
                "-T",
                "",
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n\
@CO\tReverse with: samtools fastq -n -o paired.fastq\n\
p\t77\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\tBC:Z:AAA+CCC\n\
p\t141\t*\t0\t0\t*\t*\t0\t0\tTG\t##\tBC:Z:AAA+CCC\n"
    );
}

#[test]
fn import_paired_fastq_with_index_reads() {
    let tmp = tmp_dir("imp-index");
    let r1 = tmp.join("r1.fq");
    let r2 = tmp.join("r2.fq");
    let i1 = tmp.join("i1.fq");
    let i2 = tmp.join("i2.fq");
    std::fs::write(&r1, "@one\nAC\n+\n!!\n@two\nGT\n+\n##\n").unwrap();
    std::fs::write(&r2, "@one\nTG\n+\n12\n@two\nCA\n+\n34\n").unwrap();
    std::fs::write(&i1, "@one\nAA\n+\nab\n@two\nCC\n+\ncd\n").unwrap();
    std::fs::write(&i2, "@one\nTT\n+\nef\n@two\nGG\n+\ngh\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                "--i1",
                i1.to_str().unwrap(),
                "--i2",
                i2.to_str().unwrap(),
                "--r1",
                r1.to_str().unwrap(),
                "--r2",
                r2.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n\
@CO\tReverse with: samtools fastq --i1 I1.fastq --i2 I2.fastq -1 R1.fastq -2 R2.fastq --index-format=\"i*i*\"\n\
one\t77\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\tBC:Z:AA-TT\tQT:Z:ab ef\n\
one\t141\t*\t0\t0\t*\t*\t0\t0\tTG\t12\n\
two\t77\t*\t0\t0\t*\t*\t0\t0\tGT\t##\tBC:Z:CC-GG\tQT:Z:cd gh\n\
two\t141\t*\t0\t0\t*\t*\t0\t0\tCA\t34\n"
    );
}

#[test]
fn import_index_reads_honor_custom_quality_tag_and_both_reads() {
    let tmp = tmp_dir("imp-index-tags");
    let r1 = tmp.join("r1.fq");
    let r2 = tmp.join("r2.fq");
    let i1 = tmp.join("i1.fq");
    std::fs::write(&r1, "@p\nAC\n+\n!!\n").unwrap();
    std::fs::write(&r2, "@p\nTG\n+\n##\n").unwrap();
    std::fs::write(&i1, "@p\nAA\n+\nab\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                "--i1",
                i1.to_str().unwrap(),
                "-1",
                r1.to_str().unwrap(),
                "-2",
                r2.to_str().unwrap(),
                "--barcode-tag",
                "OX",
                "--quality-tag",
                "BZ",
                "-b",
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n\
@CO\tReverse with: samtools fastq --i1 I1.fastq -1 R1.fastq -2 R2.fastq --index-format=\"i*\"\n\
p\t77\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\tOX:Z:AA\tBZ:Z:ab\n\
p\t141\t*\t0\t0\t*\t*\t0\t0\tTG\t##\tOX:Z:AA\tBZ:Z:ab\n"
    );
}

#[test]
fn import_single_fastq_with_index_reads() {
    let tmp = tmp_dir("imp-single-index");
    let reads = tmp.join("reads.fq");
    let i1 = tmp.join("i1.fq");
    std::fs::write(&reads, "@one\nAC\n+\n!!\n@two\nGT\n+\n##\n").unwrap();
    std::fs::write(&i1, "@one\nAA\n+\nab\n@two\nCC\n+\ncd\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                "--i1",
                i1.to_str().unwrap(),
                "-0",
                reads.to_str().unwrap(),
                "--quality-tag",
                "BZ",
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n\
@CO\tReverse with: samtools fastq --i1 I1.fastq -0 unpaired.fastq --index-format=\"i*\"\n\
one\t4\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\tBC:Z:AA\tBZ:Z:ab\n\
two\t4\t*\t0\t0\t*\t*\t0\t0\tGT\t##\tBC:Z:CC\tBZ:Z:cd\n"
    );
}

#[test]
fn import_interleaved_fastq_with_index_reads() {
    let tmp = tmp_dir("imp-interleaved-index");
    let reads = tmp.join("reads.fq");
    let i1 = tmp.join("i1.fq");
    let i2 = tmp.join("i2.fq");
    std::fs::write(&reads, "@p/1\nAC\n+\n!!\n@p/2\nTG\n+\n##\n").unwrap();
    std::fs::write(&i1, "@p/1\nAA\n+\nab\n@p/2\nCC\n+\ncd\n").unwrap();
    std::fs::write(&i2, "@p/1\nTT\n+\nef\n@p/2\nGG\n+\ngh\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                "--i1",
                i1.to_str().unwrap(),
                "--i2",
                i2.to_str().unwrap(),
                "-s",
                reads.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n\
@CO\tReverse with: samtools fastq --i1 I1.fastq --i2 I2.fastq -n -o paired.fastq --index-format=\"i*i*\"\n\
p\t77\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\tBC:Z:AA-TT\tQT:Z:ab ef\n\
p\t141\t*\t0\t0\t*\t*\t0\t0\tTG\t##\tBC:Z:CC-GG\tQT:Z:cd gh\n"
    );
}

#[test]
fn import_positional_interleaved_fastq_with_index_reads() {
    let tmp = tmp_dir("imp-pos-interleaved-index");
    let reads = tmp.join("reads.fq");
    let i1 = tmp.join("i1.fq");
    std::fs::write(&reads, "@p/1\nAC\n+\n!!\n@p/2\nTG\n+\n##\n").unwrap();
    std::fs::write(&i1, "@p/1\nAA\n+\nab\n@p/2\nCC\n+\ncd\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--no-PG",
                "--i1",
                i1.to_str().unwrap(),
                reads.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n\
@CO\tReverse with: samtools fastq --i1 I1.fastq -n -o paired.fastq --index-format=\"i*\"\n\
p\t77\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\tBC:Z:AA\tQT:Z:ab\n\
p\t141\t*\t0\t0\t*\t*\t0\t0\tTG\t##\tBC:Z:CC\tQT:Z:cd\n"
    );
}

#[test]
fn bedcov_basic() {
    let tmp = tmp_dir("bedcov");
    let bam = tmp.join("in.bam");
    std::fs::copy(sample_bam(), &bam).unwrap();
    assert_eq!(
        exit_to_u8(index::main(&argv("index", &[bam.to_str().unwrap()]))),
        0
    );
    let bed = tmp.join("regions.bed");
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
fn rmdup_succeeds() {
    let tmp = tmp_dir("rmdup");
    let out = tmp.join("dedup.bam");
    assert_eq!(
        exit_to_u8(rmdup::main(&argv(
            "rmdup",
            &[sample_bam().to_str().unwrap(), out.to_str().unwrap()]
        ))),
        0
    );
    assert!(out.exists());
}

#[test]
fn split_by_rg() {
    let tmp = tmp_dir("split");
    let bam = tmp.join("in.bam");
    std::fs::copy(sample_bam(), &bam).unwrap();
    let tmpl = tmp.join("out.%#.%.");
    let unk = tmp.join("unk.bam");
    assert_eq!(
        exit_to_u8(split::main(&argv(
            "split",
            &[
                "-f",
                tmpl.to_str().unwrap(),
                "-u",
                unk.to_str().unwrap(),
                bam.to_str().unwrap()
            ]
        ))),
        0
    );
}
