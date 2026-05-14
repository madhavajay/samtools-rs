//! Smoke-tests for `cat`, `reheader`, `fastq`, `samples`, `idxstats`,
//! `flagstat`, `index`, `faidx`, `import`, `bedcov`, `rmdup`, `split`.

use std::ffi::OsString;
use std::io::{BufReader, Cursor};
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::Mutex;

use samtools_rs::commands::{
    bedcov, cat, faidx, fastq, fixmate, flagstat, fqidx, idxstats, import, index, reheader, reset,
    rmdup, samples, split,
};
use samtools_rs::header_text;
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

fn write_bam_from_sam_text(path: &std::path::Path, text: &str) {
    let mut reader = htslib_rs::sam::io::Reader::new(BufReader::new(Cursor::new(text.as_bytes())));
    let header = reader.read_header().unwrap();
    let mut writer = htslib_rs::bam::io::Writer::new(std::fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();

    for result in reader.records() {
        let record = result.unwrap();
        use htslib_rs::sam::alignment::io::Write as _;
        writer.write_alignment_record(&header, &record).unwrap();
    }
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
fn flagstat_cram_uses_top_level_reference() {
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
                "flagstat",
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );
}

#[test]
fn flagstat_cram_without_reference_fails_cleanly() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let cram = htslib_fixtures_dir().join("range.cram");

    assert_ne!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &["flagstat", cram.to_str().unwrap()],
        ))),
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
fn idxstats_cram_uses_top_level_reference() {
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
                "idxstats",
                cram.to_str().unwrap(),
            ],
        ))),
        0
    );
}

#[test]
fn idxstats_cram_without_reference_fails_cleanly() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let cram = htslib_fixtures_dir().join("range.cram");

    assert_ne!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &["idxstats", cram.to_str().unwrap()],
        ))),
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
fn cat_adds_pg_by_default() {
    let tmp = tmp_dir("cat-pg");
    let out = tmp.join("cat.bam");
    let p = sample_bam();
    assert_eq!(
        exit_to_u8(cat::main(&argv(
            "cat",
            &[p.to_str().unwrap(), "-o", out.to_str().unwrap()]
        ))),
        0
    );

    let header = header_text::read_raw_header_text(&out).unwrap();
    assert!(header.contains("\tPN:samtools\tVN:"));
    assert!(header.contains("\tCL:cat "));
}

#[test]
fn cat_no_pg_suppresses_pg() {
    let tmp = tmp_dir("cat-no-pg");
    let out = tmp.join("cat.bam");
    let p = sample_bam();
    assert_eq!(
        exit_to_u8(cat::main(&argv(
            "cat",
            &["--no-PG", p.to_str().unwrap(), "-o", out.to_str().unwrap()]
        ))),
        0
    );

    let header = header_text::read_raw_header_text(&out).unwrap();
    assert!(!header.contains("\tCL:cat "));
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
fn fastq_filters_by_include_any_long_flag_alias() {
    let tmp = tmp_dir("fastq-include-any");
    let sam = tmp.join("in.sam");
    let out = tmp.join("reads.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "read1\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\n",
            "read2\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\n",
            "unpaired\t0\tchr1\t1\t60\t4M\t*\t0\t0\tNNNN\t$$$$\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "--include-flags",
                "64",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@read1/1\nACGT\n+\n!!!!\n"
    );
}

#[test]
fn fasta_filters_by_include_any_long_flag_alias() {
    let tmp = tmp_dir("fasta-include-any");
    let sam = tmp.join("in.sam");
    let out = tmp.join("reads.fa");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "read1\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\n",
            "read2\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\n",
            "unpaired\t0\tchr1\t1\t60\t4M\t*\t0\t0\tNNNN\t$$$$\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fasta",
            &[
                "--include-flags",
                "64",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(std::fs::read_to_string(out).unwrap(), ">read1/1\nACGT\n");
}

#[test]
fn fastq_excludes_secondary_and_supplementary_by_default() {
    let tmp = tmp_dir("fastq-default-exclude");
    let sam = tmp.join("in.sam");
    let default_out = tmp.join("default.fq");
    let include_out = tmp.join("include.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "primary\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "secondary\t256\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "supplementary\t2048\tchr1\t1\t60\t4M\t*\t0\t0\tNNNN\t$$$$\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &["-o", default_out.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-F",
                "0",
                "-o",
                include_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(default_out).unwrap(),
        "@primary\nACGT\n+\n!!!!\n"
    );
    assert_eq!(
        std::fs::read_to_string(include_out).unwrap(),
        "@primary\nACGT\n+\n!!!!\n@secondary\nTGCA\n+\n####\n@supplementary\nNNNN\n+\n$$$$\n"
    );
}

#[test]
fn fastq_single_sam_path_filters_by_aux_tag_value() {
    let tmp = tmp_dir("fastq-aux-filter-value");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "keep\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tNM:i:0\tMD:Z:4\n",
            "drop\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\tNM:i:1\tMD:Z:3A\n",
            "missing\t0\tchr1\t1\t60\t4M\t*\t0\t0\tNNNN\t$$$$\tMD:Z:4\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-d",
                "NM:0",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.contains("@keep\n"));
    assert!(!text.contains("@drop\n"));
    assert!(!text.contains("@missing\n"));
}

#[test]
fn fastq_single_sam_path_filters_by_aux_tag_presence() {
    let tmp = tmp_dir("fastq-aux-filter-presence");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "has\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tBC:Z:AAAA\n",
            "missing\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-d",
                "BC",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.contains("@has\n"));
    assert!(!text.contains("@missing\n"));
}

#[test]
fn fastq_single_sam_path_filters_aux_values_from_file() {
    let tmp = tmp_dir("fastq-aux-filter-file");
    let sam = tmp.join("in.sam");
    let values = tmp.join("values.txt");
    let out = tmp.join("out.fq");
    std::fs::write(&values, "AAAA\nCCCC\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "keep1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tBC:Z:AAAA\n",
            "drop\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\tBC:Z:GGGG\n",
            "keep2\t0\tchr1\t1\t60\t4M\t*\t0\t0\tNNNN\t$$$$\tBC:Z:CCCC\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-D",
                &format!("BC:{}", values.display()),
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.contains("@keep1\n"));
    assert!(text.contains("@keep2\n"));
    assert!(!text.contains("@drop\n"));
}

#[test]
fn fastq_single_bam_path_preserves_selected_aux_tags() {
    let tmp = tmp_dir("fastq-bam-aux");
    let bam = tmp.join("in.bam");
    let out = tmp.join("out.fq");
    write_bam_from_sam_text(
        &bam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:rg1\tBC:Z:ACTG\tNM:i:1\n",
        ),
    );

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-T",
                "RG,BC,NM",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@r1\tRG:Z:rg1\tBC:Z:ACTG\tNM:i:1\nACGT\n+\n!!!!\n"
    );
}

#[test]
fn fastq_single_bam_path_filters_by_aux_tag_value() {
    let tmp = tmp_dir("fastq-bam-aux-filter");
    let bam = tmp.join("in.bam");
    let out = tmp.join("out.fq");
    write_bam_from_sam_text(
        &bam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "keep\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tNM:i:0\n",
            "drop\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\tNM:i:1\n",
        ),
    );

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-d",
                "NM:0",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.contains("@keep\n"));
    assert!(!text.contains("@drop\n"));
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
fn fastq_split_sam_path_preserves_selected_aux_tags() {
    let tmp = tmp_dir("fastq-split-aux");
    let sam = tmp.join("in.sam");
    let r1 = tmp.join("r1.fq");
    let r2 = tmp.join("r2.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "pair\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\tRG:Z:rg1\tBC:Z:ACTG\tQT:Z:!!!!\tNM:i:1\n",
            "pair\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\tRG:Z:rg1\tBC:Z:TGCA\tQT:Z:####\tNM:i:2\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-N",
                "-T",
                "RG,BC,QT",
                "-1",
                r1.to_str().unwrap(),
                "-2",
                r2.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(r1).unwrap(),
        "@pair/1\tRG:Z:rg1\tBC:Z:ACTG\tQT:Z:!!!!\nACGT\n+\n!!!!\n"
    );
    assert_eq!(
        std::fs::read_to_string(r2).unwrap(),
        "@pair/2\tRG:Z:rg1\tBC:Z:TGCA\tQT:Z:####\nTGCA\n+\n####\n"
    );
}

#[test]
fn fastq_split_sam_path_t_copies_default_aux_tags() {
    let tmp = tmp_dir("fastq-split-t");
    let sam = tmp.join("in.sam");
    let r1 = tmp.join("r1.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "pair\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\tRG:Z:rg1\tBC:Z:ACTG\tQT:Z:!!!!\tNM:i:1\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-N",
                "-t",
                "-1",
                r1.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(r1).unwrap(),
        "@pair/1\tRG:Z:rg1\tBC:Z:ACTG\tQT:Z:!!!!\nACGT\n+\n!!!!\n"
    );
}

#[test]
fn fastq_split_bam_path_preserves_selected_aux_tags() {
    let tmp = tmp_dir("fastq-split-bam-aux");
    let bam = tmp.join("in.bam");
    let r1 = tmp.join("r1.fq");
    let r2 = tmp.join("r2.fq");
    write_bam_from_sam_text(
        &bam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "pair\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\tRG:Z:rg1\tBC:Z:ACTG\tQT:Z:!!!!\tNM:i:1\n",
            "pair\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\tRG:Z:rg1\tBC:Z:TGCA\tQT:Z:####\tNM:i:2\n",
        ),
    );

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-N",
                "-T",
                "RG,BC,NM",
                "-1",
                r1.to_str().unwrap(),
                "-2",
                r2.to_str().unwrap(),
                bam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(r1).unwrap(),
        "@pair/1\tRG:Z:rg1\tBC:Z:ACTG\tNM:i:1\nACGT\n+\n!!!!\n"
    );
    assert_eq!(
        std::fs::read_to_string(r2).unwrap(),
        "@pair/2\tRG:Z:rg1\tBC:Z:TGCA\tNM:i:2\nTGCA\n+\n####\n"
    );
}

#[test]
fn fastq_split_bam_path_filters_by_aux_tag_value() {
    let tmp = tmp_dir("fastq-split-bam-filter");
    let bam = tmp.join("in.bam");
    let r1 = tmp.join("r1.fq");
    let r2 = tmp.join("r2.fq");
    let singleton = tmp.join("single.fq");
    write_bam_from_sam_text(
        &bam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "keep\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\tBC:Z:KEEP\n",
            "keep\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\tBC:Z:KEEP\n",
            "drop\t65\tchr1\t1\t60\t4M\t=\t5\t8\tAAAA\t!!!!\tBC:Z:DROP\n",
            "drop\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTTTT\t####\tBC:Z:DROP\n",
            "solo\t4\t*\t0\t0\t*\t*\t0\t0\tNNNN\t$$$$\tBC:Z:KEEP\n",
        ),
    );

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-d",
                "BC:KEEP",
                "-1",
                r1.to_str().unwrap(),
                "-2",
                r2.to_str().unwrap(),
                "-s",
                singleton.to_str().unwrap(),
                bam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(r1).unwrap(),
        "@keep\nACGT\n+\n!!!!\n"
    );
    assert_eq!(
        std::fs::read_to_string(r2).unwrap(),
        "@keep\nTGCA\n+\n####\n"
    );
    assert_eq!(
        std::fs::read_to_string(singleton).unwrap(),
        "@solo\nNNNN\n+\n$$$$\n"
    );
}

#[test]
fn fastq_single_sam_path_preserves_selected_aux_tags() {
    let tmp = tmp_dir("fastq-single-aux");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tba:A:x\tbb:i:7\tbc:Z:text\tNM:i:1\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-T",
                "ba,bb,bc",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@r1\tba:A:x\tbb:i:7\tbc:Z:text\nACGT\n+\n!!!!\n"
    );
}

#[test]
fn fastq_single_sam_path_empty_t_copies_all_aux_tags() {
    let tmp = tmp_dir("fastq-single-all-aux-empty");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tMD:Z:4\tNM:i:0\tRG:Z:rg1\tBC:Z:ACTG\tba:B:c,-1,0,1\tbb:B:C,0,127,255\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &["-T", "", "-o", out.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@r1\tMD:Z:4\tNM:i:0\tRG:Z:rg1\tBC:Z:ACTG\tba:B:c,-1,0,1\tbb:B:C,0,127,255\nACGT\n+\n!!!!\n"
    );
}

#[test]
fn fastq_single_sam_path_star_t_copies_all_aux_tags_after_t() {
    let tmp = tmp_dir("fastq-single-all-aux-star");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:rg1\tBC:Z:ACTG\tQT:Z:!!!!\tMD:Z:4\tia:i:40000\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-t",
                "-T",
                "*",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@r1\tRG:Z:rg1\tBC:Z:ACTG\tQT:Z:!!!!\tMD:Z:4\tia:i:40000\nACGT\n+\n!!!!\n"
    );
}

#[test]
fn fastq_single_sam_aux_path_reverse_complements_reverse_reads() {
    let tmp = tmp_dir("fastq-single-aux-reverse");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t16\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!#$%\tNM:i:0\n",
            "r2\t16\tchr1\t1\t60\t15M\t*\t0\t0\tACGTMRWSYKVHDBN\t0123456789abcd!\tNM:i:1\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-T",
                "*",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@r1\tNM:i:0\nACGT\n+\n%$#!\n@r2\tNM:i:1\nNVHDBMRSWYKACGT\n+\n!dcba9876543210\n"
    );
}

#[test]
fn fixmate_sam_input_fills_mate_fields_to_sam_output() {
    let tmp = tmp_dir("fixmate-sam");
    let sam = tmp.join("in.sam");
    let out = tmp.join("fixed.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:queryname\n",
            "@SQ\tSN:chr1\tLN:16\n",
            "pair\t65\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "pair\t129\tchr1\t5\t60\t4M\t*\t0\t0\tTGCA\t####\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fixmate::main(&argv(
            "fixmate",
            &[
                "--output-fmt",
                "sam",
                sam.to_str().unwrap(),
                out.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    let records: Vec<Vec<_>> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0][0], "pair");
    assert_eq!(records[0][1], "65");
    assert_eq!(records[0][6], "=");
    assert_eq!(records[0][7], "5");
    assert_eq!(records[1][1], "129");
    assert_eq!(records[1][6], "=");
    assert_eq!(records[1][7], "1");
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
fn faidx_builds_index_for_bgzf_input_and_writes_gzi() {
    let tmp = tmp_dir("fai-bgzf");
    let plain = b">chr1\nACGTACGT\n>chr2\nTTTTCCCC\n";
    let encoded = htslib_rs::bgzf_compat::write_all_with_kind(
        plain,
        htslib_rs::bgzf_compat::CompressionKind::Bgzf,
    )
    .unwrap();
    let fa = tmp.join("ref.fa.gz");
    std::fs::write(&fa, encoded).unwrap();

    assert_eq!(
        exit_to_u8(faidx::main(&argv("faidx", &[fa.to_str().unwrap()]))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(tmp.join("ref.fa.gz.fai")).unwrap(),
        "chr1\t8\t6\t8\t9\nchr2\t8\t21\t8\t9\n"
    );
    assert!(tmp.join("ref.fa.gz.gzi").exists());
}

#[test]
fn faidx_extracts_from_bgzf_input_and_writes_bgzf_output() {
    let tmp = tmp_dir("fai-bgzf-region");
    let plain = b">chr1\nACGTACGTACGT\n>chr2\nTTTTCCCC\n";
    let encoded = htslib_rs::bgzf_compat::write_all_with_kind(
        plain,
        htslib_rs::bgzf_compat::CompressionKind::Bgzf,
    )
    .unwrap();
    let fa = tmp.join("ref.fa.gz");
    let out = tmp.join("out.fa.gz");
    std::fs::write(&fa, encoded).unwrap();

    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "--length",
                "4",
                "-o",
                out.to_str().unwrap(),
                fa.to_str().unwrap(),
                "chr1:3-10",
            ]
        ))),
        0
    );

    let decoded = htslib_rs::bgzf_compat::read_auto(&std::fs::read(out).unwrap()).unwrap();
    assert_eq!(
        String::from_utf8(decoded).unwrap(),
        ">chr1:3-10 length: 8\nGTAC\nGTAC\n"
    );
}

#[test]
fn faidx_accepts_equals_output_format_option_for_bgzf_output() {
    let tmp = tmp_dir("fai-output-fmt-opt");
    let fa = tmp.join("ref.fa");
    let out = tmp.join("out.fa.gz");
    std::fs::write(&fa, ">chr1\nACGTACGTACGT\n").unwrap();

    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "--length",
                "4",
                "-o",
                out.to_str().unwrap(),
                fa.to_str().unwrap(),
                "chr1:3-10",
                "--output-fmt-opt=level=4",
            ]
        ))),
        0
    );

    let decoded = htslib_rs::bgzf_compat::read_auto(&std::fs::read(out).unwrap()).unwrap();
    assert_eq!(
        String::from_utf8(decoded).unwrap(),
        ">chr1:3-10 length: 8\nGTAC\nGTAC\n"
    );
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
fn faidx_write_index_uses_default_sixty_base_output_lines() {
    let tmp = tmp_dir("fai-write-index");
    let fa = tmp.join("ref.fa");
    let regions = tmp.join("regions.txt");
    let out = tmp.join("out.fa");
    let seq = "A".repeat(61);
    std::fs::write(&fa, format!(">chr1\n{seq}\n")).unwrap();
    std::fs::write(&regions, "chr1\n").unwrap();

    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "--write-index",
                "-r",
                regions.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                fa.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        format!(">chr1\n{}\nA\n", "A".repeat(60))
    );
    assert_eq!(
        std::fs::read_to_string(tmp.join("out.fa.fai")).unwrap(),
        "chr1\t61\t6\t60\t61\n"
    );
}

#[test]
fn faidx_out_of_range_region_exits_successfully_with_empty_output() {
    let tmp = tmp_dir("fai-zero-region");
    let fa = tmp.join("ref.fa");
    let out = tmp.join("out.fa");
    std::fs::write(&fa, ">chr1\nACGTACGT\n").unwrap();

    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "-o",
                out.to_str().unwrap(),
                fa.to_str().unwrap(),
                "chr1:100-105",
            ]
        ))),
        0
    );

    assert_eq!(std::fs::read_to_string(out).unwrap(), "");
}

#[test]
fn faidx_truncated_region_exits_successfully_with_clamped_output() {
    let tmp = tmp_dir("fai-truncated-region");
    let fa = tmp.join("ref.fa");
    let out = tmp.join("out.fa");
    std::fs::write(&fa, ">chr1\nACGTACGT\n").unwrap();

    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "-o",
                out.to_str().unwrap(),
                fa.to_str().unwrap(),
                "chr1:6-12",
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        ">chr1:6-12 length: 3\nCGT\n"
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
fn fqidx_builds_index_for_bgzf_input_and_writes_gzi() {
    let tmp = tmp_dir("fqi-bgzf");
    let plain = b"@r1\nACGTACGT\n+\nabcdefgh\n";
    let encoded = htslib_rs::bgzf_compat::write_all_with_kind(
        plain,
        htslib_rs::bgzf_compat::CompressionKind::Bgzf,
    )
    .unwrap();
    let fq = tmp.join("reads.fq.gz");
    std::fs::write(&fq, encoded).unwrap();

    assert_eq!(
        exit_to_u8(fqidx::main(&argv("fqidx", &[fq.to_str().unwrap()]))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(tmp.join("reads.fq.gz.fai")).unwrap(),
        "r1\t8\t4\t8\t9\t15\n"
    );
    assert!(tmp.join("reads.fq.gz.gzi").exists());
}

#[test]
fn fqidx_extracts_from_bgzf_input_and_writes_bgzf_output() {
    let tmp = tmp_dir("fqi-bgzf-region");
    let plain = b"@r1\nACGTACGT\n+\nabcdefgh\n";
    let encoded = htslib_rs::bgzf_compat::write_all_with_kind(
        plain,
        htslib_rs::bgzf_compat::CompressionKind::Bgzf,
    )
    .unwrap();
    let fq = tmp.join("reads.fq.gz");
    let out = tmp.join("out.fq.gz");
    std::fs::write(&fq, encoded).unwrap();

    assert_eq!(
        exit_to_u8(fqidx::main(&argv(
            "fqidx",
            &[
                "--length",
                "4",
                "-o",
                out.to_str().unwrap(),
                fq.to_str().unwrap(),
                "r1:2-7",
            ]
        ))),
        0
    );

    let decoded = htslib_rs::bgzf_compat::read_auto(&std::fs::read(out).unwrap()).unwrap();
    assert_eq!(
        String::from_utf8(decoded).unwrap(),
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
fn rmdup_sam_input_keeps_best_single_end_duplicate() {
    let tmp = tmp_dir("rmdup-sam");
    let sam = tmp.join("in.sam");
    let out = tmp.join("dedup.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "low\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "reverse\t16\tchr1\t1\t30\t4M\t*\t0\t0\tCCCC\t$$$$\n",
            "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tNNNN\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(rmdup::main(&argv(
            "rmdup",
            &[sam.to_str().unwrap(), out.to_str().unwrap()]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(!text.contains("\nlow\t"));
    assert!(text.contains("\nhigh\t"));
    assert!(text.contains("\nreverse\t"));
    assert!(text.contains("\nunmapped\t"));
}

#[test]
fn reset_sam_input_clears_alignment_fields_and_default_tags() {
    let tmp = tmp_dir("reset-sam");
    let sam = tmp.join("in.sam");
    let out = tmp.join("reset.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t99\tchr1\t2\t60\t4M\t=\t6\t8\tACGT\t!!!!\tNM:i:1\tMD:Z:3A\tRG:Z:g1\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(reset::main(&argv(
            "reset",
            &[
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    let record = text.lines().find(|line| !line.starts_with('@')).unwrap();
    let fields: Vec<_> = record.split('\t').collect();
    assert_eq!(fields[1], "77");
    assert_eq!(fields[2], "*");
    assert_eq!(fields[3], "0");
    assert_eq!(fields[4], "0");
    assert_eq!(fields[5], "*");
    assert_eq!(fields[6], "*");
    assert_eq!(fields[7], "0");
    assert_eq!(fields[8], "0");
    assert!(!record.contains("\tNM:i:"));
    assert!(!record.contains("\tMD:Z:"));
    assert!(record.contains("\tRG:Z:g1"));
}

#[test]
fn reset_dupflag_preserves_duplicate_and_restores_reverse_sequence() {
    let tmp = tmp_dir("reset-dupflag-reverse");
    let sam = fixtures_dir().join("reset").join("seq.sam");
    let out = tmp.join("reset.sam");

    assert_eq!(
        exit_to_u8(reset::main(&argv(
            "reset",
            &[
                "--dupflag",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    let actual = std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('@'))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let expected =
        std::fs::read_to_string(fixtures_dir().join("reset").join("output.flg.1.expected"))
            .unwrap();

    assert_eq!(actual, expected);
}

#[test]
fn reset_no_rg_removes_read_group_headers_and_tags() {
    let tmp = tmp_dir("reset-no-rg");
    let sam = tmp.join("in.sam");
    let out = tmp.join("reset.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:g1\tSM:s1\n",
            "@RG\tID:g2\tSM:s2\n",
            "r1\t99\tchr1\t2\t60\t4M\t=\t6\t8\tACGT\t!!!!\tRG:Z:g1\tNM:i:1\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(reset::main(&argv(
            "reset",
            &[
                "--no-RG",
                "--keep-tag",
                "RG",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(!text.contains("\n@RG\t"));
    assert!(!text.contains("\tRG:Z:"));
}

#[test]
fn reset_no_pg_removes_program_headers() {
    let tmp = tmp_dir("reset-no-pg");
    let sam = tmp.join("in.sam");
    let out = tmp.join("reset.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@PG\tID:aligner\tPN:aligner\n",
            "@PG\tID:post\tPN:post\tPP:aligner\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(reset::main(&argv(
            "reset",
            &[
                "--no-PG",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(!text.contains("\n@PG\t"));
}

#[test]
fn reset_reject_pg_removes_program_chain_from_id() {
    let tmp = tmp_dir("reset-reject-pg");
    let sam = tmp.join("in.sam");
    let out = tmp.join("reset.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@PG\tID:bwa_index\tPN:bwa\n",
            "@PG\tID:bwa_aln\tPN:bwa\tPP:bwa_index\n",
            "@PG\tID:qc\tPN:qc\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(reset::main(&argv(
            "reset",
            &[
                "--reject-PG",
                "bwa_index",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(!text.contains("@PG\tID:bwa_index"));
    assert!(!text.contains("@PG\tID:bwa_aln"));
    assert!(text.contains("@PG\tID:qc"));
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

#[test]
fn split_sam_input_by_rg_to_sam_outputs() {
    let tmp = tmp_dir("split-sam");
    let sam = tmp.join("in.sam");
    let tmpl = tmp.join("out.%!.%.");
    let unk = tmp.join("unknown.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:g1\tSM:s1\n",
            "@RG\tID:g2\tSM:s2\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:g1\n",
            "r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\tRG:Z:g2\n",
            "r3\t0\tchr1\t3\t60\t4M\t*\t0\t0\tCCCC\t$$$$\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(split::main(&argv(
            "split",
            &[
                "--output-fmt",
                "sam",
                "-f",
                tmpl.to_str().unwrap(),
                "-u",
                unk.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let g1 = std::fs::read_to_string(tmp.join("out.g1.sam")).unwrap();
    let g2 = std::fs::read_to_string(tmp.join("out.g2.sam")).unwrap();
    let unknown = std::fs::read_to_string(unk).unwrap();
    assert!(g1.contains("@RG\tID:g1"));
    assert!(!g1.contains("@RG\tID:g2"));
    assert!(g2.contains("@RG\tID:g2"));
    assert!(!g2.contains("@RG\tID:g1"));
    assert!(g1.lines().any(|line| line.starts_with("r1\t")));
    assert!(g2.lines().any(|line| line.starts_with("r2\t")));
    assert!(unknown.lines().any(|line| line.starts_with("r3\t")));
}

#[test]
fn split_by_explicit_aux_tag_creates_outputs_on_demand() {
    let tmp = tmp_dir("split-aux-tag");
    let sam = tmp.join("in.sam");
    let tmpl = tmp.join("%*_%#_%!.%.");
    let unk = tmp.join("unknown.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tBC:Z:AA\n",
            "r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\tBC:Z:BB\n",
            "r3\t0\tchr1\t3\t60\t4M\t*\t0\t0\tCCCC\t$$$$\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(split::main(&argv(
            "split",
            &[
                "--output-fmt",
                "sam",
                "-d",
                "BC",
                "-f",
                tmpl.to_str().unwrap(),
                "-u",
                unk.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let aa = std::fs::read_to_string(tmp.join("in_0_AA.sam")).unwrap();
    let bb = std::fs::read_to_string(tmp.join("in_1_BB.sam")).unwrap();
    let unknown = std::fs::read_to_string(unk).unwrap();
    assert!(aa.lines().any(|line| line.starts_with("r1\t")));
    assert!(bb.lines().any(|line| line.starts_with("r2\t")));
    assert!(unknown.lines().any(|line| line.starts_with("r3\t")));
}

#[test]
fn split_by_explicit_integer_tag_honors_max_split() {
    let tmp = tmp_dir("split-aux-int-max");
    let sam = tmp.join("in.sam");
    let tmpl = tmp.join("tag_%!.%.");
    let unk = tmp.join("unknown.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tZZ:i:7\n",
            "r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\tZZ:i:8\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(split::main(&argv(
            "split",
            &[
                "--output-fmt",
                "sam",
                "-d",
                "ZZ",
                "-M",
                "1",
                "-f",
                tmpl.to_str().unwrap(),
                "-u",
                unk.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let first = std::fs::read_to_string(tmp.join("tag_7.sam")).unwrap();
    let unknown = std::fs::read_to_string(unk).unwrap();
    assert!(first.lines().any(|line| line.starts_with("r1\t")));
    assert!(unknown.lines().any(|line| line.starts_with("r2\t")));
    assert!(!tmp.join("tag_8.sam").exists());
}

#[test]
fn split_unaccounted_output_can_use_header_override() {
    let tmp = tmp_dir("split-unaccounted-header");
    let sam = tmp.join("in.sam");
    let header = tmp.join("unaccounted-header.sam");
    let tmpl = tmp.join("out.%!.%.");
    let unk = tmp.join("unknown.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:g1\tSM:s1\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:g1\n",
            "r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\n",
        ),
    )
    .unwrap();
    std::fs::write(
        &header,
        concat!(
            "@HD\tVN:1.6\tSO:unknown\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@CO\tcustom unaccounted header\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(split::main(&argv(
            "split",
            &[
                "--output-fmt",
                "sam",
                "-f",
                tmpl.to_str().unwrap(),
                "-u",
                unk.to_str().unwrap(),
                "-h",
                header.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let unknown = std::fs::read_to_string(unk).unwrap();
    assert!(unknown.contains("@CO\tcustom unaccounted header\n"));
    assert!(unknown.lines().any(|line| line.starts_with("r2\t")));
}

#[test]
fn split_explicit_rg_tag_adds_header_for_unknown_read_group() {
    let tmp = tmp_dir("split-explicit-rg-unknown");
    let sam = tmp.join("in.sam");
    let tmpl = tmp.join("out.%!.%.");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:g1\tSM:s1\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:g2\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(split::main(&argv(
            "split",
            &[
                "--output-fmt",
                "sam",
                "-d",
                "RG",
                "-f",
                tmpl.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let g2 = std::fs::read_to_string(tmp.join("out.g2.sam")).unwrap();
    assert!(g2.contains("@RG\tID:g2\n"));
    assert!(!g2.contains("@RG\tID:g1"));
    assert!(g2.lines().any(|line| line.starts_with("r1\t")));
}

#[test]
fn split_adds_pg_by_default() {
    let tmp = tmp_dir("split-default-pg");
    let sam = tmp.join("in.sam");
    let tmpl = tmp.join("out.%!.%.");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:g1\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:g1\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(split::main(&argv(
            "split",
            &[
                "--output-fmt",
                "sam",
                "-f",
                tmpl.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let g1 = std::fs::read_to_string(tmp.join("out.g1.sam")).unwrap();
    assert!(g1.contains("@PG\tID:samtools\tPN:samtools\tVN:"));
    assert!(g1.contains("\tCL:split "));
    assert!(g1.lines().any(|line| line.starts_with("r1\t")));
}

#[test]
fn split_write_index_builds_bai_for_bam_outputs() {
    let tmp = tmp_dir("split-write-index");
    let sam = tmp.join("in.sam");
    let tmpl = tmp.join("out.%!.%.");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:g1\n",
            "@RG\tID:g2\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:g1\n",
            "r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\tRG:Z:g2\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(split::main(&argv(
            "split",
            &[
                "--write-index",
                "-f",
                tmpl.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    assert!(tmp.join("out.g1.bam").exists());
    assert!(tmp.join("out.g1.bam.bai").exists());
    assert!(tmp.join("out.g2.bam").exists());
    assert!(tmp.join("out.g2.bam.bai").exists());
}

#[test]
fn split_accepts_no_pg_option() {
    let tmp = tmp_dir("split-no-pg");
    let sam = tmp.join("in.sam");
    let tmpl = tmp.join("out.%!.%.");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:g1\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:g1\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(split::main(&argv(
            "split",
            &[
                "--output-fmt",
                "sam",
                "--no-PG",
                "-f",
                tmpl.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let g1 = std::fs::read_to_string(tmp.join("out.g1.sam")).unwrap();
    assert!(!g1.contains("@PG\tID:samtools"));
    assert!(g1.lines().any(|line| line.starts_with("r1\t")));
}
