//! Smoke-tests for `cat`, `reheader`, `fastq`, `samples`, `idxstats`,
//! `flagstat`, `index`, `faidx`, `import`, `bedcov`, `rmdup`, `split`.

use std::ffi::OsString;
use std::io::{BufReader, Cursor};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Mutex;

use htslib_rs::bam;
use htslib_rs::sam;
use samtools_rs::commands::{
    addreplacerg, bedcov, calmd, cat, checksum, consensus, cram_size, depad, faidx, fastq, fixmate,
    flagstat, fqidx, idxstats, import, index, mpileup, reference, reheader, reset, rmdup, samples,
    sort, split, view,
};
use samtools_rs::header_text;
use samtools_rs::run as samtools_run;
use samtools_rs::sam_global::{SamGlobalArgs, set_current_global_args};

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

fn non_header_lines(text: &str) -> Vec<String> {
    text.lines()
        .filter(|line| !line.starts_with('@'))
        .map(str::to_string)
        .collect()
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

fn fasta_region(path: &Path, name: &str, start: usize, end: usize) -> String {
    let text = std::fs::read_to_string(path).unwrap();
    let mut current_name = None;
    let mut seq = String::new();

    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            current_name = Some(header.split_whitespace().next().unwrap_or(""));
            continue;
        }

        if current_name == Some(name) {
            seq.push_str(line);
        }
    }

    let slice = &seq[start - 1..end];
    let mut out = format!(">{name}:{start}-{end}\n");
    for chunk in slice.as_bytes().chunks(60) {
        out.push_str(std::str::from_utf8(chunk).unwrap());
        out.push('\n');
    }
    out
}

fn build_reference_embed_ref_cram(tmp: &Path) -> PathBuf {
    let d = fixtures_dir();
    let sam = d.join("dat/mpileup.1.sam");
    let refa = d.join("dat/mpileup.ref.fa");
    let cram = tmp.join("mpileup.1.tmp.cram");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "view",
                "--no-PG",
                "-e",
                "pos<1000||pos>1200",
                "-O",
                "cram,embed_ref=1",
                "-T",
                refa.to_str().unwrap(),
                "-o",
                cram.to_str().unwrap(),
                sam.to_str().unwrap(),
            ],
        ))),
        0
    );

    cram
}

fn argv(name: &str, rest: &[&str]) -> Vec<OsString> {
    std::iter::once(OsString::from(name))
        .chain(rest.iter().map(OsString::from))
        .collect()
}

fn sample_bam() -> PathBuf {
    fixtures_dir().join("checksum").join("chk1.bam")
}

fn normalize_checksum_file_line(text: &str) -> String {
    let mut lines = text.lines();
    let mut out = String::from("# Checksum 1.0 for file:\n");
    for line in lines.by_ref().skip(1) {
        out.push_str(line);
        out.push('\n');
    }
    out
}

fn checksum_all_row(text: &str) -> String {
    text.lines()
        .find(|line| line.starts_with("all        all"))
        .unwrap()
        .to_string()
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

fn without_pg_lines(text: &str) -> String {
    text.lines()
        .filter(|line| !line.starts_with("@PG\t"))
        .map(|line| {
            let mut line = line.to_string();
            line.push('\n');
            line
        })
        .collect()
}

fn without_pg_m5_lines(text: &str) -> String {
    text.lines()
        .filter(|line| !line.starts_with("@PG\t"))
        .map(|line| {
            if line.starts_with('@') {
                let fields = line
                    .split('\t')
                    .filter(|field| !field.starts_with("M5:"))
                    .collect::<Vec<_>>();
                return format!("{}\n", fields.join("\t"));
            }

            format!("{line}\n")
        })
        .collect()
}

fn without_pg_md_nm_lines(text: &str) -> String {
    text.lines()
        .filter(|line| !line.starts_with("@PG\t"))
        .map(|line| {
            if line.starts_with('@') {
                let fields = line
                    .split('\t')
                    .filter(|field| !field.starts_with("M5:"))
                    .collect::<Vec<_>>();
                return format!("{}\n", fields.join("\t"));
            }

            let fields = line
                .split('\t')
                .filter(|field| !field.starts_with("MD:Z:") && !field.starts_with("NM:i:"))
                .collect::<Vec<_>>();
            format!("{}\n", fields.join("\t"))
        })
        .collect()
}

fn write_unpadded_fasta(src: &Path, dst: &Path) {
    let text = std::fs::read_to_string(src).unwrap();
    std::fs::write(dst, text.replace('*', "")).unwrap();
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
fn depad_sam_matches_upstream_fixture() {
    let tmp = tmp_dir("depad-sam");
    let output = tmp.join("depad.sam");
    let input = fixtures_dir().join("dat").join("depad.001p.sam");
    let reference = fixtures_dir().join("dat").join("depad.001.fa");
    let expected = fixtures_dir().join("dat").join("depad.001u.sam");

    assert_eq!(
        exit_to_u8(depad::main(&argv(
            "depad",
            &[
                "-T",
                reference.to_str().unwrap(),
                "-s",
                "--no-PG",
                "-o",
                output.to_str().unwrap(),
                input.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        without_pg_m5_lines(&std::fs::read_to_string(output).unwrap()),
        without_pg_lines(&std::fs::read_to_string(expected).unwrap())
    );
}

#[test]
fn depad_bam_input_and_output_match_upstream_fixture() {
    let tmp = tmp_dir("depad-bam");
    let input_sam = fixtures_dir().join("dat").join("depad.001p.sam");
    let input_bam = tmp.join("padded.bam");
    let output_bam = tmp.join("depad.bam");
    let output_sam = tmp.join("depad-from-bam.sam");
    let reference = fixtures_dir().join("dat").join("depad.001.fa");
    let expected = without_pg_lines(
        &std::fs::read_to_string(fixtures_dir().join("dat").join("depad.001u.sam")).unwrap(),
    );

    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &input_sam,
        std::fs::File::create(&input_bam).unwrap(),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(depad::main(&argv(
            "depad",
            &[
                "-T",
                reference.to_str().unwrap(),
                "--no-PG",
                "-o",
                output_bam.to_str().unwrap(),
                input_sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    let bam_text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&output_bam, None)
            .unwrap();
    assert_eq!(without_pg_m5_lines(&bam_text), expected);

    assert_eq!(
        exit_to_u8(depad::main(&argv(
            "depad",
            &[
                "-T",
                reference.to_str().unwrap(),
                "-s",
                "--no-PG",
                "-o",
                output_sam.to_str().unwrap(),
                input_bam.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        without_pg_m5_lines(&std::fs::read_to_string(output_sam).unwrap()),
        expected
    );
}

#[test]
fn depad_cram_input_and_output_roundtrip() {
    let tmp = tmp_dir("depad-cram");
    let input_sam = tmp.join("padded.sam");
    let input_bam = tmp.join("padded.bam");
    let input_cram = tmp.join("padded.cram");
    let output_sam = tmp.join("depad-from-cram.sam");
    let output_cram = tmp.join("depad.out");
    let reference = tmp.join("padded.fa");
    let unpadded_reference = tmp.join("unpadded.fa");
    let expected = concat!(
        "@HD\tVN:1.6\n",
        "@SQ\tSN:ref\tLN:4\n",
        "r1\t0\tref\t1\t60\t2M1I2M\t*\t0\t0\tACAGT\tIIIII\n",
    );

    std::fs::write(&reference, ">ref\nAC*GT\n").unwrap();
    std::fs::write(
        &input_sam,
        "@HD\tVN:1.6\n@SQ\tSN:ref\tLN:5\nr1\t0\tref\t1\t60\t5M\t*\t0\t0\tACAGT\tIIIII\n",
    )
    .unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();
    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &input_sam,
        std::fs::File::create(&input_bam).unwrap(),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference(
        &input_bam,
        &reference,
        std::fs::File::create(&input_cram).unwrap(),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(depad::main(&argv(
            "depad",
            &[
                "-T",
                reference.to_str().unwrap(),
                "-s",
                "--no-PG",
                "-o",
                output_sam.to_str().unwrap(),
                input_cram.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        without_pg_md_nm_lines(&std::fs::read_to_string(output_sam).unwrap()),
        expected
    );

    assert_eq!(
        exit_to_u8(depad::main(&argv(
            "depad",
            &[
                "-T",
                reference.to_str().unwrap(),
                "--output-fmt=cram",
                "--no-PG",
                "-o",
                output_cram.to_str().unwrap(),
                input_sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    write_unpadded_fasta(&reference, &unpadded_reference);
    samtools_rs::reference::ensure_fai_index(&unpadded_reference, None).unwrap();
    let cram_text =
        htslib_rs::alignment_compat::view_cram_as_sam_text_from_path_with_reference_and_limit(
            &output_cram,
            &unpadded_reference,
            None,
        )
        .unwrap();
    assert_eq!(without_pg_md_nm_lines(&cram_text), expected);
}

#[test]
fn checksum_bam_matches_upstream_default_fixture() {
    let tmp = tmp_dir("checksum-default");
    let output = tmp.join("chk1.out");
    let input = sample_bam();

    assert_eq!(
        exit_to_u8(checksum::main(&argv(
            "checksum",
            &[input.to_str().unwrap(), "-o", output.to_str().unwrap()]
        ))),
        0
    );

    let actual = normalize_checksum_file_line(&std::fs::read_to_string(output).unwrap());
    let expected =
        std::fs::read_to_string(fixtures_dir().join("checksum").join("chk1.1.expected")).unwrap();
    assert_eq!(actual, expected);
}

#[test]
fn checksum_bam_qv_matches_upstream_fixture() {
    let tmp = tmp_dir("checksum-qv");
    let output = tmp.join("chk1-qv.out");
    let input = sample_bam();

    assert_eq!(
        exit_to_u8(checksum::main(&argv(
            "checksum",
            &[
                "-qv",
                input.to_str().unwrap(),
                "-o",
                output.to_str().unwrap(),
            ]
        ))),
        0
    );

    let actual = normalize_checksum_file_line(&std::fs::read_to_string(output).unwrap());
    let expected =
        std::fs::read_to_string(fixtures_dir().join("checksum").join("chk1.3.expected")).unwrap();
    assert_eq!(actual, expected);
}

#[test]
fn checksum_merge_matches_direct_checksum_for_default_reports() {
    let tmp = tmp_dir("checksum-merge");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
    let one = tmp.join("one.sam");
    let two = tmp.join("two.sam");
    let both = tmp.join("both.sam");
    std::fs::write(
        &one,
        format!("{header}r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGTT\t!!!!!\n"),
    )
    .unwrap();
    std::fs::write(
        &two,
        format!("{header}r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCAT\t#####\n"),
    )
    .unwrap();
    std::fs::write(
        &both,
        format!(
            "{header}r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGTT\t!!!!!\nr2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCAT\t#####\n"
        ),
    )
    .unwrap();

    let one_chk = tmp.join("one.chk");
    let two_chk = tmp.join("two.chk");
    let merged = tmp.join("merged.chk");
    let direct = tmp.join("direct.chk");
    for (input, output) in [(&one, &one_chk), (&two, &two_chk), (&both, &direct)] {
        assert_eq!(
            exit_to_u8(checksum::main(&argv(
                "checksum",
                &[input.to_str().unwrap(), "-o", output.to_str().unwrap()]
            ))),
            0
        );
    }
    assert_eq!(
        exit_to_u8(checksum::main(&argv(
            "checksum",
            &[
                "-m",
                one_chk.to_str().unwrap(),
                two_chk.to_str().unwrap(),
                "-o",
                merged.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        normalize_checksum_file_line(&std::fs::read_to_string(merged).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(direct).unwrap())
    );
}

#[test]
fn checksum_tabs_formats_report_as_tsv_and_merge_preserves_tsv_shape() {
    let tmp = tmp_dir("checksum-tabs");
    let input = sample_bam();
    let direct = tmp.join("direct.tsv");
    let merged = tmp.join("merged.tsv");

    assert_eq!(
        exit_to_u8(checksum::main(&argv(
            "checksum",
            &[
                "-T",
                input.to_str().unwrap(),
                "-o",
                direct.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(checksum::main(&argv(
            "checksum",
            &[
                "-m",
                "-T",
                direct.to_str().unwrap(),
                "-o",
                merged.to_str().unwrap(),
            ]
        ))),
        0
    );

    let direct_text = std::fs::read_to_string(&direct).unwrap();
    assert!(direct_text.contains("# Checksum 1.0 for file:\t"));
    assert!(direct_text.contains("# Group\tQC\tcount\tflag+seq\t+name\t+qual\t+aux\tcombined"));
    assert!(direct_text.lines().any(|line| {
        line.starts_with("all\tall\t144\t") && line.split('\t').count() == 8 && !line.contains("  ")
    }));

    assert_eq!(
        normalize_checksum_file_line(&std::fs::read_to_string(merged).unwrap()),
        normalize_checksum_file_line(&direct_text)
    );
}

#[test]
fn checksum_in_order_distinguishes_record_order() {
    let tmp = tmp_dir("checksum-order");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
    let rec_a = "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n";
    let rec_b = "b\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\n";
    let ab = tmp.join("ab.sam");
    let ba = tmp.join("ba.sam");
    std::fs::write(&ab, format!("{header}{rec_a}{rec_b}")).unwrap();
    std::fs::write(&ba, format!("{header}{rec_b}{rec_a}")).unwrap();

    let ab_default = tmp.join("ab.default.chk");
    let ba_default = tmp.join("ba.default.chk");
    let ab_ordered = tmp.join("ab.ordered.chk");
    let ba_ordered = tmp.join("ba.ordered.chk");
    for (input, output, extra) in [
        (&ab, &ab_default, &[][..]),
        (&ba, &ba_default, &[][..]),
        (&ab, &ab_ordered, &["-O"][..]),
        (&ba, &ba_ordered, &["-O"][..]),
    ] {
        let mut args = vec![input.to_str().unwrap(), "-o", output.to_str().unwrap()];
        args.extend_from_slice(extra);
        assert_eq!(exit_to_u8(checksum::main(&argv("checksum", &args))), 0);
    }

    assert_eq!(
        normalize_checksum_file_line(&std::fs::read_to_string(ab_default).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(ba_default).unwrap())
    );
    assert_ne!(
        normalize_checksum_file_line(&std::fs::read_to_string(ab_ordered).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(ba_ordered).unwrap())
    );
}

#[test]
fn checksum_check_pos_distinguishes_coordinates_and_merges() {
    let tmp = tmp_dir("checksum-check-pos");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
    let rec_at_1 = "same\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n";
    let rec_at_7 = "same\t0\tchr1\t7\t60\t4M\t*\t0\t0\tACGT\t!!!!\n";
    let one = tmp.join("one.sam");
    let two = tmp.join("two.sam");
    let both = tmp.join("both.sam");
    let shifted = tmp.join("shifted.sam");
    std::fs::write(&one, format!("{header}{rec_at_1}")).unwrap();
    std::fs::write(&two, format!("{header}{rec_at_7}")).unwrap();
    std::fs::write(&both, format!("{header}{rec_at_1}{rec_at_7}")).unwrap();
    std::fs::write(&shifted, format!("{header}{rec_at_1}{rec_at_1}")).unwrap();

    let both_default = tmp.join("both.default.chk");
    let shifted_default = tmp.join("shifted.default.chk");
    let both_pos = tmp.join("both.pos.chk");
    let shifted_pos = tmp.join("shifted.pos.chk");
    for (input, output, extra) in [
        (&both, &both_default, &[][..]),
        (&shifted, &shifted_default, &[][..]),
        (&both, &both_pos, &["-P"][..]),
        (&shifted, &shifted_pos, &["-P"][..]),
    ] {
        let mut args = vec![input.to_str().unwrap(), "-o", output.to_str().unwrap()];
        args.extend_from_slice(extra);
        assert_eq!(exit_to_u8(checksum::main(&argv("checksum", &args))), 0);
    }

    assert_eq!(
        normalize_checksum_file_line(&std::fs::read_to_string(&both_default).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(&shifted_default).unwrap())
    );
    assert_ne!(
        normalize_checksum_file_line(&std::fs::read_to_string(&both_pos).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(&shifted_pos).unwrap())
    );
    assert!(
        std::fs::read_to_string(&both_pos)
            .unwrap()
            .contains("+chr/pos  combined")
    );

    let one_pos = tmp.join("one.pos.chk");
    let two_pos = tmp.join("two.pos.chk");
    let merged_pos = tmp.join("merged.pos.chk");
    for (input, output) in [(&one, &one_pos), (&two, &two_pos)] {
        assert_eq!(
            exit_to_u8(checksum::main(&argv(
                "checksum",
                &[
                    "-P",
                    input.to_str().unwrap(),
                    "-o",
                    output.to_str().unwrap()
                ]
            ))),
            0
        );
    }
    assert_eq!(
        exit_to_u8(checksum::main(&argv(
            "checksum",
            &[
                "-m",
                one_pos.to_str().unwrap(),
                two_pos.to_str().unwrap(),
                "-o",
                merged_pos.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        normalize_checksum_file_line(&std::fs::read_to_string(merged_pos).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(both_pos).unwrap())
    );
}

#[test]
fn checksum_check_cigar_distinguishes_cigar_and_merges() {
    let tmp = tmp_dir("checksum-check-cigar");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
    let rec_4m = "same\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n";
    let rec_2m2m = "same\t0\tchr1\t1\t60\t2M2M\t*\t0\t0\tACGT\t!!!!\n";
    let one = tmp.join("one.sam");
    let two = tmp.join("two.sam");
    let both = tmp.join("both.sam");
    let merged_cigar_shape = tmp.join("merged-cigar-shape.sam");
    std::fs::write(&one, format!("{header}{rec_4m}")).unwrap();
    std::fs::write(&two, format!("{header}{rec_2m2m}")).unwrap();
    std::fs::write(&both, format!("{header}{rec_4m}{rec_2m2m}")).unwrap();
    std::fs::write(&merged_cigar_shape, format!("{header}{rec_4m}{rec_4m}")).unwrap();

    let both_default = tmp.join("both.default.chk");
    let merged_default = tmp.join("merged.default.chk");
    let both_cigar = tmp.join("both.cigar.chk");
    let merged_cigar = tmp.join("merged.cigar.chk");
    for (input, output, extra) in [
        (&both, &both_default, &[][..]),
        (&merged_cigar_shape, &merged_default, &[][..]),
        (&both, &both_cigar, &["-C"][..]),
        (&merged_cigar_shape, &merged_cigar, &["-C"][..]),
    ] {
        let mut args = vec![input.to_str().unwrap(), "-o", output.to_str().unwrap()];
        args.extend_from_slice(extra);
        assert_eq!(exit_to_u8(checksum::main(&argv("checksum", &args))), 0);
    }

    assert_eq!(
        normalize_checksum_file_line(&std::fs::read_to_string(&both_default).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(&merged_default).unwrap())
    );
    assert_ne!(
        normalize_checksum_file_line(&std::fs::read_to_string(&both_cigar).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(&merged_cigar).unwrap())
    );
    assert!(
        std::fs::read_to_string(&both_cigar)
            .unwrap()
            .contains("+cigar    combined")
    );

    let one_cigar = tmp.join("one.cigar.chk");
    let two_cigar = tmp.join("two.cigar.chk");
    let merged_report = tmp.join("merged-report.cigar.chk");
    for (input, output) in [(&one, &one_cigar), (&two, &two_cigar)] {
        assert_eq!(
            exit_to_u8(checksum::main(&argv(
                "checksum",
                &[
                    "-C",
                    input.to_str().unwrap(),
                    "-o",
                    output.to_str().unwrap(),
                ]
            ))),
            0
        );
    }
    assert_eq!(
        exit_to_u8(checksum::main(&argv(
            "checksum",
            &[
                "-m",
                one_cigar.to_str().unwrap(),
                two_cigar.to_str().unwrap(),
                "-o",
                merged_report.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        normalize_checksum_file_line(&std::fs::read_to_string(merged_report).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(both_cigar).unwrap())
    );
}

#[test]
fn checksum_check_mate_distinguishes_mate_position_and_merges() {
    let tmp = tmp_dir("checksum-check-mate");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
    let mate_at_5 = "same\t1\tchr1\t1\t60\t4M\t=\t5\t4\tACGT\t!!!!\n";
    let mate_at_9 = "same\t1\tchr1\t1\t60\t4M\t=\t9\t8\tACGT\t!!!!\n";
    let one = tmp.join("one.sam");
    let two = tmp.join("two.sam");
    let both = tmp.join("both.sam");
    let shifted = tmp.join("shifted.sam");
    std::fs::write(&one, format!("{header}{mate_at_5}")).unwrap();
    std::fs::write(&two, format!("{header}{mate_at_9}")).unwrap();
    std::fs::write(&both, format!("{header}{mate_at_5}{mate_at_9}")).unwrap();
    std::fs::write(&shifted, format!("{header}{mate_at_5}{mate_at_5}")).unwrap();

    let both_default = tmp.join("both.default.chk");
    let shifted_default = tmp.join("shifted.default.chk");
    let both_mate = tmp.join("both.mate.chk");
    let shifted_mate = tmp.join("shifted.mate.chk");
    for (input, output, extra) in [
        (&both, &both_default, &[][..]),
        (&shifted, &shifted_default, &[][..]),
        (&both, &both_mate, &["-M"][..]),
        (&shifted, &shifted_mate, &["-M"][..]),
    ] {
        let mut args = vec![input.to_str().unwrap(), "-o", output.to_str().unwrap()];
        args.extend_from_slice(extra);
        assert_eq!(exit_to_u8(checksum::main(&argv("checksum", &args))), 0);
    }

    assert_eq!(
        normalize_checksum_file_line(&std::fs::read_to_string(&both_default).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(&shifted_default).unwrap())
    );
    assert_ne!(
        normalize_checksum_file_line(&std::fs::read_to_string(&both_mate).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(&shifted_mate).unwrap())
    );
    assert!(
        std::fs::read_to_string(&both_mate)
            .unwrap()
            .contains("+mate     combined")
    );

    let one_mate = tmp.join("one.mate.chk");
    let two_mate = tmp.join("two.mate.chk");
    let merged_report = tmp.join("merged-report.mate.chk");
    for (input, output) in [(&one, &one_mate), (&two, &two_mate)] {
        assert_eq!(
            exit_to_u8(checksum::main(&argv(
                "checksum",
                &[
                    "-M",
                    input.to_str().unwrap(),
                    "-o",
                    output.to_str().unwrap(),
                ]
            ))),
            0
        );
    }
    assert_eq!(
        exit_to_u8(checksum::main(&argv(
            "checksum",
            &[
                "-m",
                one_mate.to_str().unwrap(),
                two_mate.to_str().unwrap(),
                "-o",
                merged_report.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        normalize_checksum_file_line(&std::fs::read_to_string(merged_report).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(both_mate).unwrap())
    );
}

#[test]
fn checksum_bamseqchksum_format_outputs_compat_rows() {
    let tmp = tmp_dir("checksum-bamseqchksum");
    let output = tmp.join("compat.chk");
    let input = sample_bam();

    assert_eq!(
        exit_to_u8(checksum::main(&argv(
            "checksum",
            &[
                "-B",
                input.to_str().unwrap(),
                "-o",
                output.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(output).unwrap();
    let lines = text.lines().collect::<Vec<_>>();
    assert_eq!(
        lines[0],
        "###\tset\tcount\t\tb_seq\tname_b_seq\tb_seq_qual\tb_seq_tags(BC,FI,QT,RT,TC)"
    );
    assert!(
        lines.iter().any(|line| {
            line.starts_with("all\tall\t144\t12eebb31\t7e82e134\t156fc353\t12eebb31")
        })
    );
    assert!(lines.iter().any(|line| {
        line.starts_with("ERR013140\tpass\t71\t2259e00a\t338c6a3\t7621d262\t2259e00a")
    }));
    assert!(!text.contains("\tfail\t"));
    assert!(!text.contains("combined"));
}

#[test]
fn checksum_all_expands_to_full_field_options_for_sam_and_bam() {
    let tmp = tmp_dir("checksum-all");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
    let sam = tmp.join("input.sam");
    let bam = tmp.join("input.bam");
    std::fs::write(
        &sam,
        format!(
            "{header}r1\t99\tchr1\t1\t60\t4M\t=\t10\t13\tACGT\t!!!!\tRG:Z:g1\tBC:Z:b\tNM:i:0\n"
        ),
    )
    .unwrap();
    write_bam_from_sam_text(&bam, &std::fs::read_to_string(&sam).unwrap());

    for input in [&sam, &bam] {
        let shorthand = tmp.join(format!(
            "{}.all.chk",
            input.file_stem().unwrap().to_string_lossy()
        ));
        let explicit = tmp.join(format!(
            "{}.explicit.chk",
            input.file_stem().unwrap().to_string_lossy()
        ));

        assert_eq!(
            exit_to_u8(checksum::main(&argv(
                "checksum",
                &[
                    "-a",
                    input.to_str().unwrap(),
                    "-o",
                    shorthand.to_str().unwrap()
                ]
            ))),
            0
        );
        assert_eq!(
            exit_to_u8(checksum::main(&argv(
                "checksum",
                &[
                    "-P",
                    "-C",
                    "-M",
                    "-O",
                    "-c",
                    "-b",
                    "0xfff",
                    "-f",
                    "0",
                    "-F",
                    "0",
                    "-t",
                    "*,cF,MD,NM",
                    input.to_str().unwrap(),
                    "-o",
                    explicit.to_str().unwrap(),
                ]
            ))),
            0
        );

        assert_eq!(
            normalize_checksum_file_line(&std::fs::read_to_string(shorthand).unwrap()),
            normalize_checksum_file_line(&std::fs::read_to_string(explicit).unwrap())
        );
    }
}

#[test]
fn checksum_sanitize_mutates_records_before_field_checks() {
    let tmp = tmp_dir("checksum-sanitize");
    let header = "@HD\tVN:1.6\n@SQ\tSN:x\tLN:5\n";
    let dirty = tmp.join("dirty.sam");
    let clean = tmp.join("clean.sam");
    let dirty_out = tmp.join("dirty.chk");
    let clean_out = tmp.join("clean.chk");
    std::fs::write(
        &dirty,
        format!("{header}r1\t0\tx\t6\t60\t4M\t*\t0\t0\tACGT\t!!!!\tMD:Z:4\tNM:i:0\tZZ:Z:z\n"),
    )
    .unwrap();
    std::fs::write(
        &clean,
        format!("{header}r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!\tZZ:Z:z\n"),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(checksum::main(&argv(
            "checksum",
            &[
                "-z",
                "all",
                "-P",
                "-C",
                "-t",
                "*,cF",
                dirty.to_str().unwrap(),
                "-o",
                dirty_out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(checksum::main(&argv(
            "checksum",
            &[
                "-P",
                "-C",
                "-t",
                "*,cF",
                clean.to_str().unwrap(),
                "-o",
                clean_out.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        normalize_checksum_file_line(&std::fs::read_to_string(dirty_out).unwrap()),
        normalize_checksum_file_line(&std::fs::read_to_string(clean_out).unwrap())
    );
}

#[test]
fn checksum_aux_wildcard_includes_sorted_tags() {
    let tmp = tmp_dir("checksum-aux-wildcard");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
    let input = tmp.join("input.sam");
    std::fs::write(
        &input,
        format!("{header}r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tZZ:Z:z\tBC:Z:b\tXA:B:C,1,2\tAA:i:1\n"),
    )
    .unwrap();

    let all_tags = tmp.join("all-tags.chk");
    let explicit_sorted = tmp.join("explicit-sorted.chk");
    let without_bc = tmp.join("without-bc.chk");
    let explicit_without_bc = tmp.join("explicit-without-bc.chk");
    for (tags, output) in [
        ("*", &all_tags),
        ("AA,BC,XA,ZZ", &explicit_sorted),
        ("*,BC", &without_bc),
        ("AA,XA,ZZ", &explicit_without_bc),
    ] {
        assert_eq!(
            exit_to_u8(checksum::main(&argv(
                "checksum",
                &[
                    "-t",
                    tags,
                    input.to_str().unwrap(),
                    "-o",
                    output.to_str().unwrap(),
                ]
            ))),
            0
        );
    }

    let all_row = checksum_all_row(&std::fs::read_to_string(all_tags).unwrap());
    let explicit_row = checksum_all_row(&std::fs::read_to_string(explicit_sorted).unwrap());
    let without_bc_row = checksum_all_row(&std::fs::read_to_string(without_bc).unwrap());
    let explicit_without_bc_row =
        checksum_all_row(&std::fs::read_to_string(explicit_without_bc).unwrap());

    assert_eq!(all_row, explicit_row);
    assert_eq!(without_bc_row, explicit_without_bc_row);
    assert_ne!(all_row, without_bc_row);
}

#[test]
fn reference_reconstructs_from_sam_md_tags() {
    let tmp = tmp_dir("reference-sam-md");
    let input = tmp.join("input.sam");
    let output = tmp.join("ref.fa");
    std::fs::write(
        &input,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:12\nr1\t0\tchr1\t1\t60\t4M1D4M\t*\t0\t0\tACGTTGCA\t!!!!!!!!\tMD:Z:4^A4\n",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(reference::main(&argv(
            "reference",
            &[
                "-q",
                input.to_str().unwrap(),
                "-o",
                output.to_str().unwrap()
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(output).unwrap(),
        ">chr1\nACGTATGCANNN\n"
    );
}

#[test]
fn reference_region_outputs_requested_slice_and_bam_input() {
    let tmp = tmp_dir("reference-region-bam");
    let sam = tmp.join("input.sam");
    let bam = tmp.join("input.bam");
    let sam_out = tmp.join("sam-region.fa");
    let bam_out = tmp.join("bam-region.fa");
    let text = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:12\nr1\t0\tchr1\t1\t60\t4M1D4M\t*\t0\t0\tACGTTGCA\t!!!!!!!!\tMD:Z:4^A4\n";
    std::fs::write(&sam, text).unwrap();
    write_bam_from_sam_text(&bam, text);

    for (input, output) in [(&sam, &sam_out), (&bam, &bam_out)] {
        assert_eq!(
            exit_to_u8(reference::main(&argv(
                "reference",
                &[
                    "-q",
                    "-r",
                    "chr1:3-7",
                    input.to_str().unwrap(),
                    "-o",
                    output.to_str().unwrap()
                ]
            ))),
            0
        );
    }

    let expected = ">chr1:3-7\nGTATG\n";
    assert_eq!(std::fs::read_to_string(sam_out).unwrap(), expected);
    assert_eq!(std::fs::read_to_string(bam_out).unwrap(), expected);
}

#[test]
fn reference_region_uses_indexed_bam_query_path() {
    let tmp = tmp_dir("reference-indexed-bam-region");
    let sam = tmp.join("input.sam");
    let bam = tmp.join("input.bam");
    let output = tmp.join("region.fa");
    let text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:12\nr1\t0\tchr1\t1\t60\t4M1D4M\t*\t0\t0\tACGTTGCA\t!!!!!!!!\tMD:Z:4^A4\n";
    std::fs::write(&sam, text).unwrap();
    write_bam_from_sam_text(&bam, text);
    assert_eq!(
        exit_to_u8(index::main(&argv("index", &[bam.to_str().unwrap()]))),
        0
    );

    assert_eq!(
        exit_to_u8(reference::main(&argv(
            "reference",
            &[
                "-q",
                "-r",
                "chr1:3-7",
                bam.to_str().unwrap(),
                "-o",
                output.to_str().unwrap()
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(output).unwrap(),
        ">chr1:3-7\nGTATG\n"
    );
}

/// `samtools reference` MD path on CRAM input. Builds the upstream
/// `test_reference` embed_ref CRAM in a temp dir, then checks the
/// external-reference MD2ref output byte-exactly against the committed
/// upstream `mpileup.MD.fa` fixture (whole-file + region).
#[test]
fn reference_cram_md_path_with_reference_matches_upstream() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let d = fixtures_dir();
    let refa = d.join("dat/mpileup.ref.fa");
    let tmp = tmp_dir("reference-cram-md");
    let cram = build_reference_embed_ref_cram(&tmp);
    let expected_full = d.join("reference/mpileup.MD.fa.expected");

    let cases: Vec<(&[&str], String)> = vec![
        (&[], std::fs::read_to_string(&expected_full).unwrap()),
        (
            &["-r", "17:1000-1500"],
            fasta_region(&expected_full, "17", 1000, 1500),
        ),
    ];
    for (i, (extra, expected)) in cases.into_iter().enumerate() {
        let out = tmp.join(format!("reference_md_{i}.fa"));
        let mut a: Vec<String> = vec![
            "samtools".into(),
            "--reference".into(),
            refa.to_str().unwrap().into(),
            "reference".into(),
            "-q".into(),
        ];
        a.extend(extra.iter().map(|s| s.to_string()));
        a.push(cram.to_str().unwrap().into());
        a.push("-o".into());
        a.push(out.to_str().unwrap().into());
        assert_eq!(
            exit_to_u8(samtools_run(
                a.iter().map(std::ffi::OsString::from).collect()
            )),
            0,
            "args={a:?}"
        );
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            expected,
            "reference MD case {i} must be byte-exact"
        );
    }
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
fn flagstat_cram_without_reference_succeeds() {
    // flagstat only inspects flags (reference-independent in CRAM), so
    // it now works without `--reference` via the synthesizing path,
    // matching `samtools flagstat foo.cram`.
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let cram = htslib_fixtures_dir().join("range.cram");

    assert_eq!(
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
fn idxstats_cram_without_reference_succeeds() {
    // idxstats needs only per-record reference id + flags
    // (reference-independent in CRAM), so it now works without
    // `--reference` via the synthesizing path. (Byte-exact counts vs
    // the BAM equivalent are proven by the htslib-rs unit test
    // `cram_summaries_without_reference_match_bam_flags_and_tids`.)
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let cram = htslib_fixtures_dir().join("range.cram");

    assert_eq!(
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
fn samples_custom_index_directory_reports_index_presence() {
    // Upstream `sam_index_load3` accepts a *directory* as the custom
    // index argument and finds `<dir>/<data-name>.bai` inside it. A bare
    // `.exists()` on the directory used to (wrongly) report nothing /
    // the directory itself; verify the resolver finds the relocated
    // index at this non-default location.
    let tmp = tmp_dir("samples-custom-index-dir");
    let bam = tmp.join("in.bam");
    let out = tmp.join("samples.txt");
    let idx_dir = tmp.join("indexes");
    std::fs::create_dir_all(&idx_dir).unwrap();
    std::fs::copy(sample_bam(), &bam).unwrap();

    assert_eq!(
        exit_to_u8(index::main(&argv("index", &[bam.to_str().unwrap()]))),
        0
    );
    std::fs::rename(tmp.join("in.bam.bai"), idx_dir.join("in.bam.bai")).unwrap();

    assert_eq!(
        exit_to_u8(samples::main(&argv(
            "samples",
            &[
                "-X",
                "-i",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
                idx_dir.to_str().unwrap(),
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
    let pg_line = header
        .lines()
        .find(|l| l.starts_with("@PG\t") && l.contains("\tCL:cat "))
        .expect("samtools cat @PG line present");
    assert!(pg_line.contains("\tPN:samtools"));
    assert!(pg_line.contains("\tVN:"));
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
fn cat_bam_inputs_write_bam_output_with_single_header() {
    let tmp = tmp_dir("cat-bam");
    let bam_a = tmp.join("a.bam");
    let bam_b = tmp.join("b.bam");
    let out = tmp.join("cat.bam");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
    write_bam_from_sam_text(
        &bam_a,
        &format!("{header}a1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n"),
    );
    write_bam_from_sam_text(
        &bam_b,
        &format!("{header}b1\t0\tchr1\t5\t60\t4M\t*\t0\t0\tTGCA\t####\n"),
    );

    assert_eq!(
        exit_to_u8(cat::main(&argv(
            "cat",
            &[
                "--no-PG",
                bam_a.to_str().unwrap(),
                bam_b.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    assert_eq!(text.matches("@SQ\tSN:chr1").count(), 1);
    assert!(text.lines().any(|line| line.starts_with("a1\t")));
    assert!(text.lines().any(|line| line.starts_with("b1\t")));
    assert!(!text.contains("\tCL:cat "));
}

#[test]
fn cat_reads_inputs_from_file_list_before_positionals() {
    let tmp = tmp_dir("cat-input-list");
    let bam_a = tmp.join("a.bam");
    let bam_b = tmp.join("b.bam");
    let bam_c = tmp.join("c.bam");
    let list = tmp.join("inputs.txt");
    let out = tmp.join("cat.bam");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
    write_bam_from_sam_text(
        &bam_a,
        &format!("{header}a1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n"),
    );
    write_bam_from_sam_text(
        &bam_b,
        &format!("{header}b1\t0\tchr1\t5\t60\t4M\t*\t0\t0\tTGCA\t####\n"),
    );
    write_bam_from_sam_text(
        &bam_c,
        &format!("{header}c1\t0\tchr1\t9\t60\t4M\t*\t0\t0\tAAAA\t!!!!\n"),
    );
    std::fs::write(
        &list,
        format!("{}\n\n{}\n", bam_a.display(), bam_b.display()),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(cat::main(&argv(
            "cat",
            &[
                "--no-PG",
                "-b",
                list.to_str().unwrap(),
                bam_c.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text =
        htslib_rs::alignment_compat::view_bam_as_sam_text_from_path_with_limit(&out, None).unwrap();
    let names: Vec<_> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(names, ["a1", "b1", "c1"]);
}

#[test]
fn cat_rejects_sam_input_like_upstream() {
    let sam = fixtures_dir().join("dat").join("view.001.sam");
    assert_eq!(
        exit_to_u8(cat::main(&argv("cat", &["--no-PG", sam.to_str().unwrap()]))),
        1
    );
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
fn reheader_rejects_sam_input_like_upstream() {
    let hdr = fixtures_dir().join("reheader").join("hdr.sam");
    let sam = fixtures_dir().join("dat").join("view.001.sam");
    assert_eq!(
        exit_to_u8(reheader::main(&argv(
            "reheader",
            &[hdr.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        1
    );
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
fn fastq_uses_original_quality_tag_when_requested() {
    let tmp = tmp_dir("fastq-original-quality");
    let sam = tmp.join("in.sam");
    let bam = tmp.join("in.bam");
    let sam_out = tmp.join("sam.fq");
    let bam_out = tmp.join("bam.fq");
    let text = concat!(
        "@HD\tVN:1.6\n",
        "@SQ\tSN:chr1\tLN:8\n",
        "plain\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tOQ:Z:abcd\n",
        "reverse\t16\tchr1\t1\t60\t4M\t*\t0\t0\tACGC\t####\tOQ:Z:1234\n",
        "fallback\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t$$$$\n",
    );
    std::fs::write(&sam, text).unwrap();
    write_bam_from_sam_text(&bam, text);

    for (input, output) in [(&sam, &sam_out), (&bam, &bam_out)] {
        assert_eq!(
            exit_to_u8(fastq::main(&argv(
                "fastq",
                &[
                    "-O",
                    "-o",
                    output.to_str().unwrap(),
                    input.to_str().unwrap(),
                ]
            ))),
            0
        );
    }

    let expected = concat!(
        "@plain\nACGT\n+\nabcd\n",
        "@reverse\nGCGT\n+\n4321\n",
        "@fallback\nTGCA\n+\n$$$$\n",
    );
    assert_eq!(std::fs::read_to_string(sam_out).unwrap(), expected);
    assert_eq!(std::fs::read_to_string(bam_out).unwrap(), expected);
}

#[test]
fn fastq_v_supplies_default_quality_for_missing_quality() {
    let tmp = tmp_dir("fastq-default-quality");
    let sam = tmp.join("in.sam");
    let out = tmp.join("reads.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "missing\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t*\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-v",
                "2",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@missing\nACGT\n+\n####\n"
    );
}

#[test]
fn fastq_rejects_invalid_default_quality() {
    let sam = fixtures_dir().join("dat").join("view.001.sam");
    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &["-v", "94", sam.to_str().unwrap()]
        ))),
        1
    );
}

#[test]
fn fastq_umi_appends_aux_tag_to_read_names() {
    let tmp = tmp_dir("fastq-umi");
    let sam = tmp.join("in.sam");
    let bam = tmp.join("in.bam");
    let sam_out = tmp.join("sam.fq");
    let bam_out = tmp.join("bam.fq");
    let custom_out = tmp.join("custom.fq");
    let text = concat!(
        "@HD\tVN:1.6\n",
        "@SQ\tSN:chr1\tLN:8\n",
        "rx\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRX:Z:ACG-TT\n",
        "ox\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\tOX:Z:TT\tRX:Z:AA\n",
        "hash#7\t65\tchr1\t1\t60\t4M\t=\t5\t8\tGGGG\t$$$$\tRX:Z:GG\n",
        "custom\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\t%%%%\tMI:Z:CC\tRX:Z:RR\n",
    );
    std::fs::write(&sam, text).unwrap();
    write_bam_from_sam_text(&bam, text);

    for (input, output) in [(&sam, &sam_out), (&bam, &bam_out)] {
        assert_eq!(
            exit_to_u8(fastq::main(&argv(
                "fastq",
                &[
                    "-U",
                    "-o",
                    output.to_str().unwrap(),
                    input.to_str().unwrap(),
                ]
            ))),
            0
        );
    }

    let expected = concat!(
        "@rx:ACG+TT\nACGT\n+\n!!!!\n",
        "@ox:TT\nTGCA\n+\n####\n",
        "@hash:GG#7/1\nGGGG\n+\n$$$$\n",
        "@custom:RR\nCCCC\n+\n%%%%\n",
    );
    assert_eq!(std::fs::read_to_string(sam_out).unwrap(), expected);
    assert_eq!(std::fs::read_to_string(bam_out).unwrap(), expected);

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-U",
                "--UMI-tag",
                "MI,RX",
                "-o",
                custom_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert!(
        std::fs::read_to_string(custom_out)
            .unwrap()
            .contains("@custom:CC\nCCCC\n+\n%%%%\n")
    );
}

#[test]
fn fastq_i_adds_casava_fields_from_barcode_tags() {
    let tmp = tmp_dir("fastq-casava");
    let sam = tmp.join("in.sam");
    let bam = tmp.join("in.bam");
    let sam_out = tmp.join("sam.fq");
    let bam_out = tmp.join("bam.fq");
    let custom_out = tmp.join("custom.fq");
    let fasta_out = tmp.join("reads.fa");
    let text = concat!(
        "@HD\tVN:1.6\n",
        "@SQ\tSN:chr1\tLN:8\n",
        "read1\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\tBC:Z:AAAA\n",
        "read2\t641\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\tBC:Z:CCCC\n",
        "missing\t0\tchr1\t1\t60\t4M\t*\t0\t0\tNNNN\t$$$$\n",
        "custom\t0\tchr1\t1\t60\t4M\t*\t0\t0\tGGGG\t%%%%\tXB:Z:GGGG\tBC:Z:TTTT\n",
    );
    std::fs::write(&sam, text).unwrap();
    write_bam_from_sam_text(&bam, text);

    for (input, output) in [(&sam, &sam_out), (&bam, &bam_out)] {
        assert_eq!(
            exit_to_u8(fastq::main(&argv(
                "fastq",
                &[
                    "-i",
                    "-n",
                    "-o",
                    output.to_str().unwrap(),
                    input.to_str().unwrap(),
                ]
            ))),
            0
        );
    }

    let expected = concat!(
        "@read1 1:N:0:AAAA\nACGT\n+\n!!!!\n",
        "@read2 2:Y:0:CCCC\nTGCA\n+\n####\n",
        "@missing 1:N:0:0\nNNNN\n+\n$$$$\n",
        "@custom 1:N:0:TTTT\nGGGG\n+\n%%%%\n",
    );
    assert_eq!(std::fs::read_to_string(sam_out).unwrap(), expected);
    assert_eq!(std::fs::read_to_string(bam_out).unwrap(), expected);

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-i",
                "--barcode-tag",
                "XB",
                "-n",
                "-o",
                custom_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert!(
        std::fs::read_to_string(custom_out)
            .unwrap()
            .contains("@custom 1:N:0:GGGG\nGGGG\n+\n%%%%\n")
    );

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fasta",
            &[
                "-i",
                "-n",
                "-o",
                fasta_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert!(
        std::fs::read_to_string(fasta_out)
            .unwrap()
            .starts_with(">read1 1:N:0:AAAA\nACGT\n")
    );
}

#[test]
fn fastq_index_files_extract_from_barcode_tag() {
    let tmp = tmp_dir("fastq-index-files");
    let sam = tmp.join("in.sam");
    let bam = tmp.join("in.bam");
    let single_out = tmp.join("0.fq");
    let i1_out = tmp.join("i1.fq");
    let i2_out = tmp.join("i2.fq");
    let text = concat!(
        "@HD\tVN:1.6\n",
        "foo\t4\t*\t0\t0\t*\t*\t0\t0\tACCCCCCCCCCCCCCCCCCCCT\txYYYYYYYYYYYYYYYYYYYYz\tBC:Z:AGGGGGGT-CGGGGGGT\tQT:Z:Xyyy1yyZ-Pqq1qqqR\n",
    );
    std::fs::write(&sam, text).unwrap();
    write_bam_from_sam_text(&bam, text);

    for input in [&sam, &bam] {
        std::fs::write(&single_out, "").unwrap();
        std::fs::write(&i1_out, "").unwrap();
        std::fs::write(&i2_out, "").unwrap();
        assert_eq!(
            exit_to_u8(fastq::main(&argv(
                "fastq",
                &[
                    "--index-format",
                    "i8n1i8",
                    "--i1",
                    i1_out.to_str().unwrap(),
                    "--i2",
                    i2_out.to_str().unwrap(),
                    "-0",
                    single_out.to_str().unwrap(),
                    input.to_str().unwrap(),
                ]
            ))),
            0
        );

        assert_eq!(
            std::fs::read_to_string(&single_out).unwrap(),
            "@foo\nACCCCCCCCCCCCCCCCCCCCT\n+\nxYYYYYYYYYYYYYYYYYYYYz\n"
        );
        assert_eq!(
            std::fs::read_to_string(&i1_out).unwrap(),
            "@foo\nAGGGGGGT\n+\nXyyy1yyZ\n"
        );
        assert_eq!(
            std::fs::read_to_string(&i2_out).unwrap(),
            "@foo\nCGGGGGGT\n+\nPqq1qqqR\n"
        );
    }
}

#[test]
fn fastq_index_files_accept_headerless_sam() {
    let tmp = tmp_dir("fastq-index-headerless-sam");
    let sam = tmp.join("in.sam");
    let single_out = tmp.join("0.fq");
    let i1_out = tmp.join("i1.fq");
    let i2_out = tmp.join("i2.fq");
    std::fs::write(
        &sam,
        "foo\t4\t*\t0\t0\t*\t*\t0\t0\tACCCCCCCCCCCCCCCCCCCCT\txYYYYYYYYYYYYYYYYYYYYz\tBC:Z:AGGGGGGT-CGGGGGGT\tQT:Z:Xyyy1yyZ-Pqq1qqqR\n",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "--index-format",
                "i8n1i8",
                "--i1",
                i1_out.to_str().unwrap(),
                "--i2",
                i2_out.to_str().unwrap(),
                "-0",
                single_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(&single_out).unwrap(),
        "@foo\nACCCCCCCCCCCCCCCCCCCCT\n+\nxYYYYYYYYYYYYYYYYYYYYz\n"
    );
    assert_eq!(
        std::fs::read_to_string(&i1_out).unwrap(),
        "@foo\nAGGGGGGT\n+\nXyyy1yyZ\n"
    );
    assert_eq!(
        std::fs::read_to_string(&i2_out).unwrap(),
        "@foo\nCGGGGGGT\n+\nPqq1qqqR\n"
    );
}

#[test]
fn fastq_index_emits_one_record_per_qname_group_with_casava_comment() {
    // A qname with two non-last-segment records (primary + a
    // supplementary that survives a relaxed -F) must yield exactly ONE
    // index record (upstream `flush_rec` is one-per-template). With -i
    // the index record gets the CASAVA comment, with the barcode
    // separator normalized to '+' and lower-cased bases upper-cased.
    let tmp = tmp_dir("fastq-index-group");
    let sam = tmp.join("in.sam");
    let i1_out = tmp.join("i1.fq");
    let main_out = tmp.join("0.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:16\n",
            // p1: primary R1 + supplementary R1 (same qname) — one index.
            "p1\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\tBC:Z:ac-gt\n",
            "p1\t2113\tchr1\t9\t60\t4M\t=\t5\t0\tACGT\t!!!!\tBC:Z:ac-gt\n",
            "p1\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\n",
            // p2: a separate template — its own single index record.
            "p2\t65\tchr1\t1\t60\t4M\t=\t5\t8\tGGGG\t!!!!\tBC:Z:TT+AA\n",
            "p2\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tCCCC\t####\n",
        ),
    )
    .unwrap();

    // -F 0 so the supplementary record is not filtered, proving the
    // dedup is by qname group (not by flag filtering).
    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-i",
                "-F",
                "0",
                "--index-format",
                "i*",
                "--i1",
                i1_out.to_str().unwrap(),
                "-0",
                main_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    // Exactly one index record per template, CASAVA comment present,
    // `ac-gt` → `AC+GT`, `TT+AA` stays `TT+AA`.
    // Index sequence is the raw BC segment (`ac`, `TT`); only the CASAVA
    // comment normalizes/upper-cases. Quality is the default (`"`) since
    // no QT tag is present.
    assert_eq!(
        std::fs::read_to_string(&i1_out).unwrap(),
        "@p1 1:N:0:AC+GT\nac\n+\n\"\"\n@p2 1:N:0:TT+AA\nTT\n+\n\"\"\n"
    );
}

#[test]
fn fastq_casava_barcode_propagates_from_r1_to_r2_mate() {
    // Only the R1 record carries BC; with -i, upstream copies the
    // barcode into the R2 mate's CASAVA comment within the qname group
    // (bam_fastq.c:952). So `*.2.fq` must show `2:N:0:AC+GT`, not
    // `2:N:0:0`.
    let tmp = tmp_dir("fastq-casava-propagate");
    let sam = tmp.join("in.sam");
    let r1 = tmp.join("1.fq");
    let r2 = tmp.join("2.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:16\n",
            "p\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\tBC:Z:ac-gt\n",
            "p\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-i",
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
        std::fs::read_to_string(&r1).unwrap(),
        "@p 1:N:0:AC+GT\nACGT\n+\n!!!!\n"
    );
    // R2 had no BC of its own — the group barcode is propagated in.
    assert_eq!(
        std::fs::read_to_string(&r2).unwrap(),
        "@p 2:N:0:AC+GT\nTGCA\n+\n####\n"
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
fn fastq_dash_t_and_dash_cap_t_combine_aux_tags() {
    let tmp = tmp_dir("fastq-t-T-combine");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:4\n",
            "read1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:rg1\tBC:Z:AAAA\tMD:Z:4\tNM:i:0\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-n",
                "-t",
                "-T",
                "MD",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("@read1"));
    assert!(text.contains("RG:Z:rg1"));
    assert!(text.contains("BC:Z:AAAA"));
    assert!(text.contains("MD:Z:4"));
    assert!(!text.contains("NM:i:0"));
}

#[test]
fn fastq_interleaves_read1_read2_when_paths_alias() {
    let tmp = tmp_dir("fastq-interleave");
    let sam = tmp.join("in.sam");
    let out = tmp.join("o.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:queryname\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "pair1\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\n",
            "pair1\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\n",
            "pair2\t65\tchr1\t1\t60\t4M\t=\t5\t8\tAAAA\t****\n",
            "pair2\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tCCCC\t&&&&\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-N",
                "-1",
                out.to_str().unwrap(),
                "-2",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        concat!(
            "@pair1/1\nACGT\n+\n!!!!\n",
            "@pair1/2\nTGCA\n+\n####\n",
            "@pair2/1\nAAAA\n+\n****\n",
            "@pair2/2\nCCCC\n+\n&&&&\n",
        )
    );
}

#[test]
fn fasta_reverse_strand_record_reverse_complemented_in_output() {
    let tmp = tmp_dir("fasta-revcomp");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.fa");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "fwd\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTAATT\t!!!!!!!!\n",
            "rev\t16\tchr1\t1\t60\t8M\t*\t0\t0\tACGTAATT\t!!!!!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fasta",
            &["-n", "-o", out.to_str().unwrap(), sam.to_str().unwrap(),]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    // Forward strand: as-stored.
    assert!(text.contains(">fwd\nACGTAATT\n"));
    // Reverse strand: reverse-complemented to AATTACGT.
    assert!(text.contains(">rev\nAATTACGT\n"));
}

#[test]
fn fastq_repeated_dash_d_unions_same_tag_values() {
    let tmp = tmp_dir("fastq-d-union");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:4\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tNM:i:13\n",
            "b\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\tNM:i:14\n",
            "c\t0\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\t****\tNM:i:0\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fastq::main(&argv(
            "fastq",
            &[
                "-n",
                "-d",
                "NM:13",
                "-d",
                "NM:14",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("@a\n"));
    assert!(text.contains("@b\n"));
    assert!(!text.contains("@c\n"));
}

#[test]
fn fastq_routes_r1_only_singletons_to_singleton_output() {
    let tmp = tmp_dir("fastq-r1-singleton");
    let sam = tmp.join("in.sam");
    let r1 = tmp.join("r1.fq");
    let r2 = tmp.join("r2.fq");
    let singleton = tmp.join("s.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:queryname\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "pair\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\n",
            "pair\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\n",
            "solo_r1\t73\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\t****\n",
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
        std::fs::read_to_string(&r1).unwrap(),
        "@pair\nACGT\n+\n!!!!\n"
    );
    assert_eq!(
        std::fs::read_to_string(&r2).unwrap(),
        "@pair\nTGCA\n+\n####\n"
    );
    assert_eq!(
        std::fs::read_to_string(&singleton).unwrap(),
        "@solo_r1\nAAAA\n+\n****\n"
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
fn fastq_single_default_path_reverse_complements_sam_and_bam_reads() {
    let tmp = tmp_dir("fastq-single-default-reverse");
    let sam = tmp.join("in.sam");
    let bam = tmp.join("in.bam");
    let sam_out = tmp.join("sam.fq");
    let bam_out = tmp.join("bam.fq");
    let text = concat!(
        "@HD\tVN:1.6\n",
        "@SQ\tSN:chr1\tLN:8\n",
        "forward\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!#$%\n",
        "reverse\t16\tchr1\t1\t60\t4M\t*\t0\t0\tACGC\t1234\n",
    );
    std::fs::write(&sam, text).unwrap();
    write_bam_from_sam_text(&bam, text);

    for (input, output) in [(&sam, &sam_out), (&bam, &bam_out)] {
        assert_eq!(
            exit_to_u8(fastq::main(&argv(
                "fastq",
                &["-o", output.to_str().unwrap(), input.to_str().unwrap(),]
            ))),
            0
        );
    }

    let expected = concat!("@forward\nACGT\n+\n!#$%\n", "@reverse\nGCGT\n+\n4321\n",);
    assert_eq!(std::fs::read_to_string(sam_out).unwrap(), expected);
    assert_eq!(std::fs::read_to_string(bam_out).unwrap(), expected);
}

#[test]
fn fastq_no_sc_matches_upstream_soft_clip_fixtures() {
    let fixtures = fixtures_dir();
    let input = fixtures.join("dat").join("bam2fq.sc.sam");
    let expected_dir = fixtures.join("bam2fq");
    let tmp = tmp_dir("fastq-no-sc-fixtures");
    let cases: [(&[&str], &str); 4] = [
        (
            &["-O", "--no-sc", "--no-sc-bkp", "-T", "s0"],
            "21.fq.expected",
        ),
        (&["-O", "--no-sc", "-Ts0"], "22.fq.expected"),
        (&["-O", "--no-sc"], "23.fq.expected"),
        (
            &["-O", "--no-sc", "--sc-aux", "s1", "-Ts0,s1"],
            "24.fq.expected",
        ),
    ];

    for (args, expected) in cases {
        let out = tmp.join(expected);
        let mut argv_parts = args.to_vec();
        argv_parts.extend(["-o", out.to_str().unwrap(), input.to_str().unwrap()]);
        assert_eq!(exit_to_u8(fastq::main(&argv("fastq", &argv_parts))), 0);
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            std::fs::read_to_string(expected_dir.join(expected)).unwrap(),
            "{expected}"
        );
    }
}

#[test]
fn fastq_cram_input_uses_top_level_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let fixtures = fixtures_dir();
    let tmp = tmp_dir("fastq-cram-ref-path");
    let sam = fixtures.join("dat").join("bam2fq.001.sam");
    let cram = tmp.join("in.cram");
    let sam_fq = tmp.join("sam.fq");
    let cram_fq = tmp.join("cram.fq");
    let sam_no_suffix = tmp.join("sam.n.fq");
    let cram_no_suffix = tmp.join("cram.n.fq");
    let reference = tmp.join("ref.fa");
    let ref_dir = fixtures.join("dat").join("cram_md5");
    std::fs::write(
        &reference,
        format!(
            ">ref1\n{}\n>ref2\n{}\n",
            std::fs::read_to_string(ref_dir.join("08c04d512d4797d9ba2a156c1daba468")).unwrap(),
            std::fs::read_to_string(ref_dir.join("7c35feac7036c1cdef3bee0cc4b21437")).unwrap()
        ),
    )
    .unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-C",
                "--no-PG",
                "-o",
                cram.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    for (extra, sam_out, cram_out) in [
        (&[][..], &sam_fq, &cram_fq),
        (&["-n"][..], &sam_no_suffix, &cram_no_suffix),
    ] {
        let mut sam_args = vec!["fastq"];
        sam_args.extend(extra.iter().copied());
        sam_args.extend(["-o", sam_out.to_str().unwrap(), sam.to_str().unwrap()]);
        assert_eq!(exit_to_u8(samtools_run(argv("samtools", &sam_args))), 0);

        let mut cram_args = vec!["--reference", reference.to_str().unwrap(), "fastq"];
        cram_args.extend(extra.iter().copied());
        cram_args.extend(["-o", cram_out.to_str().unwrap(), cram.to_str().unwrap()]);
        assert_eq!(exit_to_u8(samtools_run(argv("samtools", &cram_args))), 0);

        assert_eq!(
            std::fs::read_to_string(cram_out).unwrap(),
            std::fs::read_to_string(sam_out).unwrap()
        );
    }
}

#[test]
fn view_cram_roundtrip_preserves_bam2fq_setup_fields() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let fixtures = fixtures_dir();
    let tmp = tmp_dir("view-cram-bam2fq-setup");
    let sam = fixtures.join("dat").join("bam2fq.001.sam");
    let cram = tmp.join("bam2fq.001.cram");
    let out = tmp.join("roundtrip.sam");
    let reference = tmp.join("ref.fa");
    let ref_dir = fixtures.join("dat").join("cram_md5");
    std::fs::write(
        &reference,
        format!(
            ">ref1\n{}\n>ref2\n{}\n",
            std::fs::read_to_string(ref_dir.join("08c04d512d4797d9ba2a156c1daba468")).unwrap(),
            std::fs::read_to_string(ref_dir.join("7c35feac7036c1cdef3bee0cc4b21437")).unwrap()
        ),
    )
    .unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-C",
                "--no-PG",
                "-o",
                cram.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "view",
                "-h",
                "--no-PG",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    let fields = |qname: &str, flag: &str| -> Vec<String> {
        text.lines()
            .filter(|line| !line.starts_with('@'))
            .find_map(|line| {
                let fields: Vec<_> = line.split('\t').collect();

                if fields[0] == qname && fields[1] == flag {
                    Some(fields.into_iter().map(String::from).collect())
                } else {
                    None
                }
            })
            .unwrap_or_else(|| panic!("missing {qname}/{flag}"))
    };

    let primary_with_supplementary_mate = fields("ref1_grp2_p001", "99");
    assert_eq!(primary_with_supplementary_mate[7], "27");
    assert_eq!(primary_with_supplementary_mate[8], "34");

    let reversed_tlen_pair_r1 = fields("ref2_grp3_p002", "99");
    let reversed_tlen_pair_r2 = fields("ref2_grp3_p002", "147");
    assert_eq!(reversed_tlen_pair_r1[8], "-45");
    assert_eq!(reversed_tlen_pair_r2[8], "45");

    let cross_reference_pair = fields("ref12_grp1_p001", "97");
    assert_eq!(cross_reference_pair[8], "0");

    let unmapped = fields("unaligned_grp3_p001", "77");
    assert_eq!(unmapped[4], "0");
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
    assert_eq!(records[0][8], "4");
    assert!(records[0].contains(&"MC:Z:4M"));
    assert!(records[0].contains(&"MQ:i:60"));
    assert_eq!(records[1][1], "129");
    assert_eq!(records[1][6], "=");
    assert_eq!(records[1][7], "1");
    assert_eq!(records[1][8], "-4");
    assert!(records[1].contains(&"MC:Z:4M"));
    assert!(records[1].contains(&"MQ:i:60"));
}

#[test]
fn fixmate_recomputes_template_lengths_from_five_prime_positions() {
    let tmp = tmp_dir("fixmate-tlen");
    let sam = tmp.join("in.sam");
    let out = tmp.join("fixed.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:queryname\n",
            "@SQ\tSN:ref1\tLN:10000000100\n",
            "pair\t99\tref1\t10000000010\t30\t23M\t=\t10000000008\t2\tAAGTCGGCAGCGTCAGATGTGTA\t???????????????????????\n",
            "pair\t147\tref1\t10000000008\t30\t23M\t=\t10000000010\t-2\tCTGTCTCTTATACACATCTCCTT\t???????????????????????\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fixmate::main(&argv(
            "fixmate",
            &["-O", "sam", sam.to_str().unwrap(), out.to_str().unwrap(),]
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
    assert_eq!(records[0][8], "21");
    assert_eq!(records[1][8], "-21");
}

#[test]
fn fixmate_m_adds_mate_score_tags() {
    let tmp = tmp_dir("fixmate-ms");
    let sam = tmp.join("in.sam");
    let out = tmp.join("fixed.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:queryname\n",
            "@SQ\tSN:chr1\tLN:16\n",
            "pair\t65\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tms:i:1\n",
            "pair\t129\tchr1\t5\t60\t4M\t*\t0\t0\tTGCA\t5555\tms:i:2\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fixmate::main(&argv(
            "fixmate",
            &[
                "-m",
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
    assert!(records[0].contains(&"ms:i:80"));
    assert!(records[1].contains(&"ms:i:160"));
    assert!(!records[0].contains(&"ms:i:1"));
    assert!(!records[1].contains(&"ms:i:2"));
}

#[test]
fn fixmate_c_adds_template_cigar_tag_and_replaces_stale_tags() {
    let tmp = tmp_dir("fixmate-ct");
    let sam = tmp.join("in.sam");
    let out = tmp.join("fixed.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:queryname\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "pair\t65\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tct:Z:stale-a\n",
            "pair\t145\tchr1\t10\t60\t3M\t*\t0\t0\tTGA\t###\tct:Z:stale-b\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fixmate::main(&argv(
            "fixmate",
            &["-cO", "sam", sam.to_str().unwrap(), out.to_str().unwrap(),]
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
    assert!(records[0].contains(&"ct:Z:1F4M5T2R3M"));
    assert!(!records[0].contains(&"ct:Z:stale-a"));
    assert!(!records[1].contains(&"ct:Z:stale-b"));
    assert!(!records[1].iter().any(|field| field.starts_with("ct:Z:")));
}

#[test]
fn fixmate_rejects_coordinate_sorted_input() {
    let tmp = tmp_dir("fixmate-coordinate-sort");
    let sam = tmp.join("in.sam");
    let out = tmp.join("fixed.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "pair\t65\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "pair\t129\tchr1\t10\t60\t3M\t*\t0\t0\tTGA\t###\n",
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
        1
    );
    assert!(!out.exists());
}

#[test]
fn fixmate_default_sanitizer_matches_upstream_fixture() {
    let tmp = tmp_dir("fixmate-sanitize");
    let out = tmp.join("fixed.sam");
    let input = fixtures_dir().join("fixmate").join("sanitize.sam");
    let expected = fixtures_dir().join("fixmate").join("sanitize.sam.expected");

    assert_eq!(
        exit_to_u8(fixmate::main(&argv(
            "fixmate",
            &[
                "--no-PG",
                "-O",
                "sam",
                input.to_str().unwrap(),
                out.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        std::fs::read_to_string(expected).unwrap()
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
        ">chr1:3-10\nGTAC\nGTAC\n"
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
        ">chr1:3-10\nGTAC\nGTAC\n"
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
        ">chr1:3-10\nGTAC\nGTAC\n"
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
        ">chr1:1-4\nACGT\n>chr2:5-8\nCCCC\n"
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

    assert_eq!(std::fs::read_to_string(out).unwrap(), ">chr1:100-105\n");
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

    assert_eq!(std::fs::read_to_string(out).unwrap(), ">chr1:6-12\nCGT\n");
}

#[test]
fn faidx_continue_missing_region_outputs_empty_record() {
    let tmp = tmp_dir("fai-continue-missing");
    let fa = tmp.join("ref.fa");
    let out = tmp.join("out.fa");
    std::fs::write(&fa, ">chr1\nACGTACGT\n").unwrap();

    assert_eq!(
        exit_to_u8(faidx::main(&argv(
            "faidx",
            &[
                "--continue",
                "-o",
                out.to_str().unwrap(),
                fa.to_str().unwrap(),
                "chr1:1-4",
                "missing",
                "chr1:5-8",
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        ">chr1:1-4\nACGT\n>missing\n>chr1:5-8\nACGT\n"
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
        "@r1:2-7\nCGTA\nCG\n+\nbcde\nfg\n"
    );
}

#[test]
fn fqidx_out_of_range_region_exits_successfully_with_empty_record() {
    let tmp = tmp_dir("fqi-zero-region");
    let fq = tmp.join("reads.fq");
    let out = tmp.join("out.fq");
    std::fs::write(&fq, "@r1\nACGTACGT\n+\nabcdefgh\n").unwrap();

    assert_eq!(
        exit_to_u8(fqidx::main(&argv(
            "fqidx",
            &[
                "-o",
                out.to_str().unwrap(),
                fq.to_str().unwrap(),
                "r1:100-105",
            ]
        ))),
        0
    );

    assert_eq!(std::fs::read_to_string(out).unwrap(), "@r1:100-105\n+\n");
}

#[test]
fn fqidx_continue_missing_region_outputs_empty_record() {
    let tmp = tmp_dir("fqi-continue-missing");
    let fq = tmp.join("reads.fq");
    let out = tmp.join("out.fq");
    std::fs::write(&fq, "@r1\nACGTACGT\n+\nabcdefgh\n").unwrap();

    assert_eq!(
        exit_to_u8(fqidx::main(&argv(
            "fqidx",
            &[
                "--continue",
                "-o",
                out.to_str().unwrap(),
                fq.to_str().unwrap(),
                "r1:1-4",
                "missing",
                "r1:5-8",
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "@r1:1-4\nACGT\n+\nabcd\n@missing\n+\n@r1:5-8\nACGT\n+\nefgh\n"
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
        "@r1:2-7\nCGTA\nCG\n+\nbcde\nfg\n"
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
        "@r1:1-4\nACGT\n+\nabcd\n"
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
            &["-o", out.to_str().unwrap(), fq.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        htslib_rs::alignment_compat::count_bam_records_from_path(&out).unwrap(),
        2
    );
}

#[test]
fn import_fastq_to_cram() {
    let tmp = tmp_dir("imp-cram");
    let fq = tmp.join("in.fq");
    std::fs::write(&fq, "@r1\nACGT\n+\n!!!!\n@r2\nTTTT\n+\n####\n").unwrap();
    let out = tmp.join("out.cram");
    let explicit_out = tmp.join("out.byfmt");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &["-o", out.to_str().unwrap(), fq.to_str().unwrap()]
        ))),
        0
    );

    let records =
        htslib_rs::alignment_compat::summarize_cram_records_from_path_synthesizing_reference(&out)
            .unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].name_bytes(), Some(&b"r1"[..]));
    assert_eq!(records[0].flags_u16(), 4);
    assert_eq!(records[0].sequence_bytes(), b"ACGT");
    assert_eq!(records[1].name_bytes(), Some(&b"r2"[..]));
    assert_eq!(records[1].sequence_bytes(), b"TTTT");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "--output-fmt=cram",
                "-o",
                explicit_out.to_str().unwrap(),
                fq.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        htslib_rs::alignment_compat::summarize_cram_records_from_path_synthesizing_reference(
            &explicit_out
        )
        .unwrap()
        .len(),
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
fn import_paired_fastq_accepts_zero_singleton_input() {
    let tmp = tmp_dir("imp-paired-zero");
    let r1 = tmp.join("r1.fq");
    let r2 = tmp.join("r2.fq");
    let r0 = tmp.join("r0.fq");
    std::fs::write(&r1, "@p\nAC\n+\n!!\n").unwrap();
    std::fs::write(&r2, "@p\nTG\n+\n##\n").unwrap();
    std::fs::write(&r0, "@solo\nNN\n+\n$$\n").unwrap();
    let out = tmp.join("out.sam");

    assert_eq!(
        exit_to_u8(import::main(&argv(
            "import",
            &[
                "-1",
                r1.to_str().unwrap(),
                "-2",
                r2.to_str().unwrap(),
                "-0",
                r0.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "p\t77\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\n\
p\t141\t*\t0\t0\t*\t*\t0\t0\tTG\t##\n\
solo\t4\t*\t0\t0\t*\t*\t0\t0\tNN\t$$\n"
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
fn rmdup_accepts_reference_backed_cram_input_and_output() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    set_current_global_args(SamGlobalArgs::default());

    let tmp = tmp_dir("rmdup-cram");
    let reference = tmp.join("ref.fa");
    let sam = tmp.join("in.sam");
    let bam = tmp.join("in.bam");
    let cram = tmp.join("in.cram");
    let out = tmp.join("dedup.cram");

    std::fs::write(&reference, ">chr1\nACGTTGCA\n").unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "low\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "reverse\t16\tchr1\t1\t30\t4M\t*\t0\t0\tCCCC\t$$$$\n",
        ),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference(
        &bam,
        &reference,
        std::fs::File::create(&cram).unwrap(),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(rmdup::main(&argv(
            "rmdup",
            &[
                "-s",
                "--no-PG",
                "-O",
                "cram",
                "-T",
                reference.to_str().unwrap(),
                cram.to_str().unwrap(),
                out.to_str().unwrap()
            ]
        ))),
        0
    );
    set_current_global_args(SamGlobalArgs::default());

    let records = htslib_rs::alignment_compat::summarize_cram_records_from_path_with_reference(
        &out, &reference,
    )
    .unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].name_bytes(), Some(&b"high"[..]));
    assert_eq!(records[1].name_bytes(), Some(&b"reverse"[..]));
}

#[test]
fn rmdup_sam_input_removes_lower_scoring_paired_duplicate() {
    let tmp = tmp_dir("rmdup-pe-sam");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:200\n",
            "pair_a\t99\tchr1\t1\t60\t10M\t=\t91\t100\tAAAAAAAAAA\t!!!!!!!!!!\n",
            "pair_a\t147\tchr1\t91\t60\t10M\t=\t1\t-100\tTTTTTTTTTT\t!!!!!!!!!!\n",
            "pair_b\t99\tchr1\t1\t10\t10M\t=\t91\t100\tCCCCCCCCCC\t!!!!!!!!!!\n",
            "pair_b\t147\tchr1\t91\t10\t10M\t=\t1\t-100\tGGGGGGGGGG\t!!!!!!!!!!\n",
            "pair_c\t99\tchr1\t2\t10\t10M\t=\t91\t100\tACACACACAC\t!!!!!!!!!!\n",
            "pair_c\t147\tchr1\t91\t10\t10M\t=\t2\t-100\tTGTGTGTGTG\t!!!!!!!!!!\n",
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
    assert!(text.contains("\npair_a\t99\t"));
    assert!(text.contains("\npair_a\t147\t"));
    assert!(!text.contains("\npair_b\t"));
    assert!(text.contains("\npair_c\t99\t"));
    assert!(text.contains("\npair_c\t147\t"));
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
fn reset_no_pg_skips_adding_reset_program_entry() {
    let tmp = tmp_dir("reset-no-pg");
    let sam = tmp.join("in.sam");
    let no_pg_out = tmp.join("reset.no-pg.sam");
    let default_out = tmp.join("reset.default.sam");
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

    // --no-PG: existing @PG records are preserved, but no new samtools
    // @PG entry is added. This matches upstream `samtools reset` behavior
    // (`reset.c`'s `noPGentry` flag).
    assert_eq!(
        exit_to_u8(reset::main(&argv(
            "reset",
            &[
                "--no-PG",
                "-O",
                "sam",
                "-o",
                no_pg_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    let no_pg_text = std::fs::read_to_string(&no_pg_out).unwrap();
    assert!(no_pg_text.contains("@PG\tID:aligner"));
    assert!(no_pg_text.contains("@PG\tID:post"));
    assert!(!no_pg_text.contains("PN:samtools"));

    // Default: existing @PGs are preserved AND a new samtools @PG line is
    // chained onto the terminal program.
    assert_eq!(
        exit_to_u8(reset::main(&argv(
            "reset",
            &[
                "-O",
                "sam",
                "-o",
                default_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    let default_text = std::fs::read_to_string(default_out).unwrap();
    assert!(default_text.contains("@PG\tID:aligner"));
    assert!(default_text.contains("@PG\tID:post"));
    assert!(default_text.contains("PN:samtools"));
    assert!(default_text.contains("PP:post"));
}

#[test]
fn reset_reject_pg_removes_named_pg_and_all_subsequent() {
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

    // Upstream `reset.c:223`: `--reject-PG ID` keeps @PG lines until the
    // first matching ID, then drops it and *every subsequent @PG*
    // ("from this PG onwards"). bwa_index is the first @PG, so all
    // three input @PG are dropped; only the added samtools @PG remains.
    let text = std::fs::read_to_string(out).unwrap();
    assert!(!text.contains("@PG\tID:bwa_index"));
    assert!(!text.contains("@PG\tID:bwa_aln"));
    assert!(!text.contains("@PG\tID:qc"));
    assert!(text.contains("@PG\tID:samtools"));
}

#[test]
fn reset_matches_upstream_test_reset_fixtures() {
    use samtools_rs::commands::reset;
    // Byte-exact vs upstream test_reset (harness `hskip=1` drops the
    // first output line, `ignore_pg_header` strips @PG): -o SAM from a
    // SAM input, -o SAM from BAM/CRAM input, and --no-RG; plus the
    // reject.1 / reject.2 @PG-count assertions.
    let d = fixtures_dir();
    let tmp = tmp_dir("reset-fixtures");
    // hskip=1 + strip @PG
    let norm = |s: &str| -> String {
        s.lines()
            .skip(1)
            .filter(|l| !l.starts_with("@PG\t"))
            .map(|l| format!("{l}\n"))
            .collect()
    };
    let exp_norm = |s: &str| -> String {
        s.lines()
            .filter(|l| !l.starts_with("@PG\t"))
            .map(|l| format!("{l}\n"))
            .collect()
    };
    let cases: &[(&[&str], &str, &str)] = &[
        (
            &["--dupflag", "-o", "@OUT@", "dat/mpileup.1.sam"],
            "@OUT@",
            "reset/basic.output.mp.1.expected",
        ),
        (
            &["--dupflag", "-o", "@OUT@", "dat/test_input_1_a.bam"],
            "@OUT@",
            "reset/basic.bam.input.expected",
        ),
        (
            &["--dupflag", "-o", "@OUT@", "dat/test_input_1_a.cram"],
            "@OUT@",
            "reset/basic.cram.input.expected",
        ),
        (
            &[
                "--dupflag",
                "--reject-PG",
                "bwa_index",
                "dat/mpileup.1.sam",
                "--no-RG",
                "-o",
                "@OUT@",
            ],
            "@OUT@",
            "reset/output.nRG.1.expected",
        ),
        (
            &[
                "--dupflag",
                "--reject-PG",
                "bwa_index",
                "dat/mpileup.1.sam",
                "--no-RG",
                "--keep-tag",
                "RG",
                "-o",
                "@OUT@",
            ],
            "@OUT@",
            "reset/output.nRG.2.expected",
        ),
        (
            &[
                "--dupflag",
                "--reject-PG",
                "bwa_index",
                "dat/mpileup.1.sam",
                "--no-RG",
                "--keep-tag",
                "X0,MD",
                "-o",
                "@OUT@",
            ],
            "@OUT@",
            "reset/output.keep.1.expected",
        ),
        (
            &[
                "--dupflag",
                "--reject-PG",
                "bwa_index",
                "dat/mpileup.1.sam",
                "--no-RG",
                "--remove-tag",
                "X0,X1,MD",
                "--keep-tag",
                "X0,MD",
                "-o",
                "@OUT@",
            ],
            "@OUT@",
            "reset/output.keep.1.expected",
        ),
        (
            &[
                "--dupflag",
                "--reject-PG",
                "bwa_index",
                "dat/mpileup.1.sam",
                "--no-RG",
                "--remove-tag",
                "X0,X1,MD",
                "-o",
                "@OUT@",
            ],
            "@OUT@",
            "reset/output.keep.2.expected",
        ),
        (
            &[
                "--dupflag",
                "--reject-PG",
                "bwa_index",
                "dat/mpileup.1.sam",
                "--no-RG",
                "-x",
                "X0,X1,MD",
                "-o",
                "@OUT@",
            ],
            "@OUT@",
            "reset/output.keep.2.expected",
        ),
        (
            &[
                "--dupflag",
                "--reject-PG",
                "bwa_index",
                "dat/mpileup.1.sam",
                "--no-RG",
                "--remove-tag",
                "^X0,MD",
                "--keep-tag",
                "X1",
                "-o",
                "@OUT@",
            ],
            "@OUT@",
            "reset/output.keep.3.expected",
        ),
    ];
    for (i, (args, _, expected)) in cases.iter().enumerate() {
        let out = tmp.join(format!("r{i}.sam"));
        let v: Vec<OsString> = std::iter::once(OsString::from("reset"))
            .chain(args.iter().map(|a| {
                OsString::from(if *a == "@OUT@" {
                    out.to_str().unwrap().to_string()
                } else if a.starts_with("dat/") || a.starts_with("reset/") {
                    d.join(a).to_str().unwrap().to_string()
                } else {
                    a.to_string()
                })
            }))
            .collect();
        assert_eq!(exit_to_u8(reset::main(&v)), 0, "{expected}");
        assert_eq!(
            norm(&std::fs::read_to_string(&out).unwrap()),
            exp_norm(&std::fs::read_to_string(d.join(expected)).unwrap()),
            "reset {expected}"
        );
    }

    // reject.1: count of the added samtools @PG line.
    let o = tmp.join("rej1.sam");
    assert_eq!(
        exit_to_u8(reset::main(&argv(
            "reset",
            &[
                "--dupflag",
                "--reject-PG",
                "bwa_index",
                d.join("dat/mpileup.1.sam").to_str().unwrap(),
                "-o",
                o.to_str().unwrap(),
            ]
        ))),
        0
    );
    let txt = std::fs::read_to_string(&o).unwrap();
    let n = txt
        .lines()
        .filter(|l| l.starts_with("@PG\tID:samtools\tPN:samtools"))
        .count();
    assert_eq!(
        n.to_string(),
        std::fs::read_to_string(d.join("reset/reject.1.expected"))
            .unwrap()
            .trim()
    );

    // reject.2: total @PG count after positional "onwards" removal.
    let o = tmp.join("rej2.sam");
    assert_eq!(
        exit_to_u8(reset::main(&argv(
            "reset",
            &[
                "--dupflag",
                "--reject-PG",
                "sam_to_fixed_bam",
                d.join("dat/mpileup.1.sam").to_str().unwrap(),
                "-o",
                o.to_str().unwrap(),
            ]
        ))),
        0
    );
    let txt = std::fs::read_to_string(&o).unwrap();
    let n = txt.lines().filter(|l| l.starts_with("@PG\tID:")).count();
    assert_eq!(
        n.to_string(),
        std::fs::read_to_string(d.join("reset/reject.2.expected"))
            .unwrap()
            .trim()
    );
}

#[test]
fn reset_writes_cram_output() {
    let tmp = tmp_dir("reset-cram-output");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.cram");
    let explicit_out = tmp.join("out.byfmt");

    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t16\tchr1\t2\t60\t5M\t*\t0\t0\tACGTT\t!!!!!\tNM:i:1\tMD:Z:3A\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(reset::main(&argv(
            "reset",
            &[
                "--no-PG",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    let records =
        htslib_rs::alignment_compat::summarize_cram_records_from_path_synthesizing_reference(&out)
            .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].name_bytes(), Some(&b"r1"[..]));
    assert_eq!(records[0].flags_u16() & 4, 4);
    assert_eq!(records[0].sequence_bytes(), b"AACGT");

    assert_eq!(
        exit_to_u8(reset::main(&argv(
            "reset",
            &[
                "--no-PG",
                "--output-fmt=cram",
                "-o",
                explicit_out.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        htslib_rs::alignment_compat::summarize_cram_records_from_path_synthesizing_reference(
            &explicit_out
        )
        .unwrap()
        .len(),
        1
    );
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
fn split_cram_input_by_rg_to_sam_outputs() {
    let tmp = tmp_dir("split-cram");
    let cram = fixtures_dir().join("checksum").join("chk2.cram");
    let tmpl = tmp.join("out.%!.%.");
    let unk = tmp.join("unknown.sam");

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
                cram.to_str().unwrap(),
            ]
        ))),
        0
    );

    let err013140 = std::fs::read_to_string(tmp.join("out.ERR013140.sam")).unwrap();
    let err156632 = std::fs::read_to_string(tmp.join("out.ERR156632.sam")).unwrap();
    let unknown = std::fs::read_to_string(unk).unwrap();
    assert!(err013140.lines().any(|line| line.starts_with("ERR013140.")));
    assert!(err156632.lines().any(|line| line.starts_with("ERR156632.")));
    assert!(unknown.lines().any(|line| line.starts_with("ERR016352.")));
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
fn split_infers_sam_output_from_template_extension() {
    let tmp = tmp_dir("split-template-sam");
    let sam = tmp.join("in.sam");
    let tmpl = tmp.join("out.%#.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
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
                "--no-PG",
                "-f",
                tmpl.to_str().unwrap(),
                sam.to_str().unwrap()
            ]
        ))),
        0
    );

    let g1 = std::fs::read_to_string(tmp.join("out.0.sam")).unwrap();
    let g2 = std::fs::read_to_string(tmp.join("out.1.sam")).unwrap();
    assert!(g1.starts_with("@HD\t"));
    assert!(g1.lines().any(|line| line.starts_with("r1\t")));
    assert!(g2.lines().any(|line| line.starts_with("r2\t")));
}

#[test]
fn split_errors_on_missing_rg_without_unaccounted_output_like_upstream() {
    let tmp = tmp_dir("split-missing-rg-error");
    let sam = tmp.join("in.sam");
    let tmpl = tmp.join("out.%#.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:g1\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:g1\n",
            "r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\n",
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
                sam.to_str().unwrap()
            ]
        ))),
        1
    );

    let g1 = std::fs::read_to_string(tmp.join("out.0.sam")).unwrap();
    assert!(g1.lines().any(|line| line.starts_with("r1\t")));
    assert!(!g1.lines().any(|line| line.starts_with("r2\t")));
}

#[test]
fn split_sam_input_by_rg_to_cram_outputs_with_reference() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("split-cram-output");
    let sam = tmp.join("in.sam");
    let reference = tmp.join("ref.fa");
    let tmpl = tmp.join("out.%!.%.");
    let unk = tmp.join("unknown.cram");
    std::fs::write(&reference, ">chr1\nACGTACGT\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:g1\tSM:s1\n",
            "@RG\tID:g2\tSM:s2\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:g1\n",
            "r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tCGTA\t####\tRG:Z:g2\n",
            "r3\t0\tchr1\t3\t60\t4M\t*\t0\t0\tGTAC\t$$$$\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "split",
                "--output-fmt",
                "cram",
                "--no-PG",
                "-f",
                tmpl.to_str().unwrap(),
                "-u",
                unk.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let g1 = tmp.join("out.g1.cram");
    let g2 = tmp.join("out.g2.cram");
    assert!(g1.exists());
    assert!(g2.exists());
    assert!(unk.exists());
    assert_eq!(
        htslib_rs::alignment_compat::query_cram_records_all_from_path_with_reference(
            &g1, &reference
        )
        .unwrap()
        .len(),
        1
    );
    assert_eq!(
        htslib_rs::alignment_compat::query_cram_records_all_from_path_with_reference(
            &g2, &reference
        )
        .unwrap()
        .len(),
        1
    );
    assert_eq!(
        htslib_rs::alignment_compat::query_cram_records_all_from_path_with_reference(
            &unk, &reference
        )
        .unwrap()
        .len(),
        1
    );
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

#[test]
fn fixmate_adds_pg_line_by_default() {
    let tmp = tmp_dir("fixmate-pg");
    let sam = tmp.join("in.sam");
    let out = tmp.join("fixed.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:queryname\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "r1\t99\tchr1\t1\t60\t10M\t=\t91\t100\tACGTACGTAC\t!!!!!!!!!!\n",
            "r1\t147\tchr1\t91\t60\t10M\t=\t1\t-100\tACGTACGTAC\t!!!!!!!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fixmate::main(&argv(
            "fixmate",
            &["-O", "sam", sam.to_str().unwrap(), out.to_str().unwrap(),]
        ))),
        0
    );
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("PN:samtools"));
}

#[test]
fn fixmate_no_pg_omits_pg_line() {
    let tmp = tmp_dir("fixmate-no-pg");
    let sam = tmp.join("in.sam");
    let out = tmp.join("fixed.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:queryname\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "r1\t99\tchr1\t1\t60\t10M\t=\t91\t100\tACGTACGTAC\t!!!!!!!!!!\n",
            "r1\t147\tchr1\t91\t60\t10M\t=\t1\t-100\tACGTACGTAC\t!!!!!!!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fixmate::main(&argv(
            "fixmate",
            &[
                "--no-PG",
                "-O",
                "sam",
                sam.to_str().unwrap(),
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(!text.contains("PN:samtools"));
}

#[test]
fn rmdup_adds_pg_line_by_default() {
    let tmp = tmp_dir("rmdup-pg");
    let sam = tmp.join("in.sam");
    let out = tmp.join("dedup.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
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
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("PN:samtools"));
}

#[test]
fn rmdup_no_pg_omits_pg_line() {
    let tmp = tmp_dir("rmdup-no-pg");
    let sam = tmp.join("in.sam");
    let out = tmp.join("dedup.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(rmdup::main(&argv(
            "rmdup",
            &["--no-PG", sam.to_str().unwrap(), out.to_str().unwrap()]
        ))),
        0
    );
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(!text.contains("PN:samtools"));
}

#[test]
fn addreplacerg_adds_pg_line_by_default_and_omits_with_no_pg() {
    let tmp = tmp_dir("addreplacerg-pg");
    let sam = tmp.join("in.sam");
    let default_out = tmp.join("default.sam");
    let no_pg_out = tmp.join("no_pg.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(addreplacerg::main(&argv(
            "addreplacerg",
            &[
                "-r",
                "ID:g1",
                "-r",
                "SM:s1",
                "-O",
                "sam",
                "-o",
                default_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    let default_text = std::fs::read_to_string(&default_out).unwrap();
    assert!(default_text.contains("PN:samtools"));

    assert_eq!(
        exit_to_u8(addreplacerg::main(&argv(
            "addreplacerg",
            &[
                "--no-PG",
                "-r",
                "ID:g1",
                "-r",
                "SM:s1",
                "-O",
                "sam",
                "-o",
                no_pg_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    let no_pg_text = std::fs::read_to_string(&no_pg_out).unwrap();
    assert!(!no_pg_text.contains("PN:samtools"));
}

#[test]
fn addreplacerg_writes_bam_output_with_rg_header_and_tag() {
    let tmp = tmp_dir("addreplacerg-bam-output");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.bam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(addreplacerg::main(&argv(
            "addreplacerg",
            &[
                "--no-PG",
                "-r",
                "ID:g1",
                "-r",
                "SM:s1",
                "-O",
                "bam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let mut reader = bam::io::Reader::new(std::fs::File::open(&out).unwrap());
    let header = reader.read_header().unwrap();
    assert!(header.read_groups().contains_key("g1".as_bytes()));

    let mut record = sam::alignment::RecordBuf::default();
    assert_ne!(reader.read_record_buf(&header, &mut record).unwrap(), 0);
    let tag = sam::alignment::record::data::field::Tag::from([b'R', b'G']);
    let value = record.data().get(&tag).unwrap();
    assert_eq!(
        value,
        &sam::alignment::record_buf::data::field::Value::String("g1".into())
    );
}

#[test]
fn addreplacerg_writes_cram_output_with_reference() {
    let tmp = tmp_dir("addreplacerg-cram-output");
    let sam = tmp.join("in.sam");
    let reference = tmp.join("ref.fa");
    let out = tmp.join("out.cram");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();
    std::fs::write(&reference, ">chr1\nACGTACGT\n").unwrap();

    assert_eq!(
        exit_to_u8(addreplacerg::main(&argv(
            "addreplacerg",
            &[
                "--no-PG",
                "-r",
                "ID:g1",
                "-r",
                "SM:s1",
                "--output-fmt",
                "cram",
                "-T",
                reference.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    // Read the CRAM back through the shared reference-backed decoder and
    // confirm the new @RG header and per-record RG:Z: tag are present.
    let text = htslib_rs::alignment_compat::view_cram_as_sam_text_with_reference(
        std::io::Cursor::new(std::fs::read(&out).unwrap()),
        &reference,
        None,
    )
    .unwrap();
    assert!(text.contains("@RG\t") && text.contains("ID:g1"));
    assert!(text.contains("RG:Z:g1"));
}

#[test]
fn addreplacerg_accepts_reference_backed_cram_input() {
    let tmp = tmp_dir("addreplacerg-cram-input");
    let sam = tmp.join("in.sam");
    let bam = tmp.join("in.bam");
    let cram = tmp.join("in.cram");
    let reference = tmp.join("ref.fa");
    let out = tmp.join("out.sam");

    std::fs::write(&reference, ">chr1\nACGTACGT\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();
    htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference(
        &bam,
        &reference,
        std::fs::File::create(&cram).unwrap(),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(addreplacerg::main(&argv(
            "addreplacerg",
            &[
                "--no-PG",
                "-r",
                "ID:g1",
                "-T",
                reference.to_str().unwrap(),
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                cram.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.contains("@RG\tID:g1"));
    assert!(text.contains("\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:g1\n"));
}

#[test]
fn addreplacerg_bam_input_to_sam_honors_orphan_only_mode() {
    let tmp = tmp_dir("addreplacerg-bam-input");
    let sam = tmp.join("in.sam");
    let bam = tmp.join("with_rg.bam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:old\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(addreplacerg::main(&argv(
            "addreplacerg",
            &[
                "--no-PG",
                "-r",
                "ID:old",
                "-O",
                "bam",
                "-o",
                bam.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    // orphan_only must not overwrite records that already carry an RG:Z
    // tag. We add a *new* @RG via -r (so the header entry exists) and
    // confirm the existing RG:Z:old tags survive untouched.
    assert_eq!(
        exit_to_u8(addreplacerg::main(&argv(
            "addreplacerg",
            &[
                "--no-PG",
                "-r",
                "@RG\tID:new",
                "-m",
                "orphan_only",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.contains("RG:Z:old"));
    assert!(!text.contains("RG:Z:new"));
}

#[test]
fn addreplacerg_dash_cap_r_unknown_id_is_rejected() {
    let tmp = tmp_dir("addreplacerg-unknown-rg");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.4\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:present\tCN:SC\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    // -R with an ID not present in the header must fail (upstream parity).
    assert_ne!(
        exit_to_u8(addreplacerg::main(&argv(
            "addreplacerg",
            &[
                "--no-PG",
                "-R",
                "absent",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
}

#[test]
fn addreplacerg_defaults_to_first_header_rg_and_preserves_lines() {
    let tmp = tmp_dir("addreplacerg-default-rg");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.4\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:first\tCN:SC\n",
            "@RG\tID:second\tCN:SC\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\tRG:Z:second\n",
        ),
    )
    .unwrap();

    // No -r / -R: default to the first @RG ID, keep both @RG header lines,
    // overwrite all record RG tags with the first ID.
    assert_eq!(
        exit_to_u8(addreplacerg::main(&argv(
            "addreplacerg",
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
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("@RG\tID:first\tCN:SC"));
    assert!(text.contains("@RG\tID:second\tCN:SC"));
    assert!(text.contains("r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:first"));
    assert!(text.contains("r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tTGCA\t####\tRG:Z:first"));
}

#[test]
fn addreplacerg_r_overwrite_all_removes_other_header_rg_lines() {
    let tmp = tmp_dir("addreplacerg-overwrite-rg");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.4\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "@RG\tID:old\tCN:SC\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(addreplacerg::main(&argv(
            "addreplacerg",
            &[
                "--no-PG",
                "-O",
                "sam",
                "-r",
                "@RG\tID:new\tCN:SC",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("@RG\tID:new\tCN:SC"));
    assert!(!text.contains("@RG\tID:old"));
    assert!(text.contains("r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:new"));
}

#[test]
fn fixmate_r_removes_unmapped_and_secondary_records() {
    let tmp = tmp_dir("fixmate-r");
    let sam = tmp.join("in.sam");
    let out = tmp.join("fixed.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:queryname\n",
            "@SQ\tSN:chr1\tLN:100\n",
            // Pair where second mate is unmapped (FUNMAP=0x4 → flag 77 for read1 / 141 for read2)
            "r1\t77\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\t!!!!!!!!!!\n",
            "r1\t141\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\t!!!!!!!!!!\n",
            // Secondary alignment (FSECONDARY=0x100)
            "r2\t256\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t!!!!!!!!!!\n",
            // Primary mapped pair
            "r3\t99\tchr1\t1\t60\t10M\t=\t91\t100\tACGTACGTAC\t!!!!!!!!!!\n",
            "r3\t147\tchr1\t91\t60\t10M\t=\t1\t-100\tACGTACGTAC\t!!!!!!!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(fixmate::main(&argv(
            "fixmate",
            &[
                "-r",
                "-O",
                "sam",
                sam.to_str().unwrap(),
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(!text.contains("\nr1\t"), "unmapped pair must be removed");
    assert!(
        !text.contains("\nr2\t"),
        "secondary alignment must be removed"
    );
    assert!(text.contains("\nr3\t"), "primary pair must be retained");
}

#[test]
fn cat_r_region_restricts_to_indexed_bam_records() {
    use samtools_rs::commands::view;

    let tmp = tmp_dir("cat-region");
    let bam = sample_bam();
    let indexed = tmp.join("indexed.bam");
    std::fs::copy(&bam, &indexed).unwrap();
    assert_eq!(
        exit_to_u8(index::main(&argv("index", &[indexed.to_str().unwrap()]))),
        0
    );

    let full_out = tmp.join("full.bam");
    let region_out = tmp.join("region.bam");
    assert_eq!(
        exit_to_u8(cat::main(&argv(
            "cat",
            &["-o", full_out.to_str().unwrap(), indexed.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(cat::main(&argv(
            "cat",
            &[
                "-r",
                "17:1-2000",
                "-o",
                region_out.to_str().unwrap(),
                indexed.to_str().unwrap()
            ]
        ))),
        0
    );

    let full_count = tmp.join("full.count.txt");
    let region_count = tmp.join("region.count.txt");
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &[
                "-c",
                "-o",
                full_count.to_str().unwrap(),
                full_out.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &[
                "-c",
                "-o",
                region_count.to_str().unwrap(),
                region_out.to_str().unwrap()
            ]
        ))),
        0
    );
    let full: u64 = std::fs::read_to_string(&full_count)
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    let restricted: u64 = std::fs::read_to_string(&region_count)
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert!(restricted > 0);
    assert!(restricted < full);
}

#[test]
fn calmd_sam_input_recomputes_md_and_nm_tags() {
    let tmp = tmp_dir("calmd-md-nm");
    let sam = tmp.join("in.sam");
    let reference = tmp.join("ref.fa");
    let out = tmp.join("out.sam");
    std::fs::write(&reference, ">chr1\nACGTACGTACGT\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:12\n",
            "r1\t0\tchr1\t1\t60\t4M1I2M1D2M\t*\t0\t0\tACGTTCGAC\t!!!!!!!!!\tNM:i:99\tMD:Z:0A0\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(calmd::main(&argv(
            "calmd",
            &[
                "--no-PG",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
                reference.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let record = text.lines().find(|line| line.starts_with("r1\t")).unwrap();
    assert!(record.contains("\tNM:i:6"));
    assert!(record.contains("\tMD:Z:4A0C0^G0T0A0"));
    assert!(!record.contains("NM:i:99"));
    assert!(!record.contains("MD:Z:0A0"));
}

#[test]
fn calmd_bam_input_recomputes_md_and_nm_tags_to_sam_output() {
    let tmp = tmp_dir("calmd-bam-md-nm");
    let sam = tmp.join("in.sam");
    let bam = tmp.join("in.bam");
    let reference = tmp.join("ref.fa");
    let out = tmp.join("out.sam");
    std::fs::write(&reference, ">chr1\nACGTACGTACGT\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:12\n",
            "r1\t0\tchr1\t1\t60\t4M1I2M1D2M\t*\t0\t0\tACGTTCGAC\t!!!!!!!!!\tNM:i:99\tMD:Z:0A0\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &["-b", "-o", bam.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(calmd::main(&argv(
            "calmd",
            &[
                "--no-PG",
                "-o",
                out.to_str().unwrap(),
                bam.to_str().unwrap(),
                reference.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let record = text.lines().find(|line| line.starts_with("r1\t")).unwrap();
    assert!(record.contains("\tNM:i:6"));
    assert!(record.contains("\tMD:Z:4A0C0^G0T0A0"));
    assert!(!record.contains("NM:i:99"));
    assert!(!record.contains("MD:Z:0A0"));
}

#[test]
fn calmd_max_nm_masks_matching_bases_and_qualities() {
    let tmp = tmp_dir("calmd-max-nm");
    let sam = tmp.join("in.sam");
    let reference = tmp.join("ref.fa");
    let out = tmp.join("out.sam");
    std::fs::write(&reference, ">chr1\nACGTACGT\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "low\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII\tNM:i:99\tMD:Z:0A0\n",
            "high\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGA\tIIIIIIII\tNM:i:99\tMD:Z:0A0\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(calmd::main(&argv(
            "calmd",
            &[
                "--no-PG",
                "-n2",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
                reference.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let low = text.lines().find(|line| line.starts_with("low\t")).unwrap();
    let high = text
        .lines()
        .find(|line| line.starts_with("high\t"))
        .unwrap();

    assert!(low.contains("\tACGTTCGT\tIIIIIIII\t"));
    assert!(low.contains("\tNM:i:1"));
    assert!(low.contains("\tMD:Z:4A3"));
    assert!(high.contains("\tNNNNTNNA\t!!!!I!!I\t"));
    assert!(high.contains("\tNM:i:2"));
    assert!(high.contains("\tMD:Z:4A2T0"));
    assert!(!high.contains("NM:i:99"));
    assert!(!high.contains("MD:Z:0A0"));
}

#[test]
fn calmd_dash_e_changes_matching_bases_for_mapped_and_unmapped_records() {
    let tmp = tmp_dir("calmd-use-equal");
    let sam = tmp.join("in.sam");
    let reference = tmp.join("ref.fa");
    let out = tmp.join("out.sam");
    std::fs::write(&reference, ">chr1\nACGTACGT\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "mapped\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII\tNM:i:99\tMD:Z:0A0\n",
            "unmapped\t4\tchr1\t1\t0\t4M4S\t*\t0\t0\tACGTAAAA\t!!!!!!!!\tNM:i:99\tMD:Z:0A0\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(calmd::main(&argv(
            "calmd",
            &[
                "--no-PG",
                "-e",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
                reference.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let mapped = text
        .lines()
        .find(|line| line.starts_with("mapped\t"))
        .unwrap();
    let unmapped = text
        .lines()
        .find(|line| line.starts_with("unmapped\t"))
        .unwrap();

    assert!(mapped.contains("\t====T===\t"));
    assert!(mapped.contains("\tNM:i:1"));
    assert!(mapped.contains("\tMD:Z:4A3"));
    assert!(unmapped.contains("\t====AAAA\t"));
    assert!(unmapped.contains("\tNM:i:99"));
    assert!(unmapped.contains("\tMD:Z:0A0"));
}

#[test]
fn calmd_dash_q_bins_base_qualities() {
    let tmp = tmp_dir("calmd-bin-qual");
    let sam = tmp.join("in.sam");
    let reference = tmp.join("ref.fa");
    let out = tmp.join("out.sam");
    std::fs::write(&reference, ">chr1\nACGTACGT\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\t!+5?IS]g\tNM:i:99\tMD:Z:0A0\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(calmd::main(&argv(
            "calmd",
            &[
                "--no-PG",
                "-q",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
                reference.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let record = text.lines().find(|line| line.starts_with("r1\t")).unwrap();
    assert!(record.contains("\t!2<FPZdn\t"));
    assert!(record.contains("\tNM:i:1"));
    assert!(record.contains("\tMD:Z:4A3"));
}

#[test]
fn calmd_dash_cap_n_preserves_existing_md_nm_tags() {
    let tmp = tmp_dir("calmd-no-md-nm");
    let sam = tmp.join("in.sam");
    let reference = tmp.join("ref.fa");
    let out = tmp.join("out.sam");
    std::fs::write(&reference, ">chr1\nACGTACGT\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII\tNM:i:99\tMD:Z:0A0\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(calmd::main(&argv(
            "calmd",
            &[
                "--no-PG",
                "-N",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
                reference.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let record = text.lines().find(|line| line.starts_with("r1\t")).unwrap();
    assert!(record.contains("\tNM:i:99"));
    assert!(record.contains("\tMD:Z:0A0"));
}

#[test]
fn calmd_cap_mapping_quality_uses_sam_cap_mapq() {
    let tmp = tmp_dir("calmd-cap-mapq");
    let sam = tmp.join("in.sam");
    let reference = tmp.join("ref.fa");
    let out = tmp.join("out.sam");
    std::fs::write(&reference, ">sq0\nACGTACGT\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@SQ\tSN:sq0\tLN:8\n",
            "perfect\t0\tsq0\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
            "mismatch\t0\tsq0\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII\n",
            "softclip\t0\tsq0\t1\t60\t2S6M\t*\t0\t0\tTTACGTAC\tIIIIIIII\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(calmd::main(&argv(
            "calmd",
            &[
                "--no-PG",
                "-C40",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
                reference.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let records: Vec<Vec<&str>> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').collect())
        .collect();

    assert_eq!(records[0][4], "40");
    assert_eq!(records[1][4], "28");
    assert_eq!(records[2][4], "31");
    assert!(records[0].contains(&"NM:i:0"));
    assert!(records[1].contains(&"NM:i:1"));
    assert!(records[2].contains(&"NM:i:0"));
}

#[test]
fn calmd_writes_cram_output_with_reference() {
    let tmp = tmp_dir("calmd-cram-output");
    let sam = tmp.join("in.sam");
    let reference = tmp.join("ref.fa");
    let out = tmp.join("out.cram");
    let out_by_fmt = tmp.join("out.byfmt");
    std::fs::write(&reference, ">chr1\nACGTACGTACGT\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:12\n",
            "r1\t0\tchr1\t1\t60\t4M1I2M1D2M\t*\t0\t0\tACGTTCGAC\t!!!!!!!!!\tNM:i:99\tMD:Z:0A0\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(calmd::main(&argv(
            "calmd",
            &[
                "--no-PG",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
                reference.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(calmd::main(&argv(
            "calmd",
            &[
                "--no-PG",
                "--output-fmt=cram",
                "-o",
                out_by_fmt.to_str().unwrap(),
                sam.to_str().unwrap(),
                reference.to_str().unwrap(),
            ]
        ))),
        0
    );

    for cram in [&out, &out_by_fmt] {
        let text = htslib_rs::alignment_compat::view_cram_as_sam_text_with_reference(
            std::io::Cursor::new(std::fs::read(cram).unwrap()),
            &reference,
            None,
        )
        .unwrap();
        let record = text.lines().find(|line| line.starts_with("r1\t")).unwrap();
        assert!(record.contains("\tNM:i:6"));
        assert!(record.contains("\tMD:Z:4A0C0^G0T0A0"));
        assert!(!record.contains("NM:i:99"));
        assert!(!record.contains("MD:Z:0A0"));
    }
}

#[test]
fn calmd_dash_d_keeps_only_rg_aux_tag() {
    let tmp = tmp_dir("calmd-drop-bq");
    let sam = tmp.join("in.sam");
    let reference = tmp.join("ref.fa");
    let out = tmp.join("out.sam");
    std::fs::write(&reference, ">chr1\nACGTACGTACGT\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:12\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:g1\tBQ:Z:abcd\tXX:i:7\tNM:i:0\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(calmd::main(&argv(
            "calmd",
            &[
                "--no-PG",
                "-d",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
                reference.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let record = text.lines().find(|line| line.starts_with("r1\t")).unwrap();
    assert!(record.ends_with("\tRG:Z:g1"));
    assert!(!record.contains("\tBQ:Z:"));
    assert!(!record.contains("\tXX:i:"));
    assert!(!record.contains("\tNM:i:"));
}

/// Mirrors upstream `test_calmd`: `calmd -uAr mpileup.1.sam
/// mpileup.ref.fa` must emit a BGZF (BAM) stream. We additionally
/// assert the BAM round-trips with the input record count via `view`,
/// and that the glued `-uAr` cluster is split like `getopt`.
#[test]
fn calmd_dash_u_a_r_emits_bgzf_bam_like_upstream() {
    use samtools_rs::commands::view;
    let dat = fixtures_dir().join("dat");
    let sam = dat.join("mpileup.1.sam");
    let reference = dat.join("mpileup.ref.fa");
    let tmp = tmp_dir("calmd-uAr");
    let out = tmp.join("out.bam");

    assert_eq!(
        exit_to_u8(calmd::main(&argv(
            "calmd",
            &[
                "--no-PG",
                "-uAr",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
                reference.to_str().unwrap(),
            ]
        ))),
        0
    );

    // BGZF magic (gzip \x1f\x8b + BAM's FEXTRA), the exact upstream
    // `test_calmd` acceptance check.
    let bytes = std::fs::read(&out).unwrap();
    assert!(
        bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b,
        "calmd -uAr output is not BGZF-compressed"
    );

    // Record count is preserved through the SAM->BAQ->BAM path.
    let in_records = std::fs::read_to_string(&sam)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .count();
    let counted = tmp.join("count.txt");
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &["-c", "-o", counted.to_str().unwrap(), out.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&counted).unwrap().trim(),
        in_records.to_string()
    );
}

#[test]
fn calmd_baq_accepts_bam_and_cram_input() {
    let dat = fixtures_dir().join("dat");
    let sam = dat.join("mpileup.1.sam");
    let reference = dat.join("mpileup.ref.fa");
    let tmp = tmp_dir("calmd-baq-bam-cram");
    let bam = tmp.join("in.bam");
    let cram = tmp.join("in.cram");
    let sam_out = tmp.join("sam.out");
    let bam_out = tmp.join("bam.out");
    let cram_out = tmp.join("cram.out");

    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &[
                "--no-PG",
                "-b",
                "-o",
                bam.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &[
                "--no-PG",
                "-C",
                "-T",
                reference.to_str().unwrap(),
                "-o",
                cram.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    for (input, output) in [(&sam, &sam_out), (&bam, &bam_out), (&cram, &cram_out)] {
        assert_eq!(
            exit_to_u8(calmd::main(&argv(
                "calmd",
                &[
                    "--no-PG",
                    "-r",
                    "-o",
                    output.to_str().unwrap(),
                    input.to_str().unwrap(),
                    reference.to_str().unwrap(),
                ]
            ))),
            0
        );
    }

    let sam_records = non_header_lines(&std::fs::read_to_string(&sam_out).unwrap());
    assert!(!sam_records.is_empty());
    assert!(
        sam_records.iter().any(|line| line.contains("\tBQ:Z:")),
        "fixture should exercise recalculated BAQ tags"
    );

    let sam_qnames: Vec<&str> = sam_records
        .iter()
        .map(|line| line.split('\t').next().unwrap())
        .collect();
    for (label, output) in [("BAM", &bam_out), ("CRAM", &cram_out)] {
        let records = non_header_lines(&std::fs::read_to_string(output).unwrap());
        assert_eq!(records.len(), sam_records.len(), "{label} record count");
        let qnames: Vec<&str> = records
            .iter()
            .map(|line| line.split('\t').next().unwrap())
            .collect();
        assert_eq!(qnames, sam_qnames, "{label} record order");
        assert!(
            records.iter().any(|line| line.contains("\tBQ:Z:")),
            "{label} output should contain recalculated BAQ tags"
        );
        assert!(
            records.iter().any(|line| line.contains("\tNM:i:"))
                && records.iter().any(|line| line.contains("\tMD:Z:")),
            "{label} output should contain recalculated MD/NM tags"
        );
    }
}

#[test]
fn markdup_sam_input_flags_duplicates_keeping_highest_score() {
    use samtools_rs::commands::markdup;
    let tmp = tmp_dir("markdup-sam");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    // Upstream keeps the read with the higher `calc_score` (sum of base
    // quals >= 15), not the higher MAPQ. `low` quals are all phred 0;
    // `high` quals are all phred 40, so `high` wins regardless of MAPQ.
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "low\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\tIIII\n",
            "reverse\t16\tchr1\t1\t30\t4M\t*\t0\t0\tCCCC\t$$$$\n",
            "unique\t0\tchr1\t2\t30\t4M\t*\t0\t0\tGGGG\t....\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
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

    let text = std::fs::read_to_string(&out).unwrap();
    let high = text.lines().find(|l| l.starts_with("high\t")).unwrap();
    let low = text.lines().find(|l| l.starts_with("low\t")).unwrap();
    let reverse = text.lines().find(|l| l.starts_with("reverse\t")).unwrap();
    let unique = text.lines().find(|l| l.starts_with("unique\t")).unwrap();

    let flag_of = |line: &str| line.split('\t').nth(1).unwrap().parse::<u32>().unwrap();
    // 0x400 == BAM_FDUP
    assert_eq!(
        flag_of(high) & 0x400,
        0,
        "highest calc_score in group keeps primary"
    );
    assert_eq!(flag_of(low) & 0x400, 0x400, "low-score duplicate flagged");
    assert_eq!(
        flag_of(reverse) & 0x400,
        0,
        "different strand is not a duplicate"
    );
    assert_eq!(
        flag_of(unique) & 0x400,
        0,
        "different position is not a duplicate"
    );
}

#[test]
fn markdup_accepts_reference_backed_cram_input_and_output() {
    use samtools_rs::commands::markdup;

    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    set_current_global_args(SamGlobalArgs::default());

    let tmp = tmp_dir("markdup-cram");
    let reference = tmp.join("ref.fa");
    let sam = tmp.join("in.sam");
    let bam = tmp.join("in.bam");
    let cram = tmp.join("in.cram");
    let out_sam = tmp.join("out.sam");
    let out_cram = tmp.join("out.cram");

    std::fs::write(&reference, ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:32\n",
            "low\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
            "unique\t0\tchr1\t2\t30\t4M\t*\t0\t0\tCGTA\t####\n",
        ),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference(
        &bam,
        &reference,
        std::fs::File::create(&cram).unwrap(),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "--no-PG",
                "-T",
                reference.to_str().unwrap(),
                "-O",
                "sam",
                "-o",
                out_sam.to_str().unwrap(),
                cram.to_str().unwrap(),
            ]
        ))),
        0
    );
    let text = std::fs::read_to_string(&out_sam).unwrap();
    let flag_of = |text: &str, name: &str| -> u32 {
        text.lines()
            .find(|line| line.starts_with(&format!("{name}\t")))
            .unwrap()
            .split('\t')
            .nth(1)
            .unwrap()
            .parse()
            .unwrap()
    };
    assert_eq!(flag_of(&text, "high") & 0x400, 0);
    assert_eq!(flag_of(&text, "low") & 0x400, 0x400);

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "--no-PG",
                "-T",
                reference.to_str().unwrap(),
                "--output-fmt=cram",
                sam.to_str().unwrap(),
                out_cram.to_str().unwrap(),
            ]
        ))),
        0
    );
    let cram_text = htslib_rs::alignment_compat::view_cram_as_sam_text_with_reference(
        std::io::Cursor::new(std::fs::read(&out_cram).unwrap()),
        &reference,
        None,
    )
    .unwrap();
    assert_eq!(flag_of(&cram_text, "high") & 0x400, 0);
    assert_eq!(flag_of(&cram_text, "low") & 0x400, 0x400);

    set_current_global_args(SamGlobalArgs::default());
}

#[test]
fn markdup_barcode_tag_separates_duplicate_groups() {
    use samtools_rs::commands::markdup;
    let tmp = tmp_dir("markdup-barcode");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "aa_low\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\tBC:Z:AA\n",
            "aa_high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\tIIII\tBC:Z:AA\n",
            "bb\t0\tchr1\t1\t20\t4M\t*\t0\t0\tCCCC\t$$$$\tBC:Z:BB\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "-O",
                "sam",
                "-b",
                "BC",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let aa_high = text.lines().find(|l| l.starts_with("aa_high\t")).unwrap();
    let aa_low = text.lines().find(|l| l.starts_with("aa_low\t")).unwrap();
    let bb = text.lines().find(|l| l.starts_with("bb\t")).unwrap();

    let flag_of = |line: &str| line.split('\t').nth(1).unwrap().parse::<u32>().unwrap();
    assert_eq!(flag_of(aa_high) & 0x400, 0);
    assert_eq!(flag_of(aa_low) & 0x400, 0x400);
    assert_eq!(flag_of(bb) & 0x400, 0);
}

#[test]
fn markdup_r_removes_duplicates_from_output() {
    use samtools_rs::commands::markdup;
    let tmp = tmp_dir("markdup-remove");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "low\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "-r",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.lines().any(|l| l.starts_with("high\t")));
    assert!(!text.lines().any(|l| l.starts_with("low\t")));
}

#[test]
fn markdup_c_clears_existing_duplicate_marks_and_s_accepts_supplementary_mode() {
    use samtools_rs::commands::markdup;
    let tmp = tmp_dir("markdup-clear");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "previous\t1024\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tdo:Z:old\tdt:Z:LB\n",
            "unique\t0\tchr1\t10\t60\t4M\t*\t0\t0\tTGCA\t####\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "-c",
                "-S",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let previous = text
        .lines()
        .find(|line| line.starts_with("previous\t"))
        .unwrap();
    let flag = previous.split('\t').nth(1).unwrap().parse::<u32>().unwrap();
    assert_eq!(flag & 0x400, 0);
    assert!(!previous.contains("\tdo:Z:"));
    assert!(!previous.contains("\tdt:Z:"));
}

#[test]
fn markdup_t_adds_duplicate_origin_tag() {
    use samtools_rs::commands::markdup;
    let tmp = tmp_dir("markdup-origin-tag");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "low\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "-t",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let high = text
        .lines()
        .find(|line| line.starts_with("high\t"))
        .unwrap();
    let low = text.lines().find(|line| line.starts_with("low\t")).unwrap();
    assert!(!high.contains("\tdo:Z:"));
    assert!(low.contains("\tdo:Z:high"));
}

#[test]
fn markdup_d_adds_duplicate_type_tags() {
    use samtools_rs::commands::markdup;
    let tmp = tmp_dir("markdup-duplicate-type");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "INST:1:FC:1:1101:100:100\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "INST:1:FC:1:1101:105:108\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "INST:1:FC:1:1102:105:108\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "-d",
                "10",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let original = text
        .lines()
        .find(|line| line.starts_with("INST:1:FC:1:1101:100:100\t"))
        .unwrap();
    let optical = text
        .lines()
        .find(|line| line.starts_with("INST:1:FC:1:1101:105:108\t"))
        .unwrap();
    let library = text
        .lines()
        .find(|line| line.starts_with("INST:1:FC:1:1102:105:108\t"))
        .unwrap();
    assert!(!original.contains("\tdt:Z:"));
    assert!(optical.contains("\tdt:Z:SQ"));
    assert!(library.contains("\tdt:Z:LB"));
}

#[test]
fn markdup_c_removes_stale_duplicate_tags_before_t_retags_duplicates() {
    use samtools_rs::commands::markdup;
    let tmp = tmp_dir("markdup-clear-retag");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "low\t1024\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\tdo:Z:stale\tdt:Z:LB\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "-c",
                "-t",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let low = text.lines().find(|line| line.starts_with("low\t")).unwrap();
    assert!(low.contains("\tdo:Z:high"));
    assert!(!low.contains("\tdo:Z:stale"));
    assert!(!low.contains("\tdt:Z:"));
}

#[test]
fn markdup_include_fails_controls_qcfail_duplicate_marking() {
    use samtools_rs::commands::markdup;
    let tmp = tmp_dir("markdup-include-fails");
    let sam = tmp.join("in.sam");
    let default_out = tmp.join("default.sam");
    let include_out = tmp.join("include.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "pass\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "fail\t512\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "-O",
                "sam",
                "-o",
                default_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "--include-fails",
                "-O",
                "sam",
                "-o",
                include_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let default_text = std::fs::read_to_string(&default_out).unwrap();
    let include_text = std::fs::read_to_string(&include_out).unwrap();
    let default_fail = default_text
        .lines()
        .find(|line| line.starts_with("fail\t"))
        .unwrap();
    let include_fail = include_text
        .lines()
        .find(|line| line.starts_with("fail\t"))
        .unwrap();
    let default_flag = default_fail
        .split('\t')
        .nth(1)
        .unwrap()
        .parse::<u32>()
        .unwrap();
    let include_flag = include_fail
        .split('\t')
        .nth(1)
        .unwrap()
        .parse::<u32>()
        .unwrap();
    assert_eq!(default_flag & 0x400, 0);
    assert_eq!(include_flag & 0x400, 0x400);
}

#[test]
fn markdup_mode_accepts_valid_values_and_rejects_invalid() {
    use samtools_rs::commands::markdup;
    let tmp = tmp_dir("markdup-mode");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "low\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "--mode",
                "s",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "-m",
                "bad",
                "-O",
                "sam",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        1
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let low = text.lines().find(|line| line.starts_with("low\t")).unwrap();
    let flag = low.split('\t').nth(1).unwrap().parse::<u32>().unwrap();
    assert_eq!(flag & 0x400, 0x400);
}

#[test]
fn markdup_propagates_duplicate_flag_to_supplementary_records() {
    use samtools_rs::commands::markdup;
    let tmp = tmp_dir("markdup-supplementary");
    let sam = tmp.join("in.sam");
    let flagged_out = tmp.join("flagged.sam");
    let removed_out = tmp.join("removed.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:1000\n",
            "keep\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\tIIII\n",
            "dup\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\tSA:Z:chr1,20,+,4M,10,0;\n",
            "dup\t2048\tchr1\t20\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "keep\t2048\tchr1\t30\t10\t4M\t*\t0\t0\tTGCA\tIIII\n",
        ),
    )
    .unwrap();

    // Upstream only propagates to supplementary/secondary records with
    // `-S`, and only when the marked-duplicate read carries `SA`/`XA` or
    // an unmapped mate (here `dup` has an `SA` tag).
    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "-S",
                "-O",
                "sam",
                "-o",
                flagged_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    let flagged_text = std::fs::read_to_string(&flagged_out).unwrap();
    let flag_of = |text: &str, name: &str, original_flag: u32| -> u32 {
        text.lines()
            .find_map(|line| {
                let mut fields = line.split('\t');
                let qname = fields.next()?;
                let flag = fields.next()?.parse::<u32>().ok()?;
                (qname == name && flag & !0x400 == original_flag).then_some(flag)
            })
            .expect("record present")
    };

    assert_eq!(flag_of(&flagged_text, "dup", 0) & 0x400, 0x400);
    assert_eq!(flag_of(&flagged_text, "dup", 2048) & 0x400, 0x400);
    assert_eq!(flag_of(&flagged_text, "keep", 0) & 0x400, 0);
    assert_eq!(flag_of(&flagged_text, "keep", 2048) & 0x400, 0);

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
            &[
                "-S",
                "-r",
                "-O",
                "sam",
                "-o",
                removed_out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    let removed_text = std::fs::read_to_string(&removed_out).unwrap();
    assert!(!removed_text.lines().any(|line| line.starts_with("dup\t")));
    assert!(
        removed_text
            .lines()
            .any(|line| line.starts_with("keep\t0\t"))
    );
    assert!(
        removed_text
            .lines()
            .any(|line| line.starts_with("keep\t2048\t"))
    );
}

#[test]
fn ampliconstats_matches_upstream_test_ampliconstats_fixtures() {
    use samtools_rs::commands::ampliconstats;
    // Byte-exact (modulo the harness-stripped version/command-line lines)
    // vs the entire upstream test_ampliconstats harness: single-ref
    // multi-file (-S -t 50 -d 1,20,100), and the multi-ref/partial
    // single-file -c 0 cases.
    let astats = fixtures_dir().join("ampliconstats");
    let aclip = fixtures_dir().join("ampliconclip");
    let tmp = tmp_dir("ampliconstats-fixtures");
    let strip = |s: &str| -> String {
        s.lines()
            .filter(|l| !l.contains("Samtools version") && !l.contains("Command line"))
            .map(|l| format!("{l}\n"))
            .collect()
    };
    let p = |b: &std::path::Path| b.to_str().unwrap().to_string();
    let cases: Vec<(Vec<String>, &str)> = vec![
        (
            vec![
                "-S".into(),
                "-t".into(),
                "50".into(),
                "-d".into(),
                "1,20,100".into(),
                p(&aclip.join("ac_test.bed")),
                p(&aclip.join("1_hard_clipped.expected.sam")),
                p(&aclip.join("1_soft_clipped.expected.sam")),
                p(&aclip.join("1_soft_clipped_strand.expected.sam")),
                p(&aclip.join("2_both_clipped.expected.sam")),
            ],
            "stats.expected.txt",
        ),
        (
            vec![
                "-c".into(),
                "0".into(),
                p(&aclip.join("multi_ref.bed")),
                p(&astats.join("mixed_clipped.sam")),
            ],
            "stats_mixed.expected.txt",
        ),
        (
            vec![
                "-c".into(),
                "0".into(),
                p(&aclip.join("ac_test.bed")),
                p(&astats.join("mixed_clipped.sam")),
            ],
            "stats_partial.expected.txt",
        ),
    ];
    for (rest, expected) in cases {
        let out = tmp.join(expected);
        let mut v: Vec<String> = vec!["-o".into(), p(&out)];
        v.extend(rest);
        assert_eq!(
            exit_to_u8(ampliconstats::main(&argv(
                "ampliconstats",
                &v.iter().map(String::as_str).collect::<Vec<_>>()
            ))),
            0,
            "ampliconstats {expected}"
        );
        assert_eq!(
            strip(&std::fs::read_to_string(&out).unwrap()),
            strip(&std::fs::read_to_string(astats.join(expected)).unwrap()),
            "ampliconstats {expected} byte-exact",
        );
    }
}

#[test]
fn ampliconstats_use_sample_name_uses_first_read_group_sample() {
    use samtools_rs::commands::ampliconstats;

    let tmp = tmp_dir("ampliconstats-sample-name");
    let bed = tmp.join("primers.bed");
    let sam = tmp.join("input_name.sam");
    let out = tmp.join("out.txt");

    std::fs::write(
        &bed,
        concat!("chr1\t0\t10\tleft\t0\t+\n", "chr1\t90\t100\tright\t0\t-\n",),
    )
    .unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "@RG\tID:first\tSM:sample_from_header\n",
            "@RG\tID:second\tSM:ignored_sample\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(ampliconstats::main(&argv(
            "ampliconstats",
            &[
                "-s",
                "-o",
                out.to_str().unwrap(),
                bed.to_str().unwrap(),
                sam.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.contains("FSS\tsample_from_header\tchr1\traw total sequences:\t0"));
    assert!(text.contains("FREADS\tsample_from_header\t0"));
    assert!(!text.contains("FSS\tinput_name\t"));
    assert!(!text.contains("ignored_sample"));
}

#[test]
fn ampliconstats_accepts_bam_and_reference_backed_cram_input() {
    use samtools_rs::commands::ampliconstats;

    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    set_current_global_args(SamGlobalArgs::default());

    let tmp = tmp_dir("ampliconstats-bam-cram");
    let bed = tmp.join("primers.bed");
    let bam_out = tmp.join("bam.txt");
    let cram_out = tmp.join("cram.txt");
    let hts = htslib_fixtures_dir();
    let bam = hts.join("range.bam");
    let cram = hts.join("range.cram");
    let reference = hts.join("ce.fa");

    std::fs::write(
        &bed,
        concat!(
            "CHROMOSOME_II\t2960\t2970\tleft\t0\t+\n",
            "CHROMOSOME_II\t3000\t3010\tright\t0\t-\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(ampliconstats::main(&argv(
            "ampliconstats",
            &[
                "-c",
                "0",
                "-o",
                bam_out.to_str().unwrap(),
                bed.to_str().unwrap(),
                bam.to_str().unwrap(),
            ]
        ))),
        0
    );
    let bam_text = std::fs::read_to_string(&bam_out).unwrap();
    assert!(bam_text.contains("FSS\trange\tCHROMOSOME_II\traw total sequences:"));

    set_current_global_args(SamGlobalArgs {
        reference: Some(reference),
        ..SamGlobalArgs::default()
    });
    assert_eq!(
        exit_to_u8(ampliconstats::main(&argv(
            "ampliconstats",
            &[
                "-s",
                "-c",
                "0",
                "-o",
                cram_out.to_str().unwrap(),
                bed.to_str().unwrap(),
                cram.to_str().unwrap(),
            ]
        ))),
        0
    );
    set_current_global_args(SamGlobalArgs::default());

    let cram_text = std::fs::read_to_string(cram_out).unwrap();
    assert!(cram_text.contains("FSS\tERS225193\tCHROMOSOME_II\traw total sequences:"));
}

#[test]
fn stats_matches_upstream_stat_fixtures() {
    use samtools_rs::commands::stats;
    // Byte-exact end to end vs samtools' own `test/stat/*` expected
    // output (modulo the three harness-stripped header lines, i.e.
    // `| tail -n+4`): CHK, all SN lines, FFQ/LFQ, MPC, GCF/GCL,
    // GCC/GCT/FBC/FTC/LBC/LTC, IS, RL/FRL/LRL, MAPQ, the ID/IC indel
    // distribution + per-cycle rows, COV, and the GC-depth row. Covers
    // the plain map, equal/full-seq, insertion (ID/IC + bases-cigar),
    // `-i 0`, and secondary-read fixtures.
    let stat = fixtures_dir().join("stat");
    let tmp = tmp_dir("stats-fixtures");
    let reference = stat.join("test.fa");
    let p = |b: &std::path::Path| b.to_str().unwrap().to_string();
    // `tail -n+4`: drop the produced-by / contains / command-line lines.
    let strip = |s: &str| -> String { s.lines().skip(3).map(|l| format!("{l}\n")).collect() };
    // (sam, extra args, expected). Covers plain map, equal/full-seq,
    // X-cigar (MPC reference mismatch), insertion (ID/IC + bases-cigar),
    // `-i 0`, supplementary (MPC + supp aux), and secondary fixtures.
    let cases: [(&str, &[&str], &str); 8] = [
        ("1_map_cigar.sam", &[], "1.stats.expected"),
        ("2_equal_cigar_full_seq.sam", &[], "2.stats.expected"),
        ("3_map_cigar_equal_seq.sam", &[], "3.stats.expected"),
        ("4_X_cigar_full_seq.sam", &[], "4.stats.expected"),
        ("5_insert_cigar.sam", &[], "5.stats.expected"),
        ("5_insert_cigar.sam", &["-i", "0"], "6.stats.expected"),
        ("7_supp.sam", &[], "7.stats.expected"),
        ("8_secondary.sam", &[], "8.stats.expected"),
    ];
    for (sam, extra, expected) in cases {
        let out = tmp.join(expected);
        let mut v: Vec<String> = vec!["-r".into(), p(&reference), "-o".into(), p(&out)];
        v.extend(extra.iter().map(|s| s.to_string()));
        v.push(p(&stat.join(sam)));
        assert_eq!(
            exit_to_u8(stats::main(&argv(
                "stats",
                &v.iter().map(String::as_str).collect::<Vec<_>>()
            ))),
            0,
            "stats {expected}"
        );
        assert_eq!(
            strip(&std::fs::read_to_string(&out).unwrap()),
            std::fs::read_to_string(stat.join(expected)).unwrap(),
            "stats {expected} byte-exact",
        );
    }

    // stat/15: unpaired read with a 60000D deletion against the
    // mpileup `ce.fa` reference — exercises the upstream `order`
    // (unpaired => first fragment) and the `nindels` (>300) ID cap.
    {
        let ce_fa = fixtures_dir().join("mpileup").join("ce.fa");
        let out = tmp.join("15.stats.expected");
        assert_eq!(
            exit_to_u8(stats::main(&argv(
                "stats",
                &[
                    "-r",
                    &p(&ce_fa),
                    "-o",
                    &p(&out),
                    &p(&stat.join("15.big_del.sam")),
                ]
            ))),
            0,
            "stats 15"
        );
        assert_eq!(
            strip(&std::fs::read_to_string(&out).unwrap()),
            std::fs::read_to_string(stat.join("15.stats.expected")).unwrap(),
            "stats 15 byte-exact",
        );
    }

    // stat/13: barcode BCC/QTQ (and OXC/BZQ) sections, no reference.
    for (sam, expected) in [
        ("13_barcodes_ok.sam", "13.barcodes.bc.ok.expected"),
        ("13_barcodes_ok_ox_bz.sam", "13.barcodes.ox.ok.expected"),
    ] {
        let out = tmp.join(expected);
        assert_eq!(
            exit_to_u8(stats::main(&argv(
                "stats",
                &["-o", &p(&out), &p(&stat.join(sam))]
            ))),
            0,
            "stats {sam}"
        );
        assert_eq!(
            strip(&std::fs::read_to_string(&out).unwrap()),
            std::fs::read_to_string(stat.join(expected)).unwrap(),
            "stats {sam} byte-exact",
        );
    }

    // stat/12: overlapping-pairs BAM with `-t` BED. The overlap
    // variants exercise the integer-halved insert-size avg/sd; the
    // `-p`/--remove-overlaps variants the paired-overlap chunk
    // subtraction (bases-mapped-cigar + coverage) and the f32
    // error-rate cast.
    let cases12: [(&[&str], &str); 4] = [
        (&["12_3reads.bed"], "12.3reads.overlap.expected"),
        (&["-p", "12_3reads.bed"], "12.3reads.nooverlap.expected"),
        (&["12_2reads.bed"], "12.2reads.overlap.expected"),
        (&["-p", "12_2reads.bed"], "12.2reads.nooverlap.expected"),
    ];
    for (rest, expected) in cases12 {
        let (opt, bed) = if rest.len() == 2 {
            (Some("-p"), rest[1])
        } else {
            (None, rest[0])
        };
        let out = tmp.join(expected);
        let mut v: Vec<String> = vec!["-o".into(), p(&out)];
        if let Some(o) = opt {
            v.push(o.into());
        }
        v.push("-t".into());
        v.push(p(&stat.join(bed)));
        v.push(p(&stat.join("12_overlaps.bam")));
        assert_eq!(
            exit_to_u8(stats::main(&argv(
                "stats",
                &v.iter().map(String::as_str).collect::<Vec<_>>()
            ))),
            0,
            "stats {expected}"
        );
        assert_eq!(
            strip(&std::fs::read_to_string(&out).unwrap()),
            std::fs::read_to_string(stat.join(expected)).unwrap(),
            "stats {expected} byte-exact",
        );
    }

    // `--ref-stats` RFS section (upstream compares only `grep ^RFS`):
    // stat/16 no reference (GC/N = -1), stat/17 reference-backed GC/N
    // (plain and `--ref-stats-chunk -1` no-op), stat/18 positional
    // region, stat/19 `-t` target file (overlapping intervals merged,
    // including the harness' extra positional-region variant).
    let test1_fa = p(&stat.join("test1.fa"));
    let sam = p(&stat.join("11_target.sam"));
    let bam = p(&stat.join("11_target.bam"));
    let targets = p(&stat.join("11.stats.targets"));
    let rfs_cases: [(Vec<String>, &str); 6] = [
        (vec![sam.clone()], "16.stats.expected"),
        (
            vec!["-r".into(), test1_fa.clone(), sam.clone()],
            "17.stats.expected",
        ),
        (
            vec![
                "--ref-stats-chunk".into(),
                "-1".into(),
                "-r".into(),
                test1_fa.clone(),
                sam.clone(),
            ],
            "17.stats.expected",
        ),
        (
            vec![
                "-r".into(),
                test1_fa.clone(),
                bam.clone(),
                "alpha:10-20".into(),
            ],
            "18.stats.expected",
        ),
        (
            vec![
                "-r".into(),
                test1_fa.clone(),
                "-t".into(),
                targets.clone(),
                sam.clone(),
            ],
            "19.stats.expected",
        ),
        (
            vec![
                "-r".into(),
                test1_fa.clone(),
                "-t".into(),
                targets.clone(),
                bam.clone(),
                "ref1".into(),
            ],
            "19.stats.expected",
        ),
    ];
    for (rest, expected) in rfs_cases {
        let out = tmp.join(expected);
        let mut v: Vec<String> = vec!["--ref-stats".into(), "-o".into(), p(&out)];
        v.extend(rest);
        assert_eq!(
            exit_to_u8(stats::main(&argv(
                "stats",
                &v.iter().map(String::as_str).collect::<Vec<_>>()
            ))),
            0,
            "stats --ref-stats {expected}"
        );
        let rfs: String = std::fs::read_to_string(&out)
            .unwrap()
            .lines()
            .filter(|l| l.starts_with("RFS\t"))
            .map(|l| format!("{l}\n"))
            .collect();
        assert_eq!(
            rfs,
            std::fs::read_to_string(stat.join(expected)).unwrap(),
            "stats --ref-stats {expected} RFS byte-exact",
        );
    }

    // stat/11: `-t` target file on a SAM (streaming + region-clipped
    // bases-mapped-cigar with the init_regions overlap-merge), the same
    // expected via positional regions on the indexed BAM, and `-g 4`.
    let targets11 = p(&stat.join("11.stats.targets"));
    let sam11 = p(&stat.join("11_target.sam"));
    let bam11 = p(&stat.join("11_target.bam"));
    let region11_cases: [(Vec<String>, &str); 3] = [
        (
            vec!["-t".into(), targets11.clone(), sam11.clone()],
            "11.stats.expected",
        ),
        (
            vec![
                bam11.clone(),
                "ref1:10-24".into(),
                "ref1:30-46".into(),
                "ref1:39-56".into(),
            ],
            "11.stats.expected",
        ),
        (
            vec![
                "-g".into(),
                "4".into(),
                "-t".into(),
                targets11.clone(),
                sam11.clone(),
            ],
            "11.stats.g4.expected",
        ),
    ];
    for (rest, expected) in region11_cases {
        let out = tmp.join(expected);
        let mut v: Vec<String> = vec!["-o".into(), p(&out)];
        v.extend(rest);
        assert_eq!(
            exit_to_u8(stats::main(&argv(
                "stats",
                &v.iter().map(String::as_str).collect::<Vec<_>>()
            ))),
            0,
            "stats region {expected}"
        );
        assert_eq!(
            strip(&std::fs::read_to_string(&out).unwrap()),
            std::fs::read_to_string(stat.join(expected)).unwrap(),
            "stats region {expected} byte-exact",
        );
    }

    // `-S RG` split cases: stdout matches `<n>.stats.expected` and each
    // per-RG `<input>_<rg>.bamstat` matches its `.expected.bamstat`
    // (both modulo the three stripped header lines). The input is copied
    // into the temp dir so the side files land there.
    let split_cases: [(&str, &str, &[&str]); 2] = [
        ("1_map_cigar.sam", "9.stats.expected", &["s1_a_1"]),
        (
            "10_map_cigar.sam",
            "10.stats.expected",
            &["s1_a_1", "s1_b_1"],
        ),
    ];
    for (sam, expected, rgs) in split_cases {
        let local_sam = tmp.join(sam);
        std::fs::copy(stat.join(sam), &local_sam).unwrap();
        let out = tmp.join(expected);
        assert_eq!(
            exit_to_u8(stats::main(&argv(
                "stats",
                &[
                    "-S",
                    "RG",
                    "-r",
                    &p(&reference),
                    "-o",
                    &p(&out),
                    &p(&local_sam),
                ]
            ))),
            0,
            "stats -S RG {expected}"
        );
        assert_eq!(
            strip(&std::fs::read_to_string(&out).unwrap()),
            std::fs::read_to_string(stat.join(expected)).unwrap(),
            "stats -S RG {expected} stdout byte-exact",
        );
        for rg in rgs {
            let bamstat = tmp.join(format!("{sam}_{rg}.bamstat"));
            let exp = stat.join(format!("{sam}_{rg}.expected.bamstat"));
            assert_eq!(
                strip(&std::fs::read_to_string(&bamstat).unwrap()),
                std::fs::read_to_string(&exp).unwrap(),
                "stats -S RG {sam}_{rg}.bamstat byte-exact",
            );
        }
    }

    // `-I <rg>` read-group filter on an indexed BAM (stat/14), incl.
    // grp3 which matches no reads (empty FFQ/LFQ/GCF/GCL but their
    // comment headers are still emitted, as upstream).
    let bam = stat.join("11_target.bam");
    for rg in ["s1", "grp2", "grp3"] {
        let out = tmp.join(format!("14.rg.{rg}"));
        assert_eq!(
            exit_to_u8(stats::main(&argv(
                "stats",
                &["-I", rg, "-o", &p(&out), &p(&bam)]
            ))),
            0,
            "stats -I {rg}"
        );
        assert_eq!(
            strip(&std::fs::read_to_string(&out).unwrap()),
            std::fs::read_to_string(stat.join(format!("14.rg.{rg}.expected"))).unwrap(),
            "stats -I {rg} byte-exact",
        );
    }
}

#[test]
fn ampliconclip_matches_upstream_test_ampliconclip_fixtures() {
    use samtools_rs::commands::ampliconclip;
    // Byte-exact vs the entire upstream `test/ampliconclip` harness:
    // soft/hard clip, NM/MD deletion, OA tag, strand, filter-len,
    // fail-len, both-ends, multi-ref, total-hard-clip + unmap, and the
    // three primer-counts TSVs.
    let d = fixtures_dir().join("ampliconclip");
    let tmp = tmp_dir("ampliconclip-fixtures");
    let acb = d.join("ac_test.bed");
    let acb = acb.to_str().unwrap();
    let td = d.join("1_test_data.sam");
    let td = td.to_str().unwrap();
    let sam_cases: &[(&[&str], &str)] = &[
        (
            &["--no-PG", "--keep-tag", "--output-fmt=sam", "-b", acb, td],
            "1_soft_clipped",
        ),
        (
            &[
                "--no-PG",
                "--keep-tag",
                "--output-fmt=sam",
                "--hard-clip",
                "-b",
                acb,
                td,
            ],
            "1_hard_clipped",
        ),
        (
            &["--no-PG", "--output-fmt=sam", "-b", acb, td],
            "1_delete_tag",
        ),
        (
            &[
                "--no-PG",
                "--keep-tag",
                "--output-fmt=sam",
                "--original",
                "-b",
                acb,
                td,
            ],
            "1_original_tag",
        ),
        (
            &[
                "--no-PG",
                "--keep-tag",
                "--output-fmt=sam",
                "--strand",
                "-b",
                acb,
                td,
            ],
            "1_soft_clipped_strand",
        ),
        (
            &[
                "--no-PG",
                "--keep-tag",
                "--output-fmt=sam",
                "--strand",
                "--filter-len",
                "185",
                "-b",
                acb,
                td,
            ],
            "1_filter",
        ),
        (
            &[
                "--no-PG",
                "--keep-tag",
                "--output-fmt=sam",
                "--strand",
                "--fail-len",
                "185",
                "-b",
                acb,
                td,
            ],
            "1_fail",
        ),
    ];
    for (rest, expected) in sam_cases {
        let out = tmp.join(format!("{expected}.sam"));
        let mut v: Vec<String> = rest.iter().map(|s| s.to_string()).collect();
        v.extend(["-o".into(), out.to_str().unwrap().into()]);
        assert_eq!(
            exit_to_u8(ampliconclip::main(&argv(
                "ampliconclip",
                &v.iter().map(String::as_str).collect::<Vec<_>>()
            ))),
            0,
            "ampliconclip {expected}"
        );
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            std::fs::read_to_string(d.join(format!("{expected}.expected.sam"))).unwrap(),
            "ampliconclip {expected} byte-exact",
        );
    }

    // --both-ends, multi-ref, total-hard-clip.
    let extra: &[(&[&str], &str)] = &[
        (
            &[
                "--no-PG",
                "--keep-tag",
                "--output-fmt=sam",
                "--strand",
                "--both-ends",
                "-b",
                "ac_test.bed",
                "2_both_test_data.sam",
            ],
            "2_both_clipped",
        ),
        (
            &[
                "--no-PG",
                "--output-fmt=sam",
                "--keep-tag",
                "-b",
                "multi_ref.bed",
                "3_multi_ref_data.sam",
            ],
            "3_multi_ref_clip",
        ),
        (
            &[
                "--no-PG",
                "--output-fmt=sam",
                "--keep-tag",
                "--both-ends",
                "-b",
                "multi_ref.bed",
                "3_multi_ref_data.sam",
            ],
            "3_multi_ref_both_clip",
        ),
        (
            &[
                "--no-PG",
                "--output-fmt=sam",
                "--hard-clip",
                "-b",
                "ac_test2.bed",
                "4_total_hc_data.sam",
            ],
            "4_total_hc_data",
        ),
    ];
    for (rest, expected) in extra {
        let out = tmp.join(format!("{expected}.sam"));
        let mut v: Vec<String> = Vec::new();
        for tok in *rest {
            if tok.ends_with(".bed") || tok.ends_with(".sam") {
                v.push(d.join(tok).to_str().unwrap().to_string());
            } else {
                v.push((*tok).to_string());
            }
        }
        v.extend(["-o".into(), out.to_str().unwrap().into()]);
        assert_eq!(
            exit_to_u8(ampliconclip::main(&argv(
                "ampliconclip",
                &v.iter().map(String::as_str).collect::<Vec<_>>()
            ))),
            0,
            "ampliconclip {expected}"
        );
        let mut expected_text =
            std::fs::read_to_string(d.join(format!("{expected}.expected.sam"))).unwrap();
        if *expected == "3_multi_ref_both_clip" {
            // The dormant upstream expected file predates current
            // bam_ampliconclip behavior. Current C samtools clips these two
            // extra --both-ends multi-ref edges.
            expected_text = expected_text
                .replace(
                    "read_3\t161\tvir1\t132\t60\t210M\t=",
                    "read_3\t161\tvir1\t132\t60\t189M21S\t=",
                )
                .replace(
                    "read_4\t97\tvir1\t407\t60\t86S124M\t=",
                    "read_4\t97\tvir1\t411\t60\t90S120M\t=",
                );
        }
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            expected_text,
            "ampliconclip {expected} byte-exact",
        );
    }

    // primer-counts TSV outputs.
    let pc = tmp.join("pc.tsv");
    let pcs = pc.to_str().unwrap();
    let pc_cases: &[(&[&str], &str)] = &[
        (
            &[
                "--no-PG",
                "--keep-tag",
                "--primer-counts",
                pcs,
                "--output-fmt=sam",
                "-b",
                acb,
                td,
            ],
            "1_soft_clipped_primer_counts",
        ),
        (
            &[
                "--no-PG",
                "--keep-tag",
                "--output-fmt=sam",
                "--strand",
                "--primer-counts",
                pcs,
                "-b",
                acb,
                td,
            ],
            "1_soft_clipped_strand_primer_counts",
        ),
    ];
    for (rest, expected) in pc_cases {
        let out = tmp.join("oc.sam");
        let mut v: Vec<String> = rest.iter().map(|s| s.to_string()).collect();
        v.extend(["-o".into(), out.to_str().unwrap().into()]);
        assert_eq!(
            exit_to_u8(ampliconclip::main(&argv(
                "ampliconclip",
                &v.iter().map(String::as_str).collect::<Vec<_>>()
            ))),
            0,
            "ampliconclip {expected}"
        );
        assert_eq!(
            std::fs::read_to_string(&pc).unwrap(),
            std::fs::read_to_string(d.join(format!("{expected}.expected.tsv"))).unwrap(),
            "ampliconclip {expected} byte-exact",
        );
    }
}

#[test]
fn ampliconclip_accepts_reference_backed_cram_input_and_output() {
    use samtools_rs::commands::ampliconclip;

    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    set_current_global_args(SamGlobalArgs::default());

    let d = fixtures_dir().join("ampliconclip");
    let tmp = tmp_dir("ampliconclip-cram");
    let input_sam = d.join("1_test_data.sam");
    let input_bam = tmp.join("in.bam");
    let input_cram = tmp.join("in.cram");
    let reference = tmp.join("ref.fa");
    let baseline = tmp.join("baseline.sam");
    let cram_input_out = tmp.join("cram-input.sam");
    let out_cram = tmp.join("out.cram");
    let bed = d.join("ac_test.bed");

    std::fs::write(&reference, format!(">vir1\n{}\n", "N".repeat(800))).unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();
    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &input_sam,
        std::fs::File::create(&input_bam).unwrap(),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference(
        &input_bam,
        &reference,
        std::fs::File::create(&input_cram).unwrap(),
    )
    .unwrap();

    let base_args = [
        "--no-PG",
        "--keep-tag",
        "--output-fmt=sam",
        "-b",
        bed.to_str().unwrap(),
    ];
    let mut baseline_args: Vec<&str> = base_args.to_vec();
    baseline_args.extend([
        "-o",
        baseline.to_str().unwrap(),
        input_sam.to_str().unwrap(),
    ]);
    assert_eq!(
        exit_to_u8(ampliconclip::main(&argv("ampliconclip", &baseline_args))),
        0
    );

    let mut cram_input_args: Vec<&str> = base_args.to_vec();
    cram_input_args.extend([
        "-T",
        reference.to_str().unwrap(),
        "-o",
        cram_input_out.to_str().unwrap(),
        input_cram.to_str().unwrap(),
    ]);
    assert_eq!(
        exit_to_u8(ampliconclip::main(&argv("ampliconclip", &cram_input_args))),
        0
    );
    assert_eq!(
        non_header_lines(&std::fs::read_to_string(&cram_input_out).unwrap()),
        non_header_lines(&std::fs::read_to_string(&baseline).unwrap())
    );

    let without_md_nm = |text: &str| -> Vec<String> {
        non_header_lines(text)
            .into_iter()
            .map(|line| {
                line.split('\t')
                    .filter(|field| !field.starts_with("MD:Z:") && !field.starts_with("NM:i:"))
                    .collect::<Vec<_>>()
                    .join("\t")
            })
            .collect()
    };
    assert_eq!(
        exit_to_u8(ampliconclip::main(&argv(
            "ampliconclip",
            &[
                "--no-PG",
                "--keep-tag",
                "--output-fmt=cram",
                "-T",
                reference.to_str().unwrap(),
                "-b",
                bed.to_str().unwrap(),
                "-o",
                out_cram.to_str().unwrap(),
                input_sam.to_str().unwrap(),
            ]
        ))),
        0
    );
    let cram_text = htslib_rs::alignment_compat::view_cram_as_sam_text_with_reference(
        std::io::Cursor::new(std::fs::read(&out_cram).unwrap()),
        &reference,
        None,
    )
    .unwrap();
    assert_eq!(
        without_md_nm(&cram_text),
        without_md_nm(&std::fs::read_to_string(&baseline).unwrap())
    );

    set_current_global_args(SamGlobalArgs::default());
}

#[test]
fn markdup_matches_upstream_test_markdup_fixtures() {
    use samtools_rs::commands::markdup;
    // Byte-exact vs upstream `samtools/test/markdup` expected SAMs
    // (modulo @PG, suppressed by --no-PG): default template mode,
    // `-r` removal, `-S` supplementary propagation, and `--mode s`
    // sequence mode with optical-distance + barcode-tag keying.
    let d = fixtures_dir();
    let tmp = tmp_dir("markdup-fixtures");
    let cases: &[(&[&str], &str, &str)] = &[
        (&[], "5_markdup", "5_markdup.expected.sam"),
        (&["-r"], "6_remove_dups", "6_remove_dups.expected.sam"),
        (&["-S"], "7_mark_supp_dup", "7_mark_supp_dup.expected.sam"),
        (
            &[
                "-S",
                "-d",
                "100",
                "--mode",
                "s",
                "-t",
                "--barcode-tag",
                "BX",
            ],
            "13_optical_barcode_tag",
            "13_optical_barcode_tag.expected.sam",
        ),
        (
            &["-S", "-d", "100", "--mode", "s", "-t"],
            "8_optical_dup",
            "8_optical_dup.expected.sam",
        ),
        (
            &["-S", "-d", "2500", "--mode", "s", "-t", "--include-fails"],
            "9_optical_dup_qcfail",
            "9_optical_dup_qcfail.expected.sam",
        ),
        (
            &["-S", "-d", "2500", "--mode", "s", "-t", "-S"],
            "10_optical_chain",
            "10_optical_chain.expected.sam",
        ),
        (
            &[
                "--mode",
                "t",
                "-t",
                "--duplicate-count",
                "--barcode-tag",
                "BC",
                "-S",
            ],
            "18_primary_duplicate_count",
            "18_primary_duplicate_count.expected.sam",
        ),
        (
            &["-d", "100", "--mode", "s", "-t", "--use-read-groups"],
            "17_read_group",
            "17_read_group.expected.sam",
        ),
        (
            &[
                "-S",
                "-d",
                "2500",
                "--mode",
                "s",
                "-t",
                "--read-coords",
                "([[:digit:]]+):([[:digit:]]+)$",
                "--coords-order",
                "xy",
            ],
            "12_optical_chain_regex",
            "12_optical_chain_regex.expected.sam",
        ),
        (
            &["-S", "-d", "100", "--mode", "s", "-t", "--barcode-name"],
            "14_optical_barcode_name",
            "14_optical_barcode_name.expected.sam",
        ),
        (
            &[
                "-S",
                "-d",
                "100",
                "--mode",
                "s",
                "-t",
                "--barcode-rgx",
                "^([!-9;-?A-~]+):[0-9]+:",
                "--read-coords",
                "^[!-9;-?A-~]+:([0-9]+):([0-9]+)",
                "--coords-order",
                "xy",
            ],
            "15_optical_barcode_rgx_name",
            "15_optical_barcode_rgx_name.expected.sam",
        ),
        (
            &[
                "-S",
                "-d",
                "100",
                "--mode",
                "s",
                "-t",
                "--read-coords",
                "^([0-9]+):([0-9]+):([[:print:]]+)",
                "--coords-order",
                "xyt",
            ],
            "11_optical_dup_regex",
            "11_optical_dup_regex.expected.sam",
        ),
        (
            &[
                "-S",
                "-d",
                "100",
                "--mode",
                "s",
                "-t",
                "--barcode-rgx",
                "^([!-9;-?A-~]+):[0-9]+:",
                "--read-coords",
                "^[!-9;-?A-~]+:([0-9]{4})([0-9]{4})",
                "--coords-order",
                "xy",
            ],
            "16_optical_barcode_rgx_name_test_2",
            "16_optical_barcode_rgx_name_test_2.expected.sam",
        ),
    ];
    for (flags, stem, expected) in cases {
        let inp = d.join("markdup").join(format!("{stem}.sam"));
        let out = tmp.join(format!("{stem}.out.sam"));
        let mut rest: Vec<&str> = flags.to_vec();
        rest.extend(["-O", "sam", "--no-PG"]);
        let inp_s = inp.to_str().unwrap();
        let out_s = out.to_str().unwrap();
        rest.extend(["-o", out_s, inp_s]);
        assert_eq!(
            exit_to_u8(markdup::main(&argv("markdup", &rest))),
            0,
            "markdup {stem}"
        );
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            std::fs::read_to_string(d.join("markdup").join(expected)).unwrap(),
            "markdup {stem} byte-exact vs upstream"
        );
    }
}

#[test]
fn markdup_upstream_expect_fail_cases_return_one_with_expected_partial_output() {
    use samtools_rs::commands::markdup;
    let d = fixtures_dir();
    let tmp = tmp_dir("markdup-expect-fail");
    let cases = ["1_name_sort", "2_bad_order", "3_missing_mc", "4_missing_ms"];

    for stem in cases {
        let inp = d.join("markdup").join(format!("{stem}.sam"));
        let out = tmp.join(format!("{stem}.out.sam"));
        let inp_s = inp.to_str().unwrap();
        let out_s = out.to_str().unwrap();
        assert_eq!(
            exit_to_u8(markdup::main(&argv(
                "markdup",
                &["-O", "sam", "--no-PG", "-o", out_s, inp_s]
            ))),
            1,
            "markdup {stem} should fail like upstream"
        );
        let actual = std::fs::read_to_string(&out).unwrap_or_default();
        let expected =
            std::fs::read_to_string(d.join("markdup").join(format!("{stem}.expected.sam")))
                .unwrap();
        assert_eq!(actual, expected, "markdup {stem} partial stdout");
    }
}

#[test]
fn markdup_paired_end_groups_pairs_and_flags_duplicates() {
    use samtools_rs::commands::markdup;
    let tmp = tmp_dir("markdup-pe");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    // Two FR pairs at the same coordinates (positions 1 and 91). pair_a has
    // higher combined MAPQ; pair_b should be flagged as duplicate. A third
    // pair at different positions is unique.
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:1000\n",
            "pair_a\t99\tchr1\t1\t60\t10M\t=\t91\t100\tACGTACGTAC\tIIIIIIIIII\tMC:Z:10M\tms:i:400\n",
            "pair_b\t99\tchr1\t1\t10\t10M\t=\t91\t100\tACGTACGTAC\tIIIIIIIIII\tMC:Z:10M\tms:i:250\n",
            "pair_a\t147\tchr1\t91\t60\t10M\t=\t1\t-100\tACGTACGTAC\tIIIIIIIIII\tMC:Z:10M\tms:i:400\n",
            "pair_b\t147\tchr1\t91\t10\t10M\t=\t1\t-100\tACGTACGTAC\tIIIIIIIIII\tMC:Z:10M\tms:i:250\n",
            "pair_c\t99\tchr1\t200\t60\t10M\t=\t291\t100\tACGTACGTAC\tIIIIIIIIII\tMC:Z:10M\tms:i:400\n",
            "pair_c\t147\tchr1\t291\t60\t10M\t=\t200\t-100\tACGTACGTAC\tIIIIIIIIII\tMC:Z:10M\tms:i:400\n",
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(markdup::main(&argv(
            "markdup",
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

    let text = std::fs::read_to_string(&out).unwrap();
    let flag_of = |name: &str, flag_value: u32| -> u32 {
        let line = text
            .lines()
            .find(|l| {
                l.starts_with(name) && {
                    let f: u32 = l.split('\t').nth(1).unwrap().parse().unwrap();
                    f == flag_value
                }
            })
            .expect("record present");
        line.split('\t').nth(1).unwrap().parse().unwrap()
    };

    // pair_a: NOT a duplicate (highest combined MAPQ)
    assert_eq!(flag_of("pair_a\t", 99) & 0x400, 0);
    assert_eq!(flag_of("pair_a\t", 147) & 0x400, 0);
    // pair_b: BOTH reads marked as duplicates
    assert_eq!(flag_of("pair_b\t", 99 | 0x400) & 0x400, 0x400);
    assert_eq!(flag_of("pair_b\t", 147 | 0x400) & 0x400, 0x400);
    // pair_c: NOT a duplicate (different coordinates)
    assert_eq!(flag_of("pair_c\t", 99) & 0x400, 0);
    assert_eq!(flag_of("pair_c\t", 147) & 0x400, 0);
}

#[test]
fn mpileup_minus_b_ff_matches_upstream_out3() {
    let tmp = tmp_dir("mpileup-out3");
    let output = tmp.join("out.3");
    let input = fixtures_dir().join("dat").join("mpileup.1.sam");
    let reference = fixtures_dir().join("dat").join("mpileup.ref.fa");
    let expected = fixtures_dir().join("dat").join("mpileup.out.3");

    assert_eq!(
        exit_to_u8(mpileup::main(&argv(
            "mpileup",
            &[
                "-B",
                "--ff",
                "0x14",
                "-f",
                reference.to_str().unwrap(),
                "-r17:1050-1060",
                "-o",
                output.to_str().unwrap(),
                input.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(output).unwrap(),
        std::fs::read_to_string(expected).unwrap()
    );
}

#[test]
fn mpileup_default_baq_matches_upstream_out1() {
    let tmp = tmp_dir("mpileup-out1-baq");
    let list = tmp.join("inputs.list");
    let output = tmp.join("out.1");
    let dat = fixtures_dir().join("dat");
    let reference = dat.join("mpileup.ref.fa");
    let expected = dat.join("mpileup.out.1");
    std::fs::write(
        &list,
        format!(
            "{}\n{}\n{}\n",
            dat.join("mpileup.1.sam").display(),
            dat.join("mpileup.2.sam").display(),
            dat.join("mpileup.3.sam").display()
        ),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(mpileup::main(&argv(
            "mpileup",
            &[
                "-b",
                list.to_str().unwrap(),
                "-f",
                reference.to_str().unwrap(),
                "-r17:100-150",
                "-o",
                output.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        std::fs::read_to_string(output).unwrap(),
        std::fs::read_to_string(expected).unwrap()
    );
}

#[test]
fn mpileup_overlap_removal_matches_upstream_out5() {
    let tmp = tmp_dir("mpileup-out5");
    let output = tmp.join("out.5");
    let input = fixtures_dir().join("mpileup").join("overlap.bam");
    let expected = fixtures_dir().join("dat").join("mpileup.out.5");

    assert_eq!(
        exit_to_u8(mpileup::main(&argv(
            "mpileup",
            &[input.to_str().unwrap(), "-o", output.to_str().unwrap()]
        ))),
        0
    );

    let line = std::fs::read_to_string(output)
        .unwrap()
        .lines()
        .find(|l| l.contains("128814202"))
        .unwrap()
        .to_string();
    assert_eq!(line + "\n", std::fs::read_to_string(expected).unwrap());
}

#[test]
fn index_bam_without_so_coordinate_header() {
    // Upstream `samtools index` indexes coordinate-ordered BAMs whose
    // header omits `@HD SO:coordinate` (completed library batch #6).
    let tmp = tmp_dir("index-no-so");
    let bam = tmp.join("test_input_1_a.bam");
    std::fs::copy(fixtures_dir().join("dat").join("test_input_1_a.bam"), &bam).unwrap();

    assert_eq!(
        exit_to_u8(index::main(&argv("index", &[bam.to_str().unwrap()]))),
        0
    );
    assert!(bam.with_extension("bam.bai").exists());
}

#[test]
fn index_threads_use_bgzf_worker_reader() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    set_current_global_args(SamGlobalArgs::default());
    let tmp = tmp_dir("index-threads");
    let serial_bam = tmp.join("serial.bam");
    let local_threaded_bam = tmp.join("local-threaded.bam");
    let global_threaded_bam = tmp.join("global-threaded.bam");
    let fixture = fixtures_dir().join("dat").join("test_input_1_a.bam");
    for bam in [&serial_bam, &local_threaded_bam, &global_threaded_bam] {
        std::fs::copy(&fixture, bam).unwrap();
    }

    assert_eq!(
        exit_to_u8(index::main(&argv("index", &[serial_bam.to_str().unwrap()]))),
        0
    );
    assert_eq!(
        exit_to_u8(index::main(&argv(
            "index",
            &["--threads=2", local_threaded_bam.to_str().unwrap()]
        ))),
        0
    );

    set_current_global_args(SamGlobalArgs {
        threads: Some(2),
        ..SamGlobalArgs::default()
    });
    assert_eq!(
        exit_to_u8(index::main(&argv(
            "index",
            &[global_threaded_bam.to_str().unwrap()]
        ))),
        0
    );
    set_current_global_args(SamGlobalArgs::default());

    let serial_index = std::fs::read(serial_bam.with_extension("bam.bai")).unwrap();
    assert_eq!(
        serial_index,
        std::fs::read(local_threaded_bam.with_extension("bam.bai")).unwrap()
    );
    assert_eq!(
        serial_index,
        std::fs::read(global_threaded_bam.with_extension("bam.bai")).unwrap()
    );
}

#[test]
fn view_dot_region_means_whole_file() {
    // HTSlib region grammar: `.` == everything (completed library batch #10).
    let bam = fixtures_dir().join("dat").join("test_input_1_a.bam");
    let tmp = tmp_dir("view-dot");
    let with_dot = tmp.join("dot.txt");
    let no_region = tmp.join("all.txt");

    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &["-c", bam.to_str().unwrap(), "."]
        ))),
        0
    );
    // -c writes the count to stdout; compare via -o.
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &[
                "-c",
                "-o",
                with_dot.to_str().unwrap(),
                bam.to_str().unwrap(),
                "."
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &[
                "-c",
                "-o",
                no_region.to_str().unwrap(),
                bam.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(with_dot).unwrap(),
        std::fs::read_to_string(no_region).unwrap()
    );
}

#[test]
fn view_large_chrom_csi_region_matches_upstream() {
    // completed library batch #12: ref2 is 541556283 bp (> 2^29). With a header-aware
    // CSI, `view large_chrom.bam ref2` and `ref2:1-541556283` must both
    // produce dat/large_chrom.out without panicking.
    let tmp = tmp_dir("large-chrom");
    let bam = tmp.join("large_chrom.bam");
    let sam = fixtures_dir().join("dat").join("large_chrom.sam");
    let expected =
        std::fs::read_to_string(fixtures_dir().join("dat").join("large_chrom.out")).unwrap();

    // SAM -> BAM
    let bam_out = tmp.join("lc.tmp");
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &["-b", "-o", bam.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );
    // CSI index
    assert_eq!(
        exit_to_u8(index::main(&argv("index", &["-c", bam.to_str().unwrap()]))),
        0
    );

    for region in ["ref2", "ref2:1-541556283"] {
        assert_eq!(
            exit_to_u8(view::main(&argv(
                "view",
                &[
                    "-o",
                    bam_out.to_str().unwrap(),
                    bam.to_str().unwrap(),
                    region,
                ]
            ))),
            0,
            "region {region}"
        );
        assert_eq!(
            std::fs::read_to_string(&bam_out).unwrap(),
            expected,
            "region {region}"
        );
    }
}

#[test]
fn threads_flag_is_byte_identical_for_view_and_sort() {
    // completed library batch #8: `-@ N` must not change output bytes (worker-pool
    // wiring is perf-only). `--no-PG` isolates payload from the @PG CL
    // string, which legitimately embeds the thread arg.
    let tmp = tmp_dir("threads");
    let sam = fixtures_dir().join("dat").join("mpileup.1.sam");

    let v4 = tmp.join("v4.bam");
    let v1 = tmp.join("v1.bam");
    for (out, n) in [(&v4, "-@4"), (&v1, "-@1")] {
        assert_eq!(
            exit_to_u8(view::main(&argv(
                "view",
                &[
                    "--no-PG",
                    n,
                    "-b",
                    "-o",
                    out.to_str().unwrap(),
                    sam.to_str().unwrap()
                ]
            ))),
            0
        );
    }
    assert_eq!(
        std::fs::read(&v4).unwrap(),
        std::fs::read(&v1).unwrap(),
        "view -@ output must be byte-identical"
    );

    let s4 = tmp.join("s4.bam");
    let s1 = tmp.join("s1.bam");
    for (out, n) in [(&s4, "-@4"), (&s1, "-@1")] {
        assert_eq!(
            exit_to_u8(sort::main(&argv(
                "sort",
                &[
                    "--no-PG",
                    n,
                    "-o",
                    out.to_str().unwrap(),
                    v1.to_str().unwrap()
                ]
            ))),
            0
        );
    }
    assert_eq!(
        std::fs::read(&s4).unwrap(),
        std::fs::read(&s1).unwrap(),
        "sort -@ output must be byte-identical"
    );
}

#[test]
fn threads_flag_is_byte_identical_for_stats() {
    // Extends the worker-pool perf-only invariant to the reader-side
    // `stats` text path: `-@ N` must not change output bytes.
    let tmp = tmp_dir("threads-stats");
    let sam = fixtures_dir().join("dat").join("mpileup.1.sam");

    // Stage a BAM so the BGZF reader worker pool is exercised.
    let bam = tmp.join("in.bam");
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &[
                "--no-PG",
                "-b",
                "-o",
                bam.to_str().unwrap(),
                sam.to_str().unwrap(),
            ],
        ))),
        0
    );
    let bam = bam.to_str().unwrap();

    let t1 = tmp.join("stats.t1");
    let t4 = tmp.join("stats.t4");
    for (out, n) in [(&t1, "1"), (&t4, "4")] {
        assert_eq!(
            exit_to_u8(samtools_run(argv(
                "samtools",
                &["stats", "-@", n, "-o", out.to_str().unwrap(), bam],
            ))),
            0,
            "stats -@ {n}"
        );
    }
    // `#` comment lines legitimately embed the command line (which
    // includes `-@ N` and the distinct `-o` path) and the version
    // banner; the perf-only invariant is that the data body is
    // identical regardless of thread count.
    let body = |p: &std::path::Path| -> String {
        std::fs::read_to_string(p)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(
        body(&t1),
        body(&t4),
        "stats -@ data body must be identical across thread counts"
    );
}

#[test]
fn stats_cram_without_region_matches_bam_seq_quality() {
    // completed library batch #2: no-region CRAM stats now use the full-record iterator,
    // so sequence-length/quality/length SN lines match the BAM equivalent
    // (NM-derived `mismatches`/`error rate` excepted — CRAM has no NM).
    use samtools_rs::commands::stats;
    let bam = htslib_fixtures_dir().join("range.bam");
    let cram = htslib_fixtures_dir().join("range.cram");
    let reference = htslib_fixtures_dir().join("ce.fa");
    let tmp = tmp_dir("stats-cram");
    let bo = tmp.join("bam.txt");
    let co = tmp.join("cram.txt");

    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[bam.to_str().unwrap(), "-o", bo.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(stats::main(&argv(
            "stats",
            &[
                "-r",
                reference.to_str().unwrap(),
                cram.to_str().unwrap(),
                "-o",
                co.to_str().unwrap()
            ]
        ))),
        0
    );

    let pick = |s: &str| -> Vec<String> {
        s.lines()
            .filter(|l| {
                l.starts_with("SN\t") && !l.contains("mismatches:") && !l.contains("error rate:")
            })
            .map(str::to_string)
            .collect()
    };
    let b = std::fs::read_to_string(&bo).unwrap();
    let c = std::fs::read_to_string(&co).unwrap();
    assert_eq!(pick(&b), pick(&c));
    // Sanity: sequence/quality actually populated (not the old zeros).
    assert!(b.contains("SN\ttotal length:\t11200"));
    assert!(c.contains("SN\ttotal length:\t11200"));
}

#[test]
fn checksum_cram_matches_bam_via_all_record_iterator() {
    // completed library batch #2: whole-CRAM checksum via the htslib-rs all-record
    // iterator must equal the BAM checksum (checksum is order-agnostic).
    let _g = GLOBAL_ARGS_LOCK.lock().unwrap();
    let bam = htslib_fixtures_dir().join("range.bam");
    let cram = htslib_fixtures_dir().join("range.cram");
    let reference = htslib_fixtures_dir().join("ce.fa");
    let tmp = tmp_dir("checksum-cram");
    let bo = tmp.join("b.txt");
    let co = tmp.join("c.txt");

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "checksum",
                bam.to_str().unwrap(),
                "-o",
                bo.to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "--reference",
                reference.to_str().unwrap(),
                "checksum",
                cram.to_str().unwrap(),
                "-o",
                co.to_str().unwrap(),
            ]
        ))),
        0
    );

    assert_eq!(
        checksum_all_row(&std::fs::read_to_string(&bo).unwrap()),
        checksum_all_row(&std::fs::read_to_string(&co).unwrap()),
    );
}

#[test]
fn view_write_index_matches_post_pass_index() {
    // `view --write-index` BAM output follows upstream auto-indexing and writes CSI.
    let bam = htslib_fixtures_dir().join("range.bam");
    let tmp = tmp_dir("view-write-index");
    let auto = tmp.join("auto.bam");
    let post_csi = tmp.join("post.bam.csi");

    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &[
                "--write-index",
                "-b",
                "-o",
                auto.to_str().unwrap(),
                bam.to_str().unwrap()
            ]
        ))),
        0
    );
    let auto_csi = auto.with_extension("bam.csi");
    assert!(auto_csi.exists(), "view --write-index did not write a .csi");
    let auto_bytes = std::fs::read(&auto_csi).unwrap();

    assert_eq!(
        exit_to_u8(index::main(&argv(
            "index",
            &["-c", auto.to_str().unwrap(), post_csi.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        auto_bytes,
        std::fs::read(&post_csi).unwrap(),
        "--write-index CSI must equal the post-pass CSI"
    );
}

#[test]
fn view_star_region_selects_unplaced_reads() {
    // HTSlib region grammar: `*` selects only unplaced (RNAME `*`) reads
    // (completed library batch #10). Verified for SAM-text + count output.
    let tmp = tmp_dir("view-star");
    let sam = tmp.join("un.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:coordinate\n\
         @SQ\tSN:ref1\tLN:100\n\
         m1\t0\tref1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\
         u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\n\
         u2\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTT\tIIIII\n",
    )
    .unwrap();
    let bam = tmp.join("un.bam");
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &["-b", "-o", bam.to_str().unwrap(), sam.to_str().unwrap()]
        ))),
        0
    );
    assert_eq!(
        exit_to_u8(index::main(&argv("index", &[bam.to_str().unwrap()]))),
        0
    );

    let out = tmp.join("star.sam");
    assert_eq!(
        exit_to_u8(view::main(&argv(
            "view",
            &["-o", out.to_str().unwrap(), bam.to_str().unwrap(), "*"]
        ))),
        0
    );
    let body: Vec<String> = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(str::to_string)
        .collect();
    assert_eq!(body.len(), 2);
    assert!(body.iter().all(|l| l.split('\t').nth(2) == Some("*")));
}

#[test]
fn consensus_simple_matches_upstream_fixtures() {
    // completed library batch #1 wiring: `samtools consensus --mode simple` on the
    // htslib-rs pileup engine, byte-exact vs test/consensus/expected/*.
    let dir = fixtures_dir()
        .parent()
        .unwrap()
        .join("test")
        .join("consensus");
    let sam = dir.join("consen1.sam");
    let tmp = tmp_dir("consensus");

    let cases: &[(&[&str], &str)] = &[
        (&["-m", "simple", "-c", "0.6"], "1.out"),
        (&["-m", "simple", "-c", "0.6", "--show-del", "yes"], "2.out"),
        (&["-m", "simple", "-c", "0.6", "--show-ins", "no"], "3.out"),
        (
            &[
                "-m",
                "simple",
                "-c",
                "0.6",
                "--show-del",
                "yes",
                "--show-ins",
                "no",
            ],
            "4.out",
        ),
        (&["-f", "fastq", "-m", "simple", "-c", "0.6"], "1q.out"),
        (
            &[
                "-f",
                "fastq",
                "-m",
                "simple",
                "-c",
                "0.6",
                "--show-del",
                "yes",
            ],
            "2q.out",
        ),
        (
            &[
                "-f",
                "fastq",
                "-m",
                "simple",
                "-c",
                "0.6",
                "--show-ins",
                "no",
            ],
            "3q.out",
        ),
        (&["-f", "pileup", "-m", "simple", "-c", "0.6"], "1p.out"),
        (
            &["-f", "fastq", "-m", "simple", "--call-fract", "0.600"],
            "1q.out",
        ),
        (
            &["-f", "fastq", "-m", "simple", "--call-fract", "0.601"],
            "5q.out",
        ),
        (
            &["-f", "pileup", "-m", "simple", "--call-fract", "0.601"],
            "5p.out",
        ),
    ];

    for (i, (args, expected)) in cases.iter().enumerate() {
        let out = tmp.join(format!("c{i}"));
        let mut a: Vec<&str> = args.to_vec();
        a.push("-o");
        a.push(out.to_str().unwrap());
        a.push(sam.to_str().unwrap());
        assert_eq!(
            exit_to_u8(consensus::main(&argv("consensus", &a))),
            0,
            "case {expected}"
        );
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            std::fs::read_to_string(dir.join("expected").join(expected)).unwrap(),
            "case {expected} args={args:?}"
        );
    }
}

#[test]
fn coverage_matches_upstream_tabular_fixtures() {
    // completed library batch #1: exact coverage — `%g`/`%.3g` formatting, min_depth
    // gating of meandepth/meanbaseq, and pileup-arrival row ordering.
    use samtools_rs::commands::coverage;
    let d = fixtures_dir();
    let sample = d.join("dat").join("sample.sam");
    let tmp = tmp_dir("coverage-fix");
    let sample1 = tmp.join("sample1.sam");
    let text = std::fs::read_to_string(&sample).unwrap();
    let filtered: String = text
        .lines()
        .filter(|l| !l.contains("A1"))
        .map(|l| format!("{l}\n"))
        .collect();
    std::fs::write(&sample1, filtered).unwrap();

    let cases: &[(&[&str], &str)] = &[
        (&[], "1.expected"),
        (&["--min-depth", "1"], "1.expected"),
        (&["--min-depth", "2"], "2.expected"),
        (&["--min-depth", "2", "-Q", "8", "-q", "45"], "3.expected"),
    ];
    for (args, exp) in cases {
        let out = tmp.join(exp);
        let mut a: Vec<&str> = args.to_vec();
        a.push("-o");
        a.push(out.to_str().unwrap());
        a.push(sample.to_str().unwrap());
        assert_eq!(
            exit_to_u8(coverage::main(&argv("coverage", &a))),
            0,
            "{exp}"
        );
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            std::fs::read_to_string(d.join("coverage").join(exp)).unwrap(),
            "coverage {exp} args={args:?}"
        );
    }

    // Multi-input (sample.sam + sample1.sam).
    for (md, exp) in [("1", "4.expected"), ("4", "5.expected")] {
        let out = tmp.join(exp);
        assert_eq!(
            exit_to_u8(coverage::main(&argv(
                "coverage",
                &[
                    "--min-depth",
                    md,
                    "-o",
                    out.to_str().unwrap(),
                    sample.to_str().unwrap(),
                    sample1.to_str().unwrap(),
                ]
            ))),
            0
        );
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            std::fs::read_to_string(d.join("coverage").join(exp)).unwrap(),
            "coverage {exp}"
        );
    }
}

#[test]
fn depth_large_pos_matches_upstream() {
    // completed library batch #1: sparse depth (no OOM on LN:10001009800) + whitespace
    // BED parsing — upstream large_pos depth fixtures byte-exact.
    use samtools_rs::commands::depth;
    let d = fixtures_dir().join("large_pos");
    let tmp = tmp_dir("depth-largepos");

    let o1 = tmp.join("depth.out");
    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-o",
                o1.to_str().unwrap(),
                d.join("longref.sam").to_str().unwrap()
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&o1).unwrap(),
        std::fs::read_to_string(d.join("depth.expected.out")).unwrap()
    );

    let o2 = tmp.join("depth_bed.out");
    assert_eq!(
        exit_to_u8(depth::main(&argv(
            "depth",
            &[
                "-b",
                d.join("test.bed").to_str().unwrap(),
                "-o",
                o2.to_str().unwrap(),
                d.join("longref.sam").to_str().unwrap(),
            ]
        ))),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&o2).unwrap(),
        std::fs::read_to_string(d.join("depth_bed.expected.out")).unwrap()
    );
}

#[test]
fn fixmate_accepts_reference_backed_cram_input() {
    use samtools_rs::commands::fixmate;

    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    set_current_global_args(SamGlobalArgs::default());

    let tmp = tmp_dir("fixmate-cram");
    let reference = tmp.join("ref.fa");
    let input_sam = tmp.join("in.sam");
    let input_bam = tmp.join("in.bam");
    let input_cram = tmp.join("in.cram");
    let out = tmp.join("out.sam");
    let out_cram = tmp.join("out.cram");

    std::fs::write(&reference, ">chr1\nACGTTGCA\n").unwrap();
    samtools_rs::reference::ensure_fai_index(&reference, None).unwrap();
    std::fs::write(
        &input_sam,
        concat!(
            "@HD\tVN:1.6\tSO:queryname\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "r1\t65\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "r1\t129\tchr1\t5\t60\t4M\t*\t0\t0\tTGCA\t####\n",
        ),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &input_sam,
        std::fs::File::create(&input_bam).unwrap(),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_cram_from_bam_path_with_reference(
        &input_bam,
        &reference,
        std::fs::File::create(&input_cram).unwrap(),
    )
    .unwrap();

    set_current_global_args(SamGlobalArgs {
        reference: Some(reference.clone()),
        ..SamGlobalArgs::default()
    });
    assert_eq!(
        exit_to_u8(fixmate::main(&argv(
            "fixmate",
            &[
                "--no-PG",
                "-O",
                "sam",
                input_cram.to_str().unwrap(),
                out.to_str().unwrap(),
            ]
        ))),
        0
    );
    set_current_global_args(SamGlobalArgs::default());

    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.contains("r1\t65\tchr1\t1\t60\t4M\t=\t5\t4\tACGT\t!!!!"));
    assert!(text.contains("r1\t129\tchr1\t5\t60\t4M\t=\t1\t-4\tTGCA\t####"));

    set_current_global_args(SamGlobalArgs {
        reference: Some(reference.clone()),
        ..SamGlobalArgs::default()
    });
    assert_eq!(
        exit_to_u8(fixmate::main(&argv(
            "fixmate",
            &[
                "--no-PG",
                "--output-fmt=cram",
                input_sam.to_str().unwrap(),
                out_cram.to_str().unwrap(),
            ]
        ))),
        0
    );
    set_current_global_args(SamGlobalArgs::default());

    let records = htslib_rs::alignment_compat::summarize_cram_records_from_path_with_reference(
        &out_cram, &reference,
    )
    .unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].name_bytes(), Some(&b"r1"[..]));
    assert_eq!(records[0].mate_reference_sequence_id(), Some(0));
    assert_eq!(records[0].mate_alignment_start(), Some(5));
    assert_eq!(records[0].template_length(), 8);
    assert_eq!(records[1].mate_reference_sequence_id(), Some(0));
    assert_eq!(records[1].mate_alignment_start(), Some(1));
    assert_eq!(records[1].template_length(), -8);
}

#[test]
fn fixmate_matches_upstream_group() {
    // TODO.md §13 / completed library batch #7: full upstream test_fixmate group
    // byte-exact (modulo @PG, which the harness strips). Exercises
    // del-then-append aux ordering (MQ/MC/ct) and MC:Z:* semantics.
    use samtools_rs::commands::fixmate;
    let d = fixtures_dir().join("fixmate");
    let tmp = tmp_dir("fixmate-group");
    let cases: &[(&[&str], &str)] = &[
        (&["-z", "off", "-O", "sam"], "2_isize_overflow"),
        (&["-O", "sam"], "3_reverse_read_pp_lt"),
        (&["-O", "sam"], "4_reverse_read_pp_equal"),
        (&["-cO", "sam"], "5_ct"),
        (&["-cO", "sam"], "6_ct_replace"),
        (&["-z", "off", "-O", "sam"], "7_two_read_mapped"),
        (&["-z", "off", "-O", "sam"], "8_isize_overflow_64bit"),
        (&["-O", "sam"], "sanitize"),
        (&["-M", "-z", "off", "-O", "sam"], "mod_ok+"),
        (&["-M", "-z", "off", "-O", "sam"], "mod_ok-"),
        (&["-M", "-z", "off", "-O", "sam"], "mod_draft"),
        (&["-M", "-z", "off", "-O", "sam"], "mod_not_updated"),
        (&["-M", "-z", "off", "-O", "sam"], "mod_not_updated_noML"),
        (&["-M", "-z", "off", "-O", "sam"], "mod_not_updated_noMN"),
        (&["-M", "-z", "off", "-O", "sam"], "mod_bad_MN"),
        (&["-M", "-z", "off", "-O", "sam"], "mod_MN_only"),
        (&["-M", "-z", "off", "-O", "sam"], "mod_ML_only"),
        (&["-M", "-z", "off", "-O", "sam"], "mod_ML_wrong_len"),
        (&["-M", "-z", "off", "-O", "sam"], "mod_noseq"),
        (&["-M", "-z", "off", "-O", "sam"], "mod_bounds"),
    ];
    for (args, name) in cases {
        let out = tmp.join(format!("{name}.out"));
        let mut a: Vec<String> = args.iter().map(|s| s.to_string()).collect();
        a.push("--no-PG".into());
        a.push(d.join(format!("{name}.sam")).to_str().unwrap().into());
        a.push(out.to_str().unwrap().into());
        let refs: Vec<&str> = a.iter().map(String::as_str).collect();
        assert_eq!(
            exit_to_u8(fixmate::main(&argv("fixmate", &refs))),
            0,
            "{name}"
        );
        assert_eq!(
            without_pg_lines(&std::fs::read_to_string(&out).unwrap()),
            without_pg_lines(
                &std::fs::read_to_string(d.join(format!("{name}.sam.expected"))).unwrap()
            ),
            "fixmate {name}"
        );
    }
}

#[test]
fn addreplacerg_matches_upstream_group() {
    // TODO.md §13: full upstream test_addrprg group byte-exact (modulo
    // @PG). Covers -m overwrite_all/orphan_only, full `@RG\t..` -r spec
    // (escaped tabs), incremental -r ID:/-r CN:, -w edit, -R overwrite.
    let a = fixtures_dir().join("addrprg");
    let tmp = tmp_dir("addrprg-group");
    let cases: &[(&[&str], &str, &str)] = &[
        (
            &["-m", "overwrite_all"],
            "1_fixup.sam",
            "1_fixup.sam.expected",
        ),
        (
            &["-m", "orphan_only"],
            "2_fixup_orphan.sam",
            "2_fixup_orphan.sam.expected",
        ),
        (
            &["-r", "@RG\tID:1#8\tCN:SC"],
            "4_fixup_norg.sam",
            "4_fixup_norg.sam.expected",
        ),
        (
            &["-r", "ID:1#8", "-r", "CN:SC"],
            "4_fixup_norg.sam",
            "4_fixup_norg.sam.expected",
        ),
        (
            &[
                "-w",
                "-r",
                "@RG\tID:1#8\tCN:Sanger\tDS:Testing the editing code.",
            ],
            "1_fixup.sam",
            "5_editrg.sam.expected",
        ),
        (
            &["-m", "overwrite_all", "-R", "1#8"],
            "1_fixup.sam",
            "1_fixup.sam.expected",
        ),
    ];
    for (i, (args, input, expected)) in cases.iter().enumerate() {
        let out = tmp.join(format!("o{i}"));
        let mut v: Vec<String> = vec!["-O".into(), "sam".into()];
        v.extend(args.iter().map(|s| s.to_string()));
        v.push(a.join(input).to_str().unwrap().into());
        v.push("-o".into());
        v.push(out.to_str().unwrap().into());
        let refs: Vec<&str> = v.iter().map(String::as_str).collect();
        assert_eq!(
            exit_to_u8(addreplacerg::main(&argv("addreplacerg", &refs))),
            0,
            "{expected}"
        );
        assert_eq!(
            without_pg_lines(&std::fs::read_to_string(&out).unwrap()),
            without_pg_lines(&std::fs::read_to_string(a.join(expected)).unwrap()),
            "addrprg case {i} ({expected}) args={args:?}"
        );
    }
}

/// Drives the entire upstream `test/consensus/consensus.reg` harness
/// in-process: every `INIT` line builds its BAM via `view` (with the
/// `--write-index`), then every `P <name> ... consensus <args>` line is
/// run and its output compared byte-for-byte to
/// `test/consensus/expected/<name>`. Locks all 77 cases: simple +
/// bayesian/recall, fasta/fastq/pileup, -a/-aa, -r, -T/--ref-qual,
/// --min-MQ/--min-BQ, show-del/ins, and glued short options.
#[test]
fn consensus_matches_upstream_consensus_reg() {
    use samtools_rs::commands::{consensus, view};
    let dir = fixtures_dir().join("consensus");
    let tmp = tmp_dir("consensus-reg");
    // Stage the input SAM/FA(/fai) so relative names in the .reg resolve.
    for e in std::fs::read_dir(&dir).unwrap() {
        let p = e.unwrap().path();
        if matches!(
            p.extension().and_then(|s| s.to_str()),
            Some("sam") | Some("fa") | Some("fai")
        ) {
            std::fs::copy(&p, tmp.join(p.file_name().unwrap())).unwrap();
        }
    }
    // Rewrite a bare relative fixture token to its staged abs path.
    let abs = |t: &str| -> String {
        if t.ends_with(".bam") || t.ends_with(".sam") || t.ends_with(".fa") || t.ends_with(".fai") {
            tmp.join(t).to_str().unwrap().to_string()
        } else {
            t.to_string()
        }
    };
    let reg = std::fs::read_to_string(dir.join("consensus.reg")).unwrap();
    let mut n = 0usize;
    for line in reg.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("INIT ") {
            // INIT x $samtools view --write-index A.sam -o B.bam
            let toks: Vec<&str> = rest.split_whitespace().collect();
            let i = toks.iter().position(|&t| t == "view").unwrap();
            let args: Vec<String> = toks[i + 1..].iter().map(|t| abs(t)).collect();
            let a: Vec<&str> = args.iter().map(String::as_str).collect();
            assert_eq!(
                exit_to_u8(view::main(&argv("view", &a))),
                0,
                "INIT failed: {line}"
            );
        } else if let Some(rest) = line.strip_prefix("P ") {
            // P <name> $samtools consensus $PARAM <args...>
            // A few lines tee through a shell pipeline
            // (`-o cons.tmp; cat cons.tmp; rm cons.tmp`); keep only the
            // segment before the first `;` and drop its trailing
            // `-o/--output cons.tmp` so we can supply our own `-o`.
            let head = rest.split(';').next().unwrap();
            let toks: Vec<&str> = head.split_whitespace().collect();
            let name = toks[0];
            let ci = toks.iter().position(|&t| t == "consensus").unwrap();
            let out = tmp.join(format!("got.{name}"));
            let mut args: Vec<String> = toks[ci + 1..]
                .iter()
                .filter(|&&t| t != "$PARAM")
                .map(|t| abs(t))
                .collect();
            if matches!(args.last().map(String::as_str), Some(s) if s.ends_with("cons.tmp"))
                && matches!(
                    args.get(args.len().wrapping_sub(2)).map(String::as_str),
                    Some("-o") | Some("--output")
                )
            {
                args.truncate(args.len() - 2);
            }
            args.push("-o".into());
            args.push(out.to_str().unwrap().to_string());
            let a: Vec<&str> = args.iter().map(String::as_str).collect();
            assert_eq!(
                exit_to_u8(consensus::main(&argv("consensus", &a))),
                0,
                "consensus exit nonzero: {name} args={a:?}"
            );
            let got = std::fs::read(&out).unwrap();
            let exp = std::fs::read(dir.join("expected").join(name)).unwrap();
            assert_eq!(
                String::from_utf8_lossy(&got),
                String::from_utf8_lossy(&exp),
                "consensus.reg case {name} differs"
            );
            n += 1;
        }
    }
    assert_eq!(n, 77, "expected 77 consensus.reg P-cases, ran {n}");
}

/// `samtools cram-size` default and `-v` reports are byte-exact vs
/// the upstream `test/cram_size/cram_size.reg` fixtures (completed
/// library batch #3). The `-e` "Container encodings" mode is included.
#[test]
fn cram_size_matches_upstream_cram_size_reg() {
    let dir = fixtures_dir().join("cram_size");
    let cram = dir.join("mpileup.1.cram");
    let tmp = tmp_dir("cram-size");

    for (args, expected) in [
        (vec![], "normal.out"),
        (vec!["-v"], "verbose.out"),
        (vec!["-e"], "encodings.out"),
    ] {
        let out = tmp.join(expected);
        let mut a: Vec<&str> = vec!["cram-size"];
        a.extend_from_slice(&args);
        a.push(cram.to_str().unwrap());
        a.push("-o");
        a.push(out.to_str().unwrap());
        assert_eq!(
            exit_to_u8(cram_size::main(&argv("cram-size", &a[1..]))),
            0,
            "args={a:?}"
        );
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            std::fs::read_to_string(dir.join("expected").join(expected)).unwrap(),
            "cram-size {expected} must be byte-exact"
        );
    }
}

/// Full upstream `test_reference`: build an embed_ref CRAM with
/// `view -e EXPR -O cram,embed_ref=1 -T ref`, then all four
/// `samtools reference` invocations (MD path / `-e` embedded, with
/// and without `-r`) are byte-exact vs the upstream fixtures —
/// completed library batch #2 complete (embed_ref read+write + cram2ref).
#[test]
fn reference_embed_ref_full_test_reference_byte_exact() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let d = fixtures_dir();
    let tmp = tmp_dir("reference-embed");
    let cram = build_reference_embed_ref_cram(&tmp);
    let md_full = d.join("reference/mpileup.MD.fa.expected");
    let embed_full = d.join("reference/mpileup.embed.fa.expected");

    let cases: Vec<(&[&str], String)> = vec![
        (&[], std::fs::read_to_string(&md_full).unwrap()),
        (&["-e"], std::fs::read_to_string(&embed_full).unwrap()),
        (
            &["-r", "17:1000-1500"],
            fasta_region(&md_full, "17", 1000, 1500),
        ),
        (
            &["-r", "17:1000-1500", "-e"],
            fasta_region(&embed_full, "17", 1000, 1500),
        ),
    ];
    for (i, (extra, expected)) in cases.into_iter().enumerate() {
        let out = tmp.join(format!("reference_embed_{i}.fa"));
        let mut a: Vec<String> = vec!["samtools".into(), "reference".into(), "-q".into()];
        a.extend(extra.iter().map(|s| s.to_string()));
        a.push(cram.to_str().unwrap().into());
        a.push("-o".into());
        a.push(out.to_str().unwrap().into());
        assert_eq!(
            exit_to_u8(samtools_run(
                a.iter().map(std::ffi::OsString::from).collect()
            )),
            0,
            "args={a:?}"
        );
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            expected,
            "reference embed case {i} must be byte-exact"
        );
    }
}
