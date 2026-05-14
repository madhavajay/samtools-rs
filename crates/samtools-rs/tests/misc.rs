//! Smoke-tests for `cat`, `reheader`, `fastq`, `samples`, `idxstats`,
//! `flagstat`, `index`, `faidx`, `import`, `bedcov`, `rmdup`, `split`.

use std::ffi::OsString;
use std::io::{BufReader, Cursor};
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::Mutex;

use htslib_rs::bam;
use htslib_rs::sam;
use samtools_rs::commands::{
    addreplacerg, bedcov, calmd, cat, checksum, faidx, fastq, fixmate, flagstat, fqidx, idxstats,
    import, index, reference, reheader, reset, rmdup, samples, split, view,
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

#[test]
fn flagstat_succeeds() {
    let p = sample_bam();
    assert_eq!(
        exit_to_u8(flagstat::main(&argv("flagstat", &[p.to_str().unwrap()]))),
        0
    );
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

    let expected = ">chr1:3-7 length: 5\nGTATG\n";
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
        ">chr1:3-7 length: 5\nGTATG\n"
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
fn cat_sam_inputs_write_sam_output_with_single_header() {
    let tmp = tmp_dir("cat-sam");
    let sam_a = tmp.join("a.sam");
    let sam_b = tmp.join("b.sam");
    let out = tmp.join("cat.sam");
    let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
    std::fs::write(
        &sam_a,
        format!("{header}a1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n"),
    )
    .unwrap();
    std::fs::write(
        &sam_b,
        format!("{header}b1\t0\tchr1\t5\t60\t4M\t*\t0\t0\tTGCA\t####\n"),
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(cat::main(&argv(
            "cat",
            &[
                "--no-PG",
                sam_a.to_str().unwrap(),
                sam_b.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert_eq!(text.matches("@SQ\tSN:chr1").count(), 1);
    assert!(text.lines().any(|line| line.starts_with("a1\t")));
    assert!(text.lines().any(|line| line.starts_with("b1\t")));
    assert!(!text.contains("\tCL:cat "));
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
    assert_eq!(
        exit_to_u8(addreplacerg::main(&argv(
            "addreplacerg",
            &[
                "--no-PG",
                "-R",
                "new",
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
fn calmd_drop_baq_removes_bq_tags() {
    let tmp = tmp_dir("calmd-drop-bq");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "@SQ\tSN:chr1\tLN:12\n",
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tBQ:Z:abcd\tNM:i:0\n",
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
            ]
        ))),
        0
    );

    let text = std::fs::read_to_string(&out).unwrap();
    let record = text.lines().find(|line| line.starts_with("r1\t")).unwrap();
    assert!(!record.contains("\tBQ:Z:"));
    assert!(record.contains("\tNM:i:0"));
}

#[test]
fn markdup_sam_input_flags_duplicates_keeping_highest_mapq() {
    use samtools_rs::commands::markdup;
    let tmp = tmp_dir("markdup-sam");
    let sam = tmp.join("in.sam");
    let out = tmp.join("out.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "low\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
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
        "highest MAPQ in group keeps primary"
    );
    assert_eq!(flag_of(low) & 0x400, 0x400, "low MAPQ duplicate flagged");
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
            "aa_high\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\tBC:Z:AA\n",
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
            "keep\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTGCA\t####\n",
            "dup\t0\tchr1\t1\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "dup\t2048\tchr1\t20\t10\t4M\t*\t0\t0\tACGT\t!!!!\n",
            "keep\t2048\tchr1\t30\t10\t4M\t*\t0\t0\tTGCA\t####\n",
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
            "pair_a\t99\tchr1\t1\t60\t10M\t=\t91\t100\tACGTACGTAC\tIIIIIIIIII\n",
            "pair_a\t147\tchr1\t91\t60\t10M\t=\t1\t-100\tACGTACGTAC\tIIIIIIIIII\n",
            "pair_b\t99\tchr1\t1\t10\t10M\t=\t91\t100\tACGTACGTAC\tIIIIIIIIII\n",
            "pair_b\t147\tchr1\t91\t10\t10M\t=\t1\t-100\tACGTACGTAC\tIIIIIIIIII\n",
            "pair_c\t99\tchr1\t200\t60\t10M\t=\t291\t100\tACGTACGTAC\tIIIIIIIIII\n",
            "pair_c\t147\tchr1\t291\t60\t10M\t=\t200\t-100\tACGTACGTAC\tIIIIIIIIII\n",
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
