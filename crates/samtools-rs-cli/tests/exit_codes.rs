use std::path::Path;
use std::process::Command;

fn tmp_path(name: &str, extension: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!(
        "samtools-rs-cli-{name}-{}.{}",
        std::process::id(),
        extension
    ))
}

fn run_missing(command: &str, missing: &Path) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_samtools"))
        .arg(command)
        .arg(missing)
        .output()
        .unwrap_or_else(|e| panic!("run samtools {command}: {e}"))
}

fn run_args(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_samtools"))
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("run samtools {args:?}: {e}"))
}

fn hts_missing_line(path: &Path) -> String {
    format!(
        "[E::hts_open_format] Failed to open file \"{}\" : No such file or directory\n",
        path.display()
    )
}

#[test]
fn missing_input_stderr_matches_upstream_for_common_commands() {
    let missing = std::env::temp_dir().join(format!(
        "samtools-rs-cli-missing-input-{}.bam",
        std::process::id()
    ));
    let _ = std::fs::remove_file(&missing);

    let view = run_missing("view", &missing);
    assert_eq!(view.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&view.stderr),
        format!(
            "{}samtools view: failed to open \"{}\" for reading: No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let index = run_missing("index", &missing);
    assert_eq!(index.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&index.stderr),
        format!(
            "{}samtools index: failed to open \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let sort = run_missing("sort", &missing);
    assert_eq!(sort.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&sort.stderr),
        format!(
            "{}samtools sort: can't open \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let flagstat = run_missing("flagstat", &missing);
    assert_eq!(flagstat.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&flagstat.stderr),
        format!(
            "{}samtools flagstat: Cannot open input file \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let idxstats = run_missing("idxstats", &missing);
    assert_eq!(idxstats.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&idxstats.stderr),
        format!(
            "{}samtools idxstats: failed to open \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let head = run_missing("head", &missing);
    assert_eq!(head.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&head.stderr),
        format!(
            "{}samtools head: failed to open \"{}\" for reading: No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );
}

#[test]
fn dict_missing_input_stderr_matches_upstream() {
    let missing = std::env::temp_dir().join(format!(
        "samtools-rs-cli-missing-input-{}.fa",
        std::process::id()
    ));
    let _ = std::fs::remove_file(&missing);

    let output = run_missing("dict", &missing);
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        format!(
            "samtools dict: Cannot open {}: No such file or directory\n",
            missing.display()
        )
    );
}

#[test]
fn faidx_and_fqidx_missing_input_stderr_matches_upstream() {
    let missing = std::env::temp_dir().join(format!(
        "samtools-rs-cli-missing-input-{}.fa",
        std::process::id()
    ));
    let _ = std::fs::remove_file(&missing);
    let expected = format!(
        "[E::fai_build3_core] Failed to open the file {} : No such file or directory\n\
         [faidx] Could not build fai index {}.fai\n",
        missing.display(),
        missing.display()
    );

    for command in ["faidx", "fqidx"] {
        let output = run_missing(command, &missing);
        assert_eq!(output.status.code(), Some(1), "{command}");
        assert!(output.stdout.is_empty(), "{command}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            expected,
            "{command}"
        );
    }
}

#[test]
fn more_missing_input_stderr_matches_upstream() {
    let missing = std::env::temp_dir().join(format!(
        "samtools-rs-cli-missing-input-{}.bam",
        std::process::id()
    ));
    let _ = std::fs::remove_file(&missing);

    let checksum = run_missing("checksum", &missing);
    assert_eq!(checksum.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&checksum.stderr),
        format!(
            "{}samtools checksum: Cannot open input file \"{}\": No such file or directory\n\
             [checksum] Failed to process data\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let coverage = run_missing("coverage", &missing);
    assert_eq!(coverage.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&coverage.stderr),
        format!(
            "{}samtools coverage: Could not open \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let depth = run_missing("depth", &missing);
    assert_eq!(depth.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&depth.stderr),
        format!(
            "{}samtools depth: Cannot open input file \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let samples = run_missing("samples", &missing);
    assert_eq!(samples.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&samples.stderr),
        format!(
            "{}samtools samples: Failed to open \"{}\" for reading: No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let missing_arg = missing.to_string_lossy();
    let addreplacerg = run_args(&[
        "addreplacerg",
        "-r",
        r"@RG\tID:foo\tSM:bar",
        missing_arg.as_ref(),
    ]);
    assert_eq!(addreplacerg.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&addreplacerg.stderr),
        format!(
            "{}samtools addreplacerg: could not open \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let reset = run_missing("reset", &missing);
    assert_eq!(reset.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&reset.stderr),
        format!(
            "{}Could not open {}\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let consensus = run_missing("consensus", &missing);
    assert_eq!(consensus.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&consensus.stderr),
        format!(
            "{}samtools consensus: Cannot open input file \"{}\": No such file or directory\n\
             samtools consensus: failed\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let bed = tmp_path("missing-input-regions", "bed");
    std::fs::write(&bed, "").unwrap();
    let bedcov = Command::new(env!("CARGO_BIN_EXE_samtools"))
        .arg("bedcov")
        .arg(&bed)
        .arg(&missing)
        .output()
        .expect("run samtools bedcov");
    assert_eq!(bedcov.status.code(), Some(2));
    assert_eq!(
        String::from_utf8_lossy(&bedcov.stderr),
        format!(
            "{}ERROR: fail to open index BAM file '{}'\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );
    let _ = std::fs::remove_file(bed);

    let amplicon_bed = tmp_path("missing-input-amplicons", "bed");
    let amplicon_out = tmp_path("missing-input-amplicons-out", "bam");
    std::fs::write(&amplicon_bed, "chr1\t0\t10\tamp1\t0\t+\n").unwrap();
    let amplicon_bed_arg = amplicon_bed.to_string_lossy();
    let amplicon_out_arg = amplicon_out.to_string_lossy();
    let ampliconclip = run_args(&[
        "ampliconclip",
        "-b",
        amplicon_bed_arg.as_ref(),
        missing_arg.as_ref(),
        "-o",
        amplicon_out_arg.as_ref(),
    ]);
    assert_eq!(ampliconclip.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&ampliconclip.stderr),
        format!(
            "{}samtools ampliconclip: cannot open input file: No such file or directory\n",
            hts_missing_line(&missing)
        )
    );

    let ampliconstats = run_args(&[
        "ampliconstats",
        amplicon_bed_arg.as_ref(),
        missing_arg.as_ref(),
    ]);
    assert_eq!(ampliconstats.status.code(), Some(255));
    assert_eq!(
        String::from_utf8_lossy(&ampliconstats.stderr),
        format!(
            "{}samtools ampliconstats: Cannot open input file \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let missing_arg = missing.to_string_lossy();
    let collate = run_args(&["collate", "-O", missing_arg.as_ref()]);
    assert_eq!(collate.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&collate.stderr),
        format!(
            "{}samtools collate: Cannot open input file \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let cat = run_missing("cat", &missing);
    assert_eq!(cat.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&cat.stderr),
        format!(
            "{}samtools cat: failed to open file '{}': No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let split = run_missing("split", &missing);
    assert_eq!(split.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&split.stderr),
        format!(
            "{}samtools split: Could not open \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let stats = run_missing("stats", &missing);
    assert_eq!(stats.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&stats.stderr),
        format!(
            "{}samtools stats: failed to open \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let out = tmp_path("missing-input-out", "bam");
    let out_arg = out.to_string_lossy();

    let fixmate = run_args(&["fixmate", missing_arg.as_ref(), out_arg.as_ref()]);
    assert_eq!(fixmate.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&fixmate.stderr),
        format!(
            "{}samtools fixmate: cannot open input file: No such file or directory\n",
            hts_missing_line(&missing)
        )
    );

    let markdup = run_args(&["markdup", missing_arg.as_ref(), out_arg.as_ref()]);
    assert_eq!(markdup.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&markdup.stderr),
        format!(
            "{}samtools markdup: error, failed to open \"{}\" for input: No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let rmdup = run_args(&["rmdup", missing_arg.as_ref(), out_arg.as_ref()]);
    assert_eq!(rmdup.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&rmdup.stderr),
        format!(
            "{}samtools rmdup: failed to open \"{}\" for input: No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let reference = tmp_path("missing-input-ref", "fa");
    std::fs::write(&reference, ">ref\nACGT\n").unwrap();
    let reference_arg = reference.to_string_lossy();
    let calmd = run_args(&["calmd", missing_arg.as_ref(), reference_arg.as_ref()]);
    assert_eq!(calmd.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&calmd.stderr),
        format!(
            "{}samtools calmd: Failed to open input file '{}': No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    for command in ["fastq", "fasta", "bam2fq"] {
        let output = run_missing(command, &missing);
        assert_eq!(output.status.code(), Some(1), "{command}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            format!(
                "{}samtools bam2fq: Cannot read file \"{}\": No such file or directory\n",
                hts_missing_line(&missing),
                missing.display()
            ),
            "{command}"
        );
    }

    let header = tmp_path("missing-input-header", "sam");
    let missing_header = tmp_path("definitely-missing-header", "hdr");
    let _ = std::fs::remove_file(&missing_header);
    std::fs::write(&header, "@HD\tVN:1.6\n").unwrap();
    let header_arg = header.to_string_lossy();
    let missing_header_arg = missing_header.to_string_lossy();

    let reheader_missing_header = run_args(&[
        "reheader",
        missing_header_arg.as_ref(),
        missing_arg.as_ref(),
    ]);
    assert_eq!(reheader_missing_header.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&reheader_missing_header.stderr),
        format!(
            "{}samtools reheader: fail to read the header from '{}': No such file or directory\n",
            hts_missing_line(&missing_header),
            missing_header.display()
        )
    );

    let reheader_missing_input = run_args(&["reheader", header_arg.as_ref(), missing_arg.as_ref()]);
    assert_eq!(reheader_missing_input.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&reheader_missing_input.stderr),
        format!(
            "{}samtools reheader: fail to open file '{}': No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let missing_fq = tmp_path("missing-input", "fq");
    let _ = std::fs::remove_file(&missing_fq);
    let import = run_missing("import", &missing_fq);
    assert_eq!(import.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&import.stderr),
        format!(
            "{}{}: No such file or directory\n",
            hts_missing_line(&missing_fq),
            missing_fq.display()
        )
    );

    let merge = run_args(&["merge", out_arg.as_ref(), missing_arg.as_ref()]);
    assert_eq!(merge.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&merge.stderr),
        format!(
            "{}samtools merge: fail to open \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let mpileup = run_missing("mpileup", &missing);
    assert_eq!(mpileup.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&mpileup.stderr),
        format!(
            "{}[mpileup] failed to open {}: No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let reference_cmd = run_missing("reference", &missing);
    assert_eq!(reference_cmd.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&reference_cmd.stderr),
        format!(
            "{}samtools reference: failed to open file '{}': No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let missing_cram = tmp_path("missing-input", "cram");
    let _ = std::fs::remove_file(&missing_cram);
    let cram_size = run_missing("cram-size", &missing_cram);
    assert_eq!(cram_size.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&cram_size.stderr),
        format!(
            "samtools cram_size: failed to open file '{}': No such file or directory\n",
            missing_cram.display()
        )
    );

    let phase = run_missing("phase", &missing);
    assert_eq!(phase.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&phase.stderr),
        format!(
            "{}samtools phase: Couldn't open '{}': No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let targetcut = run_args(&[
        "targetcut",
        "-f",
        reference_arg.as_ref(),
        missing_arg.as_ref(),
    ]);
    assert_eq!(targetcut.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&targetcut.stderr),
        format!(
            "{}samtools targetcut: can't open \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let tview = run_missing("tview", &missing);
    assert_eq!(tview.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&tview.stderr),
        format!(
            "{}samtools tview: can't open \"{}\": No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let depad = run_missing("depad", &missing);
    assert_eq!(depad.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&depad.stderr),
        format!(
            "{}samtools depad: failed to open \"{}\" for reading: No such file or directory\n",
            hts_missing_line(&missing),
            missing.display()
        )
    );

    let _ = std::fs::remove_file(out);
    let _ = std::fs::remove_file(reference);
    let _ = std::fs::remove_file(header);
    let _ = std::fs::remove_file(amplicon_bed);
    let _ = std::fs::remove_file(amplicon_out);
}

#[test]
fn collate_requires_explicit_output_destination() {
    let dir = std::env::temp_dir().join(format!(
        "samtools-rs-cli-collate-no-output-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let input = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../samtools/test/dat/test_input_1_a.bam")
        .canonicalize()
        .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_samtools"))
        .arg("collate")
        .arg(&input)
        .current_dir(&dir)
        .output()
        .expect("run samtools collate");

    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert!(!dir.join("collated.bam").exists());
    assert!(String::from_utf8_lossy(&output.stderr).starts_with("Usage: samtools collate "));

    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn view_no_sq_header_stderr_matches_upstream() {
    let input = tmp_path("no-sq-header", "sam");
    std::fs::write(
        &input,
        "@HD\tVN:1.6\n\
         r1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\t!\n",
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_samtools"))
        .arg("view")
        .arg(&input)
        .output()
        .expect("run samtools view");

    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        format!(
            "[E::sam_parse1] no SQ lines present in the header\n\
             [W::sam_read1_sam] Parse error at line 2\n\
             samtools view: error reading file \"{}\"\n",
            input.display()
        )
    );

    let _ = std::fs::remove_file(input);
}
