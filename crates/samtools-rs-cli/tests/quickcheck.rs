use std::process::Command;

fn fixtures_dir() -> std::path::PathBuf {
    let manifest = env!("CARGO_MANIFEST_DIR");
    std::path::PathBuf::from(manifest)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("repos")
        .join("samtools")
        .join("test")
        .join("quickcheck")
}

#[test]
fn quickcheck_verbose_output_matches_upstream_all_expected() {
    let fixtures = fixtures_dir();
    let expected = std::fs::read(fixtures.join("all.expected")).unwrap();
    let inputs = [
        "1.quickcheck.badeof.bam",
        "2.quickcheck.badheader.bam",
        "3.quickcheck.ok.bam",
        "4.quickcheck.ok.bam",
        "5.quickcheck.scramble30.truncated.cram",
        "6.quickcheck.cram21.ok.cram",
        "7.quickcheck.cram30.ok.cram",
        "8.quickcheck.cram21.truncated.cram",
        "9.quickcheck.cram30.truncated.cram",
        "10.quickcheck.notargets.bam",
    ];

    let output = Command::new(env!("CARGO_BIN_EXE_samtools"))
        .current_dir(&fixtures)
        .arg("quickcheck")
        .arg("-v")
        .args(inputs)
        .output()
        .unwrap_or_else(|e| panic!("run samtools quickcheck -v: {e}"));

    assert_eq!(output.status.code(), Some(28));
    assert_eq!(output.stdout, expected);
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "\
1.quickcheck.badeof.bam was missing EOF block when one should be present.\n\
2.quickcheck.badheader.bam was not identified as sequence data.\n\
5.quickcheck.scramble30.truncated.cram was missing EOF block when one should be present.\n\
8.quickcheck.cram21.truncated.cram was missing EOF block when one should be present.\n\
9.quickcheck.cram30.truncated.cram was missing EOF block when one should be present.\n\
10.quickcheck.notargets.bam had no targets in header.\n"
    );
}
