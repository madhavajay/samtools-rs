use std::path::PathBuf;
use std::process::Command;

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("repos")
        .join("samtools")
        .join("test")
}

fn tmp_path(name: &str, extension: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "samtools-rs-cli-{name}-{}.{}",
        std::process::id(),
        extension
    ))
}

#[test]
fn fastq_dash_zero_only_keeps_paired_reads_on_suffixed_stdout() {
    let input = fixtures_dir().join("dat").join("bam2fq.001.sam");
    let other = tmp_path("fastq-dash-zero-other", "fq");
    let _ = std::fs::remove_file(&other);

    let output = Command::new(env!("CARGO_BIN_EXE_samtools"))
        .args(["fastq", "-0"])
        .arg(&other)
        .arg(input)
        .output()
        .expect("run samtools fastq");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.starts_with("@ref1_grp1_p001/1\n"));
    assert!(stdout.contains("@ref1_grp1_p001/2\n"));
    assert!(!stdout.contains("@ref1_grp1_p001\n"));
    assert!(std::fs::read(&other).unwrap().is_empty());
}
