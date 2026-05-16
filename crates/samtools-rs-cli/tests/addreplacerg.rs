use std::path::PathBuf;
use std::process::Command;

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("samtools")
        .join("test")
}

#[test]
fn unknown_dash_cap_r_error_matches_upstream_stderr() {
    let input = fixtures_dir().join("addrprg").join("1_fixup.sam");
    let expected = fixtures_dir()
        .join("addrprg")
        .join("3_fixup.sam.expected.err");

    let output = Command::new(env!("CARGO_BIN_EXE_samtools"))
        .args(["addreplacerg", "-O", "sam", "-R", "1#9"])
        .arg(input)
        .output()
        .expect("run samtools addreplacerg");

    assert!(!output.status.success());
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        std::fs::read_to_string(expected).unwrap()
    );
}
