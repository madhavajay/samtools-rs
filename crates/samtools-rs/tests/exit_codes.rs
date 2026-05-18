use std::ffi::OsString;
use std::path::PathBuf;
use std::process::ExitCode;

use samtools_rs::commands::{index, quickcheck, sort, view};
use samtools_rs::run as samtools_run;

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

fn argv(name: &str, rest: &[&str]) -> Vec<OsString> {
    std::iter::once(OsString::from(name))
        .chain(rest.iter().map(OsString::from))
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
    let p = std::env::temp_dir().join(format!(
        "samtools-rs-exit-codes-{name}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&p);
    std::fs::create_dir_all(&p).unwrap();
    p
}

#[test]
fn missing_input_paths_return_one_for_common_commands() {
    let missing = "/tmp/samtools-rs-definitely-missing-input.bam";
    assert_eq!(exit_to_u8(view::main(&argv("view", &[missing]))), 1);
    assert_eq!(exit_to_u8(index::main(&argv("index", &[missing]))), 1);
    assert_eq!(exit_to_u8(sort::main(&argv("sort", &[missing]))), 1);
}

#[test]
fn malformed_sam_cigar_returns_one() {
    let tmp = tmp_dir("malformed-cigar");
    let sam = tmp.join("bad.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\n\
         @SQ\tSN:chr1\tLN:10\n\
         r1\t0\tchr1\t1\t60\t1Z\t*\t0\t0\tA\t!\n",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(view::main(&argv("view", &[sam.to_str().unwrap()]))),
        1
    );
}

#[test]
fn mapped_sam_without_sq_header_returns_one() {
    let tmp = tmp_dir("no-sq-header");
    let sam = tmp.join("no-sq.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\n\
         r1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\t!\n",
    )
    .unwrap();

    assert_eq!(
        exit_to_u8(view::main(&argv("view", &[sam.to_str().unwrap()]))),
        1
    );
}

#[test]
fn quickcheck_preserves_bitmask_exit_codes() {
    let fixtures = fixtures_dir().join("quickcheck");
    let bad_header = fixtures.join("2.quickcheck.badheader.bam");
    let missing_eof = fixtures.join("1.quickcheck.badeof.bam");
    let truncated_cram = fixtures.join("9.quickcheck.cram30.truncated.cram");

    assert_eq!(
        exit_to_u8(quickcheck::main(&argv(
            "quickcheck",
            &[bad_header.to_str().unwrap()],
        ))),
        4
    );
    assert_eq!(
        exit_to_u8(quickcheck::main(&argv(
            "quickcheck",
            &[missing_eof.to_str().unwrap()],
        ))),
        16
    );
    assert_eq!(
        exit_to_u8(quickcheck::main(&argv(
            "quickcheck",
            &[truncated_cram.to_str().unwrap()],
        ))),
        16
    );
}

#[test]
fn unknown_top_level_command_returns_one() {
    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &["definitely-not-a-command"]
        ))),
        1
    );
}
