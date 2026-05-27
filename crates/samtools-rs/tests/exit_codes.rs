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
        .join("repos")
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

/// Drives a subcommand through the top-level dispatcher (mirrors the
/// real `samtools <sub> ...` entry point) and returns its exit code.
fn sub(name: &str, rest: &[&str]) -> u8 {
    let mut full = vec!["samtools", name];
    full.extend_from_slice(rest);
    exit_to_u8(samtools_run(
        full.iter().map(OsString::from).collect::<Vec<_>>(),
    ))
}

#[test]
fn missing_input_returns_nonzero_across_subcommands() {
    let missing = "/tmp/samtools-rs-definitely-missing-input.bam";
    // Error class: input path does not exist. Upstream samtools exits
    // non-zero (1) for these; `quickcheck` uses its own bitmask and
    // returns 2 when a listed file is missing/unopenable.
    for command in [
        "flagstat",
        "idxstats",
        "stats",
        "depth",
        "faidx",
        "fixmate",
        "markdup",
        "calmd",
        "collate",
        "addreplacerg",
    ] {
        assert_eq!(
            sub(command, &[missing]),
            1,
            "{command} on a missing input must exit 1"
        );
    }
    assert_eq!(
        sub("quickcheck", &[missing]),
        2,
        "quickcheck on a missing input must exit with bit 2"
    );
}

#[test]
fn unknown_option_returns_one() {
    let tmp = tmp_dir("unknown-option");
    let sam = tmp.join("ok.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\n@SQ\tSN:c\tLN:5\nr1\t0\tc\t1\t60\t1M\t*\t0\t0\tA\t!\n",
    )
    .unwrap();
    let sam = sam.to_str().unwrap();

    for command in ["view", "sort", "collate"] {
        assert_eq!(
            sub(command, &["--definitely-not-a-flag", sam]),
            1,
            "{command} with an unknown option must exit 1"
        );
    }
}

#[test]
fn unwritable_output_path_returns_one() {
    let tmp = tmp_dir("unwritable-out");
    let sam = tmp.join("ok.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\n@SQ\tSN:c\tLN:5\nr1\t0\tc\t1\t60\t1M\t*\t0\t0\tA\t!\n",
    )
    .unwrap();

    assert_eq!(
        sub(
            "view",
            &[
                "-b",
                "-o",
                "/no/such/directory/out.bam",
                sam.to_str().unwrap(),
            ],
        ),
        1,
        "view writing into a nonexistent directory must exit 1"
    );
}

#[test]
fn unrecognized_input_format_returns_one() {
    let tmp = tmp_dir("garbage-format");
    let garbage = tmp.join("garbage.dat");
    std::fs::write(&garbage, b"not a bam or sam or cram at all\n").unwrap();
    let garbage = garbage.to_str().unwrap();

    // Error class: input is not a recognizable alignment format.
    for command in ["view", "flagstat"] {
        assert_eq!(
            sub(command, &[garbage]),
            1,
            "{command} on an unrecognized format must exit 1"
        );
    }
}

#[test]
fn index_rejects_plain_text_sam_input() {
    // `index` needs a coordinate-sorted BGZF/CRAM file; a plain-text
    // SAM cannot be indexed and must exit non-zero rather than
    // silently producing a broken index.
    let tmp = tmp_dir("index-text-sam");
    let sam = tmp.join("plain.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:5\nr1\t0\tc\t1\t60\t1M\t*\t0\t0\tA\t!\n",
    )
    .unwrap();

    assert_eq!(sub("index", &[sam.to_str().unwrap()]), 1);
}
