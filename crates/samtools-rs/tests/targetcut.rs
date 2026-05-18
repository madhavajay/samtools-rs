use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Mutex;

use samtools_rs::commands::targetcut;
use samtools_rs::run as samtools_run;

static GLOBAL_ARGS_LOCK: Mutex<()> = Mutex::new(());

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
        "samtools-rs-targetcut-it-{name}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&p);
    std::fs::create_dir_all(&p).unwrap();
    p
}

fn write_single_read_sam(path: &Path, len: usize, qual_byte: u8) {
    let seq = "A".repeat(len);
    let qual = std::iter::repeat_n(char::from(qual_byte), len).collect::<String>();
    std::fs::write(
        path,
        format!(
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:chr1\tLN:{len}\n\
             r1\t0\tchr1\t1\t60\t{len}M\t*\t0\t0\t{seq}\t{qual}\n"
        ),
    )
    .unwrap();
}

#[test]
fn targetcut_cli_writes_supported_interval_to_output_file() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("happy");
    let sam = tmp.join("in.sam");
    let out = tmp.join("targetcut.sam");
    write_single_read_sam(&sam, 2600, b'I');

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "targetcut",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ],
        ))),
        0
    );

    let text = std::fs::read_to_string(out).unwrap();
    let fields: Vec<&str> = text.trim_end().split('\t').collect();
    assert_eq!(fields.len(), 11);
    assert!(fields[0].starts_with("chr1:"));
    assert_eq!(fields[1], "0");
    assert_eq!(fields[2], "chr1");
    assert_eq!(fields[4], "60");
    assert!(fields[5].ends_with('M'));
    assert!(fields[9].chars().all(|c| c == 'A'));
    assert_eq!(fields[9].len(), fields[10].len());
}

#[test]
fn targetcut_cli_min_base_quality_filters_low_quality_reads() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("min-baseq");
    let sam = tmp.join("in.sam");
    let out = tmp.join("targetcut.sam");
    write_single_read_sam(&sam, 2600, b'!');

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "targetcut",
                "-Q",
                "13",
                "-o",
                out.to_str().unwrap(),
                sam.to_str().unwrap(),
            ],
        ))),
        0
    );

    assert_eq!(std::fs::read_to_string(out).unwrap(), "");
}

#[test]
fn targetcut_cli_rejects_missing_input_and_unknown_option() {
    assert_ne!(exit_to_u8(targetcut::main(&argv("targetcut", &[]))), 0);
    assert_ne!(
        exit_to_u8(targetcut::main(&argv("targetcut", &["--not-an-option"]))),
        0
    );
}
