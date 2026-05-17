use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Mutex;

use samtools_rs::commands::phase;
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
        "samtools-rs-phase-it-{name}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&p);
    std::fs::create_dir_all(&p).unwrap();
    p
}

fn write_phase_sam(path: &Path) {
    let sam = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:12
r0\t0\tchr1\t1\t60\t12M\t*\t0\t0\tAAAAAACCCCCC\tFFFFFFFFFFFF
r1\t0\tchr1\t1\t60\t12M\t*\t0\t0\tAAAAAACCCCCC\tFFFFFFFFFFFF
r2\t0\tchr1\t1\t60\t12M\t*\t0\t0\tCCCCCCAAAAAA\tFFFFFFFFFFFF
r3\t0\tchr1\t1\t60\t12M\t*\t0\t0\tCCCCCCAAAAAA\tFFFFFFFFFFFF
";
    std::fs::write(path, sam).unwrap();
}

#[test]
fn phase_cli_dispatches_and_writes_split_bams() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("split");
    let input = tmp.join("in.sam");
    let prefix = tmp.join("phase-out");
    write_phase_sam(&input);

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "phase",
                "-q",
                "1",
                "-Q",
                "1",
                "-k",
                "3",
                "--no-PG",
                "-b",
                prefix.to_str().unwrap(),
                input.to_str().unwrap(),
            ],
        ))),
        0
    );

    for middle in ["0", "1", "chimera"] {
        let path = tmp.join(format!("phase-out.{middle}.bam"));
        assert!(path.exists(), "{} missing", path.display());
    }
}

#[test]
fn phase_cli_rejects_missing_input_and_invalid_k() {
    assert_ne!(exit_to_u8(phase::main(&argv("phase", &[]))), 0);
    assert_ne!(
        exit_to_u8(phase::main(&argv("phase", &["-k", "0", "in.sam"]))),
        0
    );
}
