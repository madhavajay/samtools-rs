use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Mutex;

use samtools_rs::commands::tview;
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
        "samtools-rs-tview-it-{name}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&p);
    std::fs::create_dir_all(&p).unwrap();
    p
}

fn write_tview_sam(path: &Path) {
    let sam = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:20
r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII
r2\t16\tchr1\t5\t60\t4M1I4M\t*\t0\t0\tTGCATGCAT\tIIIIIIIII
";
    std::fs::write(path, sam).unwrap();
}

#[test]
fn tview_text_mode_dispatches_for_region() {
    let _guard = GLOBAL_ARGS_LOCK.lock().unwrap();
    let tmp = tmp_dir("text");
    let sam = tmp.join("in.sam");
    write_tview_sam(&sam);

    assert_eq!(
        exit_to_u8(samtools_run(argv(
            "samtools",
            &[
                "tview",
                "-d",
                "T",
                "-p",
                "chr1:1",
                "-w",
                "20",
                sam.to_str().unwrap(),
            ],
        ))),
        0
    );
}

#[test]
fn tview_rejects_non_text_mode_and_malformed_region() {
    let tmp = tmp_dir("errors");
    let sam = tmp.join("in.sam");
    write_tview_sam(&sam);

    assert_ne!(
        exit_to_u8(tview::main(&argv(
            "tview",
            &["-p", "chr1:1", sam.to_str().unwrap()],
        ))),
        0
    );
    assert_ne!(
        exit_to_u8(tview::main(&argv(
            "tview",
            &["-d", "T", "-p", "chr1", sam.to_str().unwrap()],
        ))),
        0
    );
}

#[test]
fn tview_help_succeeds() {
    assert_eq!(exit_to_u8(tview::main(&argv("tview", &["--help"]))), 0);
}
