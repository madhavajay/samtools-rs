//! Integration tests for `samtools flags`.

use std::ffi::OsString;

use samtools_rs::commands::flags;

fn run(args: &[&str]) -> u8 {
    let argv: Vec<OsString> = std::iter::once(OsString::from("flags"))
        .chain(args.iter().map(OsString::from))
        .collect();
    let code = flags::main(&argv);
    let s = format!("{:?}", code);
    s.chars()
        .filter(|c| c.is_ascii_digit())
        .collect::<String>()
        .parse()
        .unwrap_or(255)
}

#[test]
fn flags_invalid_returns_1() {
    assert_eq!(run(&["NOPE"]), 1);
}

#[test]
fn flags_numeric_succeeds() {
    assert_eq!(run(&["0x10"]), 0);
}

#[test]
fn flags_names_succeed() {
    assert_eq!(run(&["PAIRED,READ1"]), 0);
}
