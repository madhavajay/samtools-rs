//! Diagnostic helpers matching upstream samtools' `print_error` /
//! `print_error_errno` stderr format.
//!
//! Upstream prints `[subcommand] message\n` to stderr, optionally appending
//! `: <strerror>` for the errno variant. We mirror that exact prefix because
//! some tests grep stderr for the subcommand tag.

use std::io::{self, Write};

/// Print `[subcommand] message` to stderr followed by a newline.
///
/// Matches upstream samtools' `print_error` in `sam_utils.c`.
pub fn print_error(subcommand: &str, message: impl AsRef<str>) {
    let _ = writeln!(
        io::stderr(),
        "samtools {}: {}",
        subcommand,
        message.as_ref()
    );
}

/// Print `[subcommand] message: <io-error>` to stderr followed by a newline.
///
/// Matches upstream samtools' `print_error_errno` in `sam_utils.c`.
pub fn print_error_errno(subcommand: &str, message: impl AsRef<str>, err: &io::Error) {
    let _ = writeln!(
        io::stderr(),
        "samtools {}: {}: {}",
        subcommand,
        message.as_ref(),
        err
    );
}
