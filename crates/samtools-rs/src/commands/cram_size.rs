//! `samtools cram-size` — not yet implemented in samtools-rs.

use std::ffi::OsString;
use std::process::ExitCode;

/// Entry point for `samtools cram-size`.
pub fn main(_args: &[OsString]) -> ExitCode {
    super::not_implemented("cram-size")
}
