//! `samtools targetcut` — not yet implemented in samtools-rs.

use std::ffi::OsString;
use std::process::ExitCode;

/// Entry point for `samtools targetcut`.
pub fn main(_args: &[OsString]) -> ExitCode {
    super::not_implemented("targetcut")
}
