//! `samtools reference` — not yet implemented in samtools-rs.

use std::ffi::OsString;
use std::process::ExitCode;

/// Entry point for `samtools reference`.
pub fn main(_args: &[OsString]) -> ExitCode {
    super::not_implemented("reference")
}
