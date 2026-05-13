//! `samtools mpileup` — not yet implemented in samtools-rs.

use std::ffi::OsString;
use std::process::ExitCode;

/// Entry point for `samtools mpileup`.
pub fn main(_args: &[OsString]) -> ExitCode {
    super::not_implemented("mpileup")
}
