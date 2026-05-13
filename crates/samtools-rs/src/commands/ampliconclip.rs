//! `samtools ampliconclip` — not yet implemented in samtools-rs.

use std::ffi::OsString;
use std::process::ExitCode;

/// Entry point for `samtools ampliconclip`.
pub fn main(_args: &[OsString]) -> ExitCode {
    super::not_implemented("ampliconclip")
}
