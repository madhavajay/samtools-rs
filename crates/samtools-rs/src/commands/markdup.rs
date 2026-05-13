//! `samtools markdup` — not yet implemented in samtools-rs.

use std::ffi::OsString;
use std::process::ExitCode;

/// Entry point for `samtools markdup`.
pub fn main(_args: &[OsString]) -> ExitCode {
    super::not_implemented("markdup")
}
