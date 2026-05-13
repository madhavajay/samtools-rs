//! `samtools fqidx` — FASTQ index build. Delegates to [`super::faidx::fqidx_main`].

use std::ffi::OsString;
use std::process::ExitCode;

/// Entry point for `samtools fqidx`.
pub fn main(args: &[OsString]) -> ExitCode {
    super::faidx::fqidx_main(args)
}
