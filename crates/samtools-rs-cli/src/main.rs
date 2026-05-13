//! `samtools` binary front-end.
//!
//! Forwards `argv` to [`samtools_rs::run`]. The library crate owns all
//! dispatch logic; this binary is intentionally minimal so that integration
//! tests can invoke the library directly without spawning a subprocess.

use std::env;
use std::process::ExitCode;

fn main() -> ExitCode {
    let args = env::args_os().collect();
    samtools_rs::run(args)
}
