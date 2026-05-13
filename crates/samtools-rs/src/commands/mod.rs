//! samtools subcommand implementations. One module per subcommand, each
//! exposing a `pub fn main(args: &[OsString]) -> ExitCode` entry point.
//!
//! Subcommands not yet implemented return exit code 2 after writing a
//! `not yet implemented` notice to stderr via [`not_implemented`]. The
//! dispatcher in [`crate::dispatch`] is the single source of truth for which
//! subcommands exist; modules listed there must exist here.

use std::io::{self, Write};
use std::process::ExitCode;

pub mod addreplacerg;
pub mod ampliconclip;
pub mod ampliconstats;
pub mod bedcov;
pub mod calmd;
pub mod cat;
pub mod checksum;
pub mod collate;
pub mod consensus;
pub mod coverage;
pub mod cram_size;
pub mod depad;
pub mod depth;
pub mod dict;
pub mod faidx;
pub mod fastq;
pub mod fixmate;
pub mod flags;
pub mod flagstat;
pub mod fqidx;
pub mod head;
pub mod idxstats;
pub mod import;
pub mod index;
pub mod markdup;
pub mod merge;
pub mod mpileup;
pub mod phase;
pub mod quickcheck;
pub mod reference;
pub mod reheader;
pub mod reset;
pub mod rmdup;
pub mod samples;
pub mod sort;
pub mod split;
pub mod stats;
pub mod targetcut;
pub mod view;

/// Stub used by subcommands that have not been ported yet. Writes a
/// `not yet implemented` line to stderr (prefixed with the subcommand name,
/// matching upstream's `[subcommand] message` convention) and returns
/// exit code 2.
pub(crate) fn not_implemented(name: &str) -> ExitCode {
    let _ = writeln!(
        io::stderr(),
        "samtools {}: subcommand not yet implemented in samtools-rs",
        name
    );
    ExitCode::from(2)
}
