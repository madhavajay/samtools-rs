//! samtools version string emitted in `@PG VN:` and `--version`.
//!
//! Pinned to the upstream samtools tag tracked by the `samtools/` submodule.
//! When the submodule is bumped, update [`SAMTOOLS_VERSION`] to match the
//! upstream `version.sh` output for that commit.

/// The samtools version string emitted in `@PG VN:` lines and the
/// `samtools --version` banner.
pub const SAMTOOLS_VERSION: &str = "1.23.1";

/// htslib version string emitted in the `--version` banner. Reported as
/// the htslib-rs crate version since the underlying implementation is the
/// sibling `htslib-rs` workspace, not C HTSlib.
pub const HTSLIB_VERSION: &str = "1.23.1+htslib-rs";

/// Returns [`SAMTOOLS_VERSION`].
pub fn samtools_version() -> &'static str {
    SAMTOOLS_VERSION
}

/// Returns [`HTSLIB_VERSION`].
pub fn htslib_version() -> &'static str {
    HTSLIB_VERSION
}
