//! Error type for samtools-rs.
//!
//! samtools subcommands ultimately return an exit code, so most errors are
//! reported via [`crate::diagnostics::print_error`] and exit-code conversion.
//! This type is the structured form used during processing.

use std::io;

/// Result alias used throughout samtools-rs.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors surfaced from samtools-rs subcommands.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Underlying I/O error.
    #[error("{0}")]
    Io(#[from] io::Error),

    /// A command-line argument was invalid.
    #[error("{0}")]
    InvalidArg(String),

    /// A SAM/BAM/CRAM/VCF/BCF file could not be parsed.
    #[error("{0}")]
    Parse(String),

    /// A required input file was missing or not openable.
    #[error("{0}")]
    Open(String),

    /// A generic message error.
    #[error("{0}")]
    Other(String),
}

impl Error {
    /// Construct an [`Error::InvalidArg`] from any message.
    pub fn invalid_arg(msg: impl Into<String>) -> Self {
        Error::InvalidArg(msg.into())
    }

    /// Construct an [`Error::Other`] from any message.
    pub fn other(msg: impl Into<String>) -> Self {
        Error::Other(msg.into())
    }
}
