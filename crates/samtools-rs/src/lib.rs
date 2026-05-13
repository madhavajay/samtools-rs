//! Pure Rust port of `samtools`.
//!
//! The CLI binary lives in the sibling `samtools-rs-cli` crate. This library
//! crate owns shared infrastructure (CLI dispatch, global args, `@PG` helpers,
//! diagnostics) and one module per subcommand under [`commands`].
//!
//! Format-level I/O routes through [`htslib_rs`], which itself delegates to
//! noodles. Direct use of noodles from this crate is reserved for cases that
//! have no HTSlib analogue; HTSlib-shaped helpers belong in `htslib-rs`.

pub mod aux_list;
pub mod bam_flag;
pub mod commands;
pub mod diagnostics;
pub mod dispatch;
pub mod error;
pub mod header_text;
pub mod native;
pub mod version;

pub use dispatch::run;
pub use error::{Error, Result};
