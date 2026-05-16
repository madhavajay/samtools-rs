//! `samtools cram-size` — per-Content-ID / Data-Series CRAM block size
//! and codec inventory.
//!
//! Faithful port of `samtools/cram_size.c`. Built incrementally; this
//! revision lands the compression-method detail decoder (the port's
//! crux: a Rust port of htslib `cram_expand_method` +
//! `comp_method2expanded`, with the verbatim `comp_method2char` /
//! `comp_method2str` tables) as a unit-tested module. The block walk
//! consumes `htslib_rs`'s `noodles` `Container::blocks()` +
//! `CompressionHeader` inventory surface; the aggregation/formatting
//! and the `-e` "Container encodings" dump are wired on top in
//! subsequent steps, verified byte-exact against
//! `test/cram_size/cram_size.reg`.

use std::ffi::OsString;
use std::process::ExitCode;

// Staged sub-step: the decoder is exercised by unit tests and wired
// into the block walk in the next revision, so it is not yet
// referenced by `main`.
#[allow(dead_code)]
pub(crate) mod method {
    //! CRAM compression-method detail decoding: port of htslib
    //! `cram_expand_method` + samtools `comp_method2expanded` and the
    //! `comp_method2char` / `comp_method2str` tables
    //! (`cram_size.c:136-238`, `cram_external.c:568-651`).

    /// Coarse CRAM block compression method (htslib
    /// `enum cram_block_method`), as resolved by the container reader.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub enum Coarse {
        Raw,
        Gzip,
        Bzip2,
        Lzma,
        Rans4x8,
        RansNx16,
        Arith,
        Fqzcomp,
        Tok3,
    }

    // htscodecs `RANS_ORDER_*` flag bits (rANS_static4x16.h).
    const RANS_ORDER_X32: u8 = 0x04;
    const RANS_ORDER_STRIPE: u8 = 0x08;
    const RANS_ORDER_CAT: u8 = 0x20;
    const RANS_ORDER_RLE: u8 = 0x40;
    const RANS_ORDER_PACK: u8 = 0x80;

    /// `comp_expanded` index space (`cram_size.c:55-130`). The integer
    /// value is the index into [`COMP_CHAR`] / [`COMP_STR`]. Standard
    /// CRAM methods occupy 0..=8 (matching htslib `cram_block_method`),
    /// then the localized variants follow in the exact order of the C
    /// `enum comp_expanded`.
    pub const COMP_RAW: usize = 0;
    pub const COMP_GZIP: usize = 1;
    pub const COMP_BZIP2: usize = 2;
    pub const COMP_LZMA: usize = 3;
    pub const COMP_RANS8: usize = 4;
    pub const COMP_RANS16: usize = 5;
    pub const COMP_ARITH: usize = 6;
    pub const COMP_FQZ: usize = 7;
    pub const COMP_TOK3: usize = 8;
    pub const COMP_GZIP_1: usize = 9;
    pub const COMP_GZIP_9: usize = 10;
    pub const COMP_BZIP2_1: usize = 11;
    pub const COMP_RANS4X8_O0: usize = 20;
    pub const COMP_RANS4X8_O1: usize = 21;
    pub const COMP_RANS4X16_O0: usize = 22;
    pub const COMP_RANSNX16_STRIPE: usize = 38;
    pub const COMP_RANSNX16_CAT: usize = 39;
    pub const COMP_ARITH_O0: usize = 40;
    pub const COMP_ARITH_STRIPE: usize = 48;
    pub const COMP_ARITH_CAT: usize = 49;
    pub const COMP_ARITH_EXT: usize = 50;
    pub const COMP_TOK3_RANS: usize = 51;
    pub const COMP_TOK3_ARITH: usize = 52;
    pub const COMP_MAX: usize = 53;

    /// `comp_method2char` (`cram_size.c:192-199`), verbatim.
    pub const COMP_CHAR: &[u8; COMP_MAX] = b".gblr0afn_GbbbbbbbbBrR010101014545454582aAaAaAaAaaanN";

    /// `comp_method2str` (`cram_size.c:202-238`), verbatim.
    pub const COMP_STR: [&str; COMP_MAX] = [
        "raw",
        "gzip",
        "bzip2",
        "lzma",
        "r4x8",
        "rNx16",
        "arith",
        "fqzcomp",
        "tok3",
        "gzip-min",
        "gzip-max",
        "bzip2-1",
        "bzip2-2",
        "bzip2-3",
        "bzip2-4",
        "bzip2-5",
        "bzip2-6",
        "bzip2-7",
        "bzip2-8",
        "bzip2-9",
        "r4x8-o0",
        "r4x8-o1",
        "r4x16-o0",
        "r4x16-o1",
        "r4x16-o0R",
        "r4x16-o1R",
        "r4x16-o0P",
        "r4x16-o1P",
        "r4x16-o0PR",
        "r4x16-o1PR",
        "r32x16-o0",
        "r32x16-o1",
        "r32x16-o0R",
        "r32x16-o1R",
        "r32x16-o0P",
        "r32x16-o1P",
        "r32x16-o0PR",
        "r32x16-o1PR",
        "rNx16-xo0",
        "rNx16-cat",
        "arith-o0",
        "arith-o1",
        "arith-o0R",
        "arith-o1R",
        "arith-o0P",
        "arith-o1P",
        "arith-o0PR",
        "arith-o1PR",
        "arith-stripe",
        "arith-cat",
        "arith-ext",
        "tok3-rans",
        "tok3-arith",
    ];

    /// `cram_expand_method` + `comp_method2expanded`: given the coarse
    /// method and the block's stored (compressed) bytes, return the
    /// expanded `comp_expanded` index.
    pub fn expand(coarse: Coarse, data: &[u8]) -> usize {
        let b = |i: usize| data.get(i).copied().unwrap_or(0);
        match coarse {
            Coarse::Raw => COMP_RAW,
            Coarse::Lzma => COMP_LZMA,
            Coarse::Fqzcomp => COMP_FQZ,
            Coarse::Gzip => {
                if data.len() > 8 {
                    match b(8) {
                        4 => COMP_GZIP_1, // level 1 -> gzip-min
                        2 => COMP_GZIP_9, // level 9 -> gzip-max
                        _ => COMP_GZIP,   // level 5 -> gzip
                    }
                } else {
                    COMP_GZIP
                }
            }
            Coarse::Bzip2 => {
                if data.len() > 3 && b(3) >= b'1' && b(3) <= b'9' {
                    COMP_BZIP2_1 + (b(3) - b'1') as usize
                } else {
                    COMP_BZIP2
                }
            }
            Coarse::Rans4x8 => {
                if !data.is_empty() && b(0) == 1 {
                    COMP_RANS4X8_O1
                } else {
                    COMP_RANS4X8_O0
                }
            }
            Coarse::RansNx16 => {
                if data.is_empty() {
                    return COMP_RANS16;
                }
                let f = b(0);
                if f & RANS_ORDER_STRIPE != 0 {
                    return COMP_RANSNX16_STRIPE;
                }
                if f & RANS_ORDER_CAT != 0 {
                    return COMP_RANSNX16_CAT;
                }
                // bit 0 order, bit 1 rle, bit 2 pack, bit 3 32x16.
                let mut c = COMP_RANS4X16_O0;
                c += (f & 1) as usize;
                c += if f & RANS_ORDER_RLE != 0 { 2 } else { 0 };
                c += if f & RANS_ORDER_PACK != 0 { 4 } else { 0 };
                c += if f & RANS_ORDER_X32 != 0 { 8 } else { 0 };
                c
            }
            Coarse::Arith => {
                if data.is_empty() {
                    return COMP_ARITH;
                }
                let f = b(0);
                if f & RANS_ORDER_STRIPE != 0 {
                    return COMP_ARITH_STRIPE;
                }
                if f & RANS_ORDER_CAT != 0 {
                    return COMP_ARITH_CAT;
                }
                if f & 4 != 0 {
                    return COMP_ARITH_EXT;
                }
                let mut c = COMP_ARITH_O0;
                c += (f & 3) as usize;
                c += if f & RANS_ORDER_RLE != 0 { 2 } else { 0 };
                c += if f & RANS_ORDER_PACK != 0 { 4 } else { 0 };
                c
            }
            Coarse::Tok3 => {
                if data.len() > 8 && b(8) == 1 {
                    COMP_TOK3_ARITH
                } else {
                    COMP_TOK3_RANS
                }
            }
        }
    }

    /// Compact one-char method flag (`comp_method2char[comp]`).
    pub fn flag_char(comp: usize) -> char {
        COMP_CHAR[comp] as char
    }

    /// Long method name (`comp_method2str[comp]`).
    pub fn name(comp: usize) -> &'static str {
        COMP_STR[comp]
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn tables_have_comp_max_entries_and_anchor_chars() {
            assert_eq!(COMP_CHAR.len(), COMP_MAX);
            assert_eq!(COMP_STR.len(), COMP_MAX);
            // Anchors used by the upstream cram_size fixtures.
            assert_eq!(flag_char(COMP_RAW), '.');
            assert_eq!(flag_char(COMP_GZIP), 'g');
            assert_eq!(flag_char(COMP_GZIP_1), '_');
            assert_eq!(flag_char(COMP_BZIP2_1 + 5), 'b'); // bzip2-6
            assert_eq!(flag_char(COMP_RANS4X16_O0), '0');
            assert_eq!(flag_char(COMP_RANSNX16_CAT), '2');
            assert_eq!(flag_char(COMP_RANSNX16_STRIPE), '8');
            assert_eq!(name(COMP_RANSNX16_CAT), "rNx16-cat");
            assert_eq!(name(COMP_RANSNX16_STRIPE), "rNx16-xo0");
            assert_eq!(name(COMP_BZIP2_1 + 5), "bzip2-6");
            assert_eq!(name(COMP_GZIP_1), "gzip-min");
        }

        #[test]
        fn gzip_level_from_xfl_byte() {
            let mut d = vec![0u8; 12];
            d[8] = 4;
            assert_eq!(expand(Coarse::Gzip, &d), COMP_GZIP_1);
            d[8] = 2;
            assert_eq!(expand(Coarse::Gzip, &d), COMP_GZIP_9);
            d[8] = 0;
            assert_eq!(expand(Coarse::Gzip, &d), COMP_GZIP);
        }

        #[test]
        fn bzip2_level_from_header_digit() {
            let d = b"BZh6abcd";
            assert_eq!(expand(Coarse::Bzip2, d), COMP_BZIP2_1 + 5);
            assert_eq!(name(expand(Coarse::Bzip2, d)), "bzip2-6");
        }

        #[test]
        fn rans_nx16_flag_byte() {
            assert_eq!(
                expand(Coarse::RansNx16, &[RANS_ORDER_CAT]),
                COMP_RANSNX16_CAT
            );
            assert_eq!(
                expand(Coarse::RansNx16, &[RANS_ORDER_STRIPE]),
                COMP_RANSNX16_STRIPE
            );
            assert_eq!(expand(Coarse::RansNx16, &[0]), COMP_RANS4X16_O0);
            assert_eq!(expand(Coarse::RansNx16, &[1]), COMP_RANS4X16_O0 + 1);
        }

        #[test]
        fn rans4x8_order_byte() {
            assert_eq!(expand(Coarse::Rans4x8, &[0]), COMP_RANS4X8_O0);
            assert_eq!(expand(Coarse::Rans4x8, &[1]), COMP_RANS4X8_O1);
        }
    }
}

/// Entry point for `samtools cram-size`.
pub fn main(_args: &[OsString]) -> ExitCode {
    // Aggregation/formatting + `-e` encodings dump are wired on top of
    // the `method` decoder + the noodles `Container::blocks()` /
    // `CompressionHeader` inventory surface in the next step.
    super::not_implemented("cram-size")
}
