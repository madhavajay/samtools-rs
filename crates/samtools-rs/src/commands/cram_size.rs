//! `samtools cram-size` — per-Content-ID / Data-Series CRAM block size
//! and codec inventory.
//!
//! Faithful port of `samtools/cram_size.c`. The default and `-v`
//! (verbose) reports are **byte-exact** vs the upstream
//! `test/cram_size/cram_size.reg` fixtures (`normal.out`,
//! `verbose.out`): the `method` module ports htslib
//! `cram_expand_method` + `comp_method2expanded` + the verbatim
//! `comp_method2char`/`comp_method2str` tables; the block walk
//! consumes `htslib_rs`'s `noodles` `Container::blocks()` +
//! `CompressionHeader` inventory surface, builds the
//! content_id→data-series map (`cram_cid2ds`), aggregates by
//! content_id (default) / (content_id, method) (verbose), and emits
//! the summary. The `-e` "Container encodings" dump (which depends on
//! htslib's exact internal codec-iteration order) is the remaining
//! sub-step.

use std::ffi::OsString;
use std::process::ExitCode;

// Some `comp_expanded` constants/variants are part of the faithful
// table port but not all are reachable for the test fixtures.
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

/// Builds the content-id → data-series-code map (`cram_cid2ds`).
mod cid2ds {
    use htslib_rs::cram::container::CompressionHeader;
    use htslib_rs::cram::container::compression_header::data_series_encodings::DataSeries;
    use htslib_rs::cram::container::compression_header::encoding::codec::{
        Byte, ByteArray, Integer,
    };
    use std::collections::HashMap;

    /// `DataSeries` → its CRAM 2-letter code (reverse of the noodles
    /// `TryFrom<[u8; 2]>` table).
    fn ds_code(ds: DataSeries) -> [u8; 2] {
        use DataSeries::*;
        match ds {
            BamFlags => *b"BF",
            CramFlags => *b"CF",
            ReferenceSequenceIds => *b"RI",
            ReadLengths => *b"RL",
            AlignmentStarts => *b"AP",
            ReadGroupIds => *b"RG",
            Names => *b"RN",
            MateFlags => *b"MF",
            MateReferenceSequenceIds => *b"NS",
            MateAlignmentStarts => *b"NP",
            TemplateLengths => *b"TS",
            MateDistances => *b"NF",
            TagSetIds => *b"TL",
            FeatureCounts => *b"FN",
            FeatureCodes => *b"FC",
            FeaturePositionDeltas => *b"FP",
            DeletionLengths => *b"DL",
            StretchesOfBases => *b"BB",
            StretchesOfQualityScores => *b"QQ",
            BaseSubstitutionCodes => *b"BS",
            InsertionBases => *b"IN",
            ReferenceSkipLengths => *b"RS",
            PaddingLengths => *b"PD",
            HardClipLengths => *b"HC",
            SoftClipBases => *b"SC",
            MappingQualities => *b"MQ",
            Bases => *b"BA",
            QualityScores => *b"QS",
            ReservedTc => *b"TC",
            ReservedTn => *b"TN",
        }
    }

    fn int_ids(e: &Integer, out: &mut Vec<i32>) {
        if let Integer::External { block_content_id } = e {
            out.push(*block_content_id);
        }
    }
    fn byte_ids(e: &Byte, out: &mut Vec<i32>) {
        if let Byte::External { block_content_id } = e {
            out.push(*block_content_id);
        }
    }
    fn ba_ids(e: &ByteArray, out: &mut Vec<i32>) {
        match e {
            ByteArray::ByteArrayStop {
                block_content_id, ..
            } => out.push(*block_content_id),
            ByteArray::ByteArrayLength {
                len_encoding,
                value_encoding,
            } => {
                int_ids(len_encoding.get(), out);
                byte_ids(value_encoding.get(), out);
            }
        }
    }

    /// `d` is the upstream-encoded code: 2-letter `(a<<8)|b`, or for a
    /// tag the 3-byte `(a<<16)|(b<<8)|t`. Formatted as 2 or 3 chars.
    pub fn code_str(d: i64) -> String {
        if d > 0xffff {
            format!(
                "{}{}{}",
                ((d >> 16) & 0xff) as u8 as char,
                ((d >> 8) & 0xff) as u8 as char,
                (d & 0xff) as u8 as char
            )
        } else {
            format!(
                "{}{}",
                ((d >> 8) & 0xff) as u8 as char,
                (d & 0xff) as u8 as char
            )
        }
    }

    /// content_id → ordered list of data-series codes (`int` encoded
    /// like htslib `cram_cid2ds`).
    pub fn build(map: &mut HashMap<i32, Vec<i64>>, ch: &CompressionHeader) {
        let dse = ch.data_series_encodings();
        // (DataSeries, content ids) in noodles field order.
        // `cram_update_cid2ds_map` keeps each (content_id, code) once
        // even across containers that repeat the same encodings.
        let push = |map: &mut HashMap<i32, Vec<i64>>, id: i32, code: i64| {
            let v = map.entry(id).or_default();
            if !v.contains(&code) {
                v.push(code);
            }
        };
        let mut add_int = |ds: DataSeries, e: Option<&_>| {
            if let Some(enc) = e {
                let mut ids = Vec::new();
                int_ids(
                    htslib_rs::cram::container::compression_header::Encoding::get(enc),
                    &mut ids,
                );
                let [a, b] = ds_code(ds);
                let code = ((a as i64) << 8) | b as i64;
                for id in ids {
                    push(map, id, code);
                }
            }
        };
        add_int(DataSeries::BamFlags, dse.bam_flags());
        add_int(DataSeries::CramFlags, dse.cram_flags());
        add_int(
            DataSeries::ReferenceSequenceIds,
            dse.reference_sequence_ids(),
        );
        add_int(DataSeries::ReadLengths, dse.read_lengths());
        add_int(DataSeries::AlignmentStarts, dse.alignment_starts());
        add_int(DataSeries::ReadGroupIds, dse.read_group_ids());
        add_int(DataSeries::MateFlags, dse.mate_flags());
        add_int(
            DataSeries::MateReferenceSequenceIds,
            dse.mate_reference_sequence_ids(),
        );
        add_int(DataSeries::MateAlignmentStarts, dse.mate_alignment_starts());
        add_int(DataSeries::TemplateLengths, dse.template_lengths());
        add_int(DataSeries::MateDistances, dse.mate_distances());
        add_int(DataSeries::TagSetIds, dse.tag_set_ids());
        add_int(DataSeries::FeatureCounts, dse.feature_counts());
        add_int(
            DataSeries::FeaturePositionDeltas,
            dse.feature_position_deltas(),
        );
        add_int(DataSeries::DeletionLengths, dse.deletion_lengths());
        add_int(
            DataSeries::ReferenceSkipLengths,
            dse.reference_skip_lengths(),
        );
        add_int(DataSeries::PaddingLengths, dse.padding_lengths());
        add_int(DataSeries::HardClipLengths, dse.hard_clip_lengths());
        add_int(DataSeries::MappingQualities, dse.mapping_qualities());

        let mut add_byte = |ds: DataSeries, e: Option<&_>| {
            if let Some(enc) = e {
                let mut ids = Vec::new();
                byte_ids(
                    htslib_rs::cram::container::compression_header::Encoding::get(enc),
                    &mut ids,
                );
                let [a, b] = ds_code(ds);
                let code = ((a as i64) << 8) | b as i64;
                for id in ids {
                    let v = map.entry(id).or_default();
                    if !v.contains(&code) {
                        v.push(code);
                    }
                }
            }
        };
        add_byte(DataSeries::FeatureCodes, dse.feature_codes());
        add_byte(
            DataSeries::BaseSubstitutionCodes,
            dse.base_substitution_codes(),
        );
        add_byte(DataSeries::Bases, dse.bases());
        add_byte(DataSeries::QualityScores, dse.quality_scores());

        let mut add_ba = |ds: DataSeries, e: Option<&_>| {
            if let Some(enc) = e {
                let mut ids = Vec::new();
                ba_ids(
                    htslib_rs::cram::container::compression_header::Encoding::get(enc),
                    &mut ids,
                );
                let [a, b] = ds_code(ds);
                let code = ((a as i64) << 8) | b as i64;
                for id in ids {
                    let v = map.entry(id).or_default();
                    if !v.contains(&code) {
                        v.push(code);
                    }
                }
            }
        };
        add_ba(DataSeries::Names, dse.names());
        add_ba(DataSeries::StretchesOfBases, dse.stretches_of_bases());
        add_ba(
            DataSeries::StretchesOfQualityScores,
            dse.stretches_of_quality_scores(),
        );
        add_ba(DataSeries::InsertionBases, dse.insertion_bases());
        add_ba(DataSeries::SoftClipBases, dse.soft_clip_bases());

        // Tag encodings: HashMap<tag_id, Encoding<ByteArray>>; the
        // key is the 3-byte tag code.
        let mut tags: Vec<(&i32, _)> = ch.tag_encodings().iter().collect();
        tags.sort_by_key(|(k, _)| **k);
        for (tag_id, enc) in tags {
            let mut ids = Vec::new();
            ba_ids(
                htslib_rs::cram::container::compression_header::Encoding::get(enc),
                &mut ids,
            );
            for id in ids {
                let v = map.entry(id).or_default();
                let code = *tag_id as i64;
                if !v.contains(&code) {
                    v.push(code);
                }
            }
        }
    }
}

/// `-e` "Container encodings" dump (`cram_describe_encodings` +
/// `cram_codec_describe`, `cram_codecs.c`).
mod describe {
    use htslib_rs::cram::container::CompressionHeader;
    use htslib_rs::cram::container::compression_header::data_series_encodings::DataSeries;
    use htslib_rs::cram::container::compression_header::encoding::Encoding;
    use htslib_rs::cram::container::compression_header::encoding::codec::{
        Byte, ByteArray, Integer,
    };

    fn int_desc(c: &Integer) -> String {
        match c {
            Integer::External { block_content_id } => format!("EXTERNAL(id={block_content_id})"),
            Integer::Huffman { alphabet, bit_lens } => huffman(alphabet, bit_lens),
            Integer::Beta { offset, len } => format!("BETA(offset={offset},nbits={len})"),
            Integer::Gamma { offset } => format!("GAMMA(offset={offset})"),
            Integer::Subexp { offset, k } => format!("SUBEXP(offset={offset},k={k})"),
            Integer::Golomb { offset, m } => format!("GOLOMB(offset={offset},m={m})"),
            Integer::GolombRice { offset, log2_m } => {
                format!("GOLOMB_RICE(offset={offset},log2m={log2_m})")
            }
        }
    }
    fn byte_desc(c: &Byte) -> String {
        match c {
            Byte::External { block_content_id } => format!("EXTERNAL(id={block_content_id})"),
            Byte::Huffman { alphabet, bit_lens } => huffman(alphabet, bit_lens),
        }
    }
    fn ba_desc(c: &ByteArray) -> String {
        match c {
            ByteArray::ByteArrayStop {
                stop_byte,
                block_content_id,
            } => format!("BYTE_ARRAY_STOP(stop={stop_byte},id={block_content_id})"),
            ByteArray::ByteArrayLength {
                len_encoding,
                value_encoding,
            } => format!(
                "BYTE_ARRAY_LEN(len_codec={{{}}},val_codec={{{}}}",
                int_desc(len_encoding.get()),
                byte_desc(value_encoding.get())
            ),
        }
    }
    fn huffman(alphabet: &[i32], bit_lens: &[u32]) -> String {
        let mut s = String::from("HUFFMAN(codes={");
        for (i, a) in alphabet.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&a.to_string());
        }
        s.push_str("},lengths={");
        for (i, l) in bit_lens.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&l.to_string());
        }
        s.push_str("})");
        s
    }

    /// One "Container encodings" block (no trailing blank line).
    pub fn container_encodings(ch: &CompressionHeader) -> String {
        let dse = ch.data_series_encodings();
        let mut out = String::from("Container encodings\n");
        let mut line = |code: &str, desc: String| {
            out.push('\t');
            out.push_str(code);
            out.push('\t');
            out.push_str(&desc);
            out.push('\n');
        };
        // `cram_DS_ID` order (cram_external.c:226-256).
        macro_rules! di {
            ($code:literal, $e:expr) => {
                if let Some(e) = $e {
                    line($code, int_desc(Encoding::get(e)));
                }
            };
        }
        macro_rules! db {
            ($code:literal, $e:expr) => {
                if let Some(e) = $e {
                    line($code, byte_desc(Encoding::get(e)));
                }
            };
        }
        macro_rules! da {
            ($code:literal, $e:expr) => {
                if let Some(e) = $e {
                    line($code, ba_desc(Encoding::get(e)));
                }
            };
        }
        da!("RN", dse.names());
        db!("QS", dse.quality_scores());
        da!("IN", dse.insertion_bases());
        da!("SC", dse.soft_clip_bases());
        di!("BF", dse.bam_flags());
        di!("CF", dse.cram_flags());
        di!("AP", dse.alignment_starts());
        di!("RG", dse.read_group_ids());
        di!("MQ", dse.mapping_qualities());
        di!("NS", dse.mate_reference_sequence_ids());
        di!("MF", dse.mate_flags());
        di!("TS", dse.template_lengths());
        di!("NP", dse.mate_alignment_starts());
        di!("NF", dse.mate_distances());
        di!("RL", dse.read_lengths());
        di!("FN", dse.feature_counts());
        db!("FC", dse.feature_codes());
        di!("FP", dse.feature_position_deltas());
        di!("DL", dse.deletion_lengths());
        db!("BA", dse.bases());
        db!("BS", dse.base_substitution_codes());
        di!("TL", dse.tag_set_ids());
        di!("RI", dse.reference_sequence_ids());
        di!("RS", dse.reference_skip_lengths());
        di!("PD", dse.padding_lengths());
        di!("HC", dse.hard_clip_lengths());
        da!("BB", dse.stretches_of_bases());
        da!("QQ", dse.stretches_of_quality_scores());
        let _ = (DataSeries::ReservedTc, DataSeries::ReservedTn); // legacy, unused
        // Tags follow in htslib `tag_encoding_map` order: 32 buckets,
        // `CRAM_MAP(a,b) = (a*3+b) & 31` on the two tag letters, with
        // prepend insertion (so reverse insertion order within a
        // bucket), iterated bucket 0..31. The `IndexMap` preserves
        // the compression-header insertion order this depends on.
        const N: usize = 32;
        let mut buckets: Vec<Vec<(i64, &_)>> = vec![Vec::new(); N];
        for (tag_id, enc) in ch.tag_encodings() {
            let d = *tag_id as i64;
            let c0 = (d >> 16) & 0xff;
            let c1 = (d >> 8) & 0xff;
            let b = ((c0 * 3 + c1) & (N as i64 - 1)) as usize;
            // prepend → newest first within the bucket
            buckets[b].insert(0, (d, enc));
        }
        for bucket in &buckets {
            for (d, enc) in bucket {
                line(&super::cid2ds::code_str(*d), ba_desc(Encoding::get(enc)));
            }
        }
        out
    }
}

/// Entry point for `samtools cram-size`.
pub fn main(args: &[OsString]) -> ExitCode {
    use htslib_rs::cram::container::block::{CompressionMethod, ContentType};
    use std::collections::HashMap;

    let mut verbose = false;
    let mut encodings = false;
    let mut input: Option<std::path::PathBuf> = None;
    let mut output: Option<std::path::PathBuf> = None;
    let mut it = args.iter().skip(1);
    while let Some(a) = it.next() {
        match a.to_str().unwrap_or("") {
            "-v" | "--verbose" => verbose = true,
            "-e" | "--encodings" => encodings = true,
            "-o" => output = it.next().map(std::path::PathBuf::from),
            "--help" => {
                let mut e = std::io::stderr().lock();
                let _ = std::io::Write::write_all(
                    &mut e,
                    b"Usage: samtools cram-size [-ve] [-o out] <in.cram>\n",
                );
                return ExitCode::SUCCESS;
            }
            s if s.starts_with('-') && s != "-" => {
                crate::diagnostics::print_error("cram-size", format!("unknown option {s}"));
                return ExitCode::from(1);
            }
            _ => {
                if input.is_none() {
                    input = Some(std::path::PathBuf::from(a));
                }
            }
        }
    }
    let Some(input) = input else {
        crate::diagnostics::print_error("cram-size", "an input CRAM file is required");
        return ExitCode::from(1);
    };
    // `-e`: a "Container encodings" block is emitted per container
    // (during the walk), then the normal report follows.
    let mut enc_text = String::new();

    let file_size = match std::fs::metadata(&input) {
        Ok(m) => m.len() as i64,
        Err(e) => {
            crate::diagnostics::print_error_errno("cram-size", "stat input", &e);
            return ExitCode::from(1);
        }
    };

    let reader = match std::fs::File::open(&input) {
        Ok(f) => f,
        Err(e) => {
            crate::diagnostics::print_error_errno("cram-size", "open input", &e);
            return ExitCode::from(1);
        }
    };
    let mut reader = htslib_rs::cram::io::Reader::new(std::io::BufReader::new(reader));
    if let Err(e) = reader.read_header() {
        crate::diagnostics::print_error_errno("cram-size", "read CRAM header", &e);
        return ExitCode::from(1);
    }

    // cu[content_id][comp] = (csize, usize). content_id -1 == CORE.
    let mut cu: HashMap<i32, Vec<(i64, i64)>> = HashMap::new();
    let mut cid2ds: HashMap<i32, Vec<i64>> = HashMap::new();
    let (mut ncont, mut nslice, mut nseqs, mut nbases) = (0i64, 0i64, 0i64, 0i64);
    let mut ref_seq_blk: i32 = -1;

    let coarse = |m: CompressionMethod| match m {
        CompressionMethod::None => method::Coarse::Raw,
        CompressionMethod::Gzip => method::Coarse::Gzip,
        CompressionMethod::Bzip2 => method::Coarse::Bzip2,
        CompressionMethod::Lzma => method::Coarse::Lzma,
        CompressionMethod::Rans4x8 => method::Coarse::Rans4x8,
        CompressionMethod::RansNx16 => method::Coarse::RansNx16,
        CompressionMethod::AdaptiveArithmeticCoding => method::Coarse::Arith,
        CompressionMethod::Fqzcomp => method::Coarse::Fqzcomp,
        CompressionMethod::NameTokenizer => method::Coarse::Tok3,
        _ => method::Coarse::Raw,
    };

    let mut container = htslib_rs::cram::io::reader::Container::default();
    loop {
        match reader.read_container(&mut container) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => {
                crate::diagnostics::print_error_errno("cram-size", "read container", &e);
                return ExitCode::from(1);
            }
        }
        let landmarks = container.header().landmarks();
        if landmarks.is_empty() {
            continue; // empty / EOF container
        }
        ncont += 1;
        nslice += landmarks.len() as i64;
        nseqs += container.header().record_count() as i64;
        nbases += container.header().base_count() as i64;

        if let Ok(ch) = container.compression_header() {
            cid2ds::build(&mut cid2ds, &ch);
            if encodings {
                enc_text.push_str(&describe::container_encodings(&ch));
                enc_text.push('\n');
            }
        }
        for slice in container.slices() {
            let Ok(slice) = slice else { continue };
            if let Some(id) = slice.header().embedded_reference_bases_block_content_id()
                && id >= 0
            {
                ref_seq_blk = id;
            }
        }
        let blocks = match container.blocks() {
            Ok(b) => b,
            Err(e) => {
                crate::diagnostics::print_error_errno("cram-size", "read blocks", &e);
                return ExitCode::from(1);
            }
        };
        for blk in blocks {
            if !matches!(
                blk.content_type,
                ContentType::CoreData | ContentType::ExternalData
            ) {
                continue;
            }
            let cid = if blk.content_type == ContentType::CoreData {
                -1
            } else {
                blk.content_id
            };
            let comp = method::expand(coarse(blk.compression_method), blk.src);
            let entry = cu
                .entry(cid)
                .or_insert_with(|| vec![(0, 0); method::COMP_MAX]);
            entry[comp].0 += blk.src.len() as i64;
            entry[comp].1 += blk.uncompressed_size as i64;
        }
    }

    // Report. With `-e`, the per-container "Container encodings"
    // blocks precede the normal block table + summary.
    let mut out = String::new();
    out.push_str(&enc_text);
    out.push_str(&format!(
        "#   Content_ID  Uncomp.size    Comp.size   Ratio Method{}  Data_series\n",
        if verbose { "    " } else { "" }
    ));
    let mut cids: Vec<i32> = cu.keys().copied().collect();
    cids.sort_unstable();
    let mut tot_size: i64 = 0;
    for cid in cids {
        let per = &cu[&cid];
        // Indices sorted by descending csize, tie by ascending index.
        let mut idx: Vec<usize> = (0..method::COMP_MAX).collect();
        idx.sort_by(|&a, &b| per[b].0.cmp(&per[a].0).then(a.cmp(&b)));
        let ds = cid2ds.get(&cid);
        // Data-series codes are printed on every method line; the
        // ` embedded_ref` marker is printed once, after the loop.
        let ds_codes = |s: &mut String| {
            if let Some(list) = ds {
                for &d in list {
                    s.push(' ');
                    s.push_str(&cid2ds::code_str(d));
                }
            }
        };
        let embedded = cid >= 0 && cid == ref_seq_blk;
        let cid_field = if cid < 0 {
            format!("BLOCK {:>8}", "CORE")
        } else {
            format!("BLOCK {cid:>8}")
        };
        if verbose {
            let mut first = true;
            for (c, &comp) in idx.iter().enumerate() {
                if per[comp].0 == 0 && c != 0 {
                    break;
                }
                if !first {
                    out.push('\n');
                }
                first = false;
                let (cs, us) = per[comp];
                out.push_str(&cid_field);
                out.push_str(&format!(" {us:>12} {cs:>12}"));
                let f = 100.0 * (cs as f64 + 0.0001) / (us as f64 + 0.0001);
                if f > 999.0 {
                    out.push_str(&format!("   >999% {:<11}", method::name(comp)));
                } else {
                    out.push_str(&format!(" {f:6.2}% {:<11}", method::name(comp)));
                }
                ds_codes(&mut out);
            }
            if embedded {
                out.push_str(" embedded_ref");
            }
            out.push('\n');
        } else {
            let cs: i64 = per.iter().map(|x| x.0).sum();
            let us: i64 = per.iter().map(|x| x.1).sum();
            let mut cstr = String::new();
            for &comp in &idx {
                if per[comp].0 == 0 {
                    break;
                }
                cstr.push(method::flag_char(comp));
            }
            if cstr.is_empty() {
                cstr.push('.');
            }
            out.push_str(&cid_field);
            out.push_str(&format!(" {us:>12} {cs:>12}"));
            let f = 100.0 * (cs as f64 + 0.0001) / (us as f64 + 0.0001);
            if f > 999.0 {
                out.push_str(&format!("   >999% {cstr:<7}"));
            } else {
                out.push_str(&format!(" {f:6.2}% {cstr:<7}"));
            }
            ds_codes(&mut out);
            if embedded {
                out.push_str(" embedded_ref");
            }
            out.push('\n');
        }
        tot_size += per.iter().map(|x| x.0).sum::<i64>();
    }
    out.push('\n');
    out.push_str(&format!("Number of containers  {ncont:>18}\n"));
    out.push_str(&format!("Number of slices      {nslice:>18}\n"));
    out.push_str(&format!("Number of sequences   {nseqs:>18}\n"));
    out.push_str(&format!("Number of bases       {nbases:>18}\n"));
    out.push_str(&format!("Total file size       {file_size:>18}\n"));
    out.push_str(&format!(
        "Format overhead size  {:>18}\n",
        file_size - tot_size
    ));

    let mut w: Box<dyn std::io::Write> = match output {
        Some(p) => match std::fs::File::create(&p) {
            Ok(f) => Box::new(f),
            Err(e) => {
                crate::diagnostics::print_error_errno("cram-size", "open -o", &e);
                return ExitCode::from(1);
            }
        },
        None => Box::new(std::io::stdout().lock()),
    };
    if let Err(e) = std::io::Write::write_all(&mut w, out.as_bytes())
        && e.kind() != std::io::ErrorKind::BrokenPipe
    {
        crate::diagnostics::print_error_errno("cram-size", "write output", &e);
        return ExitCode::from(1);
    }
    ExitCode::SUCCESS
}
