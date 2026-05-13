//! BAM record flag bit constants and parse/format helpers.
//!
//! Mirrors HTSlib's `BAM_F*` constants from `htslib/sam.h` plus
//! `bam_str2flag` / `bam_flag2str` from `sam.c`. Matching output is
//! required for the `samtools flags` subcommand parity.

/// Template has multiple segments in sequencing.
pub const BAM_FPAIRED: u32 = 0x1;
/// Each segment properly aligned according to the aligner.
pub const BAM_FPROPER_PAIR: u32 = 0x2;
/// Segment unmapped.
pub const BAM_FUNMAP: u32 = 0x4;
/// Next segment in template unmapped.
pub const BAM_FMUNMAP: u32 = 0x8;
/// SEQ is reverse complemented.
pub const BAM_FREVERSE: u32 = 0x10;
/// SEQ of next segment in template is reverse complemented.
pub const BAM_FMREVERSE: u32 = 0x20;
/// First segment in the template.
pub const BAM_FREAD1: u32 = 0x40;
/// Last segment in the template.
pub const BAM_FREAD2: u32 = 0x80;
/// Secondary alignment.
pub const BAM_FSECONDARY: u32 = 0x100;
/// Not passing quality controls.
pub const BAM_FQCFAIL: u32 = 0x200;
/// PCR or optical duplicate.
pub const BAM_FDUP: u32 = 0x400;
/// Supplementary alignment.
pub const BAM_FSUPPLEMENTARY: u32 = 0x800;

/// All flag bit/name pairs, in the canonical HTSlib order used by the
/// `samtools flags` usage banner and `bam_flag2str` output.
pub const FLAG_NAMES: &[(u32, &str)] = &[
    (BAM_FPAIRED, "PAIRED"),
    (BAM_FPROPER_PAIR, "PROPER_PAIR"),
    (BAM_FUNMAP, "UNMAP"),
    (BAM_FMUNMAP, "MUNMAP"),
    (BAM_FREVERSE, "REVERSE"),
    (BAM_FMREVERSE, "MREVERSE"),
    (BAM_FREAD1, "READ1"),
    (BAM_FREAD2, "READ2"),
    (BAM_FSECONDARY, "SECONDARY"),
    (BAM_FQCFAIL, "QCFAIL"),
    (BAM_FDUP, "DUP"),
    (BAM_FSUPPLEMENTARY, "SUPPLEMENTARY"),
];

/// Build a comma-separated string of flag names set in `flag`, in the
/// HTSlib canonical order. Matches `bam_flag2str` in `sam.c`.
pub fn flag_to_str(flag: u32) -> String {
    let mut out = String::new();
    for (bit, name) in FLAG_NAMES {
        if flag & *bit != 0 {
            if !out.is_empty() {
                out.push(',');
            }
            out.push_str(name);
        }
    }
    out
}

/// Parse a flag specification — either a numeric literal (decimal,
/// `0x`-prefixed hex, or `0`-prefixed octal) or a comma-separated list of
/// flag names (case-insensitive). Matches `bam_str2flag` in `sam.c`.
///
/// Returns `None` if the string cannot be parsed.
pub fn str_to_flag(s: &str) -> Option<i32> {
    // HTSlib first tries strtol(str, &end, 0). It succeeds iff at least
    // one character was consumed. We replicate by trying parsing as int
    // with auto-detected base, but only treating the *whole* string as a
    // numeric literal (mirroring how upstream lets strtol consume just
    // the prefix would let `1abc` parse as `1`, which we preserve via
    // `parse_auto_base_prefix`).
    if let Some(n) = parse_auto_base_prefix(s) {
        return Some(n);
    }

    let mut flag: u32 = 0;
    for piece in s.split(',') {
        let mut matched = false;
        for (bit, name) in FLAG_NAMES {
            if piece.eq_ignore_ascii_case(name) {
                flag |= *bit;
                matched = true;
                break;
            }
        }
        if !matched {
            return None;
        }
    }
    Some(flag as i32)
}

/// HTSlib-compatible `strtol(.., 0)` consumer: returns a signed value if
/// `s` starts with a numeric literal, even if there is trailing garbage,
/// to match `bam_str2flag`'s "the conversion was successful when end != str"
/// branch.
fn parse_auto_base_prefix(s: &str) -> Option<i32> {
    let bytes = s.as_bytes();
    let (sign, rest) = match bytes.first() {
        Some(b'-') => (-1i64, &bytes[1..]),
        Some(b'+') => (1, &bytes[1..]),
        _ => (1, bytes),
    };
    if rest.is_empty() {
        return None;
    }

    let (base, digits) = if rest.starts_with(b"0x") || rest.starts_with(b"0X") {
        (16u32, &rest[2..])
    } else if rest.first() == Some(&b'0') && rest.len() > 1 {
        (8u32, &rest[1..])
    } else {
        (10u32, rest)
    };

    let mut acc: i64 = 0;
    let mut consumed = 0usize;
    for &b in digits {
        let d = match (base, b) {
            (10, b'0'..=b'9') => (b - b'0') as i64,
            (8, b'0'..=b'7') => (b - b'0') as i64,
            (16, b'0'..=b'9') => (b - b'0') as i64,
            (16, b'a'..=b'f') => (b - b'a' + 10) as i64,
            (16, b'A'..=b'F') => (b - b'A' + 10) as i64,
            _ => break,
        };
        acc = acc.checked_mul(base as i64)?.checked_add(d)?;
        consumed += 1;
    }
    if consumed == 0 {
        return None;
    }
    let value = sign.checked_mul(acc)?;
    i32::try_from(value).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flag_to_str_basic() {
        assert_eq!(flag_to_str(0), "");
        assert_eq!(flag_to_str(BAM_FPAIRED), "PAIRED");
        assert_eq!(flag_to_str(BAM_FPAIRED | BAM_FREAD1), "PAIRED,READ1");
        assert_eq!(
            flag_to_str(0xFFF),
            "PAIRED,PROPER_PAIR,UNMAP,MUNMAP,REVERSE,MREVERSE,READ1,READ2,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY"
        );
    }

    #[test]
    fn str_to_flag_numeric() {
        assert_eq!(str_to_flag("0"), Some(0));
        assert_eq!(str_to_flag("12"), Some(12));
        assert_eq!(str_to_flag("0xff"), Some(0xff));
        assert_eq!(str_to_flag("0XFF"), Some(0xff));
        assert_eq!(str_to_flag("020"), Some(0o20));
    }

    #[test]
    fn str_to_flag_names() {
        assert_eq!(str_to_flag("PAIRED"), Some(BAM_FPAIRED as i32));
        assert_eq!(str_to_flag("paired"), Some(BAM_FPAIRED as i32));
        assert_eq!(
            str_to_flag("PAIRED,READ1"),
            Some((BAM_FPAIRED | BAM_FREAD1) as i32)
        );
    }

    #[test]
    fn str_to_flag_invalid() {
        assert_eq!(str_to_flag("NOPE"), None);
        assert_eq!(str_to_flag("PAIRED,NOPE"), None);
    }
}
