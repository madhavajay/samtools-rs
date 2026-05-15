//! Shared SAM-text rendering helpers.
//!
//! noodles' `sam::io::Writer` formats `f32` aux values as plain decimals
//! (`6.626e-34` → `0.000…0006626`, `2.9979e+09` → `2997900000`), whereas
//! htslib uses `%g`-style scientific notation with a signed two-digit
//! exponent. [`fix_sam_aux_floats`] post-processes a serialized SAM record
//! line so `:f:` scalars and `B:f,` arrays match htslib's spelling, which
//! the upstream `test.pl` fixtures expect. [`write_record`] is the drop-in
//! replacement for `sam::io::Writer::write_alignment_record` that applies
//! that fix.

use std::io::{self, Write};

use htslib_rs::sam::{self, alignment::RecordBuf};

/// Writes a SAM header to `out` (no float fields exist in `@` lines, so
/// this is a thin pass-through to noodles' header serializer; provided so
/// callers can keep a plain `Write` sink and use [`write_record`] for the
/// record stream).
pub fn write_header<W: Write>(out: &mut W, header: &sam::Header) -> io::Result<()> {
    let mut w = sam::io::Writer::new(out);
    w.write_header(header)
}

/// Writes one alignment record as a SAM text line to `out`, applying
/// [`fix_sam_aux_floats`] so float aux fields match htslib's spelling.
/// Drop-in replacement for `sam::io::Writer::write_alignment_record`
/// followed by the writer's own newline.
pub fn write_record<W: Write + ?Sized>(
    out: &mut W,
    header: &sam::Header,
    record: &RecordBuf,
) -> io::Result<()> {
    use sam::alignment::io::Write as _;

    let mut buf = Vec::new();
    {
        let mut w = sam::io::Writer::new(&mut buf);
        w.write_alignment_record(header, record)?;
    }
    // noodles terminates the record with '\n'; fix the line, re-add it.
    let line = std::str::from_utf8(&buf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
        .trim_end_matches('\n');
    out.write_all(fix_sam_aux_floats(line).as_bytes())?;
    out.write_all(b"\n")
}

/// Formats a single `f32` aux value the way htslib's `%g`/`kputd` does:
/// scientific notation with a signed, ≥2-digit exponent outside the
/// `[1e-4, 1e6)` magnitude window, plain decimal otherwise.
pub fn format_aux_float(n: f32) -> String {
    let abs = n.abs();
    if n != 0.0 && !(1e-4..1e6).contains(&abs) {
        format_htslib_exponent(n)
    } else {
        format!("{n}")
    }
}

/// Renders `n` in `<mantissa>e<+NN>` form (signed exponent, ≥2 digits),
/// matching htslib's exponent spelling.
pub fn format_htslib_exponent(n: f32) -> String {
    let raw = format!("{n:e}");
    let Some((mantissa, exponent)) = raw.split_once('e') else {
        return raw;
    };
    let value = exponent.parse::<i32>().unwrap_or(0);
    format!("{mantissa}e{value:+03}")
}

/// Applies [`fix_sam_aux_floats`] to every record line of a SAM text
/// block (multiple `\n`-terminated lines). Header lines (`@`-prefixed)
/// and empty lines pass through untouched. Use this on SAM text produced
/// from binary records (BAM/CRAM → SAM) to match htslib's float spelling.
pub fn fix_sam_text(text: &str) -> String {
    if !text.contains(":f:") && !text.contains(":B:f,") {
        return text.to_string();
    }
    let mut out = String::with_capacity(text.len());
    for piece in text.split_inclusive('\n') {
        let nl = piece.ends_with('\n');
        let body = piece.strip_suffix('\n').unwrap_or(piece);
        if body.starts_with('@') || body.is_empty() {
            out.push_str(body);
        } else {
            out.push_str(&fix_sam_aux_floats(body));
        }
        if nl {
            out.push('\n');
        }
    }
    out
}

/// Reformats float-bearing aux fields in a serialized SAM record line so
/// they match htslib's spelling. Operates only on optional fields (index
/// ≥ 11): `TAG:f:VALUE` scalars and `TAG:B:f,V1,V2,…` arrays. Mandatory
/// columns and non-float fields are passed through untouched. The input
/// must be a single record line without a trailing newline.
///
/// Float values are parsed back from the (full-precision) decimal that
/// noodles emitted; an `f32` round-trips exactly through such a string,
/// so the only change is the spelling.
pub fn fix_sam_aux_floats(line: &str) -> String {
    if !line.contains(":f:") && !line.contains(":B:f,") {
        return line.to_string();
    }

    let mut out = String::with_capacity(line.len());
    for (i, field) in line.split('\t').enumerate() {
        if i > 0 {
            out.push('\t');
        }
        if i < 11 {
            out.push_str(field);
            continue;
        }
        out.push_str(&fix_aux_field(field));
    }
    out
}

fn fix_aux_field(field: &str) -> String {
    // Scalar float: `TAG:f:VALUE`
    if let Some(value) = field.strip_prefix_aux("f") {
        let (tag_type, raw) = value;
        return match raw.parse::<f32>() {
            Ok(n) => format!("{tag_type}{}", format_aux_float(n)),
            Err(_) => field.to_string(),
        };
    }

    // Float array: `TAG:B:f,V1,V2,…`
    if field.len() >= 6 && &field[2..6] == ":B:f" {
        let prefix = &field[..6]; // `XX:B:f`
        let rest = &field[6..]; // `,V1,V2,…` (or empty)
        if let Some(values) = rest.strip_prefix(',') {
            let fixed: Vec<String> = values
                .split(',')
                .map(|v| match v.parse::<f32>() {
                    Ok(n) => format_aux_float(n),
                    Err(_) => v.to_string(),
                })
                .collect();
            return format!("{prefix},{}", fixed.join(","));
        }
    }

    field.to_string()
}

trait AuxFieldExt {
    /// For a `TAG:<ty>:VALUE` field where `<ty>` matches `ty`, returns
    /// `(("TAG:<ty>:"), "VALUE")`; otherwise `None`.
    fn strip_prefix_aux(&self, ty: &str) -> Option<(&str, &str)>;
}

impl AuxFieldExt for str {
    fn strip_prefix_aux(&self, ty: &str) -> Option<(&str, &str)> {
        // `TAG:<ty>:` is 2 + 1 + ty.len() + 1 bytes.
        let header_len = 2 + 1 + ty.len() + 1;
        if self.len() < header_len {
            return None;
        }
        let bytes = self.as_bytes();
        if bytes[2] != b':' || bytes[header_len - 1] != b':' {
            return None;
        }
        if &self[3..3 + ty.len()] != ty {
            return None;
        }
        Some((&self[..header_len], &self[header_len..]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalar_float_uses_htslib_exponent() {
        // The plain-decimal spelling noodles would emit for this f32.
        let n: f32 = 6.022e23;
        let noodles = format!("{n}");
        let line = format!("r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tfa:f:{noodles}\tNM:i:0");
        assert_eq!(
            fix_sam_aux_floats(&line),
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tfa:f:6.022e+23\tNM:i:0"
        );
    }

    #[test]
    fn plain_float_in_window_is_unchanged() {
        let line = "r1\t0\tc\t1\t60\t4M\t*\t0\t0\tA\t!\tfa:f:3.14159";
        assert_eq!(fix_sam_aux_floats(line), line);
    }

    #[test]
    fn b_float_array_each_value_fixed() {
        let got = fix_sam_aux_floats(
            "r1\t0\tc\t1\t60\t4M\t*\t0\t0\tA\t!\tbg:B:f,2.71828,0.0000000000000000000000000000000006626,2997900000",
        );
        assert_eq!(
            got,
            "r1\t0\tc\t1\t60\t4M\t*\t0\t0\tA\t!\tbg:B:f,2.71828,6.626e-34,2.9979e+09"
        );
    }

    #[test]
    fn non_float_fields_and_mandatory_columns_untouched() {
        let line = "r1\t0\tc\t1\t60\t4M\t*\t0\t0\tA\t!\tza:Z:f:not a float\tNM:i:0";
        assert_eq!(fix_sam_aux_floats(line), line);
    }

    #[test]
    fn line_without_floats_short_circuits() {
        let line = "r1\t0\tc\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tNM:i:0\tMD:Z:4";
        assert_eq!(fix_sam_aux_floats(line), line);
    }

    #[test]
    fn negative_large_exponent() {
        // -1.5e8 is outside [1e-4,1e6); htslib spells it -1.5e+08.
        assert_eq!(format_aux_float(-1.5e8), "-1.5e+08");
        assert_eq!(format_aux_float(0.0), "0");
        assert_eq!(format_aux_float(3.5), "3.5");
    }
}
