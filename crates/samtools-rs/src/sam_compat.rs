//! HTSlib-compatibility shims for SAM text that noodles' strict
//! `RecordBuf` reader does not accept but HTSlib does.

use std::io::{self, BufReader, Cursor};
use std::path::Path;

use htslib_rs::sam::{self, alignment::RecordBuf};

/// HTSlib's SAM parser accepts `c/C/s/S/I` as scalar integer-type
/// synonyms for aux fields — only `i` is in the SAM spec; the others are
/// BAM binary subtypes that HTSlib tolerates in text (e.g. upstream
/// fixtures use `AS:I:50`). noodles' `RecordBuf` reader rejects them, so
/// rewrite a scalar `TAG:[cCsSI]:` → `TAG:i:`. B-array subtypes
/// (`TAG:B:S,...`) are left untouched.
pub fn normalize_sam_aux_int_types(text: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(text.len());
    for line in text.split_inclusive(|&b| b == b'\n') {
        if line.first() == Some(&b'@') {
            out.extend_from_slice(line);
            continue;
        }
        let mut col = 0usize;
        let mut field_start = 0usize;
        let mut buf = line.to_vec();
        for i in 0..=buf.len() {
            let at_sep = i == buf.len() || buf[i] == b'\t' || buf[i] == b'\n';
            if at_sep {
                if col >= 11 {
                    let f = field_start;
                    if i >= f + 5
                        && buf[f + 2] == b':'
                        && buf[f + 4] == b':'
                        && matches!(buf[f + 3], b'c' | b'C' | b's' | b'S' | b'I')
                    {
                        buf[f + 3] = b'i';
                    }
                }
                col += 1;
                field_start = i + 1;
            }
        }
        out.extend_from_slice(&buf);
    }
    out
}

/// Reads a SAM file into `(Header, records)` tolerantly (applies
/// [`normalize_sam_aux_int_types`] first), mirroring HTSlib leniency.
pub fn read_sam_records_tolerant(input: &Path) -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let raw = std::fs::read(input)?;
    let normalized = normalize_sam_aux_int_types(&raw);
    let mut reader = sam::io::Reader::new(BufReader::new(Cursor::new(normalized)));
    let header = reader.read_header()?;
    let mut records = Vec::new();
    loop {
        let mut record = RecordBuf::default();
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        records.push(record);
    }
    Ok((header, records))
}

/// Opens a tolerant SAM reader over a file (normalized in memory).
pub fn open_sam_reader_tolerant(
    input: &Path,
) -> io::Result<sam::io::Reader<BufReader<Cursor<Vec<u8>>>>> {
    let raw = std::fs::read(input)?;
    let normalized = normalize_sam_aux_int_types(&raw);
    Ok(sam::io::Reader::new(BufReader::new(Cursor::new(
        normalized,
    ))))
}
