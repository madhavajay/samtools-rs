//! BAM sanitizer option parsing and record mutation.

use htslib_rs::sam::{
    self,
    alignment::{
        RecordBuf,
        record::{
            Flags, MappingQuality,
            cigar::{Op, op::Kind},
            data::field::Tag,
        },
        record_buf::Cigar,
    },
};

/// BAM sanitizer bit flags matching upstream `samtools.h`.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct SanitizeFlags(u16);

impl SanitizeFlags {
    pub const POS: Self = Self(2);
    pub const MQUAL: Self = Self(4);
    pub const UNMAP: Self = Self(8);
    pub const CIGAR: Self = Self(16);
    pub const AUX: Self = Self(32);
    pub const CIGDUP: Self = Self(64);
    pub const CIGARX: Self = Self(128);

    /// Upstream's `FIX_ON` default for position-sorted data.
    pub const ON: Self = Self(
        Self::MQUAL.bits()
            | Self::UNMAP.bits()
            | Self::CIGAR.bits()
            | Self::AUX.bits()
            | Self::CIGDUP.bits(),
    );

    /// Upstream's `FIX_ALL`; intentionally excludes `FIX_CIGARX`.
    pub const ALL: Self = Self(127);
    pub const ALL_WITH_CIGARX: Self = Self(Self::ALL.bits() | Self::CIGARX.bits());

    pub const fn empty() -> Self {
        Self(0)
    }

    pub const fn bits(self) -> u16 {
        self.0
    }

    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 == other.0
    }

    fn insert(&mut self, other: Self) {
        self.0 |= other.0;
    }
}

/// Parses upstream `bam_sanitize_options` comma-separated keywords.
pub fn parse_sanitize_options(raw: &str) -> Result<SanitizeFlags, String> {
    let mut flags = SanitizeFlags::empty();

    for keyword in raw.split(',').filter(|s| !s.is_empty()) {
        if keyword.starts_with("all") || keyword.starts_with('*') {
            flags = SanitizeFlags::ALL;
        } else if keyword.starts_with("none") || keyword.starts_with("off") {
            flags = SanitizeFlags::empty();
        } else if keyword.starts_with("on") {
            // Match bam_mate.c's parser exactly. This differs from FIX_ON in
            // samtools.h, which also includes FIX_CIGDUP.
            flags = SanitizeFlags::MQUAL;
            flags.insert(SanitizeFlags::UNMAP);
            flags.insert(SanitizeFlags::CIGAR);
            flags.insert(SanitizeFlags::AUX);
        } else if keyword.starts_with("pos") {
            flags.insert(SanitizeFlags::POS);
        } else if keyword.starts_with("mqual") {
            flags.insert(SanitizeFlags::MQUAL);
        } else if keyword.starts_with("unmap") {
            flags.insert(SanitizeFlags::UNMAP);
        } else if keyword.starts_with("cigdup") {
            flags.insert(SanitizeFlags::CIGDUP);
        } else if keyword.starts_with("cigarx") {
            flags.insert(SanitizeFlags::CIGARX);
            flags.insert(SanitizeFlags::CIGDUP);
        } else if keyword.starts_with("cigar") {
            flags.insert(SanitizeFlags::CIGAR);
        } else if keyword.starts_with("aux") {
            flags.insert(SanitizeFlags::AUX);
        } else {
            return Err(format!("unrecognised sanitize keyword \"{keyword}\""));
        }
    }

    Ok(flags)
}

/// Applies upstream-style `bam_sanitize` fixes to a mutable SAM/BAM record.
pub fn sanitize_record(header: &sam::Header, record: &mut RecordBuf, flags: SanitizeFlags) {
    if flags == SanitizeFlags::empty() {
        return;
    }

    if flags.contains(SanitizeFlags::POS) && record.reference_sequence_id().is_none() {
        *record.alignment_start_mut() = None;
        if flags.contains(SanitizeFlags::UNMAP) {
            let mut record_flags = record.flags();
            record_flags.insert(Flags::UNMAPPED);
            *record.flags_mut() = record_flags;
        }
    }

    if flags.contains(SanitizeFlags::CIGAR) && !record.flags().is_unmapped() {
        if record.alignment_start().is_none() {
            if flags.contains(SanitizeFlags::UNMAP) {
                let mut record_flags = record.flags();
                record_flags.insert(Flags::UNMAPPED);
                *record.flags_mut() = record_flags;
            }
        } else if let Some(reference_len) = reference_len(header, record.reference_sequence_id()) {
            let start0 = record
                .alignment_start()
                .map(|pos| pos.get().saturating_sub(1))
                .unwrap_or(usize::MAX);
            if start0 >= reference_len {
                if flags.contains(SanitizeFlags::UNMAP) {
                    let mut record_flags = record.flags();
                    record_flags.insert(Flags::UNMAPPED);
                    *record.flags_mut() = record_flags;
                    if flags.contains(SanitizeFlags::POS) {
                        *record.reference_sequence_id_mut() = None;
                        *record.alignment_start_mut() = None;
                    }
                }
            } else if let Some(end) = record.alignment_end()
                && end.get() > reference_len
            {
                trim_cigar_to_reference(record, reference_len);
            }
        }
    }

    if record.flags().is_unmapped() {
        if flags.contains(SanitizeFlags::CIGAR) && !record.cigar().as_ref().is_empty() {
            *record.cigar_mut() = Cigar::default();
        }

        if flags.contains(SanitizeFlags::MQUAL) {
            *record.mapping_quality_mut() = Some(MappingQuality::MIN);
        }

        if flags.contains(SanitizeFlags::AUX) {
            remove_sanitize_aux_tags(record);
        }
    }

    if flags.contains(SanitizeFlags::CIGARX) && !record.flags().is_unmapped() {
        for op in record.cigar_mut().as_mut() {
            if matches!(op.kind(), Kind::SequenceMatch | Kind::SequenceMismatch) {
                *op = Op::new(Kind::Match, op.len());
            }
        }
    }

    if flags.contains(SanitizeFlags::CIGDUP) && !record.flags().is_unmapped() {
        merge_adjacent_cigar_ops(record);
    }
}

fn reference_len(header: &sam::Header, tid: Option<usize>) -> Option<usize> {
    let tid = tid?;
    header
        .reference_sequences()
        .get_index(tid)
        .map(|(_, reference_sequence)| usize::from(reference_sequence.length()))
}

fn trim_cigar_to_reference(record: &mut RecordBuf, reference_len: usize) {
    let Some(start) = record.alignment_start() else {
        return;
    };

    let mut remaining_ref = reference_len.saturating_sub(start.get().saturating_sub(1));
    let mut trimmed = Vec::new();

    for op in record.cigar().as_ref() {
        let kind = op.kind();
        let len = op.len();

        if len == 0 {
            continue;
        }

        if !kind.consumes_reference() {
            trimmed.push(*op);
            continue;
        }

        if remaining_ref >= len {
            trimmed.push(*op);
            remaining_ref -= len;
            continue;
        }

        if remaining_ref > 0 {
            trimmed.push(Op::new(kind, remaining_ref));
        }

        if kind.consumes_read() {
            let clipped = len - remaining_ref;
            if clipped > 0 {
                trimmed.push(Op::new(Kind::SoftClip, clipped));
            }
        }
        remaining_ref = 0;
    }

    *record.cigar_mut() = trimmed.into_iter().collect();
}

fn remove_sanitize_aux_tags(record: &mut RecordBuf) {
    for tag in [*b"NM", *b"MD", *b"CG", *b"SM"] {
        record.data_mut().remove(&Tag::from(tag));
    }
}

fn merge_adjacent_cigar_ops(record: &mut RecordBuf) {
    let mut merged: Vec<Op> = Vec::new();

    for op in record.cigar().as_ref() {
        if op.is_empty() {
            continue;
        }

        if let Some(last) = merged.last_mut()
            && last.kind() == op.kind()
        {
            *last = Op::new(last.kind(), last.len() + op.len());
            continue;
        }

        merged.push(*op);
    }

    *record.cigar_mut() = merged.into_iter().collect();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_empty_as_no_flags() {
        assert_eq!(parse_sanitize_options("").unwrap(), SanitizeFlags::empty());
    }

    #[test]
    fn parses_named_flags() {
        let flags = parse_sanitize_options("pos,mqual,unmap,cigar,aux,cigdup").unwrap();

        assert!(flags.contains(SanitizeFlags::POS));
        assert!(flags.contains(SanitizeFlags::MQUAL));
        assert!(flags.contains(SanitizeFlags::UNMAP));
        assert!(flags.contains(SanitizeFlags::CIGAR));
        assert!(flags.contains(SanitizeFlags::AUX));
        assert!(flags.contains(SanitizeFlags::CIGDUP));
        assert!(!flags.contains(SanitizeFlags::CIGARX));
    }

    #[test]
    fn cigarx_implies_cigdup() {
        let flags = parse_sanitize_options("cigarx").unwrap();

        assert!(flags.contains(SanitizeFlags::CIGARX));
        assert!(flags.contains(SanitizeFlags::CIGDUP));
    }

    #[test]
    fn all_and_off_reset_accumulated_flags() {
        assert_eq!(
            parse_sanitize_options("pos,all").unwrap(),
            SanitizeFlags::ALL
        );
        assert_eq!(
            parse_sanitize_options("all,off").unwrap(),
            SanitizeFlags::empty()
        );
    }

    #[test]
    fn on_matches_upstream_parser_behavior() {
        let flags = parse_sanitize_options("on").unwrap();

        assert!(flags.contains(SanitizeFlags::MQUAL));
        assert!(flags.contains(SanitizeFlags::UNMAP));
        assert!(flags.contains(SanitizeFlags::CIGAR));
        assert!(flags.contains(SanitizeFlags::AUX));
        assert!(!flags.contains(SanitizeFlags::CIGDUP));
    }

    #[test]
    fn rejects_unknown_keyword() {
        assert_eq!(
            parse_sanitize_options("pos,nope").unwrap_err(),
            "unrecognised sanitize keyword \"nope\""
        );
    }
}
