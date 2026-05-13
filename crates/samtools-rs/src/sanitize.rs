//! BAM sanitizer option parsing.

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
