//! Shared parser for comma-separated SAM auxiliary tag lists.
//!
//! Upstream `parse_aux_list` accepts `TAG,TAG` strings where every tag is
//! exactly two bytes. Duplicates collapse into a set.

use std::collections::HashSet;

/// A two-byte SAM auxiliary tag.
pub type AuxTag = [u8; 2];

/// Parses a comma-separated list of two-byte SAM auxiliary tags.
pub fn parse_aux_list(raw: &str) -> Result<HashSet<AuxTag>, AuxListError> {
    let mut tags = HashSet::new();
    let mut rest = raw.as_bytes();

    while rest.len() >= 2 {
        tags.insert([rest[0], rest[1]]);
        rest = &rest[2..];

        if rest.first() == Some(&b',') {
            rest = &rest[1..];
        } else if !rest.is_empty() {
            break;
        }
    }

    if rest.is_empty() {
        Ok(tags)
    } else {
        Err(AuxListError)
    }
}

/// Error returned when an auxiliary tag list contains a malformed tag.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AuxListError;

impl std::fmt::Display for AuxListError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("auxiliary tags should be exactly two characters long")
    }
}

impl std::error::Error for AuxListError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_comma_separated_tags() {
        let tags = parse_aux_list("NM,MD,AS").unwrap();
        assert!(tags.contains(b"NM"));
        assert!(tags.contains(b"MD"));
        assert!(tags.contains(b"AS"));
        assert_eq!(tags.len(), 3);
    }

    #[test]
    fn collapses_duplicates() {
        let tags = parse_aux_list("NM,NM").unwrap();
        assert_eq!(tags.len(), 1);
        assert!(tags.contains(b"NM"));
    }

    #[test]
    fn rejects_malformed_tags() {
        assert_eq!(parse_aux_list("N").unwrap_err(), AuxListError);
        assert_eq!(parse_aux_list("NM,M").unwrap_err(), AuxListError);
        assert_eq!(parse_aux_list("NMD").unwrap_err(), AuxListError);
    }
}
