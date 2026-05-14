//! Shared BED interval parsing and lookup helpers.
//!
//! This is a small Rust analogue of samtools' `bedidx.c`: it stores 0-based
//! half-open BED intervals grouped by reference name and can emit the
//! 1-based inclusive region strings expected by HTSlib-style APIs.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// One BED interval, using BED's native 0-based half-open coordinates.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BedInterval {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl BedInterval {
    /// Returns this interval as a `chr:start-end` 1-based inclusive region.
    pub fn to_region_string(&self) -> String {
        format!("{}:{}-{}", self.chrom, self.start + 1, self.end)
    }

    /// Returns true when this BED interval overlaps `[start, end)` on `chrom`.
    pub fn overlaps(&self, chrom: &str, start: u64, end: u64) -> bool {
        self.chrom == chrom && self.start < end && start < self.end
    }
}

/// BED intervals grouped by reference name.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BedIndex {
    intervals: BTreeMap<String, Vec<BedInterval>>,
}

impl BedIndex {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn insert(&mut self, interval: BedInterval) {
        self.intervals
            .entry(interval.chrom.clone())
            .or_default()
            .push(interval);
    }

    pub fn finalize(&mut self) {
        for intervals in self.intervals.values_mut() {
            intervals.sort_by_key(|iv| (iv.start, iv.end));
        }
    }

    pub fn is_empty(&self) -> bool {
        self.intervals.values().all(Vec::is_empty)
    }

    pub fn intervals(&self) -> impl Iterator<Item = &BedInterval> {
        self.intervals.values().flat_map(|ivs| ivs.iter())
    }

    pub fn intervals_for(&self, chrom: &str) -> &[BedInterval] {
        self.intervals.get(chrom).map(Vec::as_slice).unwrap_or(&[])
    }

    pub fn overlaps(&self, chrom: &str, start: u64, end: u64) -> bool {
        self.intervals_for(chrom)
            .iter()
            .any(|iv| iv.overlaps(chrom, start, end))
    }

    pub fn to_region_strings(&self) -> Vec<String> {
        self.intervals()
            .map(BedInterval::to_region_string)
            .collect()
    }

    pub fn to_htslib_regions(&self) -> io::Result<Vec<htslib_rs::core::Region>> {
        self.to_region_strings()
            .into_iter()
            .map(|region| {
                region.parse::<htslib_rs::core::Region>().map_err(|e| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("region \"{}\": {}", region, e),
                    )
                })
            })
            .collect()
    }
}

/// Loads BED intervals from a file.
pub fn load_bed_index(path: &Path) -> io::Result<BedIndex> {
    let file = File::open(path)?;
    read_bed_index(BufReader::new(file))
}

/// Reads BED intervals from any line-oriented reader.
pub fn read_bed_index<R>(reader: R) -> io::Result<BedIndex>
where
    R: BufRead,
{
    let mut index = BedIndex::new();
    for line in reader.lines() {
        if let Some(interval) = parse_bed_line(&line?) {
            index.insert(interval);
        }
    }
    index.finalize();
    Ok(index)
}

/// Parses one BED line. Blank lines, comments, and UCSC metadata lines return
/// `None`, matching the permissive behavior used by the current commands.
pub fn parse_bed_line(line: &str) -> Option<BedInterval> {
    let s = line.trim_end();
    if s.is_empty() || s.starts_with('#') || s.starts_with("track ") || s.starts_with("browser ") {
        return None;
    }

    let mut fields = s.split('\t');
    let chrom = fields.next().unwrap_or("");
    let start: u64 = fields.next().and_then(|t| t.parse().ok()).unwrap_or(0);
    let end: u64 = fields.next().and_then(|t| t.parse().ok()).unwrap_or(0);
    if chrom.is_empty() || end <= start {
        return None;
    }

    Some(BedInterval {
        chrom: chrom.to_string(),
        start,
        end,
    })
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::{BedInterval, parse_bed_line, read_bed_index};

    #[test]
    fn parses_bed_interval_and_region_string() {
        let interval = parse_bed_line("chr1\t9\t20\tname").unwrap();

        assert_eq!(
            interval,
            BedInterval {
                chrom: "chr1".into(),
                start: 9,
                end: 20
            }
        );
        assert_eq!(interval.to_region_string(), "chr1:10-20");
    }

    #[test]
    fn skips_comments_metadata_and_invalid_intervals() {
        assert_eq!(parse_bed_line("#chrom\tstart\tend"), None);
        assert_eq!(parse_bed_line("track name=x"), None);
        assert_eq!(parse_bed_line("browser position chr1"), None);
        assert_eq!(parse_bed_line("chr1\t20\t20"), None);
        assert_eq!(parse_bed_line("\t0\t1"), None);
    }

    #[test]
    fn indexes_sorted_intervals_and_overlap_queries() {
        let data = b"chr2\t8\t10\nchr1\t20\t30\nchr1\t5\t7\n";
        let index = read_bed_index(Cursor::new(&data[..])).unwrap();

        assert_eq!(
            index.to_region_strings(),
            vec!["chr1:6-7", "chr1:21-30", "chr2:9-10"]
        );
        assert!(index.overlaps("chr1", 6, 8));
        assert!(!index.overlaps("chr1", 7, 20));
        assert!(index.overlaps("chr1", 29, 31));
        assert!(!index.overlaps("chr3", 0, 100));
    }
}
