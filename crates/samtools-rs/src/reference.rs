//! Shared FASTA reference helpers.
//!
//! These helpers cover the reference-index pieces that multiple samtools
//! commands need before the full reference-backed algorithms land.

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

/// A sequence dictionary entry represented as `(name, length)`.
pub type ReferenceDict = Vec<(String, u64)>;

/// Returns the conventional associated `.fai` path for a FASTA.
pub fn fai_path_for(fasta: &Path) -> PathBuf {
    let mut p = fasta.as_os_str().to_owned();
    p.push(".fai");
    PathBuf::from(p)
}

/// Ensures a FASTA index exists and returns its path.
pub fn ensure_fai_index(fasta: &Path, fai_path: Option<&Path>) -> io::Result<PathBuf> {
    let path = fai_path
        .map(Path::to_path_buf)
        .unwrap_or_else(|| fai_path_for(fasta));
    if path.exists() {
        return Ok(path);
    }

    let file = File::open(fasta)?;
    let index = htslib_rs::faidx_compat::build_index(BufReader::new(file))?;
    let out = File::create(&path)?;
    htslib_rs::faidx_compat::write_index(out, &index)?;
    Ok(path)
}

/// Loads a FASTA dictionary from an existing `.fai` index.
pub fn read_fai_dict(fai_path: &Path) -> io::Result<ReferenceDict> {
    let file = File::open(fai_path)?;
    read_fai_dict_from_reader(BufReader::new(file))
}

/// Ensures `<fasta>.fai` exists, then loads `(SN, LN)` entries from it.
pub fn load_fai_dict(fasta: &Path) -> io::Result<ReferenceDict> {
    let fai_path = ensure_fai_index(fasta, None)?;
    read_fai_dict(&fai_path)
}

/// Returns the FASTA path whose `.fai` dictionary exactly matches `sq_dict`.
pub fn matching_reference<P>(fa_paths: &[P], sq_dict: &[(String, u64)]) -> Option<PathBuf>
where
    P: AsRef<Path>,
{
    if sq_dict.is_empty() {
        return None;
    }

    for fa in fa_paths {
        let path = fa.as_ref();
        let dict = load_fai_dict(path).ok()?;
        if dict.len() == sq_dict.len()
            && dict
                .iter()
                .zip(sq_dict.iter())
                .all(|((n1, l1), (n2, l2))| n1 == n2 && l1 == l2)
        {
            return Some(path.to_path_buf());
        }
    }
    None
}

fn read_fai_dict_from_reader<R>(reader: R) -> io::Result<ReferenceDict>
where
    R: BufRead,
{
    let mut dict = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let mut fields = line.split('\t');
        let name = fields.next().unwrap_or("").to_string();
        let length: u64 = fields.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        if !name.is_empty() && length > 0 {
            dict.push((name, length));
        }
    }
    Ok(dict)
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io::Cursor;
    use std::path::Path;
    use std::process;

    use super::{fai_path_for, load_fai_dict, matching_reference, read_fai_dict_from_reader};

    #[test]
    fn derives_associated_fai_path() {
        assert_eq!(fai_path_for(Path::new("ref.fa")), Path::new("ref.fa.fai"));
    }

    #[test]
    fn reads_name_and_length_from_fai() {
        let data = b"chr1\t4\t6\t4\t5\nchr2\t2\t17\t2\t3\n";

        let dict = read_fai_dict_from_reader(Cursor::new(&data[..])).unwrap();

        assert_eq!(dict, vec![("chr1".into(), 4), ("chr2".into(), 2)]);
    }

    #[test]
    fn builds_index_and_matches_reference_dictionary() {
        let fasta = std::env::temp_dir().join(format!(
            "samtools-rs-reference-{}-{}.fa",
            process::id(),
            "dict"
        ));
        let fai = fai_path_for(&fasta);
        let _ = fs::remove_file(&fasta);
        let _ = fs::remove_file(&fai);

        fs::write(&fasta, b">chr1\nACGT\n>chr2\nAA\n").unwrap();

        let dict = load_fai_dict(&fasta).unwrap();
        assert_eq!(dict, vec![("chr1".into(), 4), ("chr2".into(), 2)]);

        let matched = matching_reference(std::slice::from_ref(&fasta), &dict).unwrap();
        assert_eq!(matched, fasta);

        let _ = fs::remove_file(&fasta);
        let _ = fs::remove_file(&fai);
    }
}
