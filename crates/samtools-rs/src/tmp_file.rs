//! Temporary file helpers for external sort/collate work.
//!
//! Upstream samtools has a BAM-record temp spooling layer in `tmp_file.c`.
//! This module starts with the file-lifetime pieces needed by that layer:
//! collision-resistant temp creation, path ownership, and best-effort cleanup.

use std::ffi::OsString;
use std::fs::{File, OpenOptions};
use std::io;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

static NEXT_ID: AtomicU64 = AtomicU64::new(0);

/// An owned temporary path that removes the file on drop unless persisted.
#[derive(Debug)]
pub struct TempPath {
    path: Option<PathBuf>,
}

impl TempPath {
    fn new(path: PathBuf) -> Self {
        Self { path: Some(path) }
    }

    pub fn path(&self) -> &Path {
        self.path
            .as_deref()
            .expect("temporary path already persisted")
    }

    /// Stops automatic cleanup and returns the underlying path.
    pub fn persist(mut self) -> PathBuf {
        self.path.take().expect("temporary path already persisted")
    }

    /// Removes the file immediately and consumes the owner.
    pub fn close(mut self) -> io::Result<()> {
        if let Some(path) = self.path.take()
            && path.exists()
        {
            std::fs::remove_file(path)?;
        }
        Ok(())
    }
}

impl Drop for TempPath {
    fn drop(&mut self) {
        if let Some(path) = self.path.take() {
            let _ = std::fs::remove_file(path);
        }
    }
}

/// Creates a new temporary file and returns both the open file and owned path.
pub fn create_temp_file(prefix: &str, extension: Option<&str>) -> io::Result<(File, TempPath)> {
    create_temp_file_in(std::env::temp_dir(), prefix, extension)
}

/// Creates a new temporary file in `dir`.
pub fn create_temp_file_in<P>(
    dir: P,
    prefix: &str,
    extension: Option<&str>,
) -> io::Result<(File, TempPath)>
where
    P: AsRef<Path>,
{
    let dir = dir.as_ref();
    let mut last_err = None;
    for _ in 0..128 {
        let path = candidate_path(dir, prefix, extension);
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(file) => return Ok((file, TempPath::new(path))),
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => last_err = Some(e),
            Err(e) => return Err(e),
        }
    }
    Err(last_err
        .unwrap_or_else(|| io::Error::new(io::ErrorKind::AlreadyExists, "temp path exists")))
}

fn candidate_path(dir: &Path, prefix: &str, extension: Option<&str>) -> PathBuf {
    let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);

    let mut name = OsString::from(prefix);
    name.push("-");
    name.push(std::process::id().to_string());
    name.push("-");
    name.push(nanos.to_string());
    name.push("-");
    name.push(id.to_string());
    if let Some(ext) = extension {
        name.push(".");
        name.push(ext.trim_start_matches('.'));
    }
    dir.join(name)
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::create_temp_file;

    #[test]
    fn creates_unique_temp_file_with_extension() {
        let (mut file, temp) = create_temp_file("samtools-rs-tmp-test", Some("bam")).unwrap();
        writeln!(file, "data").unwrap();

        assert!(temp.path().exists());
        assert_eq!(
            temp.path().extension().and_then(|s| s.to_str()),
            Some("bam")
        );
    }

    #[test]
    fn drop_removes_temp_file() {
        let path = {
            let (_file, temp) = create_temp_file("samtools-rs-tmp-drop", None).unwrap();
            temp.path().to_path_buf()
        };

        assert!(!path.exists());
    }

    #[test]
    fn persist_keeps_temp_file() {
        let (_file, temp) = create_temp_file("samtools-rs-tmp-persist", Some(".dat")).unwrap();
        let path = temp.persist();

        assert!(path.exists());
        assert_eq!(path.extension().and_then(|s| s.to_str()), Some("dat"));
        std::fs::remove_file(path).unwrap();
    }
}
