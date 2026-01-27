//! Memory-mapped file support for efficient I/O.
//!
//! This module provides memory-mapped file reading for large FASTA files,
//! which can be more efficient than regular I/O for certain access patterns.
//!
//! # Example
//!
//! ```rust,no_run
//! use kmerust::mmap::MmapFasta;
//!
//! let mmap = MmapFasta::open("genome.fa")?;
//! let data = mmap.as_bytes();
//! // Process the memory-mapped data...
//! # Ok::<(), std::io::Error>(())
//! ```
//!
//! # Safety
//!
//! Memory mapping relies on the underlying file not being modified while
//! the mapping is active. Modifying a mapped file leads to undefined behavior.

use memmap2::Mmap;
use std::{fs::File, io, path::Path};

/// A memory-mapped FASTA file.
///
/// Provides zero-copy access to file contents through memory mapping.
/// The file is mapped read-only and remains mapped for the lifetime of this struct.
pub struct MmapFasta {
    mmap: Mmap,
}

impl MmapFasta {
    /// Open and memory-map a FASTA file.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the FASTA file to map
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or mapped.
    ///
    /// # Safety
    ///
    /// The underlying file must not be modified while this mapping exists.
    /// Modifying a mapped file leads to undefined behavior.
    #[allow(unsafe_code)]
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        // SAFETY: We rely on the file not being modified while mapped.
        // This is documented behavior that callers must ensure.
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(Self { mmap })
    }

    /// Get a slice of the mapped file contents.
    pub fn as_bytes(&self) -> &[u8] {
        &self.mmap
    }

    /// Get the length of the mapped file in bytes.
    pub fn len(&self) -> usize {
        self.mmap.len()
    }

    /// Check if the mapped file is empty.
    pub fn is_empty(&self) -> bool {
        self.mmap.is_empty()
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn mmap_fasta_open_and_read() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, ">seq1").unwrap();
        writeln!(temp, "ACGT").unwrap();
        temp.flush().unwrap();

        let mmap = MmapFasta::open(temp.path()).unwrap();
        assert!(!mmap.is_empty());
        assert!(mmap.as_bytes().starts_with(b">seq1"));
    }

    #[test]
    fn mmap_fasta_len() {
        let mut temp = NamedTempFile::new().unwrap();
        write!(temp, "ACGT").unwrap();
        temp.flush().unwrap();

        let mmap = MmapFasta::open(temp.path()).unwrap();
        assert_eq!(mmap.len(), 4);
    }
}
