//! K-mer index serialization and deserialization.
//!
//! This module provides functionality to save k-mer counts to disk and load them
//! back without re-counting. This is useful for large genomes where counting is
//! expensive.
//!
//! # Binary Format (Version 1)
//!
//! The index file uses a simple binary format:
//!
//! ```text
//! +--------+--------+------+--------+------------------+--------+
//! | MAGIC  | VERSION|  K   | COUNT  |      DATA        | CRC32  |
//! | 4 bytes| 1 byte |1 byte| 8 bytes| 16 bytes × COUNT | 4 bytes|
//! +--------+--------+------+--------+------------------+--------+
//!
//! MAGIC:   "KMIX" (0x4B 0x4D 0x49 0x58)
//! VERSION: Format version (currently 1)
//! K:       K-mer length (1-32)
//! COUNT:   Number of distinct k-mers (little-endian u64)
//! DATA:    Array of (packed_bits: u64, count: u64) pairs (little-endian)
//! CRC32:   CRC32 checksum of all preceding bytes (little-endian)
//! ```
//!
//! # Compression
//!
//! Index files with `.gz` extension are automatically compressed/decompressed
//! using gzip (requires `gzip` feature).
//!
//! # Example
//!
//! ```rust,no_run
//! use kmerust::index::{KmerIndex, save_index, load_index};
//! use kmerust::kmer::KmerLength;
//! use std::collections::HashMap;
//!
//! // Create an index from counts
//! let mut counts = HashMap::new();
//! counts.insert(0b00_01_10_11u64, 42u64); // ACGT -> 42
//! let index = KmerIndex::new(KmerLength::new(4).unwrap(), counts);
//!
//! // Save to disk
//! save_index(&index, "kmers.kmix").unwrap();
//!
//! // Load back
//! let loaded = load_index("kmers.kmix").unwrap();
//! assert_eq!(loaded.k(), index.k());
//! ```

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use crate::error::KmeRustError;
use crate::kmer::{unpack_to_string, KmerLength};

/// Magic bytes identifying a kmerust index file.
const MAGIC: &[u8; 4] = b"KMIX";

/// Current format version.
const VERSION: u8 = 1;

/// A k-mer index containing packed k-mer counts.
///
/// The index stores k-mers in their canonical packed form (64-bit integers)
/// along with their counts. This is more compact than storing k-mer strings.
#[derive(Debug, Clone)]
pub struct KmerIndex {
    k: KmerLength,
    counts: HashMap<u64, u64>,
}

impl KmerIndex {
    /// Creates a new k-mer index.
    ///
    /// # Arguments
    ///
    /// * `k` - The k-mer length
    /// * `counts` - Map from packed canonical k-mer to count
    #[must_use]
    pub fn new(k: KmerLength, counts: HashMap<u64, u64>) -> Self {
        Self { k, counts }
    }

    /// Returns the k-mer length.
    #[must_use]
    pub fn k(&self) -> KmerLength {
        self.k
    }

    /// Returns the number of distinct k-mers in the index.
    #[must_use]
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /// Returns true if the index contains no k-mers.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// Returns a reference to the packed counts.
    #[must_use]
    pub fn counts(&self) -> &HashMap<u64, u64> {
        &self.counts
    }

    /// Consumes the index and returns the packed counts.
    #[must_use]
    pub fn into_counts(self) -> HashMap<u64, u64> {
        self.counts
    }

    /// Looks up the count for a specific k-mer (as packed bits).
    #[must_use]
    pub fn get(&self, packed_bits: u64) -> Option<u64> {
        self.counts.get(&packed_bits).copied()
    }

    /// Converts the index to string-keyed counts.
    ///
    /// This unpacks all k-mers from their 64-bit representation back to
    /// DNA strings. Useful for interoperability with text-based formats.
    #[must_use]
    pub fn to_string_counts(&self) -> HashMap<String, u64> {
        self.counts
            .iter()
            .map(|(&packed, &count)| (unpack_to_string(packed, self.k), count))
            .collect()
    }
}

/// Saves a k-mer index to a file.
///
/// The file format is detected from the extension:
/// - `.kmix` - uncompressed binary format
/// - `.kmix.gz` - gzip-compressed binary format (requires `gzip` feature)
///
/// # Errors
///
/// Returns an error if the file cannot be created or written.
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::index::{KmerIndex, save_index};
/// use kmerust::kmer::KmerLength;
/// use std::collections::HashMap;
///
/// let index = KmerIndex::new(KmerLength::new(21).unwrap(), HashMap::new());
/// save_index(&index, "output.kmix")?;
/// # Ok::<(), kmerust::error::KmeRustError>(())
/// ```
pub fn save_index<P: AsRef<Path>>(index: &KmerIndex, path: P) -> Result<(), KmeRustError> {
    let path = path.as_ref();

    #[cfg(feature = "gzip")]
    if is_gzip_path(path) {
        let file = File::create(path).map_err(|e| KmeRustError::IndexWrite {
            source: e,
            path: path.to_path_buf(),
        })?;
        let encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        let writer = BufWriter::new(encoder);
        return write_index(index, writer, path);
    }

    let file = File::create(path).map_err(|e| KmeRustError::IndexWrite {
        source: e,
        path: path.to_path_buf(),
    })?;
    let writer = BufWriter::new(file);
    write_index(index, writer, path)
}

/// Loads a k-mer index from a file.
///
/// The file format is detected from the extension:
/// - `.kmix` - uncompressed binary format
/// - `.kmix.gz` - gzip-compressed binary format (requires `gzip` feature)
///
/// # Errors
///
/// Returns an error if:
/// - The file cannot be opened
/// - The file is not a valid k-mer index (bad magic, version, or checksum)
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::index::load_index;
///
/// let index = load_index("counts.kmix")?;
/// println!("Loaded {} k-mers of length {}", index.len(), index.k().get());
/// # Ok::<(), kmerust::error::KmeRustError>(())
/// ```
pub fn load_index<P: AsRef<Path>>(path: P) -> Result<KmerIndex, KmeRustError> {
    let path = path.as_ref();

    #[cfg(feature = "gzip")]
    if is_gzip_path(path) {
        let file = File::open(path).map_err(|e| KmeRustError::IndexRead {
            source: e,
            path: path.to_path_buf(),
        })?;
        let decoder = flate2::read::GzDecoder::new(file);
        let reader = BufReader::new(decoder);
        return read_index(reader, path);
    }

    let file = File::open(path).map_err(|e| KmeRustError::IndexRead {
        source: e,
        path: path.to_path_buf(),
    })?;
    let reader = BufReader::new(file);
    read_index(reader, path)
}

/// Writes the index to a writer, computing CRC32 as we go.
fn write_index<W: Write, P: AsRef<Path>>(
    index: &KmerIndex,
    mut writer: W,
    path: P,
) -> Result<(), KmeRustError> {
    let mut crc = Crc32Writer::new(&mut writer);

    // Write header
    crc.write_all(MAGIC).map_err(|e| KmeRustError::IndexWrite {
        source: e,
        path: path.as_ref().to_path_buf(),
    })?;
    crc.write_all(&[VERSION])
        .map_err(|e| KmeRustError::IndexWrite {
            source: e,
            path: path.as_ref().to_path_buf(),
        })?;
    crc.write_all(&[index.k.as_u8()])
        .map_err(|e| KmeRustError::IndexWrite {
            source: e,
            path: path.as_ref().to_path_buf(),
        })?;
    crc.write_all(&(index.counts.len() as u64).to_le_bytes())
        .map_err(|e| KmeRustError::IndexWrite {
            source: e,
            path: path.as_ref().to_path_buf(),
        })?;

    // Write k-mer data
    for (&packed, &count) in &index.counts {
        crc.write_all(&packed.to_le_bytes())
            .map_err(|e| KmeRustError::IndexWrite {
                source: e,
                path: path.as_ref().to_path_buf(),
            })?;
        crc.write_all(&count.to_le_bytes())
            .map_err(|e| KmeRustError::IndexWrite {
                source: e,
                path: path.as_ref().to_path_buf(),
            })?;
    }

    // Write CRC32 checksum (not included in checksum itself)
    let checksum = crc.finalize();
    writer
        .write_all(&checksum.to_le_bytes())
        .map_err(|e| KmeRustError::IndexWrite {
            source: e,
            path: path.as_ref().to_path_buf(),
        })?;

    writer.flush().map_err(|e| KmeRustError::IndexWrite {
        source: e,
        path: path.as_ref().to_path_buf(),
    })?;

    Ok(())
}

/// Reads and validates an index from a reader.
fn read_index<R: Read, P: AsRef<Path>>(reader: R, path: P) -> Result<KmerIndex, KmeRustError> {
    let path = path.as_ref();

    // Read entire file into memory for CRC verification
    let mut data = Vec::new();
    let mut reader = BufReader::new(reader);
    reader
        .read_to_end(&mut data)
        .map_err(|e| KmeRustError::IndexRead {
            source: e,
            path: path.to_path_buf(),
        })?;

    // Need at least header (14 bytes) + CRC32 (4 bytes)
    if data.len() < 18 {
        return Err(KmeRustError::InvalidIndex {
            details: "file too small".into(),
            path: path.to_path_buf(),
        });
    }

    // Check magic first (before CRC) to give better error for non-index files
    if &data[..4] != MAGIC {
        return Err(KmeRustError::InvalidIndex {
            details: "invalid magic bytes (not a kmerust index file)".into(),
            path: path.to_path_buf(),
        });
    }

    // Split data and checksum
    let (content, checksum_bytes) = data.split_at(data.len() - 4);
    let stored_checksum = u32::from_le_bytes(checksum_bytes.try_into().unwrap());

    // Verify CRC32
    let computed_checksum = crc32(content);
    if computed_checksum != stored_checksum {
        return Err(KmeRustError::InvalidIndex {
            details: format!(
                "checksum mismatch (expected {stored_checksum:#x}, got {computed_checksum:#x})"
            ),
            path: path.to_path_buf(),
        });
    }

    // Parse header (magic already verified)
    let mut cursor = &content[4..];

    // Version
    if cursor.is_empty() || cursor[0] != VERSION {
        return Err(KmeRustError::InvalidIndex {
            details: format!("unsupported version {}", cursor.first().unwrap_or(&0)),
            path: path.to_path_buf(),
        });
    }
    cursor = &cursor[1..];

    // K-mer length
    if cursor.is_empty() {
        return Err(KmeRustError::InvalidIndex {
            details: "missing k-mer length".into(),
            path: path.to_path_buf(),
        });
    }
    let k_val = cursor[0];
    let k = KmerLength::new(k_val as usize).map_err(|e| KmeRustError::InvalidIndex {
        details: format!("invalid k-mer length: {e}"),
        path: path.to_path_buf(),
    })?;
    cursor = &cursor[1..];

    // Count
    if cursor.len() < 8 {
        return Err(KmeRustError::InvalidIndex {
            details: "missing k-mer count".into(),
            path: path.to_path_buf(),
        });
    }
    let count = u64::from_le_bytes(cursor[..8].try_into().unwrap());
    cursor = &cursor[8..];

    // Validate data size
    let expected_data_size = count as usize * 16; // 8 bytes packed + 8 bytes count
    if cursor.len() != expected_data_size {
        return Err(KmeRustError::InvalidIndex {
            details: format!(
                "data size mismatch (expected {expected_data_size} bytes, got {} bytes)",
                cursor.len()
            ),
            path: path.to_path_buf(),
        });
    }

    // Read k-mer data
    let mut counts = HashMap::with_capacity(count as usize);
    for _ in 0..count {
        let packed = u64::from_le_bytes(cursor[..8].try_into().unwrap());
        let kmer_count = u64::from_le_bytes(cursor[8..16].try_into().unwrap());
        counts.insert(packed, kmer_count);
        cursor = &cursor[16..];
    }

    Ok(KmerIndex { k, counts })
}

/// CRC32 (IEEE polynomial) computation.
fn crc32(data: &[u8]) -> u32 {
    // IEEE polynomial used by gzip, PNG, etc.
    const POLYNOMIAL: u32 = 0xEDB8_8320;

    // Build lookup table
    let table: [u32; 256] = {
        let mut table = [0u32; 256];
        for (i, entry) in table.iter_mut().enumerate() {
            let mut crc = i as u32;
            for _ in 0..8 {
                if crc & 1 != 0 {
                    crc = (crc >> 1) ^ POLYNOMIAL;
                } else {
                    crc >>= 1;
                }
            }
            *entry = crc;
        }
        table
    };

    let mut crc = !0u32;
    for &byte in data {
        crc = table[((crc ^ u32::from(byte)) & 0xFF) as usize] ^ (crc >> 8);
    }
    !crc
}

/// A writer that computes CRC32 as data is written.
struct Crc32Writer<W> {
    inner: W,
    data: Vec<u8>,
}

impl<W: Write> Crc32Writer<W> {
    fn new(inner: W) -> Self {
        Self {
            inner,
            data: Vec::new(),
        }
    }

    fn finalize(self) -> u32 {
        crc32(&self.data)
    }
}

impl<W: Write> Write for Crc32Writer<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let n = self.inner.write(buf)?;
        self.data.extend_from_slice(&buf[..n]);
        Ok(n)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

/// Checks if a path has a `.gz` extension.
#[cfg(feature = "gzip")]
fn is_gzip_path(path: &Path) -> bool {
    path.extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn roundtrip_empty_index() {
        let index = KmerIndex::new(KmerLength::new(21).unwrap(), HashMap::new());
        let tmp = NamedTempFile::with_suffix(".kmix").unwrap();

        save_index(&index, tmp.path()).unwrap();
        let loaded = load_index(tmp.path()).unwrap();

        assert_eq!(loaded.k(), index.k());
        assert!(loaded.is_empty());
    }

    #[test]
    fn roundtrip_with_data() {
        let mut counts = HashMap::new();
        counts.insert(0b00_01_10_11u64, 42u64); // ACGT
        counts.insert(0b11_10_01_00u64, 17u64); // TGCA
        counts.insert(0u64, 1u64); // AAAA

        let index = KmerIndex::new(KmerLength::new(4).unwrap(), counts.clone());
        let tmp = NamedTempFile::with_suffix(".kmix").unwrap();

        save_index(&index, tmp.path()).unwrap();
        let loaded = load_index(tmp.path()).unwrap();

        assert_eq!(loaded.k(), index.k());
        assert_eq!(loaded.len(), 3);
        assert_eq!(loaded.get(0b00_01_10_11), Some(42));
        assert_eq!(loaded.get(0b11_10_01_00), Some(17));
        assert_eq!(loaded.get(0), Some(1));
    }

    #[test]
    fn roundtrip_various_k_lengths() {
        for k_val in [1, 5, 16, 21, 32] {
            let mut counts = HashMap::new();
            counts.insert(1u64, 100u64);

            let k = KmerLength::new(k_val).unwrap();
            let index = KmerIndex::new(k, counts);
            let tmp = NamedTempFile::with_suffix(".kmix").unwrap();

            save_index(&index, tmp.path()).unwrap();
            let loaded = load_index(tmp.path()).unwrap();

            assert_eq!(loaded.k().get(), k_val);
        }
    }

    #[test]
    fn invalid_magic_rejected() {
        let tmp = NamedTempFile::with_suffix(".kmix").unwrap();

        // Write garbage (must be at least 18 bytes)
        std::fs::write(tmp.path(), b"GARBAGE_DATA_HERE_").unwrap();

        let result = load_index(tmp.path());
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(
            err.to_string().contains("invalid magic"),
            "expected 'invalid magic' error, got: {err}"
        );
    }

    #[test]
    fn corrupted_checksum_rejected() {
        let mut counts = HashMap::new();
        counts.insert(1u64, 1u64);
        let index = KmerIndex::new(KmerLength::new(4).unwrap(), counts);
        let tmp = NamedTempFile::with_suffix(".kmix").unwrap();

        save_index(&index, tmp.path()).unwrap();

        // Corrupt one byte in the middle
        let mut data = std::fs::read(tmp.path()).unwrap();
        if let Some(byte) = data.get_mut(10) {
            *byte ^= 0xFF;
        }
        std::fs::write(tmp.path(), data).unwrap();

        let result = load_index(tmp.path());
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("checksum"));
    }

    #[test]
    fn file_too_small_rejected() {
        let tmp = NamedTempFile::with_suffix(".kmix").unwrap();
        std::fs::write(tmp.path(), b"KMIX").unwrap(); // Only magic, no rest

        let result = load_index(tmp.path());
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("too small"));
    }

    #[test]
    fn to_string_counts() {
        let mut counts = HashMap::new();
        counts.insert(0b00_01_10_11u64, 42u64); // ACGT

        let index = KmerIndex::new(KmerLength::new(4).unwrap(), counts);
        let string_counts = index.to_string_counts();

        assert_eq!(string_counts.len(), 1);
        assert_eq!(string_counts.get("ACGT"), Some(&42));
    }

    #[test]
    fn crc32_known_values() {
        // Test against known CRC32 values
        assert_eq!(crc32(b""), 0x0000_0000);
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
    }

    #[cfg(feature = "gzip")]
    #[test]
    fn roundtrip_gzip() {
        let mut counts = HashMap::new();
        counts.insert(0b00_01_10_11u64, 42u64);

        let index = KmerIndex::new(KmerLength::new(4).unwrap(), counts);
        let tmp = NamedTempFile::with_suffix(".kmix.gz").unwrap();

        save_index(&index, tmp.path()).unwrap();
        let loaded = load_index(tmp.path()).unwrap();

        assert_eq!(loaded.k(), index.k());
        assert_eq!(loaded.get(0b00_01_10_11), Some(42));
    }
}
