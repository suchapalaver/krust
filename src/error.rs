//! Error types for kmerust.
//!
//! This module provides exhaustive, strongly-typed errors for all operations
//! in the library, enabling precise error handling and informative messages.

use std::path::PathBuf;
use thiserror::Error;

/// Errors that can occur in kmerust operations.
#[derive(Debug, Error)]
pub enum KmeRustError {
    /// K-mer length is outside the valid range (1-32).
    #[error("invalid k-mer length {k}: must be between {min} and {max}")]
    InvalidKmerLength { k: usize, min: u8, max: u8 },

    /// Encountered an invalid DNA base.
    #[error("invalid base '{base}' at position {position}")]
    InvalidBase { base: u8, position: usize },

    /// Failed to read sequence file.
    #[error("failed to read sequence file '{path}': {source}")]
    SequenceRead {
        #[source]
        source: std::io::Error,
        path: PathBuf,
    },

    /// Failed to parse sequence record.
    #[error("failed to parse sequence record: {details}")]
    SequenceParse { details: String },

    /// Failed to write output.
    #[error("failed to write output: {source}")]
    WriteError {
        #[source]
        source: std::io::Error,
    },

    /// Failed to serialize JSON output.
    #[error("failed to serialize JSON: {source}")]
    JsonError {
        #[source]
        source: serde_json::Error,
    },

    /// Failed to decompress gzip file.
    #[cfg(feature = "gzip")]
    #[error("failed to decompress gzip file '{path}': {source}")]
    GzipError {
        #[source]
        source: std::io::Error,
        path: PathBuf,
    },

    /// Failed to memory-map file.
    #[cfg(feature = "mmap")]
    #[error("failed to memory-map file '{path}': {source}")]
    MmapError {
        #[source]
        source: std::io::Error,
        path: PathBuf,
    },

    /// Failed to read index file.
    #[error("failed to read index file '{path}': {source}")]
    IndexRead {
        #[source]
        source: std::io::Error,
        path: PathBuf,
    },

    /// Failed to write index file.
    #[error("failed to write index file '{path}': {source}")]
    IndexWrite {
        #[source]
        source: std::io::Error,
        path: PathBuf,
    },

    /// Invalid or corrupted index file.
    #[error("invalid index file '{path}': {details}")]
    InvalidIndex { details: String, path: PathBuf },
}

/// Error for invalid k-mer length.
#[derive(Debug, Clone, Error, PartialEq, Eq)]
#[error("k-mer length {k} is out of range: must be between {min} and {max}")]
pub struct KmerLengthError {
    /// The invalid k value that was provided.
    pub k: usize,
    /// Minimum valid k-mer length.
    pub min: u8,
    /// Maximum valid k-mer length.
    pub max: u8,
}

/// Error for invalid DNA base.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InvalidBaseError {
    /// The invalid byte value.
    pub base: u8,
    /// Position of the invalid byte in the sequence.
    pub position: usize,
}

impl std::fmt::Display for InvalidBaseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.base.is_ascii_graphic() || self.base == b' ' {
            write!(
                f,
                "invalid base '{}' (0x{:02x}) at position {}",
                self.base as char, self.base, self.position
            )
        } else {
            write!(
                f,
                "invalid base 0x{:02x} at position {}",
                self.base, self.position
            )
        }
    }
}

impl std::error::Error for InvalidBaseError {}

impl From<std::io::Error> for KmeRustError {
    fn from(source: std::io::Error) -> Self {
        KmeRustError::WriteError { source }
    }
}

impl From<serde_json::Error> for KmeRustError {
    fn from(source: serde_json::Error) -> Self {
        KmeRustError::JsonError { source }
    }
}

impl From<KmerLengthError> for KmeRustError {
    fn from(err: KmerLengthError) -> Self {
        KmeRustError::InvalidKmerLength {
            k: err.k,
            min: err.min,
            max: err.max,
        }
    }
}

impl From<InvalidBaseError> for KmeRustError {
    fn from(err: InvalidBaseError) -> Self {
        KmeRustError::InvalidBase {
            base: err.base,
            position: err.position,
        }
    }
}

/// Errors that can occur when using the builder API.
#[derive(Debug, Error)]
pub enum BuilderError {
    /// K-mer length was not set before calling a counting method.
    #[error("k-mer length not set; call .k() first")]
    KmerLengthNotSet,

    /// Invalid k-mer length provided.
    #[error(transparent)]
    KmerLength(#[from] KmerLengthError),

    /// Error reading or parsing input file.
    #[error(transparent)]
    Kmerust(#[from] KmeRustError),

    /// I/O error during file operations.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// JSON serialization error.
    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),

    /// Error from underlying processing operations.
    #[error("{0}")]
    Process(String),
}

impl From<Box<dyn std::error::Error>> for BuilderError {
    fn from(err: Box<dyn std::error::Error>) -> Self {
        BuilderError::Process(err.to_string())
    }
}

impl From<crate::run::ProcessError> for BuilderError {
    fn from(err: crate::run::ProcessError) -> Self {
        BuilderError::Process(err.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_length_error_display() {
        let err = KmerLengthError {
            k: 50,
            min: 1,
            max: 32,
        };
        assert_eq!(
            err.to_string(),
            "k-mer length 50 is out of range: must be between 1 and 32"
        );
    }

    #[test]
    fn invalid_base_error_display() {
        let err = InvalidBaseError {
            base: b'N',
            position: 5,
        };
        assert_eq!(err.to_string(), "invalid base 'N' (0x4e) at position 5");
    }

    #[test]
    fn kmerust_error_from_kmer_length_error() {
        let err: KmeRustError = KmerLengthError {
            k: 0,
            min: 1,
            max: 32,
        }
        .into();
        assert!(matches!(err, KmeRustError::InvalidKmerLength { k: 0, .. }));
    }

    #[test]
    fn kmerust_error_from_invalid_base_error() {
        let err: KmeRustError = InvalidBaseError {
            base: b'X',
            position: 3,
        }
        .into();
        assert!(matches!(
            err,
            KmeRustError::InvalidBase {
                base: b'X',
                position: 3
            }
        ));
    }
}
