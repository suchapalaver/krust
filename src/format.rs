//! Input format detection and selection.
//!
//! This module provides types for specifying and auto-detecting input file formats
//! (FASTA or FASTQ).

use clap::ValueEnum;
use std::ffi::OsStr;
use std::path::Path;

/// Input sequence file format.
///
/// Used to specify the format of input files. When set to `Auto`, the format
/// is detected from the file extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum)]
pub enum SequenceFormat {
    /// Auto-detect format from file extension.
    ///
    /// Detection rules:
    /// - `.fq`, `.fastq`, `.fq.gz`, `.fastq.gz` -> FASTQ
    /// - `.fa`, `.fasta`, `.fna`, `.fa.gz`, `.fasta.gz`, `.fna.gz` -> FASTA
    /// - Unknown or stdin -> FASTA (default)
    #[default]
    Auto,
    /// FASTA format (`.fa`, `.fasta`, `.fna`).
    Fasta,
    /// FASTQ format (`.fq`, `.fastq`).
    Fastq,
}

impl SequenceFormat {
    /// Detects the sequence format from a file path's extension.
    ///
    /// Handles gzip-compressed files by stripping the `.gz` extension first.
    ///
    /// # Examples
    ///
    /// ```
    /// use kmerust::format::SequenceFormat;
    /// use std::path::Path;
    ///
    /// assert_eq!(SequenceFormat::from_extension(Path::new("reads.fq")), SequenceFormat::Fastq);
    /// assert_eq!(SequenceFormat::from_extension(Path::new("reads.fastq.gz")), SequenceFormat::Fastq);
    /// assert_eq!(SequenceFormat::from_extension(Path::new("genome.fa")), SequenceFormat::Fasta);
    /// assert_eq!(SequenceFormat::from_extension(Path::new("genome.fasta.gz")), SequenceFormat::Fasta);
    /// ```
    #[must_use]
    pub fn from_extension(path: &Path) -> Self {
        // Get the extension, stripping .gz if present
        let ext = path
            .extension()
            .and_then(OsStr::to_str)
            .map(|e| e.to_lowercase());

        let effective_ext = match ext.as_deref() {
            Some("gz") => {
                // Strip .gz and get the real extension
                path.file_stem()
                    .and_then(|stem| Path::new(stem).extension())
                    .and_then(OsStr::to_str)
                    .map(|e| e.to_lowercase())
            }
            other => other.map(String::from),
        };

        match effective_ext.as_deref() {
            Some("fq" | "fastq") => Self::Fastq,
            Some("fa" | "fasta" | "fna") => Self::Fasta,
            _ => Self::Fasta, // Default to FASTA for unknown extensions
        }
    }

    /// Resolves `Auto` format to a concrete format based on the file path.
    ///
    /// - If format is already `Fasta` or `Fastq`, returns it unchanged.
    /// - If format is `Auto` and a path is provided, detects from extension.
    /// - If format is `Auto` and no path is provided (stdin), defaults to `Fasta`.
    ///
    /// # Examples
    ///
    /// ```
    /// use kmerust::format::SequenceFormat;
    /// use std::path::Path;
    ///
    /// // Auto-detection from path
    /// let format = SequenceFormat::Auto.resolve(Some(Path::new("reads.fq")));
    /// assert_eq!(format, SequenceFormat::Fastq);
    ///
    /// // Explicit format is unchanged
    /// let format = SequenceFormat::Fasta.resolve(Some(Path::new("reads.fq")));
    /// assert_eq!(format, SequenceFormat::Fasta);
    ///
    /// // Stdin defaults to FASTA
    /// let format = SequenceFormat::Auto.resolve(None);
    /// assert_eq!(format, SequenceFormat::Fasta);
    /// ```
    #[must_use]
    pub fn resolve(self, path: Option<&Path>) -> Self {
        match self {
            Self::Auto => path.map_or(Self::Fasta, Self::from_extension),
            other => other,
        }
    }

    /// Returns `true` if this format is FASTQ.
    #[must_use]
    pub fn is_fastq(self) -> bool {
        matches!(self, Self::Fastq)
    }

    /// Returns `true` if this format is FASTA.
    #[must_use]
    pub fn is_fasta(self) -> bool {
        matches!(self, Self::Fasta)
    }
}

impl std::fmt::Display for SequenceFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Auto => write!(f, "auto"),
            Self::Fasta => write!(f, "fasta"),
            Self::Fastq => write!(f, "fastq"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_extension_fasta() {
        assert_eq!(
            SequenceFormat::from_extension(Path::new("test.fa")),
            SequenceFormat::Fasta
        );
        assert_eq!(
            SequenceFormat::from_extension(Path::new("test.fasta")),
            SequenceFormat::Fasta
        );
        assert_eq!(
            SequenceFormat::from_extension(Path::new("test.fna")),
            SequenceFormat::Fasta
        );
    }

    #[test]
    fn from_extension_fastq() {
        assert_eq!(
            SequenceFormat::from_extension(Path::new("test.fq")),
            SequenceFormat::Fastq
        );
        assert_eq!(
            SequenceFormat::from_extension(Path::new("test.fastq")),
            SequenceFormat::Fastq
        );
    }

    #[test]
    fn from_extension_gzipped() {
        assert_eq!(
            SequenceFormat::from_extension(Path::new("test.fa.gz")),
            SequenceFormat::Fasta
        );
        assert_eq!(
            SequenceFormat::from_extension(Path::new("test.fasta.gz")),
            SequenceFormat::Fasta
        );
        assert_eq!(
            SequenceFormat::from_extension(Path::new("test.fq.gz")),
            SequenceFormat::Fastq
        );
        assert_eq!(
            SequenceFormat::from_extension(Path::new("test.fastq.gz")),
            SequenceFormat::Fastq
        );
    }

    #[test]
    fn from_extension_unknown_defaults_to_fasta() {
        assert_eq!(
            SequenceFormat::from_extension(Path::new("test.txt")),
            SequenceFormat::Fasta
        );
        assert_eq!(
            SequenceFormat::from_extension(Path::new("test")),
            SequenceFormat::Fasta
        );
    }

    #[test]
    fn resolve_auto_with_path() {
        assert_eq!(
            SequenceFormat::Auto.resolve(Some(Path::new("test.fq"))),
            SequenceFormat::Fastq
        );
        assert_eq!(
            SequenceFormat::Auto.resolve(Some(Path::new("test.fa"))),
            SequenceFormat::Fasta
        );
    }

    #[test]
    fn resolve_auto_without_path() {
        assert_eq!(SequenceFormat::Auto.resolve(None), SequenceFormat::Fasta);
    }

    #[test]
    fn resolve_explicit_format_unchanged() {
        assert_eq!(
            SequenceFormat::Fasta.resolve(Some(Path::new("test.fq"))),
            SequenceFormat::Fasta
        );
        assert_eq!(
            SequenceFormat::Fastq.resolve(Some(Path::new("test.fa"))),
            SequenceFormat::Fastq
        );
    }

    #[test]
    fn display() {
        assert_eq!(format!("{}", SequenceFormat::Auto), "auto");
        assert_eq!(format!("{}", SequenceFormat::Fasta), "fasta");
        assert_eq!(format!("{}", SequenceFormat::Fastq), "fastq");
    }
}
