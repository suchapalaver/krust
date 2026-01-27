//! Input source abstraction for file and stdin.
//!
//! This module provides the [`Input`] enum for abstracting over different input sources,
//! enabling seamless Unix pipeline integration.
//!
//! # Example
//!
//! ```rust
//! use kmerust::input::Input;
//! use std::path::Path;
//!
//! // From a file path
//! let input = Input::from_path(Path::new("sequences.fa"));
//! assert!(matches!(input, Input::File(_)));
//!
//! // From stdin marker
//! let input = Input::from_path(Path::new("-"));
//! assert!(matches!(input, Input::Stdin));
//! ```

use std::path::{Path, PathBuf};

/// Input source for k-mer counting.
///
/// Represents either a file path or standard input, allowing the same
/// counting logic to work with both input sources.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum Input {
    /// Read from a file at the specified path.
    File(PathBuf),
    /// Read from standard input.
    #[default]
    Stdin,
}

impl Input {
    /// Creates an `Input` from a path.
    ///
    /// If the path is "-", returns [`Self::Stdin`].
    /// Otherwise, returns [`Self::File`] with the given path.
    ///
    /// # Example
    ///
    /// ```rust
    /// use kmerust::input::Input;
    /// use std::path::Path;
    ///
    /// let stdin = Input::from_path(Path::new("-"));
    /// assert!(stdin.is_stdin());
    ///
    /// let file = Input::from_path(Path::new("genome.fa"));
    /// assert!(file.is_file());
    /// ```
    #[must_use]
    pub fn from_path(path: &Path) -> Self {
        if path.as_os_str() == "-" {
            Self::Stdin
        } else {
            Self::File(path.to_path_buf())
        }
    }

    /// Creates an `Input` from an optional path.
    ///
    /// If `None` or "-", returns [`Self::Stdin`].
    /// Otherwise, returns [`Self::File`] with the given path.
    #[must_use]
    pub fn from_option(path: Option<&Path>) -> Self {
        path.map_or(Self::Stdin, Self::from_path)
    }

    /// Returns `true` if this input is stdin.
    #[must_use]
    pub const fn is_stdin(&self) -> bool {
        matches!(self, Self::Stdin)
    }

    /// Returns `true` if this input is a file.
    #[must_use]
    pub const fn is_file(&self) -> bool {
        matches!(self, Self::File(_))
    }

    /// Returns the file path if this is a file input.
    #[must_use]
    pub fn as_path(&self) -> Option<&Path> {
        match self {
            Self::File(path) => Some(path),
            Self::Stdin => None,
        }
    }
}

impl std::fmt::Display for Input {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::File(path) => write!(f, "{}", path.display()),
            Self::Stdin => write!(f, "<stdin>"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_path_stdin() {
        let input = Input::from_path(Path::new("-"));
        assert!(input.is_stdin());
        assert!(!input.is_file());
        assert!(input.as_path().is_none());
    }

    #[test]
    fn from_path_file() {
        let input = Input::from_path(Path::new("test.fa"));
        assert!(input.is_file());
        assert!(!input.is_stdin());
        assert_eq!(input.as_path(), Some(Path::new("test.fa")));
    }

    #[test]
    fn from_option_none() {
        let input = Input::from_option(None);
        assert!(input.is_stdin());
    }

    #[test]
    fn from_option_some_stdin() {
        let input = Input::from_option(Some(Path::new("-")));
        assert!(input.is_stdin());
    }

    #[test]
    fn from_option_some_file() {
        let input = Input::from_option(Some(Path::new("test.fa")));
        assert!(input.is_file());
    }

    #[test]
    fn display_stdin() {
        let input = Input::Stdin;
        assert_eq!(input.to_string(), "<stdin>");
    }

    #[test]
    fn display_file() {
        let input = Input::File(PathBuf::from("genome.fa"));
        assert_eq!(input.to_string(), "genome.fa");
    }

    #[test]
    fn default_is_stdin() {
        let input = Input::default();
        assert!(input.is_stdin());
    }
}
