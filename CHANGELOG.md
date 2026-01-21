# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Comprehensive test suite with unit and integration tests
- CI/CD pipeline with GitHub Actions
- Rustdoc documentation for public APIs
- MSRV (Minimum Supported Rust Version) set to 1.75

### Changed

- Updated all dependencies to latest stable versions
- Replaced `fxhash` with `rustc-hash` (maintained successor)
- CLI version now reads from Cargo.toml
- Improved error handling in FASTA readers

### Fixed

- Panic when sequence length is smaller than k-mer length
- Clippy warnings

## [0.1.0] - Initial Release

### Added

- Parallel k-mer counting using rayon and dashmap
- Canonical k-mer computation
- Support for rust-bio and needletail FASTA readers
- K-mer lengths from 1 to 32
- Jellyfish-compatible output format
