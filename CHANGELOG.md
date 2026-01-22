# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

#### Type Safety & API Correctness

- `KmerLength` newtype that enforces valid k-mer lengths (1-32) at construction time
- Type-state pattern for `Kmer<S>` with `Unpacked`, `Packed`, and `Canonical` states
- Exhaustive error types (`KmeRustError`, `KmerLengthError`, `InvalidBaseError`) replacing `Box<dyn Error>`
- Compile-time enforcement of k-mer operation ordering

#### Performance Hardening

- Compile-time lookup tables (`PACK_TABLE`, `COMPLEMENT_TABLE`, `UNPACK_TABLE`) for bit packing
- Zero-allocation canonical comparison using lazy iterator comparison
- Streaming API (`count_kmers_streaming`, `count_kmers_streaming_packed`) for memory-efficient processing
- `count_kmers_from_sequences` for counting from in-memory byte slices

#### Testing Rigor

- Direct library API test suite (`tests/library_tests.rs`)
- Property-based tests using proptest (`tests/property_tests.rs`)
- Jellyfish compatibility tests (`tests/jellyfish_compat.rs`)
- Fuzz targets for input handling (`fuzz/fuzz_targets/`)

#### Production Readiness

- Gzip compressed input support (`gzip` feature)
- Memory-mapped I/O for large files (`mmap` feature)
- Progress callback API (`count_kmers_with_progress`)
- Tracing instrumentation (`tracing` feature)

#### API Polish

- Builder pattern API (`KmerCounter`) for ergonomic k-mer counting
- Async API with Tokio support (`async` feature)
  - `count_kmers_async` for non-blocking k-mer counting
  - `AsyncKmerCounter` builder for async operations
- Example programs demonstrating common use cases:
  - `basic_count.rs` - Simple k-mer counting
  - `streaming_large_file.rs` - Memory-efficient processing
  - `compare_with_jellyfish.rs` - Validation against Jellyfish
  - `progress_bar.rs` - Progress reporting during counting
  - `async_count.rs` - Async k-mer counting
- `full` feature flag enabling all optional features

#### Infrastructure

- Comprehensive test suite with unit and integration tests
- CI/CD pipeline with GitHub Actions
- Rustdoc documentation for public APIs
- MSRV (Minimum Supported Rust Version) set to 1.75
- Semver checking with cargo-semver-checks in CI

### Changed

- `count_kmers` now validates k-mer length upfront and returns a clear error
- K-mer packing uses direct lookup tables instead of match statements
- Canonical computation avoids allocation when the original form is already canonical
- Updated all dependencies to latest stable versions
- Replaced `fxhash` with `rustc-hash` (maintained successor)
- CLI version now reads from Cargo.toml
- Improved error handling in FASTA readers

### Removed

- Deprecated `config.rs` module

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

[Unreleased]: https://github.com/suchapalaver/kmerust/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/suchapalaver/kmerust/releases/tag/v0.1.0
