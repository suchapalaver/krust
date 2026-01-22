# kmerust

[![Crates.io](https://img.shields.io/crates/v/kmerust.svg)](https://crates.io/crates/kmerust)
[![Documentation](https://docs.rs/kmerust/badge.svg)](https://docs.rs/kmerust)
[![CI](https://github.com/suchapalaver/kmerust/workflows/CI/badge.svg)](https://github.com/suchapalaver/kmerust/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A fast, parallel [k-mer](https://en.wikipedia.org/wiki/K-mer) counter for DNA sequences in FASTA files.

## Features

- **Fast parallel processing** using [rayon](https://docs.rs/rayon) and [dashmap](https://docs.rs/dashmap)
- **Canonical k-mers** - outputs the lexicographically smaller of each k-mer and its reverse complement
- **Flexible k-mer lengths** from 1 to 32
- **Handles N bases** by skipping invalid k-mers
- **Jellyfish-compatible output** format for easy integration with existing pipelines
- **Tested for accuracy** against [Jellyfish](https://github.com/gmarcais/Jellyfish)

## Installation

### From crates.io

```bash
cargo install kmerust
```

### From source

```bash
git clone https://github.com/suchapalaver/kmerust.git
cd kmerust
cargo install --path .
```

## Usage

```bash
kmerust <k> <path>
```

### Arguments

- `<k>` - K-mer length (1-32)
- `<path>` - Path to a FASTA file

### Options

- `-h, --help` - Print help information
- `-V, --version` - Print version information

### Examples

Count 21-mers in a FASTA file:

```bash
kmerust 21 sequences.fa > kmers.txt
```

Count 5-mers:

```bash
kmerust 5 sequences.fa > kmers.txt
```

### FASTA Readers

kmerust supports two FASTA readers via feature flags:

- `rust-bio` (default) - Uses the [rust-bio](https://docs.rs/bio) library
- `needletail` - Uses the [needletail](https://docs.rs/needletail) library

To use needletail instead:

```bash
cargo run --release --no-default-features --features needletail -- 21 sequences.fa
```

### Production Features

Enable production features for additional capabilities:

```bash
cargo build --release --features production
```

Or enable individual features:

- `gzip` - Read gzip-compressed FASTA files (`.fa.gz`)
- `mmap` - Memory-mapped I/O for large files
- `tracing` - Structured logging and diagnostics

#### Gzip Compressed Input

With the `gzip` feature, kmerust can directly read gzip-compressed files:

```bash
cargo run --release --features gzip -- 21 sequences.fa.gz
```

#### Tracing/Logging

With the `tracing` feature, use the `RUST_LOG` environment variable for diagnostic output:

```bash
RUST_LOG=kmerust=debug cargo run --features tracing -- 21 sequences.fa
```

## Output Format

Output is written to stdout in FASTA-like format:

```
>{count}
{canonical_kmer}
```

Example output:

```
>114928
ATGCC
>289495
AATCA
```

## Library Usage

kmerust can also be used as a library:

```rust
use kmerust::run::count_kmers;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = PathBuf::from("sequences.fa");
    let counts = count_kmers(&path, 21)?;
    for (kmer, count) in counts {
        println!("{kmer}: {count}");
    }
    Ok(())
}
```

### Progress Reporting

Monitor progress during long-running operations:

```rust
use kmerust::run::count_kmers_with_progress;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let counts = count_kmers_with_progress("genome.fa", 21, |progress| {
        eprintln!(
            "Processed {} sequences ({} bases)",
            progress.sequences_processed,
            progress.bases_processed
        );
    })?;
    Ok(())
}
```

### Memory-Mapped I/O

For large files, use memory-mapped I/O (requires `mmap` feature):

```rust
use kmerust::run::count_kmers_mmap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let counts = count_kmers_mmap("large_genome.fa", 21)?;
    println!("Found {} unique k-mers", counts.len());
    Ok(())
}
```

### Streaming API

For memory-efficient processing:

```rust
use kmerust::streaming::count_kmers_streaming;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let counts = count_kmers_streaming("genome.fa", 21)?;
    println!("Found {} unique k-mers", counts.len());
    Ok(())
}
```

## Performance

kmerust uses parallel processing to efficiently count k-mers:

- Sequences are processed in parallel using rayon
- A concurrent hash map (dashmap) allows lock-free updates
- FxHash provides fast hashing for 64-bit packed k-mers

## License

MIT License - see [LICENSE](LICENSE) for details.
