# krust

[![Crates.io](https://img.shields.io/crates/v/krust.svg)](https://crates.io/crates/krust)
[![Documentation](https://docs.rs/krust/badge.svg)](https://docs.rs/krust)
[![CI](https://github.com/suchapalaver/krust/workflows/CI/badge.svg)](https://github.com/suchapalaver/krust/actions)
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
cargo install krust
```

### From source

```bash
git clone https://github.com/suchapalaver/krust.git
cd krust
cargo install --path .
```

## Usage

```bash
krust <k> <path>
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
krust 21 sequences.fa > kmers.txt
```

Count 5-mers:

```bash
krust 5 sequences.fa > kmers.txt
```

### FASTA Readers

krust supports two FASTA readers via feature flags:

- `rust-bio` (default) - Uses the [rust-bio](https://docs.rs/bio) library
- `needletail` - Uses the [needletail](https://docs.rs/needletail) library

To use needletail instead:

```bash
cargo run --release --no-default-features --features needletail -- 21 sequences.fa
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

krust can also be used as a library:

```rust
use krust::run::run;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = PathBuf::from("sequences.fa");
    run(path, 21)?;
    Ok(())
}
```

## Performance

krust uses parallel processing to efficiently count k-mers:

- Sequences are processed in parallel using rayon
- A concurrent hash map (dashmap) allows lock-free updates
- FxHash provides fast hashing for 64-bit packed k-mers

## License

MIT License - see [LICENSE](LICENSE) for details.
