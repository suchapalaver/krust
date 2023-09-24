# `krust`: counts k-mers, written in rust

`krust` is a [k-mer](https://en.wikipedia.org/wiki/K-mer) counter - a bioinformatics 101 tool for counting the frequency of substrings of length `k` within strings of DNA data. `krust` is written in Rust and run from the command line. It takes a FASTA file of DNA sequences and will output all canonical k-mers (the double helix means each k-mer has a [reverse complement](https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)#DNA_and_RNA_base_pair_complementarity)) and their frequency across all records in the given data. `krust` is tested for accuracy against [jellyfish](https://github.com/gmarcais/Jellyfish).

```bash
krust: counts k-mers, written in rust

Usage: krust <k> <path>

Arguments:
  <k>     provides k length, e.g. 5
  <path>  path to a FASTA file, e.g. /home/lisa/bio/cerevisiae.pan.fa

Options:
  -h, --help     Print help information
  -V, --version  Print version information
```

`krust` supports either `rust-bio` or `needletail` to read FASTA record. Use the `--features` flag to select.  

Run `krust` with `rust-bio`'s fasta reader to count *5*-mers like this:

```bash
cargo run --release --features rust-bio -- 5 your/local/path/to/fasta_data.fa
```

or, searching for *21*-mers with `needletail` as the fasta reader, like this:  

```bash
cargo run --release --features needletail -- 21 your/local/path/to/fasta_data.fa
```

`krust` prints to `stdout`, writing, on alternate lines:

```bash
>114928
ATGCC
>289495
AATCA
...
```  
