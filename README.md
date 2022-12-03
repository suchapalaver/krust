# `krust`

`krust` is a [k-mer](https://en.wikipedia.org/wiki/K-mer) counter--a bioinformatics 101 tool for counting the frequency of substrings of length `k` within strings of DNA data. It's written in Rust and run from the command line. It takes a fasta file of DNA sequences and will output all canonical k-mers (the double helix means each k-mer has a [reverse complement](https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)#DNA_and_RNA_base_pair_complementarity)) and their frequency across all records in the given fasta file.

Run `krust` on the test data* in the [`krust` Github repo](https://github.com/suchapalaver/krust), searching for kmers of length 5, like this:

```bash
cargo run --release 5 your/local/path/to/cerevisae.pan.fa > output.tsv
```

or, searching for kmers of length 21:  

```bash
cargo run --release 21 your/local/path/to/cerevisae.pan.fa > output.tsv
```

`krust` prints to `stdout`, writing, on alternate lines:

```bash
>{frequency}  
{canonical k-mer}
>{frequency}  
{canonical k-mer}  
...
```  

`krust` uses the [`rust-bio`](https://docs.rs/bio/0.38.0/bio/), [`rayon`](https://docs.rs/rayon/1.5.1/rayon/), and [`dashmap`](https://docs.rs/crate/dashmap/4.0.2) Rust libraries.  
  
*Unusual, yes, to provide this data in the repo, but it's helped me spread word about what I'm doing.
