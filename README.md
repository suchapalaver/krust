`krust` is a [k-mer](https://en.wikipedia.org/wiki/K-mer) counter written in Rust and run from the command line that will output canonical k-mers and their frequency across the records in a fasta file. 

Run `krust` on the test data in the [`krust` Github repo](https://github.com/suchapalaver/krust), searching for kmers of length 5, like this:  
```$ cargo run --release 5 cerevisae.pan.fa > output.tsv```  
or, searching for kmers of length 21:  
```$ cargo run --release 21 cerevisae.pan.fa > output.tsv``` 

`krust` prints to `stdout`, writing, on alternate lines:  
```>{frequency}```  
```{canonical k-mer}```  
```>{frequency}
```>(canonical k-mer}```
...

`krust` uses [`rust-bio`](https://docs.rs/bio/0.38.0/bio/), [`rayon`](https://docs.rs/rayon/1.5.1/rayon/), and [`dashmap`](https://docs.rs/crate/dashmap/4.0.2).  

Future:  
A function like fn single_sequence_canonical_kmers(filepath: String, k: usize) {}    
Would returns k-mer counts for individual sequences in a fasta file.     
