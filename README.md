Run krust on the test data, searching for kmers of length 5 across all sequences, like this:

	$ cargo run --release 5 cerevisae.pan.fa > output.tsv

or, searching for kmers of length 21:

	$ cargo run --release 21 cerevisae.pan.fa > output.tsv

