To make use of k-rust's multithreaded cocurrency test with the larger fasta data with multiple records.

Run k-rust on the test data, searching for kmers of length 5, like this:

$ cargo run 5 cerevisae.pan_S288C_chrI.fa

or, searching for kmers of length 21 across multiple records:

$ cargo run 21 cerevisae.pan.fa
