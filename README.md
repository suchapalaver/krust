To make use of krust's multithreaded cocurrency test fasta data with multiple records.

Run krust on the test data, searching for kmers of length 5, like this:

	$ cargo run 5 cerevisae.pan_S288C_chrI.fa

or, searching for kmers of length 21 across multiple records:

	$ cargo run 21 cerevisae.pan.fa
	


On branches:

-main is the same as ParallelBridge

-ParallelBridge uses rayon's ParallelBridge (https://docs.rs/rayon/1.5.1/rayon/iter/trait.ParallelBridge.html)

-collect_the_iterator is able to use rayon by 'collecting' the rust-bio fasta Reader into a vector of Records

-std_threads uses Rust's standard library thread to achieve multithreadedness

-single-thread has no multithreaded processing
