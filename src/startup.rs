use super::{
    kmer::Kmer,
    reader::{Needletail, RustBio, SequenceReader},
};
use bytes::Bytes;
use dashmap::DashMap;
use fxhash::FxHasher;
use rayon::prelude::*;
use std::{
    collections::HashMap,
    error::Error,
    fmt::Debug,
    hash::BuildHasherDefault,
    io::{stdout, BufWriter, Error as IoError, Write},
    path::Path,
};

custom_error::custom_error! { pub ProcessError
    ReadError{source: Box<dyn Error>} = "Unable to read input",
    WriteError{source: IoError} = "Unable to write output",
}

pub fn run<P>(path: P, k: usize, reader: bool) -> Result<(), ProcessError>
where
    P: AsRef<Path> + Debug,
{
    let (reader, name) = match reader {
        true => (Needletail::sequence_reader(path), "needletail"),
        false => (RustBio::sequence_reader(path), "rust-bio"),
    };

    println!("\nReading fasta with {} ...", name);

    DashFx::new().build(reader?, k)?.output(k)?;

    Ok(())
}

/// A custom `DashMap` w/ `FxHasher`.
///
/// # Notes
/// Useful: [Using a Custom Hash Function in Rust](https://docs.rs/hashers/1.0.1/hashers/#using-a-custom-hash-function-in-rust)
type DashFx = DashMap<u64, i32, BuildHasherDefault<FxHasher>>;

trait KmerMap {
    fn new() -> Self;
    fn build(
        self,
        sequences: impl Iterator<Item = Bytes>,
        k: usize,
    ) -> Result<Self, Box<dyn Error>>
    where
        Self: Sized;
    fn process_sequence(&self, seq: &Bytes, k: &usize);
    fn process_valid_bytes(&self, kmer: &mut Kmer);
    fn log(&self, kmer: &Kmer);
    fn output(self, k: usize) -> Result<(), ProcessError>;
}

impl KmerMap for DashFx {
    fn new() -> Self {
        DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default())
    }

    /// Reads sequences from fasta records in parallel using [`rayon`](https://docs.rs/rayon/1.5.1/rayon/),
    /// using a customized [`dashmap`](https://docs.rs/dashmap/4.0.2/dashmap/struct.DashMap.html)
    /// with [`FxHasher`](https://docs.rs/fxhash/0.2.1/fxhash/struct.FxHasher.html) to update in parallel a
    /// hashmap of canonical k-mers (keys) and their frequency in the data (values)
    fn build(
        self,
        sequences: impl Iterator<Item = Bytes>,
        k: usize,
    ) -> Result<Self, Box<dyn Error>> {
        for seq in sequences {
            self.process_sequence(&seq, &k)
        }

        Ok(self)
    }

    /// Ignore substrings containing `N`
    ///
    /// # Notes
    /// Canonicalizes by lexicographically smaller of k-mer/reverse-complement
    fn process_sequence(&self, seq: &Bytes, k: &usize) {
        let mut i = 0;

        while i <= seq.len() - k {
            let sub = seq.slice(i..i + k);

            if let Ok(mut kmer) = Kmer::from_sub(&sub) {
                self.process_valid_bytes(&mut kmer);

                i += 1;
            } else {
                let invalid_byte_index = Kmer::find_invalid(&sub);

                i += invalid_byte_index + 1;
            }
        }
    }

    /// Convert a valid sequence substring from a bytes string to a u64
    fn process_valid_bytes(&self, kmer: &mut Kmer) {
        kmer.pack();

        // If the k-mer as found in the sequence is already a key in the `Dashmap`,
        // increment its value and move on
        if let Some(mut count) = self.get_mut(&kmer.packed_bits) {
            *count += 1;
        } else {
            kmer.canonical();

            if kmer.reverse_complement {
                // Re-initialize packed bits
                kmer.packed_bits = Default::default();
                // Compress the canonical k-mer into a 64-bit unsigned integer
                kmer.pack();
            }

            self.log(kmer);
        }
    }

    fn log(&self, kmer: &Kmer) {
        *self.entry(kmer.packed_bits).or_insert(0) += 1
    }

    fn output(self, k: usize) -> Result<(), ProcessError> {
        let mut buf = BufWriter::new(stdout());

        for (kmer, count) in self
            .into_iter()
            .par_bridge()
            .map(|(packed_bits, count)| Kmer {
                packed_bits,
                count,
                ..Default::default()
            })
            .map(|mut kmer| {
                kmer.unpack(k);
                kmer
            })
            .map(|kmer| (kmer.bytes, kmer.count))
            .map(|(kmer, count)| {
                let kmer = String::from_utf8(kmer.to_vec()).unwrap();
                (kmer, count)
            })
            .collect::<HashMap<String, i32>>()
            .into_iter()
        {
            writeln!(buf, ">{}\n{}", count, kmer)?
        }

        buf.flush()?;

        Ok(())
    }
}
