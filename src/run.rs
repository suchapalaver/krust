use super::{kmer::Kmer, reader::read};
use bytes::Bytes;
use dashmap::DashMap;
use fxhash::FxHasher;
use rayon::prelude::{ParallelBridge, ParallelIterator};
use std::{
    collections::{hash_map::IntoIter, HashMap},
    error::Error,
    fmt::Debug,
    hash::BuildHasherDefault,
    io::{stdout, BufWriter, Error as IoError, Write},
    path::Path,
};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ProcessError {
    #[error("Unable to read input: {0}")]
    ReadError(#[from] Box<dyn Error>),

    #[error("Unable to write output: {0}")]
    WriteError(#[from] IoError),
}

pub fn run<P>(path: P, k: usize) -> Result<(), ProcessError>
where
    P: AsRef<Path> + Debug,
{
    KmerMap::new().build(read(path)?, k)?.output(k)?;

    Ok(())
}

/// A custom `DashMap` w/ `FxHasher`.
///
/// # Notes
/// Useful: [Using a Custom Hash Function in Rust](https://docs.rs/hashers/1.0.1/hashers/#using-a-custom-hash-function-in-rust)
type DashFx = DashMap<u64, i32, BuildHasherDefault<FxHasher>>;

struct KmerMap(DashFx);

impl KmerMap {
    fn new() -> Self {
        Self(DashMap::with_hasher(
            BuildHasherDefault::<FxHasher>::default(),
        ))
    }

    /// Reads sequences from fasta records in parallel using [`rayon`](https://docs.rs/rayon/1.5.1/rayon/),
    /// using a customized [`dashmap`](https://docs.rs/dashmap/4.0.2/dashmap/struct.DashMap.html)
    /// with [`FxHasher`](https://docs.rs/fxhash/0.2.1/fxhash/struct.FxHasher.html) to update in parallel a
    /// hashmap of canonical k-mers (keys) and their frequency in the data (values)
    fn build(
        self,
        sequences: rayon::vec::IntoIter<Bytes>,
        k: usize,
    ) -> Result<Self, Box<dyn Error>> {
        sequences.for_each(|seq| self.process_sequence(&seq, &k));

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

            match Kmer::from_sub(sub) {
                Ok(mut kmer) => self.process_valid_bytes(&mut kmer),
                Err(invalid_byte_index) => i += invalid_byte_index,
            }

            i += 1
        }
    }

    /// Convert a valid sequence substring from a bytes string to a u64
    fn process_valid_bytes(&self, kmer: &mut Kmer) {
        kmer.pack();

        // If the k-mer as found in the sequence is already a key in the `Dashmap`,
        // increment its value and move on
        if let Some(mut count) = self.0.get_mut(&kmer.packed_bits) {
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
        *self.0.entry(kmer.packed_bits).or_insert(0) += 1
    }

    fn output(self, k: usize) -> Result<(), ProcessError> {
        let mut buf = BufWriter::new(stdout());

        for (kmer, count) in self.stream(k) {
            writeln!(buf, ">{count}\n{kmer}")?
        }

        buf.flush()?;

        Ok(())
    }

    fn stream(self, k: usize) -> IntoIter<String, i32> {
        self.0
            .into_iter()
            .par_bridge()
            .map(|(packed_bits, count)| Kmer {
                packed_bits,
                count,
                ..Default::default()
            })
            .map(|mut kmer| {
                kmer.unpack(k);
                (String::from_utf8(kmer.bytes.to_vec()).unwrap(), kmer.count)
            })
            .collect::<HashMap<String, i32>>()
            .into_iter()
    }
}
