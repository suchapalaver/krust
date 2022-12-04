use super::kmer::{Bitpack, Kmer, RevComp, Unpack};
use bio::io::fasta;
use dashmap::DashMap;
use fxhash::FxHasher;
use rayon::prelude::*;
use std::{
    collections::HashMap,
    error::Error,
    hash::BuildHasherDefault,
    io::{stdout, BufWriter, Stdout, Write},
    path::Path,
};

/// A custom `DashMap` w/ `FxHasher`.
///
/// ```use dashmap::DashMap;```
/// ```use fxhash::FxHasher;```
/// ```// skip```
/// ```let dashfx_hash: DashFx = DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default());```
///
/// # Notes
/// Useful: [Using a Custom Hash Function in Rust](https://docs.rs/hashers/1.0.1/hashers/#using-a-custom-hash-function-in-rust)
type DashFx = DashMap<u64, i32, BuildHasherDefault<FxHasher>>;

pub fn run<P: AsRef<Path> + std::fmt::Debug>(path: P, k: usize) -> Result<(), Box<dyn Error>> {
    let mut buf = BufWriter::new(stdout());

    build_map(path, k)?
        .into_iter()
        .par_bridge()
        .map(|(kmer, freq)| (Unpack::bit(kmer, k).0, freq))
        .map(|(kmer, freq)| {
            let kmer = String::from_utf8(kmer).unwrap();
            (kmer, freq)
        })
        .collect::<HashMap<String, i32>>()
        .into_iter()
        .for_each(|(kmer, count)| {
            output(&mut buf, kmer, count);
        });

    buf.flush()?;

    Ok(())
}

/// Reads sequences from fasta records in parallel using [`rayon`](https://docs.rs/rayon/1.5.1/rayon/),
/// using a customized [`dashmap`](https://docs.rs/dashmap/4.0.2/dashmap/struct.DashMap.html)
/// with [`FxHasher`](https://docs.rs/fxhash/0.2.1/fxhash/struct.FxHasher.html) to update in parallel a
/// hashmap of canonical k-mers (keys) and their frequency in the data (values)
fn build_map<P: AsRef<Path> + std::fmt::Debug>(
    path: P,
    k: usize,
) -> Result<DashFx, Box<dyn Error>> {
    let map: DashFx = DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default());

    fasta::Reader::from_file(path)?
        .records()
        .into_iter()
        .par_bridge()
        .for_each(|r| {
            let record = r.expect("Error reading fasta record.");

            let seq: &[u8] = record.seq();

            process_seq(seq, &k, &map);
        });

    Ok(map)
}

/// Ignore substrings containing `N`
///
/// # Notes
/// Canonicalizes by lexicographically smaller of k-mer/reverse-complement
fn process_seq(seq: &[u8], k: &usize, kmer_map: &DashFx) {
    let mut i = 0;

    while i <= seq.len() - k {
        let sub = &seq[i..i + k];

        if let Ok(kmer) = Kmer::from_sub(sub) {
            process_valid_bytes(kmer_map, kmer);

            i += 1;
        } else {
            let invalid_byte_index = Kmer::find_invalid(sub);

            i += invalid_byte_index + 1;
        }
    }
}

/// Convert a valid sequence substring from a bytes string to a u64
fn process_valid_bytes(kmer_map: &DashFx, kmer: Kmer) {
    let Bitpack(bitpacked_kmer) = kmer.0.iter().cloned().collect();

    // If the k-mer as found in the sequence is already a key in the `Dashmap`,
    // increment its value and move on
    if let Some(mut freq) = kmer_map.get_mut(&bitpacked_kmer) {
        *freq += 1;
    } else {
        // Initialize the reverse complement of this so-far unrecorded k-mer
        let RevComp(revcompkmer) = RevComp::from_kmer(&kmer);

        // Find the alphabetically less of the k-mer substring and its reverse complement
        let canonical_kmer = Kmer::canonical(revcompkmer, kmer.0);

        // Compress the canonical k-mer into a bitpacked 64-bit unsigned integer
        let kmer: Bitpack = canonical_kmer.collect();

        // Add k-mer key and initial value to results
        *kmer_map.entry(kmer.0).or_insert(0) += 1;
    }
}

fn output(buf: &mut BufWriter<Stdout>, kmer: String, count: i32) {
    writeln!(buf, ">{}\n{}", count, kmer).expect("Unable to write output.");
}
