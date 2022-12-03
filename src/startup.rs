use super::{
    bitpacked_kmer::BitpackedKmer, dashmaps::DashFx, kmer::Kmer, revcomp_kmer::RevCompKmer,
    unpacked_kmer::UnpackedKmer,
};
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

pub fn run<P: AsRef<Path> + std::fmt::Debug>(path: P, k: usize) -> Result<(), Box<dyn Error>> {
    let mut buf = BufWriter::new(stdout());

    build_kmer_map(path, k)?
        .into_iter()
        .par_bridge()
        .map(|(kmer, freq)| (UnpackedKmer::from_kmer_data(kmer, k).0, freq))
        .map(|(kmer, freq)| {
            let kmer = String::from_utf8(kmer).unwrap();
            (kmer, freq)
        })
        .collect::<HashMap<String, i32>>()
        .into_iter()
        .for_each(|(kmer, count)| {
            print_kmer_map(&mut buf, kmer, count);
        });

    buf.flush()?;

    Ok(())
}

/// Reads sequences from fasta records in parallel using [`rayon`](https://docs.rs/rayon/1.5.1/rayon/),
/// using a customized [`dashmap`](https://docs.rs/dashmap/4.0.2/dashmap/struct.DashMap.html)
/// with [`FxHasher`](https://docs.rs/fxhash/0.2.1/fxhash/struct.FxHasher.html) to update in parallel a
/// hashmap of canonical k-mers (keys) and their frequency in the data (values)
fn build_kmer_map<P: AsRef<Path> + std::fmt::Debug>(
    path: P,
    k: usize,
) -> Result<DashFx, Box<dyn Error>> {
    let kmer_map: DashFx = DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default());

    fasta::Reader::from_file(path)?
        .records()
        .into_iter()
        .par_bridge()
        .for_each(|r| {
            let record = r.expect("Error reading fasta record.");

            let seq: &[u8] = record.seq();

            process_seq(seq, &k, &kmer_map);
        });

    Ok(kmer_map)
}

/// Ignore substrings containing `N`
/// 
/// # Notes
/// Canonicalizes by lexicographically smaller of k-mer/reverse-complement
fn process_seq(seq: &[u8], k: &usize, kmer_map: &DashFx) {
    let mut i = 0;

    while i <= seq.len() - k {
        let sub = &seq[i..i + k];

        if let Ok(kmer) = Kmer::from_substring(sub) {
            process_valid_bytes(kmer_map, kmer);

            i += 1;
        } else {
            let invalid_byte_index = Kmer::find_invalid_byte_index(sub);

            i += invalid_byte_index + 1;
        }
    }
}

/// Convert a valid sequence substring from a bytes string to a u64
fn process_valid_bytes(kmer_map: &DashFx, kmer: Kmer) {
    let BitpackedKmer(bitpacked_kmer) = kmer.0.iter().cloned().collect();

    // If the k-mer as found in the sequence is already a key in the `Dashmap`,
    // increment its value and move on
    if let Some(mut freq) = kmer_map.get_mut(&bitpacked_kmer) {
        *freq += 1;
    } else {
        // Initialize the reverse complement of this so-far unrecorded k-mer
        let RevCompKmer(revcompkmer) = RevCompKmer::from_kmer(&kmer);

        // Find the alphabetically less of the k-mer substring and its reverse complement
        let canonical_kmer = Kmer::get_canonical_kmer(revcompkmer, kmer.0);

        // Compress the canonical k-mer into a bitpacked 64-bit unsigned integer
        let kmer: BitpackedKmer = canonical_kmer.0.into_iter().collect();

        // Add k-mer key and initial value to results
        *kmer_map.entry(kmer.0).or_insert(0) += 1;
    }
}

fn print_kmer_map(buf: &mut BufWriter<Stdout>, kmer: String, count: i32) {
    writeln!(buf, ">{}\n{}", count, kmer).expect("Unable to write output.");
}
