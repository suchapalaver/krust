use crate::bitpacked_kmer::BitpackedKmer;
use crate::canonical_kmer::CanonicalKmer;
use crate::dashmaps::DashFx;
use crate::kmer::Kmer;
use crate::revcomp_kmer::RevCompKmer;
use crate::unpacked_kmer::UnpackedKmer;
use bio::io::fasta;
use dashmap::DashMap;
use fxhash::FxHasher;
use rayon::prelude::*;
use std::{
    collections::HashMap,
    error::Error,
    hash::BuildHasherDefault,
    io::{BufWriter, Stdout, Write},
};

pub fn run(filepath: String, k: usize) -> Result<(), Box<dyn Error>> {
    let mut buf = BufWriter::new(std::io::stdout());

    let _print_results = build_kmer_map(filepath, k)?
        .into_iter()
        .par_bridge()
        .map(|(bitpacked_kmer, freq)| (UnpackedKmer::from((bitpacked_kmer, k)).0, freq))
        .map(|(unpacked_kmer, freq)| {
            let kmer_str = String::from_utf8(unpacked_kmer).unwrap();
            (kmer_str, freq)
        })
        .collect::<HashMap<String, i32>>()
        .into_iter()
        .for_each(|(kmer, count)| {
            print_kmer_map(&mut buf, kmer, count);
        });

    buf.flush()?;

    Ok(())
}

///  - Reads sequences from fasta records in parallel using [`rayon`](https://docs.rs/rayon/1.5.1/rayon/),
/// using a customized [`dashmap`](https://docs.rs/dashmap/4.0.2/dashmap/struct.DashMap.html)
/// with [`FxHasher`](https://docs.rs/fxhash/0.2.1/fxhash/struct.FxHasher.html) to update in parallel a
/// hashmap of canonical k-mers (keys) and their frequency in the data (values).
fn build_kmer_map(filepath: String, k: usize) -> Result<DashFx, Box<dyn Error>> {
    let kmer_map: DashFx = DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default());
    let _ = fasta::Reader::from_file(&filepath)?
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

///  - Ignores substrings containing `N`.
///  - Canonicalizes by lexicographically smaller of k-mer/reverse-complement.
fn process_seq(seq: &[u8], k: &usize, kmer_map: &DashFx) {
    let mut i = 0;
    while i <= seq.len() - k {
        let sub = &seq[i..i + k];
        let bytestring = Kmer::new(sub);
        if let Some(Kmer(valid_bytestring)) = bytestring {
            process_valid_bytes(kmer_map, valid_bytestring);
            i += 1;
        } else {
            let invalid_byte_index = Kmer::find_invalid(sub);
            i += invalid_byte_index + 1;
        }
    }
}

/// Converts a valid sequence substring from a bytes string to a u64.
fn process_valid_bytes(kmer_map: &DashFx, valid_bytestring: Vec<u8>) {
    let BitpackedKmer(bitpacked_kmer) = BitpackedKmer::from(&valid_bytestring);
    // If the k-mer as found in the sequence is already a key in the `Dashmap`,
    // increment its value and move on.
    if let Some(mut freq) = kmer_map.get_mut(&bitpacked_kmer) {
        *freq += 1;
    } else {
        // Initialize the reverse complement of this so-far unrecorded k-mer.
        let RevCompKmer(revcompkmer) = RevCompKmer::from(&valid_bytestring);
        // Find the alphabetically less of the k-mer substring and its reverse complement.
        let CanonicalKmer(canonical_kmer) = CanonicalKmer::from((revcompkmer, valid_bytestring));
        // Compress the canonical k-mer into a bitpacked 64-bit unsigned integer.
        let BitpackedKmer(kmer) = BitpackedKmer::from(&canonical_kmer);
        // Add k-mer key and initial value to results.
        *kmer_map.entry(kmer).or_insert(0) += 1;
    }
}

fn print_kmer_map(buf: &mut BufWriter<Stdout>, kmer: String, count: i32) {
    writeln!(buf, ">{}\n{}", count, kmer).expect("Unable to write output.");
}
