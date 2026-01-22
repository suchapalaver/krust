//! Fuzz target for pack/unpack roundtrip.
//!
//! Tests that packing and unpacking is the identity operation
//! for valid DNA sequences.

#![no_main]

use bytes::Bytes;
use kmerust::kmer::{unpack_to_bytes, Kmer, KmerLength};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // Filter to valid DNA sequences only
    if data.is_empty() || data.len() > 32 {
        return;
    }

    // Only test with valid DNA bases (we're testing pack/unpack, not validation)
    for &byte in data {
        if !matches!(byte, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
            return;
        }
    }

    let bytes = Bytes::copy_from_slice(data);
    let k = match KmerLength::new(data.len()) {
        Ok(k) => k,
        Err(_) => return,
    };

    // Pack the k-mer
    let kmer = match Kmer::from_sub(bytes.clone()) {
        Ok(k) => k,
        Err(_) => return,
    };
    let packed = kmer.pack();

    // Unpack and verify roundtrip
    let unpacked = unpack_to_bytes(packed.packed_bits(), k);

    // Normalize original to uppercase for comparison
    let normalized: Vec<u8> = data.iter().map(|b| b.to_ascii_uppercase()).collect();

    assert_eq!(
        unpacked.as_ref(),
        normalized.as_slice(),
        "Pack/unpack roundtrip failed"
    );
});
