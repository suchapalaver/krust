//! Fuzz target for canonical k-mer computation.
//!
//! Tests that canonical form has the expected properties:
//! 1. Is idempotent
//! 2. k-mer and reverse complement have same canonical form
//! 3. Canonical form is lexicographically smallest

#![no_main]

use bytes::Bytes;
use kmerust::kmer::Kmer;
use libfuzzer_sys::fuzz_target;

/// Compute reverse complement of a DNA sequence.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => unreachable!(),
        })
        .collect()
}

fuzz_target!(|data: &[u8]| {
    // Filter to valid DNA sequences only
    if data.is_empty() || data.len() > 32 {
        return;
    }

    // Only test with valid uppercase DNA bases
    for &byte in data {
        if !matches!(byte, b'A' | b'C' | b'G' | b'T') {
            return;
        }
    }

    let bytes = Bytes::copy_from_slice(data);

    let kmer = match Kmer::from_sub(bytes.clone()) {
        Ok(k) => k,
        Err(_) => return,
    };
    let canonical = kmer.pack().canonical();

    // Property 1: Canonical is idempotent
    let canonical_bytes = canonical.bytes().clone();
    let canonical2 = Kmer::from_sub(canonical_bytes)
        .unwrap()
        .pack()
        .canonical();
    assert_eq!(
        canonical.packed_bits(),
        canonical2.packed_bits(),
        "Canonical is not idempotent"
    );

    // Property 2: k-mer and RC have same canonical
    let rc = reverse_complement(data);
    let rc_canonical = Kmer::from_sub(Bytes::from(rc.clone()))
        .unwrap()
        .pack()
        .canonical();
    assert_eq!(
        canonical.packed_bits(),
        rc_canonical.packed_bits(),
        "k-mer and RC have different canonical forms"
    );

    // Property 3: Canonical is lexicographically smallest
    let canonical_str = std::str::from_utf8(canonical.bytes()).unwrap();
    let original_str = std::str::from_utf8(data).unwrap();
    let rc_str = std::str::from_utf8(&rc).unwrap();

    assert!(
        canonical_str <= original_str,
        "Canonical {} > original {}",
        canonical_str,
        original_str
    );
    assert!(
        canonical_str <= rc_str,
        "Canonical {} > RC {}",
        canonical_str,
        rc_str
    );
});
