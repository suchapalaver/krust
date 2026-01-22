//! Property-based tests using proptest.
//!
//! These tests verify invariants that should hold across all valid inputs,
//! catching edge cases that might be missed by example-based tests.

use bytes::Bytes;
use kmerust::kmer::{unpack_to_bytes, unpack_to_string, Kmer, KmerLength};
use proptest::prelude::*;

/// Strategy for generating valid DNA sequences of length 1-32.
fn dna_sequence(min_len: usize, max_len: usize) -> impl Strategy<Value = String> {
    proptest::collection::vec(
        prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
        min_len..=max_len,
    )
    .prop_map(|chars| chars.into_iter().collect())
}

/// Strategy for generating valid k-mer lengths (1-32).
fn kmer_length() -> impl Strategy<Value = usize> {
    1usize..=32
}

proptest! {
    /// Packing and unpacking a k-mer should be the identity operation.
    #[test]
    fn pack_unpack_roundtrip(seq in dna_sequence(1, 32)) {
        let bytes = Bytes::from(seq.clone());
        let k = KmerLength::new(seq.len()).unwrap();

        let kmer = Kmer::from_sub(bytes).unwrap();
        let packed = kmer.pack();
        let unpacked = unpack_to_bytes(packed.packed_bits(), k);

        prop_assert_eq!(unpacked.as_ref(), seq.as_bytes());
    }

    /// Unpacking to string should produce the same result as unpacking to bytes.
    #[test]
    fn unpack_bytes_equals_string(seq in dna_sequence(1, 32)) {
        let bytes = Bytes::from(seq.clone());
        let k = KmerLength::new(seq.len()).unwrap();

        let kmer = Kmer::from_sub(bytes).unwrap();
        let packed = kmer.pack();

        let unpacked_bytes = unpack_to_bytes(packed.packed_bits(), k);
        let unpacked_string = unpack_to_string(packed.packed_bits(), k);

        prop_assert_eq!(unpacked_bytes.as_ref(), unpacked_string.as_bytes());
    }

    /// The canonical form should be idempotent: canonical(canonical(x)) == canonical(x).
    #[test]
    fn canonical_is_idempotent(seq in dna_sequence(1, 32)) {
        let bytes = Bytes::from(seq.clone());

        let canonical1 = Kmer::from_sub(bytes.clone()).unwrap().pack().canonical();
        let canonical_bytes = canonical1.bytes().clone();

        // Taking canonical of the already-canonical form should give same result
        let canonical2 = Kmer::from_sub(canonical_bytes).unwrap().pack().canonical();

        prop_assert_eq!(canonical1.packed_bits(), canonical2.packed_bits());
    }

    /// A k-mer and its reverse complement should have the same canonical form.
    #[test]
    fn kmer_and_rc_have_same_canonical(seq in dna_sequence(1, 32)) {
        let bytes = Bytes::from(seq.clone());

        // Compute reverse complement manually
        let rc: String = seq
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _ => unreachable!(),
            })
            .collect();

        let canonical1 = Kmer::from_sub(bytes).unwrap().pack().canonical();
        let canonical2 = Kmer::from_sub(Bytes::from(rc)).unwrap().pack().canonical();

        prop_assert_eq!(canonical1.packed_bits(), canonical2.packed_bits());
    }

    /// The canonical form should be lexicographically <= both the original and its RC.
    #[test]
    fn canonical_is_lexicographically_smallest(seq in dna_sequence(1, 32)) {
        let bytes = Bytes::from(seq.clone());

        // Compute reverse complement
        let rc: String = seq
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _ => unreachable!(),
            })
            .collect();

        let canonical = Kmer::from_sub(bytes).unwrap().pack().canonical();
        let canonical_str = std::str::from_utf8(canonical.bytes()).unwrap();

        // Canonical should be <= both original and reverse complement
        prop_assert!(canonical_str <= seq.as_str());
        prop_assert!(canonical_str <= rc.as_str());
    }

    /// KmerLength should accept values 1-32 and reject others.
    #[test]
    fn kmer_length_accepts_valid_range(k in 1usize..=32) {
        let result = KmerLength::new(k);
        prop_assert!(result.is_ok());
        prop_assert_eq!(result.unwrap().get(), k);
    }

    #[test]
    fn kmer_length_rejects_zero(_dummy in Just(())) {
        let result = KmerLength::new(0);
        prop_assert!(result.is_err());
    }

    #[test]
    fn kmer_length_rejects_too_large(k in 33usize..1000) {
        let result = KmerLength::new(k);
        prop_assert!(result.is_err());
    }

    /// Soft-masked (lowercase) sequences should produce same result as uppercase.
    #[test]
    fn soft_masked_equals_uppercase(seq in dna_sequence(1, 32)) {
        let uppercase = Bytes::from(seq.clone());
        let lowercase = Bytes::from(seq.to_lowercase());

        let packed_upper = Kmer::from_sub(uppercase).unwrap().pack();
        let packed_lower = Kmer::from_sub(lowercase).unwrap().pack();

        prop_assert_eq!(packed_upper.packed_bits(), packed_lower.packed_bits());
    }

    /// Mixed case should produce same result as all uppercase.
    #[test]
    fn mixed_case_equals_uppercase(seq in dna_sequence(1, 32)) {
        let uppercase = Bytes::from(seq.clone());

        // Create mixed case version
        let mixed: String = seq
            .chars()
            .enumerate()
            .map(|(i, c)| if i % 2 == 0 { c } else { c.to_ascii_lowercase() })
            .collect();
        let mixed_bytes = Bytes::from(mixed);

        let packed_upper = Kmer::from_sub(uppercase).unwrap().pack();
        let packed_mixed = Kmer::from_sub(mixed_bytes).unwrap().pack();

        prop_assert_eq!(packed_upper.packed_bits(), packed_mixed.packed_bits());
    }

    /// Packing should be deterministic.
    #[test]
    fn packing_is_deterministic(seq in dna_sequence(1, 32)) {
        let bytes1 = Bytes::from(seq.clone());
        let bytes2 = Bytes::from(seq);

        let packed1 = Kmer::from_sub(bytes1).unwrap().pack();
        let packed2 = Kmer::from_sub(bytes2).unwrap().pack();

        prop_assert_eq!(packed1.packed_bits(), packed2.packed_bits());
    }

    /// Different sequences of the same length should have different packed representations.
    /// This tests that the packing is injective.
    #[test]
    fn different_sequences_different_packing(
        seq in dna_sequence(1, 16),
        mutation_pos in 0usize..16,
    ) {
        // Skip if mutation position is out of bounds
        prop_assume!(mutation_pos < seq.len());

        // Create a second sequence by changing one base
        let original_char = seq.chars().nth(mutation_pos).unwrap();
        let new_char = match original_char {
            'A' => 'C',
            'C' => 'G',
            'G' => 'T',
            'T' => 'A',
            _ => unreachable!(),
        };

        let mut chars: Vec<char> = seq.chars().collect();
        chars[mutation_pos] = new_char;
        let seq2: String = chars.into_iter().collect();

        let packed1 = Kmer::from_sub(Bytes::from(seq.clone())).unwrap().pack();
        let packed2 = Kmer::from_sub(Bytes::from(seq2)).unwrap().pack();

        prop_assert_ne!(packed1.packed_bits(), packed2.packed_bits());
    }

    /// The length of unpacked bytes should equal the original k-mer length.
    #[test]
    fn unpack_preserves_length(k in kmer_length(), bits in any::<u64>()) {
        let k_len = KmerLength::new(k).unwrap();
        let unpacked = unpack_to_bytes(bits, k_len);

        prop_assert_eq!(unpacked.len(), k);
    }

    /// Unpacked bytes should only contain valid bases (A, C, G, T).
    #[test]
    fn unpack_produces_valid_bases(k in kmer_length(), bits in any::<u64>()) {
        let k_len = KmerLength::new(k).unwrap();
        let unpacked = unpack_to_bytes(bits, k_len);

        for &byte in unpacked.iter() {
            prop_assert!(
                byte == b'A' || byte == b'C' || byte == b'G' || byte == b'T',
                "Invalid base: {}", byte as char
            );
        }
    }
}
