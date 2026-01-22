//! Fuzz target for `Kmer::from_sub`.
//!
//! Tests that `from_sub` handles arbitrary byte input gracefully,
//! either accepting valid DNA sequences or rejecting invalid ones.

#![no_main]

use bytes::Bytes;
use kmerust::kmer::Kmer;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // Limit input size to reasonable k-mer lengths
    if data.len() > 32 {
        return;
    }

    let bytes = Bytes::copy_from_slice(data);

    // from_sub should either succeed or fail gracefully - never panic
    match Kmer::from_sub(bytes) {
        Ok(kmer) => {
            // If successful, the bytes should only contain valid bases
            for &byte in kmer.bytes().iter() {
                assert!(
                    byte == b'A' || byte == b'C' || byte == b'G' || byte == b'T',
                    "Invalid base in accepted kmer: {}",
                    byte as char
                );
            }

            // Length should be preserved
            assert_eq!(kmer.bytes().len(), data.len());
        }
        Err(err) => {
            // Error should report a valid position
            assert!(
                err.position < data.len(),
                "Error position {} out of bounds for data len {}",
                err.position,
                data.len()
            );

            // Error should report the actual invalid byte
            assert_eq!(
                err.base, data[err.position],
                "Error byte mismatch at position {}",
                err.position
            );
        }
    }
});
