/// Unpacking compressed, bitpacked k-mer data.
#[derive(Hash, PartialEq, Eq)]
pub struct UnpackedKmer(pub Vec<u8>);

impl From<(u64, usize)> for UnpackedKmer {
    fn from(kmer_data: (u64, usize)) -> Self {
        let (kmer, k) = (kmer_data.0, kmer_data.1);
        let mut byte_string = Vec::with_capacity(k);
        for i in 0..k {
            let isolate = kmer << ((i * 2) + 64 - (k * 2));
            let base = isolate >> 62;
            let byte = UnpackedKmerByte::from(base);
            byte_string.push(byte.0);
        }
        UnpackedKmer(byte_string)
    }
}

/// Unpacking compressed, bitpacked k-mer data.
struct UnpackedKmerByte(u8);

impl From<u64> for UnpackedKmerByte {
    fn from(base: u64) -> Self {
        let unpacked_byte = match base {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => panic!("An invalid k-mer passed to here means we have a serious bug"),
        };
        UnpackedKmerByte(unpacked_byte)
    }
}
