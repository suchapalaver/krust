/// Compressing k-mers of length `0 < k < 33`, bitpacking them into unsigned integers.
pub struct BitpackedKmer(pub u64);

impl From<&Vec<u8>> for BitpackedKmer {
    fn from(sub: &Vec<u8>) -> Self {
        let bitpacked_kmer: u64 = {
            let mut k: u64 = 0;
            for byte in sub.iter() {
                k <<= 2;
                let mask = match *byte {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => panic!("`BitpackerKmer` handling an invalid k-mer bytestring is unexpected behavior"),
                };
                k |= mask;
            }
            k
        };
        BitpackedKmer(bitpacked_kmer)
    }
}
