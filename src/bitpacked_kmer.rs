/// Compressing k-mers of length `0 < k < 33`, bitpacking them into unsigned integers.
pub struct BitpackedKmer(pub u64);

impl BitpackedKmer {
    fn new() -> BitpackedKmer {
        BitpackedKmer(0)
    }

    fn pack(&mut self, elem: u8) {
        self.0 <<= 2;
        let mask = match elem {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => panic!("`BitpackerKmer` handling an invalid k-mer bytestring is unexpected behavior"),
        };
        self.0 |= mask;
    }
}

impl FromIterator<u8> for BitpackedKmer {
    fn from_iter<I: IntoIterator<Item=u8>>(iter: I) -> Self {
        let mut c = BitpackedKmer::new();

        for i in iter {
            c.pack(i);
        }
        c
    }
}

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
