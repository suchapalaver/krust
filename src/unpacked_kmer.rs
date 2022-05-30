/// Unpacking compressed, bitpacked k-mer data.
#[derive(Hash, PartialEq, Eq)]
pub struct UnpackedKmer(pub Vec<u8>);

impl UnpackedKmer {
    fn new(k: usize) -> UnpackedKmer {
        UnpackedKmer(Vec::with_capacity(k))
    }

    fn add(&mut self, elem: u8) {
        self.0.push(elem);
    }
}

impl From<(u64, usize)> for UnpackedKmer {
    fn from(kmer_data: (u64, usize)) -> Self {
        let (kmer, k) = (kmer_data.0, kmer_data.1);
        let mut unpacked_kmer = UnpackedKmer::new(k);
        for i in 0..k {
            let isolate = kmer << ((i * 2) + 64 - (k * 2));
            let base = isolate >> 62;
            let byte = UnpackedKmerByte::from(base);
            unpacked_kmer.add(byte.0);
        }
        unpacked_kmer
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
