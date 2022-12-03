/// Unpacking compressed, bitpacked k-mer data.
#[derive(Hash, PartialEq, Eq)]
pub struct UnpackedKmer(pub Vec<u8>);

impl UnpackedKmer {
    fn new() -> UnpackedKmer {
        UnpackedKmer(Vec::new())
    }

    fn add(&mut self, elem: u8) {
        self.0.push(elem);
    }

    pub fn from_kmer_data(kmer: u64, k: usize) -> Self {
        (0..k)
            .into_iter()
            .map(|i| kmer.isolate_bits(i, k).replace_bits().unpack_bits())
            .collect()
    }
}

trait Unpack {
    fn unpack_bits(self) -> u8
    where
        Self: Sized;

    fn isolate_bits(self, i: usize, k: usize) -> Self
    where
        Self: Sized;

    fn replace_bits(self) -> Self
    where
        Self: Sized;
}

impl Unpack for u64 {
    fn unpack_bits(self: u64) -> u8 {
        match self {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        }
    }
    fn isolate_bits(self: u64, i: usize, k: usize) -> Self {
        self << ((i * 2) + 64 - (k * 2))
    }

    fn replace_bits(self) -> Self {
        self >> 62
    }
}

impl FromIterator<u8> for UnpackedKmer {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        let mut c = UnpackedKmer::new();

        for i in iter {
            c.add(i)
        }
        c
    }
}
