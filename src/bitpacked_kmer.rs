/// Compressing k-mers of length `0 < k < 33`, bitpacking them into unsigned integers
pub struct BitpackedKmer(pub u64);

impl BitpackedKmer {
    fn new() -> BitpackedKmer {
        BitpackedKmer(0)
    }

    fn pack(&mut self, elem: u8) {
        self.shift();
        let mask = elem.pack_convert();
        self.0 |= mask
    }

    fn shift(&mut self) {
        self.0 <<= 2
    }
}

impl FromIterator<u8> for BitpackedKmer {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        let mut c = BitpackedKmer::new();

        for i in iter {
            c.pack(i)
        }
        c
    }
}

trait Pack {
    fn pack_convert(self) -> u64
    where
        Self: Sized;
}

impl Pack for u8 {
    fn pack_convert(self) -> u64 {
        match self {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            // can only be b'T'
            _ => 3,
        }
    }
}
