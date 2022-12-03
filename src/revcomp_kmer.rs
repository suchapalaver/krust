use crate::kmer::Kmer;

trait Complementary {
    fn parse_complement_byte(self) -> Self
    where
        Self: Sized;
}

impl Complementary for u8 {
    fn parse_complement_byte(self) -> Self {
        match self {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        }
    }
}

/// Converting a DNA string slice into its [reverse compliment](https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)#DNA_and_RNA_base_pair_complementarity).
pub struct RevCompKmer(pub Vec<u8>);

impl FromIterator<u8> for RevCompKmer {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        let mut c = RevCompKmer::new();

        for i in iter {
            c.add(i)
        }
        c
    }
}

impl RevCompKmer {
    fn new() -> Self {
        Self(vec![])
    }

    fn add(&mut self, elem: u8) {
        self.0.push(elem)
    }

    pub fn from_kmer(kmer: &Kmer) -> Self {
        kmer.0
            .iter()
            .rev()
            .map(|byte| byte.parse_complement_byte())
            .collect()
    }
}
