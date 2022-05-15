/// Converting a DNA string slice into its [reverse compliment](https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)#DNA_and_RNA_base_pair_complementarity).
pub struct RevCompKmer(pub Vec<u8>);

impl From<&Vec<u8>> for RevCompKmer {
    fn from(sub: &Vec<u8>) -> Self {
        let mut revcomp = Vec::with_capacity(sub.len());

        for byte in sub.iter().rev() {
            let comp = RevCompKmer::complement(*byte);
            revcomp.push(comp);
        }
        RevCompKmer(revcomp)
    }
}

impl RevCompKmer {
    fn complement(byte: u8) -> u8 {
        match byte {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => panic!("`RevCompKmer::from` should only be passed valid k-mers"),
        }
    }
}
