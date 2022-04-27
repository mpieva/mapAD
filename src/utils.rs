use bio::{
    alphabets::{self, RankTransform},
    data_structures::{
        bwt::{bwt, less, Occ},
        suffix_array::{suffix_array, RawSuffixArray},
    },
};

use crate::map::fmd_index::RtFmdIndex;

/// This is only used in tests and benchmarks
pub fn build_auxiliary_structures(
    mut reference: Vec<u8>,
    mut src_alphabet: alphabets::Alphabet,
) -> (RtFmdIndex, RawSuffixArray) {
    let ref_seq_revcomp = alphabets::dna::revcomp(reference.iter());
    reference.extend_from_slice(b"$");
    reference.extend_from_slice(&ref_seq_revcomp);
    drop(ref_seq_revcomp);
    reference.extend_from_slice(b"$");

    src_alphabet.insert(b'$');
    let rank_transform = RankTransform::new(&src_alphabet);
    reference = rank_transform.transform(reference);
    let rank_alphabet = alphabets::Alphabet::new(rank_transform.ranks.values());

    let sar = suffix_array(&reference);
    let bwt = bwt(&reference, &sar);
    let less = less(&bwt, &rank_alphabet);
    let occ = Occ::new(&bwt, 3, &rank_alphabet);

    (RtFmdIndex::new(bwt, less, occ, rank_transform), sar)
}
