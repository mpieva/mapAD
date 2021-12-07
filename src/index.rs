use std::fs::File;

use either::Either;
use log::info;
use rand::{
    prelude::{Rng, SeedableRng, StdRng},
    seq::SliceRandom,
};
use serde::{Deserialize, Serialize};

use bio::{
    alphabets::{dna, Alphabet, RankTransform},
    data_structures::{
        bwt::{bwt, less, Less, Occ, BWT},
        suffix_array::{suffix_array, RawSuffixArray},
    },
    io::fasta,
};

use crate::{
    errors::Result,
    map::{FastaIdPosition, FastaIdPositions},
};

// Increase this number once the on-disk index changes
pub const INDEX_VERSION: u8 = 0;

pub const DNA_UPPERCASE_ALPHABET: &[u8; 4] = b"ACGT";
// Ambiguous base symbols (which appear in stretches) can be replaced with 'X' in the index
pub const DNA_UPPERCASE_X_ALPHABET: &[u8; 5] = b"ACGTX";
const DNA_PURINE: &[u8; 2] = b"AG";
const DNA_PYRIMIDINE: &[u8; 2] = b"CT";
const DNA_KETONE: &[u8; 2] = b"GT";
const DNA_AMINO: &[u8; 2] = b"AC";
const DNA_STRONG: &[u8; 2] = b"CG";
const DNA_WEAK: &[u8; 2] = b"AT";
const DNA_NOT_A: &[u8; 3] = b"CGT";
const DNA_NOT_C: &[u8; 3] = b"AGT";
const DNA_NOT_G: &[u8; 3] = b"ACT";
const DNA_NOT_T: &[u8; 3] = b"ACG";

// Versioned index data structures
#[derive(Serialize, Deserialize)]
pub struct VersionedBwt {
    pub version: u8,
    pub data: BWT,
}

#[derive(Serialize, Deserialize)]
pub struct VersionedLess {
    pub version: u8,
    pub data: Less,
}

#[derive(Serialize, Deserialize)]
pub struct VersionedOcc {
    pub version: u8,
    pub data: Occ,
}

#[derive(Serialize, Deserialize)]
pub struct VersionedRt {
    pub version: u8,
    pub data: RankTransform,
}

#[derive(Serialize, Deserialize)]
pub struct VersionedIdPosMap {
    pub version: u8,
    pub data: FastaIdPositions,
}

#[derive(Serialize, Deserialize)]
pub struct VersionedSuffixArray {
    pub version: u8,
    pub data: RawSuffixArray,
}

/// Entry point function to launch the indexing process
pub fn run(reference_path: &str, seed: u64) -> Result<()> {
    let mut rng: StdRng = SeedableRng::seed_from_u64(seed);

    let alphabet = Alphabet::new(crate::index::DNA_UPPERCASE_X_ALPHABET);

    // Index the genome
    index(reference_path, alphabet, reference_path, &mut rng)?;

    Ok(())
}

/// Index a given reference and write the index to disk.
/// Ambiguous bases ('N') are converted to random bases unless they occur in stretches.
/// If they do, they will be converted to 'X'.
fn index<T: Rng>(
    reference_path: &str,
    mut alphabet: Alphabet,
    name: &str,
    rng: &mut T,
) -> Result<()> {
    info!("Read input reference sequence");
    let mut ref_seq = fasta::Reader::from_file(reference_path)?
        .records()
        // Convert all bases to uppercase
        .flat_map(|record| match record {
            Ok(record) => Either::Left(record.seq().to_ascii_uppercase().into_iter().map(Ok)),
            Err(e) => Either::Right(std::iter::once(e).map(|e| Err(e.into()))),
        })
        .collect::<Result<Vec<_>>>()?;

    // Replace single occurrences of ambiguous base symbols with random ones, leave runs alone
    let randomly_replace_ambiguous = |base| match base {
        b'U' => b'T',
        b'R' => *DNA_PURINE.choose(rng).unwrap(),
        b'Y' => *DNA_PYRIMIDINE.choose(rng).unwrap(),
        b'K' => *DNA_KETONE.choose(rng).unwrap(),
        b'M' => *DNA_AMINO.choose(rng).unwrap(),
        b'S' => *DNA_STRONG.choose(rng).unwrap(),
        b'W' => *DNA_WEAK.choose(rng).unwrap(),
        b'B' => *DNA_NOT_A.choose(rng).unwrap(),
        b'D' => *DNA_NOT_C.choose(rng).unwrap(),
        b'H' => *DNA_NOT_G.choose(rng).unwrap(),
        b'V' => *DNA_NOT_T.choose(rng).unwrap(),
        b'N' => *DNA_UPPERCASE_ALPHABET.choose(rng).unwrap(),
        _ => base,
    };

    let summarize_ambiguous = |base| {
        if !DNA_UPPERCASE_ALPHABET.contains(&base) {
            b'X'
        } else {
            base
        }
    };

    // Start by replacing ambiguous bases that do not occur in stretches exceeding a certain length,
    // using the closures defined above
    run_apply(
        &mut ref_seq,
        10,
        randomly_replace_ambiguous,
        summarize_ambiguous,
    );

    {
        info!("Map identifiers to positions");
        let mut end = 0;
        let identifier_position_map = FastaIdPositions::new(
            fasta::Reader::from_file(reference_path)?
                .records()
                .map(|record| {
                    let record = record.expect("Failed reading input file");
                    end += record.seq().len();
                    FastaIdPosition {
                        start: end - record.seq().len(),
                        end: end - 1,
                        identifier: record.id().to_string(),
                    }
                })
                .collect::<Vec<_>>(),
        );

        info!("Save position map to disk");
        let versioned_id_pos_map = VersionedIdPosMap {
            version: INDEX_VERSION,
            data: identifier_position_map,
        };
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{}.tpi", name))?);
        bincode::serialize_into(&mut writer, &versioned_id_pos_map)?;
    }

    info!("Add reverse complement and sentinels to reference");
    let ref_seq_rev_compl = dna::revcomp(&ref_seq);
    ref_seq.extend_from_slice(b"$");
    ref_seq.extend_from_slice(&ref_seq_rev_compl);
    drop(ref_seq_rev_compl);
    ref_seq.extend_from_slice(b"$");

    info!("Compress reference");
    alphabet.insert(b'$');
    let rank_transform = RankTransform::new(&alphabet);
    let alphabet = Alphabet::new(0..rank_transform.ranks.len() as u8);
    let ref_seq = rank_transform.transform(ref_seq);

    {
        info!("Save \"RT\" table to disk");
        let versioned_rt = VersionedRt {
            version: INDEX_VERSION,
            data: rank_transform,
        };
        let mut writer_rank_transform =
            snap::write::FrameEncoder::new(File::create(format!("{}.trt", name))?);
        bincode::serialize_into(&mut writer_rank_transform, &versioned_rt)?;
    }

    info!("Generate suffix array");
    let suffix_array = suffix_array(&ref_seq);

    info!("Generate BWT");
    let bwt = bwt(&ref_seq, &suffix_array);

    {
        info!("Save suffix array to disk");
        let versioned_suffix_array = VersionedSuffixArray {
            version: INDEX_VERSION,
            data: suffix_array,
        };
        let mut writer_suffix_array =
            snap::write::FrameEncoder::new(File::create(format!("{}.tsa", name))?);
        bincode::serialize_into(&mut writer_suffix_array, &versioned_suffix_array)?;
    }

    info!("Generate \"C\" table");
    let less = less(&bwt, &alphabet);

    info!("Generate \"Occ\" table");
    let occ = Occ::new(&bwt, 128, &alphabet);

    {
        info!("Save BWT to disk");
        let versioned_bwt = VersionedBwt {
            version: INDEX_VERSION,
            data: bwt,
        };
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{}.tbw", name))?);
        bincode::serialize_into(&mut writer, &versioned_bwt)?;
    }

    {
        info!("Save \"C\" table to disk");
        let versioned_less = VersionedLess {
            version: INDEX_VERSION,
            data: less,
        };
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{}.tle", name))?);
        bincode::serialize_into(&mut writer, &versioned_less)?;
    }

    {
        info!("Save \"Occ\" table to disk");
        let versioned_occ = VersionedOcc {
            version: INDEX_VERSION,
            data: occ,
        };
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{}.toc", name))?);
        bincode::serialize_into(&mut writer, &versioned_occ)?;
    }

    Ok(())
}

/// Calls different closures on each symbol depending on whether this symbol sits inside or outside of a run of equal symbols
fn run_apply<T, U>(ref_seq: &mut [u8], min_run_len: usize, mut non_run_fun: T, mut run_fun: U)
where
    T: FnMut(u8) -> u8,
    U: FnMut(u8) -> u8,
{
    assert!(min_run_len > 0);

    let mut run_symbol = ref_seq[0];
    let mut run_length = 1;

    let mut i = 1;
    while i < ref_seq.len() {
        while ref_seq[i] == run_symbol {
            run_length += 1;
            if i >= ref_seq.len() - 1 {
                let slice = &mut ref_seq[i - run_length..=i];
                if run_length < min_run_len {
                    for mutable_symbol in slice.iter_mut() {
                        *mutable_symbol = non_run_fun(*mutable_symbol);
                    }
                } else {
                    for mutable_symbol in slice.iter_mut() {
                        *mutable_symbol = run_fun(*mutable_symbol);
                    }
                }
                break;
            }
            i += 1;
        }
        run_symbol = ref_seq[i];
        let slice = &mut ref_seq[i - run_length..i];
        if run_length < min_run_len {
            for mutable_symbol in slice.iter_mut() {
                *mutable_symbol = non_run_fun(*mutable_symbol);
            }
        } else {
            for mutable_symbol in slice.iter_mut() {
                *mutable_symbol = run_fun(*mutable_symbol);
            }
        }

        i += 1;
        run_length = 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_replacement() {
        let ref_seq = "NNGATNTACANGATTNNACANNN".as_bytes().to_owned();

        let test_randomly_replace_ambiguous = |base| {
            if !DNA_UPPERCASE_ALPHABET.contains(&base) {
                b'A'
            } else {
                base
            }
        };
        let test_unify_ambiguous_symbols = |base| {
            if !DNA_UPPERCASE_ALPHABET.contains(&base) {
                b'X'
            } else {
                base
            }
        };

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                1,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"XXGATXTACAXGATTXXACAXXX");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                2,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"XXGATATACAAGATTXXACAXXX");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                3,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"AAGATATACAAGATTAAACAXXX");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                4,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"AAGATATACAAGATTAAACAAAA");
        }

        // Terminal singleton
        let ref_seq = "NNGATNTACANGATTNNACANNNT".as_bytes().to_owned();

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                1,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"XXGATXTACAXGATTXXACAXXXT");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                2,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"XXGATATACAAGATTXXACAXXXT");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                3,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"AAGATATACAAGATTAAACAXXXT");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                4,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"AAGATATACAAGATTAAACAAAAT");
        }

        // 5-terminal unambiguous symbol
        let ref_seq = "GNNGATNTACANGATTNNACANNNT".as_bytes().to_owned();

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                1,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"GXXGATXTACAXGATTXXACAXXXT");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                2,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"GXXGATATACAAGATTXXACAXXXT");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                3,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"GAAGATATACAAGATTAAACAXXXT");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                4,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"GAAGATATACAAGATTAAACAAAAT");
        }

        // Different ambiguous base symbol
        let ref_seq = "GNNGATNTACANGATYYYYYTNNACANNNT".as_bytes().to_owned();

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                1,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"GXXGATXTACAXGATXXXXXTXXACAXXXT");
        }

        let ref_seq = "CYNTYYNNT".as_bytes().to_owned();

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                2,
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"CAATXXXXT");
        }
    }

    #[test]
    fn char_x() {
        assert_eq!(dna::revcomp(b"GATTXACA"), "TGTXAATC".as_bytes());
    }
}
