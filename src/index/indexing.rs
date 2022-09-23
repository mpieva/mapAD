use std::{fs::File, iter, num::NonZeroUsize};

use bio::{
    alphabets::{dna, Alphabet, RankTransform},
    data_structures::{
        bwt::{bwt, less, Occ},
        suffix_array::suffix_array,
    },
    io::fasta,
};
use either::Either;
use log::info;
use rand::{
    prelude::{Rng, SeedableRng, StdRng},
    seq::SliceRandom,
};

use crate::{
    errors::Result,
    index::{
        versioned_index::VersionedIndexItem, FastaIdPosition, FastaIdPositions,
        SampledSuffixArrayOwned, DNA_AMINO, DNA_KETONE, DNA_NOT_A, DNA_NOT_C, DNA_NOT_G, DNA_NOT_T,
        DNA_PURINE, DNA_PYRIMIDINE, DNA_STRONG, DNA_UPPERCASE_ALPHABET, DNA_UPPERCASE_X_ALPHABET,
        DNA_WEAK,
    },
};

/// Entry point function to launch the indexing process
pub fn run(reference_path: &str, seed: u64) -> Result<()> {
    let mut rng: StdRng = SeedableRng::seed_from_u64(seed);

    let alphabet = Alphabet::new(DNA_UPPERCASE_X_ALPHABET);

    // Index the genome
    index(reference_path, alphabet, reference_path, &mut rng)?;

    Ok(())
}

/// Index a given reference and write the index.
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
        // `flat_map()` works by calling `into_iter` on the values returned by the closure, which
        // means that the `Result`s would get lost in the process. We use the following construct
        // to keep those.
        .flat_map(|record| match record {
            // Convert all bases to uppercase
            Ok(record) => Either::Left(record.seq().to_ascii_uppercase().into_iter().map(Ok)),
            Err(e) => Either::Right(iter::once(Err(e.into()))),
        })
        .collect::<Result<Vec<_>>>()?;

    // Replace single occurrences of ambiguous base symbols with random ones, leave runs alone
    // Calling `unwrap()` here will never panic because the slices to `choose()` from are
    // `const`.
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
        10.try_into().expect("number to be non-zero"),
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
                    end += record.seq().len() as u64;
                    FastaIdPosition {
                        start: end - record.seq().len() as u64,
                        end: end - 1,
                        identifier: record.id().to_string(),
                    }
                })
                .collect::<Vec<_>>(),
        );

        info!("Save position map");
        let versioned_id_pos_map = VersionedIndexItem::new(identifier_position_map);
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
        info!("Save \"RT\" table");
        let versioned_rt = VersionedIndexItem::new(rank_transform);
        let mut writer_rank_transform =
            snap::write::FrameEncoder::new(File::create(format!("{}.trt", name))?);
        bincode::serialize_into(&mut writer_rank_transform, &versioned_rt)?;
    }

    info!("Generate suffix array");
    let suffix_array = suffix_array(&ref_seq);

    info!("Generate BWT");
    let bwt = bwt(&ref_seq, &suffix_array);

    {
        info!("Compress suffix array");
        let owned_sampled_suffix_array =
            SampledSuffixArrayOwned::sample(&suffix_array, &ref_seq, &bwt, 32);

        info!("Save compressed suffix array");
        let versioned_suffix_array = VersionedIndexItem::new(owned_sampled_suffix_array);
        let mut writer_suffix_array =
            snap::write::FrameEncoder::new(File::create(format!("{}.tsa", name))?);
        bincode::serialize_into(&mut writer_suffix_array, &versioned_suffix_array)?;
    }

    info!("Generate \"C\" table");
    let less = less(&bwt, &alphabet);

    info!("Generate \"Occ\" table");
    let occ = Occ::new(&bwt, 128, &alphabet);

    {
        info!("Save BWT");
        let versioned_bwt = VersionedIndexItem::new(bwt);
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{}.tbw", name))?);
        bincode::serialize_into(&mut writer, &versioned_bwt)?;
    }

    {
        info!("Save \"C\" table");
        let versioned_less = VersionedIndexItem::new(less);
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{}.tle", name))?);
        bincode::serialize_into(&mut writer, &versioned_less)?;
    }

    {
        info!("Save \"Occ\" table");
        let versioned_occ = VersionedIndexItem::new(occ);
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{}.toc", name))?);
        bincode::serialize_into(&mut writer, &versioned_occ)?;
    }

    Ok(())
}

/// Calls different closures on each symbol depending on whether this symbol sits inside or outside of a run of equal symbols
fn run_apply<T, U>(
    ref_seq: &mut [u8],
    min_run_len: NonZeroUsize,
    mut non_run_fun: T,
    mut run_fun: U,
) where
    T: FnMut(u8) -> u8,
    U: FnMut(u8) -> u8,
{
    let mut i = 0;
    while i < ref_seq.len() - 1 {
        let symbol_i = ref_seq[i];
        let run_len = ref_seq[i + 1..]
            .iter()
            .enumerate()
            .find(|&(_j, &symbol_j)| symbol_j != symbol_i)
            // Correct length for starting at position i+1
            .map(|(j, _symbol_j)| j + 1)
            .unwrap_or(ref_seq.len() - i);
        // Apply closures on runs of non-alphabet symbols
        if !DNA_UPPERCASE_ALPHABET.contains(&symbol_i) {
            let run = &mut ref_seq[i..][..run_len];
            // Short run
            if run_len < min_run_len.into() {
                for base in run.iter_mut() {
                    *base = non_run_fun(*base);
                }
            // Long run
            } else {
                for base in run.iter_mut() {
                    *base = run_fun(*base);
                }
            }
        }
        i += run_len;
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
                1.try_into().unwrap(),
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"XXGATXTACAXGATTXXACAXXX");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                2.try_into().unwrap(),
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"XXGATATACAAGATTXXACAXXX");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                3.try_into().unwrap(),
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"AAGATATACAAGATTAAACAXXX");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                4.try_into().unwrap(),
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
                1.try_into().unwrap(),
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"XXGATXTACAXGATTXXACAXXXT");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                2.try_into().unwrap(),
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"XXGATATACAAGATTXXACAXXXT");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                3.try_into().unwrap(),
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"AAGATATACAAGATTAAACAXXXT");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                4.try_into().unwrap(),
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
                1.try_into().unwrap(),
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"GXXGATXTACAXGATTXXACAXXXT");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                2.try_into().unwrap(),
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"GXXGATATACAAGATTXXACAXXXT");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                3.try_into().unwrap(),
                test_randomly_replace_ambiguous,
                test_unify_ambiguous_symbols,
            );
            assert_eq!(&ref_seq, b"GAAGATATACAAGATTAAACAXXXT");
        }

        {
            let mut ref_seq = ref_seq.clone();
            run_apply(
                &mut ref_seq,
                4.try_into().unwrap(),
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
                1.try_into().unwrap(),
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
                2.try_into().unwrap(),
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
