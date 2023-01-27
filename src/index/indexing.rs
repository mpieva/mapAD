use std::{collections::BTreeMap, fs::File, iter, num::NonZeroUsize};

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
    errors::{Error, Result},
    index::{
        versioned_index::Item, FastaIdPosition, FastaIdPositions, OriginalSymbols,
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

    info!("Validate reference sequence");
    // Check if reference only contains valid IUPAC base symbols
    if !dna::iupac_alphabet().is_word(&ref_seq) {
        return Err(Error::ParseError(
            "Found non-IUPAC symbol in reference sequence".into(),
        ));
    }

    // Replace occurrences of ambiguous IUPAC DNA base symbols with random ones.
    // Calling `unwrap()` here will never panic because the slices to `choose()` from are `const`.
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
        _ => unreachable!(),
    };

    // Unconditionally return `X`
    let summarize_ambiguous = |_base| b'X';

    {
        info!("Modify reference sequence");
        // Start by replacing ambiguous bases that do not occur in stretches exceeding a certain length,
        // using the closures defined above
        let original_symbols = run_apply(
            &mut ref_seq,
            20.try_into().expect("number to be non-zero"),
            randomly_replace_ambiguous,
            summarize_ambiguous,
        );

        info!("Save original symbols");
        let versioned_original_symbols = Item::new(original_symbols);
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{name}.tos"))?);
        bincode::serialize_into(&mut writer, &versioned_original_symbols)?;
    }

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
        let versioned_id_pos_map = Item::new(identifier_position_map);
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{name}.tpi"))?);
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
    let alphabet = Alphabet::new(
        0..u8::try_from(rank_transform.ranks.len()).expect("alphabet size to be < 256"),
    );
    let ref_seq = rank_transform.transform(ref_seq);

    {
        info!("Save \"RT\" table");
        let versioned_rt = Item::new(rank_transform);
        let mut writer_rank_transform =
            snap::write::FrameEncoder::new(File::create(format!("{name}.trt"))?);
        bincode::serialize_into(&mut writer_rank_transform, &versioned_rt)?;
    }

    info!("Generate suffix array");
    let suffix_array = suffix_array(&ref_seq);

    info!("Generate BWT");
    let bwt = bwt(&ref_seq, &suffix_array);

    {
        info!("Compress suffix array");
        let owned_sampled_suffix_array = SampledSuffixArrayOwned::sample(
            &suffix_array,
            &ref_seq,
            &bwt,
            32.try_into().expect("to be non-zero"),
        );

        info!("Save compressed suffix array");
        let versioned_suffix_array = Item::new(owned_sampled_suffix_array);
        let mut writer_suffix_array =
            snap::write::FrameEncoder::new(File::create(format!("{name}.tsa"))?);
        bincode::serialize_into(&mut writer_suffix_array, &versioned_suffix_array)?;
    }

    info!("Generate \"C\" table");
    let less = less(&bwt, &alphabet);

    info!("Generate \"Occ\" table");
    let occ = Occ::new(&bwt, 128, &alphabet);

    {
        info!("Save BWT");
        let versioned_bwt = Item::new(bwt);
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{name}.tbw"))?);
        bincode::serialize_into(&mut writer, &versioned_bwt)?;
    }

    {
        info!("Save \"C\" table");
        let versioned_less = Item::new(less);
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{name}.tle"))?);
        bincode::serialize_into(&mut writer, &versioned_less)?;
    }

    {
        info!("Save \"Occ\" table");
        let versioned_occ = Item::new(occ);
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{name}.toc"))?);
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
) -> OriginalSymbols
where
    T: FnMut(u8) -> u8,
    U: FnMut(u8) -> u8,
{
    let mut original_symbols = BTreeMap::new();
    let mut i = 0;
    while let Some(&symbol_i) = ref_seq.get(i) {
        // `&array[array.len()..]` works and yields an empty slice
        let run_len = ref_seq[i + 1..]
            .iter()
            .enumerate()
            .find(|&(_j, &symbol_j)| symbol_j != symbol_i)
            // Correct length for starting at position i+1
            .map_or(ref_seq.len() - i, |(j, _symbol_j)| j + 1);
        // Apply closures on runs of non-alphabet symbols
        if !DNA_UPPERCASE_ALPHABET.contains(&symbol_i) {
            let run = &mut ref_seq[i..][..run_len];
            // Short run
            if run_len < min_run_len.into() {
                for (j, base) in run.iter_mut().enumerate() {
                    let new_symbol_tmp = non_run_fun(*base);
                    let prev_value = original_symbols.insert(i + j, *base);
                    assert!(prev_value.is_none());
                    *base = new_symbol_tmp;
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
    OriginalSymbols::new(original_symbols)
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
            let mut ref_seq = ref_seq;
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
            let mut ref_seq = ref_seq;
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
            let mut ref_seq = ref_seq;
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
            let mut ref_seq = ref_seq;
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
            let mut ref_seq = ref_seq;
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
