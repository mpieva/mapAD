use std::{collections::HashMap, fs::File, hash::BuildHasherDefault};

use bio::{
    alphabets::{dna, Alphabet, RankTransform},
    data_structures::{
        bwt::{bwt, less, Less, Occ, BWT},
        suffix_array::{suffix_array, SuffixArray},
    },
    io::fasta,
};
use either::Either;
use fxhash::FxHasher;
use log::info;
use rand::{
    prelude::{Rng, SeedableRng, StdRng},
    seq::SliceRandom,
};
use serde::{de::DeserializeOwned, Deserialize, Serialize};

use crate::{
    errors::{Error, Result},
    map::{FastaIdPosition, FastaIdPositions},
};

// Increase this number once the on-disk index changes
pub const INDEX_VERSION: u8 = 2;

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

// HashMap using a fast hasher
type HashMapFx<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

/// Owned data of sampled suffix array. The borrowed parts need to be be
/// reconstructed after deserialization.
#[derive(Serialize, Deserialize)]
pub struct SampledSuffixArrayOwned {
    sample: Vec<usize>,
    s: usize, // Rate of sampling
    extra_rows: HashMapFx<usize, usize>,
    sentinel: u8,
}

impl SampledSuffixArrayOwned {
    /// Sample the suffix array with the given sample rate.
    /// This is copied from the `bio` crate because we need more serde flexibility.
    fn sample<S>(suffix_array: &S, text: &[u8], bwt: &BWT, sampling_rate: usize) -> Self
    where
        S: SuffixArray,
    {
        let mut sample =
            Vec::with_capacity((suffix_array.len() as f32 / sampling_rate as f32).ceil() as usize);
        let mut extra_rows = HashMapFx::default();
        let sentinel = *text
            .last()
            .expect("The text should not be empty at this point");

        for (i, l_row) in bwt.iter().copied().enumerate() {
            let idx = suffix_array
                .get(i)
                .expect("BWT and suffix array have the same length");
            if (i % sampling_rate) == 0 {
                sample.push(idx);
            } else if l_row == sentinel {
                // If bwt lookup will return a sentinel
                // Text suffixes that begin right after a sentinel are always saved as extra rows
                // to help deal with FM index last to front inaccuracy when there are many sentinels
                extra_rows.insert(i, idx);
            }
        }

        Self {
            sample,
            s: sampling_rate,
            extra_rows,
            sentinel,
        }
    }

    pub fn into_sampled_suffix_array<'a, 'b, 'c>(
        self,
        bwt: &'a BWT,
        less: &'b Less,
        occ: &'c Occ,
    ) -> SampledSuffixArray<'a, 'b, 'c> {
        SampledSuffixArray {
            bwt,
            less,
            occ,
            sample: self.sample,
            s: self.s,
            extra_rows: self.extra_rows,
            sentinel: self.sentinel,
        }
    }
}

/// A sampled suffix array. The code is copied from the `bio`
/// crate because we need access to private fields.
pub struct SampledSuffixArray<'a, 'b, 'c> {
    bwt: &'a BWT,
    less: &'b Less,
    occ: &'c Occ,
    sample: Vec<usize>,
    s: usize, // Rate of sampling
    extra_rows: HashMapFx<usize, usize>,
    sentinel: u8,
}

impl<'a, 'b, 'c> SuffixArray for SampledSuffixArray<'a, 'b, 'c> {
    fn get(&self, index: usize) -> Option<usize> {
        if index < self.len() {
            let mut pos = index;
            let mut offset = 0;
            loop {
                if pos % self.s == 0 {
                    return Some(self.sample[pos / self.s] + offset);
                }

                let c = self.bwt[pos];

                if c == self.sentinel {
                    // Check if next character in the bwt is the sentinel
                    // If so, there must be a cached result to workaround FM index last to front
                    // mapping inaccuracy when there are multiple sentinels
                    // This branch should rarely be triggered so the performance impact
                    // of hashmap lookups would be low
                    return Some(self.extra_rows[&pos] + offset);
                }

                pos = self.less[c as usize] + self.occ.get(self.bwt, pos - 1, c);
                offset += 1;
            }
        } else {
            None
        }
    }

    fn len(&self) -> usize {
        self.bwt.len()
    }

    fn is_empty(&self) -> bool {
        self.bwt.is_empty()
    }
}

// Versioned index data structures
#[derive(Serialize, Deserialize)]
pub struct VersionedIndexItem<T> {
    version: u8,
    data: T,
}

impl<T> VersionedIndexItem<T> {
    pub fn new(data: T) -> Self {
        Self {
            version: INDEX_VERSION,
            data,
        }
    }

    /// Returns inner data if the version of the deserialized Struct is compatible
    pub fn try_take(self) -> Result<T> {
        if self.version == INDEX_VERSION {
            Ok(self.data)
        } else {
            Err(Error::IndexVersionMismatch)
        }
    }
}

impl<T> VersionedIndexItem<T>
where
    T: DeserializeOwned,
{
    pub fn read_from_bincode<R>(reader: R) -> Result<Self>
    where
        R: std::io::Read,
    {
        bincode::deserialize_from::<_, Self>(reader).map_err(|e| e.into())
    }
}

/// Entry point function to launch the indexing process
pub fn run(reference_path: &str, seed: u64) -> Result<()> {
    let mut rng: StdRng = SeedableRng::seed_from_u64(seed);

    let alphabet = Alphabet::new(crate::index::DNA_UPPERCASE_X_ALPHABET);

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
        // Convert all bases to uppercase
        .flat_map(|record| match record {
            Ok(record) => Either::Left(record.seq().to_ascii_uppercase().into_iter().map(Ok)),
            Err(e) => Either::Right(std::iter::once(e).map(|e| Err(e.into()))),
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
        let versioned_bwt = VersionedIndexItem::<BWT> {
            version: INDEX_VERSION,
            data: bwt,
        };
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
