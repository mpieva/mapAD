pub mod indexing;

mod versioned_index;

use std::{collections::BTreeMap, num::NonZeroUsize};

use bio::data_structures::{
    bwt::{Less, Occ, BWT},
    suffix_array::SuffixArray,
};
use log::debug;
use serde::{Deserialize, Serialize};

use crate::{errors::Result, index::versioned_index::Item, map::fmd_index::RtFmdIndex};

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

/// For multi-identifier reference sequences like the human genome (that is split by chromosome)
/// this struct is used to keep a map of IDs and positions
#[derive(Serialize, Deserialize, Debug)]
pub struct FastaIdPosition {
    pub start: u64,
    pub end: u64,
    pub identifier: String,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct FastaIdPositions {
    id_position: Vec<FastaIdPosition>,
}

impl FastaIdPositions {
    pub fn new(id_position: Vec<FastaIdPosition>) -> Self {
        Self { id_position }
    }

    pub fn iter(&self) -> impl Iterator<Item = &FastaIdPosition> {
        self.id_position.iter()
    }

    /// Find the corresponding reference identifier by position. The function
    /// returns a tuple: ("target ID", "relative position")
    pub fn get_reference_identifier(
        &self,
        position: usize,
        pattern_length: usize,
    ) -> Option<(u32, u64, &str)> {
        let position = position as u64;
        self.id_position
            .iter()
            .enumerate()
            .find(|(_, identifier)| {
                (identifier.start <= position)
                    && (position + pattern_length as u64 - 1 <= identifier.end)
            })
            .and_then(|(index, identifier)| {
                Some((
                    u32::try_from(index).ok()?,
                    position - identifier.start,
                    identifier.identifier.as_str(),
                ))
            })
    }
}

/// Owned data of sampled suffix array. The borrowed parts need to be be
/// reconstructed after deserialization.
#[derive(Serialize, Deserialize)]
pub struct SampledSuffixArrayOwned {
    sample: Vec<usize>,
    sampling_rate: NonZeroUsize,
    extra_rows: BTreeMap<usize, usize>,
    sentinel: u8,
}

impl SampledSuffixArrayOwned {
    /// Sample the suffix array with the given sample rate.
    /// This is copied from the `bio` crate because we need more serde flexibility.
    pub fn sample<S>(suffix_array: &S, text: &[u8], bwt: &BWT, sampling_rate: NonZeroUsize) -> Self
    where
        S: SuffixArray,
    {
        let mut sample = {
            let mut vec_cap = suffix_array.len() / sampling_rate;
            // Ceil
            if suffix_array.len() % sampling_rate != 0 {
                vec_cap += 1;
            }
            Vec::with_capacity(vec_cap)
        };
        let mut extra_rows = BTreeMap::default();
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
            sampling_rate,
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
            sampling_rate: self.sampling_rate,
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
    sampling_rate: NonZeroUsize,
    extra_rows: BTreeMap<usize, usize>,
    sentinel: u8,
}

impl SuffixArray for SampledSuffixArray<'_, '_, '_> {
    fn get(&self, index: usize) -> Option<usize> {
        if index < self.len() {
            let mut pos = index;
            let mut offset = 0;
            loop {
                if pos % self.sampling_rate == 0 {
                    return Some(self.sample[pos / self.sampling_rate] + offset);
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

#[derive(Debug, Serialize, Deserialize)]
pub struct OriginalSymbols(BTreeMap<usize, u8>);

impl OriginalSymbols {
    pub fn new(map_data: BTreeMap<usize, u8>) -> Self {
        Self(map_data)
    }

    /// Returns the original symbol at a _forward strand_ position
    pub fn get(&self, idx: usize) -> Option<u8> {
        self.0.get(&idx).copied()
    }
}

pub fn load_suffix_array_from_path(reference_path: &str) -> Result<SampledSuffixArrayOwned> {
    Item::read_from_path(format!("{reference_path}.tsa"))?.try_take()
}

pub fn load_id_pos_map_from_path(reference_path: &str) -> Result<FastaIdPositions> {
    Item::read_from_path(format!("{reference_path}.tpi"))?.try_take()
}

pub fn load_original_symbols_from_path(reference_path: &str) -> Result<OriginalSymbols> {
    Item::read_from_path(format!("{reference_path}.tos"))?.try_take()
}

pub fn load_index_from_path(reference_path: &str) -> Result<RtFmdIndex> {
    debug!("Load BWT");
    let bwt = Item::read_from_path(format!("{reference_path}.tbw"))?.try_take()?;

    debug!("Load \"C\" table");
    let less = Item::read_from_path(format!("{reference_path}.tle"))?.try_take()?;

    debug!("Load \"Occ\" table");
    let occ = Item::read_from_path(format!("{reference_path}.toc"))?.try_take()?;

    debug!("Load \"RT\" table");
    let rt = Item::read_from_path(format!("{reference_path}.trt"))?.try_take()?;

    debug!("Reconstruct index");
    Ok(RtFmdIndex::new(bwt, less, occ, rt))
}
