use bio::{
    alphabets::{dna, RankTransform},
    data_structures::{
        bwt::{Less, Occ, BWT},
        fmindex::{FMIndexable, Interval},
        suffix_array::SuffixArray,
    },
};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct RtFMDIndex {
    bwt: BWT,
    less: Less,
    occ: Occ,
    rank_transform: RankTransform,
    back_transform: Vec<u8>,
}

impl FMIndexable for RtFMDIndex {
    fn occ(&self, r: usize, a: u8) -> usize {
        self.occ.get(&self.bwt, r, a)
    }
    fn less(&self, a: u8) -> usize {
        self.less[a as usize]
    }
    /// Provide a reference to the underlying BWT.
    fn bwt(&self) -> &BWT {
        &self.bwt
    }
}

/// FMD-Index (Li, 2012) that operates on ranks of symbols instead of their ASCII representations
impl RtFMDIndex {
    pub fn new(bwt: BWT, less: Less, occ: Occ, rank_transform: RankTransform) -> Self {
        let mut back_transform = rank_transform
            .ranks
            .keys()
            .map(|symbol| symbol as u8)
            .collect::<Vec<_>>();
        back_transform.sort_unstable();

        Self {
            bwt,
            less,
            occ,
            rank_transform,
            back_transform,
        }
    }

    /// Initialize interval for empty pattern. The interval points at the whole suffix array.
    pub fn init_interval(&self) -> RtBiInterval {
        RtBiInterval {
            lower: 0,
            lower_rev: 0,
            size: self.bwt.len(),
            match_size: 0,
        }
    }

    /// Backward extension of given interval with given character.
    /// Takes plain, non-transformed symbol.
    pub fn backward_ext(&self, interval: &RtBiInterval, a: u8) -> RtBiInterval {
        // Non-alphabet character
        if !self.rank_transform.ranks.contains_key(a as usize) {
            return RtBiInterval {
                lower: 0,
                lower_rev: 0,
                size: 0,
                match_size: 0,
            };
        }

        self.extend_iter(interval)
            .find(|&(base, _)| base == self.rank_transform.get(a))
            .map(|(_, interval)| interval)
            .unwrap()
    }

    pub fn forward_ext(&self, interval: &RtBiInterval, a: u8) -> RtBiInterval {
        self.backward_ext(&interval.swapped(), dna::complement(a))
            .swapped()
    }

    /// Returns an Iterator over the alphabet to extend the RtBiInterval
    pub fn extend_iter<'a>(&'a self, interval: &'a RtBiInterval) -> FMDExtIterator<'a> {
        FMDExtIterator::new(interval, self)
    }

    pub fn get_rev(&self, a: u8) -> u8 {
        self.back_transform[a as usize]
    }
}

/// Extension of a RtBiInterval, implemented as an Iterator over the alphabet
pub struct FMDExtIterator<'a> {
    s: usize,
    l: usize,
    c: u8,
    input_interval: &'a RtBiInterval,
    fmd_index: &'a RtFMDIndex,
}

impl<'a> Iterator for FMDExtIterator<'a> {
    type Item = (u8, RtBiInterval);

    fn next(&mut self) -> Option<Self::Item> {
        if self.c < 2 {
            return None;
        }
        self.c -= 1;
        Some((self.c, self.extend_once_internal()))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let hint = self.c as usize - 1;
        (hint, Some(hint))
    }
}

impl<'a> ExactSizeIterator for FMDExtIterator<'a> {}

impl<'a> FMDExtIterator<'a> {
    fn new(interval: &'a RtBiInterval, fmd_index: &'a RtFMDIndex) -> Self {
        let mut initial_state = Self {
            s: 0,
            l: interval.lower_rev,
            c: 0,
            input_interval: interval,
            fmd_index,
        };
        initial_state.extend_once_internal();
        initial_state.c = 5;
        initial_state
    }

    fn extend_once_internal(&mut self) -> RtBiInterval {
        self.l += self.s;
        let o = if self.input_interval.lower == 0 {
            0
        } else {
            self.fmd_index.occ(self.input_interval.lower - 1, self.c)
        };

        // Interval size I^s
        self.s = self.fmd_index.occ(
            self.input_interval.lower + self.input_interval.size - 1,
            self.c,
        ) - o;

        RtBiInterval {
            lower: self.fmd_index.less(self.c) + o,
            lower_rev: self.l,
            size: self.s,
            match_size: self.input_interval.match_size + 1,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RtBiInterval {
    pub lower: usize,
    pub lower_rev: usize,
    pub size: usize,
    pub match_size: usize,
}

impl RtBiInterval {
    pub fn forward(&self) -> Interval {
        Interval {
            upper: self.lower + self.size,
            lower: self.lower,
        }
    }

    fn interval_to_pos<'a, S>(
        interval: Interval,
        suffix_array: &'a S,
    ) -> impl Iterator<Item = usize> + 'a
    where
        S: SuffixArray,
    {
        (interval.lower..interval.upper).map(move |pos| {
            suffix_array
                .get(pos)
                .expect("Interval out of range of suffix array")
        })
    }

    pub fn revcomp(&self) -> Interval {
        Interval {
            upper: self.lower_rev + self.size,
            lower: self.lower_rev,
        }
    }

    pub fn occ_fwd<'a, S>(&self, suffix_array: &'a S) -> impl Iterator<Item = usize> + 'a
    where
        S: SuffixArray,
    {
        Self::interval_to_pos(self.forward(), suffix_array)
    }

    pub fn occ_revcomp<'a, S>(&self, suffix_array: &'a S) -> impl Iterator<Item = usize> + 'a
    where
        S: SuffixArray,
    {
        Self::interval_to_pos(self.revcomp(), suffix_array)
    }

    pub fn swapped(&self) -> Self {
        Self {
            lower: self.lower_rev,
            lower_rev: self.lower,
            size: self.size,
            match_size: self.match_size,
        }
    }
}
