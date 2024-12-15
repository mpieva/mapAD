use std::ops::Range;

use bio::{
    alphabets::{dna, RankTransform},
    data_structures::{
        bwt::{Less, Occ, BWT},
        fmindex::{FMIndexable, Interval},
    },
};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct RtFmdIndex {
    pub bwt: BWT,
    pub less: Less,
    pub occ: Occ,
    sentinel_occ: [usize; 2],
    rank_transform: RankTransform,
    back_transform: Vec<u8>,
}

impl FMIndexable for RtFmdIndex {
    fn occ(&self, r: usize, a: u8) -> usize {
        self.occ.get_small_k(&self.bwt, r, a)
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
impl RtFmdIndex {
    pub fn new(bwt: BWT, less: Less, occ: Occ, rank_transform: RankTransform) -> Self {
        let mut sentinel_occ = [0; 2];
        sentinel_occ
            .iter_mut()
            .zip(
                bwt.iter()
                    .enumerate()
                    .filter(|(_i, &symbol)| symbol == 0)
                    .map(|(i, _symbol)| i),
            )
            .for_each(|(target, src)| *target = src);

        let mut back_transform = rank_transform
            .ranks
            .keys()
            .map(|symbol| u8::try_from(symbol).expect("alphabet size to be < 256"))
            .collect::<Vec<_>>();
        back_transform.sort_unstable();

        Self {
            bwt,
            less,
            occ,
            sentinel_occ,
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
            };
        }

        self.extend_iter(interval)
            .find(|&(base, _)| base == self.rank_transform.get(a))
            .map(|(_, interval)| interval)
            .expect("This is not expected to fail")
    }

    pub fn forward_ext(&self, interval: &RtBiInterval, a: u8) -> RtBiInterval {
        self.backward_ext(&interval.swapped(), dna::complement(a))
            .swapped()
    }

    /// Returns an Iterator over the alphabet to extend the `RtBiInterval`
    pub fn extend_iter<'a>(&'a self, interval: &'a RtBiInterval) -> FmdExtIterator<'a> {
        FmdExtIterator::new(interval, self)
    }

    pub fn get_rev(&self, a: u8) -> u8 {
        self.back_transform[a as usize]
    }
}

/// Extension of a `RtBiInterval`, implemented as an Iterator over the alphabet
pub struct FmdExtIterator<'a> {
    s: usize,
    l: usize,
    c: u8,
    input_interval: &'a RtBiInterval,
    fmd_index: &'a RtFmdIndex,
}

impl Iterator for FmdExtIterator<'_> {
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

impl ExactSizeIterator for FmdExtIterator<'_> {}

impl<'a> FmdExtIterator<'a> {
    fn new(interval: &'a RtBiInterval, fmd_index: &'a RtFmdIndex) -> Self {
        /// Finds occurrences of '$' in BWT. Since there are only two occurrences, we can cheaply
        /// store and get both positions from a special cache (avoiding the `Occ` array).
        fn sentinel_occ(interval_pos: usize, fmd_index: &RtFmdIndex) -> usize {
            fmd_index
                .sentinel_occ
                .iter()
                .position(|&sentinel_pos| interval_pos < sentinel_pos)
                .unwrap_or(fmd_index.sentinel_occ.len())
        }
        let o = if interval.lower == 0 {
            0
        } else {
            sentinel_occ(interval.lower - 1, fmd_index)
        };

        Self {
            s: sentinel_occ(interval.lower + interval.size - 1, fmd_index) - o,
            l: interval.lower_rev,
            c: 5,
            input_interval: interval,
            fmd_index,
        }
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
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RtBiInterval {
    pub lower: usize,
    pub lower_rev: usize,
    pub size: usize,
}

impl RtBiInterval {
    pub fn forward(&self) -> Interval {
        Interval {
            upper: self.lower + self.size,
            lower: self.lower,
        }
    }

    pub fn revcomp(&self) -> Interval {
        Interval {
            upper: self.lower_rev + self.size,
            lower: self.lower_rev,
        }
    }

    #[must_use]
    pub fn swapped(&self) -> Self {
        Self {
            lower: self.lower_rev,
            lower_rev: self.lower,
            size: self.size,
        }
    }

    pub fn range_fwd(&self) -> Range<usize> {
        let interval = self.forward();
        interval.lower..interval.upper
    }
}
