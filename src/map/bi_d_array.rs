use std::iter;

use either::Either;
use smallvec::SmallVec;

use crate::map::{
    fmd_index::RtFmdIndex, sequence_difference_models::SequenceDifferenceModel,
    AlignmentParameters, Direction,
};

/// Compute the lower bound of mismatches of a read per position by aligning perfectly,
/// starting from the read end and progressing to the center. As soon as the extension of
/// the alignment is not possible, we found at least one mismatch and record that per
/// read-position in the so-called D-array. The content of the array is used as "lookahead"
/// score. This allows for pruning the search tree. The values are minimal expected penalties
/// towards the respective ends of the query. Like alignment scores, these values are negative.
#[derive(Debug)]
pub struct BiDArray {
    d_composite: Vec<f32>,
    split: usize,
}

impl BiDArray {
    pub fn new<SDM>(
        pattern: &[u8],
        base_qualities: &[u8],
        split: usize,
        alignment_parameters: &AlignmentParameters,
        fmd_index: &RtFmdIndex,
        sdm: &SDM,
    ) -> Self
    where
        SDM: SequenceDifferenceModel,
    {
        // This value is never actually used as `u16` but this type
        // restricts the values to the valid range.
        const MAX_OFFSET: u16 = 15;

        // Compute a backward-D-iterator for every offset
        let mut offset_d_backward_iterators = (0..MAX_OFFSET as usize)
            .map(|offset| {
                Self::compute_part(
                    &pattern[..split],
                    &base_qualities[..split],
                    Direction::Forward,
                    pattern.len(),
                    offset as i16,
                    alignment_parameters,
                    fmd_index,
                    sdm,
                )
            })
            .collect::<SmallVec<[_; MAX_OFFSET as usize]>>();

        // Construct backward-D-array from minimal (most severe) penalties for each position
        let d_backwards = (0..split).map(|_| {
            offset_d_backward_iterators
                .iter_mut()
                .map(|offset_d_rev_iter| {
                    offset_d_rev_iter
                        .next()
                        .expect("Iterator is guaranteed to have the right length")
                })
                .fold(0_f32, |min_penalty, penalty| min_penalty.min(penalty))
        });

        // Compute a forward-D-iterator for every offset
        let mut offset_d_forward_iterators = (0..MAX_OFFSET)
            .map(|offset| {
                Self::compute_part(
                    &pattern[split..],
                    &base_qualities[split..],
                    Direction::Backward,
                    pattern.len(),
                    offset as i16,
                    alignment_parameters,
                    fmd_index,
                    sdm,
                )
            })
            .collect::<SmallVec<[_; MAX_OFFSET as usize]>>();

        // Construct forward-D-array from minimal (most severe) penalties for each position
        let d_forwards = (0..pattern.len() - split).map(|_| {
            offset_d_forward_iterators
                .iter_mut()
                .map(|offset_d_fwd_iter| {
                    offset_d_fwd_iter
                        .next()
                        .expect("Iterator is guaranteed to have the right length")
                })
                .fold(0_f32, |min_penalty, penalty| min_penalty.min(penalty))
        });

        Self {
            d_composite: d_backwards.chain(d_forwards).collect(),
            split,
        }
    }

    /// Computes either left or right part of the Bi-D-Array for the given pattern.
    /// `Direction` here is the opposite of the read alignment direction in `k_mismatch_search`.
    #[allow(clippy::too_many_arguments)]
    fn compute_part<'a, SDM>(
        pattern_part: &'a [u8],
        base_qualities_part: &'a [u8],
        direction: Direction,
        full_pattern_length: usize,
        initial_skip: i16,
        alignment_parameters: &'a AlignmentParameters,
        fmd_index: &'a RtFmdIndex,
        sdm: &'a SDM,
    ) -> impl Iterator<Item = f32> + 'a
    where
        SDM: SequenceDifferenceModel,
    {
        fn directed_index(index: usize, length: usize, direction: Direction) -> usize {
            match direction {
                Direction::Forward => index,
                Direction::Backward => length - 1 - index,
            }
        }

        iter::repeat(0.0).take(initial_skip as usize).chain(
            match direction {
                Direction::Forward => Either::Left(pattern_part.iter()),
                Direction::Backward => Either::Right(pattern_part.iter().rev()),
            }
            .enumerate()
            .skip(initial_skip as usize)
            .scan(
                (0.0, initial_skip - 1, fmd_index.init_interval()),
                move |(z, last_mismatch_pos, interval), (index, &base)| {
                    *interval = match direction {
                        Direction::Forward => fmd_index.forward_ext(interval, base),
                        Direction::Backward => fmd_index.backward_ext(interval, base),
                    };
                    if interval.size < 1 {
                        // Sub-read does not align perfectly, scan sub-sequence to find the most conservative penalty
                        *z += match direction {
                            Direction::Forward => Either::Left(pattern_part.iter()),
                            Direction::Backward => Either::Right(pattern_part.iter().rev()),
                        }
                        .enumerate()
                        .take(index + 1)
                        .skip((*last_mismatch_pos + 1) as usize)
                        .map(|(j, &base_j)| {
                            let idx_mapped_to_read =
                                directed_index(j, full_pattern_length, direction);
                            let best_penalty_mm_only = sdm.get_min_penalty(
                                idx_mapped_to_read,
                                full_pattern_length,
                                base_j,
                                base_qualities_part
                                    [directed_index(j, base_qualities_part.len(), direction)],
                                true,
                            );
                            let optimal_penalty = sdm.get_min_penalty(
                                idx_mapped_to_read,
                                full_pattern_length,
                                base_j,
                                base_qualities_part
                                    [directed_index(j, base_qualities_part.len(), direction)],
                                false,
                            );
                            // The optimal penalty, conditioning on position and base, is subtracted because we
                            // optimize the ratio `AS/optimal_AS` (in log space) to find the best mappings
                            best_penalty_mm_only - optimal_penalty
                        })
                        .fold(f32::MIN, |acc, penalty| acc.max(penalty))
                        .max(
                            alignment_parameters.penalty_gap_open
                                + alignment_parameters.penalty_gap_extend,
                        )
                        .max(alignment_parameters.penalty_gap_extend);
                        *interval = fmd_index.init_interval();
                        *last_mismatch_pos = index as i16;
                    }
                    Some(*z)
                },
            ),
        )
    }

    pub fn get(&self, backward_index: i16, forward_index: i16) -> f32 {
        let d_rev = if backward_index < 0 {
            0.0
        } else {
            self.d_composite[backward_index as usize]
        };
        let d_fwd =
            // Cover case forward_index == self.d_composite.len()
            if forward_index as usize >= self.d_composite.len() {
                0.0
            } else {
                self.d_composite[self.d_composite.len() - 1 - forward_index as usize + self.split]
            };
        d_rev + d_fwd
    }
}

#[cfg(test)]
mod tests {
    use bio::alphabets;

    use crate::{
        index::DNA_UPPERCASE_ALPHABET,
        map::{
            bi_d_array::BiDArray,
            mismatch_bounds::TestBound,
            sequence_difference_models::{
                LibraryPrep, SequenceDifferenceModel, SimpleAncientDnaModel,
            },
            AlignmentParameters,
        },
        utils::build_auxiliary_structures,
    };

    #[test]
    fn test_d() {
        let ref_seq = "GATTACA".as_bytes().to_owned(); // revcomp = TGTAATC

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let difference_model = SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.3,
                three_prime_overhang: 0.3,
            },
            0.001,
            0.8,
            0.02,
            false,
        );

        let representative_mismatch_penalty =
            difference_model.get_representative_mismatch_penalty();

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: TestBound {
                threshold: 0.0,
                representative_mm_bound: difference_model.get_representative_mismatch_penalty(),
            }
            .into(),
            penalty_gap_open: 0.00001_f32.log2(),
            penalty_gap_extend: representative_mismatch_penalty,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        let pattern = b"CCCCCCC";
        let base_qualities = [10, 40, 40, 40, 40, 10, 40];

        let center_of_read = pattern.len() / 2;

        let bi_d_array = BiDArray::new(
            pattern,
            &base_qualities,
            center_of_read,
            &parameters,
            &fmd_index,
            &difference_model,
        );

        assert_eq!(
            &*bi_d_array.d_composite,
            &[0.0, -3.6297126, -5.4444757, 0.0, -3.8959491, -3.8959491, -9.413074]
        );

        assert_eq!(
            bi_d_array.get(1, 4),
            bi_d_array.d_composite[1] + bi_d_array.d_composite[bi_d_array.split + 2]
        );
        assert_eq!(
            bi_d_array.get(2, 3),
            bi_d_array.d_composite[2] + bi_d_array.d_composite[bi_d_array.split + 3]
        );
        assert_eq!(
            bi_d_array.get(0, 6),
            bi_d_array.d_composite[0] + bi_d_array.d_composite[bi_d_array.split + 0]
        );

        assert_eq!(bi_d_array.get(0, pattern.len() as i16 - 1), 0.0,);
    }
}
