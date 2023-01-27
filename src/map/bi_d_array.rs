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
        let mut offset_d_backward_iterators = (0..MAX_OFFSET)
            .map(|offset| {
                Self::compute_part(
                    &pattern[..split],
                    &base_qualities[..split],
                    Direction::Forward,
                    pattern.len(),
                    offset,
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
                .fold(0_f32, f32::min)
        });

        // Compute a forward-D-iterator for every offset
        let mut offset_d_forward_iterators = (0..MAX_OFFSET)
            .map(|offset| {
                Self::compute_part(
                    &pattern[split..],
                    &base_qualities[split..],
                    Direction::Backward,
                    pattern.len(),
                    offset,
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
                .fold(0_f32, f32::min)
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
        initial_skip: u16,
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
            .skip(usize::from(initial_skip))
            .scan(
                (
                    0.0,
                    i16::try_from(initial_skip).unwrap_or(0) - 1,
                    fmd_index.init_interval(),
                ),
                move |(z, last_mismatch_pos, interval), (index, &base)| {
                    *interval = match direction {
                        Direction::Forward => fmd_index.forward_ext(interval, base),
                        Direction::Backward => fmd_index.backward_ext(interval, base),
                    };
                    if interval.size < 1 {
                        // Sub-read does not align perfectly, scan sub-sequence to find the most conservative penalty
                        let pattern_qual_part_iter = pattern_part.iter().zip(base_qualities_part);
                        *z += match direction {
                            Direction::Forward => Either::Left(pattern_qual_part_iter),
                            Direction::Backward => Either::Right(pattern_qual_part_iter.rev()),
                        }
                        .enumerate()
                        .take(index + 1)
                        .skip((*last_mismatch_pos + 1) as usize)
                        .map(|(j, (&base_j, &qual_j))| {
                            let idx_mapped_to_read =
                                directed_index(j, full_pattern_length, direction);
                            let best_penalty_mm_only = sdm.get_min_penalty(
                                idx_mapped_to_read,
                                full_pattern_length,
                                base_j,
                                qual_j,
                                true,
                            );

                            let optimal_penalty = sdm.get_min_penalty(
                                idx_mapped_to_read,
                                full_pattern_length,
                                base_j,
                                qual_j,
                                false,
                            );
                            // The optimal penalty, conditioning on position and base, is subtracted because we
                            // optimize the ratio `AS/optimal_AS` (in log space) to find the best mappings
                            let mm_retval = best_penalty_mm_only - optimal_penalty;

                            // Don't count gap penalties where they are not allowed.
                            // Sadly this is not expected to have much of an effect.
                            if idx_mapped_to_read.min(full_pattern_length - idx_mapped_to_read - 1)
                                > alignment_parameters.gap_dist_ends as usize
                            {
                                mm_retval.max(alignment_parameters.penalty_gap_extend)
                                // The following is not needed, since we can assume
                                // penalty_gap_extend > penalty_gap_open + penalty_gap_extend
                                // .max(
                                //     alignment_parameters.penalty_gap_open
                                //         + alignment_parameters.penalty_gap_extend,
                                // )
                            } else {
                                mm_retval
                            }
                        })
                        .fold(f32::MIN, f32::max);
                        *interval = fmd_index.init_interval();
                        *last_mismatch_pos = index as i16;
                    }
                    Some(*z)
                },
            ),
        )
    }

    pub fn get(&self, backward_index: i16, forward_index: i16) -> f32 {
        debug_assert!(backward_index < forward_index);
        debug_assert!(backward_index < self.d_composite.len() as i16);
        debug_assert!(forward_index >= 0);

        debug_assert!(backward_index >= -1);
        let d_rev = usize::try_from(backward_index)
            .ok()
            .and_then(|composite_idx| self.d_composite.get(composite_idx).copied())
            .unwrap_or(0.0);

        debug_assert!((forward_index as usize) <= self.d_composite.len());
        let d_fwd = self
            .d_composite
            .len()
            // It is important that the following subtraction underflows if
            // `forward_index >= d_composite.len()`.
            // Reordering the subtraction and addition of `self.split` will break index calculation.
            .checked_sub(1 + forward_index as usize)
            .map(|left_half_idx| left_half_idx + self.split)
            .and_then(|composite_idx| self.d_composite.get(composite_idx).copied())
            .unwrap_or(0.0);

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
            max_num_gaps_open: 2,
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
            bi_d_array.d_composite[0] + bi_d_array.d_composite[bi_d_array.split]
        );

        assert_eq!(bi_d_array.get(0, pattern.len() as i16 - 1), 0.0,);
    }
}
