pub mod backtrack_tree;
pub mod fmd_index;
pub mod input_chunk_reader;
pub mod mapping;
pub mod mismatch_bounds;
pub mod prrange;
pub mod record;
pub mod sequence_difference_models;

mod bi_d_array;

use std::cmp::Ordering;

use serde::{Deserialize, Serialize};

use crate::map::{
    backtrack_tree::NodeId, fmd_index::RtBiInterval, mismatch_bounds::MismatchBoundDispatch,
    record::EditOperationsTrack, sequence_difference_models::SequenceDifferenceModelDispatch,
};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct AlignmentParameters {
    pub difference_model: SequenceDifferenceModelDispatch,
    pub mismatch_bound: MismatchBoundDispatch,
    pub penalty_gap_open: f32,
    pub penalty_gap_extend: f32,
    pub chunk_size: usize,
    pub gap_dist_ends: u8,
    pub max_num_gaps_open: u8,
    pub stack_limit_abort: bool,
}

/// A subset of `MismatchSearchStackFrame` to store hits
#[derive(Serialize, Deserialize, Debug)]
pub struct HitInterval {
    interval: RtBiInterval,
    alignment_score: f32,
    edit_operations: EditOperationsTrack,
}

impl PartialOrd for HitInterval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HitInterval {
    fn cmp(&self, other: &Self) -> Ordering {
        self.alignment_score
            .partial_cmp(&other.alignment_score)
            .expect("This is not expected to fail")
    }
}

impl PartialEq for HitInterval {
    fn eq(&self, other: &Self) -> bool {
        self.alignment_score.eq(&other.alignment_score)
    }
}

impl Eq for HitInterval {}

/// Simple zero-cost direction enum to increase readability
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Direction {
    Forward,
    Backward,
}

impl Direction {
    /// Reverses the direction from forward to backward and vice-versa
    fn reverse(self) -> Self {
        use self::Direction::{Backward, Forward};
        match self {
            Forward => Backward,
            Backward => Forward,
        }
    }

    #[allow(dead_code)]
    fn is_forward(self) -> bool {
        use self::Direction::{Backward, Forward};
        match self {
            Forward => true,
            Backward => false,
        }
    }

    #[allow(dead_code)]
    fn is_backward(self) -> bool {
        !self.is_forward()
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
enum GapState {
    Insertion,
    Deletion,
    Closed,
}

/// Stores information about partial alignments on the priority stack.
/// There are two different measures of alignment quality:
/// `alignment_score`: Initialized with 0, penalties are simply added
/// `priority`: `alignment_score` + expected minimal amount of penalties.
/// This is used as key for the priority stack.
#[derive(Debug)]
pub struct MismatchSearchStackFrame {
    current_interval: RtBiInterval,
    backward_index: i16,
    forward_index: i16,
    direction: Direction,
    gap_forwards: GapState,
    gap_backwards: GapState,
    num_gaps_open: u8,
    alignment_score: f32,
    edit_node_id: NodeId,
}

impl PartialOrd for MismatchSearchStackFrame {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MismatchSearchStackFrame {
    fn cmp(&self, other: &Self) -> Ordering {
        self.alignment_score
            .partial_cmp(&other.alignment_score)
            .expect("This is not expected to fail")
    }
}

impl PartialEq for MismatchSearchStackFrame {
    fn eq(&self, other: &Self) -> bool {
        self.alignment_score.eq(&other.alignment_score)
    }
}

impl Eq for MismatchSearchStackFrame {}

/// Holds data used for the output of alternative alignments in a BWA-like way
#[derive(Debug)]
struct AlternativeAlignments {
    // Number of best-scoring hits
    x0: i32,
    // Number of suboptimal hits within search space
    x1: i32,
    // Alternative hits; format: (chr,pos,CIGAR,MD,NM,n_hits,AS;){0,5}
    xa: String,
    // Alignment score of the second-best hit
    xs: f32,
    // Type: _U_nique/_R_epeat/_N_ot mapped
    xt: char,
}

/// Bundles the output of `interval2coordinate()`
#[derive(Debug)]
struct IntToCoordOutput<'a> {
    tid: u32,
    contig_name: &'a str,
    relative_pos: u64,
    absolute_pos: usize,
    strand: Direction,
    num_skipped: usize,
    interval: &'a HitInterval,
}

#[cfg(test)]
mod tests {
    use std::mem;

    use super::*;

    #[test]
    fn stack_frame_size() {
        assert_eq!(40, mem::size_of::<MismatchSearchStackFrame>());
    }
}
