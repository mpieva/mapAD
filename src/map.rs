use std::fs::OpenOptions;
use std::{
    cell::RefCell,
    cmp::Ordering,
    collections::{binary_heap::BinaryHeap, BTreeMap},
    fs::File,
    io::{self, BufReader, Read, Write},
    iter::Map,
    path::Path,
    time::{Duration, Instant},
};

use bio::{alphabets::dna, data_structures::suffix_array::SuffixArray, io::fastq};
use clap::{crate_description, crate_version};
use either::Either;
use log::{debug, error, info, trace, warn};
use min_max_heap::MinMaxHeap;
use noodles::{bam, sam};
use rand::RngCore;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;

use crate::{
    backtrack_tree::{NodeId, Tree},
    errors::{Error, Result},
    fmd_index::{RtBiInterval, RtFmdIndex},
    mismatch_bounds::{MismatchBound, MismatchBoundDispatch},
    prrange::PrRange,
    sequence_difference_models::{SequenceDifferenceModel, SequenceDifferenceModelDispatch},
    utils::{
        load_id_pos_map_from_path, load_index_from_path, load_suffix_array_from_path,
        AlignmentParameters, Record,
    },
    CRATE_NAME,
};

// These settings lead to a memory consumption of ~245 MiB per thread
pub const STACK_LIMIT: u32 = 2_000_000;
pub const EDIT_TREE_LIMIT: u32 = 10_000_000;

/// A subset of MismatchSearchStackFrame to store hits
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
#[derive(Debug, Copy, Clone, PartialEq)]
enum Direction {
    Forward,
    Backward,
}

impl Direction {
    /// Reverses the direction from forward to backward and vice-versa
    fn reverse(self) -> Self {
        use self::Direction::*;
        match self {
            Forward => Backward,
            Backward => Forward,
        }
    }

    #[allow(dead_code)]
    fn is_forward(self) -> bool {
        use self::Direction::*;
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

/// Variants store position in the read and, if necessary, the reference base
#[derive(Debug, Copy, Clone, PartialEq, Serialize, Deserialize)]
pub enum EditOperation {
    Insertion(u16),
    Deletion(u16, u8),
    Match(u16),
    Mismatch(u16, u8),
}

impl Default for EditOperation {
    fn default() -> Self {
        Self::Match(0)
    }
}

impl From<&EditOperation> for sam::record::cigar::op::Kind {
    fn from(src: &EditOperation) -> Self {
        use sam::record::cigar::op::Kind;
        match src {
            EditOperation::Insertion(_) => Kind::Insertion,
            EditOperation::Deletion(_, _) => Kind::Deletion,
            EditOperation::Match(_) => Kind::Match,
            EditOperation::Mismatch(_, _) => Kind::Match,
        }
    }
}

impl From<&EditOperation> for sam::record::cigar::Op {
    fn from(src: &EditOperation) -> Self {
        use sam::record::cigar::op::{Kind, Op};
        match src {
            EditOperation::Insertion(l) => Op::new(Kind::Insertion, (*l).into()),
            EditOperation::Deletion(l, _) => Op::new(Kind::Deletion, (*l).into()),
            EditOperation::Match(l) => Op::new(Kind::Match, (*l).into()),
            EditOperation::Mismatch(l, _) => Op::new(Kind::Match, (*l).into()),
        }
    }
}

/// Contains edit operations performed in order to align the sequence
#[derive(Debug, Serialize, Deserialize)]
pub struct EditOperationsTrack(Vec<EditOperation>);

impl EditOperationsTrack {
    /// Calculates the amount of positions in the genome
    /// that are covered by this read
    pub fn effective_len(&self) -> usize {
        self.0.iter().fold(0, |acc, edit_operation| {
            acc + match edit_operation {
                EditOperation::Insertion(_) => 0,
                EditOperation::Deletion(_, _) => 1,
                EditOperation::Match(_) => 1,
                EditOperation::Mismatch(_, _) => 1,
            }
        })
    }

    /// Constructs CIGAR, MD tag, and edit distance from correctly ordered track of edit operations and yields them as a tuple
    /// The strand a read is mapped to is taken into account here.
    fn to_bam_fields(&self, strand: Direction) -> (Vec<bam::record::cigar::Op>, Vec<u8>, u16) {
        use bam::record::cigar::Op;
        // Reconstruct the order of the remaining edit operations and condense CIGAR
        let mut num_matches: u32 = 0;
        let mut num_operations = 1;
        let mut edit_distance = 0;
        let mut last_edit_operation = None;
        let mut cigar = Vec::new();
        let mut md_tag = Vec::new();

        let track = match strand {
            Direction::Forward => Either::Left(self.0.iter()),
            Direction::Backward => Either::Right(self.0.iter().rev()),
        };
        for edit_operation in track {
            edit_distance = Self::add_edit_distance(edit_operation, edit_distance);

            num_matches = Self::add_md_edit_operation(
                Some(edit_operation),
                last_edit_operation,
                strand,
                num_matches,
                &mut md_tag,
            );

            if let Some(lop) = last_edit_operation {
                // Construct CIGAR string
                match edit_operation {
                    EditOperation::Match(_) => match lop {
                        EditOperation::Match(_) | EditOperation::Mismatch(_, _) => {
                            num_operations += 1;
                        }
                        _ => {
                            cigar.push(
                                Op::new(lop.into(), num_operations)
                                    .expect("This is not expected to happen"),
                            );
                            num_operations = 1;
                            last_edit_operation = Some(edit_operation);
                        }
                    },
                    EditOperation::Mismatch(_, _) => match lop {
                        EditOperation::Mismatch(_, _) | EditOperation::Match(_) => {
                            num_operations += 1;
                        }
                        _ => {
                            cigar.push(
                                Op::new(lop.into(), num_operations)
                                    .expect("This is not expected to happen"),
                            );
                            num_operations = 1;
                            last_edit_operation = Some(edit_operation);
                        }
                    },
                    EditOperation::Insertion(_) => match lop {
                        EditOperation::Insertion(_) => {
                            num_operations += 1;
                        }
                        _ => {
                            cigar.push(
                                Op::new(lop.into(), num_operations)
                                    .expect("This is not expected to happen"),
                            );
                            num_operations = 1;
                            last_edit_operation = Some(edit_operation);
                        }
                    },
                    EditOperation::Deletion(_, _) => match lop {
                        EditOperation::Deletion(_, _) => {
                            num_operations += 1;
                        }
                        _ => {
                            cigar.push(
                                Op::new(lop.into(), num_operations)
                                    .expect("This is not expected to happen"),
                            );
                            num_operations = 1;
                            last_edit_operation = Some(edit_operation);
                        }
                    },
                }
            } else {
                last_edit_operation = Some(edit_operation);
            }
        }

        // Add remainder
        if let Some(lop) = last_edit_operation {
            cigar
                .push(Op::new(lop.into(), num_operations).expect("This is not expected to happen"));
        }
        let _ = Self::add_md_edit_operation(None, None, strand, num_matches, &mut md_tag);

        (cigar, md_tag, edit_distance)
    }

    fn add_md_edit_operation(
        edit_operation: Option<&EditOperation>,
        last_edit_operation: Option<&EditOperation>,
        strand: Direction,
        mut k: u32,
        md_tag: &mut Vec<u8>,
    ) -> u32 {
        let comp_if_necessary = |reference_base| match strand {
            Direction::Forward => reference_base,
            Direction::Backward => dna::complement(reference_base),
        };

        match edit_operation {
            Some(EditOperation::Match(_)) => k += 1,
            Some(EditOperation::Mismatch(_, reference_base)) => {
                let reference_base = comp_if_necessary(*reference_base);
                md_tag.extend_from_slice(format!("{}{}", k, reference_base as char).as_bytes());
                k = 0;
            }
            Some(EditOperation::Insertion(_)) => {
                // Insertions are ignored in MD tags
            }
            Some(EditOperation::Deletion(_, reference_base)) => {
                let reference_base = comp_if_necessary(*reference_base);
                match last_edit_operation {
                    Some(EditOperation::Deletion(_, _)) => {
                        md_tag.extend_from_slice(format!("{}", reference_base as char).as_bytes());
                    }
                    _ => {
                        md_tag.extend_from_slice(
                            format!("{}^{}", k, reference_base as char).as_bytes(),
                        );
                    }
                }
                k = 0;
            }
            None => md_tag.extend_from_slice(format!("{}", k).as_bytes()),
        }
        k
    }

    fn add_edit_distance(edit_operation: &EditOperation, distance: u16) -> u16 {
        if let EditOperation::Match(_) = edit_operation {
            distance
        } else {
            distance + 1
        }
    }

    pub fn read_len(&self) -> usize {
        self.0.iter().fold(0, |acc, edit_operation| {
            acc + match edit_operation {
                EditOperation::Insertion(_) => 1,
                EditOperation::Deletion(_, _) => 0,
                EditOperation::Match(_) => 1,
                EditOperation::Mismatch(_, _) => 1,
            }
        })
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
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
/// alignment_score: Initialized with 0, penalties are simply added
/// priority: alignment_score + expected minimal amount of penalties.
/// This is used as key for the priority stack.
#[derive(Debug)]
pub struct MismatchSearchStackFrame {
    current_interval: RtBiInterval,
    backward_index: i16,
    forward_index: i16,
    direction: Direction,
    gap_forwards: GapState,
    gap_backwards: GapState,
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
    fn get_reference_identifier(
        &self,
        position: usize,
        pattern_length: usize,
    ) -> Option<(u32, u64)> {
        let position = position as u64;
        self.id_position
            .iter()
            .enumerate()
            .find(|(_, identifier)| {
                (identifier.start <= position)
                    && (position + pattern_length as u64 - 1 <= identifier.end)
            })
            .and_then(|(index, identifier)| {
                Some((u32::try_from(index).ok()?, position - identifier.start))
            })
    }
}

/// Compute the lower bound of mismatches of a read per position by aligning perfectly,
/// starting from the read end and progressing to the center. As soon as the extension of
/// the alignment is not possible, we found at least one mismatch and record that per
/// read-position in the so-called D-array. The content of the array is used as "lookahead"
/// score. This allows for pruning the search tree. The values are minimal expected penalties
/// towards the respective ends of the query. Like alignment scores, these values are negative.
#[derive(Debug)]
struct BiDArray {
    d_composite: Vec<f32>,
    split: usize,
}

impl BiDArray {
    pub(crate) fn new<SDM>(
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

        std::iter::repeat(0.0).take(initial_skip as usize).chain(
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
                        .max(alignment_parameters.penalty_gap_open)
                        .max(alignment_parameters.penalty_gap_extend);
                        *interval = fmd_index.init_interval();
                        *last_mismatch_pos = index as i16;
                    }
                    Some(*z)
                },
            ),
        )
    }

    fn get(&self, backward_index: i16, forward_index: i16) -> f32 {
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

/// Reads chunks of configurable size from source Iterator
///
/// Very basic error checking, reporting, and recovery happens here.
struct ChunkIterator<T> {
    chunk_size: usize,
    records: T,
}

impl<T> Iterator for ChunkIterator<T>
where
    T: Iterator<Item = Result<Record>>,
{
    type Item = Vec<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let chunk = self
            .records
            .by_ref()
            .filter_map(|record| {
                if let Err(ref e) = record {
                    error!("Skip record due to an error: {}", e);
                }
                record.ok()
            })
            .filter(|record| {
                let length_check = record.sequence.len() == record.base_qualities.len();
                if !length_check {
                    error!(
                        "Skip record \"{}\" due to different length of sequence and quality strings",
                        String::from_utf8_lossy(&record.name)
                    );
                }
                length_check
            })
            .take(self.chunk_size)
            .collect::<Vec<_>>();

        // If the underlying iterator is exhausted return None, too
        if chunk.is_empty() {
            return None;
        }

        Some(chunk)
    }
}

/// Convertible to ChunkIterator
trait IntoChunkIterator<E, I, O, T>
where
    E: Into<Error>,
    I: Into<Record>,
    T: Iterator<Item = std::result::Result<I, E>>,
    O: Iterator<Item = Result<Record>>,
{
    fn into_chunks(self, chunk_size: usize) -> ChunkIterator<O>;
}

/// Adds ChunkIterator conversion method to every compatible Iterator. So when new input file types
/// are implemented, it is sufficient to impl `From<T> for Record` for the additional item `T`.
#[allow(clippy::type_complexity)]
impl<E, I, T> IntoChunkIterator<E, I, Map<T, fn(std::result::Result<I, E>) -> Result<Record>>, T>
    for T
where
    E: Into<Error>,
    I: Into<Record>,
    T: Iterator<Item = std::result::Result<I, E>>,
{
    fn into_chunks(
        self,
        chunk_size: usize,
    ) -> ChunkIterator<Map<T, fn(std::result::Result<I, E>) -> Result<Record>>> {
        ChunkIterator {
            chunk_size,
            records: self.map(|inner| inner.map(|v| v.into()).map_err(|e| e.into())),
        }
    }
}

/// Loads index files and launches the mapping process
pub fn run(
    reads_path: &str,
    reference_path: &str,
    out_file_path: &str,
    alignment_parameters: &AlignmentParameters,
) -> Result<()> {
    let reads_path = Path::new(reads_path);
    let out_file_path = Path::new(out_file_path);
    if !reads_path.exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "The given input file could not be found",
        )
        .into());
    }

    info!("Load FMD-index");
    let fmd_index = load_index_from_path(reference_path)?;

    info!("Load suffix array");
    let sampled_suffix_array_owned = load_suffix_array_from_path(reference_path)?;
    let suffix_array = sampled_suffix_array_owned.into_sampled_suffix_array(
        &fmd_index.bwt,
        &fmd_index.less,
        &fmd_index.occ,
    );

    info!("Load position map");
    let identifier_position_map = load_id_pos_map_from_path(reference_path)?;

    let mut out_file = bam::Writer::new(
        OpenOptions::new()
            .read(false)
            .write(true)
            .create_new(true)
            .open(out_file_path)?,
    );

    info!("Map reads");
    // Static dispatch of the Record type based on the filename extension
    match reads_path
        .extension()
        .ok_or(Error::InvalidInputType)?
        .to_str()
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "The provided file name contains invalid unicode",
            )
        })? {
        "bam" => {
            let mut reader = bam::Reader::new(File::open(reads_path)?);
            let header = create_bam_header(Some(&mut reader), &identifier_position_map)?;

            // Position cursor right after header
            let _ = reader.read_reference_sequences();

            out_file.write_header(&header)?;
            out_file.write_reference_sequences(header.reference_sequences())?;
            run_inner(
                reader
                    .records()
                    .into_chunks(alignment_parameters.chunk_size),
                &fmd_index,
                &suffix_array,
                alignment_parameters,
                &identifier_position_map,
                &mut out_file,
            )?
        }
        "fastq" | "fq" => {
            let reader = fastq::Reader::from_file(reads_path)?;
            let header = create_bam_header::<BufReader<File>>(None, &identifier_position_map)?;
            out_file.write_header(&header)?;
            out_file.write_reference_sequences(header.reference_sequences())?;
            run_inner(
                reader
                    .records()
                    .into_chunks(alignment_parameters.chunk_size),
                &fmd_index,
                &suffix_array,
                alignment_parameters,
                &identifier_position_map,
                &mut out_file,
            )?
        }
        _ => return Err(Error::InvalidInputType),
    }

    info!("Done");
    Ok(())
}

/// This part has been extracted from the main run() function to allow static dispatch based on the
/// input file type
fn run_inner<S, T, W>(
    records: ChunkIterator<T>,
    fmd_index: &RtFmdIndex,
    suffix_array: &S,
    alignment_parameters: &AlignmentParameters,
    identifier_position_map: &FastaIdPositions,
    out_file: &mut bam::Writer<W>,
) -> Result<()>
where
    S: SuffixArray + Send + Sync,
    T: Iterator<Item = Result<Record>>,
    W: Write,
{
    thread_local! {
        static STACK_BUF: RefCell<MinMaxHeap<MismatchSearchStackFrame>> = RefCell::new(MinMaxHeap::with_capacity(STACK_LIMIT as usize + 9));
        static TREE_BUF: RefCell<Tree<EditOperation>> = RefCell::new(Tree::with_capacity(EDIT_TREE_LIMIT + 9));
    }

    for chunk in records {
        debug!("Map chunk of reads");
        let results = chunk
            .into_par_iter()
            .map(|record| {
                STACK_BUF.with(|stack_buf| {
                    TREE_BUF.with(|tree_buf| {
                        let start = Instant::now();
                        // Here we call the instances of the generic `k_mismatch_search` function. Currently, only
                        // `SimpleAncientDnaModel` is used, the others merely serve as an example.
                        let hit_intervals = match &alignment_parameters.difference_model {
                            SequenceDifferenceModelDispatch::SimpleAncientDnaModel(sdm) => {
                                match &alignment_parameters.mismatch_bound {
                                    MismatchBoundDispatch::Discrete(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                    MismatchBoundDispatch::Continuous(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                    MismatchBoundDispatch::TestBound(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                }
                            }
                            SequenceDifferenceModelDispatch::TestDifferenceModel(sdm) => {
                                match &alignment_parameters.mismatch_bound {
                                    MismatchBoundDispatch::Discrete(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                    MismatchBoundDispatch::Continuous(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                    MismatchBoundDispatch::TestBound(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                }
                            }
                            SequenceDifferenceModelDispatch::VindijaPwm(sdm) => {
                                match &alignment_parameters.mismatch_bound {
                                    MismatchBoundDispatch::Discrete(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                    MismatchBoundDispatch::Continuous(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                    MismatchBoundDispatch::TestBound(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                }
                            }
                        };
                        let duration = start.elapsed();

                        (record, hit_intervals, duration)
                    })
                })
            })
            .map_init(
                rand::thread_rng,
                |mut rng, (record, hit_interval, duration)| -> Result<bam::Record> {
                    intervals_to_bam(
                        record,
                        hit_interval,
                        suffix_array,
                        identifier_position_map,
                        Some(&duration),
                        alignment_parameters,
                        &mut rng,
                    )
                },
            )
            .collect::<Result<Vec<_>>>()?;

        debug!("Write BAM records to output file serially");
        for record in results.iter() {
            out_file.write_record(record)?;
        }
    }
    Ok(())
}

/// Creates basic BAM header with mandatory fields pre-populated.
/// If provided, it will copy most data from an input BAM header (`@SQ` lines are removed).
pub fn create_bam_header<R>(
    reads_reader: Option<&mut bam::Reader<R>>,
    identifier_position_map: &FastaIdPositions,
) -> Result<sam::Header>
where
    R: Read,
{
    let mut header_builder = sam::Header::builder();
    let mut header_header_builder = sam::header::header::Header::builder();
    header_header_builder =
        header_header_builder.set_version(sam::header::header::Version::new(1, 6));
    header_header_builder =
        header_header_builder.set_sort_order(sam::header::header::SortOrder::Unsorted);

    let pg_id = CRATE_NAME;

    let mut program_builder = {
        let cmdline = {
            let mut out = std::env::args().fold(String::new(), |acc, part| acc + &part + " ");
            let _ = out.pop();
            out
        };
        sam::header::Program::builder()
            .set_id(pg_id)
            .set_name(CRATE_NAME)
            .set_version(crate_version!())
            .set_description(crate_description!())
            .set_command_line(cmdline)
    };

    if let Some(src_header) = reads_reader
        .and_then(|bam_reader|
            bam_reader
                .read_header()
                .ok()
                .and_then(|header| header.parse::<sam::Header>().ok())
                .or_else(|| {
                    warn!("Could not read input file header. Instead, create output header from scratch. Some metadata might be lost.");
                    None
                })
        ) {
        // We've got an header to work with!
        // Retrieve custom header (@HD) fields (version is not included):
        let header = src_header.header().map(|header| header.fields());
        if let Some(fields) = header {
            for (custom_tag, value) in fields.iter() {
                header_header_builder = header_header_builder.insert(custom_tag.to_owned(), value);
            }
        }

        // @PG chain of old entries
        for (_id, pg) in src_header.programs().iter() {
            header_builder = header_builder.add_program(pg.clone());
        }

        // Append our program line to the latest end of a chain
        for pg_id in src_header.programs().keys().rev() {
            // Ensure it's really the end of a chain
            if !src_header.programs().values().filter_map(|pg| pg.previous_id()).any(|x| x == pg_id.as_str()) {
                program_builder = program_builder.set_previous_id(pg_id);
                break;
            }
        }

        // Ensure @PG/ID is unique
        let pg_id_count = src_header.programs().keys().filter(|id| id.as_str() == pg_id || id.starts_with(&format!("{}.", pg_id))).count();
        if pg_id_count > 0 {
            program_builder = program_builder.set_id(format!("{}.{}", pg_id, pg_id_count))
        }

        for comment in src_header.comments().iter() {
            header_builder = header_builder.add_comment(comment);
        }

        for (_id, read_group) in src_header.read_groups().iter() {
            header_builder = header_builder.add_read_group(read_group.clone());
        }
    }

    // New @PG entry
    let program = program_builder
        .build()
        .expect("This is not expected to fail");
    header_builder = header_builder.add_program(program);

    // @SQ entries
    for identifier_position in identifier_position_map.iter() {
        header_builder = header_builder.add_reference_sequence(
            sam::header::ReferenceSequence::new(
                identifier_position.identifier.parse().map_err(|_e| {
                    Error::InvalidIndex(format!(
                        "Could not create header. Contig name \"{}\" can not be used as @SQ ID.",
                        identifier_position.identifier
                    ))
                })?,
                (identifier_position.end - identifier_position.start + 1) as i32,
            )
            .map_err(|_e| {
                Error::InvalidIndex("Could not create header. Contig length invalid?".into())
            })?,
        );
    }

    // @HD entries
    header_builder = header_builder.set_header(header_header_builder.build());

    Ok(header_builder.build())
}

/// Convert suffix array intervals to positions and BAM records
pub fn intervals_to_bam<R, S>(
    input_record: Record,
    mut intervals: BinaryHeap<HitInterval>,
    suffix_array: &S,
    identifier_position_map: &FastaIdPositions,
    duration: Option<&Duration>,
    alignment_parameters: &AlignmentParameters,
    rng: &mut R,
) -> Result<bam::Record>
where
    R: RngCore,
    S: SuffixArray,
{
    let mut intervals_to_coordinates =
        |interval: &HitInterval| -> Result<(usize, u32, u64, Direction)> {
            // Calculate length of reference strand ignoring sentinel characters
            let strand_len = (suffix_array.len() - 2) / 2;
            // Get amount of positions in the genome covered by the read
            let effective_read_len = interval.edit_operations.effective_len();

            // If this function returns an error that's likely because the best-scoring interval mappings
            // overlap contig boundaries. In that case, we can drop the best-scoring interval and
            // re-evaluate (coordinates & MQ) the mapping with the remainder of the hit interval heap.
            PrRange::try_from_range(&interval.interval.range_fwd(), rng.next_u32() as usize)
                .ok_or_else(|| {
                    Error::InvalidIndex("Could not enumerate possible reference positions".into())
                })?
                // This is to count the number of skipped invalid positions
                .enumerate()
                .find_map(|(i, sar_pos)| {
                    suffix_array
                        .get(sar_pos)
                        // Determine strand
                        .map(|absolute_pos| {
                            if absolute_pos < strand_len {
                                (absolute_pos, Direction::Forward)
                            } else {
                                (
                                    suffix_array.len() - absolute_pos - effective_read_len - 1,
                                    Direction::Backward,
                                )
                            }
                        })
                        // Convert to relative coordinates
                        .and_then(|(absolute_pos, strand)| {
                            if let Some((tid, rel_pos)) = identifier_position_map
                                .get_reference_identifier(absolute_pos, effective_read_len)
                            {
                                Some((i, tid, rel_pos, strand))
                            } else {
                                None
                            }
                        })
                })
                .ok_or(Error::ContigBoundaryOverlap)
        };

    let hits_found = !intervals.is_empty();

    while let Some(best_alignment) = intervals.pop() {
        // Determine relative-to-chromosome position
        match intervals_to_coordinates(&best_alignment) {
            Ok((num_failed, tid, relative_pos, strand)) => {
                let updated_best_alignment_interval_size =
                    best_alignment.interval.size - num_failed;

                let alternative_hits = if updated_best_alignment_interval_size > 1 {
                    format!("mm,{},;", updated_best_alignment_interval_size)
                } else if intervals.len() > 9 {
                    "pu,,;".into()
                } else if !intervals.is_empty() {
                    let mut buf = "pu,".to_string();
                    buf.extend(intervals.iter().map(|subopt_intv| {
                        format!(
                            "{},{:.2};",
                            subopt_intv.interval.size, subopt_intv.alignment_score
                        )
                    }));
                    buf
                } else if updated_best_alignment_interval_size == 1 {
                    "uq,,;".into()
                } else {
                    "".into()
                };
                return bam_record_helper(
                    input_record,
                    Some(relative_pos),
                    Some(&best_alignment),
                    Some(estimate_mapping_quality(
                        &best_alignment,
                        updated_best_alignment_interval_size,
                        &intervals,
                        alignment_parameters,
                    )),
                    Some(tid),
                    Some(strand),
                    duration,
                    alternative_hits,
                );
            }
            Err(e) => {
                debug!("{}. Try again with next best hit.", e);
            }
        }
    }

    if hits_found {
        debug!(
            "Hits could not be mapped to valid coordinates. Report read \"{}\" as unmapped.",
            std::str::from_utf8(&input_record.name).expect("Read names are validated elsewhere")
        );
    }

    // No match found, report unmapped read
    bam_record_helper(
        input_record,
        None,
        None,
        Some(0),
        None,
        None,
        duration,
        "".into(),
    )
}

/// Computes optimal per-base alignment scores for a read,
/// conditioned on base qualities and scoring model.
/// This function panics if `pattern` and `base_qualities` are not of equal length.
pub fn compute_optimal_scores<SDM>(
    pattern: &[u8],
    base_qualities: &[u8],
    difference_model: &SDM,
) -> Vec<f32>
where
    SDM: SequenceDifferenceModel,
{
    pattern
        .iter()
        .zip(base_qualities)
        .enumerate()
        .map(|(i, (&base, &quality))| {
            difference_model.get_min_penalty(i, pattern.len(), base, quality, false)
        })
        .collect::<Vec<_>>()
}

/// Estimate mapping quality based on the number of hits for a particular read, its alignment score,
/// and its base qualities
fn estimate_mapping_quality(
    best_alignment: &HitInterval,
    best_alignment_interval_size: usize,
    other_alignments: &BinaryHeap<HitInterval>,
    alignment_parameters: &AlignmentParameters,
) -> u8 {
    const MAX_MAPQ: u8 = 37;
    const MIN_MAPQ_UNIQ: u8 = 25;
    assert!(MIN_MAPQ_UNIQ <= MAX_MAPQ);

    let alignment_probability = {
        let ratio_best = 2_f32.powf(best_alignment.alignment_score);
        if best_alignment_interval_size > 1 {
            // Multi-mapping
            1.0 / best_alignment_interval_size as f32
        } else if other_alignments.is_empty() {
            // Unique mapping
            1.0
        } else {
            // Pseudo-unique mapping
            let weighted_suboptimal_alignments =
                other_alignments
                    .iter()
                    .fold(0.0, |acc, suboptimal_alignment| {
                        acc + 2_f32.powf(suboptimal_alignment.alignment_score)
                            * suboptimal_alignment.interval.size as f32
                    });
            ratio_best / (ratio_best + weighted_suboptimal_alignments)
        }
    }
    // Guard against rounding errors
    .clamp(0.0, 1.0);

    // Produce Phred score
    let mapping_quality = (-10.0 * (1.0 - alignment_probability).log10())
        .min(MAX_MAPQ as f32)
        .round() as u8;

    // When we found a best-scoring hit, we search up until `AS + mm_penalty` to find suboptimal
    // hits in order to estimate the mapping quality. However, we don't ever search beyond the
    // `mismatch_bound` limit. So when we are closer than `mm_penalty` to the `mismatch_bound`,
    // the mapping quality needs to be scaled down to reflect the fact that we are less likely to
    // find suboptimal alignments within these narrowed bounds.
    if mapping_quality == MAX_MAPQ {
        let scaled_mq = MIN_MAPQ_UNIQ as f32
            + ((MAX_MAPQ - MIN_MAPQ_UNIQ) as f32
                * alignment_parameters
                    .mismatch_bound
                    .remaining_frac_of_repr_mm(
                        best_alignment.alignment_score,
                        best_alignment.edit_operations.read_len(),
                    )
                    .min(1.0));
        return scaled_mq.round() as u8;
    }

    mapping_quality
}

/// Create and return a BAM record of either a hit or an unmapped read
fn bam_record_helper(
    input_record: Record,
    position: Option<u64>,
    hit_interval: Option<&HitInterval>,
    mapq: Option<u8>,
    tid: Option<u32>,
    strand: Option<Direction>,
    duration: Option<&Duration>,
    // Contains valid content for the `YA` tag
    alternative_hits: String,
) -> Result<bam::Record> {
    let mut bam_builder = bam::Record::builder();

    let (cigar, md_tag, edit_distance) = if let Some(hit_interval) = hit_interval {
        let (cigar, md_tag, edit_distance) = hit_interval
            .edit_operations
            .to_bam_fields(strand.expect("This is not expected to fail"));
        (Some(cigar), Some(md_tag), Some(edit_distance))
    } else {
        (None, None, None)
    };

    // Copy flags from input record
    let mut flags = sam::record::Flags::from(input_record.bam_flags);

    if let Some(position) = position {
        flags.remove(sam::record::Flags::UNMAPPED);
        bam_builder = bam_builder.set_position(
            i32::try_from(position + 1)
                .ok()
                .and_then(|pos32| pos32.try_into().ok())
                .ok_or_else(|| Error::InvalidIndex("Could not compute valid coordinate".into()))?,
        )
    } else {
        flags.insert(sam::record::Flags::UNMAPPED);
        flags.remove(sam::record::Flags::REVERSE_COMPLEMENTED);
        flags.remove(sam::record::Flags::PROPERLY_ALIGNED);
    }

    // Flag read that maps to reverse strand
    if let Some(Direction::Backward) = strand {
        flags.insert(sam::record::Flags::REVERSE_COMPLEMENTED);
    } else {
        flags.remove(sam::record::Flags::REVERSE_COMPLEMENTED);
    }

    // Set mandatory properties of the BAM record
    bam_builder = bam_builder
        .set_read_name(input_record.name.to_owned())
        .set_flags(flags);

    // Add optional fields (or do not)

    if let Some(tid) = tid {
        let tid = i32::try_from(tid)
            .ok()
            .and_then(|tid| tid.try_into().ok())
            .ok_or_else(|| Error::InvalidIndex("Invalid reference map".into()))?;
        bam_builder = bam_builder.set_reference_sequence_id(tid);
    }

    // Some (not all) fields need to be reversed when the read maps to the reverse strand
    match strand {
        Some(Direction::Forward) | None => {
            bam_builder = bam_builder
                .set_sequence(
                    bam::record::Sequence::try_from(input_record.sequence.as_slice()).map_err(
                        |e| Error::Hts(format!("Could not create valid sequence: {}", e)),
                    )?,
                )
                .set_quality_scores(input_record.base_qualities.into());
        }
        Some(Direction::Backward) => {
            // CIGAR strings and MD tags are reversed during generation
            bam_builder = bam_builder
                .set_sequence(
                    bam::record::Sequence::try_from(
                        dna::revcomp(&input_record.sequence).as_slice(),
                    )
                    .map_err(|e| Error::Hts(format!("Could not create valid sequence: {}", e)))?,
                )
                .set_quality_scores(
                    input_record
                        .base_qualities
                        .iter()
                        .rev()
                        .copied()
                        .collect::<Vec<_>>()
                        .into(),
                );
        }
    }

    // CIGAR strings and MD tags are reversed during generation
    if let Some(cigar) = cigar {
        bam_builder = bam_builder.set_cigar(cigar.into());
    }

    if let Some(mapq) = mapq {
        bam_builder = bam_builder
            .set_mapping_quality(mapq.try_into().expect("This is not expected to happen"));
    }

    // Add optional tags

    // Add tags that were already present in the input
    let mut aux_data = input_record
        .bam_tags
        .into_iter()
        .filter(|(tag, _v)| tag != b"AS")
        .filter(|(tag, _v)| tag != b"NM")
        .filter(|(tag, _v)| tag != b"MD")
        .filter(|(tag, _v)| tag != b"YA")
        .filter(|(tag, _v)| tag != b"XD")
        .map(|(tag, value)| {
            Ok(bam::record::data::Field::new(
                tag.as_slice()
                    .try_into()
                    .map_err(|_e| Error::ParseError("Could not read input data tag".into()))?,
                value.into(),
            ))
        })
        .collect::<Result<Vec<_>>>()?;

    if let Some(hit_interval) = hit_interval {
        aux_data.push(bam::record::data::Field::new(
            sam::record::data::field::Tag::AlignmentScore,
            bam::record::data::field::Value::Float(hit_interval.alignment_score),
        ));
    };

    if let Some(edit_distance) = edit_distance {
        aux_data.push(bam::record::data::Field::new(
            sam::record::data::field::Tag::EditDistance,
            bam::record::data::field::Value::Int32(edit_distance as i32),
        ));
    };

    // CIGAR strings and MD tags are reversed during generation
    if let Some(md_tag) = md_tag {
        aux_data.push(bam::record::data::Field::new(
            sam::record::data::field::Tag::MismatchedPositions,
            bam::record::data::field::Value::String(
                String::from_utf8(md_tag).expect("This is not expected to fail"),
            ),
        ));
    }

    // Our format differs in that we include alignment scores, so we use `YA` instead of `XA`
    if !alternative_hits.is_empty() {
        aux_data.push(bam::record::data::Field::new(
            b"YA"
                .as_slice()
                .try_into()
                .expect("This is not expected to fail"),
            bam::record::data::field::Value::String(alternative_hits),
        ));
    }

    if let Some(duration) = duration {
        // Add the time that was needed for mapping the read
        aux_data.push(bam::record::data::Field::new(
            b"XD"
                .as_slice()
                .try_into()
                .expect("This is not expected to fail"),
            bam::record::data::field::Value::Float(duration.as_secs_f32()),
        ));
    }

    bam_builder = bam_builder.set_data(
        aux_data
            .try_into()
            .map_err(|e| Error::Hts(format!("Could not create valid auxiliary data: {}", e)))?,
    );

    bam_builder
        .build()
        .map_err(|e| Error::Hts(format!("Could not create valid record: {}", e)))
}

/// Derive Cigar string from oddly-ordered tracks of edit operations.
/// Since we start aligning at the center of a read, tracks of edit operations
/// are not ordered by position in the read. Also, the track of edit operations
/// must be extracted from the (possibly huge) tree which is built during
/// backtracking for size reasons.
fn extract_edit_operations(
    end_node: NodeId,
    edit_tree: &Tree<EditOperation>,
    pattern_len: usize,
) -> EditOperationsTrack {
    // Restore outer ordering of the edit operation by the positions they carry as values.
    // Whenever there are deletions in the pattern, there is no simple rule to reconstruct the ordering.
    // So, edit operations carrying the same position are pushed onto the same bucket and dealt with later.
    let mut cigar_order_outer: BTreeMap<u16, SmallVec<[EditOperation; 8]>> = BTreeMap::new();

    edit_tree.ancestors(end_node).for_each(|&edit_operation| {
        let position = match edit_operation {
            EditOperation::Insertion(position) => position,
            EditOperation::Deletion(position, _) => position,
            EditOperation::Match(position) => position,
            EditOperation::Mismatch(position, _) => position,
        };
        cigar_order_outer
            .entry(position)
            .or_insert_with(SmallVec::new)
            .push(edit_operation);
    });

    EditOperationsTrack(
        cigar_order_outer
            .into_iter()
            .flat_map(|(i, inner_vec)| {
                if i < (pattern_len / 2) as u16 {
                    Either::Left(inner_vec.into_iter())
                } else {
                    Either::Right(inner_vec.into_iter().rev())
                }
            })
            .collect(),
    )
}

/// Checks stop-criteria of stack frames before pushing them onto the stack.
/// Since push operations on heaps are costly, this should accelerate the alignment.
fn check_and_push_stack_frame<MB>(
    mut stack_frame: MismatchSearchStackFrame,
    pattern: &[u8],
    edit_operation: EditOperation,
    edit_tree: &mut Tree<EditOperation>,
    stack: &mut MinMaxHeap<MismatchSearchStackFrame>,
    intervals: &mut BinaryHeap<HitInterval>,
    mismatch_bound: &MB,
) where
    MB: MismatchBound,
{
    // TODO: Check performance impact
    // This is technically redundant. Our micro-benchmarks suggest
    // that having this here improves the performance but it might
    // be that actually the opposite is true for large real data.
    if let Some(best_scoring_interval) = intervals.peek() {
        if mismatch_bound.reject_iterative(
            stack_frame.alignment_score,
            best_scoring_interval.alignment_score,
        ) {
            return;
        }
    }

    stack_frame.edit_node_id = edit_tree
        .add_node(edit_operation, stack_frame.edit_node_id)
        .expect("We bound the length of `edit_tree` at `STACK_LIMIT` < `u32`");

    // This route through the read graph is finished successfully, push the interval
    if stack_frame.backward_index < 0 && stack_frame.forward_index > (pattern.len() as i16 - 1) {
        let edit_operations =
            extract_edit_operations(stack_frame.edit_node_id, edit_tree, pattern.len());
        intervals.push(HitInterval {
            interval: stack_frame.current_interval,
            alignment_score: stack_frame.alignment_score,
            edit_operations,
        });
        print_debug(&stack_frame, intervals, edit_tree); // FIXME
        return;
    }

    stack.push(stack_frame);
}

/// FIXME
fn print_debug(
    stack_frame: &MismatchSearchStackFrame,
    intervals: &BinaryHeap<HitInterval>,
    edit_tree: &Tree<EditOperation>,
) {
    let switch = false; // TODO: Switch me on/off!

    if switch {
        let best_as = match intervals.peek() {
            Some(v) => v.alignment_score,
            None => 0.0,
        };

        eprintln!(
            "{}\t{}\t{}\t{}\t{}\t{:?}",
            stack_frame.alignment_score,
            stack_frame.alignment_score,
            stack_frame.backward_index,
            stack_frame.forward_index,
            best_as,
            edit_tree.ancestors(stack_frame.edit_node_id).next(),
        );
    }
}

/// Finds all suffix array intervals for the current pattern
/// w.r.t. supplied alignment parameters
pub fn k_mismatch_search<SDM, MB>(
    pattern: &[u8],
    base_qualities: &[u8],
    parameters: &AlignmentParameters,
    fmd_index: &RtFmdIndex,
    stack: &mut MinMaxHeap<MismatchSearchStackFrame>,
    edit_tree: &mut Tree<EditOperation>,
    sequence_difference_model: &SDM,
    mismatch_bound: &MB,
) -> BinaryHeap<HitInterval>
where
    SDM: SequenceDifferenceModel,
    MB: MismatchBound,
{
    let center_of_read = pattern.len() / 2;
    let bi_d_array = BiDArray::new(
        pattern,
        base_qualities,
        center_of_read,
        parameters,
        fmd_index,
        sequence_difference_model,
    );

    let optimal_penalties =
        compute_optimal_scores(pattern, base_qualities, sequence_difference_model);
    let mut hit_intervals = BinaryHeap::new();

    let mut stack_size_limit_reported = false;
    stack.clear();
    let root_node = edit_tree.clear();

    stack.push(MismatchSearchStackFrame {
        current_interval: fmd_index.init_interval(),
        backward_index: center_of_read as i16 - 1,
        forward_index: center_of_read as i16,
        direction: Direction::Forward,
        gap_backwards: GapState::Closed,
        gap_forwards: GapState::Closed,
        alignment_score: 0.0,
        edit_node_id: root_node,
    });

    while let Some(stack_frame) = stack.pop_max() {
        // Determine direction of progress for next iteration on this stack frame
        let (
            j,
            next_backward_index,
            next_forward_index,
            fmd_ext_interval,
            optimal_penalty,
            next_insertion_backward,
            next_insertion_forward,
            next_deletion_backward,
            next_deletion_forward,
            insertion_score,
            deletion_score,
            next_closed_gap_backward,
            next_closed_gap_forward,
        );
        let mut mm_scores = [0_f32; 4];
        match stack_frame.direction {
            Direction::Forward => {
                next_forward_index = stack_frame.forward_index + 1;
                next_backward_index = stack_frame.backward_index;
                j = stack_frame.forward_index;
                optimal_penalty = optimal_penalties[j as usize];
                fmd_ext_interval = stack_frame.current_interval.swapped();
                next_insertion_backward = stack_frame.gap_backwards;
                next_insertion_forward = GapState::Insertion;
                next_deletion_backward = stack_frame.gap_backwards;
                next_deletion_forward = GapState::Deletion;
                next_closed_gap_backward = stack_frame.gap_backwards;
                next_closed_gap_forward = GapState::Closed;
                insertion_score = if stack_frame.gap_forwards == GapState::Insertion {
                    parameters.penalty_gap_extend
                } else {
                    parameters.penalty_gap_open
                } + stack_frame.alignment_score;
                deletion_score = if stack_frame.gap_forwards == GapState::Deletion {
                    parameters.penalty_gap_extend
                } else {
                    parameters.penalty_gap_open
                } + stack_frame.alignment_score;
                for (score, &base) in mm_scores.iter_mut().zip(b"ACGT".iter().rev()) {
                    *score = sequence_difference_model.get(
                        j as usize,
                        pattern.len(),
                        dna::complement(base),
                        pattern[j as usize],
                        base_qualities[j as usize],
                    ) - optimal_penalty
                        + stack_frame.alignment_score;
                }
            }
            Direction::Backward => {
                next_forward_index = stack_frame.forward_index;
                next_backward_index = stack_frame.backward_index - 1;
                j = stack_frame.backward_index;
                optimal_penalty = optimal_penalties[j as usize];
                fmd_ext_interval = stack_frame.current_interval;
                next_insertion_backward = GapState::Insertion;
                next_insertion_forward = stack_frame.gap_forwards;
                next_deletion_backward = GapState::Deletion;
                next_deletion_forward = stack_frame.gap_forwards;
                next_closed_gap_backward = GapState::Closed;
                next_closed_gap_forward = stack_frame.gap_forwards;
                insertion_score = if stack_frame.gap_backwards == GapState::Insertion {
                    parameters.penalty_gap_extend
                } else {
                    parameters.penalty_gap_open
                } + stack_frame.alignment_score;
                deletion_score = if stack_frame.gap_backwards == GapState::Deletion {
                    parameters.penalty_gap_extend
                } else {
                    parameters.penalty_gap_open
                } + stack_frame.alignment_score;
                for (score, &base) in mm_scores.iter_mut().zip(b"ACGT".iter().rev()) {
                    *score = sequence_difference_model.get(
                        j as usize,
                        pattern.len(),
                        base,
                        pattern[j as usize],
                        base_qualities[j as usize],
                    ) - optimal_penalty
                        + stack_frame.alignment_score;
                }
            }
        };

        // Calculate the lower bounds for extension
        let lower_bound = bi_d_array.get(next_backward_index, next_forward_index);

        print_debug(&stack_frame, &hit_intervals, edit_tree); // FIXME

        // Since we operate on a priority stack, we can assume that there are no
        // better scoring frames on the stack, so we are going to stop the search.
        if let Some(best_scoring_interval) = hit_intervals.peek() {
            if mismatch_bound.reject_iterative(
                stack_frame.alignment_score + lower_bound,
                best_scoring_interval.alignment_score,
            ) {
                break;
            }
        }

        //
        // Insertion in read / deletion in reference
        //
        {
            if !mismatch_bound.reject(insertion_score + lower_bound, pattern.len())
                && parameters.gap_dist_ends as i16 <= j.min(pattern.len() as i16 - j)
            {
                check_and_push_stack_frame(
                    MismatchSearchStackFrame {
                        backward_index: next_backward_index,
                        forward_index: next_forward_index,
                        direction: stack_frame.direction.reverse(),
                        // Mark opened gap at the corresponding end
                        gap_backwards: next_insertion_backward,
                        gap_forwards: next_insertion_forward,
                        alignment_score: insertion_score,
                        ..stack_frame
                    },
                    pattern,
                    EditOperation::Insertion(j as u16),
                    edit_tree,
                    stack,
                    &mut hit_intervals,
                    mismatch_bound,
                );
            }
        }

        // Bidirectional extension of the (hit) interval
        for ((c, mut interval_prime), mm_score) in fmd_index
            .extend_iter(&fmd_ext_interval)
            .zip(mm_scores.iter())
        {
            if interval_prime.size < 1 {
                continue;
            }

            // Special treatment of forward extension
            let c = match stack_frame.direction {
                Direction::Forward => {
                    interval_prime = interval_prime.swapped();
                    dna::complement(fmd_index.get_rev(c))
                }
                Direction::Backward => fmd_index.get_rev(c),
            };

            //
            // Deletion in read / insertion in reference
            //
            {
                if !mismatch_bound.reject(deletion_score + lower_bound, pattern.len())
                    && parameters.gap_dist_ends as i16 <= j.min(pattern.len() as i16 - j)
                {
                    check_and_push_stack_frame(
                        MismatchSearchStackFrame {
                            current_interval: interval_prime,
                            // Mark open gap at the corresponding end
                            gap_backwards: next_deletion_backward,
                            gap_forwards: next_deletion_forward,
                            alignment_score: deletion_score,
                            ..stack_frame
                        },
                        pattern,
                        EditOperation::Deletion(j as u16, c),
                        edit_tree,
                        stack,
                        &mut hit_intervals,
                        mismatch_bound,
                    );
                }
            }

            //
            // Match/mismatch
            //
            {
                if !mismatch_bound.reject(mm_score + lower_bound, pattern.len()) {
                    check_and_push_stack_frame(
                        MismatchSearchStackFrame {
                            current_interval: interval_prime,
                            backward_index: next_backward_index,
                            forward_index: next_forward_index,
                            direction: stack_frame.direction.reverse(),
                            // Mark closed gap at the corresponding end
                            gap_backwards: next_closed_gap_backward,
                            gap_forwards: next_closed_gap_forward,
                            alignment_score: *mm_score,
                            ..stack_frame
                        },
                        pattern,
                        if c == pattern[j as usize] {
                            EditOperation::Match(j as u16)
                        } else {
                            EditOperation::Mismatch(j as u16, c)
                        },
                        edit_tree,
                        stack,
                        &mut hit_intervals,
                        mismatch_bound,
                    );
                }
            }
        }

        // Only search until we've found a multi-hit
        // FIXME: We can probably stop searching earlier
        if (hit_intervals.len() > 9)
            || (hit_intervals.len() == 1 && hit_intervals.peek().unwrap().interval.size > 1)
        {
            return hit_intervals;
        }

        // Limit stack size
        if stack.len() > STACK_LIMIT as usize || edit_tree.len() > EDIT_TREE_LIMIT as usize {
            if parameters.stack_limit_abort {
                return hit_intervals;
            } else {
                if !stack_size_limit_reported {
                    trace!(
                    "Stack size limit exceeded (read length: {} bp). Remove highly penalized partial alignments from stack (stack size: {}, edit tree size: {}).",
                    pattern.len(),
                    stack.len(),
                    edit_tree.len(),
                );
                    stack_size_limit_reported = true;
                }

                for _ in 0..(stack.len() as isize - STACK_LIMIT as isize)
                    .max(edit_tree.len() as isize - EDIT_TREE_LIMIT as isize)
                {
                    // The stack should never be empty at this point, so we boldly ignore the
                    // `None` case here
                    if let Some(min_stack_frame) = stack.pop_min() {
                        edit_tree.remove(min_stack_frame.edit_node_id);
                    }
                }
            }
        }
    }
    hit_intervals
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::{mismatch_bounds::*, sequence_difference_models::*, utils::*};
    use assert_approx_eq::assert_approx_eq;
    use bio::alphabets;
    use noodles::{bam, sam};

    #[test]
    fn test_inexact_search() {
        let difference_model = TestDifferenceModel {
            deam_score: 0.5,
            mm_score: -1.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -1.0,
            representative_mm_bound: -1.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        let ref_seq = "ACGTACGTACGTACGT".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GTTC".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -1.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_eq!(alignment_score, vec![-1.0]);

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&suffix_array))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(positions, vec![2, 6, 10, 19, 23, 27]);
    }

    #[test]
    fn test_reverse_strand_search() {
        let difference_model = TestDifferenceModel {
            deam_score: -10.0,
            mm_score: -10.0,
            match_score: 1.0,
        };
        let mmb = TestBound {
            threshold: -1.0,
            representative_mm_bound: -10.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -20.0,
            penalty_gap_extend: -10.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        let ref_seq = "GAAAAG".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "TTTT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // 1.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&suffix_array))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(positions, vec![8]);
    }

    #[test]
    fn test_d() {
        let ref_seq = "GATTACA".as_bytes().to_owned(); // revcomp = TGTAATC

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
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

    #[test]
    fn test_gapped_alignment() {
        let difference_model = TestDifferenceModel {
            deam_score: -10.0,
            mm_score: -10.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -2.0,
            representative_mm_bound: -10.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        let ref_seq = "TAT".as_bytes().to_owned(); // revcomp = "ATA"

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "TT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -2.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&suffix_array))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(positions, vec![0, 2, 5]);
    }

    #[test]
    fn test_gapped_alignment_read_end() {
        let difference_model = TestDifferenceModel {
            deam_score: -10.0,
            mm_score: -10.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -5.0,
            representative_mm_bound: -10.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 5,
            stack_limit_abort: false,
        };

        let ref_seq = "AAAAAGGGGAAAAA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        // Gap in the middle of the read (allowed)
        let pattern = "AAAAAAAAAA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&suffix_array))
            .flatten()
            .collect();
        assert_eq!(positions, vec![0]);

        // Gap near read end (not allowed)
        let pattern = "AGGGAAAA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&suffix_array))
            .flatten()
            .collect();
        assert_eq!(positions, vec![]);
    }

    #[test]
    fn test_vindija_pwm_alignment() {
        let difference_model = VindijaPwm::new();
        let mmb = TestBound {
            threshold: -30.0,
            representative_mm_bound: difference_model.get_representative_mismatch_penalty(),
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            // Disable gaps
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -200.0,
            penalty_gap_extend: -100.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        let ref_seq = "CCCCCC".as_bytes().to_owned(); // revcomp = "ATA"

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "TTCCCT".as_bytes().to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -30.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_eq!(alignment_score[0], -4.641691);

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&sar))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(positions, vec![0]);

        let pattern = "CCCCCC".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -30.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_eq!(alignment_score[0], 0.0);

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&sar))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(positions, vec![0]);

        //
        // Test "normal" mismatch
        //

        let ref_seq = "AAAAAA".as_bytes().to_owned(); // revcomp = "ATA"

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "AAGAAA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -30.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_approx_eq!(alignment_score[0], -10.965062);
    }

    #[test]
    fn test_ord_impl() {
        let mut edit_tree: Tree<usize> = Tree::new();
        let root_id = edit_tree.clear();

        let map_params_large = MismatchSearchStackFrame {
            alignment_score: -5.0,
            current_interval: RtBiInterval {
                lower: 5,
                lower_rev: 5,
                size: 5,
            },
            backward_index: 5,
            forward_index: 5,
            direction: Direction::Backward,
            gap_forwards: GapState::Closed,
            gap_backwards: GapState::Closed,
            edit_node_id: root_id,
        };
        let map_params_small = MismatchSearchStackFrame {
            alignment_score: -20.0,
            current_interval: RtBiInterval {
                lower: 5,
                lower_rev: 5,
                size: 5,
            },
            backward_index: 5,
            forward_index: 5,
            direction: Direction::Backward,
            gap_forwards: GapState::Closed,
            gap_backwards: GapState::Closed,
            edit_node_id: root_id,
        };

        assert!(map_params_large > map_params_small);
        assert!(!(map_params_large < map_params_small));
        assert_ne!(map_params_large, map_params_small);
    }

    #[test]
    fn test_corner_cases() {
        let difference_model = VindijaPwm::new();
        let repr_mm_penalty = difference_model.get_representative_mismatch_penalty();
        let mmb = Discrete::new(0.01, 0.02, repr_mm_penalty);

        let parameters = AlignmentParameters {
            penalty_gap_open: 3.5 * difference_model.get_representative_mismatch_penalty(),
            penalty_gap_extend: 1.5 * difference_model.get_representative_mismatch_penalty(),
            difference_model: difference_model.clone().into(),
            chunk_size: 1,
            mismatch_bound: mmb.clone().into(),
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        // "correct" "AAAAAAAAAAAAAAAAAAAA" (20x 'A') "incorrect"
        let ref_seq = "GTTGTATTTTTAGTAGAGACAGGGTTTCATCATGTTGGCCAGAAAAAAAAAAAAAAAAAAAATTTGTATTTTTAGTAGAGACAGGCTTTCATCATGTTGGCCAG"
            .as_bytes()
            .to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GTTGTATTTTTAGTAGAGACAGGCTTTCATCATGTTGGCCAG"
            .as_bytes()
            .to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let alignment_scores = intervals
            .iter()
            .map(|f| f.alignment_score)
            .collect::<Vec<_>>();
        assert_eq!(alignment_scores, vec![-10.936638, -38.377716, -10.965062]);

        let mut positions: Vec<usize> = intervals
            .iter()
            .map(|f| f.interval.forward().occ(&sar))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(positions, vec![0, 62, 63]);

        assert_eq!(
            intervals.peek().unwrap().interval.forward().occ(&sar),
            vec![0]
        );
    }

    #[test]
    fn test_cigar_indels() {
        let difference_model = TestDifferenceModel {
            deam_score: -10.0,
            mm_score: -10.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -3.0,
            representative_mm_bound: -10.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        //
        // Deletion
        //
        let ref_seq = "GATTAGCA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "ATTACA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit.edit_operations.to_bam_fields(Direction::Forward);

        assert_eq!(
            cigar,
            vec![
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Match, 4).unwrap(),
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Deletion, 1).unwrap(),
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Match, 2).unwrap(),
            ]
        );

        //
        // 2-base deletion
        //
        let ref_seq = "GATTACAG".as_bytes().to_owned(); // CTGTAATC

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATCAG".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit.edit_operations.to_bam_fields(Direction::Forward);

        assert_eq!(best_hit.alignment_score, -3.0);
        assert_eq!(
            cigar,
            vec![
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Match, 3).unwrap(),
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Deletion, 2).unwrap(),
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Match, 3).unwrap(),
            ]
        );

        //
        // Insertion
        //
        let ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATTAGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit.edit_operations.to_bam_fields(Direction::Forward);

        assert_eq!(best_hit.alignment_score, -2.0);
        assert_eq!(
            cigar,
            vec![
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Match, 5).unwrap(),
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Insertion, 1).unwrap(),
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Match, 2).unwrap(),
            ]
        );

        //
        // 2-base insertion
        //
        let ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATTAGGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit.edit_operations.to_bam_fields(Direction::Forward);

        assert_eq!(best_hit.alignment_score, -3.0);
        assert_eq!(
            cigar,
            vec![
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Match, 5).unwrap(),
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Insertion, 2).unwrap(),
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Match, 2).unwrap(),
            ]
        );

        //
        // 3-base insertion
        //
        let ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATTAGTGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mmb = TestBound {
            threshold: -4.0,
            representative_mm_bound: difference_model.get_representative_mismatch_penalty(),
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -4.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit.edit_operations.to_bam_fields(Direction::Forward);

        assert_eq!(best_hit.alignment_score, -4.0);
        assert_eq!(
            cigar,
            vec![
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Match, 5).unwrap(),
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Insertion, 3).unwrap(),
                bam::record::cigar::Op::new(sam::record::cigar::op::Kind::Match, 2).unwrap(),
            ]
        );
    }

    #[test]
    fn test_md_tag() {
        let difference_model = TestDifferenceModel {
            deam_score: -1.0,
            mm_score: -2.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -1.0,
            representative_mm_bound: -2.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        //
        // Mutation
        //
        let ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATTATA".as_bytes().to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -1.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit.edit_operations.to_bam_fields(Direction::Forward);

        assert_eq!(md_tag, "5C1".as_bytes());

        //
        // Deletion
        //
        let ref_seq = "GATTAGCA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "ATTACA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mmb = TestBound {
            threshold: -3.0,
            representative_mm_bound: difference_model.get_representative_mismatch_penalty(),
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit.edit_operations.to_bam_fields(Direction::Forward);

        assert_eq!(md_tag, "4^G2".as_bytes());

        //
        // 2-base deletion
        //
        let ref_seq = "GATTACAG".as_bytes().to_owned(); // CTGTAATC

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATCAG".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit.edit_operations.to_bam_fields(Direction::Forward);

        assert_eq!(md_tag, "3^TA3".as_bytes());

        //
        // Insertion
        //
        let ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATTAGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit.edit_operations.to_bam_fields(Direction::Forward);

        assert_eq!(md_tag, "7".as_bytes());

        //
        // 2-base insertion
        //
        let ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATTAGGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit.edit_operations.to_bam_fields(Direction::Forward);

        assert_eq!(md_tag, "7".as_bytes());
    }

    #[test]
    fn test_reverse_strand_search_2() {
        let difference_model = TestDifferenceModel {
            deam_score: -1.0,
            mm_score: -1.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: 0.0,
            representative_mm_bound: -1.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -3.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        let ref_seq = "AAAGCGTTTGCG".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "TTT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // 0.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let positions = if let Some(best_alignment) = intervals.pop() {
            best_alignment
                .interval
                .forward()
                .occ(&suffix_array)
                .iter()
                .filter(|&&position| position < suffix_array.len() / 2)
                .map(|&position| (position, Direction::Forward))
                .chain(
                    best_alignment
                        .interval
                        .revcomp()
                        .occ(&suffix_array)
                        .iter()
                        .filter(|&&position| position < suffix_array.len() / 2)
                        .map(|&position| (position, Direction::Backward)),
                )
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        assert_eq!(
            positions,
            vec![(6, Direction::Forward), (0, Direction::Backward),]
        );
    }

    #[test]
    fn test_edit_operations_reverse_strand() {
        let difference_model = TestDifferenceModel {
            deam_score: -1.0,
            mm_score: -1.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -1.0,
            representative_mm_bound: -1.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -3.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        let ref_seq = "GATTACA".as_bytes().to_owned(); // revcomp = "TGTAATC"

        // Reference
        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "TAGT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -1.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let best_alignment = intervals.peek().unwrap();

        let positions = best_alignment
            .interval
            .forward()
            .occ(&suffix_array)
            .iter()
            .filter(|&&position| position < suffix_array.len() / 2)
            .map(|&position| (position, Direction::Forward))
            .chain(
                best_alignment
                    .interval
                    .revcomp()
                    .occ(&suffix_array)
                    .iter()
                    .filter(|&&position| position < suffix_array.len() / 2)
                    .map(|&position| (position, Direction::Backward)),
            )
            .collect::<Vec<_>>();
        assert_eq!(positions, vec![(1, Direction::Backward),]);

        let (_cigar, md, edop) = best_alignment
            .edit_operations
            .to_bam_fields(Direction::Backward);
        assert_eq!(md, b"1T2");

        assert_eq!(edop, 1);
    }

    #[test]
    fn test_n() {
        let ref_seq = "GATTACAGATTACAGATTACA".as_bytes().to_owned();

        let difference_model = SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            },
            0.001,
            0.9,
            0.02 / 3.0,
            false,
        );

        let representative_mismatch_penalty =
            difference_model.get_representative_mismatch_penalty();

        let mismatch_bound = TestBound {
            threshold: -14.0,
            representative_mm_bound: difference_model.get_representative_mismatch_penalty(),
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mismatch_bound.clone().into(),
            penalty_gap_open: 0.001_f32.log2(),
            penalty_gap_extend: representative_mismatch_penalty,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
        };

        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        {
            let pattern = "NNNNNNNNNN".as_bytes().to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 0);
        }

        {
            let pattern = "AGATNACAG".as_bytes().to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 1);
        }
    }

    #[test]
    fn test_bench() {
        let ref_seq = "AACATCTTTTTCGGCCGTTTCATAAATACTCTTTATTCACAATTAAGGTCTCTCTGTGCTTCTACGATGGTTGGCTGGCCTGCTACTATGTAGCCCAGTT\
TTGTACTTGAAACAGATGAGTCTCTTCACACTCATACGCGTTTCAGCAAGTACCAGTGATAGCGGTGGTGCTCCTCAATATCCGCCGCTCGCAAGCGACC\
GAGAAAAGCACGCATCGAGAACCTTGTCTTACACAACTTTCAATGCTAGCGCCTAAGTTCAACCGCTAATATAACAGCTGAGTCCCCATCGCCCGAGGGT\
CGCATAGTTACTCCGTGAGGTGTATCAGAAATCCGCGATTTCGTAACTTCTACATCCCGTCCTCTCAGACCACCATTTTCCGTCATTACCAGCAAATATG\
ATCCGCCGATCGGCCGGACATCTGCCAAACCTTGCTGTCCTCCAATAAGTTTAGGTTAAGGGCCCAGTCAAGAAAAGACCCCAGTTCAAATGACGTACAA\
GGCAGCACCCCGTTTATACACCGGGGCATGCTTAGCTGCTATCGGACTTGCGCGGATAAAAATAGCCTAACTGCTAAAAGCATGTCCGCCCACGATGGAC\
TGGGGATTAACAGATAAGTTTAACATCTTAAACGTGAGTTGCGACAACAGCCAGGCTGACCGTACTGCCCGGACGTGCGGGGTCTGCTCATTGTGAATTC\
AGAGATGCCTCTATATTAGATTGGAGGGTAACTCAATAATCATGAGCGCAAGCTCACGAGAAGTACTGCTACTCCGCATATGTTGCTGAATATCTTCTGC\
CTGCGTGGCGGCCAGCATGTGATTACAATGACTCACTAGGGGGGGTAGTAGAGTCCGATTCGGTTGACTACTCGGGGTAACTAGGCATACAGTGGGCTCT\
AAGTCTGTGCGTCAAGCCCAGGGCAGGTCGTCGCTCATTAAGTTAGCCCACTTCATTAAGTATCACAGGGGCATCCTAGGTCCGTCTTAGTGTCCAACCA\
GCACCTCAATGTAATGTGGTGAGGGTAGTCACCAGTAGCGGAGTCCATCCTCTGGCCGCGTCGAGTCTTTCACGTCTATTCATGACGTATAGTCCCATTG\
GTGCGTCATGTACAGGCTTCTAGGACCACATTCTCCTTGGGCGGTTATCAAGCAAAATGCAAAGGTACCTCGGCATTCCGCCGATATTGTCTCTGAGATG\
TACGACGAGCCACACACTACGCTCGTTCTTCTGAGCCAATGGTTCGTATGGGATGCACTATAGTCGGACTAGACCCTTGCCATACTACATGCTGAAGTTT\
CAGGGATAACAATCGGTGGGCATAGACAACGCTAACGTCCGCGTTAAAGATTTAGTAACGTACGATACGCTGATTGCTGCAGGACAGAATAGACTTCAAT\
CGCGCTACCCTTCTAGCGCTGGGTAGTTTCGCGGACCGTATGCCAGGTAACAAGCACAGTAAACGCACTATTGGGTAATACGGATCTCAAAGGCTCATGC\
TCGGTGAAGGTGCCCGTAAGTGCACTTTCCTAACATAGCCGGCGTTATTAAGAGTTCCATCGGACAGCGCCCCACACCAAATCACTCAATTTCTCCTGTT\
CACCTTGCGCGTCATAATCTAACTGTGGCGAAGAACTAGGAACCTGAGAATAACCATTTCGATCGGTGATCCGTCAGATCGCTATTAGGAGCCAGAGACG\
TCTTCTACGGGGTCACCGTTCGTACCATGCTACCAGTCTTCGTCCCCTATGAGCATTTAACACAATCGTACTGTAGGTTGTCACACAACCGCTTAACGTT\
GTATCAGCTTTGCATACGCGAAGCAACGATGAGAACCGCTAAGGCTTCAATCTCCCGTTAAAGACCCTGAACTGGTCCGTGGCGGCGTGTTTTATCAGTG\
TACGGGATTCTTCGTGTACTTACTTGATAACGGCTCATACCAAGGCACGTGATTCTCTCGCAGAGGTCTCCCCCTTGGAAGTTATAGGCGCGCCGGTTTC\
CATTAGCATTGGACCACGACCGGGGGTGCTCCTATGTCCGACTCTCCATAGCGAGGAACATTGCTTCAAGCCCTTGGACGCTAGTTATTCATGCCTTGCG\
ATCTCCTAGCTTACCAAGGTTGTCGAACATCCGTTGAGAAAGTCATAAGACGTCAATAAACGTAACGGGGTTCTTTCTGCACTTGCGCTCAGACCTGCAT\
TTTAACCGTACAGGGAGCATCATACATCCCGTAACCACAGATATATCACTAACCAGGAATGACAAAGTCGGAAAGTAAGTCCCTACCGTCCGGAACCAGC\
GACGACCCAGGAAACTCCAGACTAGAAAGAGTACAGCAACCATGATATAAATAGTGCGAGCTGACAATGGTTTCGGCAGCACTGAATTTTTCGGGCTCAC\
AACAGGCGACAACCGGCCTTCGGAGAGTTTAATGCGCCATTCAGTTCGTCATCACGCTGGTACCGTCACCCTATGGGAGCTGCCAACTGTCACGATAGAA\
GATTATTTGTCACGGGACCGATTGCTTAAGACAATAACCAGGCTGGGACTACGTCATACGTCTCCAACTAACAAAAGTACATCCCCTATGAGATCAAAGG\
TACCGTATAAGAGAGTATGTCTAGGGCCTTGCATGGTTTTTACTTAACGCCGCTCCCGGAGCTTGCCCGCCATCCGGTCGAGAAAGCGTACCTCAAGTTT\
TGACCCTAACCAACACCCTAGAAGTTATTTGCTTCTAGGTTTCATACTTAAACTGCAATTGACACCGTACAAGCAGTCCGACCACTTTGGAAGCGTCGTG\
GATGCCCGGTTATACGGGCGATCGAGCTCTAGTCCTGACTCACGGCCGTGTATAACTTGATTCACGCCCCTTGCGCCGAACGCCCTATGGATCAAGGCGC\
CATCGAAGAGGTATTGTTACCATCAGAGGTTCGTGTGTCCTTTTAGCGCAACATTGATTACCGGTGCCTCCGTAAGTGGGGCCGACCAAAGTGCCGCCAG\
ATACTACCCCGCACGAAAATAGAAATGCGTTTGCGACCCTAACACTGAGCGTGCTGCCACTTGAAATTAATACTCTGAGGGCGCGGATCACTACCGGCGA\
GGCAAGATTACAATACGTACACTCCGTTAATACCATAACTCTTGGGTGTAGGAACTCCGGTAATACACCCGCAGGACGTTATTGGCCGTCGCCAGCAGCT\
TAGTAACAGATATAAAATCTAGGGACTTGTCGGGCAATGTCGGACGAATTCCGTTGTACCGTCAACTCCCTATACGCAAACTTGGTAACGTTTACTAGGT\
AGACGCGTTCAATGAACGGGTCTACTACTTCACTGGCGCCTCGTTGTGGCGTACTGCACCCGTAATACCCTATACCCGGCGGTATTCATCAGTCGCTCTG\
TCCACAGCAGCCAGGGAGGCTGTTTTCCTCGGGGCTACGGGGTCAACCCGATCCCATGGGTTTGCGTAGTACAATGAAAGGTCTCCTATCTTTCTTTTTG\
CGTACGGCTGCAGGGAAGACGCGGACATCATAGGCAACAGGGCCAGAGATATCCCCGGCTTCGAGCACTAGGGTGCTATTCCCCCAACTACGGCGCCAGG\
ACACGTGCCACTTGTTGCCACCTAACATCACTGTTGACCGTCAGATTGAAGAGTACTCAATGTCCACGATTTGCGCTGTGCCAGCGGCATCAGAACGGGA\
ATACGCAGGCATGCGTGAGAGCCAAATACCGCCCTCTGTCGGCAATCAGGTCGTAAAAGTTGAGATGTGAGGGTCTTAAGTTTGCGCCTTATGGTCTGGA\
CTCCGGTAGAACTCGGACGTTTAAGGTAACTTAGGAGGCCCTCACCTTCGAGACTGGCCGAGCAGTCAAAATCGGTCGGGCCGACGTCTTATAACTTGGC\
AAAGGTGAGGGGCACGTCTGTAATATGCAGCTGACATGTGGCTAAGAGTTGTTAAAATGTGTCGTTAAGTGTTAATACCTCACATTAGGTTTTAACGACC\
ACCTTGTTAGGTTTAAATCGGCAACTTGGCATGGACTATGTACTGGAAATCGACATAACGTGGTGAAGACTCCATACCTCGACGCATACGACCCCTCAAT\
ATCCCGCATCGTGTCCATTTATCGCTAAATAGGAAAGTGACCGACGCATACGACCGCCAGGCGCACACAAAGCTGAATCCGGCCAATAAAATAGTCTCAT\
GGTCCCGCCCGAGCCCAGTCTCGGCACCCCTAGACCCGAGTGCCGATCTCGTAGGGACACACCCTGATCGATTTTTTGCTGCATGATCTTTAAAAGGCAG\
ATGGGGGCTTAACAGTCACTACACTGACACATGAGACCTCGGAAGCGGAATCACTGCGTTCCAAACTCTCAGAGACCACGGCAAAGGAGATTTCCGGGGC\
GTAGGCCTTTTCGGGATACACTACGTATGGTTCAGCGCGACAGCCCGCATATACGGTGCATTTATATCAATCTTTAGGGAGTAACCCATCAGGTGAACAT\
GTCCAGCCCCCGGGGTCGCGAGCAATTAGCCGTGATCGGATAGGTCCAAACATGACTGATGTAGGGCCACCGGTGGAACCAGCCCGAGCAATGATTCGCC\
CCGGTGACTCCTACTCAAGCGGTCGGTGAGAGTATTGACCTTGTGACCAACGACGGGAGCAGCAAGCTGCGTGTTCCTTACCAACTCCGGGATATTCCTT\
TGATGTTCTCTAAGACAGACGTTGAGCATCAGCCCATATTAAGTCGAACGTGGTCCACTATCGATAAGTCTGCACAATTGGACGCTAGCAGGAATGGTGT\
GAGGACCCGAATCAGGCCCTCGAGTGGGACCTGGCGATCACTGCCGTCTCCGTTGCCGGTCGTTGTCCAACCGGTAAACCGGACAGCCCACATGTGCACG\
GTAAATGTTCCGTGCATAACATAATGCGAATGACCGCCTTCTTATGTCCGTTCCGAGGCCGGTCTTCTCACCCTTTTAGCATGTGGCAAGCCGGCACCGG\
AAGTGGACGACCGTCGCGTCTGTTCGTTGAGCGTGTATCTGGCGACTTAAACGATCCTACTAATCGCGTACGCTTGGTGCCTATAGAACCCGGGGAACTT\
GGCACGTCATGCCGGATACTTAATAGCCTGGCTCGTGTCGTGAGAACAATCTCCTGTGAAGTTACGTCTGAGTTTTAAAAGGCATCGTCGATCCCCTTCC\
CCTGCGCGCGGAGACGACTTAGCGACGTCGTGCACGACGAAAGCAGCCCTATTATAAGTCGTCCAGAAAGCCAGTATCTGATGTAAGCCAGACGCTAAGT\
TCAACCCAGTACATTTAAGACTAACCCATGTAAGCAAGTTATGTTGATCGCCCTCACAGCTATTCCCCGGCGGTGGACTCCATTCCATGCAGCTGATATA\
CAGATACATCACGCGTTAATTATGGTCTAGGTGTGTAAGCCGTTTCTAGGTTAGGCAGCAGTTAACGGCACTAGAGTATTGAAACACGGGAAACTAGGGG\
GCATCAATTTTCAGATGTCGTACCACAATTGCGAAATGGAACTGCAATTGTATCAAATGGCGCAATAGATATGTAGTGTCGTCCTGAGAATTCATTGGCA\
GTAAATCAACCTACTGAGCGGAGGATGTATACGCCAAATACATACATTAAAGTATCAATTATGCCCGGCCGAGTCTGGACAACCCGAGAAGCGTGATTTT\
GGGGAGTGGGGCGACGAGTCAGTAGCGTTTGGGCAGAACATGGCCTCCCTAGAGCCGCGACCGCGTGCTGGTACGGCACTCAAACTACATTTCTTTTGAG\
CGCGTAATATGAGAATTACACGGGTGCTAGATATGTAAGATACCACCAACTTCCCGCGTCCAAGGAACGGGAGACCACCTCCCGGGCCGTGCTGAGCATT\
GCCCGTTCAAGGCACCAGCAATTCTGGTGCCATAGTCGTAGAGCAGCCAACAACGGTGAAGAAGCTTACCGAAATTACAATGCGAGTGTGGCCTGCCGGC\
ACCTAAGTGTCAGCGGCCGGGGGGAGAGGATATTTTGCCGCATCCCGTCCCAAACGCGTGTTACGATCTACTGCAGATATTTGATTAGCGGTATCTTTTG\
CCCCTCCCCGCTTCAAGTTTTCATCCTAAGGGTGGGACACTTACAAGCTACAGGCCCCTTGGTAACAGGCTGTACACATCACACGTTGGAAAGCGTAATT\
ATAAAGGGTGTAGCCCTTATCTACTACTTCCAGACCGGAGTTTAACAGAGACCTGAAAGGACACTGCATTTACGAAGGGATGCCCGGGCCAACGAGAATG\
CGAACGGGTTATCGACTTAACGCGTAGTCTTCCGCTGCGGTAATTGTGTTGGATCTCCAGCTATACAGGCTGGAGCTCACTTGCCAGTCGGTAGGGTGAG\
CGGCCTTTATAGAGCGTGAGTTAGCATTGCCATTGAAATGAGATTGCTTGCCCAATACAAGTCTGTCTTCCAAGTAATAGGGCATTTGAGGCTTGAGGAA\
GCCACGCTCAGCTCTGGATGAAAAAAGTAACACGCTTCTATTAAGAAAAAGGTGAAGAAACTCCGCCTCGCAGTAGGTATCCAGCCGAACAGGGCAATCT\
CTGAGGGGACTTGTTAGACCGTCTCTGGGACTCCCATTGGCAGCAACAAAGGGATGGATCACATCTGAAGCGTCGCGGCTATTAAGAGTCTCAGCCCGGG\
CTATGCCGTATGTTTACCAAACCTCTGGTACCCATATTTCGTAGTGATGTCGCATAAAGCGCTCCTAGCAAGTTCGATCACACGTGACGAGGAGGCACGT\
TCCGCGTCATACAGACCGGGGGTGTTACTTCATTGGTACGACGGGGCACGACTCACCGAATGGGCCTGGTATTCTTCGGACGGTAGGTCCATCTACCATG\
CTCGCCACGTTTACACGGAACTTCGGGCAGCGACTCGCTAAGGGCTGCGTTAGAAACTACCTTGCAATCGGCCCAATGACTGGCAGAAATGAATCAGAAA\
CCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTAGGGTGACCGTAAGATTTAGGTCCATCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC\
TAACCGATTCTTAACTGCCAGTCCTCCGAAATTACCACCTGTTGGGGGGCATGATTCAGAATATGAAGTAGGATGGCCATGCGAAGAGGTGGCGCTTCAT\
GCTATCATTAGTACGGGATCTATAATCGAACACGATTTGTAATAGGGCGGGAGTCAAGCGCATATCCCACGTGGTTCCAACCGCTAGCCGCCTAACAAGT\
CTGTGATGACCAATATTGGAACCTTAGTACGTGGGGCCTTGAAAATCATTTGTTCGCTGAATAATTTACTGAAGATACCAGGTTTGTCATGATCAAGGTA\
GACCCTTGCCGACTCAATCCAATCTCCGGACCTAAAGGTGGCGATCTTAAGTCACCCACTTCGGCCTTGTTAATCGCGGCGCCCTGCATCCTTGGCTTTC\
AAGGAGGGGACATAAATGATTCGCTCGCACAAACGCAAGCTGGGCGTCTCAGGATTGCCGTTATTTGCAGCTTTAGCGGTACCATGCGGCGCCTAAATTT\
AACGAACTCTCGCCGGGGTCCCTCGTTCAATGTGACTATAGCTGGCTAAGGACGTTATCAGGCATCTCCAGGTCGTATCTCCTGGCTTGCTCTATTATGA\
CCTCCCGAGTTCATTACTGCCACCAAGCGTGTTATGCCAAGTAAGGACCAACGGGCCGCTTTAATGTACTCATATCGACGCGGTCTCTATGCTATCCACA\
ACGGCATTATTAGTAATCTACAACTATGCAATGAGTTGACCTTTGGTCCACGGACCAATAACCAGGGGGCTCTATACTGGCTCCACCCATCTTTTGCATG\
TTGGCCCAAATAATGCGCCCTTCCTAACTGGATACTTCTCCTTCCTTTGCCATACTGTAAACACGACTCGGCGGAGACGAGACATACCCAAGGTGGCCCA\
ATCCACAATTAACTAAAACCCAACCGCAGGTTAGCCGCCCTATATTGTTTCGTATCGCTTTAGTGTTGGAGGGTTATATCGCTGAACACATTATTATTGC\
GAACGTTCGGTTTCCTTCCGGTCGAAAAATAGTTACTTCATATGGGCCCCACGACCGATGGAGACCCACTGATACATTCTCGACTTGATCGGTCCGGGGT\
AAGGACAACAATGTGGCCAGGCCTTAATGCGATGGGACGATCATAGGGGTAGTATTGAGCAAGTTATACTTACTTACTTGCGGTCATTCGGGTTAAAGAG\
CTTTTGCTTTTGTGAGGCAACTGTTACGTGTATCGGACATGCATCCGATGGTTGCACCAATTTTAAAATTTGGACCTAATTGTGGGAGGCGGCTGGTTAT\
TCATGGAGCCCCTCTATATGGCACGTCGAACTCGTTAGTGCTTCGTCCGGTCCTAGACCAGGGAATACCGGTGAGCAAAGATACCCCACCCGATACGCCT\
CACACATCTCTTTCCGACTGTAACGAGCTGAAGCTAGAGTTCAGGTAAAGTCGGGCTAAGGACGTACATTAAAGCGGGGTACAATTCTGCACGCACCACC\
GTCTCAGATAGGCCGCCCAATCCAAGTTCGGGAAGTGATGGGATGACTACGGCTGCAATACTGTCTGGGATCTCTAACGCGATCACCCTCTCCATCCCTT\
AGGAAAGTCATCTTCTTTCCTAAACTGCCGCTGCTTATTGGACAACCGGTGACCGACAGCGAAACCATTGTTTGTGTAGTCTAGTAGACGACCAGCTGTC\
TGGGGTAATAACCAAAATCCGAATTCCACTATCCGGGGAGTCATAGCATACTGGGGTTAATCGAGCGCAGTCGAACGCAAAACTACACACCGAGTGAGGA\
CTCTCTCCAACGACCCTAATAGACTTCTGGAGCATCTTGCGAATAACGGACGGTATGTCGATGCTACCCCGCCTTCTGCTAGCCGGCATCGGCAACACGA\
TGAGCGTTCACGAAATTATTGAGAGAACAATTGCCTTAATACTAATTCGGGATATAGTTCAAAATGGACTGATATTCTCCACTATCTCTTACTAGACCCC\
ATTTCGGATAGAAGGCCATCGGAATAAAGCTGCACGAGGATCCGTCCCGACACAATCGGAGTAGAGGCTGGACATTCGATACAGGCCGTGGGCTAGCGAC\
CTCGCTCGAGCTCAAACGTTTCAGGCTGTCGGAACCCAGACGCCCCCCGAGGCCTCTGTGGTGGCAATATCCACAATAGTGCTTAGCGTCGATCCTCTGG\
CCATTCCCGTAACCCATCATAGTAACTGTGCGGATTCTCATCCTATGGGCAATAAAGCTCCAGGTCGTTCCTTATTTCCTCTCTAACTCCCCGTTGTGTC\
GTCGTTGGACGGCCTTCGTGCCGTGTTGGAATGTAAGATAAATCACGGGGCCATAATGCATTAATATTGGACGAGCTGTTCGCAAACACGTGAAGTAAAT\
AGAATTGTAATCGAAATTGTCAAGGACTATCTGATTACCAAACTACCGACTCCCCTTCACTGCGGGCAGCTCAGGTTGTCTGTACCGGGATGATGCAAGG\
CGTTCCCCATTCAGCTTGGAGGTGGACGCAGGTTGCCAATCGCCACTCGCGCTCTAGACGTTGCATCAGGTAGGGGGGGGCTTGTTTCATAACCTAGGGG\
TTCGTGAGGGTGGATTGGAATACTGACTATCACCTAGCAATGTTTTGCAGTATTATTTTAAGATCGTGTGCGGGGTGGGGAGCTATGTTTAAATAGCCTC\
GGCTCAACTTGCCAGCGTCCCATGAACTGTACAACCTAATCGTCCCTATCGTTACTGACCGACTGCGATAGAACGCCACCTTCATAGCGTCGCCGGACAA\
GCCTGTATGCAACCCATGAGTTTCCTTCGACTAGATCCAAACTCGAGGAGGTCATGGCGAGTCAAATTGTATATCTAGCGCCCACCTGATACCTAGGTTC\
".as_bytes().to_owned();

        let difference_model = SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            },
            0.001,
            0.9,
            0.02 / 3.0,
            false,
        );

        let representative_mismatch_penalty =
            difference_model.get_representative_mismatch_penalty();

        let mismatch_bound = Discrete::new(0.04, 0.02, representative_mismatch_penalty);

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mismatch_bound.clone().into(),
            penalty_gap_open: 0.00001_f32.log2(),
            penalty_gap_extend: representative_mismatch_penalty,
            chunk_size: 1,
            gap_dist_ends: 5,
            stack_limit_abort: false,
        };

        let alphabet = alphabets::Alphabet::new(crate::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _sar) = build_auxiliary_structures(ref_seq, alphabet);

        // bench_exogenous_read
        {
            let pattern = "GATATCTCGGCTGACAAACCAACAAAAAGTATCGGAACATCGCGGCGGCGTAGATGAATCTTAACCACACTCGACAGCTGTGCTTCTATACTAGCATTAC"
            .as_bytes()
            .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 0);
        }

        // bench_exogenous_read_deam
        {
            // 3' terminal C->T
            let pattern = "TTTATCTCGGCTGACAAACCAACAAAAAGTATCGGAACATCGCGGCGGCGTAGATGAATCTTAACCACACTCGACAGCTGTGCTTCTATACTAGCATTTT"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 0);
        }

        // bench_endogenous_read_perfect
        {
            let pattern = "CCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTAGGGTGACCGTAAGATTTAGGTCCATCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 1);
        }

        // bench_endogenous_read_1_mm_center
        {
            let pattern = "CCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTAGGGTGAGCGTAAGATTTAGGTCCATCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 1);
        }

        // bench_endogenous_read_2_mm_center
        {
            let pattern = "CCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTTGGGTGACCGTAAGATTTAGGTCCGTCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 1);
        }

        // bench_endogenous_read_1_deam
        {
            let pattern = "TCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTTGGGTGACCGTAAGATTTAGGTCCGTCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 1);
        }

        // bench_endogenous_read_4_deam
        {
            let pattern = "TCTGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTTGGGTGACCGTAAGATTTAGGTCCGTCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGTT"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 1);
        }
    }

    #[test]
    fn test_edop_effective_len() {
        let edop_track = EditOperationsTrack(vec![
            EditOperation::Match(0),
            EditOperation::Mismatch(1, b'C'),
            EditOperation::Match(2),
            EditOperation::Insertion(3),
            EditOperation::Match(4),
            EditOperation::Deletion(5, b'A'),
            EditOperation::Deletion(6, b'G'),
            EditOperation::Match(7),
            EditOperation::Match(8),
            EditOperation::Match(9),
            EditOperation::Match(10),
            EditOperation::Insertion(11),
            EditOperation::Mismatch(10, b'C'),
        ]);
        assert_eq!(edop_track.effective_len(), 11);

        let edop_track_2 = EditOperationsTrack(vec![
            EditOperation::Insertion(0),
            EditOperation::Insertion(1),
            EditOperation::Insertion(2),
        ]);
        assert_eq!(edop_track_2.effective_len(), 0);

        let edop_track_3 = EditOperationsTrack(vec![
            EditOperation::Deletion(0, b'A'),
            EditOperation::Deletion(1, b'C'),
            EditOperation::Deletion(2, b'G'),
            EditOperation::Deletion(3, b'T'),
        ]);
        assert_eq!(edop_track_3.effective_len(), 4);
    }

    #[test]
    fn stack_frame_size() {
        assert_eq!(40, std::mem::size_of::<MismatchSearchStackFrame>());
    }
}
