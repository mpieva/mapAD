use std::{
    cell::RefCell,
    cmp::Ordering,
    collections::{binary_heap::BinaryHeap, BTreeMap},
    error::Error,
    fs::File,
    iter::Peekable,
    time::{Duration, Instant},
};

use backtrack_tree::{NodeId, Tree};
use clap::{crate_name, crate_version};
use either::Either;
use log::{debug, trace};
use min_max_heap::MinMaxHeap;
use rand::{seq::IteratorRandom, RngCore};
use rayon::prelude::*;
use smallvec::SmallVec;

use bio::{
    alphabets::dna,
    data_structures::{
        bwt::{Less, Occ, BWT},
        fmindex::{BiInterval, FMDIndex, FMIndexable},
        suffix_array::{RawSuffixArray, SuffixArray},
    },
};

use bincode;
use rust_htslib::{bam, bam::Read};
use serde::{Deserialize, Serialize};
use snap;

use crate::{
    mismatch_bounds::MismatchBound,
    sequence_difference_models::{SequenceDifferenceModel, SequenceDifferenceModelDispatch},
    utils::{load_index_from_path, AlignmentParameters, Record},
};

pub const CRATE_NAME: &str = "mapAD";
pub const STACK_LIMIT: usize = 5_000_000;

/// Derive serde Traits for remote Struct
#[derive(Serialize, Deserialize)]
#[serde(remote = "BiInterval")]
struct BiIntervalDef {
    pub lower: usize,
    pub lower_rev: usize,
    pub size: usize,
    pub match_size: usize,
}

/// A subset of MismatchSearchStackFrame to store hits
#[derive(Serialize, Deserialize, Debug)]
pub struct HitInterval {
    #[serde(with = "BiIntervalDef")]
    interval: BiInterval,
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

    fn is_forward(self) -> bool {
        use self::Direction::*;
        match self {
            Forward => true,
            Backward => false,
        }
    }

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

/// Contains edit operations performed in order to align the sequence
#[derive(Debug, Serialize, Deserialize)]
pub struct EditOperationsTrack(Vec<EditOperation>);

impl EditOperationsTrack {
    pub fn effective_len(&self) -> usize {
        self.0
            .iter()
            .fold(0, |acc, edit_operation| {
                acc + match edit_operation {
                    EditOperation::Insertion(_) => 1,
                    EditOperation::Deletion(_, _) => -1,
                    EditOperation::Match(_) => 1,
                    EditOperation::Mismatch(_, _) => 1,
                }
            })
            .min(0) as usize
    }

    /// Constructs CIGAR, MD tag, and edit distance from correctly ordered track of edit operations and yields them as a tuple
    /// The strand a read is mapped to is taken into account here.
    fn to_bam_fields(&self, strand: Option<Direction>) -> (bam::record::CigarString, Vec<u8>, u16) {
        // Reconstruct the order of the remaining edit operations and condense CIGAR
        let mut num_matches: u32 = 0;
        let mut num_operations = 1;
        let mut edit_distance = 0;
        let mut last_edit_operation = None;
        let mut cigar = Vec::new();
        let mut md_tag = Vec::new();

        let track = match strand {
            Some(Direction::Forward) => Either::Left(self.0.iter()),
            Some(Direction::Backward) => Either::Right(self.0.iter().rev()),
            None => unreachable!(),
        };
        for edit_operation in track {
            edit_distance = Self::add_edit_distance(*edit_operation, edit_distance);

            num_matches = Self::add_md_edit_operation(
                Some(*edit_operation),
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
                            Self::add_cigar_edit_operation(lop, num_operations, &mut cigar);
                            num_operations = 1;
                            last_edit_operation = Some(*edit_operation);
                        }
                    },
                    EditOperation::Mismatch(_, _) => match lop {
                        EditOperation::Mismatch(_, _) | EditOperation::Match(_) => {
                            num_operations += 1;
                        }
                        _ => {
                            Self::add_cigar_edit_operation(lop, num_operations, &mut cigar);
                            num_operations = 1;
                            last_edit_operation = Some(*edit_operation);
                        }
                    },
                    EditOperation::Insertion(_) => match lop {
                        EditOperation::Insertion(_) => {
                            num_operations += 1;
                        }
                        _ => {
                            Self::add_cigar_edit_operation(lop, num_operations, &mut cigar);
                            num_operations = 1;
                            last_edit_operation = Some(*edit_operation);
                        }
                    },
                    EditOperation::Deletion(_, _) => match lop {
                        EditOperation::Deletion(_, _) => {
                            num_operations += 1;
                        }
                        _ => {
                            Self::add_cigar_edit_operation(lop, num_operations, &mut cigar);
                            num_operations = 1;
                            last_edit_operation = Some(*edit_operation);
                        }
                    },
                }
            } else {
                last_edit_operation = Some(*edit_operation);
            }
        }

        // Add remainder
        if let Some(lop) = last_edit_operation {
            Self::add_cigar_edit_operation(lop, num_operations, &mut cigar);
        }
        let _ = Self::add_md_edit_operation(None, None, strand, num_matches, &mut md_tag);

        (bam::record::CigarString(cigar), md_tag, edit_distance)
    }

    fn add_cigar_edit_operation(
        edit_operation: EditOperation,
        k: u32,
        cigar: &mut Vec<bam::record::Cigar>,
    ) {
        match edit_operation {
            EditOperation::Match(_) | EditOperation::Mismatch(_, _) => {
                cigar.push(bam::record::Cigar::Match(k))
            }
            EditOperation::Insertion(_) => cigar.push(bam::record::Cigar::Ins(k)),
            EditOperation::Deletion(_, _) => cigar.push(bam::record::Cigar::Del(k)),
        }
    }

    fn add_md_edit_operation(
        edit_operation: Option<EditOperation>,
        last_edit_operation: Option<EditOperation>,
        strand: Option<Direction>,
        mut k: u32,
        md_tag: &mut Vec<u8>,
    ) -> u32 {
        let comp_if_necessary = |reference_base| match strand {
            Some(Direction::Forward) => reference_base,
            Some(Direction::Backward) => dna::complement(reference_base),
            None => unreachable!(),
        };

        match edit_operation {
            Some(EditOperation::Match(_)) => k += 1,
            Some(EditOperation::Mismatch(_, reference_base)) => {
                let reference_base = comp_if_necessary(reference_base);
                md_tag.extend_from_slice(format!("{}{}", k, reference_base as char).as_bytes());
                k = 0;
            }
            Some(EditOperation::Insertion(_)) => {
                // Insertions are ignored in MD tags
            }
            Some(EditOperation::Deletion(_, reference_base)) => {
                let reference_base = comp_if_necessary(reference_base);
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

    fn add_edit_distance(edit_operation: EditOperation, distance: u16) -> u16 {
        if let EditOperation::Match(_) = edit_operation {
            distance
        } else {
            distance + 1
        }
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
    j: i16,
    current_interval: BiInterval,
    backward_index: i16,
    forward_index: i16,
    direction: Direction,
    gap_forwards: GapState,
    gap_backwards: GapState,
    alignment_score: f32,
    edit_node_id: NodeId,
    lookahead_score: f32,
}

impl PartialOrd for MismatchSearchStackFrame {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MismatchSearchStackFrame {
    fn cmp(&self, other: &Self) -> Ordering {
        self.lookahead_score
            .partial_cmp(&other.lookahead_score)
            .expect("This is not expected to fail")
    }
}

impl PartialEq for MismatchSearchStackFrame {
    fn eq(&self, other: &Self) -> bool {
        self.lookahead_score.eq(&other.lookahead_score)
    }
}

impl Eq for MismatchSearchStackFrame {}

/// For multi-identifier reference sequences like the human genome (that is split by chromosome)
/// this struct is used to keep a map of IDs and positions
#[derive(Serialize, Deserialize, Debug)]
pub struct FastaIdPosition {
    pub start: usize,
    pub end: usize,
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
    fn get_reference_identifier(&self, position: usize, pattern_length: usize) -> (i32, i64) {
        self.id_position
            .iter()
            .enumerate()
            .find(|(_, identifier)| {
                (identifier.start <= position) && (position + pattern_length - 1 <= identifier.end)
            })
            .map(|(index, identifier)| (index as i32, (position - identifier.start) as i64))
            .unwrap_or((-1, -1))
    }
}

/// A reversed FMD-index is used to compute the lower bound of mismatches of a read per position.
/// This allows for pruning the search tree. The values are minimal expected penalties towards
/// the respective ends of the query. Like alignment scores, these values are negative.
#[derive(Debug)]
struct BiDArray {
    d_backwards: SmallVec<[f32; 64]>,
    d_forwards: SmallVec<[f32; 64]>,
    split: usize,
}

impl BiDArray {
    pub(crate) fn new(
        pattern: &[u8],
        base_qualities: &[u8],
        split: usize,
        alignment_parameters: &AlignmentParameters,
        fmd_index: &FMDIndex<BWT, Less, Occ>,
    ) -> Self {
        Self {
            d_backwards: Self::compute_part(
                &pattern[..split],
                &base_qualities[..split],
                Direction::Forward,
                pattern.len(),
                alignment_parameters,
                fmd_index,
            ),
            d_forwards: Self::compute_part(
                &pattern[split..],
                &base_qualities[split..],
                Direction::Backward,
                pattern.len(),
                alignment_parameters,
                fmd_index,
            ),
            split,
        }
    }

    fn compute_part(
        pattern_part: &[u8],
        base_qualities_part: &[u8],
        direction: Direction,
        full_pattern_length: usize,
        alignment_parameters: &AlignmentParameters,
        fmd_index: &FMDIndex<BWT, Less, Occ>,
    ) -> SmallVec<[f32; 64]> {
        let directed_pattern_iterator = || match direction {
            Direction::Forward => Either::Left(pattern_part.iter()),
            Direction::Backward => Either::Right(pattern_part.iter().rev()),
        };
        let directed_index = |index, length| match direction {
            Direction::Forward => index,
            Direction::Backward => length - 1 - index,
        };

        directed_pattern_iterator()
            .enumerate()
            .scan(
                (0.0, None, fmd_index.init_interval()),
                |(z, last_mismatch_pos, interval), (index, &base)| {
                    *interval = match direction {
                        Direction::Forward => FMDIndex::forward_ext,
                        Direction::Backward => FMDIndex::backward_ext,
                    }(fmd_index, &interval, base);
                    if interval.size < 1 {
                        *z += directed_pattern_iterator()
                            .take(index + 1)
                            .skip(if let Some(lmp) = *last_mismatch_pos {
                                lmp + 1
                            } else {
                                0
                            })
                            .enumerate()
                            .map(|(j, &base_j)| {
                                let best_penalty_mm_only =
                                    alignment_parameters.difference_model.get_min_penalty(
                                        directed_index(j, full_pattern_length),
                                        full_pattern_length,
                                        base_j,
                                        base_qualities_part
                                            [directed_index(j, base_qualities_part.len())],
                                        true,
                                    );
                                let optimal_penalty =
                                    alignment_parameters.difference_model.get_min_penalty(
                                        directed_index(j, full_pattern_length),
                                        full_pattern_length,
                                        base_j,
                                        base_qualities_part
                                            [directed_index(j, base_qualities_part.len())],
                                        false,
                                    );
                                // The optimal penalty conditioned on position and base is subtracted because we
                                // optimize the ratio `AS/optimal_AS` (in log space) to find the best mappings
                                best_penalty_mm_only - optimal_penalty
                            })
                            .fold(std::f32::MIN, |acc, penalty| acc.max(penalty))
                            .max(alignment_parameters.penalty_gap_open)
                            .max(alignment_parameters.penalty_gap_extend);
                        *interval = fmd_index.init_interval();
                        *last_mismatch_pos = Some(index);
                    }
                    Some(*z)
                },
            )
            .collect()
    }

    fn get(&self, backward_index: i16, forward_index: i16) -> f32 {
        let inner = || {
            let forward_index = self
                .d_forwards
                .len()
                .checked_sub(1)?
                .checked_sub(forward_index as usize - self.split)?;
            Some(
                self.d_backwards.get(backward_index as usize)?
                    + self.d_forwards.get(forward_index)?,
            )
        };
        let out = inner().unwrap_or(0.0);
        debug_assert!(out <= 0.0);
        out
    }
}

struct ChunkIterator<'a, R>
where
    R: Read,
{
    chunk_size: usize,
    records: Peekable<bam::Records<'a, R>>,
}

impl<'a, R> Iterator for ChunkIterator<'a, R>
where
    R: Read,
{
    type Item = Vec<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let source_iterator_loan = &mut self.records;

        // If the underlying iterator is exhausted return None, too
        source_iterator_loan.peek()?;

        Some(
            source_iterator_loan
                .take(self.chunk_size)
                .map(|record| match record {
                    Ok(record) => Ok(record.into()),
                    Err(e) => Err(e),
                })
                .collect::<Result<Vec<_>, _>>()
                .expect("Input file is corrupt. Cancelling process."),
        )
    }
}

impl<'a, R> ChunkIterator<'a, R>
where
    R: Read,
{
    fn from_reader(records: bam::Records<'a, R>, chunk_size: usize) -> Self {
        Self {
            chunk_size,
            records: records.peekable(),
        }
    }
}

/// Loads index files and launches the mapping process
pub fn run(
    reads_path: &str,
    reference_path: &str,
    out_file_path: &str,
    alignment_parameters: &AlignmentParameters,
) -> Result<(), Box<dyn Error>> {
    debug!("Load FMD-index");
    let fmd_index = load_index_from_path(reference_path)?;

    debug!("Load suffix array");
    let suffix_array: RawSuffixArray = {
        let d_suffix_array =
            snap::read::FrameDecoder::new(File::open(format!("{}.tsa", reference_path))?);
        bincode::deserialize_from(d_suffix_array)?
    };

    debug!("Load position map");
    let identifier_position_map: FastaIdPositions = {
        let d_pi = snap::read::FrameDecoder::new(File::open(format!("{}.tpi", reference_path))?);
        bincode::deserialize_from(d_pi)?
    };

    let mut reads_reader = bam::Reader::from_path(reads_path)?;
    let _ = reads_reader.set_threads(4);
    let header = create_bam_header(&reads_reader, &identifier_position_map);
    let mut out_file = bam::Writer::from_path(out_file_path, &header, bam::Format::BAM)?;
    let _ = out_file.set_threads(4);

    thread_local! {
        static STACK_BUF: RefCell<MinMaxHeap<MismatchSearchStackFrame>> = RefCell::new(MinMaxHeap::with_capacity(STACK_LIMIT + 9));
        static TREE_BUF: RefCell<Tree<Option<EditOperation>>> = RefCell::new(Tree::with_capacity(STACK_LIMIT + 9));
    }

    debug!("Map reads");
    ChunkIterator::from_reader(reads_reader.records(), alignment_parameters.chunk_size)
        .map(|chunk| {
            trace!("Map chunk of reads in parallel");
            let results = chunk
                .par_iter()
                .map(|record| {
                    STACK_BUF.with(|stack_buf| {
                        TREE_BUF.with(|tree_buf| {
                            let start = Instant::now();
                            let hit_intervals = k_mismatch_search(
                                &record.sequence,
                                &record.base_qualities,
                                alignment_parameters,
                                &fmd_index,
                                &mut stack_buf.borrow_mut(),
                                &mut tree_buf.borrow_mut(),
                            );
                            let duration = start.elapsed();

                            (record, hit_intervals, duration)
                        })
                    })
                })
                .map_init(
                    rand::thread_rng,
                    |mut rng, (record, mut hit_interval, duration)| {
                        intervals_to_bam(
                            record,
                            &mut hit_interval,
                            &suffix_array,
                            &identifier_position_map,
                            Some(&duration),
                            &mut rng,
                        )
                    },
                )
                .flatten()
                .collect::<Vec<_>>();

            trace!("Write BAM records to output file serially");
            results
                .iter()
                .map(|record| out_file.write(record))
                .collect()
        })
        .collect::<Result<(), _>>()?;

    debug!("Done");
    Ok(())
}

pub fn create_bam_header(
    reads_reader: &bam::Reader,
    identifier_position_map: &FastaIdPositions,
) -> bam::Header {
    let mut header = bam::Header::from_template(reads_reader.header());
    for identifier_position in identifier_position_map.iter() {
        let mut header_record = bam::header::HeaderRecord::new(b"SQ");
        header_record.push_tag(b"SN", &identifier_position.identifier);
        header_record.push_tag(
            b"LN",
            &(identifier_position.end - identifier_position.start + 1),
        );
        header.push_record(&header_record);
    }
    {
        let mut header_record = bam::header::HeaderRecord::new(b"PG");
        header_record.push_tag(b"ID", &crate_name!());
        header_record.push_tag(b"PN", &CRATE_NAME);
        header_record.push_tag(b"VN", &crate_version!());
        let cmdline = {
            let mut tmp = std::env::args().fold(String::new(), |acc, part| acc + &part + " ");
            let _ = tmp.pop();
            tmp
        };
        header_record.push_tag(b"CL", &cmdline);
        header.push_record(&header_record);
    }
    header
}

/// Convert suffix array intervals to positions and BAM records
pub fn intervals_to_bam<R, S>(
    input_record: &Record,
    intervals: &mut BinaryHeap<HitInterval>,
    suffix_array: &S,
    identifier_position_map: &FastaIdPositions,
    duration: Option<&Duration>,
    rng: &mut R,
) -> Vec<bam::Record>
where
    R: RngCore,
    S: SuffixArray,
{
    let record = if let Some(best_alignment) = intervals.pop() {
        let mapping_quality = estimate_mapping_quality(&best_alignment, &intervals);
        best_alignment
            .interval
            .forward()
            .occ(suffix_array)
            .iter()
            .filter(|&&position| position < (suffix_array.len() - 2) / 2)
            .map(|&position| (position, Direction::Forward))
            .chain(
                best_alignment
                    .interval
                    .revcomp()
                    .occ(suffix_array)
                    .iter()
                    .filter(|&&position| position < (suffix_array.len() - 2) / 2)
                    .map(|&position| (position, Direction::Backward)),
            )
            .choose(rng)
            .map(|(position, strand)| {
                let (tid, position) = identifier_position_map.get_reference_identifier(
                    position,
                    best_alignment.edit_operations.effective_len(),
                );
                bam_record_helper(
                    input_record,
                    position,
                    Some(&best_alignment),
                    Some(mapping_quality),
                    tid,
                    Some(strand),
                    duration,
                )
            })
            .unwrap()
    } else {
        // No match found, report unmapped read
        bam_record_helper(input_record, -1, None, None, -1, None, duration)
    };
    vec![record]
}

/// Computes optimal per-base alignment scores for a read,
/// conditioned on base qualitites and scoring model.
/// This function panics if `pattern` and `base_qualities` are not of equal length.
pub fn compute_optimal_scores(
    pattern: &[u8],
    base_qualities: &[u8],
    difference_model: &SequenceDifferenceModelDispatch,
) -> SmallVec<[f32; 128]> {
    assert_eq!(pattern.len(), base_qualities.len());
    pattern
        .iter()
        .zip(base_qualities)
        .enumerate()
        .map(|(i, (&base, &quality))| {
            difference_model.get_min_penalty(i, pattern.len(), base, quality, false)
        })
        .collect::<SmallVec<_>>()
}

/// Estimate mapping quality based on the number of hits for a particular read, its alignment score,
/// and its base qualities
fn estimate_mapping_quality(
    best_alignment: &HitInterval,
    other_alignments: &BinaryHeap<HitInterval>,
) -> u8 {
    const MAX_MAPQ: f32 = 37.0;

    let alignment_probability = {
        let ratio_best = 2_f32.powf(best_alignment.alignment_score);
        if best_alignment.interval.size > 1 {
            // Multi-mapping
            1.0 / best_alignment.interval.size as f32
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
    .max(0.0)
    .min(1.0);

    // Produce Phred score
    (-10.0 * (1.0 - alignment_probability).log10())
        .min(MAX_MAPQ)
        .round() as u8
}

/// Create and return a BAM record of either a hit or an unmapped read
fn bam_record_helper(
    input_record: &Record,
    position: i64,
    hit_interval: Option<&HitInterval>,
    mapq: Option<u8>,
    tid: i32,
    strand: Option<Direction>,
    duration: Option<&Duration>,
) -> bam::Record {
    let mut bam_record = bam::record::Record::new();

    let (cigar, md_tag, edit_distance) = if let Some(hit_interval) = hit_interval {
        let (cigar, md_tag, edit_distance) = hit_interval.edit_operations.to_bam_fields(strand);
        (Some(cigar), Some(md_tag), Some(edit_distance))
    } else {
        (None, None, None)
    };

    // Set mandatory properties of the BAM record
    match strand {
        Some(Direction::Forward) | None => {
            bam_record.set(
                &input_record.name,
                cigar.as_ref(),
                &input_record.sequence,
                &input_record.base_qualities,
            );
        }
        Some(Direction::Backward) => {
            bam_record.set(
                &input_record.name,
                // CIGAR strings and MD tags are reversed during generation
                cigar.as_ref(),
                &dna::revcomp(&input_record.sequence),
                &input_record
                    .base_qualities
                    .iter()
                    .rev()
                    .copied()
                    .collect::<Vec<_>>(),
            );
        }
    }
    bam_record.set_tid(tid);
    bam_record.set_pos(position);
    if let Some(mapq) = mapq {
        bam_record.set_mapq(mapq);
    }

    // Flag read that maps to reverse strand
    if let Some(Direction::Backward) = strand {
        bam_record.set_reverse();
    }

    // Flag unmapped read
    if position == -1 {
        bam_record.set_unmapped();
    }

    // Position of mate (-1 = *)
    bam_record.set_mpos(-1);

    // Reference sequence of mate (-1 = *)
    bam_record.set_mtid(-1);

    // Add optional tags

    // Add tags that were already present in the input
    for (input_tag, value) in &input_record.bam_tags {
        bam_record.push_aux(input_tag, &value.borrow_htslib_bam_record());
    }

    if let Some(hit_interval) = hit_interval {
        bam_record.push_aux(
            b"AS",
            &bam::record::Aux::Integer(hit_interval.alignment_score.round() as i64),
        );
    };
    if let Some(edit_distance) = edit_distance {
        bam_record.push_aux(b"NM", &bam::record::Aux::Integer(i64::from(edit_distance)));
    };
    // CIGAR strings and MD tags are reversed during generation
    if let Some(md_tag) = md_tag {
        bam_record.push_aux(b"MD", &bam::record::Aux::String(&md_tag));
    }
    if let Some(duration) = duration {
        // Add the time that was needed for mapping the read
        bam_record.push_aux(
            b"XD",
            &bam::record::Aux::Integer(duration.as_micros() as i64),
        );
    }

    bam_record
}

/// Derive Cigar string from oddly-ordered tracks of edit operations.
/// Since we start aligning at the center of a read, tracks of edit operations
/// are not ordered by position in the read. Also, the track of edit operations
/// must be extracted from the (possibly huge) tree which is built during
/// backtracking for size reasons.
fn extract_edit_operations(
    end_node: NodeId,
    edit_tree: &Tree<Option<EditOperation>>,
    pattern_len: usize,
) -> EditOperationsTrack {
    // Restore outer ordering of the edit operation by the positions they carry as values.
    // Whenever there are deletions in the pattern, there is no simple rule to reconstruct the ordering.
    // So, edit operations carrying the same position are pushed onto the same bucket and dealt with later.
    let mut cigar_order_outer: BTreeMap<u16, SmallVec<[EditOperation; 8]>> = BTreeMap::new();

    edit_tree
        .ancestors(end_node)
        .filter_map(|&edit_operation| edit_operation)
        .for_each(|edit_operation| {
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
fn check_and_push_stack_frame(
    mut stack_frame: MismatchSearchStackFrame,
    pattern: &[u8],
    alignment_parameters: &AlignmentParameters,
    lower_bound: f32,
    edit_operation: EditOperation,
    edit_tree: &mut Tree<Option<EditOperation>>,
    stack: &mut MinMaxHeap<MismatchSearchStackFrame>,
    intervals: &mut BinaryHeap<HitInterval>,
) {
    // Too many mismatches
    if alignment_parameters
        .mismatch_bound
        .reject(stack_frame.lookahead_score, pattern.len())
    {
        return;
    }

    if let Some(best_scoring_interval) = intervals.peek() {
        if alignment_parameters.mismatch_bound.reject_iterative(
            stack_frame.lookahead_score,
            best_scoring_interval.alignment_score,
        ) {
            return;
        }
    }

    stack_frame.edit_node_id = edit_tree.add_node(Some(edit_operation), stack_frame.edit_node_id);

    // This route through the read graph is finished successfully, push the interval
    if stack_frame.j < 0 || stack_frame.j > (pattern.len() as i16 - 1) {
        let edit_operations =
            extract_edit_operations(stack_frame.edit_node_id, edit_tree, pattern.len());
        intervals.push(HitInterval {
            interval: stack_frame.current_interval,
            alignment_score: stack_frame.alignment_score,
            edit_operations,
        });
        print_debug(&stack_frame, intervals, lower_bound); // FIXME
        return;
    }

    stack.push(stack_frame);
}

/// FIXME
fn print_debug(
    stack_frame: &MismatchSearchStackFrame,
    intervals: &BinaryHeap<HitInterval>,
    lower_bound: f32,
) {
    let switch = false; // TODO: Switch me on/off!

    if switch {
        let best_as = match intervals.peek() {
            Some(v) => v.alignment_score,
            None => 0.0,
        };

        eprintln!(
            "{}\t{}\t{}\t{}\t{}",
            stack_frame.alignment_score,
            stack_frame.alignment_score + lower_bound,
            stack_frame.backward_index,
            stack_frame.forward_index,
            best_as,
        );
    }
}

/// Finds all suffix array intervals for the current pattern
/// w.r.t. supplied alignment parameters
pub fn k_mismatch_search(
    pattern: &[u8],
    base_qualities: &[u8],
    parameters: &AlignmentParameters,
    fmd_index: &FMDIndex<BWT, Less, Occ>,
    stack: &mut MinMaxHeap<MismatchSearchStackFrame>,
    edit_tree: &mut Tree<Option<EditOperation>>,
) -> BinaryHeap<HitInterval> {
    let center_of_read = pattern.len() / 2;
    let bi_d_array = BiDArray::new(
        pattern,
        base_qualities,
        center_of_read,
        parameters,
        fmd_index,
    );

    let optimal_penalties =
        compute_optimal_scores(pattern, base_qualities, &parameters.difference_model);
    let mut hit_intervals = BinaryHeap::new();

    let mut stack_size_limit_reported = false;
    stack.clear();
    edit_tree.clear();

    stack.push(MismatchSearchStackFrame {
        j: center_of_read as i16,
        current_interval: fmd_index.init_interval(),
        backward_index: center_of_read as i16 - 1,
        forward_index: center_of_read as i16,
        direction: Direction::Forward,
        gap_backwards: GapState::Closed,
        gap_forwards: GapState::Closed,
        alignment_score: 0.0,
        edit_node_id: edit_tree.add_node(None, None),
        lookahead_score: 0.0,
    });

    while let Some(stack_frame) = stack.pop_max() {
        // Determine direction of progress for next iteration on this stack frame
        let next_j;
        let next_backward_index;
        let next_forward_index;
        let fmd_ext_interval = match stack_frame.direction {
            Direction::Forward => {
                next_forward_index = stack_frame.forward_index + 1;
                next_backward_index = stack_frame.backward_index;
                next_j = stack_frame.backward_index;
                stack_frame.current_interval.swapped()
            }
            Direction::Backward => {
                next_forward_index = stack_frame.forward_index;
                next_backward_index = stack_frame.backward_index - 1;
                next_j = stack_frame.forward_index;
                stack_frame.current_interval
            }
        };

        let optimal_penalty = optimal_penalties[stack_frame.j as usize];

        // Calculate the lower bounds for extension
        let lower_bound = bi_d_array.get(next_backward_index, next_forward_index);

        print_debug(&stack_frame, &hit_intervals, lower_bound); // FIXME

        // Since we operate on a priority stack, we can assume that there are no
        // better scoring frames on the stack, so we are going to stop the search.
        if let Some(best_scoring_interval) = hit_intervals.peek() {
            if parameters.mismatch_bound.reject_iterative(
                stack_frame.lookahead_score,
                best_scoring_interval.alignment_score,
            ) {
                break;
            }
        }

        //
        // Insertion in read / deletion in reference
        //
        let new_alignment_score = if (stack_frame.gap_backwards == GapState::Insertion
            && stack_frame.direction.is_forward())
            || (stack_frame.gap_forwards == GapState::Insertion
                && stack_frame.direction.is_backward())
        {
            parameters.penalty_gap_extend
        } else {
            parameters.penalty_gap_open
        } - optimal_penalty
            + stack_frame.alignment_score;

        check_and_push_stack_frame(
            MismatchSearchStackFrame {
                j: next_j,
                backward_index: next_backward_index,
                forward_index: next_forward_index,
                direction: stack_frame.direction.reverse(),
                // Mark opened gap at the corresponding end
                gap_backwards: if stack_frame.direction.is_backward() {
                    stack_frame.gap_backwards
                } else {
                    GapState::Insertion
                },
                gap_forwards: if stack_frame.direction.is_forward() {
                    stack_frame.gap_forwards
                } else {
                    GapState::Insertion
                },
                alignment_score: new_alignment_score,
                lookahead_score: new_alignment_score + lower_bound,
                ..stack_frame
            },
            pattern,
            &parameters,
            lower_bound,
            EditOperation::Insertion(stack_frame.j as u16),
            edit_tree,
            stack,
            &mut hit_intervals,
        );

        // Bidirectional extension of the (hit) interval
        let mut s = 0;
        let mut o;
        let mut l = fmd_ext_interval.lower_rev;
        for &c in b"$TGCA".iter() {
            let mut interval_prime = {
                l += s;
                o = if fmd_ext_interval.lower == 0 {
                    0
                } else {
                    fmd_index.occ(fmd_ext_interval.lower - 1, c)
                };

                // Interval size I^s
                s = fmd_index.occ(fmd_ext_interval.lower + fmd_ext_interval.size - 1, c) - o;

                // No need to branch for technical characters and zero-sized intervals
                if c == b'$' || s < 1 {
                    continue;
                }

                BiInterval {
                    lower: fmd_index.less(c) + o,
                    lower_rev: l,
                    size: s,
                    match_size: fmd_ext_interval.match_size + 1,
                }
            };

            // Special treatment of forward extension
            let c = match stack_frame.direction {
                Direction::Forward => {
                    interval_prime = interval_prime.swapped();
                    dna::complement(c)
                }
                Direction::Backward => c,
            };

            //
            // Deletion in read / insertion in reference
            //
            let new_alignment_score = if (stack_frame.gap_forwards == GapState::Deletion
                && stack_frame.direction.is_forward())
                || (stack_frame.gap_backwards == GapState::Deletion
                    && stack_frame.direction.is_backward())
            {
                parameters.penalty_gap_extend
            } else {
                parameters.penalty_gap_open
            } + stack_frame.alignment_score;

            check_and_push_stack_frame(
                MismatchSearchStackFrame {
                    current_interval: interval_prime,
                    // Mark open gap at the corresponding end
                    gap_backwards: if stack_frame.direction.is_backward() {
                        GapState::Deletion
                    } else {
                        stack_frame.gap_backwards
                    },
                    gap_forwards: if stack_frame.direction.is_forward() {
                        GapState::Deletion
                    } else {
                        stack_frame.gap_forwards
                    },
                    alignment_score: new_alignment_score,
                    lookahead_score: new_alignment_score + lower_bound,
                    ..stack_frame
                },
                pattern,
                &parameters,
                lower_bound,
                EditOperation::Deletion(stack_frame.j as u16, c),
                edit_tree,
                stack,
                &mut hit_intervals,
            );

            //
            // Match/mismatch
            //
            let new_alignment_score = parameters.difference_model.get(
                stack_frame.j as usize,
                pattern.len(),
                c,
                pattern[stack_frame.j as usize],
                base_qualities[stack_frame.j as usize],
            ) - optimal_penalty
                + stack_frame.alignment_score;

            check_and_push_stack_frame(
                MismatchSearchStackFrame {
                    j: next_j,
                    current_interval: interval_prime,
                    backward_index: next_backward_index,
                    forward_index: next_forward_index,
                    direction: stack_frame.direction.reverse(),
                    // Mark closed gap at the corresponding end
                    gap_backwards: if stack_frame.direction.is_backward() {
                        stack_frame.gap_backwards
                    } else {
                        GapState::Closed
                    },
                    gap_forwards: if stack_frame.direction.is_forward() {
                        stack_frame.gap_forwards
                    } else {
                        GapState::Closed
                    },
                    alignment_score: new_alignment_score,
                    lookahead_score: new_alignment_score + lower_bound,
                    ..stack_frame
                },
                pattern,
                &parameters,
                lower_bound,
                if c == pattern[stack_frame.j as usize] {
                    EditOperation::Match(stack_frame.j as u16)
                } else {
                    EditOperation::Mismatch(stack_frame.j as u16, c)
                },
                edit_tree,
                stack,
                &mut hit_intervals,
            );
        }

        // Only search until we've found a multi-hit
        if (hit_intervals.len() > 9)
            || (hit_intervals.len() == 1 && hit_intervals.peek().unwrap().interval.size > 1)
        {
            return hit_intervals;
        }

        // Limit stack size
        if edit_tree.len() >= STACK_LIMIT {
            if !stack_size_limit_reported {
                trace!(
                    "Stack size limit reached (read length: {} bp). Remove poor partial alignments from stack (size: {}).",
                    pattern.len(),
                    stack.len(),
                );
                stack_size_limit_reported = true;
            }

            for _ in 0..(edit_tree.len() - STACK_LIMIT + 9) {
                if let Some(poor_frame) = stack.pop_min() {
                    edit_tree.remove(poor_frame.edit_node_id);
                }
            }
        }
    }
    hit_intervals
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{mismatch_bounds::*, sequence_difference_models::*};
    use assert_approx_eq::assert_approx_eq;
    use bio::{
        alphabets,
        data_structures::{
            bwt::{bwt, less},
            fmindex::FMIndex,
            suffix_array::suffix_array,
        },
    };
    use rust_htslib::bam;

    fn build_auxiliary_structures(
        reference: &mut Vec<u8>,
        alphabet: &alphabets::Alphabet,
    ) -> FMDIndex<BWT, Less, Occ> {
        let ref_seq_rev_compl = alphabets::dna::revcomp(reference.iter());
        reference.extend_from_slice(b"$");
        reference.extend_from_slice(&ref_seq_rev_compl);
        drop(ref_seq_rev_compl);
        reference.extend_from_slice(b"$");

        let sa = suffix_array(&reference);
        let bwt = bwt(&reference, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);

        FMDIndex::from(FMIndex::new(bwt, less, occ))
    }

    #[test]
    fn test_exact_search() {
        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);
        let suffix_array = suffix_array(&ref_seq);

        // Test occurring pattern
        let pattern_1 = b"TA";
        let sai_1 = fmd_index.backward_search(pattern_1.iter());
        let positions_1 = sai_1.occ(&suffix_array);
        assert_eq!(positions_1, vec![10, 3]);

        // Test non-occurring pattern
        let pattern_2 = b"GG";
        let sai_2 = fmd_index.backward_search(pattern_2.iter());
        let positions_2 = sai_2.occ(&suffix_array);
        assert_eq!(positions_2, Vec::<usize>::new());
    }

    #[test]
    fn test_inexact_search() {
        let difference_model = SequenceDifferenceModelDispatch::from(TestDifferenceModel {
            deam_score: 0.5,
            mm_score: -1.0,
            match_score: 0.0,
        });

        let parameters = AlignmentParameters {
            difference_model,
            mismatch_bound: TestBound { threshold: -1.0 }.into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
        };

        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "ACGTACGTACGTACGT".as_bytes().to_owned();

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);
        let suffix_array = suffix_array(&ref_seq);

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
        let difference_model = SequenceDifferenceModelDispatch::from(TestDifferenceModel {
            deam_score: -10.0,
            mm_score: -10.0,
            match_score: 1.0,
        });

        let parameters = AlignmentParameters {
            difference_model,
            mismatch_bound: TestBound { threshold: -1.0 }.into(),
            penalty_gap_open: -20.0,
            penalty_gap_extend: -10.0,
            chunk_size: 1,
        };

        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "GAAAAG".as_bytes().to_owned();

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);
        let suffix_array = suffix_array(&ref_seq);

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
        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "GATTACA".as_bytes().to_owned(); // revcomp = TGTAATC

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let difference_model = SequenceDifferenceModelDispatch::from(SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.3,
                three_prime_overhang: 0.3,
            },
            0.001,
            0.8,
            0.02,
        ));

        let representative_mismatch_penalty =
            difference_model.get_representative_mismatch_penalty();

        let parameters = AlignmentParameters {
            difference_model,
            mismatch_bound: TestBound { threshold: 0.0 }.into(),
            penalty_gap_open: 0.00001_f32.log2(),
            penalty_gap_extend: representative_mismatch_penalty,
            chunk_size: 1,
        };

        let pattern = b"CCCCCCC";
        let base_qualities = [9, 100, 100, 100, 100, 9, 100];

        let center_of_read = pattern.len() / 2;

        let bi_d_array = BiDArray::new(
            pattern,
            &base_qualities,
            center_of_read,
            &parameters,
            &fmd_index,
        );

        assert_eq!(&*bi_d_array.d_backwards, &[0.0, -1.7041414, -1.7041414]);
        assert_eq!(
            &*bi_d_array.d_forwards,
            &[0.0, -1.9092996, -1.9092996, -3.8185992]
        );

        assert_eq!(
            bi_d_array.get(1, 4),
            bi_d_array.d_backwards[1] + bi_d_array.d_forwards[2]
        );
        assert_eq!(
            bi_d_array.get(2, 3),
            bi_d_array.d_backwards[2] + bi_d_array.d_forwards[3]
        );
        assert_eq!(
            bi_d_array.get(0, 6),
            bi_d_array.d_backwards[0] + bi_d_array.d_forwards[0]
        );
    }

    #[test]
    fn test_gapped_alignment() {
        let difference_model = SequenceDifferenceModelDispatch::from(TestDifferenceModel {
            deam_score: -10.0,
            mm_score: -10.0,
            match_score: 0.0,
        });

        let parameters = AlignmentParameters {
            difference_model,
            mismatch_bound: TestBound { threshold: -2.0 }.into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
        };

        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "TAT".as_bytes().to_owned(); // revcomp = "ATA"

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);
        let suffix_array = suffix_array(&ref_seq);

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
    fn test_vindija_pwm_alignment() {
        let difference_model = SequenceDifferenceModelDispatch::from(VindijaPWM::new());

        let parameters = AlignmentParameters {
            difference_model,
            // Disable gaps
            mismatch_bound: TestBound { threshold: -30.0 }.into(),
            penalty_gap_open: -200.0,
            penalty_gap_extend: -100.0,
            chunk_size: 1,
        };

        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "CCCCCC".as_bytes().to_owned(); // revcomp = "ATA"

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);
        let sar = suffix_array(&ref_seq);

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

        let mut ref_seq = "AAAAAA".as_bytes().to_owned(); // revcomp = "ATA"

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

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
        );

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_approx_eq!(alignment_score[0], -10.965062);
    }

    #[test]
    fn test_ord_impl() {
        let mut edit_tree = Tree::new();
        let root_id = edit_tree.add_node(1, None);

        let map_params_large = MismatchSearchStackFrame {
            alignment_score: -5.0,
            lookahead_score: -5.0,
            j: 0,
            current_interval: BiInterval {
                lower: 5,
                lower_rev: 5,
                match_size: 5,
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
            lookahead_score: -20.0,
            j: 0,
            current_interval: BiInterval {
                lower: 5,
                lower_rev: 5,
                match_size: 5,
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
        let difference_model = SequenceDifferenceModelDispatch::from(VindijaPWM::new());
        let repr_mm_penalty = difference_model.get_representative_mismatch_penalty();

        let parameters = AlignmentParameters {
            penalty_gap_open: 3.5 * difference_model.get_representative_mismatch_penalty(),
            penalty_gap_extend: 1.5 * difference_model.get_representative_mismatch_penalty(),
            difference_model,
            chunk_size: 1,
            mismatch_bound: Discrete::new(0.01, 0.02, repr_mm_penalty).into(),
        };

        let alphabet = alphabets::dna::alphabet();

        // "correct" "AAAAAAAAAAAAAAAAAAAA" (20x 'A') "incorrect"
        let mut ref_seq = "GTTGTATTTTTAGTAGAGACAGGGTTTCATCATGTTGGCCAGAAAAAAAAAAAAAAAAAAAATTTGTATTTTTAGTAGAGACAGGCTTTCATCATGTTGGCCAG"
            .as_bytes()
            .to_owned();

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);
        let sar = suffix_array(&ref_seq);

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
        );

        let alignment_scores = intervals
            .iter()
            .map(|f| f.alignment_score)
            .collect::<Vec<_>>();
        assert_eq!(alignment_scores, vec![-10.936638, -38.379524, -10.965062]);

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

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: TestBound { threshold: -3.0 }.into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
        };
        let alphabet = alphabets::dna::alphabet();

        //
        // Deletion
        //
        let mut ref_seq = "GATTAGCA".as_bytes().to_owned();

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

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
        );
        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit
            .edit_operations
            .to_bam_fields(Some(Direction::Forward));

        assert_eq!(
            cigar,
            bam::record::CigarString(vec![
                bam::record::Cigar::Match(4),
                bam::record::Cigar::Del(1),
                bam::record::Cigar::Match(2)
            ])
        );

        //
        // 2-base deletion
        //
        let mut ref_seq = "GATTACAG".as_bytes().to_owned(); // CTGTAATC

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

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
        );

        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit
            .edit_operations
            .to_bam_fields(Some(Direction::Forward));

        assert_eq!(best_hit.alignment_score, -3.0);
        assert_eq!(
            cigar,
            bam::record::CigarString(vec![
                bam::record::Cigar::Match(3),
                bam::record::Cigar::Del(2),
                bam::record::Cigar::Match(3)
            ])
        );

        //
        // Insertion
        //
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

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
        );
        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit
            .edit_operations
            .to_bam_fields(Some(Direction::Forward));

        assert_eq!(best_hit.alignment_score, -2.0);
        assert_eq!(
            cigar,
            bam::record::CigarString(vec![
                bam::record::Cigar::Match(5),
                bam::record::Cigar::Ins(1),
                bam::record::Cigar::Match(2)
            ])
        );

        //
        // 2-base insertion
        //
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

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
        );
        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit
            .edit_operations
            .to_bam_fields(Some(Direction::Forward));

        assert_eq!(best_hit.alignment_score, -3.0);
        assert_eq!(
            cigar,
            bam::record::CigarString(vec![
                bam::record::Cigar::Match(5),
                bam::record::Cigar::Ins(2),
                bam::record::Cigar::Match(2)
            ])
        );

        //
        // 3-base insertion
        //
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let pattern = "GATTAGTGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let parameters = AlignmentParameters {
            difference_model: difference_model.into(),
            mismatch_bound: TestBound { threshold: -4.0 }.into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
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
        );
        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit
            .edit_operations
            .to_bam_fields(Some(Direction::Forward));

        assert_eq!(best_hit.alignment_score, -4.0);
        assert_eq!(
            cigar,
            bam::record::CigarString(vec![
                bam::record::Cigar::Match(5),
                bam::record::Cigar::Ins(3),
                bam::record::Cigar::Match(2)
            ])
        );
    }

    #[test]
    fn test_md_tag() {
        let difference_model = TestDifferenceModel {
            deam_score: -1.0,
            mm_score: -2.0,
            match_score: 0.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: TestBound { threshold: -1.0 }.into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
        };
        let alphabet = alphabets::dna::alphabet();

        //
        // Mutation
        //
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let pattern = "GATTATA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

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
        );
        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit
            .edit_operations
            .to_bam_fields(Some(Direction::Forward));

        assert_eq!(md_tag, "5C1".as_bytes());

        //
        // Deletion
        //
        let mut ref_seq = "GATTAGCA".as_bytes().to_owned();

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let pattern = "ATTACA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let parameters = AlignmentParameters {
            difference_model: difference_model.into(),
            mismatch_bound: TestBound { threshold: -3.0 }.into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
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
        );
        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit
            .edit_operations
            .to_bam_fields(Some(Direction::Forward));

        assert_eq!(md_tag, "4^G2".as_bytes());

        //
        // 2-base deletion
        //
        let mut ref_seq = "GATTACAG".as_bytes().to_owned(); // CTGTAATC

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

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
        );

        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit
            .edit_operations
            .to_bam_fields(Some(Direction::Forward));

        assert_eq!(md_tag, "3^TA3".as_bytes());

        //
        // Insertion
        //
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

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
        );
        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit
            .edit_operations
            .to_bam_fields(Some(Direction::Forward));

        assert_eq!(md_tag, "7".as_bytes());

        //
        // 2-base insertion
        //
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

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
        );
        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit
            .edit_operations
            .to_bam_fields(Some(Direction::Forward));

        assert_eq!(md_tag, "7".as_bytes());
    }

    #[test]
    fn test_reverse_strand_search_2() {
        let difference_model = SequenceDifferenceModelDispatch::from(TestDifferenceModel {
            deam_score: -1.0,
            mm_score: -1.0,
            match_score: 0.0,
        });

        let parameters = AlignmentParameters {
            difference_model,
            mismatch_bound: TestBound { threshold: 0.0 }.into(),
            penalty_gap_open: -3.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
        };

        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "AAAGCGTTTGCG".as_bytes().to_owned();

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);
        let suffix_array = suffix_array(&ref_seq);

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
        let difference_model = SequenceDifferenceModelDispatch::from(TestDifferenceModel {
            deam_score: -1.0,
            mm_score: -1.0,
            match_score: 0.0,
        });

        let parameters = AlignmentParameters {
            difference_model,
            mismatch_bound: TestBound { threshold: -1.0 }.into(),
            penalty_gap_open: -3.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
        };

        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "GATTACA".as_bytes().to_owned(); // revcomp = "TGTAATC"

        // Reference
        let fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);
        let suffix_array = suffix_array(&ref_seq);

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
            .to_bam_fields(Some(Direction::Backward));
        assert_eq!(md, b"1T2");

        assert_eq!(edop, 1);
    }
}
