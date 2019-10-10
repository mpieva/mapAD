use std::{
    cmp::Ordering,
    collections::{binary_heap::BinaryHeap, BTreeMap},
    error::Error,
    fs::File,
    iter::Peekable,
};

use clap::{crate_name, crate_version};
use ego_tree::{NodeId, Tree};
use either::Either;
use log::{debug, trace};
use rayon::prelude::*;
use smallvec::SmallVec;

use bio::{
    alphabets::dna,
    data_structures::{
        bwt::Occ,
        fmindex::{BiInterval, FMDIndex, FMIndex, FMIndexable},
    },
};

use bincode;
use rust_htslib::{bam, bam::Read};
use serde::{Deserialize, Serialize};
use snap;

use crate::{
    sequence_difference_models::SequenceDifferenceModel,
    utils::{AlignmentParameters, AllowedMismatches},
};

pub const CRATE_NAME: &str = "mapAD";

/// Helper struct to bundle index files
struct UnderlyingDataFMDIndex {
    bwt: Vec<u8>,
    less: Vec<usize>,
    occ: Occ,
}

impl UnderlyingDataFMDIndex {
    fn load(path: &str) -> Result<UnderlyingDataFMDIndex, bincode::Error> {
        debug!("Load BWT");
        let bwt: Vec<u8> = {
            let d_bwt = snap::Reader::new(File::open(format!("{}.tbw", path))?);
            bincode::deserialize_from(d_bwt)?
        };

        debug!("Load \"C\" table");
        let less: Vec<usize> = {
            let d_less = snap::Reader::new(File::open(format!("{}.tle", path))?);
            bincode::deserialize_from(d_less)?
        };

        debug!("Load \"Occ\" table");
        let occ: Occ = {
            let d_occ = snap::Reader::new(File::open(format!("{}.toc", path))?);
            bincode::deserialize_from(d_occ)?
        };

        Ok(UnderlyingDataFMDIndex { bwt, less, occ })
    }
}

/// A subset of MismatchSearchStackFrame to store hits
#[derive(Debug)]
pub struct HitInterval {
    interval: BiInterval,
    alignment_score: f32,
    cigar: bam::record::CigarString,
    md_tag: Vec<u8>,
}

impl PartialOrd for HitInterval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HitInterval {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.alignment_score > other.alignment_score {
            Ordering::Greater
        } else if self.alignment_score < other.alignment_score {
            Ordering::Less
        } else {
            Ordering::Equal
        }
    }
}

impl PartialEq for HitInterval {
    fn eq(&self, other: &Self) -> bool {
        (self.alignment_score - other.alignment_score).abs() < std::f32::EPSILON
    }
}

impl Eq for HitInterval {}

/// Simple zero-cost direction enum to increase readability
#[derive(Debug, Copy, Clone)]
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

#[derive(Debug, Copy, Clone, PartialEq)]
enum EditOperation {
    Insertion(usize),
    Deletion(usize, u8),
    Match(usize),
    Mismatch(usize, u8),
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
struct MismatchSearchStackFrame {
    j: i16,
    current_interval: BiInterval,
    backward_index: i16,
    forward_index: i16,
    direction: Direction,
    gap_forwards: GapState,
    gap_backwards: GapState,
    alignment_score: f32,
    priority: f32,
    edit_node_id: NodeId,
}

impl PartialOrd for MismatchSearchStackFrame {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MismatchSearchStackFrame {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.priority > other.priority {
            Ordering::Greater
        } else if self.priority < other.priority {
            Ordering::Less
        } else {
            Ordering::Equal
        }
    }
}

impl PartialEq for MismatchSearchStackFrame {
    fn eq(&self, other: &Self) -> bool {
        (self.priority - other.priority).abs() < std::f32::EPSILON
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
    fn get_reference_identifier(&self, position: usize) -> (i32, i32) {
        self.iter()
            .enumerate()
            .find(|(_, identifier)| (identifier.start <= position) && (position <= identifier.end))
            .map(|(index, identifier)| (index as i32, (position - identifier.start) as i32 + 1))
            .unwrap_or((-1, -1))
    }
}

/// A reversed FMD-index is used to compute the lower bound of mismatches of a read per position.
/// This allows for pruning the search tree. The values are minimal expected penalties towards
/// the respective ends of the query. In contrast to alignment scores, these values are _positive_.
#[derive(Debug)]
struct DArray {
    d_array: SmallVec<[f32; 64]>,
}

impl DArray {
    fn new<T: SequenceDifferenceModel + Sync>(
        pattern: &[u8],
        base_qualities: &[u8],
        pattern_length: usize,
        read_length: usize,
        direction: Direction,
        alignment_parameters: &AlignmentParameters<T>,
        fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
    ) -> Self {
        let extend_interval = match direction {
            Direction::Forward => FMDIndex::forward_ext,
            Direction::Backward => FMDIndex::backward_ext,
        };

        let directed_pattern_index = |i| match direction {
            Direction::Forward => i,
            Direction::Backward => read_length - pattern_length + i,
        };

        let mut z = 0.0;
        let mut interval = fmd_index.init_interval();
        let pattern = match direction {
            Direction::Forward => Either::Left(pattern.iter()),
            Direction::Backward => Either::Right(pattern.iter().rev()),
        };

        let pattern = pattern.enumerate().map(|(index, &base)| {
            interval = extend_interval(fmd_index, &interval, base);
            if interval.size < 1 {
                z += alignment_parameters
                    .difference_model
                    .get_min_penalty(
                        directed_pattern_index(index),
                        read_length,
                        base,
                        base_qualities[index],
                        true,
                    )
                    .max(alignment_parameters.penalty_gap_open)
                    .max(alignment_parameters.penalty_gap_extend);
                interval = fmd_index.init_interval();
            }
            z
        });

        // FIXME: It's not optimal to call collect() twice. Side-effects seem to make this necessary.
        DArray {
            d_array: match direction {
                Direction::Forward => pattern.collect(),
                Direction::Backward => pattern
                    .collect::<SmallVec<[f32; 64]>>()
                    .into_iter()
                    .rev()
                    .collect(),
            },
        }
    }

    #[inline]
    fn get(&self, index: i16) -> f32 {
        *self.d_array.get(index as usize).unwrap_or(&0.0)
    }
}

#[derive(Copy, Clone)]
struct Penalties {
    max_allowed_penalties: f32,
    lower_bound: f32,
    representative_mismatch_penalty: f32,
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
    type Item = Vec<bam::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let source_iterator_loan = &mut self.records;

        // If the underlying iterator is exhausted return None, too
        source_iterator_loan.peek()?;

        Some(
            source_iterator_loan
                .take(self.chunk_size)
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
pub fn run<T: SequenceDifferenceModel + Sync>(
    reads_path: &str,
    reference_path: &str,
    out_file_path: &str,
    alignment_parameters: &AlignmentParameters<T>,
) -> Result<(), Box<dyn Error>> {
    debug!("Load FMD-index");
    let data_fmd_index = UnderlyingDataFMDIndex::load(reference_path)?;
    debug!("Reconstruct FMD-index");
    let fm_index = FMIndex::new(
        &data_fmd_index.bwt,
        &data_fmd_index.less,
        &data_fmd_index.occ,
    );
    let fmd_index = FMDIndex::from(fm_index);

    debug!("Load suffix array");
    let suffix_array: Vec<usize> = {
        let d_suffix_array = snap::Reader::new(File::open(format!("{}.tsa", reference_path))?);
        bincode::deserialize_from(d_suffix_array)?
    };

    debug!("Load position map");
    let identifier_position_map: FastaIdPositions = {
        let d_pi = snap::Reader::new(File::open(format!("{}.tpi", reference_path))?);
        bincode::deserialize_from(d_pi)?
    };

    map_reads(
        alignment_parameters,
        reads_path,
        out_file_path,
        &fmd_index,
        &suffix_array,
        &identifier_position_map,
    )?;

    Ok(())
}

/// Transformations of the input sequence which are used throughout the program.
/// Currently, the input sequence is just converted to uppercase letters.
#[inline]
fn transform_pattern_sequence(record: &bam::Record) -> Vec<u8> {
    record.seq().as_bytes().to_ascii_uppercase()
}

/// Transformations of the base quality inputs which are used throughout the program.
/// Currently, the PHRED-scaled base qualities are assumed to be encoded with ASCII
/// codes starting from 33, so we subtract 33 to transform them to a range starting at 0.
#[inline]
fn transform_base_qualities(record: &bam::Record) -> &[u8] {
    // Subtracting offsets is not required for BAM files
    // record.qual().iter().map(|&f| f - 33).collect::<Vec<_>>()
    record.qual()
}

/// Maps reads and writes them to a file in BAM format
fn map_reads<T: SequenceDifferenceModel + Sync>(
    alignment_parameters: &AlignmentParameters<T>,
    reads_path: &str,
    out_file_path: &str,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
    suffix_array: &Vec<usize>,
    identifier_position_map: &FastaIdPositions,
) -> Result<(), Box<dyn Error>> {
    let mut reads_reader = bam::Reader::from_path(reads_path)?;
    let _ = reads_reader.set_threads(4);

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
        header.push_record(&header_record);
    }

    let allowed_mismatches = AllowedMismatches::new(&alignment_parameters);
    let mut out_file = bam::Writer::from_path(out_file_path, &header, bam::Format::BAM)?;

    debug!("Map reads");
    ChunkIterator::from_reader(reads_reader.records(), alignment_parameters.chunk_size)
        .map(|chunk| {
            trace!("Map chunk of reads in parallel");
            let results = chunk
                .par_iter()
                .map(|record| {
                    let pattern = transform_pattern_sequence(record);
                    (
                        record,
                        k_mismatch_search(
                            &pattern,
                            transform_base_qualities(record),
                            allowed_mismatches.get(pattern.len())
                                * alignment_parameters
                                    .difference_model
                                    .get_representative_mismatch_penalty(),
                            alignment_parameters,
                            fmd_index,
                        ),
                    )
                })
                .map(|(record, mut hit_interval)| {
                    intervals_to_bam(
                        record,
                        &mut hit_interval,
                        fmd_index,
                        suffix_array,
                        identifier_position_map,
                        compute_maximal_possible_score(
                            record,
                            &alignment_parameters.difference_model,
                        ),
                    )
                })
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

/// Convert suffix array intervals to positions and BAM records and write them to a BAM file eventually
fn intervals_to_bam(
    record: &bam::Record,
    intervals: &mut BinaryHeap<HitInterval>,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
    suffix_array: &Vec<usize>,
    identifier_position_map: &FastaIdPositions,
    optimal_score: f32,
) -> Vec<bam::Record> {
    let mut out = Vec::new();

    if let Some(best_alignment) = intervals.pop() {
        let mapping_quality = estimate_mapping_quality(&best_alignment, &intervals, optimal_score);

        // Max. number of alignments per read that will be reported
        let mut max_out_lines_per_read = 1;

        // Read aligns to reference strand
        for &position in best_alignment
            .interval
            .forward()
            .occ(suffix_array)
            .iter()
            .filter(|&&position| position < fmd_index.bwt().len() / 2)
            .take(max_out_lines_per_read)
        {
            max_out_lines_per_read -= 1;
            let (tid, position) = identifier_position_map.get_reference_identifier(position);
            let bam_record = create_bam_record(
                record,
                position,
                Some(&best_alignment),
                Some(mapping_quality),
                tid,
            );
            out.push(bam_record);
        }
        // Read aligns to reverse strand
        for &position in best_alignment
            .interval
            .revcomp()
            .occ(suffix_array)
            .iter()
            .filter(|&&position| position < fmd_index.bwt().len() / 2)
            .take(max_out_lines_per_read)
        {
            let (tid, position) = identifier_position_map.get_reference_identifier(position);
            let mut bam_record = create_bam_record(
                record,
                position,
                Some(&best_alignment),
                Some(mapping_quality),
                tid,
            );
            bam_record.set_reverse();
            out.push(bam_record);
        }
    } else {
        // No match found, report unmapped read
        let mut bam_record = create_bam_record(record, -1, None, None, -1);
        bam_record.set_unmapped();
        out.push(bam_record);
    }
    out
}

/// Computes theoretically maximal possible alignment score per read.
/// With this, actual alignment scores can be viewed as fraction of this
/// optimal score to overcome dependencies on read length, base composition,
/// and the used scoring function.
fn compute_maximal_possible_score<T: SequenceDifferenceModel + Sync>(
    record: &bam::Record,
    difference_model: &T,
) -> f32 {
    let pattern = transform_pattern_sequence(record);
    let base_qualities = transform_base_qualities(record);

    pattern.iter().zip(base_qualities).enumerate().fold(
        0.0,
        |optimal_score, (i, (&base, &quality))| {
            optimal_score + difference_model.get_min_penalty(i, pattern.len(), base, quality, false)
        },
    )
}

/// Estimate mapping quality based on the number of hits for a particular read, its alignment score,
/// and its base qualities
fn estimate_mapping_quality(
    best_alignment: &HitInterval,
    other_alignments: &BinaryHeap<HitInterval>,
    optimal_score: f32,
) -> u8 {
    const MAX_MAPQ: f32 = 37.0;

    let alignment_probability = {
        let ratio_best = 2_f32.powf(best_alignment.alignment_score - optimal_score);
        if best_alignment.interval.size > 1 {
            // Multi-mapping
            ratio_best / best_alignment.interval.size as f32
        } else if other_alignments.is_empty() {
            // Unique mapping
            ratio_best
        } else {
            // Pseudo-unique mapping
            let weighted_suboptimal_alignments =
                other_alignments
                    .iter()
                    .fold(0.0, |acc, suboptimal_alignment| {
                        acc + 2_f32.powf(suboptimal_alignment.alignment_score - optimal_score)
                            * suboptimal_alignment.interval.size as f32
                    });
            ratio_best.powi(2) / (ratio_best + weighted_suboptimal_alignments)
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
fn create_bam_record(
    input_record: &bam::Record,
    position: i32,
    hit_interval: Option<&HitInterval>,
    mapq: Option<u8>,
    tid: i32,
) -> bam::Record {
    let mut bam_record = bam::record::Record::new();

    let cigar = if let Some(hit_interval) = hit_interval {
        Some(&hit_interval.cigar)
    } else {
        None
    };

    // Set mandatory properties of the BAM record
    bam_record.set(
        input_record.qname(),
        cigar,
        &transform_pattern_sequence(input_record),
        &transform_base_qualities(input_record),
    );
    bam_record.set_tid(tid);
    bam_record.set_pos(position);
    if let Some(mapq) = mapq {
        bam_record.set_mapq(mapq);
    }

    // Position of mate (-1 = *)
    bam_record.set_mpos(-1);

    // Reference sequence of mate (-1 = *)
    bam_record.set_mtid(-1);

    // Add optional tags
    if let Some(input_read_group) = input_record.aux(b"RG") {
        bam_record.push_aux(b"RG", &input_read_group);
    }
    if let Some(hit_interval) = hit_interval {
        bam_record.push_aux(
            b"AS",
            &bam::record::Aux::Integer(hit_interval.alignment_score.round() as i64),
        );
        bam_record.push_aux(b"MD", &bam::record::Aux::String(&hit_interval.md_tag));
    }

    bam_record
}

fn build_edit_operation_fields(
    end_node: NodeId,
    edit_tree: &Tree<Option<EditOperation>>,
    pattern_len: usize,
) -> (bam::record::CigarString, Vec<u8>) {
    /// Derive Cigar string from oddly-ordered tracks of edit operations.
    /// Since we start aligning at the center of a read, tracks of edit operations
    /// are not ordered by position in the read.
    fn add_cigar_edit_operation(
        edit_operation: &EditOperation,
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
    };

    fn add_md_edit_operation(
        edit_operation: Option<&EditOperation>,
        last_edit_operation: Option<&EditOperation>,
        mut k: u32,
        md_tag: &mut Vec<u8>,
    ) -> u32 {
        match edit_operation {
            Some(EditOperation::Match(_)) => k += 1,
            Some(EditOperation::Mismatch(_, reference_base)) => {
                md_tag.extend_from_slice(format!("{}{}", k, *reference_base as char).as_bytes());
                k = 0;
            }
            Some(EditOperation::Insertion(_)) => {
                // Insertions are ignored in MD tags
            }
            Some(EditOperation::Deletion(_, reference_base)) => {
                match last_edit_operation {
                    Some(EditOperation::Deletion(_, _)) => {
                        md_tag.extend_from_slice(format!("{}", *reference_base as char).as_bytes());
                    }
                    _ => {
                        md_tag.extend_from_slice(
                            format!("{}^{}", k, *reference_base as char).as_bytes(),
                        );
                    }
                }
                k = 0;
            }
            None => md_tag.extend_from_slice(format!("{}", k).as_bytes()),
        }
        k
    };

    // Restore outer ordering of the edit operation by the positions they carry as values.
    // Whenever there are deletions in the pattern, there is no simple rule to reconstruct the ordering.
    // So, edit operations carrying the same position are pushed onto the same bucket and dealt with later.
    let mut cigar_order_outer: BTreeMap<usize, SmallVec<[EditOperation; 8]>> = BTreeMap::new();

    let mut insert_expanded_cigar = |edit_operation| {
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
    };

    let end_node = edit_tree
        .get(end_node)
        .expect("This is not expected to fail");
    if let Some(last_edit_operation) = end_node.value() {
        insert_expanded_cigar(*last_edit_operation);
    }

    end_node
        .ancestors()
        .map(|node| node.value())
        .filter_map(|&edit_operation| edit_operation)
        .for_each(insert_expanded_cigar);

    // Reconstruct the order of the remaining edit operations and condense CIGAR
    let mut num_matches: u32 = 0;
    let mut num_operations = 1;
    let mut last_edit_operation = None;
    let mut cigar = Vec::new();
    let mut md_tag = Vec::new();
    cigar_order_outer
        .iter()
        .flat_map(|(&i, inner_vec)| {
            if i >= pattern_len / 2 {
                Either::Left(inner_vec.iter().rev())
            } else {
                Either::Right(inner_vec.iter())
            }
        })
        .for_each(|edit_operation| {
            num_matches = add_md_edit_operation(
                Some(&edit_operation),
                last_edit_operation,
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
                            add_cigar_edit_operation(&lop, num_operations, &mut cigar);
                            num_operations = 1;
                            last_edit_operation = Some(edit_operation);
                        }
                    },
                    EditOperation::Mismatch(_, _) => match lop {
                        EditOperation::Mismatch(_, _) | EditOperation::Match(_) => {
                            num_operations += 1;
                        }
                        _ => {
                            add_cigar_edit_operation(&lop, num_operations, &mut cigar);
                            num_operations = 1;
                            last_edit_operation = Some(edit_operation);
                        }
                    },
                    EditOperation::Insertion(_) => match lop {
                        EditOperation::Insertion(_) => {
                            num_operations += 1;
                        }
                        _ => {
                            add_cigar_edit_operation(&lop, num_operations, &mut cigar);
                            num_operations = 1;
                            last_edit_operation = Some(edit_operation);
                        }
                    },
                    EditOperation::Deletion(_, _) => match lop {
                        EditOperation::Deletion(_, _) => {
                            num_operations += 1;
                        }
                        _ => {
                            add_cigar_edit_operation(&lop, num_operations, &mut cigar);
                            num_operations = 1;
                            last_edit_operation = Some(edit_operation);
                        }
                    },
                }
            } else {
                last_edit_operation = Some(edit_operation);
            }
        });
    if let Some(lop) = last_edit_operation {
        add_cigar_edit_operation(&lop, num_operations, &mut cigar);
    }
    let _ = add_md_edit_operation(None, None, num_matches, &mut md_tag);

    (bam::record::CigarString(cigar), md_tag)
}

/// Checks stop-criteria of stack frames before pushing them onto the stack.
/// Since push operations on heaps are costly, this should accelerate the alignment.
fn check_and_push(
    mut stack_frame: MismatchSearchStackFrame,
    pattern: &[u8],
    edit_operation: EditOperation,
    edit_tree: &mut Tree<Option<EditOperation>>,
    stack: &mut BinaryHeap<MismatchSearchStackFrame>,
    intervals: &mut BinaryHeap<HitInterval>,
    penalties: &Penalties,
) {
    // Empty interval
    if stack_frame.current_interval.size < 1 {
        return;
    }

    // Too many mismatches
    if stack_frame.alignment_score + penalties.lower_bound < penalties.max_allowed_penalties {
        return;
    }

    if stop_searching_suboptimal_hits(&stack_frame, intervals, penalties) {
        return;
    }

    stack_frame.edit_node_id = edit_tree
        .get_mut(stack_frame.edit_node_id)
        .expect("This is not expected to fail")
        .append(Some(edit_operation))
        .id();

    // This route through the read graph is finished successfully, push the interval
    if stack_frame.j < 0 || stack_frame.j > (pattern.len() as i16 - 1) {
        let (cigar, md_tag) =
            build_edit_operation_fields(stack_frame.edit_node_id, edit_tree, pattern.len());
        intervals.push(HitInterval {
            interval: stack_frame.current_interval,
            alignment_score: stack_frame.alignment_score,
            cigar,
            md_tag,
        });
        print_debug(&stack_frame, intervals, penalties); // FIXME
        return;
    }

    stack.push(stack_frame);
}

/// FIXME
fn print_debug(
    stack_frame: &MismatchSearchStackFrame,
    intervals: &BinaryHeap<HitInterval>,
    penalties: &Penalties,
) {
    let switch = false; // TODO: Switch me on/off!

    if switch {
        let best_as = match intervals.peek() {
            Some(v) => v.alignment_score,
            None => 0.0,
        };

        eprintln!(
            "{}\t{}\t{}\t{}",
            stack_frame.alignment_score,
            stack_frame.j,
            best_as,
            best_as + penalties.representative_mismatch_penalty,
        );
    }
}

/// If the best scoring interval has a total sum of penalties z, do not search
/// for hits with a minimal expected scored worse than z + representative_mismatch.
/// This speeds up the alignment considerably.
fn stop_searching_suboptimal_hits(
    stack_frame: &MismatchSearchStackFrame,
    hit_intervals: &BinaryHeap<HitInterval>,
    penalties: &Penalties,
) -> bool {
    if let Some(best_scoring_interval) = hit_intervals.peek() {
        if stack_frame.alignment_score + penalties.lower_bound
            < best_scoring_interval.alignment_score + penalties.representative_mismatch_penalty
        {
            return true;
        }
    }
    false
}

/// Finds all suffix array intervals for the current pattern with up to `max_allowed_penalties` mismatch penalties
pub fn k_mismatch_search<T: SequenceDifferenceModel + Sync>(
    pattern: &[u8],
    base_qualities: &[u8],
    max_allowed_penalties: f32,
    parameters: &AlignmentParameters<T>,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
) -> BinaryHeap<HitInterval> {
    let representative_mismatch_penalty = parameters
        .difference_model
        .get_representative_mismatch_penalty();
    let center_of_read = pattern.len() / 2;

    let d_part_pattern = &pattern[..center_of_read as usize];
    let d_backwards = DArray::new(
        d_part_pattern,
        &base_qualities[..center_of_read as usize],
        d_part_pattern.len(),
        pattern.len(),
        Direction::Forward,
        parameters,
        fmd_index,
    );

    let d_part_pattern = &pattern[center_of_read as usize..];
    let d_forwards = DArray::new(
        d_part_pattern,
        &base_qualities[center_of_read as usize..],
        d_part_pattern.len(),
        pattern.len(),
        Direction::Backward,
        parameters,
        fmd_index,
    );

    let mut hit_intervals: BinaryHeap<HitInterval> = BinaryHeap::new();
    let mut stack = BinaryHeap::new();
    let mut edit_tree: Tree<Option<EditOperation>> = Tree::new(None);

    stack.push(MismatchSearchStackFrame {
        j: center_of_read as i16,
        current_interval: fmd_index.init_interval(),
        backward_index: center_of_read as i16 - 1,
        forward_index: center_of_read as i16,
        direction: Direction::Forward,
        gap_backwards: GapState::Closed,
        gap_forwards: GapState::Closed,
        alignment_score: 0.0,
        priority: 0.0,
        edit_node_id: edit_tree.root().id(),
    });

    while let Some(stack_frame) = stack.pop() {
        // In the meantime, have we found a match whose score we perhaps aren't close enough to?
        let lower_bound =
            d_backwards.get(stack_frame.backward_index) + d_forwards.get(stack_frame.forward_index);

        let penalties = Penalties {
            max_allowed_penalties,
            lower_bound,
            representative_mismatch_penalty,
        };

        if stop_searching_suboptimal_hits(&stack_frame, &hit_intervals, &penalties) {
            // Since we operate on a priority stack, it's safe to assume that there are no
            // better scoring frames on the stack, so we are going to stop the search.
            break;
        }

        let next_j;
        let next_backward_index;
        let next_forward_index;

        print_debug(&stack_frame, &hit_intervals, &penalties); // FIXME

        // Determine direction of progress for next iteration on this stack frame
        let fmd_ext_interval = match stack_frame.direction {
            Direction::Forward => {
                next_forward_index = stack_frame.forward_index + 1;
                next_backward_index = stack_frame.backward_index;
                next_j = next_backward_index;
                stack_frame.current_interval.swapped()
            }
            Direction::Backward => {
                next_forward_index = stack_frame.forward_index;
                next_backward_index = stack_frame.backward_index - 1;
                next_j = next_forward_index;
                stack_frame.current_interval
            }
        };

        // Re-calculate the lower bounds for extension
        let lower_bound = d_backwards.get(next_backward_index)
            + d_forwards.get(next_forward_index - pattern.len() as i16 / 2);
        let penalties = Penalties {
            lower_bound,
            ..penalties
        };

        //
        // Insertion in read / deletion in reference
        //
        let penalty = if (stack_frame.gap_backwards == GapState::Insertion
            && stack_frame.direction.is_forward())
            || (stack_frame.gap_forwards == GapState::Insertion
                && stack_frame.direction.is_backward())
        {
            parameters.penalty_gap_extend
        } else {
            parameters.penalty_gap_open
        };

        check_and_push(
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
                alignment_score: stack_frame.alignment_score + penalty,
                priority: stack_frame.alignment_score + penalty + lower_bound,
                ..stack_frame
            },
            pattern,
            EditOperation::Insertion(stack_frame.j as usize),
            &mut edit_tree,
            &mut stack,
            &mut hit_intervals,
            &penalties,
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
            let penalty = if (stack_frame.gap_forwards == GapState::Deletion
                && stack_frame.direction.is_forward())
                || (stack_frame.gap_backwards == GapState::Deletion
                    && stack_frame.direction.is_backward())
            {
                parameters.penalty_gap_extend
            } else {
                parameters.penalty_gap_open
            };

            check_and_push(
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
                    alignment_score: stack_frame.alignment_score + penalty,
                    priority: stack_frame.alignment_score + penalty + lower_bound,
                    ..stack_frame
                },
                pattern,
                EditOperation::Deletion(stack_frame.j as usize, c),
                &mut edit_tree,
                &mut stack,
                &mut hit_intervals,
                &penalties,
            );

            //
            // Match/mismatch
            //
            let penalty = parameters.difference_model.get(
                stack_frame.j as usize,
                pattern.len(),
                c,
                pattern[stack_frame.j as usize],
                base_qualities[stack_frame.j as usize],
            );

            check_and_push(
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
                    alignment_score: stack_frame.alignment_score + penalty,
                    priority: stack_frame.alignment_score + penalty + lower_bound,
                    ..stack_frame
                },
                pattern,
                if c == pattern[stack_frame.j as usize] {
                    EditOperation::Match(stack_frame.j as usize)
                } else {
                    EditOperation::Mismatch(stack_frame.j as usize, c)
                },
                &mut edit_tree,
                &mut stack,
                &mut hit_intervals,
                &penalties,
            );
        }

        // Only search until we found a multi-hit
        if hit_intervals.len() == 1 && hit_intervals.peek().unwrap().interval.size > 1 {
            return hit_intervals;
        }
    }
    hit_intervals
}

#[cfg(test)]
mod tests {
    use bio::{
        alphabets,
        data_structures::{
            bwt::{bwt, less},
            fmindex::{FMDIndex, FMIndex},
            suffix_array::suffix_array,
        },
    };

    use super::*;

    use rust_htslib::bam;

    use crate::sequence_difference_models::VindijaPWM;

    use assert_approx_eq::assert_approx_eq;

    fn build_auxiliary_structures(
        reference: &mut Vec<u8>,
        alphabet: &alphabets::Alphabet,
    ) -> UnderlyingDataFMDIndex {
        let ref_seq_rev_compl = alphabets::dna::revcomp(reference.iter());
        reference.extend_from_slice(b"$");
        reference.extend_from_slice(&ref_seq_rev_compl);
        drop(ref_seq_rev_compl);
        reference.extend_from_slice(b"$");

        let sa = suffix_array(&reference);
        let bwt = bwt(&reference, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);

        UnderlyingDataFMDIndex { bwt, less, occ }
    }

    #[test]
    fn test_exact_search() {
        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);
        let suffix_array = suffix_array(&ref_seq);
        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

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
        struct TestDifferenceModel {}
        impl SequenceDifferenceModel for TestDifferenceModel {
            fn get(
                &self,
                _i: usize,
                _read_length: usize,
                from: u8,
                to: u8,
                _base_quality: u8,
            ) -> f32 {
                if from == b'C' && to == b'T' {
                    return 0.5;
                } else if from != to {
                    return -1.0;
                } else {
                    return 0.0;
                }
            }
        }
        let difference_model = TestDifferenceModel {};

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
        };

        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "ACGTACGTACGTACGT".as_bytes().to_owned();

        // Reference
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let suffix_array = suffix_array(&ref_seq);
        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "GTTC".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let intervals = k_mismatch_search(&pattern, &base_qualities, -1.0, &parameters, &fmd_index);

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
        struct TestDifferenceModel {}
        impl SequenceDifferenceModel for TestDifferenceModel {
            fn get(
                &self,
                _i: usize,
                _read_length: usize,
                from: u8,
                to: u8,
                _base_quality: u8,
            ) -> f32 {
                if from == b'C' && to == b'T' {
                    return -10.0;
                } else if from != to {
                    return -10.0;
                } else {
                    return 1.0;
                }
            }
        }
        let difference_model = TestDifferenceModel {};

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -20.0,
            penalty_gap_extend: -10.0,
            chunk_size: 1,
        };

        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "GAAAAG".as_bytes().to_owned();

        // Reference
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let suffix_array = suffix_array(&ref_seq);
        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "TTTT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let intervals = k_mismatch_search(&pattern, &base_qualities, 1.0, &parameters, &fmd_index);

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
        struct TestDifferenceModel {}
        impl SequenceDifferenceModel for TestDifferenceModel {
            fn get(
                &self,
                _i: usize,
                _read_length: usize,
                from: u8,
                to: u8,
                _base_quality: u8,
            ) -> f32 {
                if from == b'C' && to == b'T' {
                    return -1.0;
                } else if from != to {
                    return -1.0;
                } else {
                    return 1.0;
                }
            }
        }
        let difference_model = TestDifferenceModel {};

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -1.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
        };

        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "ACGTACGTACGTACGT".as_bytes().to_owned();

        // Reference
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "GTTC".as_bytes().to_owned();
        let base_qualities = vec![40; pattern.len()];

        let d_backward = DArray::new(
            &pattern,
            &base_qualities,
            pattern.len(),
            pattern.len(),
            Direction::Forward,
            &parameters,
            &fmd_index,
        );
        let d_forward = DArray::new(
            &pattern,
            &base_qualities,
            pattern.len(),
            pattern.len(),
            Direction::Backward,
            &parameters,
            &fmd_index,
        );

        assert_eq!(&*d_backward.d_array, &[0.0, 0.0, -1.0, -1.0]);
        assert_eq!(&*d_forward.d_array, &[-1.0, -1.0, -1.0, 0.0]);

        let pattern = "GATC".as_bytes().to_owned();

        let d_backward = DArray::new(
            &pattern,
            &base_qualities,
            pattern.len(),
            pattern.len(),
            Direction::Forward,
            &parameters,
            &fmd_index,
        );
        let d_forward = DArray::new(
            &pattern,
            &base_qualities,
            pattern.len(),
            pattern.len(),
            Direction::Backward,
            &parameters,
            &fmd_index,
        );

        assert_eq!(&*d_backward.d_array, &[0.0, -1.0, -1.0, -2.0]);
        assert_eq!(&*d_forward.d_array, &[-2.0, -1.0, -1.0, 0.0]);
    }

    #[test]
    fn test_gapped_alignment() {
        struct TestDifferenceModel {}
        impl SequenceDifferenceModel for TestDifferenceModel {
            fn get(
                &self,
                _i: usize,
                _read_length: usize,
                from: u8,
                to: u8,
                _base_quality: u8,
            ) -> f32 {
                if from == b'C' && to == b'T' {
                    return -10.0;
                } else if from != to {
                    return -10.0;
                } else {
                    return 0.0;
                }
            }
        }
        let difference_model = TestDifferenceModel {};

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
        };

        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "TAT".as_bytes().to_owned(); // revcomp = "ATA"

        // Reference
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let suffix_array = suffix_array(&ref_seq);
        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "TT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let intervals = k_mismatch_search(&pattern, &base_qualities, -2.0, &parameters, &fmd_index);

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
        let difference_model = VindijaPWM::new();

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            // Disable gaps
            penalty_gap_open: -200.0,
            penalty_gap_extend: -100.0,
            chunk_size: 1,
        };

        let alphabet = alphabets::dna::alphabet();
        let mut ref_seq = "CCCCCC".as_bytes().to_owned(); // revcomp = "ATA"

        // Reference
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let sar = suffix_array(&ref_seq);
        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "TTCCCT".as_bytes().to_owned();
        let base_qualities = vec![40; pattern.len()];

        let intervals =
            k_mismatch_search(&pattern, &base_qualities, -30.0, &parameters, &fmd_index);

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_approx_eq!(alignment_score[0], -5.3629);

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&sar))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(positions, vec![0]);

        let pattern = "CCCCCC".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let intervals =
            k_mismatch_search(&pattern, &base_qualities, -30.0, &parameters, &fmd_index);

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_approx_eq!(alignment_score[0], -2.6080124);

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
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "AAGAAA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let intervals =
            k_mismatch_search(&pattern, &base_qualities, -30.0, &parameters, &fmd_index);

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_approx_eq!(alignment_score[0], -10.969394);
    }

    #[test]
    fn test_ord_impl() {
        let map_params_large = MismatchSearchStackFrame {
            alignment_score: -5.0,
            priority: -5.0,
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
            edit_node_id: Tree::new(0).root().id(),
        };
        let map_params_small = MismatchSearchStackFrame {
            alignment_score: -20.0,
            priority: -20.0,
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
            edit_node_id: Tree::new(0).root().id(),
        };

        assert!(map_params_large > map_params_small);
        assert!(!(map_params_large < map_params_small));
        assert_ne!(map_params_large, map_params_small);
    }

    #[test]
    fn test_corner_cases() {
        let difference_model = VindijaPWM::new();

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.01,
            penalty_gap_open: 3.5 * difference_model.get_representative_mismatch_penalty(),
            penalty_gap_extend: 1.5 * difference_model.get_representative_mismatch_penalty(),
            difference_model,
            chunk_size: 1,
        };

        let alphabet = alphabets::dna::alphabet();

        // "correct" "AAAAAAAAAAAAAAAAAAAA" (20x 'A') "incorrect"
        let mut ref_seq = "GTTGTATTTTTAGTAGAGACAGGGTTTCATCATGTTGGCCAGAAAAAAAAAAAAAAAAAAAATTTGTATTTTTAGTAGAGACAGGCTTTCATCATGTTGGCCAG"
            .as_bytes()
            .to_owned();

        // Reference
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let sar = suffix_array(&ref_seq);
        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let allowed_mismatches = AllowedMismatches::new(&parameters);

        let pattern = "GTTGTATTTTTAGTAGAGACAGGCTTTCATCATGTTGGCCAG"
            .as_bytes()
            .to_owned();
        let base_qualities = vec![40; pattern.len()];

        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            allowed_mismatches.get(pattern.len())
                * parameters
                    .difference_model
                    .get_representative_mismatch_penalty(),
            &parameters,
            &fmd_index,
        );

        let alignment_scores = intervals
            .iter()
            .map(|f| f.alignment_score)
            .collect::<Vec<_>>();
        assert_eq!(alignment_scores, vec![-11.320482, -38.763355, -11.348894]);

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
        struct TestDifferenceModel {}
        impl SequenceDifferenceModel for TestDifferenceModel {
            fn get(
                &self,
                _i: usize,
                _read_length: usize,
                from: u8,
                to: u8,
                _base_quality: u8,
            ) -> f32 {
                if from == b'C' && to == b'T' {
                    return -10.0;
                } else if from != to {
                    return -10.0;
                } else {
                    return 0.0;
                }
            }
        }
        let difference_model = TestDifferenceModel {};

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
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
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "ATTACA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut intervals =
            k_mismatch_search(&pattern, &base_qualities, -3.0, &parameters, &fmd_index);
        let best_hit = intervals.pop().unwrap();

        assert_eq!(
            best_hit.cigar,
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
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "GATCAG".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut intervals =
            k_mismatch_search(&pattern, &base_qualities, -3.0, &parameters, &fmd_index);

        let best_hit = intervals.pop().unwrap();

        assert_eq!(best_hit.alignment_score, -3.0);
        assert_eq!(
            best_hit.cigar,
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
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "GATTAGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut intervals =
            k_mismatch_search(&pattern, &base_qualities, -3.0, &parameters, &fmd_index);
        let best_hit = intervals.pop().unwrap();

        assert_eq!(best_hit.alignment_score, -2.0);
        assert_eq!(
            best_hit.cigar,
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
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "GATTAGGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut intervals =
            k_mismatch_search(&pattern, &base_qualities, -3.0, &parameters, &fmd_index);
        let best_hit = intervals.pop().unwrap();

        assert_eq!(best_hit.alignment_score, -3.0);
        assert_eq!(
            best_hit.cigar,
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
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "GATTAGTGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut intervals =
            k_mismatch_search(&pattern, &base_qualities, -4.0, &parameters, &fmd_index);
        let best_hit = intervals.pop().unwrap();

        assert_eq!(best_hit.alignment_score, -4.0);
        assert_eq!(
            best_hit.cigar,
            bam::record::CigarString(vec![
                bam::record::Cigar::Match(5),
                bam::record::Cigar::Ins(3),
                bam::record::Cigar::Match(2)
            ])
        );
    }

    #[test]
    fn test_md_tag() {
        struct TestDifferenceModel {}
        impl SequenceDifferenceModel for TestDifferenceModel {
            fn get(
                &self,
                _i: usize,
                _read_length: usize,
                from: u8,
                to: u8,
                _base_quality: u8,
            ) -> f32 {
                if from == b'C' && to == b'T' {
                    return -1.0;
                } else if from != to {
                    return -2.0;
                } else {
                    return 0.0;
                }
            }
        }
        let difference_model = TestDifferenceModel {};

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
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
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "GATTATA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut intervals =
            k_mismatch_search(&pattern, &base_qualities, -1.0, &parameters, &fmd_index);
        let best_hit = intervals.pop().unwrap();

        assert_eq!(&best_hit.md_tag, &"5C1".as_bytes());

        //
        // Deletion
        //
        let mut ref_seq = "GATTAGCA".as_bytes().to_owned();

        // Reference
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "ATTACA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut intervals =
            k_mismatch_search(&pattern, &base_qualities, -3.0, &parameters, &fmd_index);
        let best_hit = intervals.pop().unwrap();

        assert_eq!(&best_hit.md_tag, &"4^G2".as_bytes());

        //
        // 2-base deletion
        //
        let mut ref_seq = "GATTACAG".as_bytes().to_owned(); // CTGTAATC

        // Reference
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "GATCAG".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut intervals =
            k_mismatch_search(&pattern, &base_qualities, -3.0, &parameters, &fmd_index);

        let best_hit = intervals.pop().unwrap();

        assert_eq!(&best_hit.md_tag, &"3^TA3".as_bytes());

        //
        // Insertion
        //
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "GATTAGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut intervals =
            k_mismatch_search(&pattern, &base_qualities, -3.0, &parameters, &fmd_index);
        let best_hit = intervals.pop().unwrap();

        assert_eq!(&best_hit.md_tag, &"7".as_bytes());

        //
        // 2-base insertion
        //
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "GATTAGGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut intervals =
            k_mismatch_search(&pattern, &base_qualities, -3.0, &parameters, &fmd_index);
        let best_hit = intervals.pop().unwrap();

        assert_eq!(&best_hit.md_tag, &"7".as_bytes());
    }
}
