use std::cmp::Ordering;
use std::collections::binary_heap::BinaryHeap;
use std::error::Error;
use std::fs::File;

use clap::{crate_name, crate_version};
use log::debug;
use smallvec::SmallVec;

use bio::alphabets::dna;
use bio::data_structures::bwt::Occ;
use bio::data_structures::fmindex::{BiInterval, FMDIndex, FMIndex, FMIndexable};

use bincode;
use bio::io::fastq;
use libflate::deflate::Decoder;
use rust_htslib::bam;
use serde::{Deserialize, Serialize};

use crate::sequence_difference_models::SequenceDifferenceModel;
use crate::utils::{AlignmentParameters, AllowedMismatches};

/// Helper struct to bundle index files
struct UnderlyingDataFMDIndex {
    bwt: Vec<u8>,
    less: Vec<usize>,
    occ: Occ,
}

impl UnderlyingDataFMDIndex {
    fn load(path: &str) -> Result<UnderlyingDataFMDIndex, Box<Error>> {
        debug!("Load BWT");
        let bwt: Vec<u8> = {
            let f_bwt = File::open(format!("{}.tbw", path))?;
            let d_bwt = Decoder::new(f_bwt);
            bincode::deserialize_from(d_bwt)?
        };

        debug!("Load \"C\" table");
        let less: Vec<usize> = {
            let f_less = File::open(format!("{}.tle", path))?;
            let d_less = Decoder::new(f_less);
            bincode::deserialize_from(d_less)?
        };

        debug!("Load \"Occ\" table");
        let occ: Occ = {
            let f_occ = File::open(format!("{}.toc", path))?;
            let d_occ = Decoder::new(f_occ);
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
    z: f32,
    edit_operations: EditOperationsTrack,
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
}

#[derive(Debug, Copy, Clone, PartialEq)]
enum EditOperation {
    Insertion,
    Deletion,
    MatchMismatch,
}

#[derive(Debug, Clone)]
struct EditOperationsTrack {
    edit_operations: SmallVec<[EditOperation; 64]>,
}

impl EditOperationsTrack {
    fn new() -> Self {
        EditOperationsTrack {
            edit_operations: SmallVec::new(),
        }
    }

    /// Derive Cigar string from oddly-ordered tracks of edit operations.
    /// Since we start aligning at the center of a read, tracks of edit operations
    /// are not ordered by position in the read.
    fn build_cigar(&self, read_length: usize) -> bam::record::CigarString {
        let mut cigar = Vec::new();

        fn add_edit_operation(
            edit_operation: EditOperation,
            k: u32,
            cigar: &mut Vec<bam::record::Cigar>,
        ) {
            match edit_operation {
                EditOperation::MatchMismatch => cigar.push(bam::record::Cigar::Match(k)),
                EditOperation::Insertion => cigar.push(bam::record::Cigar::Del(k)),
                EditOperation::Deletion => cigar.push(bam::record::Cigar::Ins(k)),
            }
        };

        let mut n = 1;
        let mut last_edit_operation = None;
        for &edit_operation in self
            .edit_operations
            .iter()
            .rev()
            .skip(if read_length % 2 == 0 { 0 } else { 1 })
            .step_by(2)
            .chain(self.edit_operations.iter().step_by(2))
        {
            last_edit_operation = match last_edit_operation {
                Some(last_edit_op) if last_edit_op == edit_operation => {
                    n += 1;
                    last_edit_operation
                }
                Some(last_edit_op) => {
                    add_edit_operation(last_edit_op, n, &mut cigar);
                    n = 1;
                    Some(edit_operation)
                }
                None => Some(edit_operation),
            };
        }
        if let Some(last_edit_operation) = last_edit_operation {
            add_edit_operation(last_edit_operation, n, &mut cigar);
        }

        bam::record::CigarString(cigar)
    }

    fn push(&mut self, edit_operation: EditOperation) {
        self.edit_operations.push(edit_operation);
    }
}

/// Stores information about partial alignments on the priority stack
#[derive(Debug)]
struct MismatchSearchStackFrame {
    j: isize,
    z: f32,
    current_interval: BiInterval,
    backward_index: isize,
    forward_index: isize,
    forward: bool, // TODO: Switch to Direction enum
    open_gap_backwards: bool,
    open_gap_forwards: bool,
    alignment_score: f32,
    edit_operations: Option<EditOperationsTrack>,
    //    debug_helper: String, // Remove this before measuring performance (it's slow)
}

impl PartialOrd for MismatchSearchStackFrame {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MismatchSearchStackFrame {
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

impl PartialEq for MismatchSearchStackFrame {
    fn eq(&self, other: &Self) -> bool {
        (self.alignment_score - other.alignment_score).abs() < std::f32::EPSILON
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

/// Loads index files and launches the mapping process
pub fn run<T: SequenceDifferenceModel>(
    reads_path: &str,
    reference_path: &str,
    out_file_path: &str,
    alignment_parameters: &AlignmentParameters<T>,
) -> Result<(), Box<Error>> {
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
        let f_suffix_array = File::open(format!("{}.tsa", reference_path))?;
        let d_suffix_array = Decoder::new(f_suffix_array);
        bincode::deserialize_from(d_suffix_array)?
    };

    debug!("Load position map");
    let identifier_position_map: FastaIdPositions = {
        let f_pi = File::open(format!("{}.tpi", reference_path))?;
        let d_pi = Decoder::new(f_pi);
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

/// Maps reads and writes them to a file in BAM format
fn map_reads<T: SequenceDifferenceModel>(
    alignment_parameters: &AlignmentParameters<T>,
    reads_path: &str,
    out_file_path: &str,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
    suffix_array: &Vec<usize>,
    identifier_position_map: &FastaIdPositions,
) -> Result<(), Box<Error>> {
    let reads_fq_reader = fastq::Reader::from_file(reads_path)?;

    let mut header = bam::Header::new();
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
        header_record.push_tag(b"PN", &crate_name!());
        header_record.push_tag(b"VN", &crate_version!());
        header.push_record(&header_record);
    }

    let mut out = bam::Writer::from_path(out_file_path, &header)?;

    let mut allowed_mismatches = AllowedMismatches::new(&alignment_parameters);

    debug!("Map reads");
    for record in reads_fq_reader.records() {
        let record = record.unwrap();
        let pattern = record.seq().to_ascii_uppercase();

        // Hardcoded value (33) that should be ok only for Illumina reads
        let base_qualities = record.qual().iter().map(|&f| f - 33).collect::<Vec<_>>();

        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            (allowed_mismatches.get(pattern.len())
                * alignment_parameters
                    .difference_model
                    .get_representative_mismatch_penalty())
            .abs(),
            alignment_parameters,
            &fmd_index,
        );

        //
        // Create BAM records
        //
        if let Some(best_alignment) = intervals.pop() {
            let mapping_quality = estimate_mapping_quality(&best_alignment, &intervals);
            let mut max_out_lines_per_read = 1;
            // Aligns to reference strand
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
                let record = create_bam_record(
                    record.id().as_bytes(),
                    record.seq(),
                    record.qual().iter(),
                    position,
                    Some(&best_alignment),
                    Some(mapping_quality),
                    tid,
                );
                out.write(&record)?;
            }
            // Aligns to reverse strand
            for &position in best_alignment
                .interval
                .revcomp()
                .occ(suffix_array)
                .iter()
                .filter(|&&position| position < fmd_index.bwt().len() / 2)
                .take(max_out_lines_per_read)
            {
                let (tid, position) = identifier_position_map.get_reference_identifier(position);
                let mut record = create_bam_record(
                    record.id().as_bytes(),
                    &dna::revcomp(record.seq()),
                    record.qual().iter().rev(),
                    position,
                    Some(&best_alignment),
                    Some(mapping_quality),
                    tid,
                );
                record.set_reverse();
                out.write(&record)?;
            }
        } else {
            // No match found
            let mut record = create_bam_record(
                record.id().as_bytes(),
                record.seq(),
                record.qual().iter(),
                -1,
                None,
                None,
                -1,
            );
            record.set_unmapped();
            out.write(&record)?;
        }
    }
    debug!("Done");
    Ok(())
}

/// Estimate mapping quality based on the number of hits for a particular read, its alignment score,
/// and its base qualities
fn estimate_mapping_quality(
    best_alignment: &HitInterval,
    other_alignments: &BinaryHeap<HitInterval>,
) -> u8 {
    // Multi-mapping
    if best_alignment.interval.size > 1 {
        (-10_f32 * (1.0 - (1.0 / best_alignment.interval.size as f32)).log10()).round() as u8
    } else {
        // "Unique" mapping
        if let Some(second_alignment) = other_alignments.peek() {
            let total_number = other_alignments
                .iter()
                .take_while(|hit_interval| {
                    hit_interval.alignment_score.round() >= second_alignment.alignment_score.round()
                })
                .fold(0, |agg, hit_interval| agg + hit_interval.interval.size);
            let nominator = 1.0;
            let denominator = total_number as f32;
            let correction =
                (best_alignment.alignment_score - second_alignment.alignment_score).abs();
            (-10_f32 * (1.0 - (nominator / denominator)).log10() + correction).round() as u8
        } else {
            37
        }
    }
}

/// Create and return a BAM record of either a hit or an unmapped read
fn create_bam_record<'a, T: Iterator<Item = &'a u8>>(
    input_name: &[u8],
    input_seq: &[u8],
    input_quality: T,
    position: i32,
    hit_interval: Option<&HitInterval>,
    mapq: Option<u8>,
    tid: i32,
) -> bam::Record {
    let mut bam_record = bam::record::Record::new();

    let cigar = match hit_interval {
        Some(hit_interval) => hit_interval.edit_operations.build_cigar(input_seq.len()),
        None => bam::record::CigarString::from_str("").unwrap(),
    };

    bam_record.set(
        input_name,
        &cigar,
        input_seq,
        input_quality
            .map(|&x| x - 33)
            .collect::<Vec<_>>()
            .as_slice(),
    );

    bam_record.set_tid(tid);
    bam_record.set_pos(position);

    if let Some(mapq) = mapq {
        bam_record.set_mapq(mapq); // Mapping quality
    }

    bam_record.set_mpos(-1); // Position of mate (-1 = *)
    bam_record.set_mtid(-1); // Reference sequence of mate (-1 = *)

    if let Some(hit_interval) = hit_interval {
        bam_record
            .push_aux(
                b"AS",
                &bam::record::Aux::Integer(hit_interval.alignment_score.round() as i64),
            )
            .unwrap();
    }

    bam_record
}

/// Checks stop-criteria of stack frames before pushing them onto the stack.
/// Since push operations on heaps are costly, this should accelerate the alignment.
fn check_and_push(
    mut stack_frame: MismatchSearchStackFrame,
    pattern: &[u8],
    edit_operation: EditOperation,
    edit_operations: &Option<EditOperationsTrack>,
    stack: &mut BinaryHeap<MismatchSearchStackFrame>,
    intervals: &mut BinaryHeap<HitInterval>,
    d_backwards: &[f32],
    d_forwards: &[f32],
    representative_mismatch_penalty: f32,
) {
    // Empty interval
    if stack_frame.current_interval.size < 1 {
        return;
    }

    // Too many mismatches
    let backwards_lower_bound = match d_backwards.get(stack_frame.backward_index as usize) {
        Some(&d_i) => d_i,
        None => 0.0,
    };
    let forwards_lower_bound =
        match d_forwards.get((stack_frame.forward_index as usize - pattern.len() / 2) as usize) {
            Some(&d_i) => d_i,
            None => 0.0,
        };
    if stack_frame.z < backwards_lower_bound + forwards_lower_bound {
        return;
    }

    if stop_searching(&stack_frame, intervals, representative_mismatch_penalty) {
        return;
    }

    let mut edit_operations = edit_operations.to_owned().unwrap();
    edit_operations.push(edit_operation);

    // This route through the read graph is finished successfully, push the interval
    if stack_frame.j < 0 || stack_frame.j > (pattern.len() as isize - 1) {
        intervals.push(HitInterval {
            interval: stack_frame.current_interval,
            alignment_score: stack_frame.alignment_score - pattern.len() as f32,
            z: stack_frame.z,
            edit_operations,
        });
        print_debug(&stack_frame, intervals, representative_mismatch_penalty); // FIXME
        return;
    }

    stack_frame.edit_operations = Some(edit_operations);
    stack.push(stack_frame);
}

/// FIXME
fn print_debug(
    stack_frame: &MismatchSearchStackFrame,
    intervals: &BinaryHeap<HitInterval>,
    representative_mismatch_penalty: f32,
) {
    let switch = false; // TODO: Switch me on/off!

    if switch {
        let best_z = match intervals.peek() {
            Some(v) => v.z,
            None => 0.0,
        };

        eprintln!(
            "{}\t{}\t{}\t{}\t{}",
            stack_frame.alignment_score,
            stack_frame.j,
            stack_frame.z,
            best_z,
            best_z + representative_mismatch_penalty,
        );
    }
}

/// If the best scoring interval has a total sum of penalties z, do not search
/// for hits scored worse than z + representative_mismatch.
/// This speeds up the alignment considerably.
fn stop_searching(
    stack_frame: &MismatchSearchStackFrame,
    hit_intervals: &BinaryHeap<HitInterval>,
    representative_mismatch_penalty: f32,
) -> bool {
    if let Some(best_scoring_interval) = hit_intervals.peek() {
        if stack_frame.z < best_scoring_interval.z + representative_mismatch_penalty {
            return true;
        }
    }
    false
}

/// Finds all suffix array intervals for the current pattern with up to z mismatch penalties
pub fn k_mismatch_search<T: SequenceDifferenceModel>(
    pattern: &[u8],
    _base_qualities: &[u8],
    z: f32,
    parameters: &AlignmentParameters<T>,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
) -> BinaryHeap<HitInterval> {
    let representative_mismatch_penalty = parameters
        .difference_model
        .get_representative_mismatch_penalty();
    let center_of_read = pattern.len() as isize / 2;

    let d_part_pattern = &pattern[..center_of_read as usize];
    let d_backwards = calculate_d(
        d_part_pattern.iter(),
        d_part_pattern.len(),
        pattern.len(),
        Direction::Forward,
        parameters,
        fmd_index,
    );

    let d_part_pattern = &pattern[center_of_read as usize..];
    let d_forwards = calculate_d(
        d_part_pattern.iter().rev(),
        d_part_pattern.len(),
        pattern.len(),
        Direction::Backward,
        parameters,
        fmd_index,
    )
    .into_iter()
    .rev()
    .collect::<Vec<_>>();

    let mut hit_intervals: BinaryHeap<HitInterval> = BinaryHeap::new();
    let mut stack = BinaryHeap::new();

    stack.push(MismatchSearchStackFrame {
        j: center_of_read,
        z,
        current_interval: fmd_index.init_interval(),
        backward_index: center_of_read - 1,
        forward_index: center_of_read,
        forward: true,
        open_gap_backwards: false,
        open_gap_forwards: false,
        alignment_score: 0.0,
        edit_operations: Some(EditOperationsTrack::new()),
        //        debug_helper: String::from("."),
    });

    while let Some(stack_frame) = stack.pop() {
        if stop_searching(
            &stack_frame,
            &hit_intervals,
            representative_mismatch_penalty,
        ) {
            continue;
        }

        let next_j;
        let next_backward_index;
        let next_forward_index;

        print_debug(
            &stack_frame,
            &hit_intervals,
            representative_mismatch_penalty,
        ); // FIXME

        // Determine direction of progress for next iteration on this stack frame
        let fmd_ext_interval = if stack_frame.forward {
            next_forward_index = stack_frame.forward_index + 1;
            next_backward_index = stack_frame.backward_index;
            next_j = next_backward_index;
            stack_frame.current_interval.swapped()
        } else {
            next_forward_index = stack_frame.forward_index;
            next_backward_index = stack_frame.backward_index - 1;
            next_j = next_forward_index;
            stack_frame.current_interval
        };

        //
        // Insertion in read / deletion in reference
        //
        let penalty = if (stack_frame.open_gap_backwards && stack_frame.forward)
            || (stack_frame.open_gap_forwards && !stack_frame.forward)
        {
            parameters.penalty_gap_extend
        } else {
            parameters.penalty_gap_open
        };

        check_and_push(
            MismatchSearchStackFrame {
                j: next_j,
                z: stack_frame.z + penalty,
                backward_index: next_backward_index,
                forward_index: next_forward_index,
                forward: !stack_frame.forward,
                // Mark opened gap at the corresponding end
                open_gap_backwards: if !stack_frame.forward {
                    stack_frame.open_gap_backwards
                } else {
                    true
                },
                open_gap_forwards: if stack_frame.forward {
                    true
                } else {
                    stack_frame.open_gap_forwards
                },
                alignment_score: stack_frame.alignment_score + penalty + 1.0,
                edit_operations: None,
                //            debug_helper: if stack_frame.forward {
                //                format!("{}(_)", stack_frame.debug_helper)
                //            } else {
                //                format!("(_){}", stack_frame.debug_helper)
                //            },
                ..stack_frame
            },
            pattern,
            EditOperation::Insertion,
            &stack_frame.edit_operations,
            &mut stack,
            &mut hit_intervals,
            &d_backwards,
            &d_forwards,
            representative_mismatch_penalty,
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
            let c = if stack_frame.forward {
                interval_prime = interval_prime.swapped();
                dna::complement(c)
            } else {
                c
            };

            //
            // Deletion in read / insertion in reference
            //
            let penalty = if (stack_frame.open_gap_backwards && stack_frame.forward)
                || (stack_frame.open_gap_forwards && !stack_frame.forward)
            {
                parameters.penalty_gap_extend
            } else {
                parameters.penalty_gap_open
            };

            check_and_push(
                MismatchSearchStackFrame {
                    z: stack_frame.z + penalty,
                    current_interval: interval_prime,
                    // Mark open gap at the corresponding end
                    open_gap_backwards: if !stack_frame.forward {
                        true
                    } else {
                        stack_frame.open_gap_backwards
                    },
                    open_gap_forwards: if stack_frame.forward {
                        true
                    } else {
                        stack_frame.open_gap_forwards
                    },
                    alignment_score: stack_frame.alignment_score + penalty + 1.0,
                    edit_operations: None,
                    //                debug_helper: if stack_frame.forward {
                    //                    format!("{}({})", stack_frame.debug_helper, c as char)
                    //                } else {
                    //                    format!("({}){}", c as char, stack_frame.debug_helper)
                    //                },
                    ..stack_frame
                },
                pattern,
                EditOperation::Deletion,
                &stack_frame.edit_operations,
                &mut stack,
                &mut hit_intervals,
                &d_backwards,
                &d_forwards,
                representative_mismatch_penalty,
            );

            //
            // Match/mismatch
            //
            let penalty = parameters.difference_model.get(
                stack_frame.j as usize,
                pattern.len(),
                c,
                pattern[stack_frame.j as usize],
            );

            check_and_push(
                MismatchSearchStackFrame {
                    j: next_j,
                    z: stack_frame.z + penalty.min(0.0),
                    current_interval: interval_prime,
                    backward_index: next_backward_index,
                    forward_index: next_forward_index,
                    forward: !stack_frame.forward,
                    // Mark closed gap at the corresponding end
                    open_gap_backwards: if !stack_frame.forward {
                        false
                    } else {
                        stack_frame.open_gap_backwards
                    },
                    open_gap_forwards: if stack_frame.forward {
                        false
                    } else {
                        stack_frame.open_gap_forwards
                    },
                    alignment_score: stack_frame.alignment_score + penalty + 1.0,
                    edit_operations: None,
                    //                    debug_helper: if stack_frame.forward {
                    //                        format!(
                    //                            "{}{}",
                    //                            stack_frame.debug_helper,
                    //                            c.to_ascii_lowercase() as char
                    //                        )
                    //                    } else {
                    //                        format!(
                    //                            "{}{}",
                    //                            c.to_ascii_lowercase() as char,
                    //                            stack_frame.debug_helper
                    //                        )
                    //                    },
                },
                pattern,
                EditOperation::MatchMismatch,
                &stack_frame.edit_operations,
                &mut stack,
                &mut hit_intervals,
                &d_backwards,
                &d_forwards,
                representative_mismatch_penalty,
            );
        }

        // Only search until we found a multi-hit (equal MAPQs) or more hits
        // with different MAPQs
        match hit_intervals.len() {
            1 if hit_intervals.peek().unwrap().interval.size > 1 => {
                return hit_intervals;
            }
            _ => {}
        }
    }
    hit_intervals
}

/// A reversed FMD-index is used to compute the lower bound of mismatches of a read per position.
/// This allows for pruning the search tree.
fn calculate_d<'a, T: Iterator<Item = &'a u8>, U: SequenceDifferenceModel>(
    pattern: T,
    pattern_length: usize,
    read_length: usize,
    direction: Direction,
    alignment_parameters: &AlignmentParameters<U>,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
) -> Vec<f32> {
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
    pattern
        .enumerate()
        .map(|(index, &base)| {
            interval = extend_interval(fmd_index, &interval, base);
            if interval.size < 1 {
                z -= alignment_parameters
                    .difference_model
                    .get_min_penalty(directed_pattern_index(index), read_length, base)
                    .max(alignment_parameters.penalty_gap_open)
                    .max(alignment_parameters.penalty_gap_extend);
                interval = fmd_index.init_interval();
            }
            z
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use bio::alphabets;
    use bio::data_structures::bwt::{bwt, less};
    use bio::data_structures::fmindex::{FMDIndex, FMIndex};
    use bio::data_structures::suffix_array::suffix_array;

    use super::*;

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
        assert_eq!(vec![10, 3], positions_1);

        // Test non-occurring pattern
        let pattern_2 = b"GG";
        let sai_2 = fmd_index.backward_search(pattern_2.iter());
        let positions_2 = sai_2.occ(&suffix_array);
        assert_eq!(Vec::<usize>::new(), positions_2);
    }

    #[test]
    fn test_inexact_search() {
        struct TestDifferenceModel {}
        impl SequenceDifferenceModel for TestDifferenceModel {
            fn new() -> Self {
                TestDifferenceModel {}
            }
            fn get(&self, _i: usize, _read_length: usize, from: u8, to: u8) -> f32 {
                if from == b'C' && to == b'T' {
                    return 0.0;
                } else if from != to {
                    return -1.0;
                } else {
                    return 1.0;
                }
            }
        }
        let difference_model = TestDifferenceModel::new();

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
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

        let intervals = k_mismatch_search(&pattern, &base_qualities, 1.0, &parameters, &fmd_index);

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_eq!(vec![2.0], alignment_score);

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&suffix_array))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(vec![2, 6, 10, 19, 23, 27], positions);
    }

    #[test]
    fn test_reverse_strand_search() {
        struct TestDifferenceModel {}
        impl SequenceDifferenceModel for TestDifferenceModel {
            fn new() -> Self {
                TestDifferenceModel {}
            }
            fn get(&self, _i: usize, _read_length: usize, from: u8, to: u8) -> f32 {
                if from == b'C' && to == b'T' {
                    return -10.0;
                } else if from != to {
                    return -10.0;
                } else {
                    return 1.0;
                }
            }
        }
        let difference_model = TestDifferenceModel::new();

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -20.0,
            penalty_gap_extend: -10.0,
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
        assert_eq!(vec![8], positions);
    }

    #[test]
    fn test_d() {
        struct TestDifferenceModel {}
        impl SequenceDifferenceModel for TestDifferenceModel {
            fn new() -> Self {
                TestDifferenceModel {}
            }
            fn get(&self, _i: usize, _read_length: usize, from: u8, to: u8) -> f32 {
                if from == b'C' && to == b'T' {
                    return -1.0;
                } else if from != to {
                    return -1.0;
                } else {
                    return 1.0;
                }
            }
        }
        let difference_model = TestDifferenceModel::new();

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -1.0,
            penalty_gap_extend: -1.0,
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

        let d_backward = calculate_d(
            pattern.iter(),
            pattern.len(),
            pattern.len(),
            Direction::Forward,
            &parameters,
            &fmd_index,
        );
        let d_forward = calculate_d(
            pattern.iter().rev(),
            pattern.len(),
            pattern.len(),
            Direction::Backward,
            &parameters,
            &fmd_index,
        )
        .into_iter()
        .rev()
        .collect::<Vec<_>>();

        assert_eq!(vec![0.0, 0.0, 1.0, 1.0], d_backward);
        assert_eq!(vec![1.0, 1.0, 1.0, 0.0], d_forward);

        let pattern = "GATC".as_bytes().to_owned();

        let d_backward = calculate_d(
            pattern.iter(),
            pattern.len(),
            pattern.len(),
            Direction::Forward,
            &parameters,
            &fmd_index,
        );
        let d_forward = calculate_d(
            pattern.iter().rev(),
            pattern.len(),
            pattern.len(),
            Direction::Backward,
            &parameters,
            &fmd_index,
        )
        .into_iter()
        .rev()
        .collect::<Vec<_>>();

        assert_eq!(vec![0.0, 1.0, 1.0, 2.0], d_backward);
        assert_eq!(vec![2.0, 1.0, 1.0, 0.0], d_forward);
    }

    #[test]
    fn test_gapped_alignment() {
        struct TestDifferenceModel {}
        impl SequenceDifferenceModel for TestDifferenceModel {
            fn new() -> Self {
                TestDifferenceModel {}
            }
            fn get(&self, _i: usize, _read_length: usize, from: u8, to: u8) -> f32 {
                if from == b'C' && to == b'T' {
                    return -10.0;
                } else if from != to {
                    return -10.0;
                } else {
                    return 1.0;
                }
            }
        }
        let difference_model = TestDifferenceModel::new();

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
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

        let intervals = k_mismatch_search(&pattern, &base_qualities, 2.0, &parameters, &fmd_index);

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&suffix_array))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(vec![0, 2, 5], positions);
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
        let base_qualities = vec![0; pattern.len()];

        let intervals = k_mismatch_search(&pattern, &base_qualities, 30.0, &parameters, &fmd_index);

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_approx_eq!(-5.3629, alignment_score[0]);

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&sar))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(vec![0], positions);

        let pattern = "CCCCCC".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let intervals = k_mismatch_search(&pattern, &base_qualities, 30.0, &parameters, &fmd_index);

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_approx_eq!(-2.6080124, alignment_score[0]);

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&sar))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(vec![0], positions);

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

        let intervals = k_mismatch_search(&pattern, &base_qualities, 30.0, &parameters, &fmd_index);

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_approx_eq!(-10.969394, alignment_score[0]);
    }

    #[test]
    fn test_ord_impl() {
        let map_params_large = MismatchSearchStackFrame {
            alignment_score: -5.0,
            j: 0,
            z: 0.0,
            current_interval: BiInterval {
                lower: 5,
                lower_rev: 5,
                match_size: 5,
                size: 5,
            },
            backward_index: 5,
            forward_index: 5,
            forward: false,
            open_gap_backwards: false,
            open_gap_forwards: false,
            edit_operations: None,
        };
        let map_params_small = MismatchSearchStackFrame {
            alignment_score: -20.0,
            j: 0,
            z: 0.0,
            current_interval: BiInterval {
                lower: 5,
                lower_rev: 5,
                match_size: 5,
                size: 5,
            },
            backward_index: 5,
            forward_index: 5,
            forward: false,
            open_gap_backwards: false,
            open_gap_forwards: false,
            edit_operations: None,
        };

        assert!(map_params_large > map_params_small);
        assert!(!(map_params_large < map_params_small));
        assert_ne!(map_params_large, map_params_small);
    }
}
