use std::cmp::Ordering;
use std::collections::binary_heap::BinaryHeap;
use std::error::Error;
use std::fs::File;

use log::debug;

use bio::alphabets::dna;
use bio::data_structures::bwt::Occ;
use bio::data_structures::fmindex::{BiInterval, FMDIndex, FMIndex, FMIndexable};

use bincode::deserialize_from;
use bio::io::fastq;
use libflate::deflate::Decoder;
use rust_htslib::bam;

use crate::sequence_difference_models::SequenceDifferenceModel;
use crate::utils::{AlignmentParameters, AllowedMismatches};

struct UnderlyingDataFMDIndex {
    bwt: Vec<u8>,
    less: Vec<usize>,
    occ: Occ,
}

impl UnderlyingDataFMDIndex {
    fn load(name: &str) -> Result<UnderlyingDataFMDIndex, Box<Error>> {
        debug!("Load BWT");
        let f_bwt = File::open(format!("{}.bwt", name))?;
        let d_bwt = Decoder::new(f_bwt);
        let bwt: Vec<u8> = deserialize_from(d_bwt)?;

        debug!("Load \"C\" table");
        let f_less = File::open(format!("{}.less", name))?;
        let d_less = Decoder::new(f_less);
        let less: Vec<usize> = deserialize_from(d_less)?;

        debug!("Load \"Occ\" table");
        let f_occ = File::open(format!("{}.occ", name))?;
        let d_occ = Decoder::new(f_occ);
        let occ: Occ = deserialize_from(d_occ)?;

        Ok(UnderlyingDataFMDIndex { bwt, less, occ })
    }
}

#[derive(Debug)]
pub struct IntervalQuality {
    interval: BiInterval,
    alignment_score: f32,
}

#[derive(Debug)]
struct MismatchSearchParameters {
    j: isize,
    z: f32,
    current_interval: BiInterval,
    backward_index: isize,
    forward_index: isize,
    forward: bool,
    open_gap_backwards: bool,
    open_gap_forwards: bool,
    alignment_score: f32,
    //    debug_helper: String, // Remove this before measuring performance (it's very slow)
}

impl PartialOrd for MismatchSearchParameters {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MismatchSearchParameters {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.alignment_score > other.alignment_score {
            Ordering::Less
        } else if self.alignment_score < other.alignment_score {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}

impl PartialEq for MismatchSearchParameters {
    fn eq(&self, other: &Self) -> bool {
        self.alignment_score == other.alignment_score
    }
}

impl Eq for MismatchSearchParameters {}

pub fn run<T: SequenceDifferenceModel>(
    reads_path: &str,
    alignment_parameters: &AlignmentParameters<T>,
) -> Result<(), Box<Error>> {
    debug!("Load FMD-index");
    let data_fmd_index = UnderlyingDataFMDIndex::load("ref")?;
    debug!("Reconstruct FMD-index");
    let fm_index = FMIndex::new(
        &data_fmd_index.bwt,
        &data_fmd_index.less,
        &data_fmd_index.occ,
    );
    let fmd_index = FMDIndex::from(fm_index);

    debug!("Load reverse FMD-index");
    let rev_data_fmd_index = UnderlyingDataFMDIndex::load("rev_ref")?;
    debug!("Reconstruct reverse FMD-index");
    let rev_fm_index = FMIndex::new(
        &rev_data_fmd_index.bwt,
        &rev_data_fmd_index.less,
        &rev_data_fmd_index.occ,
    );
    let rev_fmd_index = FMDIndex::from(rev_fm_index);

    debug!("Load suffix array");
    let f_suffix_array = File::open("ref.sar")?;
    let d_suffix_array = Decoder::new(f_suffix_array);
    let suffix_array: Vec<usize> = deserialize_from(d_suffix_array)?;

    map_reads(
        alignment_parameters,
        reads_path,
        &fmd_index,
        &rev_fmd_index,
        &suffix_array,
    )?;

    Ok(())
}

fn map_reads<T: SequenceDifferenceModel>(
    alignment_parameters: &AlignmentParameters<T>,
    reads_path: &str,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
    rev_fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
    suffix_array: &Vec<usize>,
) -> Result<(), Box<Error>> {
    let reads_fq_reader = fastq::Reader::from_file(reads_path)?;

    let mut header = bam::Header::new(); // TODO: Don't hardcode header entries
    {
        let mut header_record = bam::header::HeaderRecord::new(b"SQ");
        header_record.push_tag(b"SN", &"chr22");
        header_record.push_tag(b"LN", &"51304566");
        header.push_record(&header_record);
    }
    {
        let mut header_record = bam::header::HeaderRecord::new(b"PG");
        header_record.push_tag(b"ID", &"thrust");
        header_record.push_tag(b"PN", &"thrust");
        header_record.push_tag(b"VN", &"0.0.0");
        header.push_record(&header_record);
    }
    let mut out = bam::Writer::from_path("out.bam", &header).unwrap();

    let mut allowed_mismatches = AllowedMismatches::new(&alignment_parameters);

    debug!("Map reads");
    for record in reads_fq_reader.records() {
        let record = record.unwrap();
        let pattern = record.seq().to_ascii_uppercase();

        // Hardcoded value (33) that should be ok only for Illumina reads
        let base_qualities = record.qual().iter().map(|&f| f - 33).collect::<Vec<_>>();

        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            (allowed_mismatches.get(pattern.len())
                * alignment_parameters
                    .difference_model
                    .get_representative_mismatch_penalty())
            .abs(),
            alignment_parameters,
            &fmd_index,
            &rev_fmd_index,
        );

        const TEN_F32: f32 = 10.0;
        let mut sum_base_q_best = std::f32::MAX;
        let mut sum_base_q_all = 0.0;
        for sum_of_qualities in intervals.iter() {
            sum_base_q_all += TEN_F32.powf(-(sum_of_qualities.alignment_score));
            if sum_of_qualities.alignment_score < sum_base_q_best {
                sum_base_q_best = sum_of_qualities.alignment_score
            }
        }
        let _mapping_quality =
            -((1.0 - (TEN_F32.powf(-(sum_base_q_best)) / sum_base_q_all)).log10());

        // Create BAM records
        // Aligns to reference strand
        for imm in intervals.iter() {
            for &position in imm.interval.forward().occ(suffix_array).iter() {
                if position < fmd_index.bwt().len() / 2 {
                    out.write(&create_bam_record(
                        record.id().as_bytes(),
                        record.seq(),
                        record.qual(),
                        position,
                    ))?;
                }
            }
            // Aligns to reverse strand
            for &position in imm.interval.revcomp().occ(suffix_array).iter() {
                if position < fmd_index.bwt().len() / 2 {
                    let mut record = create_bam_record(
                        record.id().as_bytes(),
                        record.seq(),
                        record.qual(),
                        position,
                    );
                    record.set_reverse();
                    out.write(&record)?;
                }
            }
        }
    }
    debug!("Done");
    Ok(())
}

fn create_bam_record(
    input_name: &[u8],
    input_seq: &[u8],
    input_qual: &[u8],
    position: usize,
) -> bam::Record {
    let mut bam_record = bam::record::Record::new();
    let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(input_seq.len() as u32)]); // FIXME: Calculate CIGAR string appropriately
    bam_record.set(
        input_name,
        &cigar,
        input_seq,
        input_qual
            .iter()
            .map(|&x| x - 33)
            .collect::<Vec<_>>()
            .as_slice(),
    );
    bam_record.set_pos(position as i32);
    bam_record.set_mapq(37); // Mapping quality
    bam_record.set_mpos(-1); // Position of mate (-1 = *)
    bam_record.set_mtid(-1); // Reference sequence of mate (-1 = *)
    bam_record
}

/// Finds all suffix array intervals for the current pattern with up to z mismatches
pub fn k_mismatch_search<T: SequenceDifferenceModel>(
    pattern: &[u8],
    _base_qualities: &[u8],
    z: f32,
    parameters: &AlignmentParameters<T>,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
    rev_fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
) -> Vec<IntervalQuality> {
    let center_of_read = pattern.len() as isize / 2;

    let d_backwards = calculate_d(
        pattern[..center_of_read as usize].iter(),
        pattern.len(),
        parameters,
        rev_fmd_index,
    );
    let d_forwards = calculate_d(
        pattern[center_of_read as usize..].iter().rev(),
        pattern.len(),
        parameters,
        fmd_index,
    )
    .into_iter()
    .rev()
    .collect::<Vec<_>>();

    let mut intervals = Vec::new();

    let mut stack = BinaryHeap::new();
    stack.push(MismatchSearchParameters {
        j: center_of_read,
        z,
        current_interval: fmd_index.init_interval(),
        backward_index: center_of_read - 1,
        forward_index: center_of_read,
        forward: true,
        open_gap_backwards: false,
        open_gap_forwards: false,
        alignment_score: 0.0,
        //        debug_helper: String::from("."),
    });

    while let Some(stack_frame) = stack.pop() {
        // Too many mismatches
        let backwards_lower_bound = match d_backwards.get(stack_frame.backward_index as usize) {
            Some(&d_i) => d_i,
            None => 0.0,
        };
        let forwards_lower_bound =
            match d_forwards.get((stack_frame.forward_index - center_of_read) as usize) {
                Some(&d_i) => d_i,
                None => 0.0,
            };
        if stack_frame.z < backwards_lower_bound + forwards_lower_bound {
            continue;
        }

        // This route through the read graph is finished successfully, push the interval
        if stack_frame.j < 0 || stack_frame.j > (pattern.len() as isize - 1) {
            intervals.push(IntervalQuality {
                interval: stack_frame.current_interval,
                alignment_score: stack_frame.alignment_score,
            });
            continue; // FIXME
        }

        let next_j;
        let next_backward_index;
        let next_forward_index;

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

        // Insertion in read
        let penalty = if (stack_frame.open_gap_backwards && stack_frame.forward)
            || (stack_frame.open_gap_forwards && !stack_frame.forward)
        {
            parameters.penalty_gap_extend
        } else {
            parameters.penalty_gap_open
        };
        stack.push(MismatchSearchParameters {
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
            alignment_score: stack_frame.alignment_score + penalty,
            //            debug_helper: if stack_frame.forward {
            //                format!("{}(_)", stack_frame.debug_helper)
            //            } else {
            //                format!("(_){}", stack_frame.debug_helper)
            //            },
            ..stack_frame
        });

        let mut s = 0;
        let mut o;
        let mut l = fmd_ext_interval.lower_rev;
        for &c in b"$TGCNA".iter() {
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
                if c == b'$' || c == b'N' || s < 1 {
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

            // Deletion in read
            let penalty = if (stack_frame.open_gap_backwards && stack_frame.forward)
                || (stack_frame.open_gap_forwards && !stack_frame.forward)
            {
                parameters.penalty_gap_extend
            } else {
                parameters.penalty_gap_open
            };
            stack.push(MismatchSearchParameters {
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
                alignment_score: stack_frame.alignment_score + penalty,
                //                debug_helper: if stack_frame.forward {
                //                    format!("{}({})", stack_frame.debug_helper, c as char)
                //                } else {
                //                    format!("({}){}", c as char, stack_frame.debug_helper)
                //                },
                ..stack_frame
            });

            // Match
            if c == pattern[stack_frame.j as usize] {
                stack.push(MismatchSearchParameters {
                    j: next_j,
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
                    alignment_score: stack_frame.alignment_score
                        + parameters.difference_model.get(
                            stack_frame.j as usize,
                            pattern.len(),
                            c,
                            pattern[stack_frame.j as usize],
                        ),
                    //                    debug_helper: if stack_frame.forward {
                    //                        format!("{}{}", stack_frame.debug_helper, c as char)
                    //                    } else {
                    //                        format!("{}{}", c as char, stack_frame.debug_helper)
                    //                    },
                    ..stack_frame
                });

            // Mismatch
            } else {
                let penalty = parameters.difference_model.get(
                    stack_frame.j as usize,
                    pattern.len(),
                    c,
                    pattern[stack_frame.j as usize],
                );
                stack.push(MismatchSearchParameters {
                    j: next_j,
                    z: stack_frame.z + penalty,
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
                    alignment_score: stack_frame.alignment_score + penalty,
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
                });
            }
        }
    }
    intervals
}

/// A reversed FMD-index is used to compute the lower bound of mismatches of a read per position.
/// This allows for pruning the search tree. Implementation follows closely Li & Durbin (2009).
fn calculate_d<'a, T: Iterator<Item = &'a u8>, U: SequenceDifferenceModel>(
    pattern: T,
    pattern_length: usize,
    alignment_parameters: &AlignmentParameters<U>,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
) -> Vec<f32> {
    fn index_lookup<T: FMIndexable>(a: u8, l: usize, r: usize, index: &T) -> (usize, usize) {
        let less = index.less(a);
        let l = less + if l > 0 { index.occ(l - 1, a) } else { 0 };
        let r = less + index.occ(r, a) - 1;
        (l, r)
    }

    let r_upper_bound = fmd_index.bwt().len() - 1;

    let (mut l, mut r) = (0, r_upper_bound);
    let mut z = 0.0;

    pattern
        .enumerate()
        .map(|(i, &a)| {
            let tmp = index_lookup(a, l, r, fmd_index);
            l = tmp.0;
            r = tmp.1;

            // Prefix can not be found in the reference
            if l > r {
                // Allow certain transitions
                let penalty = {
                    if a == b'T' {
                        let (l_prime, r_prime) = index_lookup(b'C', l, r, fmd_index);
                        if l_prime <= r_prime {
                            return alignment_parameters.difference_model.get(
                                i,
                                pattern_length,
                                b'C',
                                a,
                            );
                        }
                    } else if a == b'A' {
                        let (l_prime, r_prime) = index_lookup(b'G', l, r, fmd_index);
                        if l_prime <= r_prime {
                            return alignment_parameters.difference_model.get(
                                i,
                                pattern_length,
                                b'G',
                                a,
                            );
                        }
                    }
                    // Return the minimum penalty
                    alignment_parameters
                        .difference_model
                        .get_min_penalty(i, pattern_length, a)
                        .max(alignment_parameters.penalty_gap_open)
                        .max(alignment_parameters.penalty_gap_extend)
                };

                l = 0;
                r = r_upper_bound;
                z -= penalty;
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
        let alphabet = alphabets::dna::n_alphabet();
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

        let alphabet = alphabets::dna::n_alphabet();
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

        // Reverse reference
        let mut reverse_reference = ref_seq.into_iter().rev().collect::<Vec<_>>();
        let rev_data_fmd_index = build_auxiliary_structures(&mut reverse_reference, &alphabet);

        let rev_fm_index = FMIndex::new(
            &rev_data_fmd_index.bwt,
            &rev_data_fmd_index.less,
            &rev_data_fmd_index.occ,
        );
        let rev_fmd_index = FMDIndex::from(rev_fm_index);

        let pattern = "GTTC".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            1.0,
            &parameters,
            &fmd_index,
            &rev_fmd_index,
        );

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

        let alphabet = alphabets::dna::n_alphabet();
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

        // Reverse reference
        let mut reverse_reference = ref_seq.into_iter().rev().collect::<Vec<_>>();
        let rev_data_fmd_index = build_auxiliary_structures(&mut reverse_reference, &alphabet);

        let rev_fm_index = FMIndex::new(
            &rev_data_fmd_index.bwt,
            &rev_data_fmd_index.less,
            &rev_data_fmd_index.occ,
        );
        let rev_fmd_index = FMDIndex::from(rev_fm_index);

        let pattern = "TTTT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            1.0,
            &parameters,
            &fmd_index,
            &rev_fmd_index,
        );

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

        let alphabet = alphabets::dna::n_alphabet();
        let mut ref_seq = "ACGTACGTACGTACGT".as_bytes().to_owned();

        // Reference
        let data_fmd_index = build_auxiliary_structures(&mut ref_seq, &alphabet);

        let fm_index = FMIndex::new(
            &data_fmd_index.bwt,
            &data_fmd_index.less,
            &data_fmd_index.occ,
        );
        let fmd_index = FMDIndex::from(fm_index);

        // Reverse reference
        let mut reverse_reference = ref_seq.into_iter().rev().collect::<Vec<_>>();
        let rev_data_fmd_index = build_auxiliary_structures(&mut reverse_reference, &alphabet);

        let rev_fm_index = FMIndex::new(
            &rev_data_fmd_index.bwt,
            &rev_data_fmd_index.less,
            &rev_data_fmd_index.occ,
        );
        let rev_fmd_index = FMDIndex::from(rev_fm_index);

        let pattern = "GTTC".as_bytes().to_owned();

        let d_backward = calculate_d(pattern.iter(), pattern.len(), &parameters, &rev_fmd_index);
        let d_forward = calculate_d(pattern.iter().rev(), pattern.len(), &parameters, &fmd_index)
            .into_iter()
            .rev()
            .collect::<Vec<_>>();

        assert_eq!(vec![0.0, 0.0, 1.0, 1.0], d_backward);
        assert_eq!(vec![1.0, 1.0, 1.0, 0.0], d_forward);

        let pattern = "GATC".as_bytes().to_owned();

        let d_backward = calculate_d(pattern.iter(), pattern.len(), &parameters, &rev_fmd_index);
        let d_forward = calculate_d(pattern.iter().rev(), pattern.len(), &parameters, &fmd_index)
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

        let alphabet = alphabets::dna::n_alphabet();
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

        // Reverse reference
        let mut reverse_reference = ref_seq.into_iter().rev().collect::<Vec<_>>();
        let rev_data_fmd_index = build_auxiliary_structures(&mut reverse_reference, &alphabet);

        let rev_fm_index = FMIndex::new(
            &rev_data_fmd_index.bwt,
            &rev_data_fmd_index.less,
            &rev_data_fmd_index.occ,
        );
        let rev_fmd_index = FMDIndex::from(rev_fm_index);

        let pattern = "TT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            2.0,
            &parameters,
            &fmd_index,
            &rev_fmd_index,
        );

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&suffix_array))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(vec![0, 0, 0, 0, 2, 2, 5, 5], positions);
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

        let alphabet = alphabets::dna::n_alphabet();
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

        // Reverse reference
        let mut reverse_reference = ref_seq.into_iter().rev().collect::<Vec<_>>();
        let rev_data_fmd_index = build_auxiliary_structures(&mut reverse_reference, &alphabet);

        let rev_fm_index = FMIndex::new(
            &rev_data_fmd_index.bwt,
            &rev_data_fmd_index.less,
            &rev_data_fmd_index.occ,
        );
        let rev_fmd_index = FMDIndex::from(rev_fm_index);

        let pattern = "TTCCCT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            30.0,
            &parameters,
            &fmd_index,
            &rev_fmd_index,
        );

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

        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            30.0,
            &parameters,
            &fmd_index,
            &rev_fmd_index,
        );

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

        // Reverse reference
        let mut reverse_reference = ref_seq.into_iter().rev().collect::<Vec<_>>();
        let rev_data_fmd_index = build_auxiliary_structures(&mut reverse_reference, &alphabet);

        let rev_fm_index = FMIndex::new(
            &rev_data_fmd_index.bwt,
            &rev_data_fmd_index.less,
            &rev_data_fmd_index.occ,
        );
        let rev_fmd_index = FMDIndex::from(rev_fm_index);

        let pattern = "AAGAAA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            30.0,
            &parameters,
            &fmd_index,
            &rev_fmd_index,
        );

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_approx_eq!(-10.969394, alignment_score[0]);
    }
}
