use std::cmp::{min, Ordering};
use std::error::Error;
use std::fs::File;

use bio::alphabets::dna;
use bio::data_structures::bwt::Occ;
use bio::data_structures::fmindex::{BiInterval, FMDIndex, FMIndex, FMIndexable};

use bincode::deserialize_from;
use bio::io::fastq;
use libflate::deflate::Decoder;
use rust_htslib::bam;

use crate::utils::{AlignmentParameters, AllowedMismatches};
use std::collections::binary_heap::BinaryHeap;

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
    alignment_score: i32,
}

#[derive(Debug)]
struct MismatchSearchParameters {
    j: isize,
    z: i32,
    interval: BiInterval,
    current_sum_base_qualities: i32,
    backward_pointer: isize,
    forward_pointer: isize,
    forward: bool,
    alignment_score: i32,
    gap: bool,
}

impl PartialOrd for MismatchSearchParameters {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MismatchSearchParameters {
    fn cmp(&self, other: &Self) -> Ordering {
        self.alignment_score.cmp(&other.alignment_score)
    }
}

impl PartialEq for MismatchSearchParameters {
    fn eq(&self, other: &Self) -> bool {
        self.alignment_score == other.alignment_score
    }
}

impl Eq for MismatchSearchParameters {}

pub fn run(reads_path: &str, alignment_parameters: &AlignmentParameters) -> Result<(), Box<Error>> {
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
        &alignment_parameters,
        reads_path,
        &fmd_index,
        &rev_fmd_index,
        &suffix_array,
    )?;

    Ok(())
}

fn map_reads(
    alignment_parameters: &AlignmentParameters,
    reads_path: &str,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
    rev_fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
    suffix_array: &Vec<usize>,
) -> Result<(), Box<Error>> {
    let header = bam::Header::new();
    let mut out = bam::Writer::from_path(&"out.bam", &header).unwrap();

    let mut allowed_mismatches = AllowedMismatches::new(&alignment_parameters);

    let reads_fq_reader = fastq::Reader::from_file(reads_path)?;

    for record in reads_fq_reader.records() {
        debug!("Map read");
        let record = record.unwrap();
        let pattern = record.seq().to_ascii_uppercase();

        // Hardcoded value (33) that should be ok only for Illumina reads
        let base_qualities = record.qual().iter().map(|&f| f - 33).collect::<Vec<_>>();

        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            allowed_mismatches.get(pattern.len()),
            &alignment_parameters,
            &fmd_index,
            &rev_fmd_index,
        );

        debug!("Estimate mapping quality");
        const TEN_F32: f32 = 10.0;
        let mut sum_base_q_best = i32::max_value();
        let mut sum_base_q_all = 0.0;
        for sum_of_qualities in intervals.iter() {
            sum_base_q_all += TEN_F32.powi(-(sum_of_qualities.alignment_score));
            if sum_of_qualities.alignment_score < sum_base_q_best {
                sum_base_q_best = sum_of_qualities.alignment_score
            }
        }
        let mapping_quality =
            -((1.0 - (TEN_F32.powi(-(sum_base_q_best)) / sum_base_q_all)).log10());

        dbg!(&record.id());
        dbg!(&mapping_quality);

        // Create SAM/BAM records
        for imm in intervals.iter() {
            for position in imm.interval.forward().occ(suffix_array).iter() {
                let mut bam_record = bam::record::Record::new();
                bam_record.set_qname(record.id().as_bytes());
                bam_record.set_pos(*position as i32);
                bam_record.set_mapq(mapping_quality as u8);
                out.write(&bam_record)?;
            }
        }
    }
    Ok(())
}

/// Finds all suffix array intervals for the current pattern with up to z mismatches
pub fn k_mismatch_search(
    pattern: &[u8],
    base_qualities: &[u8],
    z: i32,
    parameters: &AlignmentParameters,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
    rev_fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
) -> Vec<IntervalQuality> {
    let center_of_read = pattern.len() as isize / 2;

    let d_backwards = calculate_d(
        pattern[..center_of_read as usize].iter(),
        &parameters,
        rev_fmd_index,
    );
    let d_forwards = calculate_d(
        pattern[center_of_read as usize..].iter().rev(),
        &parameters,
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
        interval: fmd_index.init_interval(),
        current_sum_base_qualities: 0,
        backward_pointer: center_of_read - 1,
        forward_pointer: center_of_read,
        forward: true,
        gap: false,
        alignment_score: 0,
    });

    while let Some(stack_frame) = stack.pop() {
        // Too many mismatches
        let backwards_lower_bound = match d_backwards.get(stack_frame.backward_pointer as usize) {
            Some(&d_i) => d_i,
            None => 0,
        };
        let forwards_lower_bound =
            match d_forwards.get((stack_frame.forward_pointer - center_of_read) as usize) {
                Some(&d_i) => d_i,
                None => 0,
            };
        if stack_frame.z < backwards_lower_bound + forwards_lower_bound {
            continue;
        }

        // This route through the read graph is finished successfully, push the interval
        if stack_frame.j < 0 || stack_frame.j > (pattern.len() as isize - 1) {
            intervals.push(IntervalQuality {
                interval: stack_frame.interval,
                alignment_score: stack_frame.alignment_score,
            });
            continue;
        }

        let next_j;
        let next_backward_pointer;
        let next_forward_pointer;

        // Determine direction of progress for next iteration on this stack frame
        let fmd_ext_interval = if stack_frame.forward {
            next_forward_pointer = stack_frame.forward_pointer + 1;
            next_backward_pointer = stack_frame.backward_pointer;
            next_j = next_backward_pointer;
            stack_frame.interval.swapped()
        } else {
            next_forward_pointer = stack_frame.forward_pointer;
            next_backward_pointer = stack_frame.backward_pointer - 1;
            next_j = next_forward_pointer;
            stack_frame.interval
        };

        // Insertion in read
        let penalty = if stack_frame.gap {
            parameters.penalty_gap_extend
        } else {
            parameters.penalty_gap_open
        };
        stack.push(MismatchSearchParameters {
            j: next_j,
            z: stack_frame.z - penalty,
            current_sum_base_qualities: stack_frame.current_sum_base_qualities
                + i32::from(base_qualities[stack_frame.j as usize]),
            backward_pointer: next_backward_pointer,
            forward_pointer: next_forward_pointer,
            forward: !stack_frame.forward,
            gap: true,
            alignment_score: (stack_frame.forward_pointer as i32
                - stack_frame.backward_pointer as i32
                - (z - (stack_frame.z - penalty))),
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

                // No need to branch for technical characters and zero-intervals
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
            let penalty = if stack_frame.gap {
                parameters.penalty_gap_extend
            } else {
                parameters.penalty_gap_open
            };
            stack.push(MismatchSearchParameters {
                z: stack_frame.z - penalty,
                interval: interval_prime,
                current_sum_base_qualities: stack_frame.current_sum_base_qualities
                    + i32::from(base_qualities[stack_frame.j as usize]),
                alignment_score: (stack_frame.forward_pointer as i32
                    - stack_frame.backward_pointer as i32
                    - (z - (stack_frame.z - penalty))),
                gap: true,
                ..stack_frame
            });

            // Match
            if c == pattern[stack_frame.j as usize] {
                stack.push(MismatchSearchParameters {
                    j: next_j,
                    interval: interval_prime,
                    backward_pointer: next_backward_pointer,
                    forward_pointer: next_forward_pointer,
                    forward: !stack_frame.forward,
                    gap: false,
                    alignment_score: stack_frame.alignment_score + 1,
                    ..stack_frame
                });

            // Mismatch
            } else {
                let penalty = match (c as char, pattern[stack_frame.j as usize] as char) {
                    ('C', 'T') => parameters.penalty_c_t,
                    ('G', 'A') => parameters.penalty_g_a,
                    _ => parameters.penalty_mismatch,
                };
                stack.push(MismatchSearchParameters {
                    j: next_j,
                    z: stack_frame.z - penalty,
                    interval: interval_prime,
                    current_sum_base_qualities: stack_frame.current_sum_base_qualities
                        + i32::from(base_qualities[stack_frame.j as usize]),
                    backward_pointer: next_backward_pointer,
                    forward_pointer: next_forward_pointer,
                    forward: !stack_frame.forward,
                    gap: false,
                    alignment_score: (stack_frame.forward_pointer as i32
                        - stack_frame.backward_pointer as i32
                        - (z - (stack_frame.z - penalty))),
                });
            }
        }
    }
    intervals
}

///// A reversed FMD-index is used to compute the lower bound of mismatches of a read per position.
///// This allows for pruning the search tree. Implementation follows closely Li & Durbin (2009).
fn calculate_d<'a, T: Iterator<Item = &'a u8>>(
    pattern: T,
    alignment_parameters: &AlignmentParameters,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
) -> Vec<i32> {
    fn index_lookup<T: FMIndexable>(a: u8, l: usize, r: usize, index: &T) -> (usize, usize) {
        let less = index.less(a);
        let l = less + if l > 0 { index.occ(l - 1, a) } else { 0 };
        let r = less + index.occ(r, a) - 1;
        (l, r)
    }

    let r_upper_bound = fmd_index.bwt().len() - 1;

    let (mut l, mut r) = (0, r_upper_bound);
    let mut z = 0;

    pattern
        .map(|&a| {
            let tmp = index_lookup(a, l, r, fmd_index);
            l = tmp.0;
            r = tmp.1;

            // Prefix can not be found in the reference
            if l > r {
                // Allow certain transitions
                let penalty: i32 = {
                    if a == b'T' {
                        let (l_prime, r_prime) = index_lookup(b'C', l, r, fmd_index);
                        if l_prime <= r_prime {
                            return alignment_parameters.penalty_c_t;
                        }
                    } else if a == b'A' {
                        let (l_prime, r_prime) = index_lookup(b'G', l, r, fmd_index);
                        if l_prime <= r_prime {
                            return alignment_parameters.penalty_g_a;
                        }
                    }
                    min(
                        alignment_parameters.penalty_mismatch,
                        alignment_parameters.penalty_gap_open,
                    )
                    .min(alignment_parameters.penalty_gap_extend)
                };

                l = 0;
                r = r_upper_bound;
                z += penalty;
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
        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            penalty_mismatch: 1,
            penalty_gap_open: 1,
            penalty_gap_extend: 1,
            penalty_c_t: 0,
            penalty_g_a: 0,
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
            1,
            &parameters,
            &fmd_index,
            &rev_fmd_index,
        );

        let alignment_score: Vec<i32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_eq!(vec![3], alignment_score);

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.forward().occ(&suffix_array))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(vec![2, 6, 10, 19, 23, 27], positions);
    }

    #[test]
    fn test_d() {
        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            penalty_mismatch: 1,
            penalty_gap_open: 1,
            penalty_gap_extend: 1,
            penalty_c_t: 0,
            penalty_g_a: 0,
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

        let d_backward = calculate_d(pattern.iter(), &parameters, &rev_fmd_index);
        let d_forward = calculate_d(pattern.iter().rev(), &parameters, &fmd_index)
            .into_iter()
            .rev()
            .collect::<Vec<i32>>();

        assert_eq!(vec![0, 0, 1, 1], d_backward);
        assert_eq!(vec![1, 1, 1, 0], d_forward);

        let pattern = "GATC".as_bytes().to_owned();

        let d_backward = calculate_d(pattern.iter(), &parameters, &rev_fmd_index);
        let d_forward = calculate_d(pattern.iter().rev(), &parameters, &fmd_index)
            .into_iter()
            .rev()
            .collect::<Vec<i32>>();

        assert_eq!(vec![0, 1, 1, 2], d_backward);
        assert_eq!(vec![2, 1, 1, 0], d_forward);
    }
}
