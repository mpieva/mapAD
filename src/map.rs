use std::cmp::min;
use std::error::Error;
use std::fs::File;

use bincode::deserialize_from;
use bio::data_structures::bwt::Occ;
use bio::data_structures::fmindex::{FMDIndex, FMIndex, FMIndexable, Interval};
use bio::io::fastq;
use libflate::deflate::Decoder;

use crate::utils::{AlignmentParameters, AllowedMismatches};
use rust_htslib::bam;

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

        const TEN_F32: f32 = 10.0;
        let mut sum_base_q_best = intervals
            .iter()
            .take(1)
            .map(|alignment| alignment.sum_base_qualities)
            .sum();
        let mut sum_base_q_all = 0.0;
        for sum_of_qualities in intervals.iter() {
            sum_base_q_all += TEN_F32.powi(-(sum_of_qualities.sum_base_qualities));
            if sum_of_qualities.sum_base_qualities < sum_base_q_best {
                sum_base_q_best = sum_of_qualities.sum_base_qualities
            }
        }
        let mapping_quality =
            -((1.0 - (TEN_F32.powi(-(sum_base_q_best)) / sum_base_q_all)).log10());

        dbg!(&record.id());
        dbg!(&mapping_quality);

        // Create SAM/BAM records
        for imm in intervals.iter() {
            for position in imm.interval.occ(suffix_array).iter() {
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct IntervalQuality {
    interval: Interval,
    sum_base_qualities: i32,
}

struct MismatchSearchParameters<'a> {
    pattern: &'a [u8],
    base_qualities: &'a [u8],
    d: &'a [i32],
    parameters: &'a AlignmentParameters,
    fmd_index: &'a FMDIndex<&'a Vec<u8>, &'a Vec<usize>, &'a Occ>,
    j: isize,
    z: i32,
    interval: Interval,
    current_sum_base_qualities: i32,
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
    debug!("Calculate auxiliary array D");
    let d = calculate_d(&pattern, &parameters, rev_fmd_index);

    debug!("Map reads");
    k_mismatch_search_recursive(MismatchSearchParameters {
        pattern,
        base_qualities,
        d: &d,
        parameters,
        fmd_index,
        j: pattern.len() as isize - 1,
        z,
        interval: Interval {
            lower: 0,
            upper: fmd_index.bwt().len() - 1,
        },
        current_sum_base_qualities: 0,
    })
}

/// Follows closely the implementation of BWA-backtrack (Li & Durbin, 2009)
fn k_mismatch_search_recursive(mut par: MismatchSearchParameters) -> Vec<IntervalQuality> {
    // Too many mismatches
    if par.z < par.d[if par.j < 0 { 0 } else { par.j as usize }] {
        return Vec::new();
    }

    // This route through the read graph is finished successfully, return the interval
    if par.j < 0 {
        par.interval.upper += 1;
        let mut interval_set = Vec::new();
        interval_set.push(IntervalQuality {
            interval: par.interval,
            sum_base_qualities: par.current_sum_base_qualities,
        });
        return interval_set;
    }

    let mut interval_set = Vec::new();

    // Insertion in read
    // TODO: Adaptive penalty
    interval_set.extend(&k_mismatch_search_recursive(MismatchSearchParameters {
        j: par.j - 1,
        z: par.z - 1,
        current_sum_base_qualities: par.current_sum_base_qualities
            + i32::from(par.base_qualities[par.j as usize]),
        ..par
    }));

    for &c in b"ACGT".iter() {
        let tmp = index_lookup(c, par.interval.lower, par.interval.upper, par.fmd_index);
        let interval_prime = Interval {
            lower: tmp.0,
            upper: tmp.1,
        };

        if par.interval.lower <= par.interval.upper {
            // Deletion in read
            // TODO: Adaptive penalty
            interval_set.extend(&k_mismatch_search_recursive(MismatchSearchParameters {
                z: par.z - 1,
                interval: interval_prime,
                current_sum_base_qualities: par.current_sum_base_qualities
                    + i32::from(par.base_qualities[par.j as usize]),
                ..par
            }));

            // Match
            if c == par.pattern[par.j as usize] {
                interval_set.extend(&k_mismatch_search_recursive(MismatchSearchParameters {
                    j: par.j - 1,
                    interval: interval_prime,
                    ..par
                }));

            // Mismatch
            } else {
                let penalty = match (c as char, par.pattern[par.j as usize] as char) {
                    ('C', 'T') => par.parameters.penalty_c_t,
                    ('G', 'A') => par.parameters.penalty_g_a,
                    _ => par.parameters.penalty_mismatch,
                };
                interval_set.extend(&k_mismatch_search_recursive(MismatchSearchParameters {
                    j: par.j - 1,
                    z: par.z - penalty,
                    interval: interval_prime,
                    current_sum_base_qualities: par.current_sum_base_qualities
                        + i32::from(par.base_qualities[par.j as usize]),
                    ..par
                }));
            }
        }
    }
    interval_set
}

fn index_lookup<T: FMIndexable>(a: u8, l: usize, r: usize, index: &T) -> (usize, usize) {
    let less = index.less(a);
    let l = less + if l > 0 { index.occ(l - 1, a) } else { 0 };
    let r = less + index.occ(r, a) - 1;
    (l, r)
}

///// A reversed FMD-index is used to compute the lower bound of mismatches of a read per position.
///// This allows for pruning the search tree. Implementation follows closely Li & Durbin (2009).
fn calculate_d(
    pattern: &[u8],
    alignment_parameters: &AlignmentParameters,
    reverse_fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
) -> Vec<i32> {
    let r_upper_bound = reverse_fmd_index.bwt().len() - 1;

    let (mut l, mut r) = (0, r_upper_bound);
    let mut z = 0;

    pattern
        .iter()
        .map(|&a| {
            let tmp = index_lookup(a, l, r, reverse_fmd_index);
            l = tmp.0;
            r = tmp.1;

            // Prefix can not be found in the reference
            if l > r {
                // Allow certain transitions
                let penalty: i32 = {
                    if a == b'T' {
                        let (l_prime, r_prime) = index_lookup(b'C', l, r, reverse_fmd_index);
                        if l_prime <= r_prime {
                            return alignment_parameters.penalty_c_t;
                        }
                    } else if a == b'A' {
                        let (l_prime, r_prime) = index_lookup(b'G', l, r, reverse_fmd_index);
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
        assert_eq!(positions_1, [10, 3]);

        // Test non-occurring pattern
        let pattern_2 = b"GG";
        let sai_2 = fmd_index.backward_search(pattern_2.iter());
        let positions_2 = sai_2.occ(&suffix_array);
        assert_eq!(positions_2, []);
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
        let mut ref_seq = "GAAGAC".as_bytes().to_owned();

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

        let pattern = "CAC".as_bytes().to_owned();
        let base_qualities = [0, 0, 0];

        let d = calculate_d(&pattern, &parameters, &rev_fmd_index);
        assert_eq!(d, vec![0, 1, 1]);

        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            1,
            &parameters,
            &fmd_index,
            &rev_fmd_index,
        );
        let mut positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.interval.occ(&suffix_array))
            .flatten()
            .collect();
        positions.sort();
        assert_eq!(positions, vec![3, 4]);
    }
}
