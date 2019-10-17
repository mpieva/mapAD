use crate::sequence_difference_models::SequenceDifferenceModel;
use bio::data_structures::{
    bwt::Occ,
    fmindex::{FMDIndex, FMIndex},
};
use log::debug;
use rust_htslib::bam;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use std::fs::File;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Record {
    pub sequence: Vec<u8>,
    pub base_qualities: Vec<u8>,
    pub name: Vec<u8>,
    pub read_group: Option<Vec<u8>>,
}

impl From<bam::Record> for Record {
    fn from(input: bam::Record) -> Self {
        Self {
            sequence: input.seq().as_bytes().to_ascii_uppercase(),
            base_qualities: input.qual().to_owned(),
            name: input.qname().to_owned(),
            read_group: if let Some(rg) = input.aux(b"RG") {
                Some(rg.string().to_owned())
            } else {
                None
            },
        }
    }
}

    pub base_error_rate: f64,
    pub poisson_threshold: f64,
    pub difference_model: T,
    pub penalty_gap_open: f32,
    pub penalty_gap_extend: f32,
    pub chunk_size: usize,
}

pub struct AllowedMismatches<'a, T> {
    alignment_parameters: &'a AlignmentParameters<T>,
    cache: SmallVec<[f32; 128]>,
}

impl<'a, T: SequenceDifferenceModel + Sync> AllowedMismatches<'a, T> {
    pub fn new(alignment_parameters: &AlignmentParameters<T>) -> AllowedMismatches<T> {
        let cache = (0..128)
            .map(|read_length| {
                AllowedMismatches::<T>::calculate_max_num_mismatches(
                    read_length,
                    alignment_parameters.poisson_threshold,
                    alignment_parameters.base_error_rate,
                )
            })
            .collect();

        AllowedMismatches {
            alignment_parameters,
            cache,
        }
    }

    pub fn get(&self, read_length: usize) -> f32 {
        match self.cache.get(read_length) {
            None => AllowedMismatches::<T>::calculate_max_num_mismatches(
                read_length,
                self.alignment_parameters.poisson_threshold,
                self.alignment_parameters.base_error_rate,
            ),
            Some(v) => *v,
        }
    }

    fn calculate_max_num_mismatches(
        read_length: usize,
        poisson_threshold: f64,
        base_error_rate: f64,
    ) -> f32 {
        let lambda = read_length as f64 * base_error_rate;
        let exp_minus_lambda = (-lambda).exp();

        let mut lambda_to_the_power_of_k = 1.0;
        let mut k_factorial = 1;
        let mut sum = exp_minus_lambda;

        // k = 0. BWA allows k+1 mismatches, and so do we
        [(0 + 1, sum)]
            .iter()
            .cloned()
            // k = 1..read_length
            .chain((1..=read_length as u64).map(|k| {
                lambda_to_the_power_of_k *= lambda;
                k_factorial *= k;
                sum += lambda_to_the_power_of_k * exp_minus_lambda / k_factorial as f64;
                // BWA allows k+1 mismatches, and so do we
                (k + 1, sum)
            }))
            .take_while(|(_k, sum)| 1.0 - sum > poisson_threshold)
            .map(|(k, _sum)| k)
            .last()
            .unwrap_or(0) as f32
    }
}

/// Helper struct to bundle index files
pub struct UnderlyingDataFMDIndex {
    bwt: Vec<u8>,
    less: Vec<usize>,
    occ: Occ,
}

impl UnderlyingDataFMDIndex {
    pub fn load(path: &str) -> Result<UnderlyingDataFMDIndex, bincode::Error> {
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

impl UnderlyingDataFMDIndex {
    pub fn new(bwt: Vec<u8>, less: Vec<usize>, occ: Occ) -> Self {
        Self { bwt, less, occ }
    }
    pub fn reconstruct(&self) -> FMDIndex<&Vec<u8>, &Vec<usize>, &Occ> {
        FMDIndex::from(FMIndex::new(&self.bwt, &self.less, &self.occ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence_difference_models::VindijaPWM;

    #[test]
    fn test_allowed_mismatches() {
        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model: VindijaPWM::new(),
            penalty_gap_open: 1.0,
            penalty_gap_extend: 1.0,
            chunk_size: 1,
        };
        let allowed_mismatches = AllowedMismatches::new(&parameters);

        assert_eq!(6.0, allowed_mismatches.get(156));
        assert_eq!(6.0, allowed_mismatches.get(124));
        assert_eq!(5.0, allowed_mismatches.get(123));
        assert_eq!(5.0, allowed_mismatches.get(93));
        assert_eq!(4.0, allowed_mismatches.get(92));
        assert_eq!(4.0, allowed_mismatches.get(64));
        assert_eq!(3.0, allowed_mismatches.get(63));
        assert_eq!(3.0, allowed_mismatches.get(38));
        assert_eq!(2.0, allowed_mismatches.get(37));
        assert_eq!(2.0, allowed_mismatches.get(16));
        assert_eq!(1.0, allowed_mismatches.get(15));
        assert_eq!(1.0, allowed_mismatches.get(3));
        assert_eq!(0.0, allowed_mismatches.get(2));
        assert_eq!(0.0, allowed_mismatches.get(0));
    }

    #[test]
    fn test_allowed_mismatches_bwa_ancient_parameters() {
        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.01,
            difference_model: VindijaPWM::new(),
            penalty_gap_open: 1.0,
            penalty_gap_extend: 1.0,
            chunk_size: 1,
        };
        let allowed_mismatches = AllowedMismatches::new(&parameters);

        //        let mut old_val = 0.0;
        //        let mut count = 0;
        //        for i in 0.. {
        //            let max_diff = allowed_mismatches.get(i);
        //            if max_diff != old_val {
        //                count += 1;
        //                old_val = max_diff;
        //                println!("{}bp reads: max_diff = {}", i, max_diff);
        //            }
        //            if count >= 10 {
        //                break;
        //            }
        //        }

        assert_eq!(10.0, allowed_mismatches.get(207));
        assert_eq!(9.0, allowed_mismatches.get(176));
        assert_eq!(8.0, allowed_mismatches.get(146));
        assert_eq!(7.0, allowed_mismatches.get(117));
        assert_eq!(6.0, allowed_mismatches.get(90));
        assert_eq!(5.0, allowed_mismatches.get(64));
        assert_eq!(4.0, allowed_mismatches.get(42));
        assert_eq!(3.0, allowed_mismatches.get(22));
        assert_eq!(2.0, allowed_mismatches.get(8));
        assert_eq!(1.0, allowed_mismatches.get(1));
    }
}
