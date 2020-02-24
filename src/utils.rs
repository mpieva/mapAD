use crate::sequence_difference_models::SequenceDifferenceModel;
use bio::{
    alphabets::dna,
    data_structures::{
        bwt::Occ,
        fmindex::{FMDIndex, FMIndex},
    },
};
use log::debug;
use rust_htslib::bam;
use serde::{Deserialize, Serialize};
use std::{
    fmt::{Display, Error, Formatter},
    fs::File,
    iter::once,
};

/// Auxiliary record data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum BamTag {
    Integer(i64),
    String(Vec<u8>),
    Float(f64),
    Char(u8),
}

impl From<bam::record::Aux<'_>> for BamTag {
    fn from(input: bam::record::Aux) -> Self {
        match input {
            bam::record::Aux::Integer(v) => BamTag::Integer(v),
            bam::record::Aux::String(v) => BamTag::String(v.to_owned()),
            bam::record::Aux::Float(v) => BamTag::Float(v),
            bam::record::Aux::Char(v) => BamTag::Char(v),
        }
    }
}

impl BamTag {
    pub fn borrow_htslib_bam_record(&self) -> bam::record::Aux<'_> {
        match self {
            Self::Integer(v) => bam::record::Aux::Integer(*v),
            Self::String(v) => bam::record::Aux::String(v),
            Self::Float(v) => bam::record::Aux::Float(*v),
            Self::Char(v) => bam::record::Aux::Char(*v),
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Record {
    pub sequence: Vec<u8>,
    pub base_qualities: Vec<u8>,
    pub name: Vec<u8>,
    pub bam_tags: Vec<([u8; 2], BamTag)>,
}

impl From<bam::Record> for Record {
    fn from(input: bam::Record) -> Self {
        let (sequence, base_qualities) = {
            let sequence = input.seq().as_bytes().to_ascii_uppercase();
            let base_qualities = input.qual().to_owned();
            if input.is_reverse() {
                (
                    dna::revcomp(sequence),
                    base_qualities.into_iter().rev().collect(),
                )
            } else {
                (sequence, base_qualities)
            }
        };

        let input_tags = input
            .get_tags()
            .iter()
            .map(|&(tag, value)| ([tag[0], tag[1]], value.into()))
            .collect();

        Self {
            sequence,
            base_qualities,
            name: input.qname().to_owned(),
            bam_tags: input_tags,
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct AlignmentParameters<T> {
    pub base_error_rate: f32,
    pub poisson_threshold: f32,
    pub difference_model: T,
    pub penalty_gap_open: f32,
    pub penalty_gap_extend: f32,
    pub chunk_size: usize,
}

pub struct AllowedMismatches<'a, T> {
    alignment_parameters: &'a AlignmentParameters<T>,
    cache: [f32; 128],
}

impl<'a, T: SequenceDifferenceModel + Sync> Display for AllowedMismatches<'a, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), Error> {
        let max_length = 256;
        let width = (max_length as f32).log10().ceil() as usize;

        let text = {
            let mut tmp = (AllowedMismatches::<T>::MIN_READ_LENGTH..=max_length)
                .map(|read_length| (read_length, self.get(read_length)))
                .scan(
                    std::f32::MIN,
                    |previous, (read_length, allowed_mismatches)| {
                        if (allowed_mismatches - *previous).abs() > std::f32::EPSILON {
                            *previous = allowed_mismatches;
                            Some(Some((read_length, allowed_mismatches)))
                        } else {
                            Some(None)
                        }
                    },
                )
                .flatten()
                .map(|(read_length, allowed_mismatches)| {
                    format!(
                        "{:>width$} bp:\t{} {}\n",
                        read_length,
                        allowed_mismatches,
                        if allowed_mismatches > 1.0 + std::f32::EPSILON {
                            "mismatches"
                        } else {
                            "mismatch"
                        },
                        width = width,
                    )
                })
                .collect::<String>();
            let _ = tmp.pop();
            tmp
        };

        write!(f, "{}", text)
    }
}

impl<'a, T: SequenceDifferenceModel + Sync> AllowedMismatches<'a, T> {
    // 16 is the theoretical minimum for mapping uniquely to a random genome of human-genome-like size
    const MIN_READ_LENGTH: usize = 17;

    pub fn new(alignment_parameters: &AlignmentParameters<T>) -> AllowedMismatches<T> {
        let mut cache = [0.0; 128];
        for (read_length, value) in cache.iter_mut().enumerate() {
            *value = AllowedMismatches::<T>::calculate_max_num_mismatches(
                // Read lengths are stored with an offset to cache the "hot" lengths
                read_length + AllowedMismatches::<T>::MIN_READ_LENGTH,
                alignment_parameters.poisson_threshold,
                alignment_parameters.base_error_rate,
            );
        }

        AllowedMismatches {
            alignment_parameters,
            cache,
        }
    }

    pub fn get(&self, read_length: usize) -> f32 {
        // Reads shorter than the threshold are not allowed to differ
        if read_length < AllowedMismatches::<T>::MIN_READ_LENGTH {
            return 0.0;
        }

        match self
            .cache
            // An offset must be subtracted to point to the correct cache entries
            .get(read_length - AllowedMismatches::<T>::MIN_READ_LENGTH)
        {
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
        poisson_threshold: f32,
        base_error_rate: f32,
    ) -> f32 {
        let lambda = read_length as f32 * base_error_rate;
        let exp_minus_lambda = (-lambda).exp();

        // k = 0 (here 1, because BWA allows k+1 mismatches, and so do we)
        once((1, exp_minus_lambda))
            // k = 1..read_length
            .chain((1..=read_length as u64).scan(
                (1.0, 1, exp_minus_lambda),
                |(lambda_to_the_k, k_factorial, sum), k| {
                    *lambda_to_the_k *= lambda;
                    *k_factorial *= k;
                    *sum += *lambda_to_the_k * exp_minus_lambda / *k_factorial as f32;
                    // BWA allows k+1 mismatches, and so do we
                    Some((k + 1, *sum))
                },
            ))
            .take_while(|(_k, sum)| 1.0 - *sum > poisson_threshold)
            .last()
            .map(|(k, _sum)| k)
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
            let d_bwt = snap::read::FrameDecoder::new(File::open(format!("{}.tbw", path))?);
            bincode::deserialize_from(d_bwt)?
        };

        debug!("Load \"C\" table");
        let less: Vec<usize> = {
            let d_less = snap::read::FrameDecoder::new(File::open(format!("{}.tle", path))?);
            bincode::deserialize_from(d_less)?
        };

        debug!("Load \"Occ\" table");
        let occ: Occ = {
            let d_occ = snap::read::FrameDecoder::new(File::open(format!("{}.toc", path))?);
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
        assert_eq!(2.0, allowed_mismatches.get(17));
        assert_eq!(0.0, allowed_mismatches.get(16));
        assert_eq!(0.0, allowed_mismatches.get(15));
        assert_eq!(0.0, allowed_mismatches.get(3));
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

        assert_eq!(10.0, allowed_mismatches.get(207));
        assert_eq!(9.0, allowed_mismatches.get(176));
        assert_eq!(8.0, allowed_mismatches.get(146));
        assert_eq!(7.0, allowed_mismatches.get(117));
        assert_eq!(6.0, allowed_mismatches.get(90));
        assert_eq!(5.0, allowed_mismatches.get(64));
        assert_eq!(4.0, allowed_mismatches.get(42));
        assert_eq!(3.0, allowed_mismatches.get(22));
        assert_eq!(2.0, allowed_mismatches.get(17));
        assert_eq!(0.0, allowed_mismatches.get(8));
        assert_eq!(0.0, allowed_mismatches.get(1));
    }

    #[test]
    fn test_display() {
        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.06,
            difference_model: VindijaPWM::new(),
            penalty_gap_open: 1.0,
            penalty_gap_extend: 1.0,
            chunk_size: 1,
        };
        let allowed_mismatches = AllowedMismatches::new(&parameters);

        let comparison = " 17 bp:\t1 mismatch
 20 bp:\t2 mismatches
 45 bp:\t3 mismatches
 73 bp:\t4 mismatches
104 bp:\t5 mismatches
137 bp:\t6 mismatches
172 bp:\t7 mismatches
208 bp:\t8 mismatches
244 bp:\t9 mismatches";

        assert_eq!(comparison, format!("{}", allowed_mismatches));
    }
}
