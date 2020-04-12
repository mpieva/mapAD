use crate::{
    fmd_index::RtFMDIndex, mismatch_bounds::MismatchBoundDispatch,
    sequence_difference_models::SequenceDifferenceModelDispatch,
};
use bio::{
    alphabets,
    alphabets::{dna, RankTransform},
    data_structures::{
        bwt::{bwt, less, Less, Occ, BWT},
        suffix_array::{suffix_array, RawSuffixArray},
    },
};
use log::debug;
use rust_htslib::bam;
use serde::{Deserialize, Serialize};
use std::{fmt::Debug, fs::File};

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
pub struct AlignmentParameters {
    pub difference_model: SequenceDifferenceModelDispatch,
    pub mismatch_bound: MismatchBoundDispatch,
    pub penalty_gap_open: f32,
    pub penalty_gap_extend: f32,
    pub chunk_size: usize,
}

pub fn load_index_from_path(path: &str) -> Result<RtFMDIndex, bincode::Error> {
    debug!("Load BWT");
    let bwt: BWT = {
        let d_bwt = snap::read::FrameDecoder::new(File::open(format!("{}.tbw", path))?);
        bincode::deserialize_from(d_bwt)?
    };

    debug!("Load \"C\" table");
    let less: Less = {
        let d_less = snap::read::FrameDecoder::new(File::open(format!("{}.tle", path))?);
        bincode::deserialize_from(d_less)?
    };

    debug!("Load \"Occ\" table");
    let occ: Occ = {
        let d_occ = snap::read::FrameDecoder::new(File::open(format!("{}.toc", path))?);
        bincode::deserialize_from(d_occ)?
    };

    debug!("Load \"RT\" table");
    let rt: RankTransform = {
        let d_rt = snap::read::FrameDecoder::new(File::open(format!("{}.trt", path))?);
        bincode::deserialize_from(d_rt)?
    };

    debug!("Reconstruct index");
    Ok(RtFMDIndex::new(bwt, less, occ, rt))
}

/// This is only used in tests and benchmarks
pub fn build_auxiliary_structures(
    mut reference: Vec<u8>,
    mut src_alphabet: alphabets::Alphabet,
) -> (RtFMDIndex, RawSuffixArray) {
    let ref_seq_revcomp = alphabets::dna::revcomp(reference.iter());
    reference.extend_from_slice(b"$");
    reference.extend_from_slice(&ref_seq_revcomp);
    drop(ref_seq_revcomp);
    reference.extend_from_slice(b"$");

    src_alphabet.insert(b'$');
    let rank_transform = RankTransform::new(&src_alphabet);
    reference = rank_transform.transform(reference);
    let rank_alphabet = alphabets::Alphabet::new(rank_transform.ranks.values());

    let sar = suffix_array(&reference);
    let bwt = bwt(&reference, &sar);
    let less = less(&bwt, &rank_alphabet);
    let occ = Occ::new(&bwt, 3, &rank_alphabet);

    (RtFMDIndex::new(bwt, less, occ, rank_transform), sar)
}
