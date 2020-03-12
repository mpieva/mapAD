use crate::{
    mismatch_bounds::MismatchBoundDispatch,
    sequence_difference_models::SequenceDifferenceModelDispatch,
};
use bio::{
    alphabets::dna,
    data_structures::{
        bwt::{Less, Occ, BWT},
        fmindex::{FMDIndex, FMIndex},
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

pub fn load_index_from_path(path: &str) -> Result<FMDIndex<BWT, Less, Occ>, bincode::Error> {
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

    debug!("Reconstruct index");
    Ok(FMDIndex::from(FMIndex::new(bwt, less, occ)))
}
