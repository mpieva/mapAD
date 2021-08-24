use crate::{
    errors::{Error, Result},
    fmd_index::RtFmdIndex,
    index::{
        VersionedBwt, VersionedIdPosMap, VersionedLess, VersionedOcc, VersionedRt,
        VersionedSuffixArray, INDEX_VERSION,
    },
    map::FastaIdPositions,
    mismatch_bounds::MismatchBoundDispatch,
    sequence_difference_models::SequenceDifferenceModelDispatch,
};

use bio::{
    alphabets,
    alphabets::{dna, RankTransform},
    data_structures::{
        bwt::{bwt, less, Occ, BWT},
        suffix_array::{suffix_array, RawSuffixArray},
    },
    io::fastq,
};
use log::{debug, warn};
use rust_htslib::{bam, bam::record::Aux};
use serde::{Deserialize, Serialize};

use crate::index::SampledSuffixArrayOwned;
use std::{fmt::Debug, fs::File};

/// An owned representation of the `bam::record::Aux` data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum BamTag {
    Char(u8),
    I8(i8),
    U8(u8),
    I16(i16),
    U16(u16),
    I32(i32),
    U32(u32),
    Float(f32),
    Double(f64),
    String(String),
    HexByteArray(String),
    ArrayI8(Vec<i8>),
    ArrayU8(Vec<u8>),
    ArrayI16(Vec<i16>),
    ArrayU16(Vec<u16>),
    ArrayI32(Vec<i32>),
    ArrayU32(Vec<u32>),
    ArrayFloat(Vec<f32>),
}

impl From<Aux<'_>> for BamTag {
    fn from(input: Aux) -> Self {
        match input {
            Aux::Char(v) => BamTag::Char(v),
            Aux::I8(v) => BamTag::I8(v),
            Aux::U8(v) => BamTag::U8(v),
            Aux::I16(v) => BamTag::I16(v),
            Aux::U16(v) => BamTag::U16(v),
            Aux::I32(v) => BamTag::I32(v),
            Aux::U32(v) => BamTag::U32(v),
            Aux::Float(v) => BamTag::Float(v),
            Aux::Double(v) => BamTag::Double(v),
            Aux::String(v) => BamTag::String(v.to_owned()),
            Aux::HexByteArray(v) => BamTag::HexByteArray(v.to_owned()),
            Aux::ArrayI8(v) => BamTag::ArrayI8(v.iter().collect()),
            Aux::ArrayU8(v) => BamTag::ArrayU8(v.iter().collect()),
            Aux::ArrayI16(v) => BamTag::ArrayI16(v.iter().collect()),
            Aux::ArrayU16(v) => BamTag::ArrayU16(v.iter().collect()),
            Aux::ArrayI32(v) => BamTag::ArrayI32(v.iter().collect()),
            Aux::ArrayU32(v) => BamTag::ArrayU32(v.iter().collect()),
            Aux::ArrayFloat(v) => BamTag::ArrayFloat(v.iter().collect()),
        }
    }
}

impl<'a> From<&'a BamTag> for Aux<'a> {
    fn from(input: &'a BamTag) -> Self {
        match input {
            BamTag::Char(v) => Aux::Char(*v),
            BamTag::I8(v) => Aux::I8(*v),
            BamTag::U8(v) => Aux::U8(*v),
            BamTag::I16(v) => Aux::I16(*v),
            BamTag::U16(v) => Aux::U16(*v),
            BamTag::I32(v) => Aux::I32(*v),
            BamTag::U32(v) => Aux::U32(*v),
            BamTag::Float(v) => Aux::Float(*v),
            BamTag::Double(v) => Aux::Double(*v),
            BamTag::String(v) => Aux::String(v),
            BamTag::HexByteArray(v) => Aux::HexByteArray(v),
            BamTag::ArrayI8(v) => Aux::ArrayI8(v.into()),
            BamTag::ArrayU8(v) => Aux::ArrayU8(v.into()),
            BamTag::ArrayI16(v) => Aux::ArrayI16(v.into()),
            BamTag::ArrayU16(v) => Aux::ArrayU16(v.into()),
            BamTag::ArrayI32(v) => Aux::ArrayI32(v.into()),
            BamTag::ArrayU32(v) => Aux::ArrayU32(v.into()),
            BamTag::ArrayFloat(v) => Aux::ArrayFloat(v.into()),
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Record {
    pub sequence: Vec<u8>,
    pub base_qualities: Vec<u8>,
    pub name: Vec<u8>,
    pub bam_tags: Vec<([u8; 2], BamTag)>,
    pub bam_flags: u16,
}

impl From<bam::Record> for Record {
    fn from(input: bam::Record) -> Self {
        let (sequence, base_qualities) = {
            let sequence = input.seq().as_bytes().to_ascii_uppercase();
            // No need to subtract offsets here
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
            .aux_iter()
            .map(|maybe_tag| {
                maybe_tag
                    .map(|(tag, value)| ([tag[0], tag[1]], value.into()))
                    .map_err(|e| {
                        warn!(
                            "Error reading auxiliary data of record {}. Auxiliary data will be incomplete.",
                            String::from_utf8_lossy(input.qname())
                        );
                        e.into()
                    })
            })
            .filter_map(Result::ok)
            .collect::<Vec<_>>();

        Self {
            sequence,
            base_qualities,
            name: input.qname().to_owned(),
            bam_tags: input_tags,
            bam_flags: input.flags(),
        }
    }
}

impl From<fastq::Record> for Record {
    fn from(fq_record: fastq::Record) -> Self {
        let sequence = fq_record.seq().to_ascii_uppercase();
        // Subtract offset
        let base_qualities = fq_record.qual().iter().map(|qual| qual - 33).collect();
        let name = fq_record.id().as_bytes().to_owned();

        Self {
            sequence,
            base_qualities,
            name,
            bam_tags: Vec::new(),
            bam_flags: 0,
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
    pub gap_dist_ends: u8,
    pub stack_limit_abort: bool,
}

pub fn load_suffix_array_from_path(reference_path: &str) -> Result<SampledSuffixArrayOwned> {
    let d_suffix_array =
        snap::read::FrameDecoder::new(File::open(format!("{}.tsa", reference_path))?);
    let versioned_suffix_array: VersionedSuffixArray = bincode::deserialize_from(d_suffix_array)?;
    if versioned_suffix_array.version == INDEX_VERSION {
        Ok(versioned_suffix_array.data)
    } else {
        Err(Error::IndexVersionMismatch)
    }
}

pub fn load_id_pos_map_from_path(reference_path: &str) -> Result<FastaIdPositions> {
    let d_pi = snap::read::FrameDecoder::new(File::open(format!("{}.tpi", reference_path))?);
    let versioned_identifier_position_map: VersionedIdPosMap = bincode::deserialize_from(d_pi)?;
    if versioned_identifier_position_map.version == INDEX_VERSION {
        Ok(versioned_identifier_position_map.data)
    } else {
        Err(Error::IndexVersionMismatch)
    }
}

pub fn load_index_from_path(reference_path: &str) -> Result<RtFmdIndex> {
    debug!("Load BWT");
    let bwt: BWT = {
        let d_bwt = snap::read::FrameDecoder::new(File::open(format!("{}.tbw", reference_path))?);
        let versioned_bwt: VersionedBwt = bincode::deserialize_from(d_bwt)?;
        if versioned_bwt.version == INDEX_VERSION {
            versioned_bwt.data
        } else {
            return Err(Error::IndexVersionMismatch);
        }
    };

    debug!("Load \"C\" table");
    let less = {
        let d_less = snap::read::FrameDecoder::new(File::open(format!("{}.tle", reference_path))?);
        let versioned_less: VersionedLess = bincode::deserialize_from(d_less)?;
        if versioned_less.version == INDEX_VERSION {
            versioned_less.data
        } else {
            return Err(Error::IndexVersionMismatch);
        }
    };

    debug!("Load \"Occ\" table");
    let occ = {
        let d_occ = snap::read::FrameDecoder::new(File::open(format!("{}.toc", reference_path))?);
        let versioned_occ: VersionedOcc = bincode::deserialize_from(d_occ)?;
        if versioned_occ.version == INDEX_VERSION {
            versioned_occ.data
        } else {
            return Err(Error::IndexVersionMismatch);
        }
    };

    debug!("Load \"RT\" table");
    let rt = {
        let d_rt = snap::read::FrameDecoder::new(File::open(format!("{}.trt", reference_path))?);
        let versioned_rt: VersionedRt = bincode::deserialize_from(d_rt)?;
        if versioned_rt.version == INDEX_VERSION {
            versioned_rt.data
        } else {
            return Err(Error::IndexVersionMismatch);
        }
    };

    debug!("Reconstruct index");
    Ok(RtFmdIndex::new(bwt, less, occ, rt))
}

/// This is only used in tests and benchmarks
pub fn build_auxiliary_structures(
    mut reference: Vec<u8>,
    mut src_alphabet: alphabets::Alphabet,
) -> (RtFmdIndex, RawSuffixArray) {
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

    (RtFmdIndex::new(bwt, less, occ, rank_transform), sar)
}
