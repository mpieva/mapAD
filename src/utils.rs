use std::fmt::{self, Debug};

use bio::{
    alphabets,
    alphabets::{dna, RankTransform},
    data_structures::{
        bwt::{bwt, less, Occ},
        suffix_array::{suffix_array, RawSuffixArray},
    },
    io::fastq,
};
use log::debug;
use noodles::{
    bam,
    sam::{self, AlignmentRecord},
};
use serde::{Deserialize, Serialize};

use crate::{
    errors::Result,
    fmd_index::RtFmdIndex,
    index::{SampledSuffixArrayOwned, VersionedIndexItem},
    map::FastaIdPositions,
    mismatch_bounds::MismatchBoundDispatch,
    sequence_difference_models::SequenceDifferenceModelDispatch,
};

/// An owned representation of the `bam::record::Aux` data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum BamAuxField {
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

// We (currently) only get references to the internal fields of `noodles::bam::Record`s,
// so we have to copy/clone data over
impl From<&sam::record::data::field::Value> for BamAuxField {
    fn from(input: &sam::record::data::field::Value) -> Self {
        use sam::record::data::field::Value;
        match input {
            Value::Char(v) => BamAuxField::Char(*v as u8),
            Value::Int8(v) => BamAuxField::I8(*v),
            Value::UInt8(v) => BamAuxField::U8(*v),
            Value::Int16(v) => BamAuxField::I16(*v),
            Value::UInt16(v) => BamAuxField::U16(*v),
            Value::Int32(v) => BamAuxField::I32(*v),
            Value::UInt32(v) => BamAuxField::U32(*v),
            Value::Float(v) => BamAuxField::Float(*v),
            //Value::Double(v) => BamTag::Double(*v),
            Value::String(v) => BamAuxField::String(v.to_owned()),
            Value::Hex(v) => BamAuxField::HexByteArray(v.to_owned()),
            Value::Int8Array(v) => BamAuxField::ArrayI8(v.to_owned()),
            Value::UInt8Array(v) => BamAuxField::ArrayU8(v.to_owned()),
            Value::Int16Array(v) => BamAuxField::ArrayI16(v.to_owned()),
            Value::UInt16Array(v) => BamAuxField::ArrayU16(v.to_owned()),
            Value::Int32Array(v) => BamAuxField::ArrayI32(v.to_owned()),
            Value::UInt32Array(v) => BamAuxField::ArrayU32(v.to_owned()),
            Value::FloatArray(v) => BamAuxField::ArrayFloat(v.to_owned()),
        }
    }
}

impl From<BamAuxField> for sam::record::data::field::Value {
    fn from(input: BamAuxField) -> Self {
        use sam::record::data::field::Value;
        match input {
            BamAuxField::Char(v) => Value::Char(v.into()),
            BamAuxField::I8(v) => Value::Int8(v),
            BamAuxField::U8(v) => Value::UInt8(v),
            BamAuxField::I16(v) => Value::Int16(v),
            BamAuxField::U16(v) => Value::UInt16(v),
            BamAuxField::I32(v) => Value::Int32(v),
            BamAuxField::U32(v) => Value::UInt32(v),
            BamAuxField::Float(v) => Value::Float(v),
            BamAuxField::Double(v) => Value::Float(v as f32), // FIXME
            BamAuxField::String(v) => Value::String(v),
            BamAuxField::HexByteArray(v) => Value::Hex(v),
            BamAuxField::ArrayI8(v) => Value::Int8Array(v),
            BamAuxField::ArrayU8(v) => Value::UInt8Array(v),
            BamAuxField::ArrayI16(v) => Value::Int16Array(v),
            BamAuxField::ArrayU16(v) => Value::UInt16Array(v),
            BamAuxField::ArrayI32(v) => Value::Int32Array(v),
            BamAuxField::ArrayU32(v) => Value::UInt32Array(v),
            BamAuxField::ArrayFloat(v) => Value::FloatArray(v),
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Record {
    pub sequence: Vec<u8>,
    pub base_qualities: Vec<u8>,
    pub name: Option<Vec<u8>>,
    pub bam_tags: Vec<([u8; 2], BamAuxField)>,
    pub bam_flags: u16,
}

impl From<bam::Record> for Record {
    fn from(input: bam::Record) -> Self {
        let mut sequence = input.sequence().to_string().into_bytes();

        let mut base_qualities = input
            .quality_scores()
            .as_ref()
            .iter()
            .copied()
            .map(u8::from)
            .collect::<Vec<_>>();

        if input.flags().is_reverse_complemented() {
            base_qualities.reverse();
            sequence = dna::revcomp(sequence);
        };

        let input_tags = input
            .data()
            .values()
            .map(|field| (field.tag().as_ref().to_owned(), field.value().into()))
            .collect::<Vec<_>>();

        let read_name = input.read_name().map(|name| {
            let name: &[u8] = name.as_ref();
            name.to_vec()
        });

        Self {
            sequence,
            base_qualities,
            name: read_name,
            bam_tags: input_tags,
            bam_flags: input.flags().bits(),
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
            name: Some(name),
            bam_tags: Vec::new(),
            bam_flags: 0,
        }
    }
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let read_name = self
            .name
            .as_ref()
            .map(|name| name.as_ref())
            .unwrap_or(sam::record::read_name::MISSING);
        write!(f, "{}", String::from_utf8_lossy(read_name))
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
    VersionedIndexItem::read_from_path(format!("{}.tsa", reference_path))?.try_take()
}

pub fn load_id_pos_map_from_path(reference_path: &str) -> Result<FastaIdPositions> {
    VersionedIndexItem::read_from_path(format!("{}.tpi", reference_path))?.try_take()
}

pub fn load_index_from_path(reference_path: &str) -> Result<RtFmdIndex> {
    debug!("Load BWT");
    let bwt = VersionedIndexItem::read_from_path(format!("{}.tbw", reference_path))?.try_take()?;

    debug!("Load \"C\" table");
    let less = VersionedIndexItem::read_from_path(format!("{}.tle", reference_path))?.try_take()?;

    debug!("Load \"Occ\" table");
    let occ = VersionedIndexItem::read_from_path(format!("{}.toc", reference_path))?.try_take()?;

    debug!("Load \"RT\" table");
    let rt = VersionedIndexItem::read_from_path(format!("{}.trt", reference_path))?.try_take()?;

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
