use std::{collections::BTreeMap, fmt};

use bio::alphabets::dna;
use bstr::{BString, ByteSlice};
use either::Either;
use log::warn;
use noodles::{cram, fastq, sam};
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;

use crate::{
    errors::{Error, Result},
    index::OriginalSymbols,
    map::{
        backtrack_tree::{NodeId, Tree},
        Direction,
    },
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
    String(BString),
    HexByteArray(BString),
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
impl From<&sam::alignment::record::data::field::Value<'_>> for BamAuxField {
    fn from(value: &sam::alignment::record::data::field::Value<'_>) -> Self {
        use sam::alignment::record::data::field::value::{Array, Value};
        match value {
            Value::Character(v) => Self::Char(*v),
            Value::Int8(v) => Self::I8(*v),
            Value::UInt8(v) => Self::U8(*v),
            Value::Int16(v) => Self::I16(*v),
            Value::UInt16(v) => Self::U16(*v),
            Value::Int32(v) => Self::I32(*v),
            Value::UInt32(v) => Self::U32(*v),
            Value::Float(v) => Self::Float(*v),
            //Value::Double(v) => Ok(Self::Double(*v)),
            Value::String(v) => Self::String((*v).to_owned()),
            Value::Hex(v) => Self::HexByteArray((*v).to_owned()),
            Value::Array(Array::Int8(v)) => Self::ArrayI8(
                v.iter()
                    .collect::<std::io::Result<_>>()
                    .unwrap_or_else(|_e| {
                        warn!("Dropped unreadable array");
                        Vec::new()
                    }),
            ),
            Value::Array(Array::UInt8(v)) => Self::ArrayU8(
                v.iter()
                    .collect::<std::io::Result<_>>()
                    .unwrap_or_else(|_e| {
                        warn!("Dropped unreadable array");
                        Vec::new()
                    }),
            ),
            Value::Array(Array::Int16(v)) => Self::ArrayI16(
                v.iter()
                    .collect::<std::io::Result<_>>()
                    .unwrap_or_else(|_e| {
                        warn!("Dropped unreadable array");
                        Vec::new()
                    }),
            ),
            Value::Array(Array::UInt16(v)) => Self::ArrayU16(
                v.iter()
                    .collect::<std::io::Result<_>>()
                    .unwrap_or_else(|_e| {
                        warn!("Dropped unreadable array");
                        Vec::new()
                    }),
            ),
            Value::Array(Array::Int32(v)) => Self::ArrayI32(
                v.iter()
                    .collect::<std::io::Result<_>>()
                    .unwrap_or_else(|_e| {
                        warn!("Dropped unreadable array");
                        Vec::new()
                    }),
            ),
            Value::Array(Array::UInt32(v)) => Self::ArrayU32(
                v.iter()
                    .collect::<std::io::Result<_>>()
                    .unwrap_or_else(|_e| {
                        warn!("Dropped unreadable array");
                        Vec::new()
                    }),
            ),
            Value::Array(Array::Float(v)) => Self::ArrayFloat(
                v.iter()
                    .collect::<std::io::Result<_>>()
                    .unwrap_or_else(|_e| {
                        warn!("Dropped unreadable array");
                        Vec::new()
                    }),
            ),
        }
    }
}

impl From<&sam::alignment::record_buf::data::field::Value> for BamAuxField {
    fn from(value: &sam::alignment::record_buf::data::field::Value) -> Self {
        value.into()
    }
}

impl From<BamAuxField> for sam::alignment::record_buf::data::field::Value {
    fn from(input: BamAuxField) -> Self {
        use sam::alignment::record_buf::data::field::value::Array;
        match input {
            BamAuxField::Char(v) => Self::Character(v),
            BamAuxField::I8(v) => Self::Int8(v),
            BamAuxField::U8(v) => Self::UInt8(v),
            BamAuxField::I16(v) => Self::Int16(v),
            BamAuxField::U16(v) => Self::UInt16(v),
            BamAuxField::I32(v) => Self::Int32(v),
            BamAuxField::U32(v) => Self::UInt32(v),
            BamAuxField::Float(v) => Self::Float(v),
            BamAuxField::Double(v) => {
                warn!("Lost precision: Converted 8-byte to 4-byte floating point number array");
                Self::Float(v as f32)
            }
            BamAuxField::String(v) => Self::String(v.as_bstr().to_owned()),
            BamAuxField::HexByteArray(v) => Self::Hex(v.as_bstr().to_owned()),
            BamAuxField::ArrayI8(v) => Self::Array(Array::Int8(v)),
            BamAuxField::ArrayU8(v) => Self::Array(Array::UInt8(v)),
            BamAuxField::ArrayI16(v) => Self::Array(Array::Int16(v)),
            BamAuxField::ArrayU16(v) => Self::Array(Array::UInt16(v)),
            BamAuxField::ArrayI32(v) => Self::Array(Array::Int32(v)),
            BamAuxField::ArrayU32(v) => Self::Array(Array::UInt32(v)),
            BamAuxField::ArrayFloat(v) => Self::Array(Array::Float(v)),
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

impl TryFrom<&dyn sam::alignment::Record> for Record {
    type Error = Error;

    fn try_from(input: &dyn sam::alignment::Record) -> Result<Self> {
        let mut sequence = input.sequence().iter().collect::<Vec<_>>();

        // Reads can not be longer than `i16::MAX`
        if i16::try_from(sequence.len()).is_err() {
            return Err(Error::SeqLenError(input.name().map_or_else(
                || String::from("unnamed record"),
                |maybe_record| String::from_utf8_lossy(maybe_record.as_bytes()).into_owned(),
            )));
        }

        let mut base_qualities = input.quality_scores().as_ref().iter().collect::<Vec<_>>();

        if let Ok(flags) = input.flags() {
            if flags.is_reverse_complemented() {
                base_qualities.reverse();
                sequence = dna::revcomp(sequence);
            }
        } else {
            warn!("Dropped unreadable flags");
        };

        let input_tags = input
            .data()
            .iter()
            .filter_map(|maybe_tv| {
                maybe_tv
                    .ok()
                    .as_ref()
                    .map(|(tag, value)| (tag.as_ref().to_owned(), value.into()))
            })
            .collect::<Vec<_>>();

        let read_name = input.name().map(|name| name.as_bytes().to_owned());

        Ok(Self {
            sequence,
            base_qualities,
            name: read_name,
            bam_tags: input_tags,
            bam_flags: input.flags().map(|flags| flags.bits()).unwrap_or(0),
        })
    }
}

impl TryFrom<fastq::Record> for Record {
    type Error = Error;

    fn try_from(fq_record: fastq::Record) -> Result<Self> {
        let sequence = fq_record.sequence().to_ascii_uppercase();

        // Reads can not be longer than `i16::MAX`
        if i16::try_from(sequence.len()).is_err() {
            return Err(Error::SeqLenError(
                String::from_utf8_lossy(fq_record.name()).into_owned(),
            ));
        }

        // Subtract offset
        let base_qualities = fq_record
            .quality_scores()
            .iter()
            .map(|qual| qual - 33)
            .collect();
        let name = fq_record.name().to_owned();

        Ok(Self {
            sequence,
            base_qualities,
            name: Some(name),
            bam_tags: Vec::new(),
            bam_flags: 0,
        })
    }
}

impl TryFrom<cram::record::Record> for Record {
    type Error = Error;

    fn try_from(input: cram::Record) -> Result<Self> {
        use sam::alignment::record::{Name, Sequence};

        let mut sequence = input.sequence().iter().collect::<Vec<_>>();

        // Reads can not be longer than `i16::MAX`
        if i16::try_from(sequence.len()).is_err() {
            return Err(Error::SeqLenError(input.name().map_or_else(
                || String::from("unnamed record"),
                |maybe_record| String::from_utf8_lossy(maybe_record.as_bytes()).into_owned(),
            )));
        }

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
            .iter()
            .map(|(tag, value)| (tag.as_ref().to_owned(), value.into()))
            .collect::<Vec<_>>();

        let read_name = input.name().map(|name| name.as_bytes().to_owned());

        Ok(Self {
            sequence,
            base_qualities,
            name: read_name,
            bam_tags: input_tags,
            bam_flags: input.flags().bits(),
        })
    }
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let read_name = self.name.as_ref().map_or(b"*".as_slice(), AsRef::as_ref);
        write!(f, "{}", String::from_utf8_lossy(read_name))
    }
}

/// Variants store position in the read and, if necessary, the reference base
#[derive(Debug, Copy, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum EditOperation {
    Insertion(u16),
    Deletion(u16, u8),
    Match(u16),
    Mismatch(u16, u8),
}

impl Default for EditOperation {
    fn default() -> Self {
        Self::Match(0)
    }
}

impl From<EditOperation> for sam::alignment::record::cigar::op::Kind {
    fn from(src: EditOperation) -> Self {
        match src {
            EditOperation::Insertion(_) => Self::Insertion,
            EditOperation::Deletion(_, _) => Self::Deletion,
            EditOperation::Match(_) | EditOperation::Mismatch(_, _) => Self::Match,
        }
    }
}

impl From<EditOperation> for sam::alignment::record::cigar::Op {
    fn from(src: EditOperation) -> Self {
        use sam::alignment::record::cigar::op::Kind;
        match src {
            EditOperation::Insertion(l) => Self::new(Kind::Insertion, l.into()),
            EditOperation::Deletion(l, _) => Self::new(Kind::Deletion, l.into()),
            EditOperation::Match(l) | EditOperation::Mismatch(l, _) => {
                Self::new(Kind::Match, l.into())
            }
        }
    }
}

/// Contains edit operations performed in order to align the sequence
#[derive(Debug, Serialize, Deserialize)]
pub struct EditOperationsTrack(Vec<EditOperation>);

impl EditOperationsTrack {
    /// Calculates the amount of positions in the genome
    /// that are covered by this read
    pub fn effective_len(&self) -> usize {
        self.0.iter().fold(0, |acc, edit_operation| {
            acc + match edit_operation {
                EditOperation::Insertion(_) => 0,
                EditOperation::Deletion(_, _) => 1,
                EditOperation::Match(_) => 1,
                EditOperation::Mismatch(_, _) => 1,
            }
        })
    }

    /// Constructs CIGAR, MD tag, and edit distance from correctly ordered track of edit operations and yields them as a tuple
    /// The strand a read is mapped to is taken into account here.
    pub fn to_bam_fields(
        &self,
        strand: Direction,
        absolute_pos: usize,
        original_symbols: &OriginalSymbols,
    ) -> (Vec<sam::alignment::record::cigar::Op>, Vec<u8>, u16) {
        use sam::alignment::record::cigar::Op;
        // Reconstruct the order of the remaining edit operations and condense CIGAR
        let mut num_matches: u32 = 0;
        let mut num_operations = 1;
        let mut edit_distance = 0;
        let mut last_edit_operation = None;
        let mut cigar = Vec::new();
        let mut md_tag = Vec::new();

        let track = match strand {
            Direction::Forward => Either::Left(self.0.iter()),
            Direction::Backward => Either::Right(self.0.iter().rev()),
        };
        for (i, edit_operation) in track.enumerate() {
            let edit_operation = match edit_operation {
                EditOperation::Insertion(_) => *edit_operation,
                EditOperation::Match(j) => original_symbols
                    .get(absolute_pos + i)
                    .map_or(*edit_operation, |original_symbol| {
                        EditOperation::Mismatch(*j, original_symbol)
                    }),
                EditOperation::Deletion(j, reference_base) => EditOperation::Deletion(
                    *j,
                    original_symbols
                        .get(absolute_pos + i)
                        .unwrap_or(*reference_base),
                ),
                EditOperation::Mismatch(j, reference_base) => EditOperation::Mismatch(
                    *j,
                    original_symbols
                        .get(absolute_pos + i)
                        .unwrap_or(*reference_base),
                ),
            };

            edit_distance = Self::add_edit_distance(edit_operation, edit_distance);

            num_matches = Self::add_md_edit_operation(
                Some(edit_operation),
                last_edit_operation,
                strand,
                num_matches,
                &mut md_tag,
            );

            if let Some(lop) = last_edit_operation {
                // Construct CIGAR string
                match edit_operation {
                    EditOperation::Match(_) => match lop {
                        EditOperation::Match(_) | EditOperation::Mismatch(_, _) => {
                            num_operations += 1;
                        }
                        _ => {
                            cigar.push(Op::new(lop.into(), num_operations));
                            num_operations = 1;
                            last_edit_operation = Some(edit_operation);
                        }
                    },
                    EditOperation::Mismatch(_, _) => match lop {
                        EditOperation::Mismatch(_, _) | EditOperation::Match(_) => {
                            num_operations += 1;
                        }
                        _ => {
                            cigar.push(Op::new(lop.into(), num_operations));
                            num_operations = 1;
                            last_edit_operation = Some(edit_operation);
                        }
                    },
                    EditOperation::Insertion(_) => match lop {
                        EditOperation::Insertion(_) => {
                            num_operations += 1;
                        }
                        _ => {
                            cigar.push(Op::new(lop.into(), num_operations));
                            num_operations = 1;
                            last_edit_operation = Some(edit_operation);
                        }
                    },
                    EditOperation::Deletion(_, _) => match lop {
                        EditOperation::Deletion(_, _) => {
                            num_operations += 1;
                        }
                        _ => {
                            cigar.push(Op::new(lop.into(), num_operations));
                            num_operations = 1;
                            last_edit_operation = Some(edit_operation);
                        }
                    },
                }
            } else {
                last_edit_operation = Some(edit_operation);
            }
        }

        // Add remainder
        if let Some(lop) = last_edit_operation {
            cigar.push(Op::new(lop.into(), num_operations));
        }
        let _ = Self::add_md_edit_operation(None, None, strand, num_matches, &mut md_tag);

        (cigar, md_tag, edit_distance)
    }

    fn add_md_edit_operation(
        edit_operation: Option<EditOperation>,
        last_edit_operation: Option<EditOperation>,
        strand: Direction,
        mut k: u32,
        md_tag: &mut Vec<u8>,
    ) -> u32 {
        let comp_if_necessary = |reference_base| match strand {
            Direction::Forward => reference_base,
            Direction::Backward => dna::complement(reference_base),
        };

        match edit_operation {
            Some(EditOperation::Match(_)) => k += 1,
            Some(EditOperation::Mismatch(_, reference_base)) => {
                let reference_base = comp_if_necessary(reference_base);
                md_tag.extend_from_slice(format!("{}{}", k, reference_base as char).as_bytes());
                k = 0;
            }
            Some(EditOperation::Insertion(_)) => {
                // Insertions are ignored in MD tags
            }
            Some(EditOperation::Deletion(_, reference_base)) => {
                let reference_base = comp_if_necessary(reference_base);
                match last_edit_operation {
                    Some(EditOperation::Deletion(_, _)) => {
                        md_tag.extend_from_slice(format!("{}", reference_base as char).as_bytes());
                    }
                    _ => {
                        md_tag.extend_from_slice(
                            format!("{}^{}", k, reference_base as char).as_bytes(),
                        );
                    }
                }
                k = 0;
            }
            None => md_tag.extend_from_slice(format!("{k}").as_bytes()),
        }
        k
    }

    fn add_edit_distance(edit_operation: EditOperation, distance: u16) -> u16 {
        if let EditOperation::Match(_) = edit_operation {
            distance
        } else {
            distance + 1
        }
    }

    pub fn read_len(&self) -> usize {
        self.0.iter().fold(0, |acc, edit_operation| {
            acc + match edit_operation {
                EditOperation::Insertion(_) => 1,
                EditOperation::Deletion(_, _) => 0,
                EditOperation::Match(_) => 1,
                EditOperation::Mismatch(_, _) => 1,
            }
        })
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

/// Derive Cigar string from oddly-ordered tracks of edit operations.
/// Since we start aligning at the center of a read, tracks of edit operations
/// are not ordered by position in the read. Also, the track of edit operations
/// must be extracted from the (possibly huge) tree which is built during
/// backtracking for size reasons.
pub fn extract_edit_operations(
    end_node: NodeId,
    edit_tree: &Tree<EditOperation>,
    pattern_len: usize,
) -> EditOperationsTrack {
    // Restore outer ordering of the edit operation by the positions they carry as values.
    // Whenever there are deletions in the query, there is no simple rule to reconstruct the ordering.
    // So, edit operations carrying the same position are pushed onto the same bucket and dealt with later.
    let mut cigar_order_outer: BTreeMap<u16, SmallVec<[EditOperation; 8]>> = BTreeMap::new();

    for &edit_operation in edit_tree.ancestors(end_node) {
        let position = match edit_operation {
            EditOperation::Insertion(position) => position,
            EditOperation::Deletion(position, _) => position,
            EditOperation::Match(position) => position,
            EditOperation::Mismatch(position, _) => position,
        };
        cigar_order_outer
            .entry(position)
            .or_default()
            .push(edit_operation);
    }

    EditOperationsTrack(
        cigar_order_outer
            .into_iter()
            .flat_map(|(i, inner_vec)| {
                if i < (pattern_len / 2) as u16 {
                    Either::Left(inner_vec.into_iter())
                } else {
                    Either::Right(inner_vec.into_iter().rev())
                }
            })
            .collect(),
    )
}

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn test_edop_effective_len() {
        let edop_track = EditOperationsTrack(vec![
            EditOperation::Match(0),
            EditOperation::Mismatch(1, b'C'),
            EditOperation::Match(2),
            EditOperation::Insertion(3),
            EditOperation::Match(4),
            EditOperation::Deletion(5, b'A'),
            EditOperation::Deletion(6, b'G'),
            EditOperation::Match(7),
            EditOperation::Match(8),
            EditOperation::Match(9),
            EditOperation::Match(10),
            EditOperation::Insertion(11),
            EditOperation::Mismatch(10, b'C'),
        ]);
        assert_eq!(edop_track.effective_len(), 11);

        let edop_track_2 = EditOperationsTrack(vec![
            EditOperation::Insertion(0),
            EditOperation::Insertion(1),
            EditOperation::Insertion(2),
        ]);
        assert_eq!(edop_track_2.effective_len(), 0);

        let edop_track_3 = EditOperationsTrack(vec![
            EditOperation::Deletion(0, b'A'),
            EditOperation::Deletion(1, b'C'),
            EditOperation::Deletion(2, b'G'),
            EditOperation::Deletion(3, b'T'),
        ]);
        assert_eq!(edop_track_3.effective_len(), 4);
    }
}
