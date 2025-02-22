use std::{
    fmt::{Display, Formatter},
    fs::File,
    io::{stdin, BufRead, BufReader, Read},
    path::Path,
};

use flate2::read::MultiGzDecoder;
use log::{debug, error, info};
use noodles::{
    bam, bgzf, cram, fastq,
    sam::{self, alignment::io::Read as AlignmentRead},
};
use serde::{Deserialize, Serialize};

use crate::{
    distributed::Message,
    errors::{Error, Result},
    map::{record::Record, AlignmentParameters},
};

pub enum InputSource {
    Bam(
        bam::io::reader::Reader<bgzf::Reader<Box<dyn Read>>>,
        Box<sam::Header>,
    ),
    Cram(cram::io::reader::Reader<Box<dyn Read>>, Box<sam::Header>),
    //Sam(sam::reader::Reader<BR>, Box<sam::Header>),
    // .fastq and fastq.gz
    Fastq(fastq::reader::Reader<Box<dyn BufRead>>),
}

enum Format {
    Bam,
    Cram,
    //Sam,
    Fastq,
    FastqGz,
}

impl InputSource {
    pub fn from_path<P>(path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let path = path.as_ref();

        // Input file options: `.bam`, `.cram`, `.fastq.gz`, and `.fastq` as fallback

        // Use `BufRead` function to peek into the file/stdin, then convert to `Read` trait object.
        // Copying the buffer is expensive, but only done once.

        // Stdin
        let (file_handle, magic_bytes): (Box<dyn Read>, _) = if path.to_str() == Some("-") {
            // BGZIP blocks are <= 64kb, so in order to be able to peek into the file and read the
            // magic bytes we might need a buffer that big
            let mut stdin_guard = BufReader::with_capacity(65536, stdin().lock());
            let buf_copy = stdin_guard.fill_buf()?.to_owned();
            (Box::new(stdin_guard), buf_copy)
        // File
        } else {
            // BGZIP blocks are <= 64kb, so in order to be able to peek into the file and read the
            // magic bytes we might need a buffer that big
            let buf_copy = BufReader::with_capacity(65536, File::open(path)?)
                .fill_buf()?
                .to_owned();
            (Box::new(File::open(path)?), buf_copy)
        };

        match Self::detect_format(&magic_bytes)? {
            Format::Bam => {
                debug!("Try reading input in BAM format");
                let mut reader = bam::io::Reader::new(file_handle);
                let header = Box::new(
                    reader
                        .read_alignment_header()
                        .map_err(Into::<Error>::into)?,
                );
                Ok(Self::Bam(reader, header))
            }
            Format::Cram => {
                debug!("Try reading input in CRAM format");
                let mut reader = cram::io::Reader::new(file_handle);
                let header = Box::new(
                    reader
                        .read_alignment_header()
                        .map_err(Into::<Error>::into)?,
                );
                Ok(Self::Cram(reader, header))
            }
            Format::Fastq => {
                debug!("Try reading input in FASTQ format");
                let reader: fastq::Reader<Box<dyn BufRead>> =
                    fastq::Reader::new(Box::new(BufReader::new(file_handle)));
                Ok(Self::Fastq(reader))
            }
            Format::FastqGz => {
                debug!("Try reading input in gzip compressed FASTQ format");
                let reader: fastq::Reader<Box<dyn BufRead>> =
                    fastq::Reader::new(Box::new(BufReader::new(MultiGzDecoder::new(file_handle))));
                Ok(Self::Fastq(reader))
            }
        }
    }

    // Inspired by `noodles-utils`
    fn detect_format(src: &[u8]) -> Result<Format> {
        const CRAM_MAGIC_NUMBER: [u8; 4] = [b'C', b'R', b'A', b'M'];
        const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];
        const BAM_MAGIC_NUMBER: [u8; 4] = [b'B', b'A', b'M', 0x01];

        let buf = src.get(..4).ok_or(Error::InvalidInputType)?;
        if buf == CRAM_MAGIC_NUMBER {
            return Ok(Format::Cram);
        }

        if buf[..2] == GZIP_MAGIC_NUMBER {
            let mut reader = bgzf::Reader::new(src);
            let mut buf = [0; 4];

            // Read exactly four bytes to see whether they contain the BAM magic numbers
            if reader.read_exact(&mut buf).is_err() {
                // We are going to ignore this error, since we probably got GZ instead of BGZF data.
                // If that is not the case, potential errors are properly handled somewhere else.
            }

            return if buf == BAM_MAGIC_NUMBER {
                Ok(Format::Bam)
            } else {
                Ok(Format::FastqGz)
            };
        }

        Ok(Format::Fastq)
    }

    pub fn header(&self) -> Option<&sam::Header> {
        match self {
            Self::Bam(_, header) | Self::Cram(_, header) => Some(header),
            Self::Fastq(_) => None,
        }
    }

    pub fn task_queue(
        &mut self,
        chunk_size: usize,
    ) -> TaskQueue<Box<dyn Iterator<Item = Result<Record>> + '_>> {
        match self {
            Self::Bam(reader, header) => TaskQueue::new(
                Box::new(reader.alignment_records(header).map(|maybe_record| {
                    maybe_record
                        .map_err(Into::into)
                        .and_then(|record| record.as_ref().try_into())
                })),
                chunk_size,
            ),
            Self::Cram(reader, header) => TaskQueue::new(
                Box::new(reader.alignment_records(header).map(|maybe_record| {
                    maybe_record
                        .map_err(Into::into)
                        .and_then(|record| record.as_ref().try_into())
                })),
                chunk_size,
            ),
            Self::Fastq(reader) => TaskQueue::new(
                Box::new(reader.records().map(|maybe_record| {
                    maybe_record.map_err(Into::into).and_then(TryInto::try_into)
                })),
                chunk_size,
            ),
        }
    }
}

/// Keeps track of the processing state of chunks of reads
///
/// Very basic error checking, reporting, and recovery happens here.
pub struct TaskQueue<T> {
    chunk_id: usize,
    chunk_size: usize,
    records: T,
    requeried_tasks: Vec<TaskSheet>,
}

impl<T> Iterator for TaskQueue<T>
where
    T: Iterator<Item = Result<Record>>,
{
    type Item = TaskSheet;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(task) = self.requeried_tasks.pop() {
            info!("Retrying previously failed task {}", task);
            return Some(task);
        }

        let chunk = self
            .records
            .by_ref()
            .filter_map(|record| {
                if let Err(ref e) = record {
                    error!("Skip record due to an error: {}", e);
                }
                record.ok()
            })
            .filter(|record| {
                let length_check_ok = record.sequence.len() == record.base_qualities.len();
                if !length_check_ok {
                    error!(
                        "Skip record \"{}\" due to different length of sequence and quality strings", record
                    );
                }
                length_check_ok
            })
            .take(self.chunk_size)
            .collect::<Vec<_>>();
        self.chunk_id += 1;

        if chunk.is_empty() {
            None
        } else {
            Some(TaskSheet::from_records(self.chunk_id - 1, chunk))
        }
    }
}

impl<T> TaskQueue<T>
where
    T: Iterator<Item = Result<Record>>,
{
    pub fn new(inner: T, chunk_size: usize) -> Self {
        Self {
            chunk_id: 0,
            chunk_size,
            records: inner,
            requeried_tasks: Vec::new(),
        }
    }

    pub fn requery_task(&mut self, task_sheet: TaskSheet) {
        self.requeried_tasks.push(task_sheet);
    }
}

/// Task wrapped in a Struct for sending to a worker
#[derive(Serialize, Deserialize, Debug)]
pub struct TaskSheet {
    encoded_size: u64, // must be of type `u64` and the first element
    chunk_id: usize,
    records: Vec<Record>,
    reference_path: Option<String>,
    alignment_parameters: Option<AlignmentParameters>,
}

impl Display for TaskSheet {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.chunk_id)
    }
}

impl TaskSheet {
    pub fn from_records(read_set: usize, records: Vec<Record>) -> Self {
        Self {
            encoded_size: Default::default(),
            chunk_id: read_set,
            records,
            reference_path: None,
            alignment_parameters: None,
        }
    }

    pub fn get_reference_path(&self) -> Option<&str> {
        self.reference_path.as_deref()
    }

    pub fn get_alignment_parameters(&self) -> Option<&AlignmentParameters> {
        self.alignment_parameters.as_ref()
    }

    pub fn get_id(&self) -> usize {
        self.chunk_id
    }

    pub fn set_reference_path<T>(&mut self, reference_path: T)
    where
        T: Into<String>,
    {
        self.reference_path = Some(reference_path.into());
    }

    pub fn set_alignment_parameters(&mut self, alignment_parameters: &AlignmentParameters) {
        self.alignment_parameters = Some(alignment_parameters.clone());
    }

    pub fn get_records(self) -> Vec<Record> {
        self.records
    }
}

impl Message for TaskSheet {
    fn encode(&mut self) -> Vec<u8> {
        let encoded_size = bincode::serialized_size(&self).expect("This is not expected to fail");
        self.encoded_size = encoded_size;
        bincode::serialize(&self).expect("This is not expected to fail")
    }
}
