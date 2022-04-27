use std::{
    fmt::{Display, Formatter},
    iter::Map,
};

use log::{error, info};
use serde::{Deserialize, Serialize};

use crate::{
    distributed::Message,
    errors::{Error, Result},
    map::{record::Record, AlignmentParameters},
};

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
    pub fn requery_task(&mut self, task_sheet: TaskSheet) {
        self.requeried_tasks.push(task_sheet);
    }
}

/// Convertible to TaskQueue
pub trait IntoTaskQueue<E, I, O, T>
where
    E: Into<Error>,
    I: Into<Record>,
    T: Iterator<Item = std::result::Result<I, E>>,
    O: Iterator<Item = Result<Record>>,
{
    fn into_tasks(self, chunk_size: usize) -> TaskQueue<O>;
}

/// Adds ChunkIterator conversion method to every compatible Iterator. So when new input file types
/// are implemented, it is sufficient to impl `From<T> for Record` for the additional item.
#[allow(clippy::type_complexity)]
impl<E, I, T> IntoTaskQueue<E, I, Map<T, fn(std::result::Result<I, E>) -> Result<Record>>, T> for T
where
    I: Into<Record>,
    E: Into<Error>,
    T: Iterator<Item = std::result::Result<I, E>>,
{
    fn into_tasks(
        self,
        chunk_size: usize,
    ) -> TaskQueue<Map<T, fn(std::result::Result<I, E>) -> Result<Record>>> {
        TaskQueue {
            chunk_id: 0,
            chunk_size,
            records: self.map(|inner| inner.map(|v| v.into()).map_err(|e| e.into())),
            requeried_tasks: Vec::new(),
        }
    }
}

/// Task wrapped in a Struct for sending to a worker
#[derive(Serialize, Deserialize, Debug)]
pub struct TaskSheet {
    encoded_size: u64,
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
            encoded_size: 0,
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
