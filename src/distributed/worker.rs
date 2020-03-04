use crate::{
    distributed::{ResultSheet, TaskRxBuffer, TaskSheet},
    map,
    sequence_difference_models::SequenceDifferenceModel,
    utils::{AlignmentParameters, AllowedMismatches, Record, UnderlyingDataFMDIndex},
};
use log::debug;
use rayon::prelude::*;
use serde::{de::DeserializeOwned, Serialize};
use std::{
    cell::RefCell,
    collections::BinaryHeap,
    error::Error,
    fmt::Debug,
    io,
    io::{Read, Write},
    net::TcpStream,
};

pub struct Worker<T> {
    network_buffer: TaskRxBuffer<T>,
    connection: TcpStream,
    fmd_data: Option<UnderlyingDataFMDIndex>,
    alignment_parameters: Option<AlignmentParameters<T>>,
}

impl<T> Worker<T>
where
    T: SequenceDifferenceModel + Serialize + DeserializeOwned + Sync + Debug,
{
    pub fn new(host: &str, port: &str) -> Result<Self, io::Error> {
        Ok(Self {
            network_buffer: TaskRxBuffer::new(),
            connection: TcpStream::connect(format!("{}:{}", host, port))?,
            fmd_data: None,
            alignment_parameters: None,
        })
    }

    pub fn run(&mut self) -> Result<(), Box<dyn Error>> {
        loop {
            // Receive task
            match self.read_message() {
                Ok(mut task) => {
                    // Load FMD index if necessary
                    if self.fmd_data.is_none() {
                        if let Some(reference_path) = task.reference_path {
                            debug!("Load FMD-index");
                            self.fmd_data = Some(UnderlyingDataFMDIndex::load(&reference_path)?);
                        }
                    }

                    // Extract alignment parameters if necessary
                    if self.alignment_parameters.is_none() {
                        if let Some(alignment_parameters) = task.alignment_parameters {
                            self.alignment_parameters = Some(alignment_parameters);
                        }
                    }

                    // If requirements are met, run task
                    if let Some(data_fmd_index) = &self.fmd_data {
                        if let Some(alignment_parameters) = &self.alignment_parameters {
                            debug!("Reconstruct FMD-index");
                            let fmd_index = data_fmd_index.reconstruct();
                            let allowed_mismatches = AllowedMismatches::new(alignment_parameters);

                            debug!("Map reads");
                            thread_local! {
                                static STACK_BUF: RefCell<BinaryHeap<map::MismatchSearchStackFrame>> = RefCell::new(BinaryHeap::with_capacity(map::STACK_LIMIT + 9))
                            }
                            let results = std::mem::replace(&mut task.records, Vec::new())
                                .into_par_iter()
                                .map(|record: Record| {
                                    let seq_len = record.sequence.len();
                                    let allowed_number_of_mismatches =
                                        allowed_mismatches.get(seq_len);
                                    let representative_mismatch_penalty = alignment_parameters
                                        .difference_model
                                        .get_representative_mismatch_penalty();

                                    STACK_BUF.with(|stack_buf| {
                                        let hit_intervals = map::k_mismatch_search(
                                            &record.sequence,
                                            &record.base_qualities,
                                            allowed_number_of_mismatches
                                                * representative_mismatch_penalty,
                                            alignment_parameters,
                                            &fmd_index,
                                            &mut stack_buf.borrow_mut(),
                                        );

                                        (record, hit_intervals)
                                    })
                                })
                                .collect::<Vec<_>>();

                            // Return results
                            self.connection
                                .write_all(&ResultSheet::new(task.chunk_id, results).encode())?;
                        }
                    }
                }
                Err(ref e) if e.kind() == io::ErrorKind::ConnectionAborted => {
                    debug!("The dispatcher has dropped the connection, shutting down gracefully");
                    return Ok(());
                }
                Err(e) => {
                    return Err(Box::new(e));
                }
            }
        }
    }

    /// Reads task sheet completely from connection in a blocking way
    /// and decodes it eventually
    fn read_message(&mut self) -> Result<TaskSheet<T>, io::Error> {
        // Read and decode first bytes in which the message size is stored
        if self
            .connection
            .read_exact(self.network_buffer.buf_mut_unfilled())
            .is_err()
        {
            // Apparently, the dispatcher has dropped the connection
            return Err(io::Error::new(
                io::ErrorKind::ConnectionAborted,
                "Connection aborted",
            ));
        }

        // Enlarge buffer to fit the entire message
        if self.network_buffer.decode_header().is_err() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Could not decode task header",
            ));
        }

        // Read (blocking) and decode message
        self.connection
            .read_exact(&mut self.network_buffer.buf_mut_unfilled())?;
        match self.network_buffer.decode_and_reset() {
            Ok(task) => Ok(task),
            Err(_) => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Could not decode task message",
            )),
        }
    }
}
