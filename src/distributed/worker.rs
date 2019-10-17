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
                    dbg!(&task);

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
                            let out_buffer = std::mem::replace(&mut task.records, Vec::new())
                                .into_par_iter()
                                .map(|record: Record| {
                                    let hit_intervals = map::k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        allowed_mismatches.get(record.sequence.len())
                                            * alignment_parameters
                                                .difference_model
                                                .get_representative_mismatch_penalty(),
                                        alignment_parameters,
                                        &fmd_index,
                                    );

                                    (record, hit_intervals)
                                })
                                .collect::<Vec<_>>();

                            let mut results = ResultSheet::new(task.chunk_id, out_buffer);

                            // Return results
                            let message = results.encode();
                            self.connection.write_all(&message)?;
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
