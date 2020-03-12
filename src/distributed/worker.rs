use crate::{
    distributed::{ResultSheet, TaskRxBuffer, TaskSheet},
    map,
    utils::{load_index_from_path, AlignmentParameters},
};
use backtrack_tree::Tree;
use bio::data_structures::{
    bwt::{Less, Occ, BWT},
    fmindex::FMDIndex,
};
use log::debug;
use min_max_heap::MinMaxHeap;
use rayon::prelude::*;
use std::{
    cell::RefCell,
    error::Error,
    io,
    io::{Read, Write},
    net::TcpStream,
};

pub struct Worker {
    network_buffer: TaskRxBuffer,
    connection: TcpStream,
    fmd_index: Option<FMDIndex<BWT, Less, Occ>>,
    alignment_parameters: Option<AlignmentParameters>,
}

impl Worker {
    pub fn new(host: &str, port: &str) -> Result<Self, io::Error> {
        Ok(Self {
            network_buffer: TaskRxBuffer::new(),
            connection: TcpStream::connect(format!("{}:{}", host, port))?,
            fmd_index: None,
            alignment_parameters: None,
        })
    }

    pub fn run(&mut self) -> Result<(), Box<dyn Error>> {
        thread_local! {
            static STACK_BUF: RefCell<MinMaxHeap<map::MismatchSearchStackFrame>> = RefCell::new(MinMaxHeap::with_capacity(map::STACK_LIMIT + 9));
            static TREE_BUF: RefCell<Tree<Option<map::EditOperation>>> = RefCell::new(Tree::with_capacity(map::STACK_LIMIT + 9));
        }

        loop {
            // Receive task
            match self.read_message() {
                Ok(mut task) => {
                    // Load FMD index if necessary
                    if self.fmd_index.is_none() {
                        if let Some(reference_path) = task.reference_path {
                            debug!("Load FMD-index");
                            self.fmd_index = Some(load_index_from_path(&reference_path)?);
                        }
                    }

                    // Extract alignment parameters if necessary
                    if self.alignment_parameters.is_none() {
                        debug!("Extract alignment parameters");
                        if let Some(alignment_parameters) = task.alignment_parameters {
                            self.alignment_parameters = Some(alignment_parameters);
                        }
                    }

                    // If requirements are met, run task
                    if let Some(fmd_index) = &self.fmd_index {
                        if let Some(alignment_parameters) = &self.alignment_parameters {
                            debug!("Map reads");
                            let results = std::mem::replace(&mut task.records, Vec::new())
                                .into_par_iter()
                                .map(|record| {
                                    STACK_BUF.with(|stack_buf| {
                                        TREE_BUF.with(|tree_buf| {
                                            let hit_intervals = map::k_mismatch_search(
                                                &record.sequence,
                                                &record.base_qualities,
                                                alignment_parameters,
                                                fmd_index,
                                                &mut stack_buf.borrow_mut(),
                                                &mut tree_buf.borrow_mut(),
                                            );
                                            (record, hit_intervals)
                                        })
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
    fn read_message(&mut self) -> Result<TaskSheet, io::Error> {
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
