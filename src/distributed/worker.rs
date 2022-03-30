use std::{
    cell::RefCell,
    io::{self, Read, Write},
    mem,
    net::TcpStream,
};

use log::{debug, info};
use min_max_heap::MinMaxHeap;
use rayon::prelude::*;

use crate::{
    backtrack_tree::Tree,
    distributed::{ResultSheet, TaskRxBuffer, TaskSheet},
    errors::{Error, Result},
    fmd_index::RtFmdIndex,
    map,
    mismatch_bounds::MismatchBoundDispatch,
    sequence_difference_models::SequenceDifferenceModelDispatch,
    utils::{load_index_from_path, AlignmentParameters},
};

pub struct Worker {
    network_buffer: TaskRxBuffer,
    connection: TcpStream,
    fmd_index: Option<RtFmdIndex>,
    alignment_parameters: Option<AlignmentParameters>,
}

impl Worker {
    pub fn new(host: &str, port: u16) -> Result<Self> {
        info!("Wait for dispatcher to respond");
        Ok(Self {
            network_buffer: TaskRxBuffer::new(),
            connection: TcpStream::connect((host, port))?,
            fmd_index: None,
            alignment_parameters: None,
        })
    }

    pub fn run(&mut self) -> Result<()> {
        thread_local! {
            static STACK_BUF: RefCell<MinMaxHeap<map::MismatchSearchStackFrame>> = RefCell::new(MinMaxHeap::with_capacity(map::STACK_LIMIT as usize + 9));
            static TREE_BUF: RefCell<Tree<map::EditOperation>> = RefCell::new(Tree::with_capacity(map::EDIT_TREE_LIMIT + 9));
        }

        loop {
            // Receive task
            match self.read_message() {
                Ok(mut task) => {
                    debug!("Received task");
                    // Load FMD index if necessary
                    if self.fmd_index.is_none() {
                        if let Some(reference_path) = task.reference_path {
                            info!("Load FMD-index");
                            self.fmd_index = Some(load_index_from_path(&reference_path)?);
                        }
                    }

                    // Extract alignment parameters if necessary
                    if self.alignment_parameters.is_none() {
                        info!("Extract alignment parameters");
                        if let Some(alignment_parameters) = task.alignment_parameters {
                            self.alignment_parameters = Some(alignment_parameters);
                        }
                    }

                    // If requirements are met, run task
                    if let Some(fmd_index) = &self.fmd_index {
                        if let Some(alignment_parameters) = &self.alignment_parameters {
                            info!("Map reads");
                            let results = mem::take(&mut task.records)
                                .into_par_iter()
                                .map(|record| {
                                    STACK_BUF.with(|stack_buf| {
                                        TREE_BUF.with(|tree_buf| {
                                            // Here we call the instances of the generic `k_mismatch_search` function. Currently, only
                                            // `SimpleAncientDnaModel` is used, as well as `Discrete` and `Continuuous` `MismatchBound`s.
                                            // The others merely serve as an example.
                                            // TODO: I wonder if there's a less-boilerplate way of dispatching.
                                            let hit_intervals = match &alignment_parameters.difference_model {
                                                SequenceDifferenceModelDispatch::SimpleAncientDnaModel(sdm) => {
                                                    match &alignment_parameters.mismatch_bound {
                                                        MismatchBoundDispatch::Discrete(mb) => map::k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                        MismatchBoundDispatch::Continuous(mb) => map::k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                        MismatchBoundDispatch::TestBound(mb) => map::k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                    }
                                                }
                                                SequenceDifferenceModelDispatch::TestDifferenceModel(sdm) => {
                                                    match &alignment_parameters.mismatch_bound {
                                                        MismatchBoundDispatch::Discrete(mb) => map::k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                        MismatchBoundDispatch::Continuous(mb) => map::k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                        MismatchBoundDispatch::TestBound(mb) => map::k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                    }
                                                },
                                                SequenceDifferenceModelDispatch::VindijaPwm(sdm) => {
                                                    match &alignment_parameters.mismatch_bound {
                                                        MismatchBoundDispatch::Discrete(mb) => map::k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                        MismatchBoundDispatch::Continuous(mb) => map::k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                        MismatchBoundDispatch::TestBound(mb) => map::k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                    }
                                                },
                                            };
                                            (record, hit_intervals)})
                                    })
                                })
                                .collect::<Vec<_>>();

                            info!("Return results");
                            self.connection
                                .write_all(&ResultSheet::new(task.chunk_id, results).encode())?;
                        }
                    }
                }
                Err(Error::Io(e)) if e.kind() == io::ErrorKind::ConnectionAborted => {
                    info!("The dispatcher has dropped the connection, shutting down");
                    return Ok(());
                }
                Err(e) => {
                    return Err(e);
                }
            }
        }
    }

    /// Reads task sheet completely from connection in a blocking way
    /// and decodes it eventually
    fn read_message(&mut self) -> Result<TaskSheet> {
        // Read and decode first bytes in which the message size is stored
        self.connection
            .read_exact(self.network_buffer.buf_mut_unfilled())
            .map_err(|_e| io::Error::new(io::ErrorKind::ConnectionAborted, "Connection aborted"))?;

        // Enlarge buffer to fit the entire message
        self.network_buffer.decode_header().map_err(|_e| {
            io::Error::new(io::ErrorKind::InvalidData, "Could not decode task header")
        })?;

        // Read (blocking) and decode message
        self.connection
            .read_exact(self.network_buffer.buf_mut_unfilled())?;
        self.network_buffer.decode_and_reset().map_err(|_e| {
            io::Error::new(io::ErrorKind::InvalidData, "Could not decode task message").into()
        })
    }
}
