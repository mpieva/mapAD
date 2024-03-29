use std::{
    cell::RefCell,
    io::{self, Read, Write},
    net::TcpStream,
    time::Instant,
};

use log::{debug, info};
use min_max_heap::MinMaxHeap;
use rayon::prelude::*;

use crate::{
    distributed::{comm_buffers::TaskRxBuffer, Message, ResultSheet, TaskSheet},
    errors::{Error, Result},
    index::load_index_from_path,
    map::{
        backtrack_tree::Tree,
        fmd_index::RtFmdIndex,
        mapping::{k_mismatch_search, EDIT_TREE_LIMIT, STACK_LIMIT},
        mismatch_bounds::MismatchBoundDispatch,
        record::EditOperation,
        sequence_difference_models::SequenceDifferenceModelDispatch,
        AlignmentParameters, MismatchSearchStackFrame,
    },
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
            static STACK_BUF: RefCell<MinMaxHeap<MismatchSearchStackFrame>> = RefCell::new(MinMaxHeap::with_capacity(STACK_LIMIT as usize + 9));
            static TREE_BUF: RefCell<Tree<EditOperation>> = RefCell::new(Tree::with_capacity(EDIT_TREE_LIMIT + 9));
        }

        loop {
            // Receive task
            match self.read_message() {
                Ok(task) => {
                    debug!("Received task");
                    // Load FMD index if necessary
                    if self.fmd_index.is_none() {
                        info!("Load FMD-index");
                        self.fmd_index = Some(load_index_from_path(
                            task.get_reference_path()
                                .expect("This is not expected to fail"),
                        )?);
                    }

                    // Extract alignment parameters if necessary
                    if self.alignment_parameters.is_none() {
                        info!("Extract alignment parameters");
                        self.alignment_parameters = Some(
                            task.get_alignment_parameters()
                                .expect("This is not expected to fail")
                                .clone(),
                        );
                    }

                    // If requirements are met, run task
                    if let Some(fmd_index) = &self.fmd_index {
                        if let Some(alignment_parameters) = &self.alignment_parameters {
                            info!("Map reads");
                            let task_id = task.get_id();
                            let results = task.get_records()
                                .into_par_iter()
                                .map(|record| {
                                    STACK_BUF.with(|stack_buf| {
                                        TREE_BUF.with(|tree_buf| {
                                            let start = Instant::now();
                                            // Here we call the instances of the generic `k_mismatch_search` function. Currently, only
                                            // `SimpleAncientDnaModel` is used, as well as `Discrete` and `Continuuous` `MismatchBound`s.
                                            // The others merely serve as an example.
                                            // TODO: I wonder if there's a less-boilerplate way of dispatching.
                                            let hit_intervals = match &alignment_parameters.difference_model {
                                                SequenceDifferenceModelDispatch::SimpleAncientDnaModel(sdm) => {
                                                    match &alignment_parameters.mismatch_bound {
                                                        MismatchBoundDispatch::Discrete(mb) => k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                        MismatchBoundDispatch::Continuous(mb) => k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                        MismatchBoundDispatch::TestBound(mb) => k_mismatch_search(
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
                                                        MismatchBoundDispatch::Discrete(mb) => k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                        MismatchBoundDispatch::Continuous(mb) => k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                        MismatchBoundDispatch::TestBound(mb) => k_mismatch_search(
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
                                                        MismatchBoundDispatch::Discrete(mb) => k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                        MismatchBoundDispatch::Continuous(mb) => k_mismatch_search(
                                                            &record.sequence,
                                                            &record.base_qualities,
                                                            alignment_parameters,
                                                            fmd_index,
                                                            &mut stack_buf.borrow_mut(),
                                                            &mut tree_buf.borrow_mut(),
                                                            sdm,
                                                            mb,
                                                        ),
                                                        MismatchBoundDispatch::TestBound(mb) => k_mismatch_search(
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
                                            let duration = start.elapsed();
                                            (record, hit_intervals, duration)})
                                    })
                                })
                                .collect::<Vec<_>>();

                            info!("Return results");
                            self.connection
                                .write_all(&ResultSheet::new(task_id, results).encode())?;
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
