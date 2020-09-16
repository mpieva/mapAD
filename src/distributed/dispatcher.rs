use crate::{distributed::*, map, utils::*};

use bio::data_structures::suffix_array::SuffixArray;
use log::debug;
use mio::*;
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read as BamRead};

use std::{
    collections::{BinaryHeap, HashMap},
    error::Error,
    fs::File,
    io,
    io::{ErrorKind::WouldBlock, Read, Write},
    iter::Peekable,
    net::{IpAddr, Ipv4Addr, SocketAddr},
    path::Path,
};

enum TransportState<E>
where
    E: Error,
{
    // The message has been successfully transferred
    Finished,
    // This operation would block
    Stalled,
    // An error has occurred
    Error(E),
    // There are no tasks left in the queue
    Complete,
}

/// Keeps track of the processing state of chunks of reads
struct TaskQueue<'a> {
    chunk_id: usize,
    chunk_size: usize,
    records: Peekable<bam::Records<'a, bam::Reader>>,
    requeried_tasks: Vec<TaskSheet>,
}

impl<'a> Iterator for TaskQueue<'a> {
    type Item = TaskSheet;

    fn next(&mut self) -> Option<TaskSheet> {
        let source_iterator_loan = &mut self.records;

        // If the underlying iterator is exhausted, resend chunks that were previously assigned to now disconnected workers
        if source_iterator_loan.peek().is_some() {
            self.chunk_id += 1;
            let chunk = source_iterator_loan
                .take(self.chunk_size)
                .map(|record| match record {
                    Ok(record) => Ok(record.into()),
                    Err(e) => Err(e),
                })
                .collect::<Result<Vec<_>, _>>()
                .expect("Input file is corrupt. Cancelled process");

            Some(TaskSheet::from_records(self.chunk_id - 1, chunk, None))
        } else if let Some(task) = self.requeried_tasks.pop() {
            debug!("Retrying previously failed task {}", task.chunk_id);
            Some(task)
        } else {
            None
        }
    }
}

impl<'a> TaskQueue<'a> {
    fn from_reader(reader: &'a mut bam::Reader, chunk_size: usize) -> Self {
        Self {
            chunk_id: 0,
            chunk_size,
            records: reader.records().peekable(),
            requeried_tasks: Vec::new(),
        }
    }

    fn requery_task(&mut self, task_sheet: TaskSheet) {
        self.requeried_tasks.push(task_sheet);
    }
}

// Central place to store data relevant for one dispatcher <-> worker connection
struct Connection {
    pub stream: net::TcpStream,
    pub assigned_task: Option<TaskSheet>,
    pub rx_buffer: ResultRxBuffer,
    pub tx_buffer: TaskTxBuffer,
}

pub struct Dispatcher<'a, 'b> {
    reads_path: &'b str,
    reference_path: &'b str,
    out_file_path: &'b str,
    alignment_parameters: &'a AlignmentParameters,
    connections: HashMap<Token, Connection>,
    accept_connections: bool,
}

impl<'a, 'b> Dispatcher<'a, 'b> {
    // Token numbering starts with '1' because of weird behaviour with Token(0) on BSDs
    const DISPATCHER_TOKEN: Token = Token(1);

    pub fn new(
        reads_path: &'b str,
        reference_path: &'b str,
        out_file_path: &'b str,
        alignment_parameters: &'a AlignmentParameters,
    ) -> Result<Self, bam::Error> {
        Ok(Self {
            reads_path,
            reference_path,
            out_file_path,
            alignment_parameters,
            connections: HashMap::new(),
            accept_connections: true,
        })
    }

    pub fn run(&mut self, port: u16) -> Result<(), Box<dyn Error>> {
        // Set up networking
        let mut max_token = Self::DISPATCHER_TOKEN.0;
        let addr = SocketAddr::new(IpAddr::V4(Ipv4Addr::new(0, 0, 0, 0)), port);
        let mut listener = net::TcpListener::bind(addr)?;

        // Create a poll instance
        let mut poll = Poll::new()?;

        // Start listening for incoming connection attempts and store connections
        poll.registry()
            .register(&mut listener, Self::DISPATCHER_TOKEN, Interest::READABLE)?;

        debug!("Load position map");
        let identifier_position_map: map::FastaIdPositions = {
            let d_pi =
                snap::read::FrameDecoder::new(File::open(format!("{}.tpi", &self.reference_path))?);
            bincode::deserialize_from(d_pi)?
        };

        // Set up input and output files
        let mut bam_reader = bam::Reader::from_path(self.reads_path)?;
        let _ = bam_reader.set_threads(4);
        let header = map::create_bam_header(&bam_reader, &identifier_position_map);
        let mut task_queue =
            TaskQueue::from_reader(&mut bam_reader, self.alignment_parameters.chunk_size);
        let out_file_path = Path::new(&self.out_file_path);
        if out_file_path.exists() {
            return Err(Box::new(io::Error::new(
                io::ErrorKind::AlreadyExists,
                "The given output file already exists",
            )));
        }
        let mut out_file = bam::Writer::from_path(out_file_path, &header, bam::Format::BAM)?;
        let _ = out_file.set_threads(4);

        debug!("Load suffix array");
        let suffix_array: Vec<usize> = {
            let d_suffix_array =
                snap::read::FrameDecoder::new(File::open(format!("{}.tsa", &self.reference_path))?);
            bincode::deserialize_from(d_suffix_array)?
        };

        // Create storage for events
        let mut events = Events::with_capacity(1024);

        debug!("Ready to distribute work");
        loop {
            poll.poll(&mut events, None)?;
            for event in events.iter() {
                match event.token() {
                    // Workers of the world, register!
                    Self::DISPATCHER_TOKEN => {
                        if self.accept_connections {
                            while let Ok((mut remote_stream, remote_addr)) = listener.accept() {
                                debug!("Connection established ({:?})", remote_addr);
                                max_token += 1;
                                let remote_token = Token(max_token);
                                poll.registry().register(
                                    &mut remote_stream,
                                    remote_token,
                                    Interest::READABLE | Interest::WRITABLE,
                                )?;
                                let connection = Connection {
                                    stream: remote_stream,
                                    assigned_task: None,
                                    rx_buffer: ResultRxBuffer::new(),
                                    tx_buffer: TaskTxBuffer::new(),
                                };
                                self.connections.insert(remote_token, connection);
                            }
                        } else {
                            debug!("Task queue is empty: declined connection attempt");
                        }
                    }

                    //
                    // Communication with existing workers
                    //
                    _ => {
                        // Receive results from workers
                        if event.is_readable() {
                            // After finishing previous tasks, the worker is ready to receive fresh jobs
                            match self.read_rx_buffer(event.token()) {
                                TransportState::Finished => {
                                    if let Ok(results) = self
                                        .connections
                                        .get_mut(&event.token())
                                        .expect("This is not expected to fail")
                                        .rx_buffer
                                        .decode_and_reset()
                                    {
                                        debug!(
                                            "Worker {} sent results of task {}",
                                            event.token().0,
                                            results.chunk_id,
                                        );
                                        self.write_results(
                                            results.results,
                                            &suffix_array,
                                            &identifier_position_map,
                                            &mut out_file,
                                        )?;

                                        // Remove completed task from assignments
                                        self.connections
                                            .get_mut(&event.token())
                                            .expect("This is not expected to fail")
                                            .assigned_task = None;
                                    }

                                    poll.registry().reregister(
                                        &mut self
                                            .connections
                                            .get_mut(&event.token())
                                            .expect("This is not expected to fail")
                                            .stream,
                                        event.token(),
                                        Interest::READABLE | Interest::WRITABLE,
                                    )?;
                                }
                                TransportState::Stalled => {
                                    // We're patient!!1!
                                }
                                TransportState::Error(_) => {
                                    debug!("Connection is no longer valid, removed worker {} from pool", event.token().0);
                                    self.release_worker(event.token(), &mut task_queue);
                                }
                                TransportState::Complete => {
                                    // TransportState::Complete means the input file has been processed completely
                                    // which workers don't know. That's why this match arm here is
                                    unreachable!();
                                }
                            }
                        }

                        // Distribution of work
                        if event.is_writable() {
                            match self.write_tx_buffer(event.token(), &mut task_queue) {
                                TransportState::Finished => {
                                    poll.registry().reregister(
                                        &mut self
                                            .connections
                                            .get_mut(&event.token())
                                            .expect("This is not expected to fail")
                                            .stream,
                                        event.token(),
                                        Interest::READABLE,
                                    )?;
                                }
                                TransportState::Stalled => {
                                    // We're patient!!1!
                                }
                                TransportState::Error(_) => {
                                    debug!("Connection is no longer valid, removed worker {} from pool", event.token().0);
                                    self.release_worker(event.token(), &mut task_queue);
                                }
                                TransportState::Complete => {
                                    self.accept_connections = false;
                                    self.release_worker(event.token(), &mut task_queue);
                                    if self.connections.is_empty() {
                                        debug!("All tasks have been completed, shutting down gracefully");
                                        return Ok(());
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fn write_results<S>(
        &self,
        mut hits: Vec<(Record, BinaryHeap<map::HitInterval>)>,
        suffix_array: &S,
        identifier_position_map: &map::FastaIdPositions,
        out_file: &mut bam::Writer,
    ) -> Result<(), bam::Error>
    where
        S: SuffixArray + Sync,
    {
        let bam_records = hits
            .par_iter_mut()
            .map_init(rand::thread_rng, |mut rng, (record, hit_interval)| {
                map::intervals_to_bam(
                    record,
                    hit_interval,
                    suffix_array,
                    identifier_position_map,
                    None,
                    &mut rng,
                )
            })
            .flatten()
            .collect::<Vec<_>>();

        bam_records
            .into_iter()
            .map(|record| out_file.write(&record))
            .collect()
    }

    /// Removes worker_to_be_removed from the event queue and terminates the connection.
    /// This causes the worker_to_be_removed process to terminate.
    /// Also, this worker's task will be requeried.
    fn release_worker(&mut self, worker_to_be_removed: Token, task_queue: &mut TaskQueue) {
        if let Some(assigned_task) = self
            .connections
            .get_mut(&worker_to_be_removed)
            .expect("This is not expected to fail")
            .assigned_task
            .take()
        {
            debug!("Requeried task {} of failed worker", assigned_task.chunk_id);
            task_queue.requery_task(assigned_task);
        }
        self.connections
            .remove(&worker_to_be_removed)
            .expect("This is not expected to fail");
    }

    fn read_rx_buffer(&mut self, worker: Token) -> TransportState<io::Error> {
        let connection = self
            .connections
            .get_mut(&worker)
            .expect("This is not expected to fail");
        let stream = &mut connection.stream;
        let result_buffer = &mut connection.rx_buffer;

        // After we have finished reading the header, the buffer gets enlarged and remaining
        // bytes from the message body will be read to the buffer in subsequent iterations
        loop {
            let read_results = stream.read(result_buffer.buf_mut_unfilled());
            match read_results {
                Ok(0) => {
                    return TransportState::Error(io::Error::new(
                        io::ErrorKind::ConnectionAborted,
                        "Connection aborted",
                    ));
                }
                Ok(bytes_read) => {
                    result_buffer.update_bytes_read(bytes_read);
                }
                // When errors are returned, it's guaranteed that nothing was read during this iteration,
                // so we don't need to check here if we're perhaps finished
                Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {
                    // Retry...
                }
                Err(ref e) if e.kind() == WouldBlock => {
                    // The underlying OS socket is empty, wait for another event to occur
                    return TransportState::Stalled;
                }
                Err(e) => {
                    return TransportState::Error(e);
                }
            }

            if result_buffer.is_finished_reading_header() {
                if result_buffer.decode_header().is_err() {
                    return TransportState::Error(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Could not decode header",
                    ));
                }
            } else if result_buffer.is_finished() {
                return TransportState::Finished;
            }
        }
    }

    fn write_tx_buffer(
        &mut self,
        worker: Token,
        task_queue: &mut TaskQueue,
    ) -> TransportState<io::Error> {
        let connection = self
            .connections
            .get_mut(&worker)
            .expect("This is not expected to fail");
        let stream = &mut connection.stream;
        let send_buffer = &mut connection.tx_buffer;

        if send_buffer.is_ready() {
            if let Some(mut task) = task_queue.next() {
                task.reference_path = Some(self.reference_path.to_string());
                task.alignment_parameters = Some(self.alignment_parameters.to_owned());

                // The task gets stored in the assignment record and send buffer
                connection.assigned_task = Some(task);
                send_buffer.reload(connection.assigned_task.as_mut().unwrap());
            } else {
                return TransportState::Complete;
            }
        }

        loop {
            let write_results = stream.write(send_buffer.buf_unsent());
            match write_results {
                Ok(0) => {
                    return TransportState::Error(io::Error::new(
                        io::ErrorKind::ConnectionAborted,
                        "Connection aborted",
                    ));
                }
                Ok(bytes_sent) => {
                    send_buffer.update_bytes_sent(bytes_sent);
                }
                // When errors are returned, it's guaranteed that nothing was read during this iteration,
                // so we don't need to check here if we're perhaps finished
                Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {
                    // Retry...
                }
                Err(ref e) if e.kind() == WouldBlock => {
                    // The underlying OS socket is empty, wait for another event to occur
                    return TransportState::Stalled;
                }
                Err(e) => {
                    return TransportState::Error(e);
                }
            }

            if send_buffer.is_ready() {
                // The contents of the buffer have been sent successfully
                return TransportState::Finished;
            }
        }
    }
}
