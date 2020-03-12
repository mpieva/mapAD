use crate::{distributed::*, map, utils::*};
use bio::data_structures::suffix_array::SuffixArray;
use log::debug;
use mio::{
    net::{TcpListener, TcpStream},
    *,
};
use rand;
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read as BamRead};
use std::{
    collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap},
    error::Error,
    fs::File,
    io::{ErrorKind::WouldBlock, Read, Write},
    iter::Peekable,
    net::{IpAddr, Ipv4Addr, SocketAddr},
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
}

impl<'a> Iterator for TaskQueue<'a> {
    type Item = TaskSheet;

    fn next(&mut self) -> Option<TaskSheet> {
        let source_iterator_loan = &mut self.records;

        // If the underlying iterator is exhausted return None, too
        source_iterator_loan.peek()?;

        self.chunk_id += 1;
        let chunk = source_iterator_loan
            .take(self.chunk_size)
            .map(|record| match record {
                Ok(record) => Ok(record.into()),
                Err(e) => Err(e),
            })
            .collect::<Result<Vec<_>, _>>()
            .expect("Input file is corrupt. Cancelling process.");

        Some(TaskSheet::from_records(self.chunk_id - 1, chunk, None))
    }
}

impl<'a> TaskQueue<'a> {
    fn from_reader(reader: &'a mut bam::Reader, chunk_size: usize) -> Self {
        Self {
            chunk_id: 0,
            chunk_size,
            records: reader.records().peekable(),
        }
    }
}

pub struct Dispatcher<'a, 'b> {
    reads_path: &'b str,
    reference_path: &'b str,
    out_file_path: &'b str,
    alignment_parameters: &'a AlignmentParameters,
    connections: HashMap<Token, TcpStream>,
    result_buffers: HashMap<Token, ResultRxBuffer>,
    send_buffers: HashMap<Token, TaskTxBuffer>,
    checklist: BTreeMap<usize, Token>,
    failed_tasks: BTreeSet<usize>,
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
            result_buffers: HashMap::new(),
            send_buffers: HashMap::new(),
            checklist: BTreeMap::new(),
            failed_tasks: BTreeSet::new(),
        })
    }

    pub fn run(&mut self, port: u16) -> Result<(), Box<dyn Error>> {
        // Set up networking
        let mut max_token = Self::DISPATCHER_TOKEN.0;
        let addr = SocketAddr::new(IpAddr::V4(Ipv4Addr::new(0, 0, 0, 0)), port);
        let listener = TcpListener::bind(&addr)?;

        // Create a poll instance
        let poll = Poll::new()?;

        // Start listening for incoming connection attempts and store connections
        poll.register(
            &listener,
            Self::DISPATCHER_TOKEN,
            Ready::readable(),
            PollOpt::edge(),
        )?;

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
            TaskQueue::from_reader(&mut bam_reader, self.alignment_parameters.chunk_size)
                .peekable();
        let mut out_file = bam::Writer::from_path(&self.out_file_path, &header, bam::Format::BAM)?;
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
                        while let Ok((remote_stream, remote_addr)) = listener.accept() {
                            debug!("Connection established ({:?})", remote_addr);
                            max_token += 1;
                            let remote_token = Token(max_token);
                            poll.register(
                                &remote_stream,
                                remote_token,
                                Ready::all(),
                                PollOpt::edge(),
                            )?;
                            self.connections.insert(remote_token, remote_stream);
                        }
                    }

                    //
                    // Communication with existing workers
                    //
                    _ => {
                        // Receive results from workers
                        if event.readiness().is_readable() {
                            // After finishing previous tasks, the worker is ready to receive fresh jobs
                            match self.read_rx_buffer(event.token()) {
                                TransportState::Finished => {
                                    if let Ok(results) = self
                                        .result_buffers
                                        .get_mut(&event.token())
                                        .expect("This is not expected to fail")
                                        .decode_and_reset()
                                    {
                                        debug!(
                                            "Worker has sent data: {:?}, writing it down.",
                                            event
                                        );
                                        self.write_results(
                                            results.results,
                                            &suffix_array,
                                            &identifier_position_map,
                                            &mut out_file,
                                        )?;

                                        self.mark_task_completed(results.chunk_id);
                                    }

                                    poll.reregister(
                                        self.connections
                                            .get(&event.token())
                                            .expect("This is not expected to fail"),
                                        event.token(),
                                        Ready::all(),
                                        PollOpt::edge(),
                                    )?;
                                }
                                TransportState::Stalled => {
                                    poll.reregister(
                                        self.connections
                                            .get(&event.token())
                                            .expect("This is not expected to fail"),
                                        event.token(),
                                        Ready::readable(),
                                        PollOpt::edge(),
                                    )?;
                                }
                                TransportState::Error(_) => {
                                    debug!("Connection is no longer valid, removing worker from the pool.");
                                    self.release_worker(event.token());
                                }
                                TransportState::Complete => {
                                    // TransportState::Complete means the input file has been processed completely
                                    // which workers don't know. That's why this match arm here is
                                    unreachable!();
                                }
                            }
                        }

                        // Distribution of work
                        if event.readiness().is_writable() {
                            match self.write_tx_buffer(event.token(), &mut task_queue) {
                                TransportState::Finished => {
                                    poll.reregister(
                                        self.connections
                                            .get(&event.token())
                                            .expect("This is not expected to fail"),
                                        event.token(),
                                        Ready::readable(),
                                        PollOpt::edge(),
                                    )?;
                                }
                                TransportState::Stalled => {
                                    poll.reregister(
                                        self.connections
                                            .get(&event.token())
                                            .expect("This is not expected to fail"),
                                        event.token(),
                                        Ready::all(),
                                        PollOpt::edge(),
                                    )?;
                                }
                                TransportState::Error(_) => {
                                    debug!("Connection is no longer valid, removing worker from the pool.");
                                    self.release_worker(event.token());
                                }
                                TransportState::Complete => {
                                    self.release_worker(event.token());
                                }
                            }
                        }
                    }
                }
            }
            if self.all_tasks_finished(&mut task_queue) {
                debug!("All tasks have been completed, shutting down gracefully");
                return Ok(());
            }
        }
    }

    fn all_tasks_finished(&mut self, task_queue: &mut Peekable<TaskQueue>) -> bool {
        task_queue.peek().is_none() && self.checklist.is_empty() && self.failed_tasks.is_empty()
    }

    fn mark_task_completed(&mut self, chunk_id: usize) {
        self.checklist.remove(&chunk_id);
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
    fn release_worker(&mut self, worker_to_be_removed: Token) {
        // `None` variant is non-fatal here
        let _ = self.connections.remove(&worker_to_be_removed);

        // Try to restart pending tasks that are assigned to this worker
        let mut additional_failed_tasks = self
            .checklist
            .iter()
            .filter(|(_, &worker)| worker == worker_to_be_removed)
            .map(|(&chunk_id, _)| chunk_id)
            .collect::<BTreeSet<_>>();

        self.failed_tasks.append(&mut additional_failed_tasks);

        self.checklist = self
            .checklist
            .iter()
            .filter(|(k, _)| !additional_failed_tasks.contains(k))
            .map(|(&k, &v)| (k, v))
            .collect();
    }

    fn read_rx_buffer(&mut self, worker: Token) -> TransportState<std::io::Error> {
        let connection = self
            .connections
            .get_mut(&worker)
            .expect("This is not expected to fail");
        let result_buffer = self
            .result_buffers
            .entry(worker)
            .or_insert_with(ResultRxBuffer::new);

        // After we have finished reading the header, the buffer gets enlarged and remaining
        // bytes from the message body will be read to the buffer in subsequent iterations
        loop {
            let read_results = connection.read(result_buffer.buf_mut_unfilled());
            match read_results {
                Ok(0) => {
                    return TransportState::Error(std::io::Error::new(
                        std::io::ErrorKind::ConnectionAborted,
                        "Connection aborted",
                    ));
                }
                Ok(bytes_read) => {
                    result_buffer.update_bytes_read(bytes_read);
                }
                // When errors are returned, it's guaranteed that nothing was read during this iteration,
                // so we don't need to check here if we're perhaps finished
                Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => {
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
                    return TransportState::Error(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
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
        task_queue: &mut Peekable<TaskQueue>,
    ) -> TransportState<std::io::Error> {
        let connection = self
            .connections
            .get_mut(&worker)
            .expect("This is not expected to fail");
        let send_buffer = self
            .send_buffers
            .entry(worker)
            .or_insert_with(TaskTxBuffer::new);

        if send_buffer.is_ready() {
            if let Some(mut task) = task_queue.next() {
                // TODO (optimization): Don't send parameters every time
                task.reference_path = Some(self.reference_path.to_string());
                task.alignment_parameters = Some(self.alignment_parameters.to_owned());

                // Mark task pending
                self.checklist.insert(task.chunk_id, worker);

                send_buffer.reload(task);
            } else {
                return TransportState::Complete;
            }
        }

        loop {
            let write_results = connection.write(send_buffer.buf_unsent());
            match write_results {
                Ok(0) => {
                    return TransportState::Error(std::io::Error::new(
                        std::io::ErrorKind::ConnectionAborted,
                        "Connection aborted",
                    ));
                }
                Ok(bytes_sent) => {
                    send_buffer.update_bytes_sent(bytes_sent);
                }
                // When errors are returned, it's guaranteed that nothing was read during this iteration,
                // so we don't need to check here if we're perhaps finished
                Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => {
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
