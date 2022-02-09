use std::{
    collections::{BinaryHeap, HashMap},
    io::{self, Read, Write},
    iter::Map,
    net::{IpAddr, Ipv4Addr, SocketAddr},
    path::Path,
};

use bio::{data_structures::suffix_array::SuffixArray, io::fastq};
use log::{debug, info, warn};
use mio::{net::TcpListener, *};
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read as BamRead};

use crate::{
    distributed::*,
    errors::{Error, Result},
    map,
    utils::*,
};

enum TransportState {
    // The message has been successfully transferred
    Finished,
    // This operation would block
    Stalled,
    // An error has occurred
    Error(std::io::Error),
    // There are no tasks left in the queue
    Complete,
}

/// Keeps track of the processing state of chunks of reads
struct TaskQueue<T> {
    chunk_id: usize,
    chunk_size: usize,
    records: T,
    requeried_tasks: Vec<TaskSheet>,
}

impl<T> Iterator for TaskQueue<T>
where
    T: Iterator<Item = Result<Record>>,
{
    type Item = Result<TaskSheet>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(task) = self.requeried_tasks.pop() {
            info!("Retrying previously failed task {}", task.chunk_id);
            return Some(Ok(task));
        }

        let chunk = self
            .records
            .by_ref()
            .take(self.chunk_size)
            .collect::<Result<Vec<_>>>();
        self.chunk_id += 1;

        match chunk {
            Ok(inner) if !inner.is_empty() => {
                Some(Ok(TaskSheet::from_records(self.chunk_id - 1, inner, None)))
            }
            Ok(_) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

impl<T> TaskQueue<T>
where
    T: Iterator<Item = Result<Record>>,
{
    fn requery_task(&mut self, task_sheet: TaskSheet) {
        self.requeried_tasks.push(task_sheet);
    }
}

/// Convertible to TaskQueue
trait IntoTaskQueue<E, I, O, T>
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

// Central place to store data relevant for one dispatcher <-> worker connection
struct Connection {
    pub stream: net::TcpStream,
    pub assigned_task: Option<TaskSheet>,
    pub rx_buffer: ResultRxBuffer,
    pub tx_buffer: TaskTxBuffer,
}

pub struct Dispatcher<'a, 'b> {
    reads_path: &'b Path,
    reference_path: &'b Path,
    out_file_path: &'b Path,
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
    ) -> Result<Self> {
        let reads_path = Path::new(reads_path);
        let out_file_path = Path::new(out_file_path);
        if !reads_path.exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "The given input file could not be found",
            )
            .into());
        }
        if out_file_path.exists() {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                "The given output file already exists",
            )
            .into());
        }
        Ok(Self {
            reads_path,
            reference_path: Path::new(reference_path),
            out_file_path,
            alignment_parameters,
            connections: HashMap::new(),
            accept_connections: true,
        })
    }

    pub fn run(&mut self, port: u16) -> Result<()> {
        // Set up networking
        let addr = SocketAddr::new(IpAddr::V4(Ipv4Addr::new(0, 0, 0, 0)), port);
        let mut listener = net::TcpListener::bind(addr)?;

        // Things that determine the concrete values of the generics used in `run_inner(...)`
        // are set up here to allow static dispatch
        info!("Load position map");
        let identifier_position_map =
            load_id_pos_map_from_path(self.reference_path.to_str().ok_or(Error::InvalidIndex(
                "Cannot access the index (file paths contain invalid UTF-8 unicode)",
            ))?)?;

        info!("Load suffix array");
        let sampled_suffix_array_owned = load_suffix_array_from_path(
            self.reference_path.to_str().ok_or(Error::InvalidIndex(
                "Cannot access the index (file paths contain invalid UTF-8 unicode)",
            ))?,
        )?;
        let fmd_index =
            load_index_from_path(self.reference_path.to_str().ok_or(Error::InvalidIndex(
                "Cannot access the index (file paths contain invalid UTF-8 unicode)",
            ))?)?;
        let suffix_array = sampled_suffix_array_owned.into_sampled_suffix_array(
            &fmd_index.bwt,
            &fmd_index.less,
            &fmd_index.occ,
        );

        // Static dispatch of the Record type based on the filename extension
        match self
            .reads_path
            .extension()
            .ok_or(Error::InvalidInputType)?
            .to_str()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Cannot access the input file (file path contains invalid UTF-8 unicode)",
                )
            })? {
            "bam" => {
                let mut reader = bam::Reader::from_path(self.reads_path)?;
                let _ = reader.set_threads(4);
                let header = map::create_bam_header(Some(&reader), &identifier_position_map);
                let mut out_file =
                    bam::Writer::from_path(self.out_file_path, &header, bam::Format::Bam)?;
                let _ = out_file.set_threads(4);
                let mut task_queue = reader
                    .records()
                    .into_tasks(self.alignment_parameters.chunk_size);
                self.run_inner(
                    &mut task_queue,
                    &suffix_array,
                    &identifier_position_map,
                    &mut listener,
                    &mut out_file,
                )
            }
            "fastq" | "fq" => {
                let reader = fastq::Reader::from_file(self.reads_path)?;
                let header = map::create_bam_header(None, &identifier_position_map);
                let mut out_file =
                    bam::Writer::from_path(self.out_file_path, &header, bam::Format::Bam)?;
                let _ = out_file.set_threads(4);
                let mut task_queue = reader
                    .records()
                    .into_tasks(self.alignment_parameters.chunk_size);
                self.run_inner(
                    &mut task_queue,
                    &suffix_array,
                    &identifier_position_map,
                    &mut listener,
                    &mut out_file,
                )
            }
            _ => Err(Error::InvalidInputType),
        }
    }

    /// This part has been extracted from the main run() function to allow static dispatch based on the
    /// input file type
    fn run_inner<S, T>(
        &mut self,
        task_queue: &mut TaskQueue<T>,
        suffix_array: &S,
        identifier_position_map: &map::FastaIdPositions,
        listener: &mut TcpListener,
        out_file: &mut bam::Writer,
    ) -> Result<()>
    where
        S: SuffixArray + Send + Sync,
        T: Iterator<Item = Result<Record>>,
    {
        let mut max_token = Self::DISPATCHER_TOKEN.0;

        // Create a poll instance
        let mut poll = Poll::new()?;

        // Start listening for incoming connection attempts and store connections
        poll.registry()
            .register(listener, Self::DISPATCHER_TOKEN, Interest::READABLE)?;

        // Create storage for events
        let mut events = Events::with_capacity(1024);
        info!("Ready to distribute work");
        loop {
            poll.poll(&mut events, None)?;
            for event in events.iter() {
                match event.token() {
                    // Workers of the world, register!
                    Self::DISPATCHER_TOKEN => {
                        if self.accept_connections {
                            while let Ok((mut remote_stream, remote_addr)) = listener.accept() {
                                info!("Connection established ({:?})", remote_addr);
                                max_token += 1;
                                let remote_token = Token(max_token);
                                poll.registry().register(
                                    &mut remote_stream,
                                    remote_token,
                                    Interest::WRITABLE,
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
                            debug!("Task queue is empty: decline connection attempt");
                        }
                    }

                    //
                    // Communication with existing workers
                    //
                    _ => {
                        // Receive results from workers
                        if event.is_readable() {
                            match self.read_rx_buffer(event.token()) {
                                TransportState::Finished => {
                                    let results = self
                                        .connections
                                        .get_mut(&event.token())
                                        .expect("This is not expected to fail")
                                        .rx_buffer
                                        .decode_and_reset()?;

                                    debug!(
                                        "Worker {} sent results of task {}",
                                        event.token().0,
                                        results.chunk_id,
                                    );
                                    self.write_results(
                                        results.results,
                                        suffix_array,
                                        identifier_position_map,
                                        self.alignment_parameters,
                                        out_file,
                                    )?;
                                    debug!("Finished task {}", results.chunk_id,);

                                    // Remove completed task from assignments
                                    self.connections
                                        .get_mut(&event.token())
                                        .expect("This is not expected to fail")
                                        .assigned_task = None;

                                    poll.registry().reregister(
                                        &mut self
                                            .connections
                                            .get_mut(&event.token())
                                            .expect("This is not expected to fail")
                                            .stream,
                                        event.token(),
                                        Interest::WRITABLE,
                                    )?;
                                }
                                TransportState::Stalled => {
                                    // We're patient!!1!
                                }
                                TransportState::Error(_) => {
                                    warn!(
                                        "Connection is no longer valid, remove worker {} from pool",
                                        event.token().0
                                    );
                                    self.release_worker(event.token(), task_queue);
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
                            // `TransportState` is wrapped in a `Result` to indicate whether or
                            // not there was an error reading the input file so we bubble it up.
                            // `TransportState::Err(e)`, however, means that there was a communication
                            // error with one of the workers which we can recover from by dropping
                            // the connection and requerying its tasks.
                            match self.write_tx_buffer(event.token(), task_queue)? {
                                TransportState::Finished => {
                                    debug!(
                                        "Assigned and sent task {} to worker {}",
                                        self.connections
                                            .get(&event.token())
                                            .expect("This is not expected to fail")
                                            .assigned_task
                                            .as_ref()
                                            .expect("This is not expected to fail")
                                            .chunk_id,
                                        event.token().0,
                                    );
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
                                    warn!(
                                        "Connection is no longer valid, remove worker {} from pool",
                                        event.token().0
                                    );
                                    self.release_worker(event.token(), task_queue);
                                }
                                TransportState::Complete => {
                                    self.accept_connections = false;
                                    self.release_worker(event.token(), task_queue);
                                    if self.connections.is_empty() {
                                        info!("All tasks have been completed, shutting down gracefully");
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
        alignment_parameters: &AlignmentParameters,
        out_file: &mut bam::Writer,
    ) -> Result<()>
    where
        S: SuffixArray + Sync,
    {
        debug!("Translate suffix array intervals to genomic positions");
        let bam_records = hits
            .par_iter_mut()
            .map_init(rand::thread_rng, |mut rng, (record, hit_interval)| {
                map::intervals_to_bam(
                    record,
                    hit_interval,
                    suffix_array,
                    identifier_position_map,
                    None,
                    alignment_parameters,
                    &mut rng,
                )
            })
            .collect::<Result<Vec<_>>>()?;

        debug!("Write to output file");
        for record in bam_records.iter() {
            out_file.write(record)?;
        }

        Ok(())
    }

    /// Removes worker_to_be_removed from the event queue and terminates the connection.
    /// This causes the worker_to_be_removed process to terminate.
    /// Also, this worker's task will be requeried.
    fn release_worker<T>(&mut self, worker_to_be_removed: Token, task_queue: &mut TaskQueue<T>)
    where
        T: Iterator<Item = Result<Record>>,
    {
        if let Some(assigned_task) = self
            .connections
            .get_mut(&worker_to_be_removed)
            .expect("This is not expected to fail")
            .assigned_task
            .take()
        {
            warn!("Requeried task {} of failed worker", assigned_task.chunk_id);
            task_queue.requery_task(assigned_task);
        }
        self.connections
            .remove(&worker_to_be_removed)
            .expect("This is not expected to fail");
    }

    fn read_rx_buffer(&mut self, worker: Token) -> TransportState {
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
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {
                    // Retry...
                }
                Err(e) if e.kind() == io::ErrorKind::WouldBlock => {
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

    /// The resulting `TransportState` is wrapped in a `Result` to indicate
    /// whether or not there was an error reading the input file so we bubble
    /// it up. `TransportState::Err(e)`, however, means that there was a
    /// communication problem with one of the workers, so we can react and
    /// recover accordingly.
    fn write_tx_buffer<T>(
        &mut self,
        worker: Token,
        task_queue: &mut TaskQueue<T>,
    ) -> Result<TransportState>
    where
        T: Iterator<Item = Result<Record>>,
    {
        let connection = self
            .connections
            .get_mut(&worker)
            .expect("This is not expected to fail");
        let stream = &mut connection.stream;
        let send_buffer = &mut connection.tx_buffer;

        if send_buffer.is_ready() {
            if let Some(task) = task_queue.next() {
                let mut task = task?;
                task.reference_path = Some(self.reference_path.to_string_lossy().into_owned());
                task.alignment_parameters = Some(self.alignment_parameters.to_owned());

                // The task gets stored in the assignment record and send buffer
                connection.assigned_task = Some(task);
                send_buffer.reload(
                    connection
                        .assigned_task
                        .as_mut()
                        .expect("This is not expected to fail"),
                );
            } else {
                return Ok(TransportState::Complete);
            }
        }

        loop {
            match stream.write(send_buffer.buf_unsent()) {
                Ok(0) => {
                    return Ok(TransportState::Error(io::Error::new(
                        io::ErrorKind::ConnectionAborted,
                        "Connection aborted",
                    )));
                }
                Ok(bytes_sent) => {
                    send_buffer.update_bytes_sent(bytes_sent);
                }
                // When errors are returned, it's guaranteed that nothing was read during this iteration,
                // so we don't need to check here if we're perhaps finished
                Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {
                    // Retry...
                }
                Err(ref e) if e.kind() == io::ErrorKind::WouldBlock => {
                    // The underlying OS socket is empty, wait for another event to occur
                    return Ok(TransportState::Stalled);
                }
                Err(e) => {
                    return Ok(TransportState::Error(e));
                }
            }

            if send_buffer.is_ready() {
                // The contents of the buffer have been sent successfully
                return Ok(TransportState::Finished);
            }
        }
    }
}
