use std::{
    collections::{BinaryHeap, HashMap},
    fs::OpenOptions,
    io::{self, Read, Write},
    net::{IpAddr, Ipv4Addr, SocketAddr},
    path::Path,
    time::Duration,
};

use bio::data_structures::suffix_array::SuffixArray;
use log::{debug, info, warn};
use mio::{net, Events, Interest, Poll, Token};
use noodles::{
    bam,
    sam::{self, alignment::io::Write as AlignmentWrite},
};
use rayon::prelude::*;

use crate::{
    distributed::comm_buffers::{ResultRxBuffer, TaskTxBuffer},
    errors::{Error, Result},
    index::{
        load_id_pos_map_from_path, load_index_from_path, load_original_symbols_from_path,
        load_suffix_array_from_path, FastaIdPositions, OriginalSymbols,
    },
    map::{
        input_chunk_reader::{InputSource, TaskQueue, TaskSheet},
        mapping::{create_bam_header, intervals_to_bam},
        record::Record,
        AlignmentParameters, HitInterval,
    },
};

enum TransportState {
    // The message has been successfully transferred
    Finished,
    // This operation would block
    Stalled,
    // An error has occurred
    #[allow(dead_code)]
    Error(io::Error),
    // There are no tasks left in the queue
    Complete,
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
        if !reads_path.exists() && reads_path.to_str() != Some("-") {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "The given input file could not be found",
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

        let reference_path = self.reference_path.to_str().ok_or_else(|| {
            Error::InvalidIndex(
                "Cannot access the index (file paths contain invalid UTF-8 unicode)".into(),
            )
        })?;

        info!("Load position map");
        let identifier_position_map = load_id_pos_map_from_path(reference_path)?;

        info!("Load original symbols");
        let original_symbols = load_original_symbols_from_path(reference_path)?;

        info!("Load suffix array");
        let sampled_suffix_array_owned = load_suffix_array_from_path(reference_path)?;
        let fmd_index = load_index_from_path(reference_path)?;
        let suffix_array = sampled_suffix_array_owned.into_sampled_suffix_array(
            &fmd_index.bwt,
            &fmd_index.less,
            &fmd_index.occ,
        );

        let mut out_file = bam::io::Writer::new(
            OpenOptions::new()
                .read(false)
                .write(true)
                .create_new(true)
                .open(self.out_file_path)?,
        );

        let mut input_source = InputSource::from_path(self.reads_path)?;
        let out_header = create_bam_header(input_source.header(), &identifier_position_map)?;
        out_file.write_header(&out_header)?;
        let mut task_queue = input_source.task_queue(self.alignment_parameters.chunk_size);
        self.run_inner(
            &mut task_queue,
            &suffix_array,
            &identifier_position_map,
            &original_symbols,
            &mut listener,
            &out_header,
            &mut out_file,
        )
    }

    /// This part has been extracted from the main run() function to allow static dispatch based on the
    /// input file type
    #[allow(clippy::too_many_arguments)]
    fn run_inner<S, T, W>(
        &mut self,
        task_queue: &mut TaskQueue<T>,
        suffix_array: &S,
        identifier_position_map: &FastaIdPositions,
        original_symbols: &OriginalSymbols,
        listener: &mut net::TcpListener,
        out_header: &sam::Header,
        out_file: &mut bam::io::Writer<W>,
    ) -> Result<()>
    where
        S: SuffixArray + Send + Sync,
        T: Iterator<Item = Result<Record>>,
        W: Write,
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
                if event.token() == Self::DISPATCHER_TOKEN {
                    // Workers of the world, register!
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
                } else {
                    //
                    // Communication with existing workers
                    //
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
                                Self::write_results(
                                    results.results,
                                    suffix_array,
                                    identifier_position_map,
                                    original_symbols,
                                    self.alignment_parameters,
                                    out_header,
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
                                        .expect("This is not expected to fail"),
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
                                    info!(
                                        "All tasks have been completed, shutting down gracefully"
                                    );
                                    return Ok(());
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn write_results<S, W>(
        hits: Vec<(Record, BinaryHeap<HitInterval>, Duration)>,
        suffix_array: &S,
        identifier_position_map: &FastaIdPositions,
        original_symbols: &OriginalSymbols,
        alignment_parameters: &AlignmentParameters,
        out_header: &sam::Header,
        out_file: &mut bam::io::Writer<W>,
    ) -> Result<()>
    where
        S: SuffixArray + Sync,
        W: Write,
    {
        debug!("Translate suffix array intervals to genomic positions");
        let bam_records = hits
            .into_par_iter()
            .map_init(
                rand::thread_rng,
                |mut rng, (record, hit_interval, duration)| {
                    intervals_to_bam(
                        record,
                        hit_interval,
                        suffix_array,
                        identifier_position_map,
                        original_symbols,
                        Some(&duration),
                        alignment_parameters,
                        &mut rng,
                    )
                },
            )
            .collect::<Result<Vec<_>>>()?;

        debug!("Write chunk of BAM records to output file");
        for record in &bam_records {
            bam::io::Writer::write_alignment_record(out_file, out_header, record)?;
        }

        Ok(())
    }

    /// Removes `worker_to_be_removed` from the event queue and terminates the connection.
    /// This causes the `worker_to_be_removed` process to terminate.
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
            warn!("Requery task {} of failed worker", assigned_task);
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
            if let Some(mut task) = task_queue.next() {
                task.set_reference_path(
                    self.reference_path
                        .to_str()
                        .expect("This is not expected to fail"),
                );
                task.set_alignment_parameters(self.alignment_parameters);
                // The task gets stored in the assignment record and send buffer
                send_buffer.reload(&mut task);
                connection.assigned_task = Some(task);
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
