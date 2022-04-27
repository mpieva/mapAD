pub mod dispatcher;
pub mod worker;

use std::{collections::BinaryHeap, iter};

use serde::{Deserialize, Serialize};

use crate::{errors::Result, io::TaskSheet, map::HitInterval, utils::Record};

pub trait Message {
    const PROTO_LEN: usize = 8;
    fn encode(&mut self) -> Vec<u8>;
}

/// Results wrapped in a Struct for sending back to the dispatcher
#[derive(Serialize, Deserialize, Debug)]
struct ResultSheet {
    encoded_size: u64,
    chunk_id: usize,
    results: Vec<(Record, BinaryHeap<HitInterval>)>,
}

impl Message for ResultSheet {
    fn encode(&mut self) -> Vec<u8> {
        let encoded_size = bincode::serialized_size(self).expect("This is not expected to fail");
        self.encoded_size = encoded_size;
        bincode::serialize(self).expect("This is not expected to fail")
    }
}

impl ResultSheet {
    fn new(read_set: usize, results: Vec<(Record, BinaryHeap<HitInterval>)>) -> Self {
        Self {
            encoded_size: 0,
            chunk_id: read_set,
            results,
        }
    }
}

#[derive(Debug)]
struct TaskRxBuffer {
    pub expected_size: usize,
    pub buf: Vec<u8>,
}

impl TaskRxBuffer {
    fn new() -> Self {
        Self {
            buf: vec![0; TaskSheet::PROTO_LEN],
            expected_size: TaskSheet::PROTO_LEN,
        }
    }

    fn reset(&mut self) {
        self.buf.clear();
        self.buf.extend(iter::repeat(0).take(TaskSheet::PROTO_LEN));
        self.expected_size = TaskSheet::PROTO_LEN;
    }

    fn buf_mut_unfilled(&mut self) -> &mut [u8] {
        if self.expected_size == TaskSheet::PROTO_LEN {
            &mut self.buf
        } else {
            &mut self.buf[TaskSheet::PROTO_LEN..]
        }
    }

    fn decode_and_reset(&mut self) -> Result<TaskSheet> {
        let out = bincode::deserialize(&self.buf)?;
        self.reset();
        Ok(out)
    }

    fn decode_header(&mut self) -> Result<()> {
        let message_size: u64 = bincode::deserialize(&self.buf)?;
        self.expected_size = message_size as usize;
        self.buf
            .extend(iter::repeat(0).take(self.expected_size - TaskSheet::PROTO_LEN));
        Ok(())
    }
}

#[derive(Debug)]
struct ResultRxBuffer {
    pub expected_size: usize,
    pub buf: Vec<u8>,
    pub already_read: usize,
}

impl ResultRxBuffer {
    fn new() -> Self {
        Self {
            buf: vec![0; ResultSheet::PROTO_LEN],
            // Actually, we expect a size of `ResultSheet::PROTO_LEN`,
            // but the value is not read before it gets replaced with
            // the decoded message size
            expected_size: 0,
            already_read: 0,
        }
    }

    fn reset(&mut self) {
        self.buf.clear();
        self.buf
            .extend(iter::repeat(0).take(ResultSheet::PROTO_LEN));
        // Actually, we expect a size of `ResultSheet::PROTO_LEN`,
        // but the value is not read before it gets replaced with
        // the decoded message size
        self.already_read = 0;
    }

    fn buf_mut_unfilled(&mut self) -> &mut [u8] {
        &mut self.buf[self.already_read..]
    }

    fn decode_and_reset(&mut self) -> Result<ResultSheet> {
        let out = bincode::deserialize(&self.buf)?;
        self.reset();
        Ok(out)
    }

    fn decode_header(&mut self) -> Result<()> {
        let message_size: u64 = bincode::deserialize(&self.buf)?;
        self.expected_size = message_size as usize;
        self.buf
            .extend(iter::repeat(0).take(self.expected_size - ResultSheet::PROTO_LEN));
        Ok(())
    }

    fn is_finished_reading_header(&self) -> bool {
        self.already_read == ResultSheet::PROTO_LEN
    }
    fn is_finished(&self) -> bool {
        self.expected_size > 0 && self.expected_size == self.already_read
    }
    fn update_bytes_read(&mut self, bytes_read: usize) {
        self.already_read += bytes_read;
    }
}

struct TaskTxBuffer {
    bytes_sent: usize,
    buf: Vec<u8>,
}

impl TaskTxBuffer {
    fn new() -> Self {
        Self {
            bytes_sent: 0,
            buf: Vec::new(),
        }
    }

    fn reload(&mut self, task_sheet: &mut TaskSheet) {
        self.bytes_sent = 0;
        self.buf = task_sheet.encode();
    }

    fn update_bytes_sent(&mut self, bytes_sent: usize) {
        self.bytes_sent += bytes_sent;
    }

    fn is_ready(&self) -> bool {
        self.bytes_sent == self.buf.len()
    }

    fn buf_unsent(&self) -> &[u8] {
        &self.buf[self.bytes_sent..]
    }
}
