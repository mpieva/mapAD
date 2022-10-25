use std::iter;

use crate::{
    distributed::{Message, ResultSheet},
    errors::{Error, Result},
    map::input_chunk_reader::TaskSheet,
};

#[derive(Debug)]
pub struct TaskRxBuffer {
    pub expected_size: usize,
    pub buf: Vec<u8>,
}

impl TaskRxBuffer {
    pub fn new() -> Self {
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

    pub fn buf_mut_unfilled(&mut self) -> &mut [u8] {
        if self.expected_size > TaskSheet::PROTO_LEN {
            &mut self.buf[TaskSheet::PROTO_LEN..]
        } else {
            &mut self.buf[..TaskSheet::PROTO_LEN]
        }
    }

    pub fn decode_and_reset(&mut self) -> Result<TaskSheet> {
        let out = bincode::deserialize(&self.buf)?;
        self.reset();
        Ok(out)
    }

    pub fn decode_header(&mut self) -> Result<()> {
        let message_size: u64 = bincode::deserialize(&self.buf)?;
        self.expected_size = message_size
            .try_into()
            .map_err(|_e| Error::ArchitectureError)?;
        self.buf
            .extend(iter::repeat(0).take(self.expected_size - TaskSheet::PROTO_LEN));
        Ok(())
    }
}

#[derive(Debug)]
pub struct ResultRxBuffer {
    pub expected_size: usize,
    pub buf: Vec<u8>,
    pub already_read: usize,
}

impl ResultRxBuffer {
    pub fn new() -> Self {
        Self {
            buf: vec![0; ResultSheet::PROTO_LEN],
            expected_size: ResultSheet::PROTO_LEN,
            already_read: 0,
        }
    }

    fn reset(&mut self) {
        self.buf.clear();
        self.buf
            .extend(iter::repeat(0).take(ResultSheet::PROTO_LEN));
        self.already_read = 0;
    }

    pub fn buf_mut_unfilled(&mut self) -> &mut [u8] {
        &mut self.buf[self.already_read..]
    }

    pub fn decode_and_reset(&mut self) -> Result<ResultSheet> {
        let out = bincode::deserialize(&self.buf)?;
        self.reset();
        Ok(out)
    }

    pub fn decode_header(&mut self) -> Result<()> {
        let message_size: u64 = bincode::deserialize(&self.buf)?;
        self.expected_size = message_size
            .try_into()
            .map_err(|_e| Error::ArchitectureError)?;
        self.buf
            .extend(iter::repeat(0).take(self.expected_size - ResultSheet::PROTO_LEN));
        Ok(())
    }

    pub fn is_finished_reading_header(&self) -> bool {
        self.already_read == ResultSheet::PROTO_LEN
    }

    pub fn is_finished(&self) -> bool {
        self.expected_size > ResultSheet::PROTO_LEN && self.already_read == self.expected_size
    }

    pub fn update_bytes_read(&mut self, bytes_read: usize) {
        self.already_read += bytes_read;
    }
}

pub struct TaskTxBuffer {
    bytes_sent: usize,
    buf: Vec<u8>,
}

impl TaskTxBuffer {
    pub fn new() -> Self {
        Self {
            bytes_sent: 0,
            buf: Vec::new(),
        }
    }

    pub fn reload(&mut self, task_sheet: &mut TaskSheet) {
        self.bytes_sent = 0;
        self.buf = task_sheet.encode();
    }

    pub fn update_bytes_sent(&mut self, bytes_sent: usize) {
        self.bytes_sent += bytes_sent;
    }

    /// Is the buffer ready to be used?
    pub fn is_ready(&self) -> bool {
        self.bytes_sent == self.buf.len()
    }

    pub fn buf_unsent(&self) -> &[u8] {
        &self.buf[self.bytes_sent..]
    }
}
