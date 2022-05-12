pub mod dispatcher;
pub mod worker;

mod comm_buffers;

use std::{collections::BinaryHeap, mem};

use serde::{Deserialize, Serialize};

use crate::map::{input_chunk_reader::TaskSheet, record::Record, HitInterval};

/// Implementors must ensure that the first element in the struct is the
/// size of the serialized struct, stored in an `u64`
pub trait Message {
    // Size of the first field (`encoded_size: u64`) in bytes
    const PROTO_LEN: usize = mem::size_of::<u64>();
    fn encode(&mut self) -> Vec<u8>;
}

/// Results wrapped in a Struct for sending back to the dispatcher
#[derive(Serialize, Deserialize, Debug)]
pub struct ResultSheet {
    encoded_size: u64, // must be of type `u64` and the first element
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
