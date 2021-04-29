use clap::crate_version;

use core::fmt;
use std::{error, io, result};

/// Internally, we only use this Error type and a newtype wrapper around `std::result::Result<T, E>`
/// where `E` is fixed. When, for example, an additional input file type machinery is added,
/// additional `From<E> for Error` impls might be needed to allow to plug in an `Record`-yielding
/// `Iterator` that returns `Result<T, E>` on calls to its `next()` method.
#[derive(Debug)]
pub enum Error {
    Hts(rust_htslib::errors::Error),
    Io(io::Error),
    InvalidInputType,
    InvalidIndex(String),
    IndexVersionMismatch,
}

impl fmt::Display for Error {
    #[cold]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::Hts(_err) => write!(f, "HTS error: Details unknown"),
            Error::Io(err) => write!(f, "IO error: {}", err),
            Error::InvalidInputType => write!(f, "Please specify a path to an input file that ends either with \".bam\", \".fq\", or \".fastq\""),
            Error::InvalidIndex(err) => write!(f, "Index is invalid: {}", err),
            Error::IndexVersionMismatch => write!(f, "The provided index is incompatible with version {} of {}. Please re-create the index.", crate_version!(), crate::map::CRATE_NAME),
        }
    }
}

impl From<io::Error> for Error {
    #[cold]
    fn from(e: io::Error) -> Self {
        Error::Io(e)
    }
}

impl From<bincode::Error> for Error {
    #[cold]
    fn from(e: bincode::Error) -> Self {
        match *e {
            bincode::ErrorKind::Io(e) => Error::Io(e),
            _ => Error::InvalidIndex("Invalid encoding".to_string()),
        }
    }
}

impl From<rust_htslib::errors::Error> for Error {
    #[cold]
    fn from(e: rust_htslib::errors::Error) -> Self {
        Self::Hts(e)
    }
}

impl error::Error for Error {}

pub type Result<T> = result::Result<T, Error>;
