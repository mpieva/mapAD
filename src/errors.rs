use core::fmt;
use std::{error, io, result, str};

use anyhow;
use clap::crate_version;

use crate::CRATE_NAME;

/// Internally, we only use this Error type and a newtype wrapper around `std::result::Result<T, E>`
/// where `E` is fixed. When, for example, an additional input file type machinery is added,
/// additional `From<E> for Error` impls might be needed to allow to plug in an `Record`-yielding
/// `Iterator` that returns `Result<T, E>` on calls to its `next()` method.
#[derive(Debug)]
pub enum Error {
    FastQ(bio::io::fastq::Error),
    Hts(String),
    Io(io::Error),
    Utf8Error(str::Utf8Error),
    ParseError(String),
    InvalidInputType,
    InvalidIndex(String),
    IndexVersionMismatch,
    AnyhowError(String),
    ContigBoundaryOverlap,
}

impl fmt::Display for Error {
    #[cold]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::FastQ(_err) => write!(f, "Error reading FASTQ file"),
            Error::Hts(err) => write!(f, "Error reading/writing HTS file format: {}", err),
            Error::Io(err) => write!(f, "IO error: {}", err),
            Error::Utf8Error(err) => write!(f, "UTF-8 error: {}", err),
            Error::ParseError(err) => write!(f, "Parse error: {}", err),
            Error::InvalidInputType => write!(f, "Please specify a path to an input file that ends either with \".bam\", \".fq\", or \".fastq\""),
            Error::InvalidIndex(err) => write!(f, "Index is invalid: {}", err),
            Error::IndexVersionMismatch => write!(f, "The provided index is incompatible with version {} of {}. Please re-create the index.", crate_version!(), CRATE_NAME),
            Error::AnyhowError(err) => write!(f, "Internal error: {}", err),
            Error::ContigBoundaryOverlap => write!(f, "Mapped coordinate overlaps contig boundary"),
        }
    }
}

impl From<io::Error> for Error {
    #[cold]
    fn from(e: io::Error) -> Self {
        Error::Io(e)
    }
}

impl From<str::Utf8Error> for Error {
    #[cold]
    fn from(e: str::Utf8Error) -> Self {
        Error::Utf8Error(e)
    }
}

impl From<bincode::Error> for Error {
    #[cold]
    fn from(e: bincode::Error) -> Self {
        match *e {
            bincode::ErrorKind::Io(e) => Error::Io(e),
            _ => Error::InvalidIndex("Invalid encoding".into()),
        }
    }
}

impl From<anyhow::Error> for Error {
    #[cold]
    fn from(e: anyhow::Error) -> Self {
        Error::AnyhowError(e.to_string())
    }
}

impl From<bio::io::fastq::Error> for Error {
    #[cold]
    fn from(e: bio::io::fastq::Error) -> Self {
        Error::FastQ(e)
    }
}

impl error::Error for Error {}

pub type Result<T> = result::Result<T, Error>;
