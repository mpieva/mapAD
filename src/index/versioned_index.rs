use std::{fs::File, path::Path};

use serde::{de::DeserializeOwned, Deserialize, Serialize};
use snap::read::FrameDecoder;

use crate::errors::{Error, Result};

/// Versioned index data. A version number is attached to all on-disk index files to ensure we can
/// detect an incompatible on-disk index. The index data can only be extracted if the version
/// matches exactly the currently defined `INDEX_VERSION`.
#[derive(Serialize, Deserialize)]
pub struct VersionedIndexItem<T> {
    version: u8,
    data: T,
}

impl<T> VersionedIndexItem<T> {
    /// Increase this number once the on-disk index changes
    const INDEX_VERSION: u8 = 5;

    /// Creates a new versioned index item with the version set to the current value of
    /// `Self::INDEX_VERSION`
    pub fn new(data: T) -> Self {
        Self {
            version: Self::INDEX_VERSION,
            data,
        }
    }

    /// Returns inner data only if the version of the deserialized Struct is compatible
    pub fn try_take(self) -> Result<T> {
        if self.version == Self::INDEX_VERSION {
            Ok(self.data)
        } else {
            Err(Error::IndexVersionMismatch {
                running: Self::INDEX_VERSION,
                on_disk: self.version,
            })
        }
    }
}

impl<T> VersionedIndexItem<T>
where
    T: DeserializeOwned,
{
    /// Reads a versioned index item from a given `Path`
    pub fn read_from_path<P>(path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        Ok(File::open(path)
            .map(FrameDecoder::new)
            .map(bincode::deserialize_from)??)
    }
}
