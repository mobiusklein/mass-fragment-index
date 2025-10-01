#![cfg(feature = "binary_storage")]

mod fragment_parquet;
mod peak_parquet;
mod split;
mod util;

pub use fragment_parquet::{read_fragment_index, write_fragment_index};
pub use peak_parquet::{read_peak_index, write_peak_index};
pub use split::{
    BinStorageStrategy, SearchIndexOnDisk, SplitArrowStorage, SplitBand, SplitIndexBinaryStorage,
    SplitStorageOptions,
};
pub use util::{ArrowStorage, IndexBinaryStorage, IndexMetadata};

#[doc(hidden)]
pub use parquet::basic::{BrotliLevel, Compression, GzipLevel, ZstdLevel};
