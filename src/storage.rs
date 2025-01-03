#![cfg(feature = "binary_storage")]

mod peak_parquet;
mod fragment_parquet;
mod util;
mod split;

pub use peak_parquet::{read_peak_index, write_peak_index};
pub use fragment_parquet::{read_fragment_index, write_fragment_index};
pub use util::{ArrowStorage, IndexMetadata, IndexBinaryStorage};
pub use split::{SplitIndexBinaryStorage, SplitBand};

#[doc(hidden)]
pub use parquet::basic::{Compression, ZstdLevel, GzipLevel, BrotliLevel};