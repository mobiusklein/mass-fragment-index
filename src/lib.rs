pub mod fragment;
pub mod index;
pub mod interval;
pub mod r#match;
pub mod parent;
pub mod peak;
pub mod sort;

#[cfg(feature = "binary_storage")]
pub mod storage;

pub use crate::fragment::{Fragment, FragmentSeriesParsingError};
pub use crate::index::SearchIndex;
pub use crate::interval::Interval;
pub use crate::parent::{ParentMolecule, Peptide, Spectrum};
pub use crate::peak::{DeconvolutedPeak, MZPeak};
pub use crate::sort::{IndexSortable, MassType, Tolerance, ToleranceParsingError};

pub type PeptideFragmentIndex = SearchIndex<Fragment, Peptide>;
pub type SpectrumIndex = SearchIndex<MZPeak, Spectrum>;
pub type DeconvolutedSpectrumIndex = SearchIndex<DeconvolutedPeak, Spectrum>;
