use std::collections::HashMap;
use std::io;
use std::iter::FusedIterator;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[cfg(feature = "parallelism")]
use rayon::prelude::*;

#[cfg(feature = "binary_storage")]
use crate::storage::{ArrowStorage, IndexBinaryStorage, IndexMetadata, SplitIndexBinaryStorage};

use crate::interval::Interval;
use crate::sort::{
    IndexBin, IndexSortable, MassType, ParentSortedIndexBinSearchIter, SortType, Tolerance,
};

#[derive(Debug, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SearchIndex<T: IndexSortable + Default, P: IndexSortable + Default> {
    pub bins: Vec<IndexBin<T>>,
    pub parents: IndexBin<P>,
    pub(crate) bins_per_dalton: u32,
    pub(crate) max_item_mass: MassType,
    pub(crate) sort_type: SortType,
}

impl<T: IndexSortable + Default, P: IndexSortable + Default> SearchIndex<T, P> {
    pub fn empty(bins_per_dalton: u32, max_fragment_size: MassType) -> Self {
        let mut inst = Self {
            bins_per_dalton,
            max_item_mass: max_fragment_size,
            ..Self::default()
        };
        inst.initialize_bins();
        inst
    }

    fn initialize_bins(&mut self) {
        let mut mass_step: MassType = 0.0;
        self.bins = Vec::new();
        while mass_step < self.max_item_mass {
            let bin: IndexBin<T> = IndexBin::default();
            self.bins.push(bin);
            mass_step += 1 as MassType / self.bins_per_dalton as MassType;
        }
        let bin: IndexBin<T> = IndexBin::default();
        self.bins.push(bin);
    }

    pub fn num_bins(&self) -> usize {
        self.bins.len()
    }

    pub fn num_entries(&self) -> usize {
        self.bins.iter().map(|b| b.len()).sum()
    }

    pub fn iter_bins(&self) -> std::slice::Iter<'_, IndexBin<T>> {
        self.bins.iter()
    }

    pub fn new(
        bins: Vec<IndexBin<T>>,
        parents: IndexBin<P>,
        bins_per_dalton: u32,
        max_item_mass: MassType,
        sort_type: SortType,
    ) -> Self {
        Self {
            bins,
            parents: parents,
            bins_per_dalton,
            max_item_mass,
            sort_type,
        }
    }

    pub fn total_bins_for_mass(&self) -> u32 {
        self.bins_per_dalton * (self.max_item_mass.round() as u32)
    }

    pub fn bin_for_mass(&self, mass: MassType) -> usize {
        let i = (mass * self.bins_per_dalton as MassType).round() as usize;
        let bin_index = if i >= self.bins.len() {
            self.bins.len() - 1
        } else {
            i
        };
        return bin_index;
    }

    pub fn sort(&mut self, ordering: SortType) {
        for bin in self.bins.iter_mut() {
            bin.sort(ordering)
        }
        self.sort_type = ordering
    }

    #[cfg(feature = "parallelism")]
    pub fn par_sort(&mut self, ordering: SortType)
    where
        T: Send,
    {
        self.bins.par_iter_mut().for_each(|bin| bin.sort(ordering));
        self.sort_type = ordering;
    }

    pub fn add_parent(&mut self, parent_molecule: P) {
        self.parents.push(parent_molecule)
    }

    pub fn add(&mut self, entry: T) -> usize {
        let mass = entry.mass();
        let bin_index = self.bin_for_mass(mass);
        self.bins[bin_index].push(entry);
        bin_index
    }

    pub fn parents_for(&self, mass: MassType, error_tolerance: Tolerance) -> Interval {
        let iv = self.parents.search_mass(mass, error_tolerance);
        iv
    }

    pub fn parents_for_range(
        &self,
        low: MassType,
        high: MassType,
        error_tolerance: Tolerance,
    ) -> Interval {
        let mut out = Interval::default();
        out.start = self.parents_for(low, error_tolerance).start;
        out.end = self.parents_for(high, error_tolerance).end;
        out
    }

    pub fn search(
        &self,
        query: MassType,
        error_tolerance: Tolerance,
        parent_interval: Option<Interval>,
    ) -> SearchIndexSearchIter<'_, T, P> {
        SearchIndexSearchIter::new(
            self,
            query,
            error_tolerance,
            parent_interval.unwrap_or_else(|| Interval::new(0, self.parents.len())),
        )
    }

    pub fn bins_per_dalton(&self) -> u32 {
        self.bins_per_dalton
    }

    pub fn max_item_mass(&self) -> f32 {
        self.max_item_mass
    }

    pub fn sort_type(&self) -> SortType {
        self.sort_type
    }
}

#[cfg(feature = "binary_storage")]
mod storage {
    use crate::storage::{SplitArrowStorage, SplitStorageOptions};

    use super::*;

    impl<
            'a,
            T: IndexSortable + Default + ArrowStorage + 'a,
            P: IndexSortable + Default + ArrowStorage + 'a,
        > IndexBinaryStorage<'a, T, P, IndexMetadata> for SearchIndex<T, P>
    {
        fn parents(&self) -> &[P] {
            self.parents.as_slice()
        }

        fn iter_entries(&'a self) -> impl Iterator<Item = &'a [T]> + 'a {
            self.bins.iter().map(|b| b.as_slice())
        }

        fn to_metadata(&self) -> IndexMetadata {
            IndexMetadata {
                bins_per_dalton: self.bins_per_dalton,
                max_item_mass: self.max_item_mass,
            }
        }

        fn from_components(
            metadata: IndexMetadata,
            parents: Vec<P>,
            entries: HashMap<u64, Vec<T>>,
        ) -> Self {
            let mut parents = IndexBin::from(parents);
            parents.assume_sorted(SortType::ByMass);
            let mut this = Self::empty(metadata.bins_per_dalton, metadata.max_item_mass);
            this.parents = parents;
            entries.into_iter().for_each(|(k, b)| {
                let mut bin = IndexBin::from(b);
                bin.assume_sorted(SortType::ByParentId);
                this.bins[k as usize] = bin;
            });
            this
        }
    }

    impl<
            'a,
            T: IndexSortable + Default + ArrowStorage + 'a + Clone + SplitArrowStorage,
            P: IndexSortable + Default + ArrowStorage + 'a,
        > SplitIndexBinaryStorage<'a, T, P, IndexMetadata> for SearchIndex<T, P>
    {
    }

    impl<
            'a,
            T: IndexSortable + Default + ArrowStorage + 'a,
            P: IndexSortable + Default + ArrowStorage + 'a,
        > SearchIndex<T, P>
    {
        pub fn write_parquet<D: AsRef<std::path::Path>>(
            &'a self,
            directory: &D,
            compression_level: Option<parquet::basic::Compression>,
        ) -> io::Result<()> {
            self.write(directory, compression_level)
        }

        pub fn read_parquet<D: AsRef<std::path::Path>>(directory: &D) -> io::Result<Self> {
            Self::read(directory)
        }
    }

    impl<
            'a,
            T: IndexSortable + Default + ArrowStorage + 'a + Clone + SplitArrowStorage,
            P: IndexSortable + Default + ArrowStorage + 'a,
        > SearchIndex<T, P>
    {
        pub fn write_banded_parquet<D: AsRef<std::path::Path>>(
            &'a self,
            directory: &D,
            bin_options: SplitStorageOptions,
            compression_level: Option<parquet::basic::Compression>,
        ) -> io::Result<()> {
            self.write_split(directory, bin_options, compression_level)
        }

        pub fn read_banded_parquet<D: AsRef<std::path::Path>>(directory: &D) -> io::Result<Self> {
            Self::read_split(directory)
        }
    }
}

pub struct SearchIndexBinIter<'a, T: IndexSortable + Default, P: IndexSortable + Default> {
    index: &'a SearchIndex<T, P>,
    pub query: MassType,
    pub error_tolerance: Tolerance,
    pub low_bin: usize,
    pub high_bin: usize,
    pub current_bin: usize,
}

impl<'a, T: IndexSortable + Default, P: IndexSortable + Default> ExactSizeIterator
    for SearchIndexBinIter<'a, T, P>
{
    fn len(&self) -> usize {
        self.high_bin.saturating_sub(self.low_bin)
    }
}

impl<'a, T: IndexSortable + Default, P: IndexSortable + Default> FusedIterator
    for SearchIndexBinIter<'a, T, P>
{
}

impl<'a, T: IndexSortable + Default, P: IndexSortable + Default> Iterator
    for SearchIndexBinIter<'a, T, P>
{
    type Item = &'a IndexBin<T>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_bin()
    }
}

impl<'a, T: IndexSortable + Default, P: IndexSortable + Default> SearchIndexBinIter<'a, T, P> {
    pub fn new(index: &'a SearchIndex<T, P>, query: MassType, error_tolerance: Tolerance) -> Self {
        let (low_mass, high_mass) = error_tolerance.bounds(query);
        let low_bin = index.bin_for_mass(low_mass);
        let high_bin = (index.bin_for_mass(high_mass) + 1).min(index.bins.len());
        Self {
            index,
            query,
            error_tolerance,
            low_bin,
            high_bin,
            current_bin: low_bin,
        }
    }

    #[inline(always)]
    fn next_bin(&mut self) -> Option<&'a IndexBin<T>> {
        if self.current_bin < self.high_bin {
            let bin = unsafe { self.index.bins.get_unchecked(self.current_bin) };
            self.current_bin += 1;
            Some(bin)
        } else {
            None
        }
    }
}

pub struct SearchIndexSearchIter<'a, T: IndexSortable + Default, P: IndexSortable + Default> {
    pub query: MassType,
    pub error_tolerance: Tolerance,
    pub parent_range: Interval,
    bin_iter: SearchIndexBinIter<'a, T, P>,
    item_iter: Option<ParentSortedIndexBinSearchIter<'a, T>>,
}

impl<'a, T: IndexSortable + Default, P: IndexSortable + Default> Iterator
    for SearchIndexSearchIter<'a, T, P>
{
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_entry()
    }
}

impl<'a, T: IndexSortable + Default, P: IndexSortable + Default> SearchIndexSearchIter<'a, T, P> {
    pub fn new(
        index: &'a SearchIndex<T, P>,
        query: MassType,
        error_tolerance: Tolerance,
        parent_range: Interval,
    ) -> Self {
        let bin_iter = SearchIndexBinIter::new(index, query, error_tolerance);
        let mut this = Self {
            query,
            error_tolerance,
            parent_range,
            bin_iter,
            item_iter: None,
        };
        this.next_bin_iterator();
        this
    }

    #[inline(always)]
    fn next_bin_iterator(&mut self) -> bool {
        self.item_iter = self.bin_iter.next().map(|bin| {
            ParentSortedIndexBinSearchIter::new(
                bin,
                self.parent_range,
                self.query,
                self.error_tolerance,
            )
        });
        self.item_iter.is_some()
    }

    #[inline(always)]
    fn next_entry(&mut self) -> Option<&'a T> {
        loop {
            if let Some(t) = self.item_iter.as_mut().and_then(|it| it.next()) {
                return Some(t);
            } else {
                if !self.next_bin_iterator() {
                    return None;
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::parent::Spectrum;
    use crate::peak::DeconvolutedPeak;

    #[test]
    fn test_build() {
        let spectra = vec![
            Spectrum::new(2300.0, 2, 0, 0, 0),
            Spectrum::new(2301.0, 2, 0, 1, 1),
            Spectrum::new(2401.0, 2, 1, 0, 2),
            Spectrum::new(4100.0, 4, 0, 2, 3),
        ];
        let mut parent_list = IndexBin::new(spectra, SortType::Unsorted, 0.0, 0.0);
        parent_list.sort(SortType::ByMass);

        let peaks = vec![
            DeconvolutedPeak::new(251.5, 1, 0.0, 0),
            DeconvolutedPeak::new(251.5, 1, 0.0, 2),
            DeconvolutedPeak::new(251.6, 1, 0.0, 1),
            DeconvolutedPeak::new(303.7, 1, 0.0, 0),
            DeconvolutedPeak::new(501.2, 1, 0.0, 1),
        ];
        let mut index: SearchIndex<DeconvolutedPeak, Spectrum> = SearchIndex::empty(10, 1000.0);
        for peak in peaks {
            index.add(peak);
        }

        let idx = index.bin_for_mass(251.5);
        assert!(idx == 2515);
    }
}
