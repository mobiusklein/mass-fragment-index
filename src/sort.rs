use std::{
    error::Error,
    fmt::Display,
    iter::FusedIterator,
    ops::{Index, Mul},
    str::FromStr,
};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::interval::Interval;

pub type ParentID = u32;
pub type MassType = f32;

pub fn _isclose(x: MassType, y: MassType, rtol: MassType, atol: MassType) -> bool {
    (x - y).abs() <= (atol + rtol * y.abs())
}

pub fn isclose(x: MassType, y: MassType) -> bool {
    _isclose(x, y, 1e-5, 1e-8)
}

pub fn aboutzero(x: MassType) -> bool {
    isclose(x, 0.0)
}

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum SortType {
    ByMass,
    ByParentId,
    Unsorted,
}

impl Default for SortType {
    fn default() -> Self {
        Self::Unsorted
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Tolerance {
    PPM(MassType),
    Da(MassType),
}

#[derive(Debug, PartialEq, Eq)]
pub enum ToleranceParsingError {
    UnknownUnit,
    InvalidMagnitude,
}

impl Display for ToleranceParsingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{:?}", self))
    }
}

impl Error for ToleranceParsingError {}

impl FromStr for Tolerance {
    type Err = ToleranceParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let n = s.len();
        if n <= 2 {
            return Err(ToleranceParsingError::InvalidMagnitude);
        }
        let s = s.to_lowercase();
        if s.ends_with("da") {
            if let Ok(magnitude) = s[0..n - 2].parse::<MassType>() {
                Ok(Self::Da(magnitude))
            } else {
                Err(ToleranceParsingError::InvalidMagnitude)
            }
        } else if s.ends_with("ppm") {
            if let Ok(magnitude) = s[0..n - 3].parse::<MassType>() {
                Ok(Self::PPM(magnitude))
            } else {
                Err(ToleranceParsingError::InvalidMagnitude)
            }
        } else {
            Err(ToleranceParsingError::UnknownUnit)
        }
    }
}

impl Tolerance {
    pub fn bounds(&self, query: MassType) -> (MassType, MassType) {
        match self {
            Tolerance::PPM(tol) => {
                let width = query * *tol / 1e6;
                (query - width, query + width)
            }
            Tolerance::Da(tol) => (query - *tol, query + *tol),
        }
    }

    pub fn test(&self, query: MassType, reference: MassType) -> bool {
        let (lower_bound, upper_bound) = self.bounds(reference);
        query >= lower_bound && query <= upper_bound
    }

    pub fn format_error(&self, query: MassType, reference: MassType) -> String {
        match self {
            Self::PPM(_tol) => {
                let magnitude = (query - reference) / reference * 1e6;
                format!("{}PPM", magnitude).to_string()
            }
            Self::Da(_tol) => {
                let magnitude = query - reference;
                format!("{}Da", magnitude).to_string()
            }
        }
    }
}

impl Mul<MassType> for Tolerance {
    type Output = Tolerance;

    fn mul(self, rhs: MassType) -> Self::Output {
        match self {
            Self::Da(val) => Self::Da(rhs * val),
            Self::PPM(val) => Self::PPM(rhs * val),
        }
    }
}

pub trait IndexSortable {
    fn mass(&self) -> MassType;
    fn parent_id(&self) -> ParentID;
}

impl<T: IndexSortable> IndexSortable for &T {
    fn mass(&self) -> MassType {
        (*self).mass()
    }

    fn parent_id(&self) -> ParentID {
        (*self).parent_id()
    }
}

#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct IndexBin<T: IndexSortable> {
    pub(crate) entries: Vec<T>,
    pub(crate) sort_type: SortType,
    pub(crate) min_mass: MassType,
    pub(crate) max_mass: MassType,
}

impl<T: IndexSortable> IndexBin<T> {
    pub fn new(
        entries: Vec<T>,
        sort_type: SortType,
        min_mass: MassType,
        max_mass: MassType,
    ) -> Self {
        Self {
            entries,
            sort_type,
            min_mass,
            max_mass,
        }
    }

    pub fn push(&mut self, entry: T) {
        self.entries.push(entry);
        self.sort_type = SortType::Unsorted;
    }

    pub fn find_min_max_masses(&self) -> (MassType, MassType) {
        let mut min_mass = MassType::INFINITY;
        let mut max_mass = 0.0 as MassType;

        for f in self.entries.iter() {
            if f.mass() < min_mass {
                min_mass = f.mass();
            }
            if f.mass() > max_mass {
                max_mass = f.mass();
            }
        }
        return (min_mass, max_mass);
    }

    pub fn sort(&mut self, ordering: SortType) {
        match ordering {
            SortType::ByMass => {
                self.entries
                    .sort_by(|a, b| a.mass().partial_cmp(&b.mass()).unwrap());
                if let Some(f) = self.entries.first() {
                    self.min_mass = f.mass()
                }
                if let Some(f) = self.entries.last() {
                    self.max_mass = f.mass()
                }
            }
            SortType::ByParentId => {
                self.entries
                    .sort_by(|a, b| a.parent_id().cmp(&b.parent_id()));
                (self.min_mass, self.max_mass) = self.find_min_max_masses();
            }
            SortType::Unsorted => {
                (self.min_mass, self.max_mass) = self.find_min_max_masses();
            }
        }
        self.sort_type = ordering;
    }

    pub fn assume_sorted(&mut self, sort_type: SortType) {
        (self.min_mass, self.max_mass) = self.find_min_max_masses();
        self.sort_type = sort_type;
    }

    pub fn len(&self) -> usize {
        return self.entries.len();
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<T> {
        self.entries.iter()
    }

    #[allow(unused)]
    pub(crate) fn iter_mut(&mut self) -> std::slice::IterMut<T> {
        self.entries.iter_mut()
    }

    pub fn as_slice(&self) -> &[T] {
        &self.entries
    }

    pub fn first(&self) -> Option<&T> {
        self.entries.first()
    }

    pub fn get(&self, index: usize) -> Option<&T> {
        self.entries.get(index)
    }

    pub fn last(&self) -> Option<&T> {
        self.entries.last()
    }

    pub fn search_mass(&self, query: MassType, error_tolerance: Tolerance) -> Interval {
        let (lower_bound, upper_bound) = error_tolerance.bounds(query);

        let mut lower_i = self
            .entries
            .partition_point(|entry| entry.mass() <= lower_bound);
        let mut upper_i = self.entries[lower_i..self.len()]
            .partition_point(|entry| entry.mass() <= upper_bound)
            + lower_i;

        while lower_i > 0 {
            if error_tolerance.test(query, self.entries[lower_i - 1].mass()) {
                lower_i -= 1;
            } else {
                break;
            }
        }

        while upper_i + 1 < self.len() {
            if error_tolerance.test(query, self.entries[upper_i + 1].mass()) {
                upper_i += 1;
            } else {
                break;
            }
        }

        return Interval::new(lower_i, upper_i);
    }

    pub fn search_parent_id(&self, parent_id_range: Interval) -> Interval {
        let mut result = Interval::default();
        if self.is_empty() {
            return result
        }
        let start_idx = self.entries.partition_point(|i| {
            parent_id_range.start > i.parent_id() as usize
        });

        let end_idx = self.entries.partition_point(|i| {
            parent_id_range.end > i.parent_id() as usize
        });
        result.start = start_idx;
        result.end = end_idx;
        result
    }

    pub fn min_mass(&self) -> f32 {
        self.min_mass
    }

    pub fn max_mass(&self) -> f32 {
        self.max_mass
    }

    pub fn sort_type(&self) -> SortType {
        self.sort_type
    }
}

impl<T: IndexSortable + Default> Index<usize> for IndexBin<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.entries[index]
    }
}

impl<I: IndexSortable> FromIterator<I> for IndexBin<I> {
    fn from_iter<T: IntoIterator<Item = I>>(iter: T) -> Self {
        let entries = Vec::from_iter(iter);
        let bin = entries.into();
        bin
    }
}

impl<T: IndexSortable + PartialEq> PartialEq for IndexBin<T> {
    fn eq(&self, other: &Self) -> bool {
        self.entries == other.entries
            && self.sort_type == other.sort_type
            && self.min_mass == other.min_mass
            && self.max_mass == other.max_mass
    }
}

impl<T: IndexSortable> From<Vec<T>> for IndexBin<T> {
    fn from(value: Vec<T>) -> Self {
        let mut this = Self::new(value, SortType::Unsorted, 0.0, 0.0);
        let (min, max) = this.find_min_max_masses();
        this.min_mass = min;
        this.max_mass = max;
        this
    }
}

#[derive(Debug)]
pub struct ParentSortedIndexBinSearchIter<'a, T: IndexSortable> {
    bin_iter: std::slice::Iter<'a, T>,
    parent_range: Interval,
    query: f32,
    error_tolerance: Tolerance,
    spanned: bool,
}

impl<'a, T: IndexSortable> FusedIterator for ParentSortedIndexBinSearchIter<'a, T> {}

impl<'a, T: IndexSortable> Iterator for ParentSortedIndexBinSearchIter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_entry()
    }
}

impl<'a, T: IndexSortable> ParentSortedIndexBinSearchIter<'a, T> {
    pub fn new(
        bin: &'a IndexBin<T>,
        parent_range: Interval,
        query: f32,
        error_tolerance: Tolerance,
    ) -> Self {
        let bin_iter = bin.iter();
        let (lo, hi) = error_tolerance.bounds(query);
        let spanned = lo <= bin.min_mass && hi >= bin.max_mass;
        Self {
            bin_iter,
            parent_range,
            query,
            error_tolerance,
            spanned,
        }
    }

    fn next_entry(&mut self) -> Option<&'a T> {
        while let Some(t) = self.bin_iter.next() {
            if self.spanned && self.parent_range.contains(t.parent_id() as usize) {
                return Some(t);
            }
            if self.error_tolerance.test(self.query, t.mass())
                && self.parent_range.contains(t.parent_id() as usize)
            {
                return Some(t);
            }
        }
        None
    }
}

mod soa_bin {
    use std::{marker::PhantomData, ops::RangeBounds};

    use soa_derive::*;

    use super::*;

    pub trait SoAIndexSortable<T: StructOfArray>: SoAVec<T>
    where
        for<'t> Self::Ref<'t>: IndexSortable,
    {
        fn mass(&self) -> &[MassType];

        fn parent_id(&self) -> &[ParentID];

        fn convert_ref(val_ref: Self::Ref<'_>) -> T;
    }

    pub trait SoAIndexSortableSlice<'a, T: StructOfArray>: SoASlice<T> + 'a
    where
        Self::Ref<'a>: IndexSortable
    {
        fn mass(&self) -> &[MassType];

        fn parent_id(&self) -> &[ParentID];
    }


    #[derive(Debug, Clone)]
    #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
    pub struct SoAIndexBin<T: StructOfArray, V: SoAVec<T> + SoAIndexSortable<T> = <T as StructOfArray>::Type>
    where
        for<'t> V::Ref<'t>: IndexSortable,
    {
        entries: V,
        pub(crate) sort_type: SortType,
        pub(crate) min_mass: MassType,
        pub(crate) max_mass: MassType,
        _t: PhantomData<T>,
    }

    impl<T: StructOfArray, V: SoAVec<T> + SoAIndexSortable<T>> From<V> for SoAIndexBin<T, V>
    where
        for<'t> V::Ref<'t>: IndexSortable,
    {
        fn from(value: V) -> Self {
            let mut this = Self::new(value, SortType::Unsorted, 0.0, 0.0);
            (this.min_mass, this.max_mass) = this.find_min_max_masses();
            this
        }
    }

    impl<T: StructOfArray, V: SoAVec<T> + SoAIndexSortable<T>> Default for SoAIndexBin<T, V>
    where
        for<'t> V::Ref<'t>: IndexSortable,
    {
        fn default() -> Self {
            Self {
                entries: V::new(),
                sort_type: SortType::Unsorted,
                min_mass: Default::default(),
                max_mass: Default::default(),
                _t: Default::default(),
            }
        }
    }

    impl<T: StructOfArray, V: SoAVec<T> + SoAIndexSortable<T>> SoAIndexBin<T, V>
    where
        for<'t> V::Ref<'t>: IndexSortable,
    {
        pub fn new(
            entries: V,
            sort_type: SortType,
            min_mass: MassType,
            max_mass: MassType,
        ) -> Self {
            Self {
                entries,
                sort_type,
                min_mass,
                max_mass,
                _t: PhantomData,
            }
        }

        pub fn push(&mut self, entry: T) {
            self.entries.push(entry);
            self.sort_type = SortType::Unsorted;
        }

        pub fn find_min_max_masses(&self) -> (MassType, MassType) {
            let mut min_mass = MassType::INFINITY;
            let mut max_mass = 0.0 as MassType;

            for f in self.entries.iter() {
                if f.mass() < min_mass {
                    min_mass = f.mass();
                }
                if f.mass() > max_mass {
                    max_mass = f.mass();
                }
            }
            return (min_mass, max_mass);
        }

        pub fn sort(&mut self, ordering: SortType) {
            match ordering {
                SortType::ByMass => {
                    self.entries
                        .sort_by(|a, b| a.mass().partial_cmp(&b.mass()).unwrap());
                    if let Some(f) = self.entries.first() {
                        self.min_mass = f.mass()
                    }
                    if let Some(f) = self.entries.last() {
                        self.max_mass = f.mass()
                    }
                }
                SortType::ByParentId => {
                    self.entries
                        .sort_by(|a, b| a.parent_id().cmp(&b.parent_id()));
                    (self.min_mass, self.max_mass) = self.find_min_max_masses();
                }
                SortType::Unsorted => {
                    (self.min_mass, self.max_mass) = self.find_min_max_masses();
                }
            }
            self.sort_type = ordering;
        }

        pub fn assume_sorted(&mut self, sort_type: SortType) {
            (self.min_mass, self.max_mass) = self.find_min_max_masses();
            self.sort_type = sort_type;
        }

        pub fn len(&self) -> usize {
            return self.entries.len();
        }

        pub fn is_empty(&self) -> bool {
            self.entries.is_empty()
        }

        pub fn iter(&self) -> <V as SoAVec<T>>::Iter<'_> {
            self.entries.iter()
        }

        #[allow(unused)]
        pub(crate) fn iter_mut(&mut self) -> <V as SoAVec<T>>::IterMut<'_> {
            self.entries.iter_mut()
        }

        pub fn as_slice(&self) -> <V as SoAVec<T>>::Slice<'_> {
            self.entries.as_slice()
        }

        pub fn first(&self) -> Option<<V as SoAVec<T>>::Ref<'_>> {
            self.entries.first()
        }

        pub fn get(&self, index: usize) -> Option<<V as SoAVec<T>>::Ref<'_>> {
            self.entries.get(index)
        }

        pub fn last(&self) -> Option<<V as SoAVec<T>>::Ref<'_>> {
            self.entries.last()
        }

        pub fn search_mass(&self, query: MassType, error_tolerance: Tolerance) -> Interval {
            let (lower_bound, upper_bound) = error_tolerance.bounds(query);

            let mut lower_i = self
                .entries
                .mass()
                .partition_point(|entry| *entry <= lower_bound);
            let mut upper_i = self.entries.mass()[lower_i..self.len()]
                .partition_point(|entry| *entry <= upper_bound)
                + lower_i;

            while lower_i > 0 {
                if error_tolerance.test(query, self.entries.mass()[lower_i - 1]) {
                    lower_i -= 1;
                } else {
                    break;
                }
            }

            while upper_i + 1 < self.len() {
                if error_tolerance.test(query, self.entries.mass()[upper_i + 1]) {
                    upper_i += 1;
                } else {
                    break;
                }
            }

            return Interval::new(lower_i, upper_i);
        }

        pub fn search_parent_id(&self, parent_id_range: Interval) -> Interval {
            let mut result = Interval::default();
            if self.is_empty() {
                return result
            }
            let parent_ids = self.entries.parent_id();
            let start_idx = parent_ids.partition_point(|i| {
                parent_id_range.start > *i as usize
            });

            let end_idx = parent_ids.partition_point(|i| {
                parent_id_range.end > *i as usize
            });
            result.start = start_idx;
            result.end = end_idx;
            result
        }

        pub fn select_parent_id(&self, parent_id_range: Interval) -> <V as SoAVec<T>>::Slice<'_> {
            if self.is_empty() {
                return self.as_slice()
            }

            let idx = self.search_parent_id(parent_id_range);

            // let parent_ids = self.entries.parent_id();
            // let slc = Index::index(parent_ids, idx.start..idx.end);
            // let contained = slc.iter().all(|i| parent_id_range.contains(*i as usize));
            // debug_assert!(contained, "{slc:?} {idx:?} {parent_id_range:?}");

            let slc: <V as SoAVec<T>>::Slice<'_> = self.slice(idx);
            return slc
        }

        pub fn min_mass(&self) -> f32 {
            self.min_mass
        }

        pub fn max_mass(&self) -> f32 {
            self.max_mass
        }

        pub fn sort_type(&self) -> SortType {
            self.sort_type
        }

        pub fn slice(&self, index: impl RangeBounds<usize>) -> <V as SoAVec<T>>::Slice<'_> {
            self.entries.slice(index)
        }

        pub fn entries(&self) -> &V {
            &self.entries
        }
    }

    pub struct SoAParentSortedIndexBinSearchIter<
        'a,
        T: StructOfArray,
        V: SoAVec<T> + SoAIndexSortable<T> + 'a,
    >
    where
        for<'t> V::Ref<'t>: IndexSortable,
    {
        bin_iter: V::Iter<'a>,
        parent_range: Interval,
        query: f32,
        error_tolerance: Tolerance,
        spanned: bool,
    }

    impl<'a, T: StructOfArray, V: SoAVec<T> + SoAIndexSortable<T> + 'a> Iterator
        for SoAParentSortedIndexBinSearchIter<'a, T, V>
    where
        for<'t> V::Ref<'t>: IndexSortable,
    {
        type Item = V::Ref<'a>;

        fn next(&mut self) -> Option<Self::Item> {
            self.next_entry()
        }
    }

    impl<'a, T: StructOfArray, V: SoAVec<T> + SoAIndexSortable<T> + 'a>
        SoAParentSortedIndexBinSearchIter<'a, T, V>
    where
        for<'t> V::Ref<'t>: IndexSortable,
    {
        pub fn new(
            bin: &'a SoAIndexBin<T, V>,
            parent_range: Interval,
            query: f32,
            error_tolerance: Tolerance,
        ) -> Self {
            let bin_iter = bin.iter();
            let (lo, hi) = error_tolerance.bounds(query);
            let spanned = lo <= bin.min_mass() && hi >= bin.max_mass();
            Self {
                bin_iter,
                parent_range,
                query,
                error_tolerance,
                spanned,
            }
        }

        pub fn next_entry(&mut self) -> Option<V::Ref<'a>> {
            while let Some(t) = self.bin_iter.next() {
                let parent_id = IndexSortable::parent_id(&t) as usize;
                if self.spanned && self.parent_range.contains(parent_id) {
                    return Some(t);
                }
                if self.error_tolerance.test(self.query, t.mass())
                    && self.parent_range.contains(parent_id)
                {
                    return Some(t);
                }
            }
            None
        }
    }
}

pub use soa_bin::{SoAIndexBin, SoAIndexSortable, SoAIndexSortableSlice, SoAParentSortedIndexBinSearchIter};

#[cfg(test)]
mod test {
    use super::*;
    use crate::parent::Spectrum;

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
        assert!(parent_list.len() == 4);

        let search_out = parent_list.search_mass(2300.01, "5ppm".parse().unwrap());
        assert!(search_out.start == 0);
        assert!(search_out.end == 1);

        let search_out = parent_list.search_mass(2401.0, "5ppm".parse().unwrap());
        assert!(search_out.start == 2);
        assert!(search_out.end == 3);
    }
}
