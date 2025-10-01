use std::collections::hash_map::Entry;
use std::collections::HashMap;

use crate::fragment::{Fragment, FragmentSeries};
use crate::index::SearchIndex;
use crate::sort::{IndexSortable, MassType, ParentID, Tolerance};

#[derive(Debug, Clone, Copy)]
pub struct IndexMatcher {
    pub tolerance: Tolerance,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct HyperscoreMatcher {
    pub nt_intensity: f32,
    pub ct_intensity: f32,
    pub total: f32,
    pub n_nt_hits: u16,
    pub n_ct_hits: u16,
}

fn factorial(value: u16) -> f64 {
    let mut acc = value as f64;
    let mut next = acc - 1.0;

    while next > 0.0 {
        acc *= next;
        next -= 1.0;
    }
    acc
}

impl From<HyperscoreMatcher> for f64 {
    fn from(value: HyperscoreMatcher) -> Self {
        let score = (value.nt_intensity as f64).ln()
            + factorial(value.n_nt_hits)
            + (value.ct_intensity as f64).ln()
            + factorial(value.n_ct_hits);
        score
    }
}

impl HyperscoreMatcher {
    pub fn add(&mut self, fragment: Fragment, intensity: f32) {
        match fragment.series {
            FragmentSeries::a | FragmentSeries::b | FragmentSeries::c => {
                self.nt_intensity += intensity;
                self.n_nt_hits += 1;
            }
            FragmentSeries::x | FragmentSeries::y | FragmentSeries::z => {
                self.ct_intensity += intensity;
                self.n_ct_hits += 1;
            }
            _ => {}
        }
    }
}

impl IndexMatcher {
    pub fn new(tolerance: Tolerance) -> Self {
        Self { tolerance }
    }

    pub fn search_index<'a, T: IndexSortable + Default + Clone, P: IndexSortable + Default>(
        &self,
        queries: &[MassType],
        parent_range: (MassType, MassType),
        index: &'a SearchIndex<T, P>,
    ) -> HashMap<ParentID, Vec<&'a T>> {
        let parent_interval =
            index.parents_for_range(parent_range.0, parent_range.1, self.tolerance);

        let mut matches: HashMap<ParentID, Vec<&T>> = HashMap::new();

        for query in queries {
            for hit in index.search(*query, self.tolerance, Some(parent_interval)) {
                match matches.entry(hit.parent_id()) {
                    Entry::Occupied(mut entry) => {
                        entry.get_mut().push(hit);
                    }
                    Entry::Vacant(entry) => {
                        entry.insert(vec![hit]);
                    }
                }
            }
        }
        matches
    }

    pub fn search_index_hyperscore<P: IndexSortable + Default>(
        &self,
        queries: &[(MassType, f32)],
        parent_range: (MassType, MassType),
        index: &SearchIndex<Fragment, P>,
    ) -> HashMap<ParentID, f64> {
        let parent_interval =
            index.parents_for_range(parent_range.0, parent_range.1, self.tolerance);

        let mut matches: HashMap<ParentID, HyperscoreMatcher> = HashMap::new();
        for (mass, intensity) in queries {
            for hit in index.search(*mass, self.tolerance, Some(parent_interval)) {
                match matches.entry(hit.parent_id()) {
                    Entry::Occupied(mut entry) => {
                        entry.get_mut().add(*hit, *intensity);
                    }
                    Entry::Vacant(entry) => {
                        let mut state = HyperscoreMatcher::default();
                        state.add(*hit, *intensity);
                        entry.insert(state);
                    }
                }
            }
        }

        matches.into_iter().map(|(k, v)| (k, v.into())).collect()
    }
}
