use std::iter::IntoIterator;
use std::ops::Range;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Interval {
    pub start: usize,
    pub end: usize,
}

impl Interval {
    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    pub fn contains(&self, point: usize) -> bool {
        self.start <= point && self.end > point
    }

    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }

    pub fn len(&self) -> usize {
        self.end - self.start
    }

    pub fn overlaps(&self, interval: &Interval) -> bool {
        self.end > interval.start && interval.end > self.start
    }
}

impl Into<Range<usize>> for Interval {
    fn into(self) -> Range<usize> {
        Range {
            start: self.start,
            end: self.end,
        }
    }
}

impl From<(usize, usize)> for Interval {
    fn from(value: (usize, usize)) -> Self {
        Self::new(value.0, value.1)
    }
}

impl From<Interval> for (usize, usize) {
    fn from(value: Interval) -> Self {
        (value.start, value.end)
    }
}

impl IntoIterator for Interval {
    type Item = usize;

    type IntoIter = <Range<usize> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        let range: Range<usize> = self.into();
        range.into_iter()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_contains() {
        let iv = Interval::new(5, 10);
        assert!(iv.contains(5));
        assert!(iv.contains(9));
        assert!(!iv.contains(10));
    }
}
