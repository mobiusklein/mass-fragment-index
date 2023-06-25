
use std::{str::FromStr, error::Error, fmt::Display};

use crate::sort::{IndexSortable};

#[allow(non_snake_case, non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FragmentSeries {
    b,
    y,
    c,
    z,
    a,
    x,
    Precursor,
    PeptideY,
    Oxonium,
    Internal,
    Unknown
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum FragmentSeriesParsingError {
    Empty,
    UnknownSeries(String),
    InvalidOrdinal(String),
}

impl Display for FragmentSeriesParsingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let text = match &self {
            Self::Empty => "Fragment name cannot be an empty string".to_string(),
            Self::UnknownSeries(series_label) => format!("Unknown series label \"{}\"", series_label),
            Self::InvalidOrdinal(ordinal_label) => format!("Invalid ordinal value \"{}\", should be an integer", ordinal_label),
        };
        f.write_str(&text)
    }
}

impl Error for FragmentSeriesParsingError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FragmentName(pub FragmentSeries, pub u16);

impl FromStr for FragmentName {
    type Err = FragmentSeriesParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() == 0 {
            return Err(FragmentSeriesParsingError::Empty)
        }
        let series = match &s[0..1] {
            "b" => FragmentSeries::b,
            "y" => FragmentSeries::y,
            "c" => FragmentSeries::c,
            "z" => FragmentSeries::z,
            "a" => FragmentSeries::a,
            "x" => FragmentSeries::x,
            _ => {
                return Err(FragmentSeriesParsingError::UnknownSeries(s[0..1].to_string()))
            }
        };
        let ordinal = match s[1..s.len()].parse() {
            Ok(size) => size,
            Err(_) => {
                return Err(FragmentSeriesParsingError::InvalidOrdinal(s[1..].to_string()))
            }
        };
        Ok(FragmentName(series, ordinal))
    }
}




impl Default for FragmentSeries {
    fn default() -> Self {
        Self::Unknown
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Fragment {
    pub mass: f32,
    pub parent_id: usize,
    pub series: FragmentSeries,
    pub ordinal: u16,
}


impl IndexSortable for Fragment {
    fn mass(&self) -> f32 {
        self.mass
    }

    fn parent_id(&self) -> usize {
        self.parent_id
    }
}

impl Fragment {
    pub fn new(mass: f32, parent_id: usize, series: FragmentSeries, ordinal: u16) -> Self {
        Self {
            mass,
            parent_id,
            series,
            ordinal,
        }
    }
}