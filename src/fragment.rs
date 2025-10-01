use std::{error::Error, fmt::Display, str::FromStr};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::sort::{IndexSortable, MassType, ParentID};

#[allow(non_snake_case, non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
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
    Unknown,
}

impl FragmentSeries {
    pub const fn series_name(&self) -> &'static str {
        match self {
            FragmentSeries::b => "b",
            FragmentSeries::y => "y",
            FragmentSeries::c => "c",
            FragmentSeries::z => "z",
            FragmentSeries::a => "a",
            FragmentSeries::x => "x",
            FragmentSeries::Precursor => "Precursor",
            FragmentSeries::PeptideY => "PeptideY",
            FragmentSeries::Oxonium => "Oxonium",
            FragmentSeries::Internal => "Internal",
            FragmentSeries::Unknown => "Unknown",
        }
    }
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
            Self::UnknownSeries(series_label) => {
                format!("Unknown series label \"{}\"", series_label)
            }
            Self::InvalidOrdinal(ordinal_label) => format!(
                "Invalid ordinal value \"{}\", should be an integer",
                ordinal_label
            ),
        };
        f.write_str(&text)
    }
}

impl Error for FragmentSeriesParsingError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FragmentName(pub FragmentSeries, pub u16);

impl FromStr for FragmentSeries {
    type Err = FragmentSeriesParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() == 0 {
            return Err(FragmentSeriesParsingError::Empty);
        }
        let series = match s {
            "b" => FragmentSeries::b,
            "y" => FragmentSeries::y,
            "c" => FragmentSeries::c,
            "z" => FragmentSeries::z,
            "a" => FragmentSeries::a,
            "x" => FragmentSeries::x,
            "Precursor" => FragmentSeries::Precursor,
            "PeptideY" => FragmentSeries::PeptideY,
            "Oxonium" => FragmentSeries::Oxonium,
            "Internal" => FragmentSeries::Internal,
            "Unknown" => FragmentSeries::Unknown,
            _ => {
                return Err(FragmentSeriesParsingError::UnknownSeries(
                    s[0..1].to_string(),
                ))
            }
        };
        Ok(series)
    }
}

impl FromStr for FragmentName {
    type Err = FragmentSeriesParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() == 0 {
            return Err(FragmentSeriesParsingError::Empty);
        }
        let series = match &s[0..1] {
            "b" => FragmentSeries::b,
            "y" => FragmentSeries::y,
            "c" => FragmentSeries::c,
            "z" => FragmentSeries::z,
            "a" => FragmentSeries::a,
            "x" => FragmentSeries::x,
            _ => {
                return Err(FragmentSeriesParsingError::UnknownSeries(
                    s[0..1].to_string(),
                ))
            }
        };
        let ordinal = match s[1..s.len()].parse() {
            Ok(size) => size,
            Err(_) => {
                return Err(FragmentSeriesParsingError::InvalidOrdinal(
                    s[1..].to_string(),
                ))
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
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Fragment {
    pub mass: MassType,
    pub parent_id: ParentID,
    pub series: FragmentSeries,
    pub ordinal: u16,
}

impl IndexSortable for Fragment {
    fn mass(&self) -> MassType {
        self.mass
    }

    fn parent_id(&self) -> ParentID {
        self.parent_id
    }
}

impl Fragment {
    pub fn new(mass: MassType, parent_id: ParentID, series: FragmentSeries, ordinal: u16) -> Self {
        Self {
            mass,
            parent_id,
            series,
            ordinal,
        }
    }
}
