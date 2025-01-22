#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use soa_derive::StructOfArray;

use crate::sort::{MassType, ParentID};

#[derive(Debug, Clone, Copy, PartialEq, Default, StructOfArray)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[soa_derive(Debug, Clone)]
#[generate_traits]
pub struct DeconvolutedPeak {
    pub mass: MassType,
    pub charge: i16,
    pub intensity: f32,
    pub scan_ref: ParentID,
}

impl PartialOrd for DeconvolutedPeak {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.mass.partial_cmp(&other.mass)
    }
}

impl DeconvolutedPeak {
    pub fn new(mass: MassType, charge: i16, intensity: f32, scan_ref: ParentID) -> Self {
        Self {
            mass,
            charge,
            intensity,
            scan_ref,
        }
    }
}

crate::generate_index_sortable!(
    DeconvolutedPeak,
    mass,
    scan_ref,
    DeconvolutedPeakVec,
    DeconvolutedPeakRef<'t>
);

#[derive(Debug, Clone, Copy, PartialEq, Default, StructOfArray)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[soa_derive(Debug, Clone)]
#[generate_traits]
pub struct MZPeak {
    pub mz: MassType,
    pub intensity: f32,
    pub scan_ref: ParentID,
}

impl PartialOrd for MZPeak {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.mz.partial_cmp(&other.mz)
    }
}

impl MZPeak {
    pub fn new(mz: MassType, intensity: f32, scan_ref: ParentID) -> Self {
        Self {
            mz,
            intensity,
            scan_ref,
        }
    }
}

crate::generate_index_sortable!(MZPeak, mz, scan_ref, MZPeakVec, MZPeakRef<'t>);

#[cfg(test)]
mod test {
    use super::*;
    use crate::IndexSortable;

    #[test]
    fn test_creation() {
        let peak = DeconvolutedPeak::new(256.03, 1, 0.0, 300);
        assert!(peak.mass == 256.03);
        assert!(peak.mass() == 256.03);
        assert!(peak.scan_ref == 300);
        assert!(peak.parent_id() == 300);
    }
}
