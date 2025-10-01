#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::sort::{IndexSortable, MassType, ParentID};

#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ParentMolecule {
    pub mass: MassType,
    pub id: ParentID,
    pub source_id: ParentID,
    pub start_position: u16,
    pub size: u16,
}

impl ParentMolecule {
    pub fn new(
        mass: MassType,
        id: ParentID,
        parent_id: ParentID,
        start_position: u16,
        size: u16,
    ) -> Self {
        Self {
            mass,
            id,
            source_id: parent_id,
            start_position,
            size,
        }
    }
}

impl IndexSortable for ParentMolecule {
    fn mass(&self) -> MassType {
        self.mass
    }

    fn parent_id(&self) -> ParentID {
        self.source_id
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Peptide {
    pub mass: MassType,
    pub id: ParentID,
    pub protein_id: ParentID,
    pub start_position: u16,
    pub sequence: String,
}

impl Peptide {
    pub fn new(
        mass: MassType,
        id: ParentID,
        protein_id: ParentID,
        start_position: u16,
        sequence: String,
    ) -> Self {
        Self {
            mass,
            id,
            protein_id,
            start_position,
            sequence,
        }
    }
}

impl IndexSortable for Peptide {
    fn mass(&self) -> MassType {
        self.mass
    }

    fn parent_id(&self) -> ParentID {
        self.protein_id as ParentID
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Spectrum {
    pub precursor_mass: MassType,
    pub precursor_charge: i32,
    pub source_file_id: ParentID,
    pub scan_number: ParentID,
    pub sort_id: ParentID,
}

impl Spectrum {
    pub fn new(
        precursor_mass: MassType,
        precursor_charge: i32,
        source_file_id: ParentID,
        scan_number: ParentID,
        sort_id: ParentID,
    ) -> Self {
        Self {
            precursor_mass,
            precursor_charge,
            source_file_id,
            scan_number,
            sort_id,
        }
    }
}

impl IndexSortable for Spectrum {
    fn mass(&self) -> MassType {
        self.precursor_mass
    }

    fn parent_id(&self) -> ParentID {
        self.source_file_id
    }
}
