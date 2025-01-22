#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use soa_derive::StructOfArray;

use crate::sort::{MassType, ParentID};

#[derive(Debug, Clone, Copy, Default, PartialEq, StructOfArray)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[soa_derive(Debug, Clone)]
#[generate_traits]
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

crate::generate_index_sortable!(
    ParentMolecule,
    mass,
    source_id,
    ParentMoleculeVec,
    ParentMoleculeRef<'t>
);

#[derive(Debug, Clone, Default, PartialEq, StructOfArray)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[soa_derive(Debug, Clone)]
#[generate_traits]
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

crate::generate_index_sortable!(Peptide, mass, protein_id, PeptideVec, PeptideRef<'t>);

#[derive(Debug, Clone, Copy, Default, PartialEq, StructOfArray)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[soa_derive(Debug, Clone)]
#[generate_traits]
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

crate::generate_index_sortable!(
    Spectrum,
    precursor_mass,
    source_file_id,
    SpectrumVec,
    SpectrumRef<'t>
);
