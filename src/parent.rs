#[cfg(feature = "serde")]
use serde::{Serialize, Deserialize};
use soa_derive::StructOfArray;

use crate::sort::{IndexSortable, MassType, ParentID, SoAIndexSortable, SoAIndexSortableSlice};

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
    pub fn new(mass: MassType, id: ParentID, parent_id: ParentID, start_position: u16, size: u16) -> Self {
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
    pub fn new(mass: MassType, id: ParentID, protein_id: ParentID, start_position: u16, sequence: String) -> Self { Self { mass, id, protein_id, start_position, sequence } }
}


impl IndexSortable for Peptide {
    fn mass(&self) -> MassType {
        self.mass
    }

    fn parent_id(&self) -> ParentID {
        self.protein_id as ParentID
    }
}


#[derive(Debug, Clone, Copy, Default, PartialEq, StructOfArray)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[soa_derive(Debug, Clone)]
#[generate_traits]
pub struct Spectrum {
    pub precursor_mass: MassType,
    pub precursor_charge: i32,
    pub source_file_id: ParentID,
    pub scan_number: ParentID,
    pub sort_id: ParentID
}

impl Spectrum {
    pub fn new(
        precursor_mass: MassType,
        precursor_charge: i32,
        source_file_id: ParentID,
        scan_number: ParentID,
        sort_id: ParentID
    ) -> Self {
        Self {
            precursor_mass,
            precursor_charge,
            source_file_id,
            scan_number,
            sort_id
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

impl<'t> IndexSortable for ParentMoleculeRef<'t> {
    fn mass(&self) -> MassType {
        *self.mass
    }

    fn parent_id(&self) -> ParentID {
        *self.source_id
    }
}

impl SoAIndexSortable<ParentMolecule> for ParentMoleculeVec {
    fn mass(&self) -> &[MassType] {
        &self.mass
    }

    fn parent_id(&self) -> &[ParentID] {
        &self.source_id
    }

    fn convert_ref(val_ref: Self::Ref<'_>) -> ParentMolecule {
        val_ref.to_owned()
    }
}

impl<'t> IndexSortable for PeptideRef<'t> {
    fn mass(&self) -> MassType {
        *self.mass
    }

    fn parent_id(&self) -> ParentID {
        *self.protein_id
    }
}

impl SoAIndexSortable<Peptide> for PeptideVec {
    fn mass(&self) -> &[MassType] {
        &self.mass
    }

    fn parent_id(&self) -> &[ParentID] {
        &self.protein_id
    }

    fn convert_ref(val_ref: Self::Ref<'_>) -> Peptide {
        val_ref.to_owned()
    }
}

impl<'t> IndexSortable for SpectrumRef<'t> {
    fn mass(&self) -> MassType {
        *self.precursor_mass
    }

    fn parent_id(&self) -> ParentID {
        *self.source_file_id
    }
}

impl SoAIndexSortable<Spectrum> for SpectrumVec {
    fn mass(&self) -> &[MassType] {
        &self.precursor_mass
    }

    fn parent_id(&self) -> &[ParentID] {
        &self.source_file_id
    }

    fn convert_ref(val_ref: Self::Ref<'_>) -> Spectrum {
        val_ref.to_owned()
    }
}

impl<'t> SoAIndexSortableSlice<'t, Spectrum> for SpectrumSlice<'t> {
    fn mass(&self) -> &[MassType] {
        &self.precursor_mass
    }

    fn parent_id(&self) -> &[ParentID] {
        &self.source_file_id
    }
}