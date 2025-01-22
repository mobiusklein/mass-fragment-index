
#[macro_export]
macro_rules! generate_index_sortable {
    ($tname:ty, $mass_name:tt, $parent_id:tt) => {
        impl $crate::IndexSortable for $tname {
            fn mass(&self) -> $crate::MassType {
                self.$mass_name
            }

            fn parent_id(&self) -> $crate::ParentID {
                self.$parent_id
            }
        }
    };

    ($tname:ty, $mass_name:tt, $parent_id:tt, $vec_name:ty, $ref_name:ty) => {
        impl $crate::IndexSortable for $tname {
            fn mass(&self) -> $crate::MassType {
                self.$mass_name
            }

            fn parent_id(&self) -> $crate::ParentID {
                self.$parent_id
            }
        }

        impl $crate::SoAIndexSortable<$tname> for $vec_name {
            fn mass(&self) -> &[$crate::MassType] {
                &self.$mass_name
            }

            fn parent_id(&self) -> &[$crate::ParentID] {
                &self.$parent_id
            }

            fn convert_ref(val_ref: Self::Ref<'_>) -> $tname {
                val_ref.to_owned()
            }
        }

        impl<'t> $crate::IndexSortable for $ref_name {
            fn mass(&self) -> $crate::MassType {
                *self.$mass_name
            }

            fn parent_id(&self) -> $crate::ParentID {
                *self.$parent_id
            }
        }
    };
}
