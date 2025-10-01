use std::{collections::HashMap, fs, io, path::Path, sync::Arc};

use arrow::{
    array::{ArrayRef, AsArray, Float32Array, RecordBatch, UInt32Array},
    datatypes::{DataType, Field, Float32Type, Schema, SchemaRef, UInt32Type},
    error::ArrowError,
    json::{LineDelimitedWriter, ReaderBuilder as JSONReaderBuilder},
};
use parquet::{
    arrow::{arrow_reader::ArrowReaderBuilder, ArrowWriter},
    basic::{Compression, ZstdLevel},
    file::properties::{WriterProperties, WriterPropertiesBuilder},
};

use crate::MassType;

pub trait ArrowStorage: Sized {
    fn schema() -> SchemaRef;

    fn from_batch<'a>(
        batch: &'a RecordBatch,
        schema: SchemaRef,
    ) -> impl Iterator<Item = (Self, u64)> + 'a;

    fn to_batch(
        batch: &[Self],
        schema: SchemaRef,
        segment_id: u64,
    ) -> Result<RecordBatch, ArrowError>;

    fn archive_name() -> String;

    fn writer_properties() -> WriterPropertiesBuilder;

    fn mass_column() -> Option<usize> {
        None
    }

    fn parent_id_column() -> Option<usize> {
        None
    }

    fn sort_id_column() -> Option<usize> {
        None
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub struct IndexMetadata {
    pub bins_per_dalton: u32,
    pub max_item_mass: MassType,
}

macro_rules! afield {
    ($name:expr, $ctype:expr) => {
        Arc::new(Field::new($name, $ctype, false))
    };
}

macro_rules! as_array_ref {
    ($a:expr) => {
        Arc::new($a.finish()) as ArrayRef
    };
}

macro_rules! field_of {
    ($batch:expr, $name:expr) => {
        $batch.column_by_name($name).unwrap()
    };
}

pub(crate) use afield;
pub(crate) use as_array_ref;
pub(crate) use field_of;

impl ArrowStorage for IndexMetadata {
    fn schema() -> SchemaRef {
        let bins_per_dalton = afield!("bins_per_dalton", DataType::UInt32);
        let max_mass = afield!("max_item_mass", DataType::Float32);
        Arc::new(Schema::new(vec![bins_per_dalton, max_mass]))
    }

    fn from_batch<'a>(
        batch: &'a RecordBatch,
        _schema: SchemaRef,
    ) -> impl Iterator<Item = (Self, u64)> + 'a {
        let max_item_mass = batch
            .column_by_name("max_item_mass")
            .unwrap()
            .as_primitive::<Float32Type>()
            .into_iter()
            .flatten()
            .next()
            .unwrap();
        let bins_per_dalton = batch
            .column_by_name("bins_per_dalton")
            .unwrap()
            .as_primitive::<UInt32Type>()
            .into_iter()
            .flatten()
            .next()
            .unwrap();
        let this = Self {
            max_item_mass,
            bins_per_dalton,
        };
        [(this, 0)].into_iter()
    }

    fn to_batch(
        batch: &[Self],
        schema: SchemaRef,
        _segment_id: u64,
    ) -> Result<RecordBatch, ArrowError> {
        let this = batch.first().unwrap();
        let bins_per_dalton = UInt32Array::from(vec![this.bins_per_dalton]);
        let max_item_mass = Float32Array::from(vec![this.max_item_mass]);
        RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(bins_per_dalton) as ArrayRef,
                Arc::new(max_item_mass) as ArrayRef,
            ],
        )
    }

    fn archive_name() -> String {
        "meta.json".into()
    }

    fn writer_properties() -> WriterPropertiesBuilder {
        WriterProperties::builder()
    }
}

pub trait IndexBinaryStorage<'a, T: ArrowStorage + 'a, P: ArrowStorage, M: ArrowStorage> {
    fn write_metadata(&self, directory: &Path) -> io::Result<()> {
        let metadata = self.to_metadata();
        let meta_path = directory.join(M::archive_name());
        let meta_schema = M::schema();
        let meta_fh = io::BufWriter::new(fs::File::create(meta_path)?);
        let mut writer = LineDelimitedWriter::new(meta_fh);
        let metadata = M::to_batch(&[metadata], meta_schema, 0).unwrap();

        writer.write(&metadata).unwrap();
        writer.finish().unwrap();
        Ok(())
    }

    fn write_parents(&self, directory: &Path, compression_level: &Compression) -> io::Result<()> {
        let parent_path = directory.join(P::archive_name());
        let parent_schema = P::schema();
        let props = P::writer_properties()
            .set_compression(compression_level.clone())
            .build();
        let mut writer = ArrowWriter::try_new(
            fs::File::create(parent_path)?,
            parent_schema.clone(),
            Some(props),
        )?;
        let batch = P::to_batch(self.parents(), parent_schema.clone(), 0).unwrap();
        writer.write(&batch)?;
        writer.close()?;
        Ok(())
    }

    fn write_entries(
        &'a self,
        directory: &Path,
        compression_level: &Compression,
    ) -> io::Result<()> {
        let entries_path = directory.join(T::archive_name());
        let entries_schema = T::schema();
        let props = T::writer_properties()
            .set_compression(compression_level.clone())
            .build();
        let mut writer = ArrowWriter::try_new(
            fs::File::create(entries_path)?,
            entries_schema.clone(),
            Some(props),
        )?;
        for (i, bin) in self.iter_entries().enumerate() {
            let batch = T::to_batch(bin, entries_schema.clone(), i as u64).unwrap();
            writer.write(&batch)?;
        }
        writer.close()?;
        Ok(())
    }

    fn write<D: AsRef<Path>>(
        &'a self,
        directory: &D,
        compression_level: Option<Compression>,
    ) -> io::Result<()> {
        let directory = directory.as_ref();

        let compression_level =
            compression_level.unwrap_or_else(|| Compression::ZSTD(ZstdLevel::try_new(9).unwrap()));

        self.write_metadata(directory)?;
        self.write_parents(directory, &compression_level)?;
        self.write_entries(directory, &compression_level)?;

        Ok(())
    }

    fn parents(&self) -> &[P];

    fn iter_entries(&'a self) -> impl Iterator<Item = &'a [T]> + 'a;

    fn to_metadata(&self) -> M;

    fn from_components(metadata: M, parents: Vec<P>, entries: HashMap<u64, Vec<T>>) -> Self;

    fn read<D: AsRef<Path>>(directory: &D) -> io::Result<Self>
    where
        Self: Sized,
    {
        let parents_path = directory.as_ref().join(P::archive_name());
        let entries_path = directory.as_ref().join(T::archive_name());
        let meta_path = directory.as_ref().join(M::archive_name());

        let metadata = {
            let meta_schema = M::schema();
            let meta_fh = io::BufReader::new(fs::File::open(meta_path)?);
            let meta_rec = JSONReaderBuilder::new(meta_schema.clone())
                .build(meta_fh)
                .unwrap()
                .next()
                .unwrap()
                .unwrap();

            let (metadata, _) = M::from_batch(&meta_rec, meta_schema.clone())
                .next()
                .unwrap();
            metadata
        };

        let parents = {
            let parent_schema = P::schema();
            let parents_fh = fs::File::open(parents_path)?;

            let reader = ArrowReaderBuilder::try_new(parents_fh)?.build()?;
            let mut parents = Vec::new();
            for batch in reader {
                parents
                    .extend(P::from_batch(&batch.unwrap(), parent_schema.clone()).map(|(p, _)| p));
            }

            parents
        };

        let entries = {
            let mut bin_collector: HashMap<u64, Vec<T>> = HashMap::default();
            let entries_fh = fs::File::open(entries_path)?;
            let reader = ArrowReaderBuilder::try_new(entries_fh)?.build()?;
            let entry_schema = T::schema();

            for batch in reader {
                for (entry, segment_id) in T::from_batch(&batch.unwrap(), entry_schema.clone()) {
                    bin_collector.entry(segment_id).or_default().push(entry);
                }
            }

            bin_collector
        };

        let this = Self::from_components(metadata, parents, entries);
        Ok(this)
    }
}
