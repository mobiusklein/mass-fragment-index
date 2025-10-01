#![allow(unused)]

use std::collections::HashMap;
use std::fs;
use std::io;
use std::{path::Path, sync::Arc};

use arrow::array::Array;
use arrow::array::DictionaryArray;
use arrow::array::StringArray;
use arrow::array::UInt16Array;
use arrow::array::{
    ArrayRef, AsArray, Float32Array, Float32Builder, Int16Array, Int16Builder, Int32Builder,
    StringBuilder, StringDictionaryBuilder, UInt16Builder, UInt32Array, UInt32Builder, UInt64Array,
    UInt64Builder,
};
use arrow::datatypes::SchemaRef;
use arrow::datatypes::Utf8Type;
use arrow::datatypes::{
    DataType, Field, Float32Type, Int32Type, Schema, UInt16Type, UInt32Type, UInt8Type,
};
use arrow::json::{LineDelimitedWriter as JSONArrayLineWriter, ReaderBuilder as JSONReaderBuilder};
use arrow::record_batch::RecordBatch;

use itertools::izip;

use parquet::arrow::arrow_reader::ArrowReaderBuilder;
use parquet::basic::Compression;
use parquet::basic::ZstdLevel;
use parquet::{arrow::ArrowWriter, file::properties::*};

use super::util::{afield, as_array_ref, field_of, ArrowStorage};
use super::SplitArrowStorage;
use crate::index::SearchIndex;
use crate::sort::IndexBin;
use crate::sort::SortType;
use crate::Fragment;
use crate::Peptide;

pub fn make_fragment_schema() -> Arc<Schema> {
    let mass = afield!("mass", DataType::Float32);
    let parent_id = afield!("parent_id", DataType::UInt32);
    let series = afield!(
        "series",
        DataType::Dictionary(Box::new(DataType::UInt8), Box::new(DataType::Utf8))
    );
    let ordinal = afield!("ordinal", DataType::UInt16);
    let segment_id = afield!("segment_id", DataType::UInt64);
    Arc::new(Schema::new(vec![
        mass, parent_id, series, ordinal, segment_id,
    ]))
}

pub fn make_peptide_schema() -> Arc<Schema> {
    let mass = afield!("mass", DataType::Float32);
    let id = afield!("id", DataType::UInt32);
    let protein_id = afield!("protein_id", DataType::UInt32);
    let start_position = afield!("start_position", DataType::UInt16);
    let sequence = afield!("sequence", DataType::Utf8);

    Arc::new(Schema::new(vec![
        mass,
        id,
        protein_id,
        start_position,
        sequence,
    ]))
}

pub fn make_meta_schema() -> Arc<Schema> {
    let bins_per_dalton = afield!("bins_per_dalton", DataType::UInt32);
    let max_mass = afield!("max_item_mass", DataType::Float32);
    Arc::new(Schema::new(vec![bins_per_dalton, max_mass]))
}

impl ArrowStorage for Peptide {
    fn schema() -> SchemaRef {
        make_peptide_schema()
    }

    fn mass_column() -> Option<usize> {
        Self::schema().column_with_name("mass").map(|(i, field)| i)
    }

    fn parent_id_column() -> Option<usize> {
        Self::schema()
            .column_with_name("protein_id")
            .map(|(i, field)| i)
    }

    fn sort_id_column() -> Option<usize> {
        Self::schema().column_with_name("id").map(|(i, field)| i)
    }

    fn from_batch<'a>(
        batch: &'a RecordBatch,
        schema: SchemaRef,
    ) -> impl Iterator<Item = (Self, u64)> + 'a {
        let mass = batch
            .column_by_name("mass")
            .unwrap()
            .as_primitive::<Float32Type>();
        let start_position = batch
            .column_by_name("start_position")
            .unwrap()
            .as_primitive::<UInt16Type>();
        let protein_id = batch
            .column_by_name("protein_id")
            .unwrap()
            .as_primitive::<UInt32Type>();
        let sequence = batch.column_by_name("sequence").unwrap().as_string::<i32>();
        let id = batch
            .column_by_name("id")
            .unwrap()
            .as_primitive::<UInt32Type>();
        izip!(mass, start_position, protein_id, sequence, id).map(
            |(mass, start_position, protein_id, sequence, id)| {
                (
                    Peptide::new(
                        mass.unwrap(),
                        id.unwrap(),
                        protein_id.unwrap(),
                        start_position.unwrap(),
                        sequence.unwrap().to_string(),
                    ),
                    0,
                )
            },
        )
    }

    fn to_batch(
        batch: &[Self],
        schema: SchemaRef,
        segment_id: u64,
    ) -> Result<RecordBatch, arrow::error::ArrowError> {
        peptide_to_arrow(batch, schema)
    }

    fn archive_name() -> String {
        "peptides.parquet".to_string()
    }

    fn writer_properties() -> WriterPropertiesBuilder {
        WriterProperties::builder()
            .set_column_encoding("mass".into(), parquet::basic::Encoding::BYTE_STREAM_SPLIT)
    }
}

impl ArrowStorage for Fragment {
    fn schema() -> SchemaRef {
        make_fragment_schema()
    }

    fn mass_column() -> Option<usize> {
        Self::schema().column_with_name("mass").map(|(i, field)| i)
    }

    fn parent_id_column() -> Option<usize> {
        Self::schema()
            .column_with_name("parent_id")
            .map(|(i, field)| i)
    }

    fn archive_name() -> String {
        "fragments.parquet".into()
    }

    fn to_batch(
        batch: &[Self],
        schema: SchemaRef,
        segment_id: u64,
    ) -> Result<RecordBatch, arrow::error::ArrowError> {
        fragment_to_arrow(batch, schema, segment_id)
    }

    fn from_batch<'a>(
        batch: &'a RecordBatch,
        schema: SchemaRef,
    ) -> impl Iterator<Item = (Self, u64)> + 'a {
        let mass = field_of!(batch, "mass")
            .as_any()
            .downcast_ref::<Float32Array>()
            .unwrap();
        let ordinal = field_of!(batch, "ordinal")
            .as_any()
            .downcast_ref::<UInt16Array>()
            .unwrap();
        let series = field_of!(batch, "series")
            .as_dictionary::<UInt8Type>()
            .downcast_dict::<StringArray>()
            .unwrap();
        let parent_id = field_of!(batch, "parent_id")
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let segment_id = field_of!(batch, "segment_id")
            .as_any()
            .downcast_ref::<UInt64Array>()
            .unwrap();

        izip!(mass, series, ordinal, parent_id, segment_id).map(
            |(mass, series, ordinal, parent_id, segment_id)| {
                let peak = Fragment::new(
                    mass.unwrap(),
                    parent_id.unwrap(),
                    series.unwrap().parse().unwrap(),
                    ordinal.unwrap() as u16,
                );
                (peak, segment_id.unwrap())
            },
        )
    }

    fn writer_properties() -> WriterPropertiesBuilder {
        WriterProperties::builder()
            .set_column_encoding("mass".into(), parquet::basic::Encoding::BYTE_STREAM_SPLIT)
            .set_column_encoding("segment_id".into(), parquet::basic::Encoding::RLE)
    }
}

pub fn peptide_to_arrow(
    peptides: &[Peptide],
    schema: Arc<Schema>,
) -> Result<RecordBatch, arrow::error::ArrowError> {
    let mut mass_builder = Float32Builder::new();
    let mut sequence_builder = StringBuilder::new();
    let mut protein_id_builder = UInt32Builder::new();
    let mut start_position_builder = UInt16Builder::new();
    let mut id_builder = UInt32Builder::new();

    peptides.iter().for_each(|s| {
        mass_builder.append_value(s.mass);
        id_builder.append_value(s.id);
        protein_id_builder.append_value(s.protein_id);
        sequence_builder.append_value(s.sequence.clone());
        start_position_builder.append_value(s.start_position);
    });

    let columns = vec![
        as_array_ref!(mass_builder),
        as_array_ref!(id_builder),
        as_array_ref!(protein_id_builder),
        as_array_ref!(start_position_builder),
        as_array_ref!(sequence_builder),
    ];

    let batch = RecordBatch::try_new(schema, columns);
    batch
}

pub fn fragment_to_arrow(
    fragments: &[Fragment],
    schema: Arc<Schema>,
    segment_id: u64,
) -> Result<RecordBatch, arrow::error::ArrowError> {
    let mut mass_builder = Float32Builder::new();
    let mut series_builder = StringDictionaryBuilder::<UInt8Type>::new();
    let mut ordinal_builder = UInt16Builder::new();
    let mut parent_id_builder = UInt32Builder::new();
    let mut segment_id_builder = UInt64Builder::new();

    fragments.iter().for_each(|p| {
        mass_builder.append_value(p.mass);
        series_builder.append_value(p.series.series_name());
        ordinal_builder.append_value(p.ordinal);
        parent_id_builder.append_value(p.parent_id);
        segment_id_builder.append_value(segment_id);
    });

    RecordBatch::try_new(
        schema,
        vec![
            as_array_ref!(mass_builder),
            as_array_ref!(parent_id_builder),
            as_array_ref!(series_builder),
            as_array_ref!(ordinal_builder),
            as_array_ref!(segment_id_builder),
        ],
    )
}

pub fn write_fragment_index<P: AsRef<Path>>(
    index: &SearchIndex<Fragment, Peptide>,
    directory: &P,
    compression_level: Option<Compression>,
) -> io::Result<()> {
    let peptides_path = directory.as_ref().join("peptides.parquet");
    let fragments_path = directory.as_ref().join("fragments.parquet");
    let meta_path = directory.as_ref().join("meta.json");

    let peptides_fh = fs::File::create(peptides_path)?;
    let peptides_schema = make_peptide_schema();
    let compression =
        compression_level.unwrap_or(Compression::ZSTD(ZstdLevel::try_new(20).unwrap()));

    {
        let props = WriterProperties::builder()
            .set_compression(compression.clone())
            .set_column_encoding("mass".into(), parquet::basic::Encoding::BYTE_STREAM_SPLIT)
            .build();
        let mut writer =
            ArrowWriter::try_new(peptides_fh, peptides_schema.clone(), Some(props.clone()))?;
        let batch = peptide_to_arrow(index.parents.as_slice(), peptides_schema.clone()).unwrap();
        writer.write(&batch)?;
        writer.close()?;
    }

    let fragments_fh = fs::File::create(fragments_path)?;
    let fragments_schema = make_fragment_schema();

    {
        let props = WriterProperties::builder()
            .set_compression(compression.clone())
            .set_column_encoding("mass".into(), parquet::basic::Encoding::BYTE_STREAM_SPLIT)
            .set_column_encoding("segment_id".into(), parquet::basic::Encoding::RLE)
            .build();
        let mut writer =
            ArrowWriter::try_new(fragments_fh, fragments_schema.clone(), Some(props.clone()))?;
        for (i, bin) in index.bins.iter().enumerate() {
            let batch =
                fragment_to_arrow(bin.as_slice(), fragments_schema.clone(), i as u64).unwrap();
            writer.write(&batch)?;
        }
        writer.close()?;
    }

    let meta_schema = make_meta_schema();
    {
        let meta_fh = io::BufWriter::new(fs::File::create(meta_path)?);
        let bins_per_dalton = UInt32Array::from(vec![index.bins_per_dalton]);
        let max_item_mass = Float32Array::from(vec![index.max_item_mass]);
        let mut writer = JSONArrayLineWriter::new(meta_fh);

        writer
            .write_batches(&vec![&RecordBatch::try_new(
                meta_schema.clone(),
                vec![
                    Arc::new(bins_per_dalton) as ArrayRef,
                    Arc::new(max_item_mass) as ArrayRef,
                ],
            )
            .unwrap()])
            .unwrap();
        writer.finish().unwrap();
    }

    Ok(())
}

pub fn read_fragment_index<P: AsRef<Path>>(
    directory: &P,
) -> io::Result<SearchIndex<Fragment, Peptide>> {
    let peptides_path = directory.as_ref().join("peptides.parquet");
    let fragments_path = directory.as_ref().join("fragments.parquet");
    let meta_path = directory.as_ref().join("meta.json");

    let meta_fh = io::BufReader::new(fs::File::open(meta_path)?);
    let meta_rec = JSONReaderBuilder::new(make_meta_schema())
        .build(meta_fh)
        .unwrap()
        .next()
        .unwrap()
        .unwrap();
    let max_item_mass = meta_rec
        .column_by_name("max_item_mass")
        .unwrap()
        .as_primitive::<Float32Type>()
        .into_iter()
        .flatten()
        .next()
        .unwrap();
    let bins_per_dalton = meta_rec
        .column_by_name("bins_per_dalton")
        .unwrap()
        .as_primitive::<UInt32Type>()
        .into_iter()
        .flatten()
        .next()
        .unwrap();

    let peptides_fh = fs::File::open(peptides_path)?;

    let reader = ArrowReaderBuilder::try_new(peptides_fh)?.build()?;
    let mut peptides: IndexBin<_> = reader
        .map(|b| {
            let b = b.unwrap();
            let mass = b
                .column_by_name("mass")
                .unwrap()
                .as_primitive::<Float32Type>();
            let start_position = b
                .column_by_name("start_position")
                .unwrap()
                .as_primitive::<UInt16Type>();
            let protein_id = b
                .column_by_name("protein_id")
                .unwrap()
                .as_primitive::<UInt32Type>();
            let sequence = b.column_by_name("sequence").unwrap().as_string::<i32>();
            let id = b.column_by_name("id").unwrap().as_primitive::<UInt32Type>();
            let b: Vec<_> = izip!(mass, start_position, protein_id, sequence, id)
                .map(|(mass, start_position, protein_id, sequence, id)| {
                    Peptide::new(
                        mass.unwrap(),
                        id.unwrap(),
                        protein_id.unwrap(),
                        start_position.unwrap(),
                        sequence.unwrap().to_string(),
                    )
                })
                .collect();
            b.into_iter()
        })
        .flatten()
        .collect();

    peptides.sort_type = SortType::ByMass;

    let mut bin_collector: HashMap<u64, Vec<Fragment>> = HashMap::default();
    let fragments_fh = fs::File::open(fragments_path)?;
    let reader = ArrowReaderBuilder::try_new(fragments_fh)?.build()?;

    reader.enumerate().for_each(|(_, b)| {
        let b = b.unwrap();
        let mass = field_of!(b, "mass")
            .as_any()
            .downcast_ref::<Float32Array>()
            .unwrap();
        let ordinal = field_of!(b, "ordinal")
            .as_any()
            .downcast_ref::<UInt16Array>()
            .unwrap();
        let series = field_of!(b, "series")
            .as_dictionary::<UInt8Type>()
            .values()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let parent_id = field_of!(b, "parent_id")
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let segment_id = field_of!(b, "segment_id")
            .as_any()
            .downcast_ref::<UInt64Array>()
            .unwrap();

        izip!(mass, series, ordinal, parent_id, segment_id).for_each(
            |(mass, series, ordinal, parent_id, segment_id)| {
                let peak = Fragment::new(
                    mass.unwrap(),
                    parent_id.unwrap(),
                    series.unwrap().parse().unwrap(),
                    ordinal.unwrap() as u16,
                );
                bin_collector
                    .entry(segment_id.unwrap())
                    .or_default()
                    .push(peak);
            },
        )
    });

    let mut index = SearchIndex::empty(bins_per_dalton, max_item_mass);
    index.parents = peptides;
    for (key, bin) in bin_collector.into_iter() {
        let mut bin = IndexBin::new(bin, SortType::ByParentId, 0.0, 0.0);
        (bin.min_mass, bin.max_mass) = bin.find_min_max_masses();
        index.bins[key as usize] = bin;
    }

    index.sort(SortType::ByParentId);
    Ok(index)
}

impl SplitArrowStorage for Fragment {}
impl SplitArrowStorage for Peptide {}
