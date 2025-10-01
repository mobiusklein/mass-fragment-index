use std::collections::HashMap;
use std::fs;
use std::io;
use std::{path::Path, sync::Arc};

use arrow::array::Float32Array;
use arrow::array::Int16Array;
use arrow::array::UInt32Array;
use arrow::array::UInt64Array;
use arrow::array::UInt64Builder;
use arrow::array::{ArrayRef, AsArray, Float32Builder, Int16Builder, Int32Builder, UInt32Builder};
use arrow::datatypes::{DataType, Field, Float32Type, Int32Type, Schema, UInt32Type};
use arrow::json::{LineDelimitedWriter as JSONArrayLineWriter, ReaderBuilder as JSONReaderBuilder};
use arrow::record_batch::RecordBatch;

use itertools::izip;

use parquet::arrow::arrow_reader::ArrowReaderBuilder;
use parquet::basic::Compression;
use parquet::basic::Encoding;
use parquet::basic::ZstdLevel;
use parquet::{arrow::ArrowWriter, file::properties::*};

use crate::index::SearchIndex;
use crate::sort::IndexBin;
use crate::sort::SortType;
use crate::{parent::Spectrum, peak::DeconvolutedPeak};

use super::util::{afield, as_array_ref, field_of};
use super::{ArrowStorage, SplitArrowStorage};

pub fn make_peak_schema() -> Arc<Schema> {
    let mass = afield!("mass", DataType::Float32);
    let charge = afield!("charge", DataType::Int16);
    let intensity = afield!("intensity", DataType::Float32);
    let scan_ref = afield!("scan_ref", DataType::UInt32);
    let segment_id = afield!("segment_id", DataType::UInt64);
    Arc::new(Schema::new(vec![
        mass, charge, intensity, scan_ref, segment_id,
    ]))
}

pub fn make_spectrum_schema() -> Arc<Schema> {
    let mass = afield!("precursor_mass", DataType::Float32);
    let charge = afield!("precursor_charge", DataType::Int32);
    let source_file_id = afield!("source_file_id", DataType::UInt32);
    let scan_number = afield!("scan_number", DataType::UInt32);
    let sort_id = afield!("sort_id", DataType::UInt32);

    Arc::new(Schema::new(vec![
        mass,
        charge,
        source_file_id,
        scan_number,
        sort_id,
    ]))
}

pub fn make_meta_schema() -> Arc<Schema> {
    let bins_per_dalton = afield!("bins_per_dalton", DataType::UInt32);
    let max_mass = afield!("max_item_mass", DataType::Float32);
    Arc::new(Schema::new(vec![bins_per_dalton, max_mass]))
}

pub fn spectra_to_arrow(
    spectra: &[Spectrum],
    schema: Arc<Schema>,
) -> Result<RecordBatch, arrow::error::ArrowError> {
    let mut mass_builder = Float32Builder::new();
    let mut charge_builder = Int32Builder::new();
    let mut source_file_id_builder = UInt32Builder::new();
    let mut scan_number_builder = UInt32Builder::new();
    let mut sort_id_builder = UInt32Builder::new();

    spectra.iter().for_each(|s| {
        mass_builder.append_value(s.precursor_mass);
        charge_builder.append_value(s.precursor_charge);
        source_file_id_builder.append_value(s.source_file_id);
        scan_number_builder.append_value(s.scan_number);
        sort_id_builder.append_value(s.sort_id);
    });

    let columns = vec![
        as_array_ref!(mass_builder),
        as_array_ref!(charge_builder),
        as_array_ref!(source_file_id_builder),
        as_array_ref!(scan_number_builder),
        as_array_ref!(sort_id_builder),
    ];

    let batch = RecordBatch::try_new(schema, columns);
    batch
}

pub fn peaks_to_arrow(
    peaks: &[DeconvolutedPeak],
    schema: Arc<Schema>,
    segment_id: u64,
) -> Result<RecordBatch, arrow::error::ArrowError> {
    let mut mass_builder = Float32Builder::new();
    let mut charge_builder = Int16Builder::new();
    let mut intensity_builder = Float32Builder::new();
    let mut scan_ref_builder = UInt32Builder::new();
    let mut segment_id_builder = UInt64Builder::new();

    peaks.iter().for_each(|p| {
        mass_builder.append_value(p.mass);
        charge_builder.append_value(p.charge);
        intensity_builder.append_value(p.intensity);
        scan_ref_builder.append_value(p.scan_ref);
        segment_id_builder.append_value(segment_id);
    });

    RecordBatch::try_new(
        schema,
        vec![
            as_array_ref!(mass_builder),
            as_array_ref!(charge_builder),
            as_array_ref!(intensity_builder),
            as_array_ref!(scan_ref_builder),
            as_array_ref!(segment_id_builder),
        ],
    )
}

impl ArrowStorage for DeconvolutedPeak {
    fn schema() -> arrow::datatypes::SchemaRef {
        make_peak_schema()
    }

    fn mass_column() -> Option<usize> {
        Self::schema().column_with_name("mass").map(|(i, _)| i)
    }

    fn parent_id_column() -> Option<usize> {
        Self::schema().column_with_name("scan_ref").map(|(i, _)| i)
    }

    fn from_batch<'a>(
        batch: &'a RecordBatch,
        _schema: arrow::datatypes::SchemaRef,
    ) -> impl Iterator<Item = (Self, u64)> + 'a {
        let mass = field_of!(batch, "mass")
            .as_any()
            .downcast_ref::<Float32Array>()
            .unwrap();
        let intensity = field_of!(batch, "intensity")
            .as_any()
            .downcast_ref::<Float32Array>()
            .unwrap();
        let charge = field_of!(batch, "charge")
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();
        let scan_ref = field_of!(batch, "scan_ref")
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let segment_id = field_of!(batch, "segment_id")
            .as_any()
            .downcast_ref::<UInt64Array>()
            .unwrap();
        izip!(mass, charge, intensity, scan_ref, segment_id).map(
            |(mass, charge, intensity, scan_ref, segment_id)| {
                let peak = DeconvolutedPeak::new(
                    mass.unwrap(),
                    charge.unwrap() as i16,
                    intensity.unwrap(),
                    scan_ref.unwrap(),
                );
                (peak, segment_id.unwrap())
            },
        )
    }

    fn to_batch(
        batch: &[Self],
        schema: arrow::datatypes::SchemaRef,
        segment_id: u64,
    ) -> Result<RecordBatch, arrow::error::ArrowError> {
        peaks_to_arrow(batch, schema, segment_id)
    }

    fn archive_name() -> String {
        "peaks.parquet".into()
    }

    fn writer_properties() -> WriterPropertiesBuilder {
        WriterProperties::builder()
            .set_column_encoding("segment_id".into(), parquet::basic::Encoding::RLE)
            .set_column_encoding("mass".into(), parquet::basic::Encoding::BYTE_STREAM_SPLIT)
    }
}

impl ArrowStorage for Spectrum {
    fn schema() -> arrow::datatypes::SchemaRef {
        make_spectrum_schema()
    }

    fn mass_column() -> Option<usize> {
        Self::schema()
            .column_with_name("precursor_mass")
            .map(|(i, _)| i)
    }

    fn parent_id_column() -> Option<usize> {
        Self::schema()
            .column_with_name("source_file_id")
            .map(|(i, _)| i)
    }

    fn sort_id_column() -> Option<usize> {
        Self::schema().column_with_name("sort_id").map(|(i, _)| i)
    }

    fn from_batch<'a>(
        batch: &'a RecordBatch,
        _schema: arrow::datatypes::SchemaRef,
    ) -> impl Iterator<Item = (Self, u64)> + 'a {
        let mass = batch
            .column_by_name("precursor_mass")
            .unwrap()
            .as_primitive::<Float32Type>();
        let charge = batch
            .column_by_name("precursor_charge")
            .unwrap()
            .as_primitive::<Int32Type>();
        let source_file_id = batch
            .column_by_name("source_file_id")
            .unwrap()
            .as_primitive::<UInt32Type>();
        let scan_number = batch
            .column_by_name("scan_number")
            .unwrap()
            .as_primitive::<UInt32Type>();
        let sort_id = batch
            .column_by_name("sort_id")
            .unwrap()
            .as_primitive::<UInt32Type>();
        let b: Vec<_> = izip!(mass, charge, source_file_id, scan_number, sort_id)
            .map(
                |(precursor_mass, precursor_charge, source_file_id, scan_number, sort_id)| {
                    Spectrum::new(
                        precursor_mass.unwrap(),
                        precursor_charge.unwrap(),
                        source_file_id.unwrap(),
                        scan_number.unwrap(),
                        sort_id.unwrap(),
                    )
                },
            )
            .collect();
        b.into_iter().map(|p| (p, 0))
    }

    fn to_batch(
        batch: &[Self],
        schema: arrow::datatypes::SchemaRef,
        _segment_id: u64,
    ) -> Result<RecordBatch, arrow::error::ArrowError> {
        spectra_to_arrow(batch, schema)
    }

    fn archive_name() -> String {
        "spectra.parquet".into()
    }

    fn writer_properties() -> WriterPropertiesBuilder {
        WriterProperties::builder().set_column_encoding(
            "precursor_mass".into(),
            parquet::basic::Encoding::BYTE_STREAM_SPLIT,
        )
    }
}

pub fn write_peak_index<P: AsRef<Path>>(
    index: &SearchIndex<DeconvolutedPeak, Spectrum>,
    directory: &P,
    compression_level: Option<Compression>,
) -> io::Result<()> {
    let spectra_path = directory.as_ref().join("spectra.parquet");
    let peaks_path = directory.as_ref().join("peaks.parquet");
    let meta_path = directory.as_ref().join("meta.json");

    let spectra_fh = fs::File::create(spectra_path)?;
    let spectra_schema = make_spectrum_schema();
    let compression =
        compression_level.unwrap_or(Compression::ZSTD(ZstdLevel::try_new(20).unwrap()));

    {
        let props = WriterProperties::builder()
            .set_compression(compression.clone())
            .set_column_encoding("precursor_mass".into(), Encoding::BYTE_STREAM_SPLIT)
            .build();
        let mut writer =
            ArrowWriter::try_new(spectra_fh, spectra_schema.clone(), Some(props.clone()))?;
        let batch = spectra_to_arrow(index.parents.as_slice(), spectra_schema.clone()).unwrap();
        writer.write(&batch)?;
        writer.close()?;
    }

    let peaks_fh = fs::File::create(peaks_path)?;
    let peaks_schema = make_peak_schema();

    {
        let props = WriterProperties::builder()
            .set_compression(compression.clone())
            .set_column_encoding("segment_id".into(), parquet::basic::Encoding::RLE)
            .set_column_encoding("mass".into(), parquet::basic::Encoding::BYTE_STREAM_SPLIT)
            .build();
        let mut writer = ArrowWriter::try_new(peaks_fh, peaks_schema.clone(), Some(props.clone()))?;
        for (i, bin) in index.bins.iter().enumerate() {
            let batch = peaks_to_arrow(bin.as_slice(), peaks_schema.clone(), i as u64).unwrap();
            writer.write(&batch)?;
        }
        writer.close()?;
    }

    let meta_schema = make_meta_schema();
    {
        let meta_fh = io::BufWriter::new(fs::File::create(meta_path)?);
        let mut writer = JSONArrayLineWriter::new(meta_fh);
        let bins_per_dalton = UInt32Array::from(vec![index.bins_per_dalton]);
        let max_item_mass = Float32Array::from(vec![index.max_item_mass]);

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

pub fn read_peak_index<P: AsRef<Path>>(
    directory: &P,
) -> io::Result<SearchIndex<DeconvolutedPeak, Spectrum>> {
    let spectra_path = directory.as_ref().join("spectra.parquet");
    let peaks_path = directory.as_ref().join("peaks.parquet");
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

    let spectra_fh = fs::File::open(spectra_path)?;

    let reader = ArrowReaderBuilder::try_new(spectra_fh)?.build()?;
    let mut spectra: IndexBin<_> = reader
        .map(|b| {
            let b = b.unwrap();
            let mass = b
                .column_by_name("precursor_mass")
                .unwrap()
                .as_primitive::<Float32Type>();
            let charge = b
                .column_by_name("precursor_charge")
                .unwrap()
                .as_primitive::<Int32Type>();
            let source_file_id = b
                .column_by_name("source_file_id")
                .unwrap()
                .as_primitive::<UInt32Type>();
            let scan_number = b
                .column_by_name("scan_number")
                .unwrap()
                .as_primitive::<UInt32Type>();
            let sort_id = b
                .column_by_name("sort_id")
                .unwrap()
                .as_primitive::<UInt32Type>();
            let b: Vec<_> = izip!(mass, charge, source_file_id, scan_number, sort_id)
                .map(
                    |(precursor_mass, precursor_charge, source_file_id, scan_number, sort_id)| {
                        Spectrum::new(
                            precursor_mass.unwrap(),
                            precursor_charge.unwrap(),
                            source_file_id.unwrap(),
                            scan_number.unwrap(),
                            sort_id.unwrap(),
                        )
                    },
                )
                .collect();
            b.into_iter()
        })
        .flatten()
        .collect();
    spectra.sort_type = SortType::ByMass;
    let mut bin_collector: HashMap<u64, Vec<DeconvolutedPeak>> = HashMap::default();
    let peaks_fh = fs::File::open(peaks_path)?;
    let reader = ArrowReaderBuilder::try_new(peaks_fh)?.build()?;
    reader.enumerate().for_each(|(_, b)| {
        let b = b.unwrap();
        let mass = field_of!(b, "mass")
            .as_any()
            .downcast_ref::<Float32Array>()
            .unwrap();
        let intensity = field_of!(b, "intensity")
            .as_any()
            .downcast_ref::<Float32Array>()
            .unwrap();
        let charge = field_of!(b, "charge")
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();
        let scan_ref = field_of!(b, "scan_ref")
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let segment_id = field_of!(b, "segment_id")
            .as_any()
            .downcast_ref::<UInt64Array>()
            .unwrap();
        izip!(mass, charge, intensity, scan_ref, segment_id).for_each(
            |(mass, charge, intensity, scan_ref, segment_id)| {
                let peak = DeconvolutedPeak::new(
                    mass.unwrap(),
                    charge.unwrap() as i16,
                    intensity.unwrap(),
                    scan_ref.unwrap(),
                );
                bin_collector
                    .entry(segment_id.unwrap())
                    .or_default()
                    .push(peak);
            },
        )
    });

    let mut index = SearchIndex::empty(bins_per_dalton, max_item_mass);
    index.parents = spectra;
    for (key, bin) in bin_collector.into_iter() {
        let mut bin = IndexBin::new(bin, SortType::ByParentId, 0.0, 0.0);
        (bin.min_mass, bin.max_mass) = bin.find_min_max_masses();
        index.bins[key as usize] = bin;
    }

    index.sort(SortType::ByParentId);
    Ok(index)
}

impl SplitArrowStorage for DeconvolutedPeak {}
impl SplitArrowStorage for Spectrum {}
