use std::{
    collections::{HashMap, HashSet},
    fs, io,
    marker::PhantomData,
    ops::Range,
    path::Path,
    sync::Arc,
};

use arrow::{
    array::{
        ArrayRef, AsArray, Float32Array, Int16Array, Int32Array, RecordBatch, StringArray, StringDictionaryBuilder, StructArray, UInt16Array, UInt32Array, UInt64Array
    },
    datatypes::{SchemaRef, UInt8Type},
    json::{LineDelimitedWriter, ReaderBuilder},
};
use parquet::{
    arrow::{arrow_reader::ArrowReaderBuilder, ArrowWriter},
    basic::{Compression, ZstdLevel},
};
use soa_derive::prelude::*;

use crate::{
    fragment::FragmentVec,
    parent::{PeptideVec, Spectrum, SpectrumVec},
    peak::{DeconvolutedPeak, DeconvolutedPeakVec},
    sort::{ParentID, SoAIndexBin, SoAIndexSortable},
    storage::{ArrowStorage, SplitArrowStorage},
    Fragment, IndexSortable, Interval, MassType, Peptide,
};

use super::{BinStorageStrategy, SplitBand, SplitStorageOptions};

macro_rules! col_to_vec {
    ($batch:ident, $col:literal, $arr:ty) => {{
        let x: &$arr = $batch.column($col).as_any().downcast_ref().unwrap();
        x.values().to_vec()
    }};
}

pub trait SoAArrowStorage<V: SoAVec<Self> + SoAAppendVec<Self>> : ArrowStorage + StructOfArray {
    fn to_batch_soa<'a>(
        batch: <V as SoAVec<Self>>::Slice<'a>,
        schema: SchemaRef,
        segment_id: u64,
    ) -> RecordBatch;

    fn from_batch_soa(batch: &RecordBatch, schema: SchemaRef) -> (V, Vec<u64>);

    fn partition_by_segments(batch: V, segments: Vec<u64>) -> HashMap<u64, V> {
        let mut segment_map = HashMap::default();
        if segments.is_empty() {
            return segment_map;
        }

        let mut start = 0;
        let mut current_segment: Option<u64> = None;
        let n_segs = segments.len();

        for (i, seg_id) in segments.into_iter().enumerate() {
            if Some(seg_id) == current_segment {
                continue;
            } else {
                match current_segment {
                    Some(active_seg) => {
                        let chunk = batch.slice(start..i);
                        segment_map
                            .entry(active_seg)
                            .or_insert_with(|| V::new())
                            .extend_from_slice(chunk);
                        start = i;
                        current_segment = Some(seg_id);
                    }
                    None => {
                        current_segment = Some(seg_id);
                    }
                }
            }
        }

        if let Some(current_segment) = current_segment {
            let chunk = batch.slice(start..n_segs);
            segment_map
                .entry(current_segment)
                .or_insert_with(|| V::new())
                .extend_from_slice(chunk);
        }

        segment_map
    }
}

pub trait SoAIndexBinaryStorage<
    'a,
    T: SoAArrowStorage<TV> + 'a,
    TV: SoAVec<T> + SoAIndexSortable<T> + SoAAppendVec<T>,
    P: SoAArrowStorage<PV>,
    PV: SoAVec<P> + SoAIndexSortable<P> + SoAAppendVec<P>,
    M: ArrowStorage,
> where
    for<'t> <TV as soa_derive::SoAVec<T>>::Ref<'t>: IndexSortable,
    for<'t> <PV as soa_derive::SoAVec<P>>::Ref<'t>: IndexSortable,
{
    fn parents(&self) -> &PV;

    fn iter_entries<'b>(&'b self) -> impl Iterator<Item = &'b SoAIndexBin<T, TV>> where T: 'b, TV: 'b;

    fn to_metadata(&self) -> M;

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
        let batch = P::to_batch_soa(self.parents().as_slice(), parent_schema.clone(), 0);
        writer.write(&batch)?;
        writer.close()?;
        Ok(())
    }

    fn write_entries(&self, directory: &Path, compression_level: &Compression) -> io::Result<()> {
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
            let batch = T::to_batch_soa(bin.as_slice(), entries_schema.clone(), i as u64);
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

    fn from_components(metadata: M, parents: PV, entries: HashMap<u64, TV>) -> Self;

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
            let meta_rec = ReaderBuilder::new(meta_schema.clone())
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
            let mut parents = PV::new();
            for batch in reader {
                let (mut p, _) = P::from_batch_soa(&batch.unwrap(), parent_schema.clone());
                parents.append(&mut p);
            }

            parents
        };

        let entries = {
            let mut bin_collector: HashMap<u64, TV> = HashMap::default();
            let entries_fh = fs::File::open(entries_path)?;
            let reader = ArrowReaderBuilder::try_new(entries_fh)?.build()?;
            let entry_schema = T::schema();

            for batch in reader {
                let (entries, segments) = T::from_batch_soa(&batch.unwrap(), entry_schema.clone());
                let mut parts: HashMap<u64, TV> =
                    T::partition_by_segments(entries, segments);
                for (k, v) in parts.iter_mut() {
                    bin_collector
                        .entry(*k)
                        .or_insert_with(|| TV::new())
                        .append(v);
                }
            }

            bin_collector
        };

        let this = Self::from_components(metadata, parents, entries);
        Ok(this)
    }
}

impl SoAArrowStorage<PeptideVec> for Peptide {
    fn to_batch_soa<'a>(
        batch: <PeptideVec as SoAVec<Self>>::Slice<'a>,
        schema: SchemaRef,
        _segment_id: u64,
    ) -> RecordBatch {
        let masses = Float32Array::from(batch.mass.to_vec());
        let protein_ids = UInt32Array::from(batch.protein_id.to_vec());
        let start_positions = UInt16Array::from(batch.start_position.to_vec());
        let ids = UInt32Array::from(batch.id.to_vec());
        let sequences = StringArray::from(batch.sequence.to_vec());

        let columns = vec![
            Arc::new(masses) as ArrayRef,
            Arc::new(ids) as ArrayRef,
            Arc::new(protein_ids),
            Arc::new(start_positions),
            Arc::new(sequences),
        ];

        let batch = RecordBatch::try_new(schema, columns);
        batch.unwrap()
    }

    fn from_batch_soa(batch: &RecordBatch, _schema: SchemaRef) -> (PeptideVec, Vec<u64>) {
        let mut items = PeptideVec::new();

        let masses: &Float32Array = batch.column(0).as_any().downcast_ref().unwrap();
        items.mass = masses.values().to_vec();

        let ids: &UInt32Array = batch.column(1).as_any().downcast_ref().unwrap();
        items.id = ids.values().to_vec();

        let protein_ids: &UInt32Array = batch.column(2).as_any().downcast_ref().unwrap();
        items.protein_id = protein_ids.values().to_vec();

        let start_positions: &UInt16Array = batch.column(3).as_any().downcast_ref().unwrap();
        items.start_position = start_positions.values().to_vec();

        let sequences: &StringArray = batch.column(4).as_string();
        items.sequence = sequences.iter().map(|s| s.unwrap().to_string()).collect();

        let segments = Vec::new();

        (items, segments)
    }
}

impl SoAArrowStorage<FragmentVec> for Fragment {

    fn to_batch_soa<'a>(
        batch: <FragmentVec as soa_derive::SoAVec<Fragment>>::Slice<'a>,
        schema: SchemaRef,
        segment_id: u64,
    ) -> RecordBatch {
        let masses = Float32Array::from(batch.mass.to_vec());
        let parent_ids = UInt32Array::from(batch.parent_id.to_vec());
        let mut series_builder = StringDictionaryBuilder::<UInt8Type>::new();
        series_builder.extend(batch.series.iter().map(|s| Some(s.series_name())));
        let series = series_builder.finish();
        let ordinals = UInt16Array::from(batch.ordinal.to_vec());
        let segment_ids = UInt64Array::from_value(segment_id, batch.len());

        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(masses) as ArrayRef,
                Arc::new(parent_ids),
                Arc::new(series),
                Arc::new(ordinals),
                Arc::new(segment_ids),
            ],
        )
        .unwrap()
    }

    fn from_batch_soa(batch: &RecordBatch, _schema: SchemaRef) -> (FragmentVec, Vec<u64>) {
        let mut items = FragmentVec::new();

        let segments = col_to_vec!(batch, 4, UInt64Array);

        let series_strs = batch
            .column(2)
            .as_dictionary::<UInt8Type>()
            .downcast_dict::<StringArray>()
            .unwrap();

        items.mass = col_to_vec!(batch, 0, Float32Array);
        items.parent_id = col_to_vec!(batch, 1, UInt32Array);
        items.series = series_strs
            .into_iter()
            .map(|s| s.unwrap().parse().unwrap())
            .collect();
        items.ordinal = col_to_vec!(batch, 3, UInt16Array);

        (items, segments)
    }
}

impl SoAArrowStorage<SpectrumVec> for Spectrum {
    fn to_batch_soa<'a>(
        batch: <SpectrumVec as SoAVec<Self>>::Slice<'a>,
        schema: SchemaRef,
        _segment_id: u64,
    ) -> RecordBatch {
        let mass = Float32Array::from(batch.precursor_mass.to_vec());
        let charge = Int32Array::from(batch.precursor_charge.to_vec());
        let source_file_id = UInt32Array::from(batch.source_file_id.to_vec());
        let scan_number = UInt32Array::from(batch.scan_number.to_vec());
        let sort_id = UInt32Array::from(batch.sort_id.to_vec());

        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(mass) as ArrayRef,
                Arc::new(charge),
                Arc::new(source_file_id),
                Arc::new(scan_number),
                Arc::new(sort_id),
            ],
        )
        .unwrap()
    }

    fn from_batch_soa(batch: &RecordBatch, _schema: SchemaRef) -> (SpectrumVec, Vec<u64>) {
        let mut items = SpectrumVec::new();
        items.precursor_mass = col_to_vec!(batch, 0, Float32Array);
        items.precursor_charge = col_to_vec!(batch, 1, Int32Array);
        items.source_file_id = col_to_vec!(batch, 2, UInt32Array);
        items.scan_number = col_to_vec!(batch, 3, UInt32Array);
        items.sort_id = col_to_vec!(batch, 4, UInt32Array);

        (items, Vec::new())
    }
}

impl SoAArrowStorage<DeconvolutedPeakVec> for DeconvolutedPeak {
    fn to_batch_soa<'a>(
        batch: <DeconvolutedPeakVec as SoAVec<Self>>::Slice<'a>,
        schema: SchemaRef,
        segment_id: u64,
    ) -> RecordBatch {
        let mass = Float32Array::from(batch.mass.to_vec());
        let charge = Int16Array::from(batch.charge.to_vec());
        let intensity = Float32Array::from(batch.intensity.to_vec());
        let scan_ref = UInt32Array::from(batch.scan_ref.to_vec());
        let segment_ids = UInt64Array::from_value(segment_id, batch.len());

        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(mass) as ArrayRef,
                Arc::new(charge),
                Arc::new(intensity),
                Arc::new(scan_ref),
                Arc::new(segment_ids),
            ],
        )
        .unwrap()
    }

    fn from_batch_soa(batch: &RecordBatch, _schema: SchemaRef) -> (DeconvolutedPeakVec, Vec<u64>) {
        let mut this = DeconvolutedPeakVec::new();
        this.mass = col_to_vec!(batch, 0, Float32Array);
        this.charge = col_to_vec!(batch, 1, Int16Array);
        this.intensity = col_to_vec!(batch, 2, Float32Array);
        this.scan_ref = col_to_vec!(batch, 3, UInt32Array);

        let segment_ids = col_to_vec!(batch, 4, UInt64Array);
        (this, segment_ids)
    }
}

pub struct SoASplitIndexBinaryStorageWriter<
    'a,
    T: ArrowStorage + 'a + IndexSortable + Clone + SplitArrowStorage + SoAArrowStorage<TV>,
    TV: SoAVec<T> + SoAIndexSortable<T> + SoAAppendVec<T>,
    P: ArrowStorage + IndexSortable + SoAArrowStorage<PV>,
    PV: SoAVec<P> + SoAIndexSortable<P> + SoAAppendVec<P>,
    M: ArrowStorage,
    I: SoASplitIndexBinaryStorage<'a, T, TV, P, PV, M>,
> where
    for<'t> TV::Ref<'t>: IndexSortable,
    for<'t> PV::Ref<'t>: IndexSortable,
{
    config: SplitStorageOptions,
    index: &'a I,
    bands: Vec<SplitBand>,
    _t: PhantomData<T>,
    _tv: PhantomData<TV>,
    _p: PhantomData<P>,
    _pv: PhantomData<PV>,
    _m: PhantomData<M>,
}

impl<
        'a,
        T: ArrowStorage + 'a + IndexSortable + Clone + SplitArrowStorage + SoAArrowStorage<TV>,
        TV: SoAVec<T> + SoAIndexSortable<T> + SoAAppendVec<T>,
        P: ArrowStorage + IndexSortable + SoAArrowStorage<PV>,
        PV: SoAVec<P> + SoAIndexSortable<P> + SoAAppendVec<P>,
        M: ArrowStorage,
        I: SoASplitIndexBinaryStorage<'a, T, TV, P, PV, M>,
    > SoASplitIndexBinaryStorageWriter<'a, T, TV, P, PV, M, I>
where
    for<'t> TV::Ref<'t>: IndexSortable,
    for<'t> PV::Ref<'t>: IndexSortable,
{
    pub fn new(config: SplitStorageOptions, index: &'a I, bands: Vec<SplitBand>) -> Self {
        Self {
            config,
            index,
            bands,
            _t: PhantomData,
            _p: PhantomData,
            _m: PhantomData,
            _tv: PhantomData,
            _pv: PhantomData,
        }
    }

    pub fn write(&mut self, directory: &Path, compression_level: &Compression) -> io::Result<()> {
        match self.config.bin_storage_strategy {
            BinStorageStrategy::SingleFile => self.write_single_file(directory, compression_level),
            BinStorageStrategy::FilePerBin => self.write_bin_per_file(directory, compression_level),
            BinStorageStrategy::NBinsPerFile(_) => {
                self.write_n_bands_per_file(directory, compression_level)
            }
            BinStorageStrategy::NEntriesPerFile(_) => {
                // self.write_n_entries_per_file(directory, compression_level)
                todo!()
            }
        }
    }

    pub fn write_single_file(
        &mut self,
        directory: &Path,
        compression_level: &Compression,
    ) -> io::Result<()> {
        self.write_bands(0..self.bands.len(), 0, directory, compression_level)?;
        Ok(())
    }

    pub fn write_bin_per_file(
        &mut self,
        directory: &Path,
        compression_level: &Compression,
    ) -> io::Result<()> {
        let idxes = 0..self.bands.len();
        for band_i in idxes {
            self.write_bands(band_i..(band_i + 1), band_i, directory, compression_level)?;
        }
        Ok(())
    }

    pub fn write_n_bands_per_file(
        &mut self,
        directory: &Path,
        compression_level: &Compression,
    ) -> io::Result<()> {
        if let BinStorageStrategy::NBinsPerFile(n_per_file) =
            self.config.bin_storage_strategy.clone()
        {
            let mut start = 0usize;
            let n = self.bands.len();
            let mut bands_i = 0;
            while start < n {
                let end = if (start + n_per_file) > n {
                    n
                } else {
                    start + n_per_file
                };
                let indices = start..end;
                self.write_bands(indices, bands_i, directory, compression_level)?;
                bands_i += 1;
                start += n_per_file;
            }
        }
        Ok(())
    }

    pub fn write_bands(
        &mut self,
        band_indices: Range<usize>,
        bands_i: usize,
        directory: &Path,
        compression_level: &Compression,
    ) -> io::Result<()> {
        let archive_name = self
            .config
            .bin_storage_strategy
            .make_file_name::<T>(&self.bands, bands_i);

        let entries_path = directory.join(archive_name.clone());
        let entries_schema = T::schema();
        let props = T::writer_properties()
            .set_compression(compression_level.clone())
            // .set_column_encoding("band_id".into(), parquet::basic::Encoding::RLE)
            .set_writer_version(parquet::file::properties::WriterVersion::PARQUET_2_0)
            .set_statistics_enabled(parquet::file::properties::EnabledStatistics::Page)
            .build();

        log::debug!(
            "Opening {archive_name} for bands {}..{}",
            band_indices.start,
            band_indices.end
        );
        let ext_schema = I::make_item_schema();
        let mut writer = ArrowWriter::try_new(
            fs::File::create(entries_path)?,
            ext_schema.clone(),
            Some(props),
        )?;

        let n_bins = self.index.iter_entries().count();
        let mut bin_counts = vec![0usize; n_bins];

        for j in band_indices {
            let band: &mut SplitBand = self.bands.get_mut(j).unwrap();
            // log::debug!(
            //     "Writing band {j}: {:0.2}-{:0.2}",
            //     band.start_mass,
            //     band.end_mass
            // );
            let interval = Interval::new(band.start_id as usize, band.end_id as usize + 1);
            // log::debug!("Indices: {interval:?}");
            band.file_name = Some(archive_name.clone());

            for (i, bin ) in self.index.iter_entries().enumerate() {
                // let idx = SoAIndexBin::search_parent_id(bin, interval);
                // let entries_of: <TV as SoAVec<T>>::Slice<'_> = SoAIndexBin::slice(bin, idx);
                let entries_of: <TV as SoAVec<T>>::Slice<'_> = SoAIndexBin::select_parent_id(bin, interval);

                let n_entries_of = entries_of.len();
                log::debug!("Band {j} Bin {i} Count {n_entries_of}");
                if n_entries_of == 0 {
                    continue;
                }
                bin_counts[i as usize] += n_entries_of;

                let batch = T::to_batch_soa(entries_of, entries_schema.clone(), i as u64);
                let (_fields, mut arrays, _null_buffer) = StructArray::from(batch).into_parts();
                let band_id_col: Vec<u32> = vec![band.band_id; n_entries_of];
                let band_id_col = Arc::new(UInt32Array::from(band_id_col));
                arrays.push(band_id_col);
                let batch = RecordBatch::try_new(ext_schema.clone(), arrays).unwrap();
                writer.write(&batch)?;
            }
        }
        log::debug!("{} items in bins", bin_counts[11904]);
        writer.close()?;
        Ok(())
    }
}

pub trait SoASplitIndexBinaryStorage<
    'a,
    T: ArrowStorage + 'a + IndexSortable + Clone + SplitArrowStorage + SoAArrowStorage<TV>,
    TV: SoAVec<T> + SoAIndexSortable<T> + SoAAppendVec<T>,
    P: ArrowStorage + IndexSortable + SoAArrowStorage<PV>,
    PV: SoAVec<P> + SoAIndexSortable<P> + SoAAppendVec<P>,
    M: ArrowStorage,
>: SoAIndexBinaryStorage<'a, T, TV, P, PV, M> + Sized where
    for<'t> TV::Ref<'t>: IndexSortable,
    for<'t> PV::Ref<'t>: IndexSortable,
{
    fn make_item_schema() -> SchemaRef {
        T::split_schema()
    }

    fn compute_parent_bands(&self, bin_width: MassType) -> Vec<SplitBand> {
        let parents = self.parents();
        let mut min_value = parents.mass().first().copied().unwrap_or_default();
        let mut max_value = min_value + bin_width;
        let mut min_index = 0usize;
        let mut bands = Vec::new();
        for (i, val) in parents.mass().iter().copied().enumerate() {
            if val > max_value {
                bands.push(SplitBand::new(
                    bands.len() as ParentID,
                    min_index as ParentID,
                    i.saturating_sub(1) as ParentID,
                    min_value,
                    max_value,
                    None,
                ));
                min_index = i;
                min_value = val;
                max_value = min_value + bin_width;
            }
        }
        bands.push(SplitBand::new(
            bands.len() as ParentID,
            min_index as ParentID,
            parents.len() as ParentID,
            min_value,
            max_value,
            None,
        ));
        bands
    }

    fn write_entries_split(
        &'a self,
        directory: &Path,
        mut bands: Vec<SplitBand>,
        compression_level: &Compression,
    ) -> io::Result<Vec<SplitBand>> {
        let opts = SplitStorageOptions::default();
        let mut writer = SoASplitIndexBinaryStorageWriter::new(opts, self, bands);
        writer.write(directory, compression_level)?;
        bands = writer.bands;
        Ok(bands)
    }

    fn write_split<D: AsRef<Path>>(
        &'a self,
        directory: &D,
        bin_options: SplitStorageOptions,
        compression_level: Option<Compression>,
    ) -> io::Result<()> {
        let directory = directory.as_ref();

        let compression_level =
            compression_level.unwrap_or_else(|| Compression::ZSTD(ZstdLevel::try_new(9).unwrap()));
        let mut bands = self.compute_parent_bands(bin_options.bin_width);

        self.write_metadata(directory)?;
        self.write_parents(directory, &compression_level)?;
        bands = self.write_entries_split(directory, bands, &compression_level)?;
        self.write_split_log(directory, &bands)?;
        Ok(())
    }

    fn band_log_name() -> String {
        SplitBand::archive_name()
    }

    fn write_split_log(&self, directory: &Path, split_log: &[SplitBand]) -> io::Result<()> {
        let split_log_path = directory.join(Self::band_log_name());

        let split_log_fh = io::BufWriter::new(fs::File::create(split_log_path)?);
        let mut writer = LineDelimitedWriter::new(split_log_fh);
        let split_log = SplitBand::to_batch(split_log, SplitBand::schema(), 0).unwrap();

        writer.write(&split_log).unwrap();
        writer.finish().unwrap();
        Ok(())
    }

    fn read_split<D: AsRef<Path>>(directory: &D) -> io::Result<Self>
    where
        Self: Sized,
    {
        let root = directory.as_ref();

        let metadata = Self::read_metadata(root)?;
        let parents = Self::read_parents(root)?;
        let split_log = Self::read_split_log(root)?;

        let file_names: HashSet<&str> = split_log
            .iter()
            .flat_map(|s| s.file_name.as_deref())
            .collect();
        let mut file_names: Vec<_> = file_names.into_iter().collect();
        file_names.sort();

        let entries = {
            let mut bin_collector: HashMap<u64, TV> = HashMap::default();
            for archive_name in file_names.iter() {
                let entries_fh = fs::File::open(root.join(archive_name))?;
                let reader = ArrowReaderBuilder::try_new(entries_fh)?.build()?;
                let entry_schema = T::schema();

                for batch in reader {
                    let (batch, segments) =
                        T::from_batch_soa(&batch.unwrap(), entry_schema.clone());
                    let parts = T::partition_by_segments(batch, segments);
                    for (seg, mut part) in parts.into_iter() {
                        match bin_collector.entry(seg) {
                            std::collections::hash_map::Entry::Occupied(mut occupied_entry) => {
                                occupied_entry.get_mut().append(&mut part);
                            }
                            std::collections::hash_map::Entry::Vacant(vacant_entry) => {
                                vacant_entry.insert(part);
                            }
                        }
                    }
                }
            }

            bin_collector
        };

        let this = Self::from_components(metadata, parents, entries);
        Ok(this)
    }

    fn read_split_log(directory: &Path) -> io::Result<Vec<SplitBand>> {
        let split_log_path = directory.join(Self::band_log_name());
        let split_log_fh = io::BufReader::new(fs::File::open(split_log_path)?);

        let split_log = ReaderBuilder::new(SplitBand::schema())
            .build(split_log_fh)
            .unwrap()
            .map(|batch| {
                let batch: Vec<_> = SplitBand::from_batch(&batch.unwrap(), SplitBand::schema())
                    .map(|(x, _)| x)
                    .collect();
                batch
            })
            .flatten()
            .collect();
        Ok(split_log)
    }

    fn read_parents(directory: &Path) -> io::Result<PV> {
        let parents_path = directory.join(P::archive_name());
        let parent_schema = P::schema();
        let parents_fh = fs::File::open(parents_path)?;

        let reader = ArrowReaderBuilder::try_new(parents_fh)?.build()?;
        let mut parents = PV::new();
        for batch in reader {
            parents.append(&mut P::from_batch_soa(&batch.unwrap(), parent_schema.clone()).0);
        }
        Ok(parents)
    }

    fn read_metadata(directory: &Path) -> io::Result<M> {
        let meta_path = directory.join(M::archive_name());

        let metadata = {
            let meta_schema = M::schema();
            let meta_fh = io::BufReader::new(fs::File::open(meta_path)?);
            let meta_rec = ReaderBuilder::new(meta_schema.clone())
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
        Ok(metadata)
    }
}
