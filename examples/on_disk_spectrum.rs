use std::{env, io, path::PathBuf};

use mass_fragment_index::{
    storage::{IndexMetadata, SearchIndexOnDisk},
    DeconvolutedPeak, Spectrum,
};

fn main() -> io::Result<()> {
    pretty_env_logger::init_timed();

    let mut args = env::args().skip(1);

    let storage_dir: PathBuf = args
        .next()
        .unwrap_or_else(|| panic!("Please provide a path to read index from"))
        .into();

    let index: SearchIndexOnDisk<DeconvolutedPeak, Spectrum, IndexMetadata> =
        SearchIndexOnDisk::new(storage_dir)?;

    let iv = index.parents_for_range(1200.0, 2300.0, mass_fragment_index::Tolerance::Da(0.2));
    eprintln!("{iv:?}");
    let peaks = index.load_for_parents(iv)?;

    let query = 1200.0;
    let tol = mass_fragment_index::Tolerance::Da(0.2);

    let mut query_hit = 0;
    for peak in peaks.iter() {
        assert!(iv.contains(peak.scan_ref as usize));
        query_hit += tol.test(peak.mass, query) as i32;
        if tol.test(peak.mass, query) {
            eprintln!("{:?}", peak);
        }
    }
    eprintln!("{query_hit}");
    Ok(())
}
