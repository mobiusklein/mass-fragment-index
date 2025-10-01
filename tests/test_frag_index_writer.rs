#![cfg(test)]
#![allow(unused)]

use std::{fs, io, mem};

use csv;

use mass_fragment_index::fragment::{Fragment, FragmentName};
use mass_fragment_index::index::SearchIndex;
use mass_fragment_index::parent::Peptide;
use mass_fragment_index::sort::{MassType, ParentID, SortType};

use mass_fragment_index::storage::{read_fragment_index, write_fragment_index};
use mass_fragment_index::Tolerance;

fn parse_csv<R: io::BufRead>(reader: R) -> io::Result<Vec<(Peptide, Vec<Fragment>)>> {
    let mut csv_reader = csv::Reader::from_reader(reader);
    let mut parent_i: i32 = -1;

    let mut accumulator = Vec::new();
    let mut peprec: Option<Peptide> = None;
    let mut fragments = Vec::new();

    for result in csv_reader.records() {
        let record = result?;
        match record.get(0).expect("Missing record type") {
            "PEPTIDE" => {
                parent_i += 1;
                if let Some(ref peptide) = peprec {
                    accumulator.push((peptide.clone(), mem::take(&mut fragments)));
                }

                peprec = Some(Peptide::new(
                    record.get(2).unwrap().parse::<MassType>().unwrap(),
                    parent_i as ParentID,
                    0,
                    0,
                    record.get(1).unwrap().to_string(),
                ));
            }
            "FRAGMENT" => {
                let name: FragmentName = record.get(1).unwrap().parse().unwrap();
                let frag = Fragment::new(
                    record.get(2).unwrap().parse::<MassType>().unwrap(),
                    parent_i as ParentID,
                    name.0,
                    name.1,
                );
                fragments.push(frag);
            }
            field => {
                panic!("Unknown record type {}", field)
            }
        }
    }
    if let Some(peptide) = peprec {
        accumulator.push((peptide, mem::take(&mut fragments)));
    }
    Ok(accumulator)
}

fn build_index<R: io::BufRead>(reader: R) -> io::Result<SearchIndex<Fragment, Peptide>> {
    let pepfrags = parse_csv(reader)?;
    let mut search_index: SearchIndex<Fragment, Peptide> = SearchIndex::empty(100, 10000.0);
    pepfrags.into_iter().for_each(|(pep, frags)| {
        search_index.add_parent(pep);
        frags.into_iter().for_each(|frag| {
            search_index.add(frag);
        });
    });

    search_index.sort(SortType::ByParentId);
    Ok(search_index)
}

#[test]
fn test_index_build_traversal() -> io::Result<()> {
    let reader = io::BufReader::new(fs::File::open("tests/data/test_data.csv")?);
    let search_index: SearchIndex<Fragment, Peptide> = build_index(reader)?;

    let tmpdir = tempfile::tempdir()?;
    let tmpdir_path = tmpdir.path();
    search_index.write_parquet(&tmpdir_path, None)?;

    let duplicate_index = SearchIndex::<Fragment, Peptide>::read_parquet(&tmpdir_path)?;
    assert_eq!(duplicate_index.parents.len(), search_index.parents.len());
    assert_eq!(duplicate_index.num_entries(), search_index.num_entries());

    let parent_interval = search_index.parents_for_range(200.0, 1200.0, Tolerance::PPM(10.0));
    let search: Vec<_> = search_index
        .search(
            113.08406397713001,
            Tolerance::PPM(10.0),
            Some(parent_interval),
        )
        .collect();
    assert_eq!(search.len(), 504);

    let parent_interval = duplicate_index.parents_for_range(200.0, 1200.0, Tolerance::PPM(10.0));
    let dup_search: Vec<_> = duplicate_index
        .search(
            113.08406397713001,
            Tolerance::PPM(10.0),
            Some(parent_interval),
        )
        .collect();
    assert_eq!(dup_search.len(), 504);

    Ok(())
}
