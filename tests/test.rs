#[cfg(test)]
use std::{fs, io, mem};

use csv;
use rand::Rng;

// #[cfg(feature = "parallelism")]
use rayon::prelude::*;

use mass_fragment_index::fragment::{Fragment, FragmentName};
use mass_fragment_index::index::SearchIndex;
use mass_fragment_index::parent::Peptide;
use mass_fragment_index::sort::{IndexSortable, MassType, ParentID, SortType, Tolerance};

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

fn test_exact(
    search_index: &SearchIndex<Fragment, Peptide>,
    entries: &Vec<(Peptide, Vec<Fragment>)>,
    precursor_error_tolerance: MassType,
    product_error_tolerance: MassType,
) -> bool {
    let precursor_error_tolerance = Tolerance::PPM(precursor_error_tolerance);
    let product_error_tolerance = Tolerance::PPM(product_error_tolerance);

    let it = {
        // #[cfg(feature = "parallelism")]
        {
            entries.par_iter()
        }
        // #[cfg(not(feature = "parallelism"))]
        // {
        //     entries.iter()
        // }
    };

    it.for_each(|(pept, frags)| {
        let parent_interval = search_index.parents_for(pept.mass, precursor_error_tolerance);
        for i in parent_interval {
            assert!(parent_interval.contains(i));
            let parent = &search_index.parents[i];
            assert!(
                precursor_error_tolerance.test(parent.mass, pept.mass),
                "Parent molecule {}/{:?} {:?} did not pass precursor error tolerance {} for {:?}",
                i,
                parent_interval,
                parent,
                ((parent.mass - pept.mass) / pept.mass).abs() * 1e6,
                pept,
            );
        }

        let n_expected_parents = entries
            .iter()
            .filter(|(alt_pept, _)| precursor_error_tolerance.test(alt_pept.mass, pept.mass))
            .count();

        let expected_parents: Vec<_> = entries
            .iter()
            .map(|(alt_pept, _)| alt_pept)
            .filter(|alt_pept| precursor_error_tolerance.test(alt_pept.mass, pept.mass))
            .collect();

        let n_parents = parent_interval.clone().into_iter().count();

        assert_eq!(
            n_expected_parents, n_parents,
            "Expected {} parent matches, found {} for query {}: {:?}/{:?}/{:?}",
            n_expected_parents, n_parents, pept.mass, parent_interval, pept, expected_parents
        );

        for expected_frag in frags.iter() {
            let search_iter = search_index.search(
                expected_frag.mass,
                product_error_tolerance,
                Some(parent_interval),
            );

            let n_frags = search_iter
                .enumerate()
                .map(|(i, frag)| {
                    assert!(
                        parent_interval.contains(frag.parent_id as usize),
                        "{}th hit {:?} did not fall within parent interval {:?} {:?} of batch {}",
                        i,
                        frag,
                        parent_interval,
                        expected_frag,
                        pept.id,
                    );
                    assert!(
                        product_error_tolerance.test(frag.mass, expected_frag.mass),
                        "{}th {:?} did not fall within mass error tolerance for {:?} {}",
                        i,
                        frag,
                        expected_frag,
                        (frag.mass - expected_frag.mass).abs() / frag.mass * 1e6,
                    );
                    frag
                })
                .count();

            let n_expected_frags = entries
                .iter()
                .filter(|(alt_pept, _)| precursor_error_tolerance.test(alt_pept.mass, pept.mass))
                .map(|(_, frags)| {
                    frags
                        .iter()
                        .filter(|frag| product_error_tolerance.test(frag.mass, expected_frag.mass))
                })
                .flatten()
                .count();
            assert_eq!(n_expected_frags, n_frags);
        }
    });
    true
}

fn test_permuted(
    search_index: &SearchIndex<Fragment, Peptide>,
    entries: &Vec<(Peptide, Vec<Fragment>)>,
    precursor_error_tolerance: MassType,
    product_error_tolerance: MassType,
) -> bool {
    let precursor_error_tolerance = Tolerance::PPM(precursor_error_tolerance);
    let product_error_tolerance = Tolerance::PPM(product_error_tolerance);
    let it = {
        // #[cfg(feature = "parallelism")]
        {
            entries.par_iter()
        }
        // #[cfg(not(feature = "parallelism"))]
        // {
        //     entries.iter()
        // }
    };
    it.for_each(|(pept, frags)| {
        let mut rng = rand::thread_rng();
        let mass_error_scale = ((precursor_error_tolerance) * (1.0)).bounds(pept.mass());

        // When the random offset is very close to the outer limit of the tolerance window, things
        // seem to fall apart. Try using a larger dataset?
        let precursor = rng.gen_range(mass_error_scale.0..=mass_error_scale.1);

        let parent_interval =
            search_index.parents_for(precursor, precursor_error_tolerance);
        for i in parent_interval {
            assert!(parent_interval.contains(i));
            let parent = &search_index.parents[i];
            if parent == pept {
                assert!(
                    precursor_error_tolerance.test(parent.mass, precursor)
                )
            }
            assert!(
                precursor_error_tolerance.test(parent.mass, precursor),
                "Parent molecule {}/{:?} {:?} did not pass precursor error tolerance {} for {:?}  with random offset {} ({}% of {})",
                i,
                parent_interval,
                parent,
                ((parent.mass - (precursor)) / (precursor)).abs(),
                pept,
                precursor - pept.mass,
                (precursor - pept.mass).abs() * 100.0 / mass_error_scale.1,
                mass_error_scale.1
            );
        }

        let n_expected_parents = entries
            .iter()
            .filter(|(alt_pept, _)| {
                precursor_error_tolerance.test(alt_pept.mass, precursor)
            })
            .count();

        let n_parents = parent_interval.clone().into_iter().count();

        assert_eq!(
            n_expected_parents,
            n_parents,
            "Expected parent count differed for {} with random offset {} ({}% of {})",
            pept.mass,
            precursor - pept.mass,
            (precursor - pept.mass).abs() * 100.0 / mass_error_scale.1,
            mass_error_scale.1
        );

        for expected_frag in frags.iter() {
            let product_mass_error_scale = ((product_error_tolerance) * (1.0 - 1e-6)).bounds(expected_frag.mass);
            let product = rng.gen_range(product_mass_error_scale.0..=product_mass_error_scale.1);

            let search_iter = search_index.search(
                product,
                product_error_tolerance,
                Some(parent_interval),
            );

            let n_frags =
                search_iter
                    .enumerate()
                    .map(|(i, frag)| {
                        assert!(
                            parent_interval.contains(frag.parent_id as usize),
                            "{}th hit {:?} did not fall within parent interval {:?} {:?} of batch {}",
                            i, frag, parent_interval, expected_frag, pept.sequence
                        );
                        assert!(
                            product_error_tolerance.test(product, frag.mass()),
                            "{}th hit {} PPM exceeds error tolerance for {:?}",
                            i,
                            (frag.mass() - product).abs() * 1e6 / product,
                            product,
                        );
                        frag
                    })
                    .count();

            let n_expected_frags = entries
                .iter()
                .filter(|(alt_pept, _)| {
                    precursor_error_tolerance.test(alt_pept.mass(), precursor)
                })
                .map(|(_, frags)| {
                    frags.iter().filter(|frag| {
                        product_error_tolerance.test(product, frag.mass())
                    })
                })
                .flatten()
                .count();
            assert!(n_expected_frags == n_frags);
        }
    });
    true
}

#[test]
fn test_index_build_traversal() -> io::Result<()> {
    let reader = io::BufReader::new(fs::File::open("tests/data/test_data.csv")?);
    let search_index = build_index(reader)?;
    let reader = io::BufReader::new(fs::File::open("tests/data/test_data.csv")?);
    let entries = parse_csv(reader)?;

    if MassType::DIGITS == f32::DIGITS {
        assert_eq!(1013803, search_index.bins.len());
    } else if MassType::DIGITS == f64::DIGITS {
        assert_eq!(1000001, search_index.bins.len());
    }
    assert_eq!(search_index.parents.len(), entries.len());
    assert_eq!(
        search_index.bins.iter().map(|b| b.len()).sum::<usize>(),
        entries.iter().map(|(_, frags)| frags.len()).sum::<usize>()
    );

    let precursor_error_tolerance = 5f64 as MassType;
    let product_error_tolerance = 10f64 as MassType;

    assert!(test_exact(
        &search_index,
        &entries,
        precursor_error_tolerance,
        product_error_tolerance
    ));

    (0..5).into_iter().for_each(|_| {
        assert!(test_permuted(
            &search_index,
            &entries,
            precursor_error_tolerance,
            product_error_tolerance
        ));
    });

    Ok(())
}
