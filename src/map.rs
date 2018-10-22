use bincode::deserialize_from;
use bio::data_structures::bwt::Occ;
use bio::data_structures::fmindex::{FMDIndex, FMIndex, FMIndexable, Interval};
use bio::io::fastq;
use libflate::deflate::Decoder;
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;

pub fn run(reads_path: &str) -> Result<(), Box<Error>> {
    let intervals = calculate_intervals(reads_path)?;

    debug!("Load suffix array");
    let f_suffix_array = File::open("reference.sar")?;
    let d_suffix_array = Decoder::new(f_suffix_array);
    let suffix_array: Vec<usize> = deserialize_from(d_suffix_array)?;

    debug!("Print results");
    // Loop through the results, extracting the positions array for each pattern
    for interval_calculator in intervals {
        let positions = interval_calculator.occ(&suffix_array);
        println!("{:#?}", positions);
    }
    Ok(())
}

fn calculate_intervals(reads_path: &str) -> Result<Vec<Interval>, Box<Error>> {
    debug!("Load BWT");
    let f_bwt = File::open("reference.bwt")?;
    let d_bwt = Decoder::new(f_bwt);
    let bwt: Vec<u8> = deserialize_from(d_bwt)?;
    debug!("Load \"C\" table");
    let f_less = File::open("reference.less")?;
    let d_less = Decoder::new(f_less);
    let less: Vec<usize> = deserialize_from(d_less)?;
    debug!("Load \"Occ\" table");
    let f_occ = File::open("reference.occ")?;
    let d_occ = Decoder::new(f_occ);
    let occ: Occ = deserialize_from(d_occ)?;

    debug!("Reconstruct FMD-index");
    let fm_index = FMIndex::new(&bwt, &less, &occ);
    let fmd_index = FMDIndex::from(fm_index);

    debug!("Map reads");
    // TODO: Load reads in batches to memory to be able to process them in parallel
    let reads_fq_reader = fastq::Reader::from_file(reads_path)?;
    let intervals = reads_fq_reader
        .records()
        .flat_map(|pattern| {
            let pattern = pattern.unwrap().seq().to_ascii_uppercase();
            k_mismatch_search(
                &pattern,
                // These values are borrowed from BWA
                match pattern.len() {
                    0..=14 => 0,
                    15..=37 => 2,
                    38..=63 => 3,
                    64..=92 => 4,
                    93..=123 => 5,
                    124..=156 => 6,
                    _ => 7,
                },
                &fmd_index,
            )
        }).collect::<Vec<_>>();
    Ok(intervals)
}

/// Finds all suffix array intervals for the current pattern with up to z mismatches
pub fn k_mismatch_search(
    pattern: &[u8],
    z: u32,
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
) -> HashSet<Interval> {
    let d = calculate_d(&pattern, fmd_index);
    let interval = Interval {
        lower: 0,
        upper: fmd_index.bwt().len() - 1,
    };
    let mut interval_set = k_mismatch_search_recursive(
        pattern,
        &d,
        fmd_index,
        pattern.len() as isize - 1,
        z as i32,
        interval,
    );

    interval_set = interval_set
        .iter()
        .map(|i| Interval {
            lower: i.lower,
            upper: i.upper + 1,
        }).collect();

    interval_set
}

/// Follows closely the implementation of BWA-backtrack (Li & Durbin, 2009)
fn k_mismatch_search_recursive(
    pattern: &[u8],
    d: &[i32],
    fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>,
    j: isize,
    z: i32,
    interval: Interval,
) -> HashSet<Interval> {
    // Too many mismatches
    if z < d[if j < 0 { 0 } else { j } as usize] {
        return HashSet::new();
    }

    // This route through the read graph is finished successfully, return the interval
    if j < 0 {
        return [interval].iter().cloned().collect();
    }

    let mut interval_set = HashSet::new();

    // Insertion
    interval_set = interval_set
        .union(&k_mismatch_search_recursive(
            pattern,
            d,
            fmd_index,
            j - 1,
            z - 1,
            interval,
        )).cloned()
        .collect();

    for &c in b"ACGT".iter() {
        let less = fmd_index.less(c);

        let interval_prime = Interval {
            lower: less + if interval.lower > 0 {
                fmd_index.occ(interval.lower - 1, c)
            } else {
                0
            },
            upper: less + fmd_index.occ(interval.upper, c) - 1,
        };

        if interval.lower <= interval.upper {
            // Deletion
            interval_set = interval_set
                .union(&k_mismatch_search_recursive(
                    pattern,
                    d,
                    fmd_index,
                    j,
                    z - 1,
                    interval_prime,
                )).cloned()
                .collect();

            // Match
            if c == pattern[j as usize] {
                interval_set = interval_set
                    .union(&k_mismatch_search_recursive(
                        pattern,
                        d,
                        fmd_index,
                        j - 1,
                        z,
                        interval_prime,
                    )).cloned()
                    .collect();

            // Mismatch
            } else {
                interval_set = interval_set
                    .union(&k_mismatch_search_recursive(
                        pattern,
                        d,
                        fmd_index,
                        j - 1,
                        z - 1,
                        interval_prime,
                    )).cloned()
                    .collect();
            }
        }
    }
    interval_set
}

/// Calculate lower bound of mismatches for pruning of the search tree
fn calculate_d(pattern: &[u8], fmd_index: &FMDIndex<&Vec<u8>, &Vec<usize>, &Occ>) -> Vec<i32> {
    let mut d = vec![0; pattern.len()];

    let mut z = 0;
    let mut j = 0;
    for i in 0..pattern.len() {
        // TODO: It's inefficient to search for the whole part of the pattern again and again
        let interval = fmd_index.backward_search(pattern[j..=i].iter());
        if interval.lower >= interval.upper {
            z += 1;
            j = i + 1;
        }
        d[i] = z;
    }
    d
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::alphabets;
    use bio::data_structures::bwt::{bwt, less};
    use bio::data_structures::fmindex::{FMDIndex, FMIndex};
    use bio::data_structures::suffix_array::suffix_array;

    #[test]
    fn test_exact_search() {
        // Prepare reference for usage with FM(D)-index by
        // making it double stranded and adding sentinels
        let mut ref_seq = "GATTACA".as_bytes().to_owned();
        let ref_seq_rev_compl = alphabets::dna::revcomp(ref_seq.iter());
        ref_seq.extend_from_slice(b"$");
        ref_seq.extend_from_slice(&ref_seq_rev_compl);
        drop(ref_seq_rev_compl);
        ref_seq.extend_from_slice(b"$");

        let alphabet = alphabets::dna::n_alphabet();

        let sa = suffix_array(&ref_seq);
        let bwt = bwt(&ref_seq, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);

        let fm_index = FMIndex::new(&bwt, &less, &occ);
        let fmd_index = FMDIndex::from(fm_index);

        // Test occurring pattern
        let pattern_1 = b"TA";
        let sai_1 = fmd_index.backward_search(pattern_1.iter());
        let positions_1 = sai_1.occ(&sa);
        assert_eq!(positions_1, [10, 3]);

        // Test non-occurring pattern
        let pattern_2 = b"GG";
        let sai_2 = fmd_index.backward_search(pattern_2.iter());
        let positions_2 = sai_2.occ(&sa);
        assert_eq!(positions_2, []);
    }

    #[test]
    fn test_inexact_search() {
        let mut ref_seq = "GATTACA".as_bytes().to_owned();
        let ref_seq_rev_compl = alphabets::dna::revcomp(ref_seq.iter());
        ref_seq.extend_from_slice(b"$");
        ref_seq.extend_from_slice(&ref_seq_rev_compl);
        drop(ref_seq_rev_compl);
        ref_seq.extend_from_slice(b"$");

        let alphabet = alphabets::dna::n_alphabet();

        let sa = suffix_array(&ref_seq);
        let bwt = bwt(&ref_seq, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);

        let fm_index = FMIndex::new(&bwt, &less, &occ);
        let fmd_index = FMDIndex::from(fm_index);

        let pattern = "GTTT".as_bytes().to_owned();

        let d = calculate_d(&pattern, &fmd_index);
        assert_eq!(d, vec![0, 0, 1, 1]);

        let intervals = k_mismatch_search(&pattern, 1, &fmd_index);
        let positions: Vec<usize> = intervals
            .into_iter()
            .map(|f| f.occ(&sa))
            .flatten()
            .collect();
        assert_eq!(positions, vec![0]);
    }

}
