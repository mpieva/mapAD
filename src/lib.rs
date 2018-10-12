extern crate bincode;
extern crate bio;
#[macro_use]
extern crate log;

use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMDIndex, FMIndex, FMIndexable, Interval};
use bio::data_structures::suffix_array::suffix_array;
use bio::io::{fasta, fastq};
use std::error::Error;
use std::fs::File;

pub fn index(
    reference_path: &str,
    bwt_alphabet: &alphabets::Alphabet,
    rank_alphabet: &alphabets::Alphabet,
) -> Result<(), Box<Error>> {
    // Handle on reference FASTA
    // TODO: Would this iterator only return one result?
    for record in fasta::Reader::from_file(reference_path).unwrap().records() {
        debug!("Convert reference to uppercase letters");
        let mut ref_seq = record.unwrap().seq().to_ascii_uppercase();

        debug!("Add reverse complement and sentinels to reference");
        let ref_seq_rev_compl = alphabets::dna::revcomp(&ref_seq);
        ref_seq.extend_from_slice(b"$");
        ref_seq.extend_from_slice(&ref_seq_rev_compl);
        drop(ref_seq_rev_compl);
        ref_seq.extend_from_slice(b"$");

        let rank_transform = alphabets::RankTransform::new(&bwt_alphabet);
        let ref_seq = rank_transform.transform(&ref_seq);

        debug!("Generate suffix array");
        let suffix_array = suffix_array(&ref_seq);

        debug!("Generate BWT");
        let bwt = bwt(&ref_seq, &suffix_array);

        debug!("Save suffix array to disk");
        let mut f_suffix_array = File::create("reference.sar")?;
        bincode::serialize_into(f_suffix_array, &suffix_array)?;

        debug!("Drop suffix array");
        drop(suffix_array);

        debug!("Drop source sequence");
        drop(ref_seq);

        debug!("Generate \"C\" table");
        let less = less(&bwt, &rank_alphabet);

        debug!("Generate \"Occ\" table");
        let occ = Occ::new(&bwt, 128, &rank_alphabet);

        debug!("Save FMD-index to disk");
        let mut f_less = File::create("reference.less")?;
        bincode::serialize_into(f_less, &less)?;
        let mut f_bwt = File::create("reference.bwt")?;
        bincode::serialize_into(f_bwt, &bwt)?;
        let mut f_occ = File::create("reference.occ")?;
        bincode::serialize_into(f_occ, &occ)?;

        // There is actually only one sequence in the reference genome fasta file for now
        break;
    }
    Ok(())
}

pub fn map(reads_path: &str, bwt_alphabet: &alphabets::Alphabet) -> Result<(), Box<Error>> {
    let intervals = calculate_intervals(reads_path, bwt_alphabet)?;

    debug!("Load suffix array");
    let f_suffix_array = File::open("reference.sar")?;
    let suffix_array: Vec<usize> = bincode::deserialize_from(f_suffix_array)?;

    debug!("Print results");
    // Loop through the results, extracting the positions array for each pattern
    for interval_calculator in intervals {
        let positions = interval_calculator.occ(&suffix_array);
        println!("{:#?}", positions);
    }
    Ok(())
}

fn calculate_intervals(
    reads_path: &str,
    bwt_alphabet: &alphabets::Alphabet,
) -> Result<Vec<Interval>, Box<Error>> {
    debug!("Load FMD-index files");
    let f_less = File::open("reference.less")?;
    let less: Vec<usize> = bincode::deserialize_from(f_less)?;
    let f_bwt = File::open("reference.bwt")?;
    let bwt: Vec<u8> = bincode::deserialize_from(f_bwt)?;
    let f_occ = File::open("reference.occ")?;
    let occ: Occ = bincode::deserialize_from(f_occ)?;

    debug!("Reconstruct FMD-index");
    let rank_transform = alphabets::RankTransform::new(&bwt_alphabet);
    let fm_index = FMIndex::new(
        &bwt,
        &less,
        &occ,
        alphabets::RankTransform::new(&bwt_alphabet),
    );
    let fmd_index = FMDIndex::from(fm_index);

    debug!("Map reads");
    // TODO: Load reads in batches to memory to be able to process them in parallel
    let reads_fq_reader = fastq::Reader::from_file(reads_path)?;
    let interval_calculators = reads_fq_reader
        .records()
        .map(|pattern| {
            fmd_index.backward_search(
                rank_transform
                    .transform(&pattern.unwrap().seq().to_ascii_uppercase())
                    .iter(),
            )
        }).collect::<Vec<_>>();
    Ok(interval_calculators)
}
