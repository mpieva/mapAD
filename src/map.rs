use bincode::deserialize_from;
use bio::alphabets;
use bio::data_structures::bwt::Occ;
use bio::data_structures::fmindex::{FMDIndex, FMIndex, FMIndexable, Interval};
use bio::io::fastq;
use libflate::deflate::Decoder;
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
    let interval_calculators = reads_fq_reader
        .records()
        .map(|pattern| {
            fmd_index.backward_search(pattern.unwrap().seq().to_ascii_uppercase().iter())
        }).collect::<Vec<_>>();
    Ok(interval_calculators)
}
