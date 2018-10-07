extern crate bio;
extern crate clap;
#[macro_use]
extern crate log;
extern crate simple_logger;

use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;
use bio::io::{fasta, fastq};
use clap::{App, Arg};

fn main() {
    simple_logger::init().unwrap();

    let matches = App::new("Thrust")
        .about("An ancient aware short-read mapper")
        .arg(
            Arg::with_name("v")
                .short("v")
                .multiple(true)
                .help("Sets the level of verbosity"),
        ).arg(
            Arg::with_name("reference")
                .required(true)
                .long("reference")
                .help("FASTA file containing the genome we are about to map against")
                .value_name("FASTA file"),
        ).arg(
            Arg::with_name("reads")
                .required(true)
                .long("reads")
                .help("FASTQ file containing adapter-trimmed and quality-controlled reads")
                .value_name("FASTQ file"),
        ).get_matches();

    // Handle on reference FASTA
    // TODO: Would this iterator only return one result?
    for record in fasta::Reader::from_file(matches.value_of("reference").unwrap())
        .unwrap()
        .records()
    {
        debug!("Convert reference to uppercase letters");
        let mut ref_seq = record.unwrap().seq().to_ascii_uppercase();

        debug!("Add sentinel character to reference");
        ref_seq.extend_from_slice(b"$");

        // Handle on HT-sequencing reads in FASTQ format
        // TODO: Load reads in batches to memory
        // in order to be able to process them in parallel
        let reads_fq_reader = fastq::Reader::from_file(matches.value_of("reads").unwrap()).unwrap();

        let alphabet = alphabets::dna::n_alphabet();

        debug!("Generate suffix array");
        let sa = suffix_array(&ref_seq);

        debug!("Generate BWT");
        let bwt = bwt(&ref_seq, &sa);

        debug!("Drop source sequence");
        drop(ref_seq);

        debug!("Generate \"C\" table");
        let less = less(&bwt, &alphabet);

        debug!("Generate \"Occ\" table");
        let occ = Occ::new(&bwt, 3, &alphabet);

        debug!("Generate FM index");
        let fm_index = FMIndex::new(&bwt, &less, &occ);

        debug!("Generate FMD index");
        let fmd_index = FMDIndex::from(fm_index);

        debug!("Map reads");
        let interval_calculators = reads_fq_reader
            .records()
            .map(|pattern| {
                fmd_index.backward_search(pattern.unwrap().seq().to_ascii_uppercase().iter())
            }).collect::<Vec<_>>();

        debug!("Print results");
        // Loop through the results, extracting the positions array for each pattern
        for interval_calculator in interval_calculators {
            let positions = interval_calculator.occ(&sa);
            println!("{:#?}", positions);
        }

        // There is actually only one sequence in the reference genome fasta file for now
        break;
    }
}
