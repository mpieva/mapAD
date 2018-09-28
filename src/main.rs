extern crate bio;
extern crate clap;
#[macro_use]
extern crate log;
extern crate simple_logger;

use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;
use bio::io::fasta;
use bio::io::fastq;
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
                .short("r")
                .long("reference")
                .help("FASTA file containing the genome we are about to map against")
                .value_name("FILE"),
        ).get_matches();

    // Handle on reference FASTA
    // TODO: Would this iterator only return one result?
    let mut ref_seq = vec![];
    for record in fasta::Reader::from_file(matches.value_of("reference").unwrap())
        .unwrap()
        .records()
    {
        let record = record.unwrap();
        ref_seq.extend(record.seq().iter().cloned());

        debug!("Convert reference to uppercase");

        debug!("Add sentinel to reference");
        ref_seq.extend_from_slice(b"$");

        // Handle on HT-sequencing reads in FASTQ format
        // TODO: Load reads in batches to memory to be able to process them in parallel
        let reads_fq_reader =
            fastq::Reader::from_file("example/simulated_reads/test.bwa.read1.fastq").unwrap();

        debug!("Rank-transform input sequence");
        // Create an FM-Index for the reference genome
        // TODO: Use FMD-index instead, to not have to search two indices
        // TODO: Read about FMD index
        // TODO: It's perhaps worth to do a rank transform (memory reduction by 10* (?))
        let symbols = b"$ACGTN";
        let rank_symbols = symbols
            .iter()
            .enumerate()
            .map(|(i, _v)| i as u8)
            .collect::<Vec<_>>();

        // Define two alphabets, one with actual chars and
        // another one containing its corresponding ranks
        let bwt_alphabet = alphabets::Alphabet::new(symbols.iter());
        let rank_alphabet = alphabets::Alphabet::new(&rank_symbols);

        let rank_transform = alphabets::RankTransform::new(&bwt_alphabet);
        let ref_seq = rank_transform.transform(&ref_seq);

        debug!("Generate suffix array");
        let sa = suffix_array(&ref_seq);

        debug!("Generate BWT");
        let bwt = bwt(&ref_seq, &sa);

        debug!("Drop source sequence");
        drop(ref_seq);

        debug!("Generate less");
        let less = less(&bwt, &rank_alphabet);

        debug!("Generate OCC");
        let occ = Occ::new(&bwt, 3, &rank_alphabet);

        debug!("Generate FM index");
        let fmindex = FMIndex::new(&bwt, &less, &occ);

        debug!("Map reads");
        let interval_calculators = reads_fq_reader
            .records()
            .map(|pattern| fmindex.backward_search(pattern.unwrap().seq().iter()))
            .collect::<Vec<_>>();

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
