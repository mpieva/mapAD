extern crate bio;
extern crate clap;
#[macro_use]
extern crate log;
extern crate simple_logger;

use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMDIndex, FMIndex, FMIndexable};
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

        debug!("Add reverse complement and sentinels to reference");
        let ref_seq_rev_compl = alphabets::dna::revcomp(&ref_seq);
        ref_seq.extend_from_slice(b"$");
        ref_seq.extend_from_slice(&ref_seq_rev_compl);
        drop(ref_seq_rev_compl);
        ref_seq.extend_from_slice(b"$");

        // TODO: Load reads in batches to memory to be able to process them in parallel
        let reads_fq_reader = fastq::Reader::from_file(matches.value_of("reads").unwrap()).unwrap();

        debug!("Rank-transform input sequence");
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
        let suffix_array = suffix_array(&ref_seq);

        debug!("Generate BWT");
        let bwt = bwt(&ref_seq, &suffix_array);

        debug!("Drop source sequence");
        drop(ref_seq);

        debug!("Generate \"C\" table");
        let less = less(&bwt, &rank_alphabet);

        debug!("Generate \"Occ\" table");
        let occ = Occ::new(&bwt, 128, &rank_alphabet);

        debug!("Generate FM-index");
        let fm_index = FMIndex::new(
            &bwt,
            &less,
            &occ,
            alphabets::RankTransform::new(&bwt_alphabet),
        );

        debug!("Generate FMD-index");
        let fmd_index = FMDIndex::from(fm_index);

        debug!("Map reads");
        let interval_calculators = reads_fq_reader
            .records()
            .map(|pattern| {
                let pattern = pattern.unwrap().seq().to_ascii_uppercase();
                fmd_index.backward_search(rank_transform.transform(&pattern).iter())
            }).collect::<Vec<_>>();

        debug!("Print results");
        // Loop through the results, extracting the positions array for each pattern
        for interval_calculator in interval_calculators {
            let positions = interval_calculator.occ(&suffix_array);
            println!("{:#?}", positions);
        }

        // There is actually only one sequence in the reference genome fasta file for now
        break;
    }
}
