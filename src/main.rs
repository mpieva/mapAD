extern crate bio;

use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;
use bio::io::fasta;
use bio::io::fastq;

fn main() {
    // Handle on reference FASTA
    // TODO: Would this iterator only return one result?
    let mut ref_seq = vec![];
    for record in fasta::Reader::from_file("example/chr22.fa")
        .unwrap()
        .records()
    {
        let record = record.unwrap();
        ref_seq.extend(record.seq().iter().cloned());
        println!("Add sentinel to reference");
        ref_seq.extend_from_slice(b"$");

        // Handle on HT-sequencing reads in FASTQ format
        // TODO: Load reads in batches to memory to be able to process them in parallel
        let reads_fq_reader =
            fastq::Reader::from_file("example/simulated_reads/test.bwa.read1.fastq").unwrap();

        // Create an FM-Index for the reference genome
        // TODO: Use FMD-index instead, to not have to search two indices
        // TODO: Read about FMD index
        // TODO: It's perhaps worth to do a rank transform (memory reduction by 10* (?))
        let alphabet = alphabets::dna::n_alphabet();

        println!("Generate suffix array");
        let sa = suffix_array(&ref_seq);

        println!("Generate BWT");
        let bwt = bwt(&ref_seq, &sa);

        println!("Drop source sequence");
        drop(ref_seq);

        println!("Generate less");
        let less = less(&bwt, &alphabet);

        println!("Generate OCC");
        let occ = Occ::new(&bwt, 18, &alphabet);

        println!("Generate FM index");
        let fmindex = FMIndex::new(&bwt, &less, &occ);

        println!("Map reads");
        let interval_calculators = reads_fq_reader
            .records()
            .map(|pattern| fmindex.backward_search(pattern.unwrap().seq().iter()))
            .collect::<Vec<_>>();

        // Loop through the results, extracting the positions array for each pattern
        for interval_calculator in interval_calculators {
            let positions = interval_calculator.occ(&sa);
            println!("{:#?}", positions);
        }

        // There is actually only one sequence in the reference genome fasta file for now
        break;
    }
}
