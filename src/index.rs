use std::error::Error;
use std::fs::File;

use log::debug;

use bio::alphabets::dna::{n_alphabet, revcomp};
use bio::alphabets::Alphabet;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::suffix_array::suffix_array;

use bincode::serialize_into;
use bio::io::fasta;
use libflate::deflate::Encoder;

pub fn run(reference_path: &str) -> Result<(), Box<Error>> {
    let alphabet = n_alphabet();

    debug!("Read input reference sequence");
    let mut ref_seq = fasta::Reader::from_file(reference_path)
        .unwrap()
        .records()
        .flat_map(|record| record.unwrap().seq().to_ascii_uppercase())
        .collect::<Vec<_>>();

    // Index the genome
    index(&mut ref_seq, &alphabet, "ref", true)?;

    Ok(())
}

fn index(
    ref_seq: &mut Vec<u8>,
    alphabet: &Alphabet,
    name: &str,
    save_suffix_array: bool,
) -> Result<(), Box<Error>> {
    debug!("Add reverse complement and sentinels to reference");
    let ref_seq_rev_compl = revcomp(ref_seq.as_slice());
    ref_seq.extend_from_slice(b"$");
    ref_seq.extend_from_slice(&ref_seq_rev_compl);
    drop(ref_seq_rev_compl);
    ref_seq.extend_from_slice(b"$");

    debug!("Generate suffix array");
    let suffix_array = suffix_array(ref_seq);

    debug!("Generate BWT");
    let bwt = bwt(&ref_seq, &suffix_array);

    if save_suffix_array {
        debug!("Save suffix array to disk");
        let f_suffix_array = File::create(format!("{}.sar", name))?;
        let mut e_suffix_array = Encoder::new(f_suffix_array);
        serialize_into(&mut e_suffix_array, &suffix_array)?;
        e_suffix_array.finish();
    }

    debug!("Drop suffix array");
    drop(suffix_array);

    debug!("Generate \"C\" table");
    let less = less(&bwt, &alphabet);

    debug!("Generate \"Occ\" table");
    let occ = Occ::new(&bwt, 128, &alphabet);

    debug!("Save BWT to disk");
    let f_bwt = File::create(format!("{}.bwt", name))?;
    let mut e_bwt = Encoder::new(f_bwt);
    serialize_into(&mut e_bwt, &bwt)?;
    e_bwt.finish();

    debug!("Save \"C\" table to disk");
    let f_less = File::create(format!("{}.less", name))?;
    let mut e_less = Encoder::new(f_less);
    serialize_into(&mut e_less, &less)?;
    e_less.finish();

    debug!("Save \"Occ\" table to disk");
    let f_occ = File::create(format!("{}.occ", name))?;
    let mut e_occ = Encoder::new(f_occ);
    serialize_into(&mut e_occ, &occ)?;
    e_occ.finish();

    Ok(())
}
