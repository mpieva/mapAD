use bincode::serialize_into;
use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::suffix_array::suffix_array;
use bio::io::fasta;
use libflate::deflate::Encoder;
use std::error::Error;
use std::fs::File;

pub fn run(
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
        let f_suffix_array = File::create("reference.sar")?;
        let mut e_suffix_array = Encoder::new(f_suffix_array);
        serialize_into(&mut e_suffix_array, &suffix_array)?;
        e_suffix_array.finish();

        debug!("Drop suffix array");
        drop(suffix_array);

        debug!("Drop source sequence");
        drop(ref_seq);

        debug!("Generate \"C\" table");
        let less = less(&bwt, &rank_alphabet);

        debug!("Generate \"Occ\" table");
        let occ = Occ::new(&bwt, 128, &rank_alphabet);

        debug!("Save BWT to disk");
        let f_bwt = File::create("reference.bwt")?;
        let mut e_bwt = Encoder::new(f_bwt);
        serialize_into(&mut e_bwt, &bwt)?;
        e_bwt.finish();
        debug!("Save \"C\" table to disk");
        let f_less = File::create("reference.less")?;
        let mut e_less = Encoder::new(f_less);
        serialize_into(&mut e_less, &less)?;
        e_less.finish();
        debug!("Save \"Occ\" table to disk");
        let f_occ = File::create("reference.occ")?;
        let mut e_occ = Encoder::new(f_occ);
        serialize_into(&mut e_occ, &occ)?;
        e_occ.finish();

        // There is actually only one sequence in the reference genome fasta file for now
        break;
    }
    Ok(())
}
