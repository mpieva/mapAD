use bincode::serialize_into;
use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::suffix_array::suffix_array;
use bio::io::fasta;
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
        let mut f_suffix_array = File::create("reference.sar")?;
        serialize_into(f_suffix_array, &suffix_array)?;

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
        serialize_into(f_less, &less)?;
        let mut f_bwt = File::create("reference.bwt")?;
        serialize_into(f_bwt, &bwt)?;
        let mut f_occ = File::create("reference.occ")?;
        serialize_into(f_occ, &occ)?;

        // There is actually only one sequence in the reference genome fasta file for now
        break;
    }
    Ok(())
}
