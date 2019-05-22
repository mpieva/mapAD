use std::error::Error;
use std::fs::File;

use log::debug;
use rand::prelude::{Rng, SeedableRng, StdRng};
use rand::seq::SliceRandom;

use bio::alphabets::{dna, Alphabet};
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::suffix_array::suffix_array;

use bincode;
use bio::io::fasta;
use libflate::deflate::Encoder;

use crate::map::{FastaIdPosition, FastaIdPositions};

const DNA_UPPERCASE_ALPHABET: &[u8; 4] = b"ACGT";

pub fn run(reference_path: &str) -> Result<(), Box<Error>> {
    let mut rng: StdRng = SeedableRng::seed_from_u64(1234);

    let alphabet = dna::alphabet();

    // Index the genome
    index(reference_path, &alphabet, reference_path, &mut rng)?;

    Ok(())
}

fn index<T: Rng>(
    reference_path: &str,
    alphabet: &Alphabet,
    name: &str,
    rng: &mut T,
) -> Result<(), Box<Error>> {
    debug!("Read input reference sequence");
    let mut ref_seq = fasta::Reader::from_file(reference_path)?
        .records()
        .filter_map(Result::ok)
        .flat_map(|record| record.seq().to_ascii_uppercase())
        .map(|c| match c {
            b'N' => *DNA_UPPERCASE_ALPHABET.choose(rng).unwrap(),
            _ => c,
        })
        .collect::<Vec<_>>();

    {
        debug!("Map identifiers to positions");
        let mut end = 0;
        let identifier_position_map = FastaIdPositions::new(
            fasta::Reader::from_file(reference_path)?
                .records()
                .filter_map(Result::ok)
                .map(|record| {
                    end += record.seq().len();
                    FastaIdPosition {
                        start: end - record.seq().len() + 1,
                        end,
                        identifier: record.id().to_string(),
                    }
                })
                .collect::<Vec<_>>(),
        );

        debug!("Save position map to disk");
        let f_pi = File::create(format!("{}.tpi", name))?;
        let mut e_pi = Encoder::new(f_pi);
        bincode::serialize_into(&mut e_pi, &identifier_position_map)?;
        e_pi.finish();
    }

    debug!("Add reverse complement and sentinels to reference");
    ref_seq.extend_from_slice(b"$");
    let ref_seq_rev_compl = dna::revcomp(&ref_seq);
    ref_seq.extend_from_slice(&ref_seq_rev_compl);
    drop(ref_seq_rev_compl);
    ref_seq.extend_from_slice(b"$");

    debug!("Generate suffix array");
    let suffix_array = suffix_array(&ref_seq);

    debug!("Generate BWT");
    let bwt = bwt(&ref_seq, &suffix_array);

    {
        debug!("Save suffix array to disk");
        let f_suffix_array = File::create(format!("{}.tsa", name))?;
        let mut e_suffix_array = Encoder::new(f_suffix_array);
        bincode::serialize_into(&mut e_suffix_array, &suffix_array)?;
        e_suffix_array.finish();
    }

    debug!("Drop suffix array");
    drop(suffix_array);

    debug!("Generate \"C\" table");
    let less = less(&bwt, &alphabet);

    debug!("Generate \"Occ\" table");
    let occ = Occ::new(&bwt, 128, &alphabet);

    debug!("Save BWT to disk");
    let f_bwt = File::create(format!("{}.tbw", name))?;
    let mut e_bwt = Encoder::new(f_bwt);
    bincode::serialize_into(&mut e_bwt, &bwt)?;
    e_bwt.finish();

    debug!("Save \"C\" table to disk");
    let f_less = File::create(format!("{}.tle", name))?;
    let mut e_less = Encoder::new(f_less);
    bincode::serialize_into(&mut e_less, &less)?;
    e_less.finish();

    debug!("Save \"Occ\" table to disk");
    let f_occ = File::create(format!("{}.toc", name))?;
    let mut e_occ = Encoder::new(f_occ);
    bincode::serialize_into(&mut e_occ, &occ)?;
    e_occ.finish();

    Ok(())
}
