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
use snap;

use crate::map::{FastaIdPosition, FastaIdPositions};

const DNA_UPPERCASE_ALPHABET: &[u8; 4] = b"ACGT";

/// Entry point function to launch the indexing process
pub fn run(reference_path: &str) -> Result<(), Box<Error>> {
    let mut rng: StdRng = SeedableRng::seed_from_u64(1234);

    let alphabet = dna::alphabet();

    // Index the genome
    index(reference_path, &alphabet, reference_path, &mut rng)?;

    Ok(())
}

/// Index a given reference and write the index to disk.
/// Ambiguous bases ('N') are converted to random bases.
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
        // Convert all bases to uppercase
        .flat_map(|record| record.seq().to_ascii_uppercase())
        // Replace ambiguous bases with random ones
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
        let mut writer = snap::Writer::new(File::create(format!("{}.tpi", name))?);
        bincode::serialize_into(&mut writer, &identifier_position_map)?;
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
        let mut writer_suffix_array = snap::Writer::new(File::create(format!("{}.tsa", name))?);
        bincode::serialize_into(&mut writer_suffix_array, &suffix_array)?;
    }

    debug!("Drop suffix array");
    drop(suffix_array);

    debug!("Generate \"C\" table");
    let less = less(&bwt, &alphabet);

    debug!("Generate \"Occ\" table");
    let occ = Occ::new(&bwt, 128, &alphabet);

    {
        debug!("Save BWT to disk");
        let mut writer = snap::Writer::new(File::create(format!("{}.tbw", name))?);
        bincode::serialize_into(&mut writer, &bwt)?;
    }

    {
        debug!("Save \"C\" table to disk");
        let mut writer = snap::Writer::new(File::create(format!("{}.tle", name))?);
        bincode::serialize_into(&mut writer, &less)?;
    }

    {
        debug!("Save \"Occ\" table to disk");
        let mut writer = snap::Writer::new(File::create(format!("{}.toc", name))?);
        bincode::serialize_into(&mut writer, &occ)?;
    }

    Ok(())
}
