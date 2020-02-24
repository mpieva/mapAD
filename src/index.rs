use std::{error::Error, fs::File};

use log::debug;
use rand::{
    prelude::{Rng, SeedableRng, StdRng},
    seq::SliceRandom,
};

use bio::{
    alphabets::{dna, Alphabet},
    data_structures::{
        bwt::{bwt, less, Occ},
        suffix_array::suffix_array,
    },
    io::fasta,
};

use bincode;
use snap;

use crate::map::{FastaIdPosition, FastaIdPositions};

const DNA_UPPERCASE_ALPHABET: &[u8; 4] = b"ACGT";
const DNA_PURINES: &[u8; 2] = b"AG";
const DNA_PYRIMIDINES: &[u8; 2] = b"CT";
const DNA_KETONE: &[u8; 2] = b"GT";
const DNA_AMINO: &[u8; 2] = b"AC";
const DNA_STRONG: &[u8; 2] = b"CG";
const DNA_WEAK: &[u8; 2] = b"AT";
const DNA_NOT_A: &[u8; 3] = b"CGT";
const DNA_NOT_C: &[u8; 3] = b"AGT";
const DNA_NOT_G: &[u8; 3] = b"ACT";
const DNA_NOT_T: &[u8; 3] = b"ACG";

/// Entry point function to launch the indexing process
pub fn run(reference_path: &str, seed: u64) -> Result<(), Box<dyn Error>> {
    let mut rng: StdRng = SeedableRng::seed_from_u64(seed);

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
) -> Result<(), Box<dyn Error>> {
    debug!("Read input reference sequence");
    let mut ref_seq = fasta::Reader::from_file(reference_path)?
        .records()
        // Convert all bases to uppercase
        .flat_map(|record| {
            record
                .expect("Failed reading input file")
                .seq()
                .to_ascii_uppercase()
        })
        // Replace ambiguous bases with random ones
        .map(|c| match c {
            b'U' => b'T',
            b'R' => *DNA_PURINES.choose(rng).unwrap(),
            b'Y' => *DNA_PYRIMIDINES.choose(rng).unwrap(),
            b'K' => *DNA_KETONE.choose(rng).unwrap(),
            b'M' => *DNA_AMINO.choose(rng).unwrap(),
            b'S' => *DNA_STRONG.choose(rng).unwrap(),
            b'W' => *DNA_WEAK.choose(rng).unwrap(),
            b'B' => *DNA_NOT_A.choose(rng).unwrap(),
            b'D' => *DNA_NOT_C.choose(rng).unwrap(),
            b'H' => *DNA_NOT_G.choose(rng).unwrap(),
            b'V' => *DNA_NOT_T.choose(rng).unwrap(),
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
                .map(|record| {
                    let record = record.expect("Failed reading input file");
                    end += record.seq().len();
                    FastaIdPosition {
                        start: end - record.seq().len(),
                        end: end - 1,
                        identifier: record.id().to_string(),
                    }
                })
                .collect::<Vec<_>>(),
        );

        debug!("Save position map to disk");
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{}.tpi", name))?);
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
        let mut writer_suffix_array =
            snap::write::FrameEncoder::new(File::create(format!("{}.tsa", name))?);
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
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{}.tbw", name))?);
        bincode::serialize_into(&mut writer, &bwt)?;
    }

    {
        debug!("Save \"C\" table to disk");
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{}.tle", name))?);
        bincode::serialize_into(&mut writer, &less)?;
    }

    {
        debug!("Save \"Occ\" table to disk");
        let mut writer = snap::write::FrameEncoder::new(File::create(format!("{}.toc", name))?);
        bincode::serialize_into(&mut writer, &occ)?;
    }

    Ok(())
}
