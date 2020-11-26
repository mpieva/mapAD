use mapad::{
    index, map,
    mismatch_bounds::Discrete,
    sequence_difference_models::{LibraryPrep, SequenceDifferenceModel, SimpleAncientDnaModel},
    utils,
};

use rust_htslib::{bam, bam::Read};
use std::{fs::File, io::Write};
use tempfile::tempdir;

#[test]
fn integration_1() {
    let temp_dir = tempdir().unwrap();

    let test_genome_path = temp_dir.path().clone().join("test_genome.fa");
    let input_bam_path = temp_dir.path().clone().join("input_reads.bam");
    let output_bam_path = temp_dir.path().clone().join("output_reads.bam");

    // Create test genome and index it
    {
        let mut file = File::create(&test_genome_path).unwrap();
        let fasta_content = ">chr1
TGTACTCGGGTGCCGAAGCCTACAGCTGGACCACCCGATGGCGTGCCTCTATCGGCACTC
GGCAGAATTGTTCCGGACGTATTGCAACTCCTCCGTACTTTGGTCCGTAAACTCACTTAG
CTACCCTGTCACCCCTGCGGTATTTAAAAGGCCTAAGCTGATCTTGCACGTGAGAGCCTC
GCGTCTTGTGAGAAAAAGGTCCGGAAGTAATGGTTTGACACGATCAACGCCCGTCACGCC
GTATGGTCTGCTTAGCCCAACTAGAGTTTTAACAATGAACTTAGGGAACGACCAGGGAAC
ATATGCGACGTAAGAATGTTTGCCAGCCTCAGTAATTTGCAGGGGATAGTCTCCATTAGA
GCTTCCGGGTGGACATTTTTCGTGTCACTTGCCCCGACAAGCGACTAGCGTGTAGAGGGA
CAAAAGTCACAGGATTCCCAGGCATCTCTACTCCATAAGACTTTGTCACGAACTCATTAG
ACCTATGTCGCGACTACCCATGTATGGGCTCGCACCCTTCATGATTCTGCGCTGACCCTA
GGATGCCGAGTAGCACTTCCGCTGTGTATGTGGGGTTAGACCGAACACTAAGACCTTCAG
>Chromosome_02
CAGTGATGAAATGCCAAAGTCTAGGTTGGGGGAATAGGGCCGCGCCCTCTCCAGCGGCTC
TATGGCCGGACAATTTCGGACAGGCCTCATACAGGGTTCAAAGGTCAGGCCACGCGGGCT
GATCTTCCCTTCTGAGGCCCTCATGTATGTACTAAATAGCTAACGCTATGACTCGGCGTT
TAATACTTCAAGAATCCGTAGACTCTGATCGATCATGCTAAAAATCGATCGAGCATCAAC
TCCAATTGGAGGTCTTTACATTAGGACCTGACTCACTACGTACGCTGTGGTACATAATAG
CGATACTCATCGTCCAAGTTCAACGTGGGTAACAACCCTACTGGCTCCCCCGAATAGTAG
TACCAGGACGGGCTCAACAATACTGGAAGTAACGGAATTTTTTGCCGTAATTCTCAAAAT
AAAGAGGTAATTGACCGAAAACCCTGTAACTCACCAATATGGGTTGGCAATCTTACCAAA
ATTCTGATGACGAAGTGTATACCCTGGCGTGCTAGTCCCTCGGCGTTGGATATCCTAGAT
TGAGAATCCTGTCGCGGGACCTCGTTTAGGAAGCGAATGGTTGCACATCCGTCTAAACTA";
        writeln!(file, "{}", fasta_content).unwrap();

        index::run(test_genome_path.to_str().unwrap(), 1234).unwrap();
    }

    // Create test reads and store them in a BAM file
    {
        let sam_content = b"\
        A00123_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00234_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t589\t*\t0\t0\t*\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00345_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t16\t*\t0\t0\t*\t*\t0\t0\tCCTGGTCGTTCCCTAAGTTCATTGTTAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00456_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t16\t*\t0\t0\t*\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00567_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tCCTGGTCGTTCCCAAGTTCATTGTTAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00678_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tCCTGGTCGTTCCCTTAAGTTCATTGTTAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00789_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t0\tchr1\t269\t37\t28M\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\tAS:i:0\tNM:i:0\tMD:Z:28\tXD:i:195";

        let input_bam_header = bam::Header::from_template(&bam::HeaderView::from_bytes(
            b"@HD\tVN:1.0\n@RG\tID:A12345\n@SQ\tSN:chr1\tLN:600\n",
        ));
        let mut input_bam_file =
            bam::Writer::from_path(&input_bam_path, &input_bam_header, bam::Format::BAM).unwrap();
        for sam_line in sam_content.split(|&byte| byte == b'\n') {
            let bam_record = bam::Record::from_sam(input_bam_file.header(), sam_line).unwrap();
            input_bam_file.write(&bam_record).unwrap();
        }
    }

    // Start mapping
    {
        let base_error_rate = 0.02;
        let adna_scoring_model = SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.6,
                three_prime_overhang: 0.55,
            },
            1.0,
            0.01,
            base_error_rate,
            false,
        );
        let representative_mm_penalty = adna_scoring_model.get_representative_mismatch_penalty();
        let alignment_parameters = utils::AlignmentParameters {
            difference_model: adna_scoring_model.into(),
            mismatch_bound: Discrete::new(0.03, base_error_rate, representative_mm_penalty).into(),
            penalty_gap_open: representative_mm_penalty * 3.0,
            penalty_gap_extend: representative_mm_penalty,
            chunk_size: 1,
        };
        map::run(
            input_bam_path.to_str().unwrap(),
            test_genome_path.to_str().unwrap(),
            output_bam_path.to_str().unwrap(),
            &alignment_parameters,
        )
        .unwrap();
    }

    // Check mappings
    {
        let mut bam_reader = bam::Reader::from_path(&output_bam_path).unwrap();

        // println!(
        //     "{}",
        //     String::from_utf8(bam_reader.header().as_bytes().to_owned()).unwrap()
        // );

        // for processed_read in bam_reader.records() {
        //     let processed_read = processed_read.unwrap();
        //     print!("Chr: {}\t", processed_read.tid());
        //     print!("Pos: {}\t", processed_read.pos());
        //     print!("Flags: {}\t", processed_read.flags());
        //     println!(
        //         "Seq: {}",
        //         String::from_utf8(processed_read.seq().as_bytes()).unwrap()
        //     );
        // }

        let result_sample = bam_reader
            .records()
            .map(|record| record.unwrap())
            .map(|record| {
                (
                    record.tid(),
                    record.pos(),
                    record.seq().len(),
                    record.flags(),
                    record.seq().as_bytes(),
                )
            })
            .collect::<Vec<_>>();

        let comp = vec![
            (0, 268, 28, 0, b"TTAACAATGAACTTAGGGAACGACCAGG".to_vec()),
            (0, 268, 28, 585, b"TTAACAATGAACTTAGGGAACGACCAGG".to_vec()),
            (0, 268, 28, 0, b"TTAACAATGAACTTAGGGAACGACCAGG".to_vec()),
            (0, 268, 28, 16, b"TTAACAATGAACTTAGGGAACGACCAGG".to_vec()),
            (0, 268, 27, 16, b"TTAACAATGAACTTGGGAACGACCAGG".to_vec()),
            (0, 268, 29, 16, b"TTAACAATGAACTTAAGGGAACGACCAGG".to_vec()),
            (0, 268, 28, 0, b"TTAACAATGAACTTAGGGAACGACCAGG".to_vec()),
        ];
        assert_eq!(comp, result_sample);
    }

    temp_dir.close().unwrap();
}
