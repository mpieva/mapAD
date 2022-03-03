use std::{
    fs::File,
    io::{self, Write},
    str::FromStr,
};

use noodles::{bam, sam};
use tempfile::tempdir;

use mapad::{
    index, map,
    mismatch_bounds::Discrete,
    sequence_difference_models::{LibraryPrep, SequenceDifferenceModel, SimpleAncientDnaModel},
    utils,
};

#[test]
fn integration_1() {
    let temp_dir = tempdir().unwrap();

    let test_genome_path = temp_dir.path().join("test_genome.fa");
    let input_bam_path = temp_dir.path().join("input_reads.bam");
    let output_bam_path = temp_dir.path().join("output_reads.bam");

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
        let sam_content = "\
        A00123_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t28M\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00234_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t589\t*\t0\t0\t28M\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00345_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t16\t*\t0\t0\t28M\t*\t0\t0\tCCTGGTCGTTCCCTAAGTTCATTGTTAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00456_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t16\t*\t0\t0\t28M\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00567_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t28M\t*\t0\t0\tCCTGGTCGTTCCCAAGTTCATTGTTAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00678_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t28M\t*\t0\t0\tCCTGGTCGTTCCCTTAAGTTCATTGTTAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00789_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t0\tchr1\t269\t37\t28M\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\tAS:i:0\tNM:i:0\tMD:Z:28\tXD:i:195";

        let input_bam_header =
            sam::Header::from_str("\
            @HD\tVN:1.0\n\
            @RG\tID:A12345\tSM:Sample1\n\
            @SQ\tSN:chr1\tLN:600\n\
            @PG\tID:samtools\tPN:samtools\tVN:1.13\tCL:samtools view -h interesting_specimen.bam -o input_reads.bam\n\
            @PG\tID:mapad\tPN:mapAD\tCL:mapad map\tPP:samtools\tDS:An aDNA aware short-read mapper\tVN:0.0.33\n\
            @PG\tID:mapad.1\tPN:mapAD\tCL:mapad map\tPP:mapad\tDS:An aDNA aware short-read mapper\tVN:0.0.33\n\
            ").unwrap();

        let mut input_bam_file = bam::Writer::new(File::create(&input_bam_path).unwrap());
        input_bam_file.write_header(&input_bam_header).unwrap();
        input_bam_file
            .write_reference_sequences(input_bam_header.reference_sequences())
            .unwrap();
        for sam_line in sam_content.lines() {
            let sam_record = sam::Record::from_str(sam_line).unwrap();
            input_bam_file
                .write_sam_record(input_bam_header.reference_sequences(), &sam_record)
                .unwrap();
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
            0.01,
            1.0,
            base_error_rate,
            false,
        );
        let representative_mm_penalty = adna_scoring_model.get_representative_mismatch_penalty();
        let mismatch_bound = Discrete::new(0.03, base_error_rate, representative_mm_penalty);
        let alignment_parameters = utils::AlignmentParameters {
            difference_model: adna_scoring_model.into(),
            mismatch_bound: mismatch_bound.into(),
            penalty_gap_open: representative_mm_penalty * 2.0,
            penalty_gap_extend: representative_mm_penalty,
            chunk_size: 1,
            gap_dist_ends: 5,
            stack_limit_abort: false,
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
        let mut bam_reader = bam::Reader::new(File::open(&output_bam_path).unwrap());

        // Move cursor to the right place
        let _header = bam_reader.read_header().unwrap();
        let _header_reference_sequences = bam_reader.read_reference_sequences().unwrap();

        let result_sample = bam_reader
            .records()
            .map(|maybe_record| {
                maybe_record.map(|record| {
                    (
                        record.reference_sequence_id().unwrap().into(),
                        record.position().unwrap().into(),
                        record.sequence().len(),
                        record.flags().bits(),
                        record.sequence().to_string(),
                        record.cigar().to_string(),
                    )
                })
            })
            .collect::<io::Result<Vec<_>>>()
            .unwrap();

        let comp = vec![
            (
                0,
                269,
                28,
                0,
                "TTAACAATGAACTTAGGGAACGACCAGG".into(),
                "28M".into(),
            ),
            (
                0,
                269,
                28,
                585,
                "TTAACAATGAACTTAGGGAACGACCAGG".into(),
                "28M".into(),
            ),
            (
                0,
                269,
                28,
                0,
                "TTAACAATGAACTTAGGGAACGACCAGG".into(),
                "28M".into(),
            ),
            (
                0,
                269,
                28,
                16,
                "TTAACAATGAACTTAGGGAACGACCAGG".into(),
                "28M".into(),
            ),
            (
                0,
                269,
                27,
                16,
                "TTAACAATGAACTTGGGAACGACCAGG".into(),
                "14M1D13M".into(),
            ),
            (
                0,
                269,
                29,
                16,
                "TTAACAATGAACTTAAGGGAACGACCAGG".into(),
                "15M1I13M".into(),
            ),
            (
                0,
                269,
                28,
                0,
                "TTAACAATGAACTTAGGGAACGACCAGG".into(),
                "28M".into(),
            ),
        ];
        assert_eq!(comp, result_sample);
    }

    temp_dir.close().unwrap();
}
