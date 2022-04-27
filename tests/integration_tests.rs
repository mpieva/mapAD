use std::{
    fs::File,
    io::{self, Write},
};

use noodles::{
    bam,
    sam::{self, AlignmentRecord, AlignmentWriter},
};
use tempfile::tempdir;

use mapad::{
    index::indexing,
    map::{
        mapping,
        mismatch_bounds::Discrete,
        sequence_difference_models::{LibraryPrep, SequenceDifferenceModel, SimpleAncientDnaModel},
        AlignmentParameters,
    },
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

        indexing::run(test_genome_path.to_str().unwrap(), 1234).unwrap();
    }

    // Create test reads and store them in a BAM file
    {
        let sam_content = "\
        A00123_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00234_0124_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t589\t*\t0\t0\t*\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00345_0125_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t16\t*\t0\t0\t28M\t*\t0\t0\tCCTGGTCGTTCCCTAAGTTCATTGTTAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00456_0126_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t16\t*\t0\t0\t28M\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00567_0127_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tCCTGGTCGTTCCCAAGTTCATTGTTAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00678_0128_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tCCTGGTCGTTCCCTTAAGTTCATTGTTAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00789_0129_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t0\tchr1\t269\t37\t28M\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\tAS:i:0\tNM:i:0\tMD:Z:28\tXD:i:195\n\
        A00789_0130_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tGATTGGTGCACGGACGCGCGTTGAAAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]";

        let input_bam_header =
            "\
            @HD\tVN:1.0\n\
            @RG\tID:A12345\tSM:Sample1\n\
            @SQ\tSN:chr1\tLN:600\n\
            @PG\tID:samtools\tPN:samtools\tVN:1.13\tCL:samtools view -h interesting_specimen.bam -o input_reads.bam\n\
            @PG\tID:mapAD\tPN:mapAD\tCL:mapad map\tPP:samtools\tDS:An aDNA aware short-read mapper\tVN:0.0.33\n\
            @PG\tID:mapAD.1\tPN:mapAD\tCL:mapad map\tPP:mapAD\tDS:An aDNA aware short-read mapper\tVN:0.0.33\n\
            ".parse().unwrap();

        let mut input_bam_file = bam::Writer::new(File::create(&input_bam_path).unwrap());
        input_bam_file.write_header(&input_bam_header).unwrap();
        input_bam_file
            .write_reference_sequences(input_bam_header.reference_sequences())
            .unwrap();
        for sam_line in sam_content.lines() {
            let sam_record = sam_line.parse::<sam::Record>().unwrap();
            input_bam_file
                .write_alignment_record(&input_bam_header, &sam_record)
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
        let alignment_parameters = AlignmentParameters {
            difference_model: adna_scoring_model.into(),
            mismatch_bound: mismatch_bound.into(),
            penalty_gap_open: representative_mm_penalty * 1.5,
            penalty_gap_extend: representative_mm_penalty * 0.5,
            chunk_size: 1,
            gap_dist_ends: 5,
            stack_limit_abort: false,
        };
        mapping::run(
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
                        record.read_name().cloned(),
                        record.flags().to_owned(),
                        record.reference_sequence_id(),
                        record.position(),
                        record.mapping_quality(),
                        record.cigar().to_owned(),
                        record.sequence().len(),
                        record.sequence().to_owned(),
                        record.quality_scores().to_owned(),
                    )
                })
            })
            .collect::<io::Result<Vec<_>>>()
            .unwrap();

        let comp = vec![
            (
                Some(
                    "A00123_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678"
                        .parse()
                        .unwrap(),
                ),
                0.into(),
                Some(0_i32.try_into().unwrap()),
                Some(269.try_into().unwrap()),
                Some(37_u8.try_into().unwrap()),
                "28M".parse().unwrap(),
                28,
                "TTAACAATGAACTTAGGGAACGACCAGG".parse().unwrap(),
                "]]]]]]]]]]]]]]]]]]]]]]]]]]]]".parse().unwrap(),
            ),
            (
                Some(
                    "A00234_0124_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678"
                        .parse()
                        .unwrap(),
                ),
                585.into(),
                Some(0_i32.try_into().unwrap()),
                Some(269.try_into().unwrap()),
                Some(37_u8.try_into().unwrap()),
                "28M".parse().unwrap(),
                28,
                "TTAACAATGAACTTAGGGAACGACCAGG".parse().unwrap(),
                "]]]]]]]]]]]]]]]]]]]]]]]]]]]]".parse().unwrap(),
            ),
            (
                Some(
                    "A00345_0125_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678"
                        .parse()
                        .unwrap(),
                ),
                0.into(),
                Some(0_i32.try_into().unwrap()),
                Some(269.try_into().unwrap()),
                Some(37_u8.try_into().unwrap()),
                "28M".parse().unwrap(),
                28,
                "TTAACAATGAACTTAGGGAACGACCAGG".parse().unwrap(),
                "]]]]]]]]]]]]]]]]]]]]]]]]]]]]".parse().unwrap(),
            ),
            (
                Some(
                    "A00456_0126_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678"
                        .parse()
                        .unwrap(),
                ),
                16.into(),
                Some(0_i32.try_into().unwrap()),
                Some(269.try_into().unwrap()),
                Some(37_u8.try_into().unwrap()),
                "28M".parse().unwrap(),
                28,
                "TTAACAATGAACTTAGGGAACGACCAGG".parse().unwrap(),
                "]]]]]]]]]]]]]]]]]]]]]]]]]]]]".parse().unwrap(),
            ),
            (
                Some(
                    "A00567_0127_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678"
                        .parse()
                        .unwrap(),
                ),
                16.into(),
                Some(0_i32.try_into().unwrap()),
                Some(269.try_into().unwrap()),
                Some(3_u8.try_into().unwrap()),
                "14M1D13M".parse().unwrap(),
                27,
                "TTAACAATGAACTTGGGAACGACCAGG".parse().unwrap(),
                "]]]]]]]]]]]]]]]]]]]]]]]]]]]".parse().unwrap(),
            ),
            (
                Some(
                    "A00678_0128_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678"
                        .parse()
                        .unwrap(),
                ),
                16.into(),
                Some(0_i32.try_into().unwrap()),
                Some(269.try_into().unwrap()),
                Some(3_u8.try_into().unwrap()),
                "15M1I13M".parse().unwrap(),
                29,
                "TTAACAATGAACTTAAGGGAACGACCAGG".parse().unwrap(),
                "]]]]]]]]]]]]]]]]]]]]]]]]]]]]]".parse().unwrap(),
            ),
            (
                Some(
                    "A00789_0129_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678"
                        .parse()
                        .unwrap(),
                ),
                0.into(),
                Some(0_i32.try_into().unwrap()),
                Some(269.try_into().unwrap()),
                Some(37_u8.try_into().unwrap()),
                "28M".parse().unwrap(),
                28,
                "TTAACAATGAACTTAGGGAACGACCAGG".parse().unwrap(),
                "]]]]]]]]]]]]]]]]]]]]]]]]]]]]".parse().unwrap(),
            ),
            (
                Some(
                    "A00789_0130_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678"
                        .parse()
                        .unwrap(),
                ),
                4.into(),
                None,
                None,
                Some(0_u8.try_into().unwrap()),
                sam::record::Cigar::default(),
                28,
                "GATTGGTGCACGGACGCGCGTTGAAAGG".parse().unwrap(),
                "]]]]]]]]]]]]]]]]]]]]]]]]]]]]".parse().unwrap(),
            ),
        ];
        assert_eq!(comp, result_sample);
    }

    temp_dir.close().unwrap();
}
