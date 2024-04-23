use std::{
    fs::File,
    io::{self, Write},
    path::{Path, PathBuf},
    thread,
};

use noodles::{
    bam,
    core::Position,
    sam::{
        self,
        alignment::{io::Write as AlignmentWrite, record::cigar, record_buf::Cigar},
    },
};
use tempfile::{tempdir, TempDir};

use mapad::{
    distributed::{dispatcher::Dispatcher, worker::Worker},
    index::indexing,
    map::{
        mapping,
        mismatch_bounds::Discrete,
        sequence_difference_models::{LibraryPrep, SequenceDifferenceModel, SimpleAncientDnaModel},
        AlignmentParameters,
    },
};

#[derive(Debug, Clone, PartialEq)]
struct BamFieldSubset {
    name: Option<sam::alignment::record_buf::Name>,
    flags: sam::alignment::record::Flags,
    tid: Option<usize>,
    pos: Option<Position>,
    mq: Option<sam::alignment::record::MappingQuality>,
    cigar: Cigar,
    seq_len: usize,
    seq: sam::alignment::record_buf::Sequence,
    qual: sam::alignment::record_buf::QualityScores,
    md: Option<String>,
    x0: Option<i32>,
    x1: Option<i32>,
    xa: Option<String>,
    xs: Option<f32>,
    xt: Option<char>,
}

#[derive(Debug)]
struct Environment {
    temp_dir: TempDir,
    alignment_parameters: AlignmentParameters,
    test_genome_path: PathBuf,
    input_bam_path: PathBuf,
}

fn prepare() -> Environment {
    let temp_dir = tempdir().unwrap();

    let test_genome_path = temp_dir.path().join("test_genome.fa");
    let input_bam_path = temp_dir.path().join("input_reads.bam");

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
ATTCTGATGACGAAGTGTATACCCTGGCGTGCTNGTCCCTCGGCGTTGGATATCCTAGAT
TGAGAATCCTGTCGCGGGACCTCGTTTAGGAAGCGAATGGTTGCACATCCGTCTAAACTA
>Chromosome_03
CCAAGAATCCGTAGACTCTGATCGATCATGCTAAAAATCGACCCAAGAATCCGTAGACTC
TGATCGATCATGCTAAAAATCGAT";
        writeln!(file, "{fasta_content}").unwrap();

        indexing::run(test_genome_path.to_str().unwrap(), 1234).unwrap();
    }

    // Create test reads and store them in a BAM file
    {
        let sam_content = b"\
        @HD\tVN:1.0\n\
        @RG\tID:A12345\tSM:Sample1\n\
        @SQ\tSN:chr1\tLN:600\n\
        @PG\tID:samtools\tPN:samtools\tVN:1.13\tCL:samtools view -h interesting_specimen.bam -o input_reads.bam\n\
        @PG\tID:mapAD\tPN:mapAD\tCL:mapad map\tPP:samtools\tDS:An aDNA aware short-read mapper\tVN:0.0.33\n\
        @PG\tID:mapAD.1\tPN:mapAD\tCL:mapad map\tPP:mapAD\tDS:An aDNA aware short-read mapper\tVN:0.0.33\n\
        A00123_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00234_0124_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t589\t*\t0\t0\t*\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00345_0125_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t16\t*\t0\t0\t28M\t*\t0\t0\tCCTGGTCGTTCCCTAAGTTCATTGTTAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00456_0126_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t16\t*\t0\t0\t28M\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00567_0127_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tCCTGGTCGTTCCCAAGTTCATTGTTAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00678_0128_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tCCTGGTCGTTCCCTTAAGTTCATTGTTAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\n\
        A00789_0129_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t0\tchr1\t269\t37\t28M\t*\t0\t0\tTTAACAATGAACTTAGGGAACGACCAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\tXI:Z:ACGACGT\tYI:Z::BBBBGG\tXJ:Z:TGCTGCA\tYJ:Z:AAAAABB\tFF:i:3\tZ0:i:0\tRG:Z:A12345\tAS:i:0\tNM:i:0\tMD:Z:28\tXD:i:195\n\
        A00789_0130_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tGATTGGTGCACGGACGCGCGTTGAAAGG\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]\n\
        A00791_0131_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tCCTCAT\t]]]]]]\n\
        A00792_0132_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tTCAAGAATCCGTAGACTCTGATCGATCATGCTAAAAATCGAT\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\n\
        A00793_0133_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tCTGGCGTGCTAGTCCCTCGGCG\t]]]]]]]]]]]]]]]]]]]]]]\n\
        A00794_0134_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tCGCCGAGGGACTAGCACGCCAG\t]]]]]]]]]]]]]]]]]]]]]]\n\
        A00795_0135_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678\t4\t*\t0\t0\t*\t*\t0\t0\tCGCCGAGGGACTAGCACCCCAG\t]]]]]]]]]]]]]]]]]]]]]]";

        let mut sam_reader = sam::io::Reader::new(&sam_content[..]);
        let input_sam_header = sam_reader.read_header().unwrap();

        let mut input_bam_file = bam::io::Writer::new(File::create(&input_bam_path).unwrap());
        input_bam_file.write_header(&input_sam_header).unwrap();
        for sam_record in sam_reader.records() {
            let sam_record = sam_record.unwrap();
            input_bam_file
                .write_alignment_record(&input_sam_header, &sam_record)
                .unwrap();
        }
    }

    // Mapping parameters
    {
        let base_error_rate = 0.02;
        let adna_scoring_model = SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.6,
                three_prime_overhang: 0.55,
            },
            0.01,
            1.0,
            base_error_rate / 3.0,
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
            max_num_gaps_open: 2,
        };

        Environment {
            temp_dir,
            alignment_parameters,
            test_genome_path,
            input_bam_path,
        }
    }
}

#[test]
fn integration_1_local() {
    let env = prepare();
    let output_bam_path_local = env.temp_dir.path().join("output_reads_local.bam");

    mapping::run(
        env.input_bam_path.to_str().unwrap(),
        env.test_genome_path.to_str().unwrap(),
        output_bam_path_local.to_str().unwrap(),
        &env.alignment_parameters,
    )
    .unwrap();

    check_results(&output_bam_path_local);
}

/// This test is disabled by default since it is flaky. Sometimes one of the workers can not
/// connect to the dispatcher because the mapping process is finished already, causing the test to
/// fail.
#[test]
#[ignore]
fn integration_1_distributed() {
    let env = prepare();
    let output_bam_path_distr = env.temp_dir.path().join("output_reads_distr.bam");
    let port = 4321;

    let dispatcher_handle = {
        let output_bam_path_distr_clone = output_bam_path_distr.clone();
        thread::spawn(move || {
            let mut dispatcher = Dispatcher::new(
                env.input_bam_path.to_str().unwrap(),
                env.test_genome_path.to_str().unwrap(),
                output_bam_path_distr_clone.to_str().unwrap(),
                &env.alignment_parameters,
            )
            .unwrap();
            dispatcher.run(port).unwrap();
        })
    };

    let w0_handle = thread::spawn(move || {
        let mut worker = Worker::new("127.0.0.1", port).unwrap();
        worker.run().unwrap();
    });
    let w1_handle = thread::spawn(move || {
        let mut worker = Worker::new("127.0.0.1", port).unwrap();
        worker.run().unwrap();
    });

    w0_handle.join().unwrap();
    w1_handle.join().unwrap();
    dispatcher_handle.join().unwrap();

    check_results(&output_bam_path_distr);
}

fn check_results<P>(bam_path: P)
where
    P: AsRef<Path>,
{
    let mut bam_reader = bam::io::Reader::new(File::open(bam_path).unwrap());

    // Check header
    {
        let header = bam_reader.read_header().unwrap();

        let header_prefix = "\
        @HD\tVN:1.6\tSO:unsorted\n\
        @SQ\tSN:chr1\tLN:600\n\
        @SQ\tSN:Chromosome_02\tLN:600\n\
        @SQ\tSN:Chromosome_03\tLN:84\n\
        @RG\tID:A12345\tSM:Sample1\n\
        @PG\tID:samtools\tPN:samtools\tVN:1.13\tCL:samtools view -h interesting_specimen.bam -o input_reads.bam\n\
        @PG\tID:mapAD\tPN:mapAD\tCL:mapad map\tPP:samtools\tDS:An aDNA aware short-read mapper";
        let mut header_writer = sam::io::Writer::new(Vec::new());
        header_writer.write_header(&header).unwrap();
        assert!(std::str::from_utf8(header_writer.get_ref())
            .unwrap()
            .starts_with(header_prefix));
    }

    let mut result_sample = bam_reader
        .records()
        .map(|maybe_record| {
            maybe_record.map(|record| BamFieldSubset {
                name: record
                    .name()
                    .map(|v| v.as_bytes().to_owned())
                    .map(Into::into),
                flags: record.flags().to_owned(),
                tid: record.reference_sequence_id().map(|v| v.unwrap()),
                pos: record.alignment_start().map(|v| v.unwrap()),
                mq: record.mapping_quality(),
                cigar: record.cigar().try_into().unwrap(),
                seq_len: record.sequence().len(),
                seq: record.sequence().iter().collect::<Vec<_>>().into(),
                qual: record.quality_scores().into(),
                md: record
                    .data()
                    .get(&sam::alignment::record::data::field::tag::Tag::MISMATCHED_POSITIONS)
                    .map(|v| match v {
                        Ok(sam::alignment::record::data::field::Value::String(value)) => {
                            value.to_string()
                        }
                        _ => {
                            panic!();
                        }
                    }),
                x0: record.data().get(b"X0").map(|v| match v {
                    Ok(sam::alignment::record::data::field::Value::Int32(value)) => value,
                    _ => {
                        panic!();
                    }
                }),
                x1: record
                    .data()
                    .get(b"X1")
                    .map(|value| value.unwrap().as_int().map(|v| v as i32).unwrap()),
                xa: record.data().get(b"XA").map(|v| match v {
                    Ok(sam::alignment::record::data::field::Value::String(value)) => {
                        value.to_string()
                    }
                    _ => {
                        panic!();
                    }
                }),
                xs: record.data().get(b"XS").map(|v| match v {
                    Ok(sam::alignment::record::data::field::Value::Float(value)) => value,
                    _ => {
                        panic!();
                    }
                }),
                xt: record.data().get(b"XT").map(|v| match v {
                    Ok(sam::alignment::record::data::field::Value::Character(value)) => {
                        value as char
                    }
                    _ => {
                        panic!();
                    }
                }),
            })
        })
        .collect::<io::Result<Vec<_>>>()
        .unwrap();

    let comp = vec![
        BamFieldSubset {
            name: Some(b"A00123_0123_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 0.into(),
            tid: Some(0_i32.try_into().unwrap()),
            pos: Some(269.try_into().unwrap()),
            mq: Some(37_u8.try_into().unwrap()),
            cigar: Cigar::from_iter([cigar::Op::new(cigar::op::Kind::Match, 28)]),
            seq_len: 28,
            seq: b"TTAACAATGAACTTAGGGAACGACCAGG".into(),
            qual: b"]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: Some("28".into()),
            x0: Some(1),
            x1: Some(0),
            xa: None,
            xs: None,
            xt: Some('U'),
        },
        BamFieldSubset {
            name: Some(b"A00234_0124_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 577.into(),
            tid: Some(0_i32.try_into().unwrap()),
            pos: Some(269.try_into().unwrap()),
            mq: Some(37_u8.try_into().unwrap()),
            cigar: Cigar::from_iter([cigar::Op::new(cigar::op::Kind::Match, 28)]),
            seq_len: 28,
            seq: b"TTAACAATGAACTTAGGGAACGACCAGG".into(),
            qual: b"]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: Some("28".into()),
            x0: Some(1),
            x1: Some(0),
            xa: None,
            xs: None,
            xt: Some('U'),
        },
        BamFieldSubset {
            name: Some(b"A00345_0125_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 0.into(),
            tid: Some(0_i32.try_into().unwrap()),
            pos: Some(269.try_into().unwrap()),
            mq: Some(37_u8.try_into().unwrap()),
            cigar: Cigar::from_iter([cigar::Op::new(cigar::op::Kind::Match, 28)]),
            seq_len: 28,
            seq: b"TTAACAATGAACTTAGGGAACGACCAGG".into(),
            qual: b"]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: Some("28".into()),
            x0: Some(1),
            x1: Some(0),
            xa: None,
            xs: None,
            xt: Some('U'),
        },
        BamFieldSubset {
            name: Some(b"A00456_0126_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 16.into(),
            tid: Some(0_i32.try_into().unwrap()),
            pos: Some(269.try_into().unwrap()),
            mq: Some(37_u8.try_into().unwrap()),
            cigar: Cigar::from_iter([cigar::Op::new(cigar::op::Kind::Match, 28)]),
            seq_len: 28,
            seq: b"TTAACAATGAACTTAGGGAACGACCAGG".into(),
            qual: b"]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: Some("28".into()),
            x0: Some(1),
            x1: Some(0),
            xa: None,
            xs: None,
            xt: Some('U'),
        },
        BamFieldSubset {
            name: Some(b"A00567_0127_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 16.into(),
            tid: Some(0_i32.try_into().unwrap()),
            pos: Some(269.try_into().unwrap()),
            mq: Some(20_u8.try_into().unwrap()),
            cigar: Cigar::from_iter([
                cigar::Op::new(cigar::op::Kind::Match, 14),
                cigar::Op::new(cigar::op::Kind::Deletion, 1),
                cigar::Op::new(cigar::op::Kind::Match, 13),
            ]),
            seq_len: 27,
            seq: b"TTAACAATGAACTTGGGAACGACCAGG".into(),
            qual: b"]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: Some("14^A13".into()),
            x0: Some(1),
            x1: Some(0),
            xa: None,
            xs: None,
            xt: Some('U'),
        },
        BamFieldSubset {
            name: Some(b"A00678_0128_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 16.into(),
            tid: Some(0_i32.try_into().unwrap()),
            pos: Some(269.try_into().unwrap()),
            mq: Some(20_u8.try_into().unwrap()),
            cigar: Cigar::from_iter([
                cigar::Op::new(cigar::op::Kind::Match, 15),
                cigar::Op::new(cigar::op::Kind::Insertion, 1),
                cigar::Op::new(cigar::op::Kind::Match, 13),
            ]),
            seq_len: 29,
            seq: b"TTAACAATGAACTTAAGGGAACGACCAGG".into(),
            qual: b"]]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: Some("28".into()),
            x0: Some(1),
            x1: Some(0),
            xa: None,
            xs: None,
            xt: Some('U'),
        },
        BamFieldSubset {
            name: Some(b"A00789_0129_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 0.into(),
            tid: Some(0_i32.try_into().unwrap()),
            pos: Some(269.try_into().unwrap()),
            mq: Some(37_u8.try_into().unwrap()),
            cigar: Cigar::from_iter([cigar::Op::new(cigar::op::Kind::Match, 28)]),
            seq_len: 28,
            seq: b"TTAACAATGAACTTAGGGAACGACCAGG".into(),
            qual: b"]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: Some("28".into()),
            x0: Some(1),
            x1: Some(0),
            xa: None,
            xs: None,
            xt: Some('U'),
        },
        BamFieldSubset {
            name: Some(b"A00789_0130_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 4.into(),
            tid: None,
            pos: None,
            mq: Some(0_u8.try_into().unwrap()),
            cigar: Cigar::default(),
            seq_len: 28,
            seq: b"GATTGGTGCACGGACGCGCGTTGAAAGG".into(),
            qual: b"]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: None,
            x0: None,
            x1: None,
            xa: None,
            xs: None,
            xt: None,
        },
        BamFieldSubset {
            name: Some(b"A00791_0131_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 0.into(),
            tid: Some(1),
            pos: Some(85.try_into().unwrap()),
            mq: Some(3_u8.try_into().unwrap()),
            cigar: Cigar::from_iter([cigar::Op::new(cigar::op::Kind::Match, 6)]),
            seq_len: 6,
            seq: b"CCTCAT".into(),
            qual: b"]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: Some("6".into()),
            x0: Some(2),
            x1: Some(0),
            xa: Some("Chromosome_02,+139,6M,6,0,2,0.00;".into()),
            xs: None,
            xt: Some('R'),
        },
        BamFieldSubset {
            name: Some(b"A00792_0132_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 0.into(),
            tid: Some(1),
            pos: Some(188.try_into().unwrap()),
            mq: Some(3_u8.try_into().unwrap()),
            cigar: Cigar::from_iter([cigar::Op::new(cigar::op::Kind::Match, 42)]),
            seq_len: 42,
            seq: b"TCAAGAATCCGTAGACTCTGATCGATCATGCTAAAAATCGAT".into(),
            qual: b"]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: Some("42".into()),
            x0: Some(1),
            x1: Some(2),
            xa: Some(
                "Chromosome_03,+43,42M,0C41,1,1,-0.72;Chromosome_03,+1,42M,0C40C0,2,1,-1.56;"
                    .into(),
            ),
            xs: Some(-0.7209588),
            xt: Some('U'),
        },
        BamFieldSubset {
            name: Some(b"A00793_0133_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 0.into(),
            tid: Some(1),
            pos: Some(504.try_into().unwrap()),
            mq: Some(37_u8.try_into().unwrap()),
            cigar: Cigar::from_iter([cigar::Op::new(cigar::op::Kind::Match, 22)]),
            seq_len: 22,
            seq: b"CTGGCGTGCTAGTCCCTCGGCG".into(),
            qual: b"]]]]]]]]]]]]]]]]]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: Some("10N11".into()),
            x0: Some(1),
            x1: Some(0),
            xa: None,
            xs: None,
            xt: Some('U'),
        },
        BamFieldSubset {
            name: Some(b"A00794_0134_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 16.into(),
            tid: Some(1),
            pos: Some(504.try_into().unwrap()),
            mq: Some(37_u8.try_into().unwrap()),
            cigar: Cigar::from_iter([cigar::Op::new(cigar::op::Kind::Match, 22)]),
            seq_len: 22,
            seq: b"CTGGCGTGCTAGTCCCTCGGCG".into(),
            qual: b"]]]]]]]]]]]]]]]]]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: Some("10N11".into()),
            x0: Some(1),
            x1: Some(0),
            xa: None,
            xs: None,
            xt: Some('U'),
        },
        BamFieldSubset {
            name: Some(b"A00795_0135_ABC12XXXXX_ABcd_AB_CC_DE:1:2345:1234:5678".into()),
            flags: 16.into(),
            tid: Some(1),
            pos: Some(504.try_into().unwrap()),
            mq: Some(37_u8.try_into().unwrap()),
            cigar: Cigar::from_iter([cigar::Op::new(cigar::op::Kind::Match, 22)]),
            seq_len: 22,
            seq: b"CTGGGGTGCTAGTCCCTCGGCG".into(),
            qual: b"]]]]]]]]]]]]]]]]]]]]]]"
                .iter()
                .map(|&encoded| encoded - 33)
                .collect::<Vec<_>>()
                .into(),
            md: Some("4C5N11".into()),
            x0: Some(1),
            x1: Some(0),
            xa: None,
            xs: None,
            xt: Some('U'),
        },
    ];
    result_sample.sort_by_key(|k| k.name.clone().map(|v| v.as_ref().to_vec()));

    assert_eq!(result_sample, comp);
}
