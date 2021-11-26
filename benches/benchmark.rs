use bio::alphabets;
use criterion::{criterion_group, criterion_main, Criterion};
use min_max_heap::MinMaxHeap;

use mapad::{
    backtrack_tree::Tree,
    map::k_mismatch_search,
    mismatch_bounds::*,
    sequence_difference_models::*,
    utils::{build_auxiliary_structures, AlignmentParameters},
};

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("3_mismatch_search", |b| {
        let ref_seq = "GATTACA".as_bytes().to_owned();

        let difference_model = SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            },
            0.001,
            0.9,
            0.02 / 3.0,
            false,
        );

        let representative_mismatch_penalty =
            difference_model.get_representative_mismatch_penalty();

        let mismatch_bound =
            MismatchBoundDispatch::from(Discrete::new(0.02, 0.02, representative_mismatch_penalty));

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound,
            penalty_gap_open: 0.00001_f32.log2(),
            penalty_gap_extend: representative_mismatch_penalty,
            chunk_size: 1,
            gap_dist_ends: 5,
        };

        let alphabet = alphabets::Alphabet::new(mapad::index::DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GTTT".as_bytes().to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack = MinMaxHeap::new();
        let mut tree = Tree::new();
        b.iter(|| {
            k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
            );
        })
    });
}

fn bench_exo_endo(c: &mut Criterion) {
    let ref_seq = "AACATCTTTTTCGGCCGTTTCATAAATACTCTTTATTCACAATTAAGGTCTCTCTGTGCTTCTACGATGGTTGGCTGGCCTGCTACTATGTAGCCCAGTT\
TTGTACTTGAAACAGATGAGTCTCTTCACACTCATACGCGTTTCAGCAAGTACCAGTGATAGCGGTGGTGCTCCTCAATATCCGCCGCTCGCAAGCGACC\
GAGAAAAGCACGCATCGAGAACCTTGTCTTACACAACTTTCAATGCTAGCGCCTAAGTTCAACCGCTAATATAACAGCTGAGTCCCCATCGCCCGAGGGT\
CGCATAGTTACTCCGTGAGGTGTATCAGAAATCCGCGATTTCGTAACTTCTACATCCCGTCCTCTCAGACCACCATTTTCCGTCATTACCAGCAAATATG\
ATCCGCCGATCGGCCGGACATCTGCCAAACCTTGCTGTCCTCCAATAAGTTTAGGTTAAGGGCCCAGTCAAGAAAAGACCCCAGTTCAAATGACGTACAA\
GGCAGCACCCCGTTTATACACCGGGGCATGCTTAGCTGCTATCGGACTTGCGCGGATAAAAATAGCCTAACTGCTAAAAGCATGTCCGCCCACGATGGAC\
TGGGGATTAACAGATAAGTTTAACATCTTAAACGTGAGTTGCGACAACAGCCAGGCTGACCGTACTGCCCGGACGTGCGGGGTCTGCTCATTGTGAATTC\
AGAGATGCCTCTATATTAGATTGGAGGGTAACTCAATAATCATGAGCGCAAGCTCACGAGAAGTACTGCTACTCCGCATATGTTGCTGAATATCTTCTGC\
CTGCGTGGCGGCCAGCATGTGATTACAATGACTCACTAGGGGGGGTAGTAGAGTCCGATTCGGTTGACTACTCGGGGTAACTAGGCATACAGTGGGCTCT\
AAGTCTGTGCGTCAAGCCCAGGGCAGGTCGTCGCTCATTAAGTTAGCCCACTTCATTAAGTATCACAGGGGCATCCTAGGTCCGTCTTAGTGTCCAACCA\
GCACCTCAATGTAATGTGGTGAGGGTAGTCACCAGTAGCGGAGTCCATCCTCTGGCCGCGTCGAGTCTTTCACGTCTATTCATGACGTATAGTCCCATTG\
GTGCGTCATGTACAGGCTTCTAGGACCACATTCTCCTTGGGCGGTTATCAAGCAAAATGCAAAGGTACCTCGGCATTCCGCCGATATTGTCTCTGAGATG\
TACGACGAGCCACACACTACGCTCGTTCTTCTGAGCCAATGGTTCGTATGGGATGCACTATAGTCGGACTAGACCCTTGCCATACTACATGCTGAAGTTT\
CAGGGATAACAATCGGTGGGCATAGACAACGCTAACGTCCGCGTTAAAGATTTAGTAACGTACGATACGCTGATTGCTGCAGGACAGAATAGACTTCAAT\
CGCGCTACCCTTCTAGCGCTGGGTAGTTTCGCGGACCGTATGCCAGGTAACAAGCACAGTAAACGCACTATTGGGTAATACGGATCTCAAAGGCTCATGC\
TCGGTGAAGGTGCCCGTAAGTGCACTTTCCTAACATAGCCGGCGTTATTAAGAGTTCCATCGGACAGCGCCCCACACCAAATCACTCAATTTCTCCTGTT\
CACCTTGCGCGTCATAATCTAACTGTGGCGAAGAACTAGGAACCTGAGAATAACCATTTCGATCGGTGATCCGTCAGATCGCTATTAGGAGCCAGAGACG\
TCTTCTACGGGGTCACCGTTCGTACCATGCTACCAGTCTTCGTCCCCTATGAGCATTTAACACAATCGTACTGTAGGTTGTCACACAACCGCTTAACGTT\
GTATCAGCTTTGCATACGCGAAGCAACGATGAGAACCGCTAAGGCTTCAATCTCCCGTTAAAGACCCTGAACTGGTCCGTGGCGGCGTGTTTTATCAGTG\
TACGGGATTCTTCGTGTACTTACTTGATAACGGCTCATACCAAGGCACGTGATTCTCTCGCAGAGGTCTCCCCCTTGGAAGTTATAGGCGCGCCGGTTTC\
CATTAGCATTGGACCACGACCGGGGGTGCTCCTATGTCCGACTCTCCATAGCGAGGAACATTGCTTCAAGCCCTTGGACGCTAGTTATTCATGCCTTGCG\
ATCTCCTAGCTTACCAAGGTTGTCGAACATCCGTTGAGAAAGTCATAAGACGTCAATAAACGTAACGGGGTTCTTTCTGCACTTGCGCTCAGACCTGCAT\
TTTAACCGTACAGGGAGCATCATACATCCCGTAACCACAGATATATCACTAACCAGGAATGACAAAGTCGGAAAGTAAGTCCCTACCGTCCGGAACCAGC\
GACGACCCAGGAAACTCCAGACTAGAAAGAGTACAGCAACCATGATATAAATAGTGCGAGCTGACAATGGTTTCGGCAGCACTGAATTTTTCGGGCTCAC\
AACAGGCGACAACCGGCCTTCGGAGAGTTTAATGCGCCATTCAGTTCGTCATCACGCTGGTACCGTCACCCTATGGGAGCTGCCAACTGTCACGATAGAA\
GATTATTTGTCACGGGACCGATTGCTTAAGACAATAACCAGGCTGGGACTACGTCATACGTCTCCAACTAACAAAAGTACATCCCCTATGAGATCAAAGG\
TACCGTATAAGAGAGTATGTCTAGGGCCTTGCATGGTTTTTACTTAACGCCGCTCCCGGAGCTTGCCCGCCATCCGGTCGAGAAAGCGTACCTCAAGTTT\
TGACCCTAACCAACACCCTAGAAGTTATTTGCTTCTAGGTTTCATACTTAAACTGCAATTGACACCGTACAAGCAGTCCGACCACTTTGGAAGCGTCGTG\
GATGCCCGGTTATACGGGCGATCGAGCTCTAGTCCTGACTCACGGCCGTGTATAACTTGATTCACGCCCCTTGCGCCGAACGCCCTATGGATCAAGGCGC\
CATCGAAGAGGTATTGTTACCATCAGAGGTTCGTGTGTCCTTTTAGCGCAACATTGATTACCGGTGCCTCCGTAAGTGGGGCCGACCAAAGTGCCGCCAG\
ATACTACCCCGCACGAAAATAGAAATGCGTTTGCGACCCTAACACTGAGCGTGCTGCCACTTGAAATTAATACTCTGAGGGCGCGGATCACTACCGGCGA\
GGCAAGATTACAATACGTACACTCCGTTAATACCATAACTCTTGGGTGTAGGAACTCCGGTAATACACCCGCAGGACGTTATTGGCCGTCGCCAGCAGCT\
TAGTAACAGATATAAAATCTAGGGACTTGTCGGGCAATGTCGGACGAATTCCGTTGTACCGTCAACTCCCTATACGCAAACTTGGTAACGTTTACTAGGT\
AGACGCGTTCAATGAACGGGTCTACTACTTCACTGGCGCCTCGTTGTGGCGTACTGCACCCGTAATACCCTATACCCGGCGGTATTCATCAGTCGCTCTG\
TCCACAGCAGCCAGGGAGGCTGTTTTCCTCGGGGCTACGGGGTCAACCCGATCCCATGGGTTTGCGTAGTACAATGAAAGGTCTCCTATCTTTCTTTTTG\
CGTACGGCTGCAGGGAAGACGCGGACATCATAGGCAACAGGGCCAGAGATATCCCCGGCTTCGAGCACTAGGGTGCTATTCCCCCAACTACGGCGCCAGG\
ACACGTGCCACTTGTTGCCACCTAACATCACTGTTGACCGTCAGATTGAAGAGTACTCAATGTCCACGATTTGCGCTGTGCCAGCGGCATCAGAACGGGA\
ATACGCAGGCATGCGTGAGAGCCAAATACCGCCCTCTGTCGGCAATCAGGTCGTAAAAGTTGAGATGTGAGGGTCTTAAGTTTGCGCCTTATGGTCTGGA\
CTCCGGTAGAACTCGGACGTTTAAGGTAACTTAGGAGGCCCTCACCTTCGAGACTGGCCGAGCAGTCAAAATCGGTCGGGCCGACGTCTTATAACTTGGC\
AAAGGTGAGGGGCACGTCTGTAATATGCAGCTGACATGTGGCTAAGAGTTGTTAAAATGTGTCGTTAAGTGTTAATACCTCACATTAGGTTTTAACGACC\
ACCTTGTTAGGTTTAAATCGGCAACTTGGCATGGACTATGTACTGGAAATCGACATAACGTGGTGAAGACTCCATACCTCGACGCATACGACCCCTCAAT\
ATCCCGCATCGTGTCCATTTATCGCTAAATAGGAAAGTGACCGACGCATACGACCGCCAGGCGCACACAAAGCTGAATCCGGCCAATAAAATAGTCTCAT\
GGTCCCGCCCGAGCCCAGTCTCGGCACCCCTAGACCCGAGTGCCGATCTCGTAGGGACACACCCTGATCGATTTTTTGCTGCATGATCTTTAAAAGGCAG\
ATGGGGGCTTAACAGTCACTACACTGACACATGAGACCTCGGAAGCGGAATCACTGCGTTCCAAACTCTCAGAGACCACGGCAAAGGAGATTTCCGGGGC\
GTAGGCCTTTTCGGGATACACTACGTATGGTTCAGCGCGACAGCCCGCATATACGGTGCATTTATATCAATCTTTAGGGAGTAACCCATCAGGTGAACAT\
GTCCAGCCCCCGGGGTCGCGAGCAATTAGCCGTGATCGGATAGGTCCAAACATGACTGATGTAGGGCCACCGGTGGAACCAGCCCGAGCAATGATTCGCC\
CCGGTGACTCCTACTCAAGCGGTCGGTGAGAGTATTGACCTTGTGACCAACGACGGGAGCAGCAAGCTGCGTGTTCCTTACCAACTCCGGGATATTCCTT\
TGATGTTCTCTAAGACAGACGTTGAGCATCAGCCCATATTAAGTCGAACGTGGTCCACTATCGATAAGTCTGCACAATTGGACGCTAGCAGGAATGGTGT\
GAGGACCCGAATCAGGCCCTCGAGTGGGACCTGGCGATCACTGCCGTCTCCGTTGCCGGTCGTTGTCCAACCGGTAAACCGGACAGCCCACATGTGCACG\
GTAAATGTTCCGTGCATAACATAATGCGAATGACCGCCTTCTTATGTCCGTTCCGAGGCCGGTCTTCTCACCCTTTTAGCATGTGGCAAGCCGGCACCGG\
AAGTGGACGACCGTCGCGTCTGTTCGTTGAGCGTGTATCTGGCGACTTAAACGATCCTACTAATCGCGTACGCTTGGTGCCTATAGAACCCGGGGAACTT\
GGCACGTCATGCCGGATACTTAATAGCCTGGCTCGTGTCGTGAGAACAATCTCCTGTGAAGTTACGTCTGAGTTTTAAAAGGCATCGTCGATCCCCTTCC\
CCTGCGCGCGGAGACGACTTAGCGACGTCGTGCACGACGAAAGCAGCCCTATTATAAGTCGTCCAGAAAGCCAGTATCTGATGTAAGCCAGACGCTAAGT\
TCAACCCAGTACATTTAAGACTAACCCATGTAAGCAAGTTATGTTGATCGCCCTCACAGCTATTCCCCGGCGGTGGACTCCATTCCATGCAGCTGATATA\
CAGATACATCACGCGTTAATTATGGTCTAGGTGTGTAAGCCGTTTCTAGGTTAGGCAGCAGTTAACGGCACTAGAGTATTGAAACACGGGAAACTAGGGG\
GCATCAATTTTCAGATGTCGTACCACAATTGCGAAATGGAACTGCAATTGTATCAAATGGCGCAATAGATATGTAGTGTCGTCCTGAGAATTCATTGGCA\
GTAAATCAACCTACTGAGCGGAGGATGTATACGCCAAATACATACATTAAAGTATCAATTATGCCCGGCCGAGTCTGGACAACCCGAGAAGCGTGATTTT\
GGGGAGTGGGGCGACGAGTCAGTAGCGTTTGGGCAGAACATGGCCTCCCTAGAGCCGCGACCGCGTGCTGGTACGGCACTCAAACTACATTTCTTTTGAG\
CGCGTAATATGAGAATTACACGGGTGCTAGATATGTAAGATACCACCAACTTCCCGCGTCCAAGGAACGGGAGACCACCTCCCGGGCCGTGCTGAGCATT\
GCCCGTTCAAGGCACCAGCAATTCTGGTGCCATAGTCGTAGAGCAGCCAACAACGGTGAAGAAGCTTACCGAAATTACAATGCGAGTGTGGCCTGCCGGC\
ACCTAAGTGTCAGCGGCCGGGGGGAGAGGATATTTTGCCGCATCCCGTCCCAAACGCGTGTTACGATCTACTGCAGATATTTGATTAGCGGTATCTTTTG\
CCCCTCCCCGCTTCAAGTTTTCATCCTAAGGGTGGGACACTTACAAGCTACAGGCCCCTTGGTAACAGGCTGTACACATCACACGTTGGAAAGCGTAATT\
ATAAAGGGTGTAGCCCTTATCTACTACTTCCAGACCGGAGTTTAACAGAGACCTGAAAGGACACTGCATTTACGAAGGGATGCCCGGGCCAACGAGAATG\
CGAACGGGTTATCGACTTAACGCGTAGTCTTCCGCTGCGGTAATTGTGTTGGATCTCCAGCTATACAGGCTGGAGCTCACTTGCCAGTCGGTAGGGTGAG\
CGGCCTTTATAGAGCGTGAGTTAGCATTGCCATTGAAATGAGATTGCTTGCCCAATACAAGTCTGTCTTCCAAGTAATAGGGCATTTGAGGCTTGAGGAA\
GCCACGCTCAGCTCTGGATGAAAAAAGTAACACGCTTCTATTAAGAAAAAGGTGAAGAAACTCCGCCTCGCAGTAGGTATCCAGCCGAACAGGGCAATCT\
CTGAGGGGACTTGTTAGACCGTCTCTGGGACTCCCATTGGCAGCAACAAAGGGATGGATCACATCTGAAGCGTCGCGGCTATTAAGAGTCTCAGCCCGGG\
CTATGCCGTATGTTTACCAAACCTCTGGTACCCATATTTCGTAGTGATGTCGCATAAAGCGCTCCTAGCAAGTTCGATCACACGTGACGAGGAGGCACGT\
TCCGCGTCATACAGACCGGGGGTGTTACTTCATTGGTACGACGGGGCACGACTCACCGAATGGGCCTGGTATTCTTCGGACGGTAGGTCCATCTACCATG\
CTCGCCACGTTTACACGGAACTTCGGGCAGCGACTCGCTAAGGGCTGCGTTAGAAACTACCTTGCAATCGGCCCAATGACTGGCAGAAATGAATCAGAAA\
CCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTAGGGTGACCGTAAGATTTAGGTCCATCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC\
TAACCGATTCTTAACTGCCAGTCCTCCGAAATTACCACCTGTTGGGGGGCATGATTCAGAATATGAAGTAGGATGGCCATGCGAAGAGGTGGCGCTTCAT\
GCTATCATTAGTACGGGATCTATAATCGAACACGATTTGTAATAGGGCGGGAGTCAAGCGCATATCCCACGTGGTTCCAACCGCTAGCCGCCTAACAAGT\
CTGTGATGACCAATATTGGAACCTTAGTACGTGGGGCCTTGAAAATCATTTGTTCGCTGAATAATTTACTGAAGATACCAGGTTTGTCATGATCAAGGTA\
GACCCTTGCCGACTCAATCCAATCTCCGGACCTAAAGGTGGCGATCTTAAGTCACCCACTTCGGCCTTGTTAATCGCGGCGCCCTGCATCCTTGGCTTTC\
AAGGAGGGGACATAAATGATTCGCTCGCACAAACGCAAGCTGGGCGTCTCAGGATTGCCGTTATTTGCAGCTTTAGCGGTACCATGCGGCGCCTAAATTT\
AACGAACTCTCGCCGGGGTCCCTCGTTCAATGTGACTATAGCTGGCTAAGGACGTTATCAGGCATCTCCAGGTCGTATCTCCTGGCTTGCTCTATTATGA\
CCTCCCGAGTTCATTACTGCCACCAAGCGTGTTATGCCAAGTAAGGACCAACGGGCCGCTTTAATGTACTCATATCGACGCGGTCTCTATGCTATCCACA\
ACGGCATTATTAGTAATCTACAACTATGCAATGAGTTGACCTTTGGTCCACGGACCAATAACCAGGGGGCTCTATACTGGCTCCACCCATCTTTTGCATG\
TTGGCCCAAATAATGCGCCCTTCCTAACTGGATACTTCTCCTTCCTTTGCCATACTGTAAACACGACTCGGCGGAGACGAGACATACCCAAGGTGGCCCA\
ATCCACAATTAACTAAAACCCAACCGCAGGTTAGCCGCCCTATATTGTTTCGTATCGCTTTAGTGTTGGAGGGTTATATCGCTGAACACATTATTATTGC\
GAACGTTCGGTTTCCTTCCGGTCGAAAAATAGTTACTTCATATGGGCCCCACGACCGATGGAGACCCACTGATACATTCTCGACTTGATCGGTCCGGGGT\
AAGGACAACAATGTGGCCAGGCCTTAATGCGATGGGACGATCATAGGGGTAGTATTGAGCAAGTTATACTTACTTACTTGCGGTCATTCGGGTTAAAGAG\
CTTTTGCTTTTGTGAGGCAACTGTTACGTGTATCGGACATGCATCCGATGGTTGCACCAATTTTAAAATTTGGACCTAATTGTGGGAGGCGGCTGGTTAT\
TCATGGAGCCCCTCTATATGGCACGTCGAACTCGTTAGTGCTTCGTCCGGTCCTAGACCAGGGAATACCGGTGAGCAAAGATACCCCACCCGATACGCCT\
CACACATCTCTTTCCGACTGTAACGAGCTGAAGCTAGAGTTCAGGTAAAGTCGGGCTAAGGACGTACATTAAAGCGGGGTACAATTCTGCACGCACCACC\
GTCTCAGATAGGCCGCCCAATCCAAGTTCGGGAAGTGATGGGATGACTACGGCTGCAATACTGTCTGGGATCTCTAACGCGATCACCCTCTCCATCCCTT\
AGGAAAGTCATCTTCTTTCCTAAACTGCCGCTGCTTATTGGACAACCGGTGACCGACAGCGAAACCATTGTTTGTGTAGTCTAGTAGACGACCAGCTGTC\
TGGGGTAATAACCAAAATCCGAATTCCACTATCCGGGGAGTCATAGCATACTGGGGTTAATCGAGCGCAGTCGAACGCAAAACTACACACCGAGTGAGGA\
CTCTCTCCAACGACCCTAATAGACTTCTGGAGCATCTTGCGAATAACGGACGGTATGTCGATGCTACCCCGCCTTCTGCTAGCCGGCATCGGCAACACGA\
TGAGCGTTCACGAAATTATTGAGAGAACAATTGCCTTAATACTAATTCGGGATATAGTTCAAAATGGACTGATATTCTCCACTATCTCTTACTAGACCCC\
ATTTCGGATAGAAGGCCATCGGAATAAAGCTGCACGAGGATCCGTCCCGACACAATCGGAGTAGAGGCTGGACATTCGATACAGGCCGTGGGCTAGCGAC\
CTCGCTCGAGCTCAAACGTTTCAGGCTGTCGGAACCCAGACGCCCCCCGAGGCCTCTGTGGTGGCAATATCCACAATAGTGCTTAGCGTCGATCCTCTGG\
CCATTCCCGTAACCCATCATAGTAACTGTGCGGATTCTCATCCTATGGGCAATAAAGCTCCAGGTCGTTCCTTATTTCCTCTCTAACTCCCCGTTGTGTC\
GTCGTTGGACGGCCTTCGTGCCGTGTTGGAATGTAAGATAAATCACGGGGCCATAATGCATTAATATTGGACGAGCTGTTCGCAAACACGTGAAGTAAAT\
AGAATTGTAATCGAAATTGTCAAGGACTATCTGATTACCAAACTACCGACTCCCCTTCACTGCGGGCAGCTCAGGTTGTCTGTACCGGGATGATGCAAGG\
CGTTCCCCATTCAGCTTGGAGGTGGACGCAGGTTGCCAATCGCCACTCGCGCTCTAGACGTTGCATCAGGTAGGGGGGGGCTTGTTTCATAACCTAGGGG\
TTCGTGAGGGTGGATTGGAATACTGACTATCACCTAGCAATGTTTTGCAGTATTATTTTAAGATCGTGTGCGGGGTGGGGAGCTATGTTTAAATAGCCTC\
GGCTCAACTTGCCAGCGTCCCATGAACTGTACAACCTAATCGTCCCTATCGTTACTGACCGACTGCGATAGAACGCCACCTTCATAGCGTCGCCGGACAA\
GCCTGTATGCAACCCATGAGTTTCCTTCGACTAGATCCAAACTCGAGGAGGTCATGGCGAGTCAAATTGTATATCTAGCGCCCACCTGATACCTAGGTTC\
".as_bytes().to_owned();

    let difference_model = SimpleAncientDnaModel::new(
        LibraryPrep::SingleStranded {
            five_prime_overhang: 0.475,
            three_prime_overhang: 0.475,
        },
        0.001,
        0.9,
        0.02 / 3.0,
        false,
    );

    let representative_mismatch_penalty = difference_model.get_representative_mismatch_penalty();

    let mismatch_bound =
        MismatchBoundDispatch::from(Discrete::new(0.02, 0.02, representative_mismatch_penalty));

    let parameters = AlignmentParameters {
        difference_model: difference_model.clone().into(),
        mismatch_bound,
        penalty_gap_open: 0.00001_f32.log2(),
        penalty_gap_extend: representative_mismatch_penalty,
        chunk_size: 1,
        gap_dist_ends: 5,
    };

    let alphabet = alphabets::Alphabet::new(mapad::index::DNA_UPPERCASE_ALPHABET);
    let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

    c.bench_function("bench_exogenous_read", |b| {
        let pattern = "GATATCTCGGCTGACAAACCAACAAAAAGTATCGGAACATCGCGGCGGCGTAGATGAATCTTAACCACACTCGACAGCTGTGCTTCTATACTAGCATTAC"
            .as_bytes()
            .to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack = MinMaxHeap::new();
        let mut tree = Tree::new();
        b.iter(|| {
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
            );
            assert_eq!(intervals.len(), 0);
        })
    });

    c.bench_function("bench_exogenous_read_deam", |b| {
        // 3' terminal C->T
        let pattern = "TTTATCTCGGCTGACAAACCAACAAAAAGTATCGGAACATCGCGGCGGCGTAGATGAATCTTAACCACACTCGACAGCTGTGCTTCTATACTAGCATTTT"
            .as_bytes()
            .to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack = MinMaxHeap::new();
        let mut tree = Tree::new();
        b.iter(|| {
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model
            );
            assert_eq!(intervals.len(), 0);
        })
    });

    c.bench_function("bench_endogenous_read_perfect", |b| {
        let pattern = "CCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTAGGGTGACCGTAAGATTTAGGTCCATCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC"
            .as_bytes()
            .to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack = MinMaxHeap::new();
        let mut tree = Tree::new();
        b.iter(|| {
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model
            );
            assert_eq!(intervals.len(), 1);
        })
    });

    c.bench_function("bench_endogenous_read_1_mm_center", |b| {
        let pattern = "CCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTAGGGTGAGCGTAAGATTTAGGTCCATCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC"
            .as_bytes()
            .to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack = MinMaxHeap::new();
        let mut tree = Tree::new();
        b.iter(|| {
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model
            );
            assert_eq!(intervals.len(), 1);
        })
    });

    c.bench_function("bench_endogenous_read_2_mm_center", |b| {
        let pattern = "CCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTTGGGTGACCGTAAGATTTAGGTCCGTCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC"
            .as_bytes()
            .to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack = MinMaxHeap::new();
        let mut tree = Tree::new();
        b.iter(|| {
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model
            );
            assert_eq!(intervals.len(), 1);
        })
    });

    c.bench_function("bench_endogenous_read_1_deam", |b| {
        let pattern = "TCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTTGGGTGACCGTAAGATTTAGGTCCGTCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC"
            .as_bytes()
            .to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack = MinMaxHeap::new();
        let mut tree = Tree::new();
        b.iter(|| {
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model
            );
            assert_eq!(intervals.len(), 1);
        })
    });

    c.bench_function("bench_endogenous_read_4_deam", |b| {
        let pattern = "TCTGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTTGGGTGACCGTAAGATTTAGGTCCGTCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGTT"
            .as_bytes()
            .to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack = MinMaxHeap::new();
        let mut tree = Tree::new();
        b.iter(|| {
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model
            );
            assert_eq!(intervals.len(), 1);
        })
    });
}

criterion_group!(benches, criterion_benchmark, bench_exo_endo,);
criterion_main!(benches);
