use bio::{
    alphabets,
    data_structures::{
        bwt::{bwt, less, Occ},
        fmindex::{FMDIndex, FMIndex},
        suffix_array::suffix_array,
    },
};
use criterion::{criterion_group, criterion_main, Criterion};
use min_max_heap::MinMaxHeap;

use backtrack_tree::Tree;
use mapad::{
    map::k_mismatch_search, mismatch_bounds::*, sequence_difference_models::*,
    utils::AlignmentParameters,
};

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("3_mismatch_search", |b| {
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        let difference_model = SequenceDifferenceModelDispatch::from(SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            },
            0.001,
            0.9,
            0.02 / 3.0,
        ));

        let representative_mismatch_penalty =
            difference_model.get_representative_mismatch_penalty();

        let mismatch_bound =
            MismatchBoundDispatch::from(Discrete::new(0.02, 0.02, representative_mismatch_penalty));

        let parameters = AlignmentParameters {
            difference_model,
            mismatch_bound,
            penalty_gap_open: 0.00001_f32.log2(),
            penalty_gap_extend: representative_mismatch_penalty,
            chunk_size: 1,
        };

        // Reference
        let ref_seq_rev_compl = alphabets::dna::revcomp(ref_seq.iter());
        ref_seq.extend_from_slice(b"$");
        ref_seq.extend_from_slice(&ref_seq_rev_compl);
        drop(ref_seq_rev_compl);
        ref_seq.extend_from_slice(b"$");

        let alphabet = alphabets::dna::alphabet();

        let sa = suffix_array(&ref_seq);
        let bwtr = bwt(&ref_seq, &sa);
        let lessa = less(&bwtr, &alphabet);
        let occ = Occ::new(&bwtr, 3, &alphabet);

        let fmd_index = FMDIndex::from(FMIndex::new(bwtr, lessa, occ));

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
            );
        })
    });
}

fn bench_mapping_read(c: &mut Criterion) {
    c.bench_function("bench_mapping_read", |b| {
        let mut ref_seq = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGA\
                           TTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTA\
                           CAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAG\
                           ATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT\
                           ACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA\
                           GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGAT\
                           TACAGATTACAGATTACAGATTACAGATTACA"
            .as_bytes()
            .to_owned();

        let difference_model = SequenceDifferenceModelDispatch::from(SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            },
            0.001,
            0.9,
            0.02 / 3.0,
        ));

        let representative_mismatch_penalty =
            difference_model.get_representative_mismatch_penalty();

        let mismatch_bound =
            MismatchBoundDispatch::from(Discrete::new(0.02, 0.02, representative_mismatch_penalty));

        let parameters = AlignmentParameters {
            difference_model,
            mismatch_bound,
            penalty_gap_open: 0.00001_f32.log2(),
            penalty_gap_extend: representative_mismatch_penalty,
            chunk_size: 1,
        };

        // Reference
        let ref_seq_rev_compl = alphabets::dna::revcomp(ref_seq.iter());
        ref_seq.extend_from_slice(b"$");
        ref_seq.extend_from_slice(&ref_seq_rev_compl);
        drop(ref_seq_rev_compl);
        ref_seq.extend_from_slice(b"$");

        let alphabet = alphabets::dna::alphabet();

        let sa = suffix_array(&ref_seq);
        let bwtr = bwt(&ref_seq, &sa);
        let lessa = less(&bwtr, &alphabet);
        let occ = Occ::new(&bwtr, 3, &alphabet);

        let fmd_index = FMDIndex::from(FMIndex::new(bwtr, lessa, occ));

        let pattern = "TAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAG"
            .as_bytes()
            .to_owned();
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
            );
        })
    });
}

fn bench_exogenous_read(c: &mut Criterion) {
    c.bench_function("bench_exogenous_read", |b| {
        let mut ref_seq = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGA\
                           TTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTA\
                           CAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAG\
                           ATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT\
                           ACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA\
                           GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGAT\
                           TACAGATTACAGATTACAGATTACAGATTACA"
            .as_bytes()
            .to_owned();

        let difference_model = SequenceDifferenceModelDispatch::from(SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            },
            0.001,
            0.9,
            0.02 / 3.0,
        ));

        let representative_mismatch_penalty =
            difference_model.get_representative_mismatch_penalty();

        let mismatch_bound =
            MismatchBoundDispatch::from(Discrete::new(0.02, 0.02, representative_mismatch_penalty));

        let parameters = AlignmentParameters {
            difference_model,
            mismatch_bound,
            penalty_gap_open: 0.00001_f32.log2(),
            penalty_gap_extend: representative_mismatch_penalty,
            chunk_size: 1,
        };

        // Reference
        let ref_seq_rev_compl = alphabets::dna::revcomp(ref_seq.iter());
        ref_seq.extend_from_slice(b"$");
        ref_seq.extend_from_slice(&ref_seq_rev_compl);
        drop(ref_seq_rev_compl);
        ref_seq.extend_from_slice(b"$");

        let alphabet = alphabets::dna::alphabet();

        let sa = suffix_array(&ref_seq);
        let bwtr = bwt(&ref_seq, &sa);
        let lessa = less(&bwtr, &alphabet);
        let occ = Occ::new(&bwtr, 3, &alphabet);

        let fmd_index = FMDIndex::from(FMIndex::new(bwtr, lessa, occ));

        let pattern = "TTTTTTTTTTGGGGGTTACAGATTACAGATTACAGGGGGGTTTTTTTTTT"
            .as_bytes()
            .to_owned();
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
            );
        })
    });
}

fn bench_long_read(c: &mut Criterion) {
    c.bench_function("bench_long_read", |b| {
        let mut ref_seq = "GTCTGCATCCCCAGGACCACCATGGGTGGGGAGGGCAGAGATTGGGGAGCACCTATAGAGGCTCTAATGCTCTAAGGTGACAGTGATGAGGACCTGGGTGCACCCATGAGTGGAGAAGCTAGGCCTGTCCAGAGAAGCAAGACAAACACACACATACACACTCACACACACACAGGCACATATGCATACACAAATACATTGCAT".as_bytes().to_owned();
        let pattern = "ACTAAAGAGCTTCTGCACAGGAAAAGAAACTACCATCAGAACCACCAGGCAACCTACAACATGGGATAAAATTTTCACAACCTACTCATCTGACAAAGGGCCAATATCCAGAATCTACAATGAACTCCAACAAATTTACAAGAAAAAAACAAACAACCCCATCAAAAAGTGGGCAAAGGACATGAACAGACACTTCTCAAAAGAAGATATTTATGCAGCCAAGAAAACACATAAAAA".as_bytes().to_owned();

        let difference_model = SequenceDifferenceModelDispatch::from(SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            },
            0.001,
            0.9,
            0.02 / 3.0,
        ));

        let representative_mismatch_penalty =
            difference_model.get_representative_mismatch_penalty();

        let mismatch_bound = MismatchBoundDispatch::from(Discrete::new(0.04, 0.02, representative_mismatch_penalty));

        let parameters = AlignmentParameters {
            difference_model,
            mismatch_bound,
            penalty_gap_open: 0.00001_f32.log2(),
            penalty_gap_extend: representative_mismatch_penalty,
            chunk_size: 1,
        };

        // Reference
        let ref_seq_rev_compl = alphabets::dna::revcomp(ref_seq.iter());
        ref_seq.extend_from_slice(b"$");
        ref_seq.extend_from_slice(&ref_seq_rev_compl);
        drop(ref_seq_rev_compl);
        ref_seq.extend_from_slice(b"$");

        let alphabet = alphabets::dna::alphabet();

        let sar = suffix_array(&ref_seq);
        let bwtr = bwt(&ref_seq, &sar);
        let lessa = less(&bwtr, &alphabet);
        let occ = Occ::new(&bwtr, 3, &alphabet);

        let fmd_index = FMDIndex::from(FMIndex::new(bwtr, lessa, occ));

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
                    );
            })
    });
}

criterion_group!(
    benches,
    criterion_benchmark,
    bench_mapping_read,
    bench_exogenous_read,
    bench_long_read,
);
criterion_main!(benches);
