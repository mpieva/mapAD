use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMDIndex, FMIndex};
use bio::data_structures::suffix_array::suffix_array;

use criterion::{criterion_group, criterion_main, Criterion};

use thrust::sequence_difference_models::SequenceDifferenceModel;
use thrust::map::k_mismatch_search;
use thrust::utils::{AlignmentParameters, AllowedMismatches};

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("3_mismatch_search", |b| {
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            penalty_gap_open: -1.0,
            penalty_gap_extend: -1.0,
        };

        struct TestDifferenceModel {}
        impl SequenceDifferenceModel for TestDifferenceModel {
            fn new_default() -> Self {
                TestDifferenceModel {}
            }
            fn get(&self, _i: usize, _read_length: usize, from: u8, to: u8) -> f32 {
                if from == b'C' && to == b'T' {
                    return 0.0;
                } else if from != to {
                    return -1.0;
                } else {
                    return 1.0;
                }
            }
        }
        let difference_model = TestDifferenceModel::new_default();

        // Reference
        let ref_seq_rev_compl = alphabets::dna::revcomp(ref_seq.iter());
        ref_seq.extend_from_slice(b"$");
        ref_seq.extend_from_slice(&ref_seq_rev_compl);
        drop(ref_seq_rev_compl);
        ref_seq.extend_from_slice(b"$");

        let alphabet = alphabets::dna::n_alphabet();

        let sa = suffix_array(&ref_seq);
        let bwtr = bwt(&ref_seq, &sa);
        let lessa = less(&bwtr, &alphabet);
        let occ = Occ::new(&bwtr, 3, &alphabet);

        let fm_index = FMIndex::new(&bwtr, &lessa, &occ);
        let fmd_index = FMDIndex::from(fm_index);

        // Reverse reference
        let mut rev_ref_seq = ref_seq.into_iter().rev().collect::<Vec<_>>();
        let rev_ref_seq_rev_compl = alphabets::dna::revcomp(rev_ref_seq.iter());
        rev_ref_seq.extend_from_slice(b"$");
        rev_ref_seq.extend_from_slice(&rev_ref_seq_rev_compl);
        drop(rev_ref_seq_rev_compl);
        rev_ref_seq.extend_from_slice(b"$");

        let rev_sa = suffix_array(&rev_ref_seq);
        let rev_bwtr = bwt(&rev_ref_seq, &rev_sa);
        let rev_less = less(&rev_bwtr, &alphabet);
        let rev_occ = Occ::new(&rev_bwtr, 3, &alphabet);

        let rev_fm_index = FMIndex::new(&rev_bwtr, &rev_less, &rev_occ);
        let rev_fmd_index = FMDIndex::from(rev_fm_index);

        let mut allowed_mismatches = AllowedMismatches::new(&parameters);

        let pattern = "GTTT".as_bytes().to_owned();
        let base_qualities = [0, 0, 0, 0];

        b.iter(|| {
            k_mismatch_search(
                &pattern,
                &base_qualities,
                allowed_mismatches.get(pattern.len()),
                &parameters,
                &difference_model,
                &fmd_index,
                &rev_fmd_index,
            )
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
