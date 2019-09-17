use bio::{
    alphabets,
    data_structures::{
        bwt::{bwt, less, Occ},
        fmindex::{FMDIndex, FMIndex},
        suffix_array::suffix_array,
    },
};

use criterion::{criterion_group, criterion_main, Criterion};

use mapad::{
    map::k_mismatch_search,
    sequence_difference_models::{LibraryPrep, SimpleAncientDnaModel},
    utils::{AlignmentParameters, AllowedMismatches},
};

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("3_mismatch_search", |b| {
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        let difference_model = SimpleAncientDnaModel {
            library_prep: (LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            }),
            ds_deamination_rate: 0.001,
            ss_deamination_rate: 0.9,
            divergence: 0.02 / 3.0,
        };

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -1.0,
            penalty_gap_extend: -1.0,
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

        let fm_index = FMIndex::new(&bwtr, &lessa, &occ);
        let fmd_index = FMDIndex::from(fm_index);

        let allowed_mismatches = AllowedMismatches::new(&parameters);

        let pattern = "GTTT".as_bytes().to_owned();
        let base_qualities = vec![40; pattern.len()];

        b.iter(|| {
            k_mismatch_search(
                &pattern,
                &base_qualities,
                -allowed_mismatches.get(pattern.len()),
                &parameters,
                &fmd_index,
            )
        })
    });
}

fn bench_multiple_reads(c: &mut Criterion) {
    c.bench_function("bench_multiple_reads", |b| {
        let mut ref_seq = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGA\
                           TTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTA\
                           CAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAG\
                           ATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT\
                           ACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA\
                           GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGAT\
                           TACAGATTACAGATTACAGATTACAGATTACA"
            .as_bytes()
            .to_owned();

        let difference_model = SimpleAncientDnaModel {
            library_prep: (LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            }),
            ds_deamination_rate: 0.001,
            ss_deamination_rate: 0.9,
            divergence: 0.02 / 3.0,
        };

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -1.0,
            penalty_gap_extend: -1.0,
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

        let fm_index = FMIndex::new(&bwtr, &lessa, &occ);
        let fmd_index = FMDIndex::from(fm_index);

        let allowed_mismatches = AllowedMismatches::new(&parameters);

        let patterns = vec![
            "TAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAG"
                .as_bytes()
                .to_owned();
            100
        ];
        let base_qualities = vec![40; patterns[0].len()];

        for pattern in patterns.iter() {
            b.iter(|| {
                k_mismatch_search(
                    &pattern,
                    &base_qualities,
                    -allowed_mismatches.get(pattern.len()),
                    &parameters,
                    &fmd_index,
                )
            })
        }
    });
}

fn bench_exogenous_reads(c: &mut Criterion) {
    c.bench_function("bench_exogenous_reads", |b| {
        let mut ref_seq = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGA\
                           TTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTA\
                           CAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAG\
                           ATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT\
                           ACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA\
                           GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGAT\
                           TACAGATTACAGATTACAGATTACAGATTACA"
            .as_bytes()
            .to_owned();

        let difference_model = SimpleAncientDnaModel {
            library_prep: (LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            }),
            ds_deamination_rate: 0.001,
            ss_deamination_rate: 0.9,
            divergence: 0.02 / 3.0,
        };

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -1.0,
            penalty_gap_extend: -1.0,
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

        let fm_index = FMIndex::new(&bwtr, &lessa, &occ);
        let fmd_index = FMDIndex::from(fm_index);

        let allowed_mismatches = AllowedMismatches::new(&parameters);

        let patterns = vec![
            "TTTTTTTTTTGGGGGTTACAGATTACAGATTACAGGGGGGTTTTTTTTTT"
                .as_bytes()
                .to_owned();
            100
        ];
        let base_qualities = vec![40; patterns[0].len()];

        for pattern in patterns.iter() {
            b.iter(|| {
                k_mismatch_search(
                    &pattern,
                    &base_qualities,
                    -allowed_mismatches.get(pattern.len()),
                    &parameters,
                    &fmd_index,
                )
            })
        }
    });
}

criterion_group!(
    benches,
    criterion_benchmark,
    bench_multiple_reads,
    bench_exogenous_reads
);
criterion_main!(benches);
