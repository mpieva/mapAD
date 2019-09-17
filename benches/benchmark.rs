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
    sequence_difference_models::SequenceDifferenceModel,
    utils::{AlignmentParameters, AllowedMismatches},
};

struct TestDifferenceModel {}
impl SequenceDifferenceModel for TestDifferenceModel {
    fn get(&self, _i: usize, _read_length: usize, from: u8, to: u8, _base_quality: u8) -> f32 {
        if from == b'C' && to == b'T' {
            return -0.5;
        } else if from != to {
            return -1.0;
        } else {
            return 0.0;
        }
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("3_mismatch_search", |b| {
        let mut ref_seq = "GATTACA".as_bytes().to_owned();

        let difference_model = TestDifferenceModel {};

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -1.0,
            penalty_gap_extend: -1.0,
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

        let difference_model = TestDifferenceModel {};

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -1.0,
            penalty_gap_extend: -1.0,
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

        let difference_model = TestDifferenceModel {};

        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            difference_model,
            penalty_gap_open: -1.0,
            penalty_gap_extend: -1.0,
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
