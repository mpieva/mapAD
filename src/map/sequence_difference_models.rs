use std::fmt::{self, Display};

use either::Either;
use log::info;
use serde::{Deserialize, Serialize};

use crate::index::DNA_UPPERCASE_ALPHABET;

const MAX_ENCODED_BASE_QUALITY: u8 = u8::MAX;

/// New models must impl this trait and be included in the `SequenceDifferenceModelDispatch` enum in order to
/// be usable directly from functions calling `k_mismatch_search::<SDM>()`.
/// Sequence difference models are expected to yield only non-positive values (and 0.0), for example log probabilities.
pub trait SequenceDifferenceModel {
    fn get(&self, i: usize, read_length: usize, from: u8, to: u8, base_quality: u8) -> f32;
    fn get_representative_mismatch_penalty(&self) -> f32 {
        let read_length = 80;
        self.get(
            read_length / 2,
            read_length,
            b'T',
            b'A',
            MAX_ENCODED_BASE_QUALITY,
        ) - self.get(
            read_length / 2,
            read_length,
            b'T',
            b'T',
            MAX_ENCODED_BASE_QUALITY,
        )
    }

    /// Needed for the calculation of D arrays
    fn get_min_penalty(
        &self,
        i: usize,
        read_length: usize,
        to: u8,
        base_quality: u8,
        only_mismatches: bool,
    ) -> f32 {
        let iterator = DNA_UPPERCASE_ALPHABET.iter();

        if only_mismatches {
            Either::Left(iterator.filter(|&&base| base != to))
        } else {
            // Is the reference symbol ambiguous (not part of the DNA alphabet)? Then we don't want
            // to subtract an optimal score since this would make the ambiguous symbol case
            // penalty-free.
            if !DNA_UPPERCASE_ALPHABET.contains(&to) {
                return 0.0;
            }
            Either::Right(iterator)
        }
        .map(|&base| self.get(i, read_length, base, to, base_quality))
        .fold(f32::MIN, |acc, v| acc.max(v))
    }
}

/// Used to allow static dispatch. No trait objects needed! Method call speed is not negatively
/// affected by vtable lookups. Every type implementing `SequenceDifferenceModel` also has to be
/// added as variant to this enum.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum SequenceDifferenceModelDispatch {
    SimpleAncientDnaModel(SimpleAncientDnaModel),
    VindijaPwm(VindijaPwm),
    TestDifferenceModel(TestDifferenceModel),
}

impl From<SimpleAncientDnaModel> for SequenceDifferenceModelDispatch {
    fn from(sdm: SimpleAncientDnaModel) -> Self {
        Self::SimpleAncientDnaModel(sdm)
    }
}

impl From<TestDifferenceModel> for SequenceDifferenceModelDispatch {
    fn from(sdm: TestDifferenceModel) -> Self {
        Self::TestDifferenceModel(sdm)
    }
}

impl From<VindijaPwm> for SequenceDifferenceModelDispatch {
    fn from(sdm: VindijaPwm) -> Self {
        Self::VindijaPwm(sdm)
    }
}

/// Library preparation methods commonly used for ancient DNA. Values are overhang probabilities.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum LibraryPrep {
    SingleStranded {
        five_prime_overhang: f32,
        three_prime_overhang: f32,
    },
    DoubleStranded(f32),
}

/// Model of deamination (aDNA degradation), divergence, and sequencing error
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SimpleAncientDnaModel {
    library_prep: LibraryPrep,
    // Deamination rate in double-stranded stems
    ds_deamination_rate: f32,
    // Deamination rate in single-stranded overhangs
    ss_deamination_rate: f32,
    divergence: f32,
    use_default_base_quality: Option<f32>,
    cache: Vec<f32>,
}

impl SequenceDifferenceModel for SimpleAncientDnaModel {
    fn get(&self, i: usize, read_length: usize, from: u8, to: u8, base_quality: u8) -> f32 {
        let fp_dist = i;
        let tp_dist = read_length - 1 - i;

        // Only C->T and G->A substitutions are part of the deamination model.
        // Therefore, the following calculations are offloaded into a closure that
        // is only called when one of those substitutions is observed.
        let compute_deamination_part = || {
            let (p_fwd, p_rev) = match self.library_prep {
                LibraryPrep::SingleStranded {
                    five_prime_overhang,
                    three_prime_overhang,
                } => {
                    let five_prime_overhang = five_prime_overhang.powi(fp_dist as i32 + 1);
                    let three_prime_overhang = three_prime_overhang.powi(tp_dist as i32 + 1);
                    (
                        (five_prime_overhang + three_prime_overhang)
                            - (five_prime_overhang * three_prime_overhang),
                        0.0,
                    )
                }
                LibraryPrep::DoubleStranded(overhang) => (
                    overhang.powi(fp_dist as i32 + 1),
                    overhang.powi(tp_dist as i32 + 1),
                ),
            };

            // Probabilities of seeing C->T or G->A substitutions
            let c_to_t =
                self.ss_deamination_rate * p_fwd + self.ds_deamination_rate * (1.0 - p_fwd);
            let g_to_a =
                self.ss_deamination_rate * p_rev + self.ds_deamination_rate * (1.0 - p_rev);

            (c_to_t, g_to_a)
        };

        let sequencing_error = self.use_default_base_quality.unwrap_or_else(|| {
            // We need the specific value
            self.cache
                .get(base_quality as usize)
                .copied()
                // We need to actually compute the value
                .unwrap_or_else(|| Self::qual2prob(base_quality))
        });

        // Probability of seeing a mutation or sequencing error
        let independent_error =
            sequencing_error + self.divergence - sequencing_error * self.divergence;

        match from {
            b'A' => match to {
                b'A' => 1.0 - 3.0 * independent_error,
                _ => independent_error,
            },
            b'C' => match to {
                b'C' => {
                    let (c_to_t, _) = compute_deamination_part();
                    1.0 - 3.0 * independent_error - c_to_t + 4.0 * independent_error * c_to_t
                }
                b'T' => {
                    let (c_to_t, _) = compute_deamination_part();
                    independent_error + c_to_t - 4.0 * independent_error * c_to_t
                }
                _ => independent_error,
            },
            b'G' => match to {
                b'A' => {
                    let (_, g_to_a) = compute_deamination_part();
                    independent_error + g_to_a - 4.0 * independent_error * g_to_a
                }
                b'G' => {
                    let (_, g_to_a) = compute_deamination_part();
                    1.0 - 3.0 * independent_error - g_to_a + 4.0 * independent_error * g_to_a
                }
                _ => independent_error,
            },
            b'T' => match to {
                b'T' => 1.0 - 3.0 * independent_error,
                _ => independent_error,
            },
            _ => independent_error,
        }
        .max(f32::EPSILON)
        .log2()
    }
}

impl Display for SimpleAncientDnaModel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        const BASE_QUALITY: u8 = 37;
        const READ_LEN: usize = 50;

        // "Ordinary" mismatch penalty
        writeln!(
            f,
            "\"Ordinary\" mismatch: {:.2}",
            self.get_representative_mismatch_penalty()
        )?;

        // Central cheap mismatches
        writeln!(
            f,
            "Central C->T / G->A: {:.2}",
            self.get(READ_LEN / 2, READ_LEN, b'C', b'T', BASE_QUALITY)
        )?;

        // 5' end
        write!(f, "5' C->T: ")?;
        (0..10).try_for_each(|pos| {
            write!(
                f,
                "{:.2} ",
                self.get(pos, READ_LEN, b'C', b'T', BASE_QUALITY)
            )
        })?;
        write!(f, "...")?;

        // 3' end
        let three_prime_end = (READ_LEN - 10)..READ_LEN;
        match &self.library_prep {
            LibraryPrep::SingleStranded { .. } => {
                write!(f, "\n3' C->T: ")?;
                three_prime_end.rev().try_for_each(|pos| {
                    write!(
                        f,
                        "{:.2} ",
                        self.get(pos, READ_LEN, b'C', b'T', BASE_QUALITY)
                    )
                })?;
            }
            LibraryPrep::DoubleStranded(_) => {
                write!(f, "\n3' G->A: ")?;
                three_prime_end.rev().try_for_each(|pos| {
                    write!(
                        f,
                        "{:.2} ",
                        self.get(pos, READ_LEN, b'G', b'A', BASE_QUALITY)
                    )
                })?;
            }
        }
        write!(f, "...")?;

        Ok(())
    }
}

impl SimpleAncientDnaModel {
    fn qual2prob(encoded_base_quality: u8) -> f32 {
        10_f32.powf(-1.0 * encoded_base_quality as f32 / 10.0) / 3.0
    }

    pub fn new(
        library_prep: LibraryPrep,
        ds_deamination_rate: f32,
        ss_deamination_rate: f32,
        divergence: f32,
        ignore_base_qualities: bool,
    ) -> Self {
        let use_default_base_quality =
            ignore_base_qualities.then_some(Self::qual2prob(MAX_ENCODED_BASE_QUALITY));

        let cache = if use_default_base_quality.is_some() {
            // No cache needed
            Vec::new()
        } else {
            (0..=MAX_ENCODED_BASE_QUALITY)
                .map(Self::qual2prob)
                .collect::<Vec<_>>()
        };

        let out = Self {
            library_prep,
            ds_deamination_rate,
            ss_deamination_rate,
            divergence,
            use_default_base_quality,
            cache,
        };

        // Requiring `Display` as a supertrait to `SequenceDifferenceModel` would be a better
        // solution, but due to the static dispatch tricks it's not that simple.
        info!("{}", out);

        out
    }
}

/// Very simple model of ancient DNA degradation for starters.
/// It only takes C->T deaminations into account and assumes
/// symmetry between 5' and 3' ends of the deamination pattern.
#[derive(Default, Serialize, Deserialize, Debug, Clone)]
pub struct VindijaPwm {
    // "Sparse" position probability matrix
    ppm_read_ends_symmetric_ct: [f32; 7],
    position_probability_ct_default: f32,
    observed_substitution_probability_default: f32,
}

impl Display for VindijaPwm {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "This is only for testing purposes.")
    }
}

impl SequenceDifferenceModel for VindijaPwm {
    fn get(&self, i: usize, read_length: usize, from: u8, to: u8, _base_quality: u8) -> f32 {
        let position_probability = match from {
            b'C' => {
                let i = i.min(read_length - (i + 1));
                let position_probability_ct = *self
                    .ppm_read_ends_symmetric_ct
                    .get(i)
                    .unwrap_or(&self.position_probability_ct_default);
                match to {
                    // C->T or C->C
                    b'T' => position_probability_ct,
                    b'C' => 1.0 - position_probability_ct,
                    // "Normal" mismatch
                    _ => self.observed_substitution_probability_default,
                }
            }
            _ => {
                if from == to {
                    // "Normal" match
                    1.0 - self.observed_substitution_probability_default
                } else {
                    // "Normal" mismatch
                    self.observed_substitution_probability_default
                }
            }
        };
        position_probability.log2()
    }
}

impl VindijaPwm {
    pub fn new() -> Self {
        VindijaPwm {
            // The following values are roughly derived with
            // the naked eye from Prüfer et al. (2017), Fig. S3.
            ppm_read_ends_symmetric_ct: [0.4, 0.25, 0.1, 0.06, 0.05, 0.04, 0.03],
            position_probability_ct_default: 0.02,
            observed_substitution_probability_default: 0.0005,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestDifferenceModel {
    pub deam_score: f32,
    pub mm_score: f32,
    pub match_score: f32,
}

impl Display for TestDifferenceModel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "This is only for testing purposes.")
    }
}

impl SequenceDifferenceModel for TestDifferenceModel {
    fn get(&self, _i: usize, _read_len: usize, from: u8, to: u8, _base_quality: u8) -> f32 {
        if from == b'C' && to == b'T' {
            self.deam_score
        } else if from == to {
            self.match_score
        } else {
            self.mm_score
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_vindija_pwm() {
        let vindija_pwm = VindijaPwm::new();
        let read_length = 35;

        //        println!("Vindija PWM");
        //        for i in 0..read_length {
        //            println!(
        //                "{i})\tC->T: {c_t}\t\tC->C: {c_c}\t\tA->A: {a_a}\t\t G->A: {g_a}",
        //                i = i,
        //                c_t = vindija_pwm.get(i, read_length, b'C', b'T', 40),
        //                c_c = vindija_pwm.get(i, read_length, b'C', b'C', 40),
        //                a_a = vindija_pwm.get(i, read_length, b'A', b'A', 40),
        //                g_a = vindija_pwm.get(i, read_length, b'G', b'A', 40),
        //            );
        //        }

        assert_approx_eq!(-1.321928, vindija_pwm.get(0, read_length, b'C', b'T', 40));
        assert_approx_eq!(-0.736965, vindija_pwm.get(0, read_length, b'C', b'C', 40));
        assert_approx_eq!(-5.643856, vindija_pwm.get(15, read_length, b'C', b'T', 40));
        assert_approx_eq!(-10.965784, vindija_pwm.get(15, read_length, b'G', b'C', 40));
        assert_approx_eq!(-0.000721, vindija_pwm.get(15, read_length, b'A', b'A', 40));
    }

    #[test]
    fn test_simple_adna_model() {
        // Model heavily damaged data
        let adna_model = SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.6,
                three_prime_overhang: 0.55,
            },
            0.01,
            1.0,
            0.02 / 3.0,
            false,
        );

        assert_approx_eq!(-0.183332, adna_model.get(0, 25, b'A', b'A', 10));
        assert_approx_eq!(-0.029293, adna_model.get(1, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.183332, adna_model.get(2, 25, b'A', b'A', 10));
        assert_approx_eq!(-0.029293, adna_model.get(3, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(4, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.183332, adna_model.get(5, 25, b'A', b'A', 10));
        assert_approx_eq!(-0.029293, adna_model.get(6, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(7, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(8, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(9, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.183332, adna_model.get(10, 25, b'A', b'A', 10));
        assert_approx_eq!(-0.029293, adna_model.get(11, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.183332, adna_model.get(12, 25, b'A', b'A', 10));
        assert_approx_eq!(-0.029293, adna_model.get(13, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(14, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(15, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(16, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(17, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(18, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(19, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(20, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(21, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(22, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(23, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.183332, adna_model.get(24, 25, b'A', b'A', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'A', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'A', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'A', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'A', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'A', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'A', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'A', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'A', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'A', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'A', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'A', b'C', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'A', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'A', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'A', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'A', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'A', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'A', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'A', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'A', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'A', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'A', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'A', b'G', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'A', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'A', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'A', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'A', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'A', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'A', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'A', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'A', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'A', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'A', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'A', b'T', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'C', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'C', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'C', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'C', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'C', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'C', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'C', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'C', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'C', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'C', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'C', b'A', 10));
        assert_approx_eq!(-1.423644, adna_model.get(0, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.681956, adna_model.get(1, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.530236, adna_model.get(2, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.242462, adna_model.get(3, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.159644, adna_model.get(4, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.262897, adna_model.get(5, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.084385, adna_model.get(6, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.067990, adna_model.get(7, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.058259, adna_model.get(8, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.052482, adna_model.get(9, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.202353, adna_model.get(10, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.047147, adna_model.get(11, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.199553, adna_model.get(12, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.045914, adna_model.get(13, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.046364, adna_model.get(14, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.047730, adna_model.get(15, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.050548, adna_model.get(16, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.055885, adna_model.get(17, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.065759, adna_model.get(18, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.083959, adna_model.get(19, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.117695, adna_model.get(20, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.181160, adna_model.get(21, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.304246, adna_model.get(22, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.559120, adna_model.get(23, 25, b'C', b'C', 40));
        assert_approx_eq!(-1.270929, adna_model.get(24, 25, b'C', b'C', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'C', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'C', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'C', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'C', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'C', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'C', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'C', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'C', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'C', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'C', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'C', b'G', 10));
        assert_approx_eq!(-0.868609, adna_model.get(0, 25, b'C', b'T', 10));
        assert_approx_eq!(-1.460842, adna_model.get(1, 25, b'C', b'T', 40));
        assert_approx_eq!(-2.132875, adna_model.get(2, 25, b'C', b'T', 10));
        assert_approx_eq!(-2.823177, adna_model.get(3, 25, b'C', b'T', 40));
        assert_approx_eq!(-3.452384, adna_model.get(4, 25, b'C', b'T', 40));
        assert_approx_eq!(-3.522311, adna_model.get(5, 25, b'C', b'T', 10));
        assert_approx_eq!(-4.525707, adna_model.get(6, 25, b'C', b'T', 40));
        assert_approx_eq!(-4.937460, adna_model.get(7, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.255495, adna_model.get(8, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.485218, adna_model.get(9, 25, b'C', b'T', 40));
        assert_approx_eq!(-4.284543, adna_model.get(10, 25, b'C', b'T', 10));
        assert_approx_eq!(-5.736821, adna_model.get(11, 25, b'C', b'T', 40));
        assert_approx_eq!(-4.332809, adna_model.get(12, 25, b'C', b'T', 10));
        assert_approx_eq!(-5.801927, adna_model.get(13, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.777827, adna_model.get(14, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.707015, adna_model.get(15, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.571322, adna_model.get(16, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.345414, adna_model.get(17, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.004263, adna_model.get(18, 25, b'C', b'T', 40));
        assert_approx_eq!(-4.534981, adna_model.get(19, 25, b'C', b'T', 40));
        assert_approx_eq!(-3.944710, adna_model.get(20, 25, b'C', b'T', 40));
        assert_approx_eq!(-3.256952, adna_model.get(21, 25, b'C', b'T', 40));
        assert_approx_eq!(-2.500338, adna_model.get(22, 25, b'C', b'T', 40));
        assert_approx_eq!(-1.699540, adna_model.get(23, 25, b'C', b'T', 40));
        assert_approx_eq!(-0.982643, adna_model.get(24, 25, b'C', b'T', 10));
        assert_approx_eq!(-4.375222, adna_model.get(0, 25, b'G', b'A', 10));
        assert_approx_eq!(-5.927367, adna_model.get(1, 25, b'G', b'A', 40));
        assert_approx_eq!(-4.375222, adna_model.get(2, 25, b'G', b'A', 10));
        assert_approx_eq!(-5.927367, adna_model.get(3, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(4, 25, b'G', b'A', 40));
        assert_approx_eq!(-4.375222, adna_model.get(5, 25, b'G', b'A', 10));
        assert_approx_eq!(-5.927367, adna_model.get(6, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(7, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(8, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(9, 25, b'G', b'A', 40));
        assert_approx_eq!(-4.375222, adna_model.get(10, 25, b'G', b'A', 10));
        assert_approx_eq!(-5.927367, adna_model.get(11, 25, b'G', b'A', 40));
        assert_approx_eq!(-4.375222, adna_model.get(12, 25, b'G', b'A', 10));
        assert_approx_eq!(-5.927367, adna_model.get(13, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(14, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(15, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(16, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(17, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(18, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(19, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(20, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(21, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(22, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927367, adna_model.get(23, 25, b'G', b'A', 40));
        assert_approx_eq!(-4.375222, adna_model.get(24, 25, b'G', b'A', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'G', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'G', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'G', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'G', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'G', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'G', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'G', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'G', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'G', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'G', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'G', b'C', 10));
        assert_approx_eq!(-0.197174, adna_model.get(0, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.043693, adna_model.get(1, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.197174, adna_model.get(2, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.043693, adna_model.get(3, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(4, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.197174, adna_model.get(5, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.043693, adna_model.get(6, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(7, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(8, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(9, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.197174, adna_model.get(10, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.043693, adna_model.get(11, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.197174, adna_model.get(12, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.043693, adna_model.get(13, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(14, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(15, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(16, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(17, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(18, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(19, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(20, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(21, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(22, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(23, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.197174, adna_model.get(24, 25, b'G', b'G', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'G', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'G', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'G', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'G', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'G', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'G', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'G', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'G', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'G', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'G', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'G', b'T', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'T', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'T', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'T', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'T', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'T', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'T', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'T', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'T', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'T', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'T', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'T', b'A', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'T', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'T', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'T', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'T', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'T', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'T', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'T', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'T', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'T', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'T', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'T', b'C', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'T', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'T', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'T', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'T', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'T', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'T', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'T', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'T', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'T', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'T', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'T', b'G', 10));
        assert_approx_eq!(-0.183332, adna_model.get(0, 25, b'T', b'T', 10));
        assert_approx_eq!(-0.029293, adna_model.get(1, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.183332, adna_model.get(2, 25, b'T', b'T', 10));
        assert_approx_eq!(-0.029293, adna_model.get(3, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(4, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.183332, adna_model.get(5, 25, b'T', b'T', 10));
        assert_approx_eq!(-0.029293, adna_model.get(6, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(7, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(8, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(9, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.183332, adna_model.get(10, 25, b'T', b'T', 10));
        assert_approx_eq!(-0.029293, adna_model.get(11, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.183332, adna_model.get(12, 25, b'T', b'T', 10));
        assert_approx_eq!(-0.029293, adna_model.get(13, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(14, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(15, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(16, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(17, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(18, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(19, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(20, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(21, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(22, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(23, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.183332, adna_model.get(24, 25, b'T', b'T', 10));
    }

    #[test]
    fn test_simple_adna_model_ds() {
        let adna_model = SimpleAncientDnaModel::new(
            LibraryPrep::DoubleStranded(0.475),
            0.01,
            0.9,
            0.02 / 3.0,
            false,
        );

        assert_approx_eq!(-0.183332, adna_model.get(0, 25, b'A', b'A', 10));
        assert_approx_eq!(-0.029293, adna_model.get(1, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.183332, adna_model.get(2, 25, b'A', b'A', 10));
        assert_approx_eq!(-0.029293, adna_model.get(3, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(4, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.183332, adna_model.get(5, 25, b'A', b'A', 10));
        assert_approx_eq!(-0.029293, adna_model.get(6, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(7, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(8, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(9, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.183332, adna_model.get(10, 25, b'A', b'A', 10));
        assert_approx_eq!(-0.029293, adna_model.get(11, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.183332, adna_model.get(12, 25, b'A', b'A', 10));
        assert_approx_eq!(-0.029293, adna_model.get(13, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(14, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(15, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(16, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(17, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(18, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(19, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(20, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(21, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(22, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029293, adna_model.get(23, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.183332, adna_model.get(24, 25, b'A', b'A', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'A', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'A', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'A', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'A', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'A', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'A', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'A', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'A', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'A', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'A', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'A', b'C', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'A', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'A', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'A', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'A', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'A', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'A', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'A', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'A', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'A', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'A', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'A', b'G', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'A', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'A', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'A', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'A', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'A', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'A', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'A', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'A', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'A', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'A', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'A', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'A', b'T', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'C', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'C', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'C', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'C', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'C', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'C', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'C', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'C', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'C', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'C', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'C', b'A', 10));
        assert_approx_eq!(-0.952400, adna_model.get(0, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.368209, adna_model.get(1, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.336334, adna_model.get(2, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.110798, adna_model.get(3, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.075179, adna_model.get(4, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.211461, adna_model.get(5, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.050737, adna_model.get(6, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.047034, adna_model.get(7, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.045279, adna_model.get(8, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.044446, adna_model.get(9, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.197517, adna_model.get(10, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.043862, adna_model.get(11, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.197251, adna_model.get(12, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.043731, adna_model.get(13, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.043711, adna_model.get(14, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.043701, adna_model.get(15, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.043697, adna_model.get(16, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.043694, adna_model.get(17, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.043693, adna_model.get(18, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.043693, adna_model.get(19, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.043693, adna_model.get(20, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.043693, adna_model.get(21, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.043693, adna_model.get(22, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.043693, adna_model.get(23, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.197174, adna_model.get(24, 25, b'C', b'C', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'C', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'C', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'C', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'C', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'C', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'C', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'C', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'C', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'C', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'C', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'C', b'G', 10));
        assert_approx_eq!(-1.308743, adna_model.get(0, 25, b'C', b'T', 10));
        assert_approx_eq!(-2.238840, adna_model.get(1, 25, b'C', b'T', 40));
        assert_approx_eq!(-2.961360, adna_model.get(2, 25, b'C', b'T', 10));
        assert_approx_eq!(-4.046337, adna_model.get(3, 25, b'C', b'T', 40));
        assert_approx_eq!(-4.741751, adna_model.get(4, 25, b'C', b'T', 40));
        assert_approx_eq!(-4.138409, adna_model.get(5, 25, b'C', b'T', 10));
        assert_approx_eq!(-5.562702, adna_model.get(6, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.742640, adna_model.get(7, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.836668, adna_model.get(8, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.883573, adna_model.get(9, 25, b'C', b'T', 40));
        assert_approx_eq!(-4.369012, adna_model.get(10, 25, b'C', b'T', 10));
        assert_approx_eq!(-5.917369, adna_model.get(11, 25, b'C', b'T', 40));
        assert_approx_eq!(-4.373819, adna_model.get(12, 25, b'C', b'T', 10));
        assert_approx_eq!(-5.925105, adna_model.get(13, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.926292, adna_model.get(14, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.926856, adna_model.get(15, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.927124, adna_model.get(16, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.927252, adna_model.get(17, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.927312, adna_model.get(18, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.927341, adna_model.get(19, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.927354, adna_model.get(20, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.927361, adna_model.get(21, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.927364, adna_model.get(22, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.927366, adna_model.get(23, 25, b'C', b'T', 40));
        assert_approx_eq!(-4.375222, adna_model.get(24, 25, b'C', b'T', 10));
        assert_approx_eq!(-4.375222, adna_model.get(0, 25, b'G', b'A', 10));
        assert_approx_eq!(-5.927366, adna_model.get(1, 25, b'G', b'A', 40));
        assert_approx_eq!(-4.375221, adna_model.get(2, 25, b'G', b'A', 10));
        assert_approx_eq!(-5.927361, adna_model.get(3, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927354, adna_model.get(4, 25, b'G', b'A', 40));
        assert_approx_eq!(-4.375215, adna_model.get(5, 25, b'G', b'A', 10));
        assert_approx_eq!(-5.927312, adna_model.get(6, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927252, adna_model.get(7, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.927124, adna_model.get(8, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.926856, adna_model.get(9, 25, b'G', b'A', 40));
        assert_approx_eq!(-4.374905, adna_model.get(10, 25, b'G', b'A', 10));
        assert_approx_eq!(-5.925105, adna_model.get(11, 25, b'G', b'A', 40));
        assert_approx_eq!(-4.373819, adna_model.get(12, 25, b'G', b'A', 10));
        assert_approx_eq!(-5.917369, adna_model.get(13, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.906399, adna_model.get(14, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.883573, adna_model.get(15, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.836668, adna_model.get(16, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.742640, adna_model.get(17, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.562702, adna_model.get(18, 25, b'G', b'A', 40));
        assert_approx_eq!(-5.244400, adna_model.get(19, 25, b'G', b'A', 40));
        assert_approx_eq!(-4.741751, adna_model.get(20, 25, b'G', b'A', 40));
        assert_approx_eq!(-4.046337, adna_model.get(21, 25, b'G', b'A', 40));
        assert_approx_eq!(-3.194182, adna_model.get(22, 25, b'G', b'A', 40));
        assert_approx_eq!(-2.238840, adna_model.get(23, 25, b'G', b'A', 40));
        assert_approx_eq!(-1.308743, adna_model.get(24, 25, b'G', b'A', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'G', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'G', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'G', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'G', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'G', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'G', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'G', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'G', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'G', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'G', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'G', b'C', 10));
        assert_approx_eq!(-0.197174, adna_model.get(0, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.043693, adna_model.get(1, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.197174, adna_model.get(2, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.043693, adna_model.get(3, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043693, adna_model.get(4, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.197174, adna_model.get(5, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.043693, adna_model.get(6, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043694, adna_model.get(7, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043697, adna_model.get(8, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.043701, adna_model.get(9, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.197191, adna_model.get(10, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.043731, adna_model.get(11, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.197251, adna_model.get(12, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.043862, adna_model.get(13, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.044050, adna_model.get(14, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.044446, adna_model.get(15, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.045279, adna_model.get(16, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.047034, adna_model.get(17, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.050737, adna_model.get(18, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.058563, adna_model.get(19, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.075179, adna_model.get(20, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.110798, adna_model.get(21, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.188789, adna_model.get(22, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.368209, adna_model.get(23, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.952400, adna_model.get(24, 25, b'G', b'G', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'G', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'G', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'G', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'G', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'G', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'G', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'G', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'G', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'G', b'T', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'G', b'T', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'G', b'T', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'T', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'T', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'T', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'T', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'T', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'T', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'T', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'T', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'T', b'A', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'T', b'A', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'T', b'A', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'T', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'T', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'T', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'T', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'T', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'T', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'T', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'T', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'T', b'C', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'T', b'C', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'T', b'C', 10));
        assert_approx_eq!(-4.651894, adna_model.get(0, 25, b'T', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(1, 25, b'T', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(2, 25, b'T', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(3, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(4, 25, b'T', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(5, 25, b'T', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(6, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(7, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(8, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(9, 25, b'T', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(10, 25, b'T', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(11, 25, b'T', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(12, 25, b'T', b'G', 10));
        assert_approx_eq!(-7.221671, adna_model.get(13, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(14, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(15, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(16, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(17, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(18, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(19, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(20, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(21, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(22, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.221671, adna_model.get(23, 25, b'T', b'G', 40));
        assert_approx_eq!(-4.651894, adna_model.get(24, 25, b'T', b'G', 10));
        assert_approx_eq!(-0.183332, adna_model.get(0, 25, b'T', b'T', 10));
        assert_approx_eq!(-0.029293, adna_model.get(1, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.183332, adna_model.get(2, 25, b'T', b'T', 10));
        assert_approx_eq!(-0.029293, adna_model.get(3, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(4, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.183332, adna_model.get(5, 25, b'T', b'T', 10));
        assert_approx_eq!(-0.029293, adna_model.get(6, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(7, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(8, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(9, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.183332, adna_model.get(10, 25, b'T', b'T', 10));
        assert_approx_eq!(-0.029293, adna_model.get(11, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.183332, adna_model.get(12, 25, b'T', b'T', 10));
        assert_approx_eq!(-0.029293, adna_model.get(13, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(14, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(15, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(16, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(17, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(18, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(19, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(20, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(21, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(22, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.029293, adna_model.get(23, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.183332, adna_model.get(24, 25, b'T', b'T', 10));
    }

    #[test]
    fn test_simple_adna_wo_deam() {
        let adna_model = SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.0,
                three_prime_overhang: 0.0,
            },
            0.0,
            0.0,
            0.02 / 3.0,
            false,
        );

        assert_eq!(
            adna_model.get(0, 25, b'C', b'T', 40),
            adna_model.get(13, 25, b'T', b'A', 40)
        );
        assert_eq!(
            adna_model.get(24, 25, b'C', b'T', 40),
            adna_model.get(13, 25, b'T', b'A', 40)
        );
        assert_eq!(
            adna_model.get(13, 25, b'C', b'C', 40),
            adna_model.get(0, 25, b'C', b'C', 40)
        );
    }

    #[test]
    fn display_simple_adna_model() {
        let adna_model_sist = SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.4,
                three_prime_overhang: 0.3,
            },
            0.02,
            1.0,
            0.02 / 3.0,
            false,
        );

        let sist_display = "\"Ordinary\" mismatch: -7.20\n\
Central C->T / G->A: -5.25\n\
5' C->T: -1.29 -2.48 -3.52 -4.30 -4.80 -5.05 -5.17 -5.22 -5.24 -5.25 ...\n\
3' C->T: -1.68 -3.16 -4.27 -4.88 -5.13 -5.22 -5.24 -5.25 -5.25 -5.25 ...";

        assert_eq!(sist_display, format!("{}", adna_model_sist));

        let adna_model_ds = SimpleAncientDnaModel::new(
            LibraryPrep::DoubleStranded(0.4),
            0.02,
            1.0,
            0.02 / 3.0,
            false,
        );

        let ds_display = "\"Ordinary\" mismatch: -7.20\n\
Central C->T / G->A: -5.25\n\
5' C->T: -1.29 -2.48 -3.52 -4.30 -4.80 -5.05 -5.17 -5.22 -5.24 -5.25 ...\n\
3' G->A: -1.29 -2.48 -3.52 -4.30 -4.80 -5.05 -5.17 -5.22 -5.24 -5.25 ...";

        assert_eq!(ds_display, format!("{}", adna_model_ds));
    }
}
