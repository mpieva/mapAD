use either::Either;

// TODO: Use lookup tables instead of computing scores on-demand

/// Sequence difference models are expected to yield only non-positive values (and 0.0), for example log probabilities.
pub trait SequenceDifferenceModel {
    fn get(&self, i: usize, read_length: usize, from: u8, to: u8, base_quality: u8) -> f32;
    fn get_representative_mismatch_penalty(&self) -> f32 {
        self.get(40, 80, b'T', b'A', 40)
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
        let iterator = b"ACGT".iter();

        if only_mismatches {
            Either::Left(iterator.filter(|&&base| base != to))
        } else {
            Either::Right(iterator)
        }
        .map(|&base| self.get(i, read_length, base, to, base_quality))
        .filter(|&penalty| penalty < 0.0)
        .fold(std::f32::MIN, |acc, v| acc.max(v))
    }
}

/// Library preparation methods commonly used for ancient DNA. Values are overhang probabilities.
pub enum LibraryPrep {
    SingleStranded {
        five_prime_overhang: f32,
        three_prime_overhang: f32,
    },
    DoubleStranded(f32),
}

/// Model of deamination (aDNA degradation), divergence, and sequencing error
pub struct SimpleAncientDnaModel {
    pub library_prep: LibraryPrep,
    // Deamination rate in double-stranded stems
    pub ds_deamination_rate: f32,
    // Deamination rate in single-stranded overhangs
    pub ss_deamination_rate: f32,
    pub divergence: f32,
}

impl SequenceDifferenceModel for SimpleAncientDnaModel {
    fn get(&self, i: usize, read_length: usize, from: u8, to: u8, base_quality: u8) -> f32 {
        let (p_fwd, p_rev) = match self.library_prep {
            LibraryPrep::SingleStranded {
                five_prime_overhang,
                three_prime_overhang,
            } => {
                let five_prime_overhang = five_prime_overhang.powi((i + 1) as i32);
                let three_prime_overhang = three_prime_overhang.powi((read_length - i) as i32);
                (
                    (five_prime_overhang + three_prime_overhang)
                        - (five_prime_overhang * three_prime_overhang),
                    0.0,
                )
            }
            LibraryPrep::DoubleStranded(overhang) => (
                overhang.powi(i as i32 + 1),
                overhang.powi((read_length - i) as i32),
            ),
        };

        let sequencing_error = 10_f32.powf(-1.0 * f32::from(base_quality) / 10.0);

        // Probability of seeing a mutation or sequencing error
        // Artificial boundary at 0.25 to ensure the model's output won't be negative
        let independent_error =
            (sequencing_error + self.divergence - sequencing_error * self.divergence).min(0.25);

        // Probabilities of seeing C->T or G->A substitutions
        let c_to_t = self.ss_deamination_rate * p_fwd + self.ds_deamination_rate * (1.0 - p_fwd);
        let g_to_a = self.ss_deamination_rate * p_rev + self.ds_deamination_rate * (1.0 - p_rev);

        match from {
            b'A' => match to {
                b'A' => 1.0 - 3.0 * independent_error,
                _ => independent_error,
            },
            b'C' => match to {
                b'C' => 1.0 - 3.0 * independent_error - c_to_t + 4.0 * independent_error * c_to_t,
                b'T' => independent_error + c_to_t - 4.0 * independent_error * c_to_t,
                _ => independent_error,
            },
            b'G' => match to {
                b'A' => independent_error + g_to_a - 4.0 * independent_error * g_to_a,
                b'G' => 1.0 - 3.0 * independent_error - g_to_a + 4.0 * independent_error * g_to_a,
                _ => independent_error,
            },
            b'T' => match to {
                b'T' => 1.0 - 3.0 * independent_error,
                _ => independent_error,
            },
            _ => independent_error,
        }
        .log2()
    }
}

/// Very simple model of ancient DNA degradation for starters.
/// It only takes C->T deaminations into account and assumes
/// symmetry between 5' and 3' ends of the deamination pattern.
#[derive(Default)]
pub struct VindijaPWM {
    // "Sparse" position probability matrix
    ppm_read_ends_symmetric_ct: [f32; 7],
    position_probability_ct_default: f32,
    observed_substitution_probability_default: f32,
}

impl VindijaPWM {
    pub fn new() -> Self {
        VindijaPWM {
            // The following values are roughly derived with
            // the naked eye from Prüfer et al. (2017), Fig. S3.
            ppm_read_ends_symmetric_ct: [0.4, 0.25, 0.1, 0.06, 0.05, 0.04, 0.03],
            position_probability_ct_default: 0.02,
            observed_substitution_probability_default: 0.0005,
        }
    }
}

impl SequenceDifferenceModel for VindijaPWM {
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

#[cfg(test)]
mod tests {
    use super::*;

    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_vindija_pwm() {
        let vindija_pwm = VindijaPWM::new();
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
        let adna_model = SimpleAncientDnaModel {
            library_prep: (LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            }),
            ds_deamination_rate: 0.001,
            ss_deamination_rate: 0.9,
            divergence: 0.02 / 3.0,
        };

        // 'C' -> 'T'
        assert_approx_eq!(-1.504131, adna_model.get(0, 25, b'C', b'T', 10));
        assert_approx_eq!(-2.285697, adna_model.get(1, 25, b'C', b'T', 40));
        assert_approx_eq!(-2.625292, adna_model.get(2, 25, b'C', b'T', 10));
        assert_approx_eq!(-4.257998, adna_model.get(3, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.113328, adna_model.get(4, 25, b'C', b'T', 40));
        assert_approx_eq!(-3.151698, adna_model.get(5, 25, b'C', b'T', 10));
        assert_approx_eq!(-6.320596, adna_model.get(6, 25, b'C', b'T', 40));
        assert_approx_eq!(-6.642854, adna_model.get(7, 25, b'C', b'T', 40));
        assert_approx_eq!(-6.825267, adna_model.get(8, 25, b'C', b'T', 40));
        assert_approx_eq!(-6.920301, adna_model.get(9, 25, b'C', b'T', 40));
        assert_approx_eq!(-3.228001, adna_model.get(10, 25, b'C', b'T', 10));
        assert_approx_eq!(-6.987523, adna_model.get(11, 25, b'C', b'T', 40));
        assert_approx_eq!(-3.229167, adna_model.get(12, 25, b'C', b'T', 10));
        assert_approx_eq!(-6.987523, adna_model.get(13, 25, b'C', b'T', 40));
        assert_approx_eq!(-6.966826, adna_model.get(14, 25, b'C', b'T', 40));
        assert_approx_eq!(-6.920301, adna_model.get(15, 25, b'C', b'T', 40));
        assert_approx_eq!(-6.825267, adna_model.get(16, 25, b'C', b'T', 40));
        assert_approx_eq!(-6.642854, adna_model.get(17, 25, b'C', b'T', 40));
        assert_approx_eq!(-6.320596, adna_model.get(18, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.813153, adna_model.get(19, 25, b'C', b'T', 40));
        assert_approx_eq!(-5.113328, adna_model.get(20, 25, b'C', b'T', 40));
        assert_approx_eq!(-4.257998, adna_model.get(21, 25, b'C', b'T', 40));
        assert_approx_eq!(-3.300748, adna_model.get(22, 25, b'C', b'T', 40));
        assert_approx_eq!(-2.285697, adna_model.get(23, 25, b'C', b'T', 40));
        assert_approx_eq!(-1.504131, adna_model.get(24, 25, b'C', b'T', 10));

        // Reference 'C', not substituted
        assert_approx_eq!(-1.199396, adna_model.get(0, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.355901, adna_model.get(1, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.675932, adna_model.get(2, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.098193, adna_model.get(3, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.062537, adna_model.get(4, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.566023, adna_model.get(5, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.038071, adna_model.get(6, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.034366, adna_model.get(7, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.032611, adna_model.get(8, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.031781, adna_model.get(9, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.553695, adna_model.get(10, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.031227, adna_model.get(11, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.553513, adna_model.get(12, 25, b'C', b'C', 10));
        assert_approx_eq!(-0.031227, adna_model.get(13, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.031395, adna_model.get(14, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.031781, adna_model.get(15, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.032611, adna_model.get(16, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.034366, adna_model.get(17, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.038071, adna_model.get(18, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.045904, adna_model.get(19, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.062537, adna_model.get(20, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.098193, adna_model.get(21, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.176268, adna_model.get(22, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.355901, adna_model.get(23, 25, b'C', b'C', 40));
        assert_approx_eq!(-1.199396, adna_model.get(24, 25, b'C', b'C', 10));

        // 'G' -> 'A'
        assert_approx_eq!(-3.230046, adna_model.get(0, 25, b'G', b'A', 10));
        assert_approx_eq!(-7.013649, adna_model.get(1, 25, b'G', b'A', 40));
        assert_approx_eq!(-3.230046, adna_model.get(2, 25, b'G', b'A', 10));
        assert_approx_eq!(-7.013649, adna_model.get(3, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(4, 25, b'G', b'A', 40));
        assert_approx_eq!(-3.230046, adna_model.get(5, 25, b'G', b'A', 10));
        assert_approx_eq!(-7.013649, adna_model.get(6, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(7, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(8, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(9, 25, b'G', b'A', 40));
        assert_approx_eq!(-3.230046, adna_model.get(10, 25, b'G', b'A', 10));
        assert_approx_eq!(-7.013649, adna_model.get(11, 25, b'G', b'A', 40));
        assert_approx_eq!(-3.230046, adna_model.get(12, 25, b'G', b'A', 10));
        assert_approx_eq!(-7.013649, adna_model.get(13, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(14, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(15, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(16, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(17, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(18, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(19, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(20, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(21, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(22, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.013649, adna_model.get(23, 25, b'G', b'A', 40));
        assert_approx_eq!(-3.230046, adna_model.get(24, 25, b'G', b'A', 10));

        // Reference 'G', not substituted
        assert_approx_eq!(-0.553375, adna_model.get(0, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.031019, adna_model.get(1, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.553375, adna_model.get(2, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.031019, adna_model.get(3, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(4, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.553375, adna_model.get(5, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.031019, adna_model.get(6, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(7, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(8, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(9, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.553375, adna_model.get(10, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.031019, adna_model.get(11, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.553375, adna_model.get(12, 25, b'G', b'G', 10));
        assert_approx_eq!(-0.031019, adna_model.get(13, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(14, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(15, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(16, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(17, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(18, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(19, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(20, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(21, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(22, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(23, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.553375, adna_model.get(24, 25, b'G', b'G', 10));

        // Random samples
        assert_approx_eq!(-3.237864, adna_model.get(5, 25, b'A', b'C', 10));
        assert_approx_eq!(-3.230046, adna_model.get(10, 25, b'G', b'A', 10));
        assert_approx_eq!(-7.207481, adna_model.get(23, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.207481, adna_model.get(9, 25, b'A', b'T', 40));
        assert_approx_eq!(-3.237864, adna_model.get(10, 25, b'A', b'C', 10));
        assert_approx_eq!(-7.207481, adna_model.get(18, 25, b'C', b'G', 40));
        assert_approx_eq!(-2.285697, adna_model.get(23, 25, b'C', b'T', 40));
        assert_approx_eq!(-7.207481, adna_model.get(21, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.207481, adna_model.get(21, 25, b'T', b'G', 40));
        assert_approx_eq!(-0.552156, adna_model.get(5, 25, b'A', b'A', 10));
        assert_approx_eq!(-7.207481, adna_model.get(9, 25, b'G', b'C', 40));
        assert_approx_eq!(-7.013649, adna_model.get(9, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(17, 25, b'G', b'C', 40));
        assert_approx_eq!(-3.237864, adna_model.get(5, 25, b'T', b'G', 10));
        assert_approx_eq!(-7.207481, adna_model.get(21, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.207481, adna_model.get(21, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(23, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.207481, adna_model.get(19, 25, b'T', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(15, 25, b'G', b'G', 40));
        assert_approx_eq!(-7.207481, adna_model.get(3, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(15, 25, b'C', b'A', 40));
        assert_approx_eq!(-0.031227, adna_model.get(11, 25, b'C', b'C', 40));
        assert_approx_eq!(-7.207481, adna_model.get(20, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.207481, adna_model.get(17, 25, b'C', b'G', 40));
        assert_approx_eq!(-3.237864, adna_model.get(12, 25, b'A', b'C', 10));
        assert_approx_eq!(-7.207481, adna_model.get(18, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.207481, adna_model.get(7, 25, b'C', b'A', 40));
        assert_approx_eq!(-0.553375, adna_model.get(0, 25, b'G', b'G', 10));
        assert_approx_eq!(-3.237864, adna_model.get(24, 25, b'A', b'T', 10));
        assert_approx_eq!(-7.207481, adna_model.get(6, 25, b'A', b'T', 40));
        assert_approx_eq!(-0.029585, adna_model.get(19, 25, b'A', b'A', 40));
        assert_approx_eq!(-3.237864, adna_model.get(2, 25, b'A', b'G', 10));
        assert_approx_eq!(-7.207481, adna_model.get(20, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(4, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(1, 25, b'C', b'A', 40));
        assert_approx_eq!(-3.230046, adna_model.get(0, 25, b'G', b'A', 10));
        assert_approx_eq!(-0.029585, adna_model.get(7, 25, b'T', b'T', 40));
        assert_approx_eq!(-7.207481, adna_model.get(15, 25, b'G', b'C', 40));
        assert_approx_eq!(-3.237864, adna_model.get(2, 25, b'G', b'C', 10));
        assert_approx_eq!(-3.237864, adna_model.get(0, 25, b'C', b'G', 10));
        assert_approx_eq!(-6.825267, adna_model.get(8, 25, b'C', b'T', 40));
        assert_approx_eq!(-7.207481, adna_model.get(1, 25, b'G', b'T', 40));
        assert_approx_eq!(-7.207481, adna_model.get(19, 25, b'T', b'C', 40));
        assert_approx_eq!(-7.207481, adna_model.get(14, 25, b'C', b'A', 40));
        assert_approx_eq!(-3.237864, adna_model.get(10, 25, b'A', b'G', 10));
        assert_approx_eq!(-3.237864, adna_model.get(12, 25, b'T', b'C', 10));
        assert_approx_eq!(-0.029585, adna_model.get(14, 25, b'A', b'A', 40));
        assert_approx_eq!(-0.029585, adna_model.get(1, 25, b'A', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(7, 25, b'G', b'C', 40));
        assert_approx_eq!(-3.237864, adna_model.get(2, 25, b'T', b'A', 10));
        assert_approx_eq!(-7.207481, adna_model.get(11, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(19, 25, b'C', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(6, 25, b'T', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(6, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.029585, adna_model.get(11, 25, b'A', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(19, 25, b'T', b'A', 40));
        assert_approx_eq!(-0.029585, adna_model.get(3, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.031019, adna_model.get(3, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.029585, adna_model.get(1, 25, b'T', b'T', 40));
        assert_approx_eq!(-0.031019, adna_model.get(22, 25, b'G', b'G', 40));
        assert_approx_eq!(-7.207481, adna_model.get(17, 25, b'A', b'C', 40));
        assert_approx_eq!(-7.207481, adna_model.get(1, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.207481, adna_model.get(7, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(13, 25, b'T', b'G', 40));
        assert_approx_eq!(-3.237864, adna_model.get(12, 25, b'T', b'A', 10));
        assert_approx_eq!(-0.553695, adna_model.get(10, 25, b'C', b'C', 10));
        assert_approx_eq!(-7.207481, adna_model.get(7, 25, b'A', b'C', 40));
        assert_approx_eq!(-3.237864, adna_model.get(5, 25, b'G', b'T', 10));
        assert_approx_eq!(-0.355901, adna_model.get(1, 25, b'C', b'C', 40));
        assert_approx_eq!(-0.032611, adna_model.get(16, 25, b'C', b'C', 40));
        assert_approx_eq!(-7.207481, adna_model.get(9, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.207481, adna_model.get(4, 25, b'T', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(8, 25, b'A', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(1, 25, b'G', b'G', 40));
        assert_approx_eq!(-6.825267, adna_model.get(16, 25, b'C', b'T', 40));
        assert_approx_eq!(-7.207481, adna_model.get(16, 25, b'A', b'T', 40));
        assert_approx_eq!(-0.031019, adna_model.get(11, 25, b'G', b'G', 40));
        assert_approx_eq!(-0.031019, adna_model.get(8, 25, b'G', b'G', 40));
        assert_approx_eq!(-3.151698, adna_model.get(5, 25, b'C', b'T', 10));
        assert_approx_eq!(-7.013649, adna_model.get(21, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(8, 25, b'C', b'G', 40));
        assert_approx_eq!(-3.230046, adna_model.get(24, 25, b'G', b'A', 10));
        assert_approx_eq!(-7.207481, adna_model.get(7, 25, b'C', b'G', 40));
        assert_approx_eq!(-7.013649, adna_model.get(18, 25, b'G', b'A', 40));
        assert_approx_eq!(-7.207481, adna_model.get(1, 25, b'A', b'T', 40));
        assert_approx_eq!(-0.029585, adna_model.get(14, 25, b'T', b'T', 40));
        assert_approx_eq!(-6.642854, adna_model.get(17, 25, b'C', b'T', 40));
        assert_approx_eq!(-6.320596, adna_model.get(18, 25, b'C', b'T', 40));
        assert_approx_eq!(-3.237864, adna_model.get(0, 25, b'A', b'G', 10));
        assert_approx_eq!(-0.031019, adna_model.get(4, 25, b'G', b'G', 40));
        assert_approx_eq!(-3.237864, adna_model.get(24, 25, b'T', b'A', 10));
        assert_approx_eq!(-7.207481, adna_model.get(14, 25, b'T', b'G', 40));
        assert_approx_eq!(-7.207481, adna_model.get(16, 25, b'T', b'A', 40));
        assert_approx_eq!(-0.566023, adna_model.get(5, 25, b'C', b'C', 10));
        assert_approx_eq!(-3.237864, adna_model.get(12, 25, b'G', b'C', 10));
        assert_approx_eq!(-7.207481, adna_model.get(14, 25, b'A', b'G', 40));
        assert_approx_eq!(-7.207481, adna_model.get(17, 25, b'T', b'C', 40));
        assert_approx_eq!(-1.504131, adna_model.get(0, 25, b'C', b'T', 10));
        assert_approx_eq!(-7.207481, adna_model.get(15, 25, b'T', b'A', 40));
        assert_approx_eq!(-3.237864, adna_model.get(5, 25, b'C', b'G', 10));
    }
}
