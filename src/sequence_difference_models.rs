// TODO: Use lookup tables instead of computing scores on-demand

/// Sequence difference models are expected to yield only non-positive values (and 0.0), for example log probabilities.
pub trait SequenceDifferenceModel {
    fn get(&self, i: usize, read_length: usize, from: u8, to: u8, base_quality: u8) -> f32;
    fn get_representative_mismatch_penalty(&self) -> f32 {
        self.get(40, 80, b'T', b'A', 40)
    }

    /// Needed for the calculation of D arrays
    fn get_min_penalty(&self, i: usize, read_length: usize, to: u8, base_quality: u8) -> f32 {
        b"ACGT"
            .iter()
            .filter(|&&base| base != to)
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

/// Briggs, A. W. et al. Patterns of damage in genomic DNA sequences from a Neandertal. PNAS 104, 14616–14621 (2007)
pub struct BriggsEtAl2007aDNA {
    pub library_prep: LibraryPrep,
    // Deamination rate in double-stranded stems
    pub ds_deamination_rate: f32,
    // Deamination rate in single-stranded overhangs
    pub ss_deamination_rate: f32,
    pub divergence: f32,
}

impl SequenceDifferenceModel for BriggsEtAl2007aDNA {
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
        let divergence_and_seq_err =
            sequencing_error + self.divergence - sequencing_error * self.divergence;

        let c_to_t = self.ss_deamination_rate * p_fwd + self.ds_deamination_rate * (1.0 - p_fwd);
        let g_to_a = self.ss_deamination_rate * p_rev + self.ds_deamination_rate * (1.0 - p_rev);
        match from {
            b'A' => match to {
                b'A' => 1.0 - 3.0 * divergence_and_seq_err,
                _ => divergence_and_seq_err,
            },
            b'C' => match to {
                b'C' => {
                    1.0 - 3.0 * divergence_and_seq_err - c_to_t
                        + 4.0 * divergence_and_seq_err * c_to_t
                }
                b'T' => divergence_and_seq_err + c_to_t - 4.0 * divergence_and_seq_err * c_to_t,
                _ => divergence_and_seq_err,
            },
            b'G' => match to {
                b'A' => divergence_and_seq_err + g_to_a - 4.0 * divergence_and_seq_err * g_to_a,
                b'G' => {
                    1.0 - 3.0 * divergence_and_seq_err - g_to_a
                        + 4.0 * divergence_and_seq_err * g_to_a
                }
                _ => divergence_and_seq_err,
            },
            b'T' => match to {
                b'T' => 1.0 - 3.0 * divergence_and_seq_err,
                _ => divergence_and_seq_err,
            },
            _ => divergence_and_seq_err,
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
    fn test_briggs_model() {
        let briggs_model = BriggsEtAl2007aDNA {
            library_prep: (LibraryPrep::SingleStranded {
                five_prime_overhang: 0.63,
                three_prime_overhang: 0.8,
            }),
            ds_deamination_rate: 0.07,
            ss_deamination_rate: 0.6,
            divergence: 0.001,
        };
        let read_length = 35;

        //        println!("Briggs model");
        //        for i in 0..read_length {
        //            println!(
        //                "{i})\tC->T: {c_t}\t\tC->C: {c_c}\t\tA->A: {a_a}\t\t G->A: {g_a}",
        //                i = i,
        //                c_t = briggs_model.get(i, read_length, b'C', b'T', 40),
        //                c_c = briggs_model.get(i, read_length, b'C', b'C', 40),
        //                a_a = briggs_model.get(i, read_length, b'A', b'A', 40),
        //                g_a = briggs_model.get(i, read_length, b'G', b'A', 40),
        //            );
        //        }

        assert_approx_eq!(-1.310067, briggs_model.get(0, read_length, b'C', b'T', 40));
        assert_approx_eq!(-0.750255, briggs_model.get(0, read_length, b'C', b'C', 40));
        assert_approx_eq!(-3.695316, briggs_model.get(15, read_length, b'C', b'T', 40));
        assert_approx_eq!(-9.828412, briggs_model.get(15, read_length, b'G', b'C', 40));
        assert_approx_eq!(-0.004768, briggs_model.get(15, read_length, b'A', b'A', 40));
    }
}
