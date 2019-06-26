// TODO: Use lookup tables instead of computing scores on-demand

/// Sequence difference models are expected to yield only non-positive values (and 0.0), for example log probabilities.
pub trait SequenceDifferenceModel {
    fn get(&self, i: usize, read_length: usize, from: u8, to: u8, base_quality: u8) -> f32;
    fn get_representative_mismatch_penalty(&self) -> f32 {
        self.get(40, 80, b'T', b'A', 40)
    }

    /// Needed for the calculation of D arrays
    fn get_min_penalty(&self, i: usize, read_length: usize, to: u8) -> f32 {
        b"ACGT"
            .iter()
            .filter(|&&base| base != to)
            .map(|&base| self.get(i, read_length, base, to, 40))
            .filter(|&penalty| penalty < 0.0)
            .fold(std::f32::MIN, |acc, v| acc.max(v))
    }
}

///// Briggs, A. W. et al. Patterns of damage in genomic DNA sequences from a Neandertal. PNAS 104, 14616–14621 (2007)
//struct BriggsEtAl2007aDNA {
//    limit: u32,
//    lparam: Option<u32>,
//    rparam: Option<u32>,
//    ds_deam: f32,
//    ss_deam: f32,
//    div: f32,
//    gap: u32,
//    hit_limit: u32,
//    cnt_limit: u32,
//}
//
//impl DnaDamageModel for BriggsEtAl2007aDNA {
//    fn get(&self, i: usize, from: u8, to: u8) -> f32 {}
//
//    fn default() -> Self {
//        BriggsEtAl2007aDNA {
//            limit: 12, // about two mismatches?
//            lparam: None,
//            rparam: None,
//            ds_deam: 0.02,
//            ss_deam: 0.45,
//            div: 0.015 / 3.0, // approx. human-chimp
//            gap: 11,
//            hit_limit: 64,
//            cnt_limit: 256,
//        }
//    }
//}

/// Very simple model of ancient DNA degradation for starters.
/// It only takes C->T deaminations into account and assumes
/// symmetry between 5' and 3' ends of the deamination pattern.
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
        //        for i in 0..read_length {
        //            println!(
        //                "{i})\tC->T: {c_t}\t\tC->C: {c_c}\t\tA->A: {a_a}\t\t G->A: {g_a}",
        //                i = i,
        //                c_t = vindija_pwm.get(i, read_length, b'C', b'T'),
        //                c_c = vindija_pwm.get(i, read_length, b'C', b'C'),
        //                a_a = vindija_pwm.get(i, read_length, b'A', b'A'),
        //                g_a = vindija_pwm.get(i, read_length, b'G', b'A'),
        //            );
        //        }

        assert_approx_eq!(-1.321928, vindija_pwm.get(0, read_length, b'C', b'T'));
        assert_approx_eq!(-0.736965, vindija_pwm.get(0, read_length, b'C', b'C'));
        assert_approx_eq!(-5.643856, vindija_pwm.get(15, read_length, b'C', b'T'));
        assert_approx_eq!(-10.965784, vindija_pwm.get(15, read_length, b'G', b'C'));
        assert_approx_eq!(-0.000721, vindija_pwm.get(15, read_length, b'A', b'A'));
    }
}
