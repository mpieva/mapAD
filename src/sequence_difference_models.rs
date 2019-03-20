pub trait SequenceDifferenceModel {
    fn new() -> Self;
    fn get(&self, i: usize, read_length: usize, from: u8, to: u8) -> f32;

    // TODO: Cache results
    fn get_min_penalty(&self, i: usize, read_length: usize, to: u8) -> f32 {
        b"ACGT"
            .iter()
            .filter(|&&base| base != to) // TODO: Don't filter for mismatches
            .map(|&base| self.get(i, read_length, base, to))
            .filter(|&penalty| penalty <= (0.0 + std::f32::EPSILON))
            .fold(std::f32::MIN, |acc: f32, v| acc.max(v))
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
    background_substitution_probability: f32,
    observed_substitution_probability_default: f32,
}

impl SequenceDifferenceModel for VindijaPWM {
    fn new() -> Self {
        VindijaPWM {
            ppm_read_ends_symmetric_ct: [0.4, 0.25, 0.1, 0.06, 0.05, 0.04, 0.03],
            position_probability_ct_default: 0.02,
            background_substitution_probability: 0.25,
            observed_substitution_probability_default: 0.0005,
        }
    }
    fn get(&self, i: usize, read_length: usize, from: u8, to: u8) -> f32 {
        let position_probability = match from {
            b'C' => {
                let i = i.min(read_length - (i + 1));
                let position_probability_ct = *self
                    .ppm_read_ends_symmetric_ct
                    .get(i)
                    .unwrap_or(&self.position_probability_ct_default);
                match to {
                    b'T' => position_probability_ct,
                    b'C' => 1.0 - position_probability_ct,
                    // "Normal" mismatch
                    _ => self.observed_substitution_probability_default,
                }
            }
            // "Normal" mismatch
            _ => {
                if from == to {
                    1.0
                } else {
                    self.observed_substitution_probability_default
                }
            }
        };
        (position_probability / self.background_substitution_probability).log2()
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

        //        for i in 0..read_length {
        //            let c_t = vindija_pwm.get(i, read_length, b'C', b'T');
        //            let c_c = vindija_pwm.get(i, read_length, b'C', b'C');
        //            let a_a = vindija_pwm.get(i, read_length, b'A', b'A');
        //            let g_a = vindija_pwm.get(i, read_length, b'G', b'A');
        //            println!(
        //                "{})\tC->T: {}\t\tC->C: {}\t\tA->A: {}\t\t G->A: {}",
        //                i, c_t, c_c, a_a, g_a
        //            );
        //        }

        assert_approx_eq!(0.678072, vindija_pwm.get(0, read_length, b'C', b'T'));
        assert_approx_eq!(1.2630345, vindija_pwm.get(0, read_length, b'C', b'C'));
        assert_approx_eq!(-3.6438563, vindija_pwm.get(15, read_length, b'C', b'T'));
    }
}
