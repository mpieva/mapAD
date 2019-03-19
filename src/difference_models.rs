use std::collections::HashMap;

pub trait DifferenceModel {
    fn new_default() -> Self;
    fn get(&self, i: usize, from: u8, to: u8) -> f32;
    fn get_min_penalty(&self, i: usize, to: u8) -> f32 {
        b"ACGT"
            .iter()
            .filter(|&&base| base != to)
            .map(|&base| self.get(i, base, to))
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

pub struct SimplisticVindijaPattern {
    /// substitution_matrix[&to, &from]
    substitution_matrix: HashMap<u8, HashMap<u8, f32>>,
    read_end: [f32; 7],
}

impl DifferenceModel for SimplisticVindijaPattern {
    fn new_default() -> Self {
        let substitution_matrix = hashmap!(
                b'A' => hashmap!(
                    b'A' => 1.0,
                    b'C' => -1.0,
                    b'G' => -1.0,
                    b'T' => -1.0,
                ),
                b'C' => hashmap!(
                    b'A' => -1.0,
                    b'C' => 1.0,
                    b'G' => -1.0,
                    b'T' => -1.0,
                ),
                b'G' => hashmap!(
                    b'A' => -1.0,
                    b'C' => -1.0,
                    b'G' => 1.0,
                    b'T' => -1.0,
                ),
                b'T' => hashmap!(
                    b'A' => -1.0,
                    b'C' => -1.0,
                    b'G' => -1.0,
                    b'T' => 1.0,
                ),
        );
        let read_end = [0.4, 0.25, 0.1, 0.06, 0.05, 0.04, 0.03];
        SimplisticVindijaPattern {
            substitution_matrix,
            read_end,
        }
    }
    fn get(&self, i: usize, from: u8, to: u8) -> f32 {
        match from {
            b'C' => match to {
                b'T' => self.deamination_helper(i),
                b'C' => 1.0 - self.deamination_helper(i),
                _ => self.substitution_matrix[&to][&from],
            },
            _ => self.substitution_matrix[&to][&from],
        }
    }
}

impl SimplisticVindijaPattern {
    fn deamination_helper(&self, i: usize) -> f32 {
        if i >= self.read_end.len() {
            0.02
        } else {
            self.read_end[i]
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::difference_models::DifferenceModel;
    use crate::difference_models::SimplisticVindijaPattern;

    #[test]
    fn test_simplistic_adna_damage_model() {
        let simplistic_model = SimplisticVindijaPattern::new_default();

        assert_eq!(0.4, simplistic_model.get(0, b'C', b'T'));
        assert_eq!(0.6, simplistic_model.get(0, b'C', b'C'));
        assert_eq!(0.02, simplistic_model.get(15, b'C', b'T'));
    }
}
