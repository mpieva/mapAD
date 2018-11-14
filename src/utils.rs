use std::collections::HashMap;

pub struct AlignmentParameters {
    pub base_error_rate: f32,
    pub poisson_threshold: f32,
    pub penalty_mismatch: i32,
    pub penalty_gap_open: i32,
    pub penalty_gap_extend: i32,
    pub penalty_c_t: i32,
    pub penalty_g_a: i32,
}

pub struct AllowedMismatches<'a> {
    alignment_parameters: &'a AlignmentParameters,
    cache: HashMap<usize, i32>,
}

impl<'a> AllowedMismatches<'a> {
    pub fn new(alignment_parameters: &AlignmentParameters) -> AllowedMismatches {
        AllowedMismatches {
            alignment_parameters,
            cache: HashMap::new(),
        }
    }

    pub fn get(&mut self, read_length: usize) -> i32 {
        match self.cache.get(&read_length) {
            None => {
                let max_num_mismatches = self.calculate_max_num_mismatches(read_length);
                self.cache.insert(read_length, max_num_mismatches);
                max_num_mismatches
            }
            Some(v) => *v,
        }
    }

    fn factorial(n: u32) -> u32 {
        (2..=n).product()
    }

    fn ppoisson(lambda: f32, k: u32) -> f32 {
        let lower_tail: f32 = (0..=k)
            .map(|i| {
                lambda.powi(i as i32) * (-lambda).exp() / AllowedMismatches::factorial(i) as f32
            })
            .sum();
        1.0 - lower_tail
    }

    fn calculate_max_num_mismatches(&self, read_length: usize) -> i32 {
        (1..=read_length as i32)
            .take_while(|&i| {
                AllowedMismatches::ppoisson(
                    self.alignment_parameters.base_error_rate * read_length as f32,
                    i as u32,
                ) > self.alignment_parameters.poisson_threshold
            })
            .last()
            .unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_allowed_mismatches() {
        let parameters = AlignmentParameters {
            base_error_rate: 0.02,
            poisson_threshold: 0.04,
            penalty_mismatch: 1,
            penalty_gap_open: 1,
            penalty_gap_extend: 1,
            penalty_c_t: 0,
            penalty_g_a: 0,
        };
        let mut allowed_mismatches = AllowedMismatches::new(&parameters);

        assert_eq!(4, allowed_mismatches.get(100));
        assert_eq!(2, allowed_mismatches.get(38));
        assert_eq!(0, allowed_mismatches.get(1));
        assert_eq!(0, allowed_mismatches.get(0));
    }
}
