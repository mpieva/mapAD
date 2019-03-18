use std::collections::HashMap;

pub struct AlignmentParameters {
    pub base_error_rate: f64,
    pub poisson_threshold: f64,
    pub penalty_c_t: i32,
    pub penalty_g_a: i32,
    pub penalty_mismatch: f32,
    pub penalty_gap_open: f32,
    pub penalty_gap_extend: f32,
}

pub struct AllowedMismatches<'a> {
    alignment_parameters: &'a AlignmentParameters,
    cache: HashMap<usize, f32>,
}

impl<'a> AllowedMismatches<'a> {
    pub fn new(alignment_parameters: &AlignmentParameters) -> AllowedMismatches {
        AllowedMismatches {
            alignment_parameters,
            cache: HashMap::new(),
        }
    }

    pub fn get(&mut self, read_length: usize) -> f32 {
        match self.cache.get(&read_length) {
            None => {
                let max_num_mismatches = self.calculate_max_num_mismatches(read_length);
                self.cache.insert(read_length, max_num_mismatches);
                max_num_mismatches
            }
            Some(v) => *v,
        }
    }

    fn calculate_max_num_mismatches(&self, read_length: usize) -> f32 {
        let lambda = read_length as f64 * self.alignment_parameters.base_error_rate;
        let exp_minus_lambda = (-lambda).exp();

        let mut lambda_to_the_power_of_k = 1.0;
        let mut k_factorial = 1;
        let mut sum = exp_minus_lambda;

        // k = 0. BWA allows k+1 mismatches, and so do we
        [(0 + 1, sum)]
            .iter()
            .cloned()
            // k = 1..read_length
            .chain((1..=read_length as u64).map(|k| {
                lambda_to_the_power_of_k *= lambda;
                k_factorial *= k;
                sum += lambda_to_the_power_of_k * exp_minus_lambda / k_factorial as f64;
                // BWA allows k+1 mismatches, and so do we
                (k + 1, sum)
            }))
            .take_while(|(_k, sum)| 1.0 - sum > self.alignment_parameters.poisson_threshold)
            .map(|(k, _sum)| k)
            .last()
            .unwrap_or(0) as f32
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
            penalty_c_t: 0,
            penalty_g_a: 0,
            penalty_mismatch: 1.0,
            penalty_gap_open: 1.0,
            penalty_gap_extend: 1.0,
        };
        let mut allowed_mismatches = AllowedMismatches::new(&parameters);

        assert_eq!(6.0, allowed_mismatches.get(156));
        assert_eq!(6.0, allowed_mismatches.get(124));
        assert_eq!(5.0, allowed_mismatches.get(123));
        assert_eq!(5.0, allowed_mismatches.get(93));
        assert_eq!(4.0, allowed_mismatches.get(92));
        assert_eq!(4.0, allowed_mismatches.get(64));
        assert_eq!(3.0, allowed_mismatches.get(63));
        assert_eq!(3.0, allowed_mismatches.get(38));
        assert_eq!(2.0, allowed_mismatches.get(37));
        assert_eq!(2.0, allowed_mismatches.get(16));
        assert_eq!(1.0, allowed_mismatches.get(15));
        assert_eq!(1.0, allowed_mismatches.get(3));
        assert_eq!(0.0, allowed_mismatches.get(2));
        assert_eq!(0.0, allowed_mismatches.get(0));
    }
}
