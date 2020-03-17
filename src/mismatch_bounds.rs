use enum_dispatch::enum_dispatch;
use serde::{Deserialize, Serialize};
use std::{
    fmt,
    fmt::{Display, Formatter},
    iter::once,
};

const MAX_CACHED_READ_LENGTH: usize = 256;

#[enum_dispatch]
pub trait MismatchBound {
    fn reject(&self, value: f32, read_length: usize) -> bool;

    /// If the best scoring interval has a total sum of penalties z, do not search
    /// for hits with a minimal expected scored worse than z + representative_mismatch
    /// to speed up the alignment of endogenous reads.
    fn reject_iterative(&self, value: f32, reference: f32) -> bool;
}

#[enum_dispatch(MismatchBound)]
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum MismatchBoundDispatch {
    Continuous,
    Discrete,
    TestBound,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Continuous {
    cutoff: f32,
    exponent: f32,
    representative_mismatch_penalty: f32,
    cache: Vec<f32>,
}

impl MismatchBound for Continuous {
    fn reject(&self, value: f32, read_length: usize) -> bool {
        let scaled_read_length = match self.cache.get(read_length) {
            Some(&v) => v,
            None => (read_length as f32).powf(self.exponent),
        };
        (value / scaled_read_length) < self.cutoff
    }

    fn reject_iterative(&self, value: f32, reference: f32) -> bool {
        value < reference + self.representative_mismatch_penalty
    }
}

impl Continuous {
    pub fn new(cutoff: f32, exponent: f32, representative_mismatch_penalty: f32) -> Self {
        let cache = (0..MAX_CACHED_READ_LENGTH)
            .map(|read_length| (read_length as f32).powf(exponent))
            .collect();

        Self {
            cutoff,
            exponent,
            representative_mismatch_penalty,
            cache,
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Discrete {
    poisson_threshold: f32,
    base_error_rate: f32,
    representative_mismatch_penalty: f32,
    cache: Vec<f32>,
}

impl MismatchBound for Discrete {
    fn reject(&self, value: f32, read_length: usize) -> bool {
        // Poisson-modelled number (BWA-like) times representative mismatch penalty
        value < self.get(read_length) * self.representative_mismatch_penalty
    }

    fn reject_iterative(&self, value: f32, reference: f32) -> bool {
        value < reference + self.representative_mismatch_penalty
    }
}

impl Display for Discrete {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let width = (MAX_CACHED_READ_LENGTH as f32).log10().ceil() as usize;

        let text = {
            let mut tmp = (Self::MIN_READ_LENGTH..=MAX_CACHED_READ_LENGTH)
                .map(|read_length| (read_length, self.get(read_length)))
                .scan(
                    std::f32::MIN,
                    |previous, (read_length, allowed_mismatches)| {
                        if (allowed_mismatches - *previous).abs() > std::f32::EPSILON {
                            *previous = allowed_mismatches;
                            Some(Some((read_length, allowed_mismatches)))
                        } else {
                            Some(None)
                        }
                    },
                )
                .flatten()
                .map(|(read_length, allowed_mismatches)| {
                    format!(
                        "{:>width$} bp:\t{} {}\n",
                        read_length,
                        allowed_mismatches,
                        if allowed_mismatches > 1.0 + std::f32::EPSILON {
                            "mismatches"
                        } else {
                            "mismatch"
                        },
                        width = width,
                    )
                })
                .collect::<String>();
            let _ = tmp.pop();
            tmp
        };

        write!(f, "{}", text)
    }
}

impl Discrete {
    // 16 is the theoretical minimum for mapping uniquely to a random genome of human-genome-like size
    const MIN_READ_LENGTH: usize = 17;

    pub fn new(
        poisson_threshold: f32,
        base_error_rate: f32,
        representative_mismatch_penalty: f32,
    ) -> Discrete {
        let cache = (0..MAX_CACHED_READ_LENGTH)
            .map(|read_length| {
                Self::calculate_max_num_mismatches(
                    // Read lengths are stored with an offset to cache the "hot" lengths
                    read_length + Self::MIN_READ_LENGTH,
                    poisson_threshold,
                    base_error_rate,
                )
            })
            .collect();

        Self {
            poisson_threshold,
            base_error_rate,
            representative_mismatch_penalty,
            cache,
        }
    }

    fn calculate_max_num_mismatches(
        read_length: usize,
        poisson_threshold: f32,
        base_error_rate: f32,
    ) -> f32 {
        let lambda = read_length as f32 * base_error_rate;
        let exp_minus_lambda = (-lambda).exp();

        // k = 0 (here 1, because BWA allows k+1 mismatches, and so do we)
        once((1, exp_minus_lambda))
            // k = 1..read_length
            .chain((1..=read_length as u64).scan(
                (1.0, 1, exp_minus_lambda),
                |(lambda_to_the_k, k_factorial, sum), k| {
                    *lambda_to_the_k *= lambda;
                    *k_factorial *= k;
                    *sum += *lambda_to_the_k * exp_minus_lambda / *k_factorial as f32;
                    // BWA allows k+1 mismatches, and so do we
                    Some((k + 1, *sum))
                },
            ))
            .take_while(|(_k, sum)| 1.0 - *sum > poisson_threshold)
            .last()
            .map(|(k, _sum)| k)
            .unwrap_or(0) as f32
    }

    fn get(&self, read_length: usize) -> f32 {
        // Reads shorter than the threshold are not allowed to differ
        if read_length < Self::MIN_READ_LENGTH {
            return 0.0;
        }

        match self
            .cache
            // An offset must be subtracted to point to the correct cache entries
            .get(read_length - Self::MIN_READ_LENGTH)
        {
            None => Self::calculate_max_num_mismatches(
                read_length,
                self.poisson_threshold,
                self.base_error_rate,
            ),
            Some(v) => *v,
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct TestBound {
    pub threshold: f32,
}

impl MismatchBound for TestBound {
    fn reject(&self, value: f32, _read_length: usize) -> bool {
        value < self.threshold
    }

    fn reject_iterative(&self, _value: f32, _reference: f32) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence_difference_models::*;

    #[test]
    fn test_allowed_mismatches() {
        let difference_model = VindijaPWM::new();
        let repr_mm_penalty = difference_model.get_representative_mismatch_penalty();
        let mismatch_bound = Discrete::new(0.04, 0.02, repr_mm_penalty);

        assert_eq!(mismatch_bound.get(156), 6.0);
        assert_eq!(mismatch_bound.get(124), 6.0);
        assert_eq!(mismatch_bound.get(123), 5.0);
        assert_eq!(mismatch_bound.get(93), 5.0);
        assert_eq!(mismatch_bound.get(92), 4.0);
        assert_eq!(mismatch_bound.get(64), 4.0);
        assert_eq!(mismatch_bound.get(63), 3.0);
        assert_eq!(mismatch_bound.get(38), 3.0);
        assert_eq!(mismatch_bound.get(37), 2.0);
        assert_eq!(mismatch_bound.get(17), 2.0);
        assert_eq!(mismatch_bound.get(16), 0.0);
        assert_eq!(mismatch_bound.get(15), 0.0);
        assert_eq!(mismatch_bound.get(3), 0.0);
        assert_eq!(mismatch_bound.get(2), 0.0);
        assert_eq!(mismatch_bound.get(0), 0.0);
    }

    #[test]
    fn test_allowed_mismatches_bwa_ancient_parameters() {
        let difference_model = VindijaPWM::new();
        let repr_mm_penalty = difference_model.get_representative_mismatch_penalty();
        let mismatch_bound = Discrete::new(0.01, 0.02, repr_mm_penalty);

        assert_eq!(mismatch_bound.get(207), 10.0);
        assert_eq!(mismatch_bound.get(176), 9.0);
        assert_eq!(mismatch_bound.get(146), 8.0);
        assert_eq!(mismatch_bound.get(117), 7.0);
        assert_eq!(mismatch_bound.get(90), 6.0);
        assert_eq!(mismatch_bound.get(64), 5.0);
        assert_eq!(mismatch_bound.get(42), 4.0);
        assert_eq!(mismatch_bound.get(22), 3.0);
        assert_eq!(mismatch_bound.get(17), 2.0);
        assert_eq!(mismatch_bound.get(8), 0.0);
        assert_eq!(mismatch_bound.get(1), 0.0);
    }

    #[test]
    fn test_display() {
        let difference_model = VindijaPWM::new();
        let representative_mismatch_boundary =
            difference_model.get_representative_mismatch_penalty();
        let mismatch_bound = Discrete::new(0.06, 0.02, representative_mismatch_boundary);

        let comparison = " 17 bp:\t1 mismatch
 20 bp:\t2 mismatches
 45 bp:\t3 mismatches
 73 bp:\t4 mismatches
104 bp:\t5 mismatches
137 bp:\t6 mismatches
172 bp:\t7 mismatches
208 bp:\t8 mismatches
244 bp:\t9 mismatches";

        assert_eq!(comparison, format!("{}", mismatch_bound));
    }
}
