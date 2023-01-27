/// <https://chrismaclellan.com/blog/lazy-shuffled-list-generator>
/// <https://en.wikipedia.org/wiki/Lehmer_random_number_generator>
/// <https://en.wikipedia.org/wiki/Linear_congruential_generator>
use std::ops::Range;

pub struct PrRange {
    start: usize,
    l: usize,
    m: usize,
    a: usize,
    x: usize,
    seed: usize,
    count: usize,
}

impl Iterator for PrRange {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.count == 0 && self.l == 1 {
            self.count += 1;
            return Some(self.start);
        }
        loop {
            let prev_x = self.x;
            self.x = (self.a * self.x) % self.m;
            if self.count > 0 && prev_x == self.seed {
                return None;
            }
            if prev_x <= self.l {
                self.count += 1;
                return Some(prev_x - 1 + self.start);
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.l - self.count, Some(self.l - self.count))
    }
}

impl PrRange {
    pub fn try_new(start: usize, end: usize, seed: usize) -> Option<Self> {
        // We don't allow empty or negative ranges
        let l = end.saturating_sub(start);
        if l == 0 {
            return None;
        }

        let m = next_prime(l);

        let a = {
            let mut a = 2;
            while !is_primitive_root(a, m)? {
                a += 1;
            }
            a
        };

        let seed = (seed % l).max(1);

        Some(Self {
            start,
            l,
            m,
            a,
            x: seed,
            seed,
            count: 0,
        })
    }

    pub fn try_from_range<Idx>(range: &Range<Idx>, seed: usize) -> Option<Self>
    where
        Idx: Into<usize> + Copy,
    {
        Self::try_new(range.start.into(), range.end.into(), seed)
    }
}

fn next_prime(n: usize) -> usize {
    let mut p = n + 1;
    if p <= 2 {
        return 2;
    }
    if p % 2 == 0 {
        p += 1;
    }
    while !is_prime(p) {
        p += 2;
    }
    p
}

fn is_prime(n: usize) -> bool {
    if n <= 1 {
        return false;
    } else if n <= 3 {
        return true;
    } else if n % 2 == 0 || n % 3 == 0 {
        return false;
    }
    let mut i = 5;
    while i * i <= n {
        if n % i == 0 || n % (i + 2) == 0 {
            return false;
        }
        i += 6;
    }
    true
}

fn is_primitive_root(a: usize, n: usize) -> Option<bool> {
    let phi = n - 1;
    for p in PrimeFactorIterator::new(phi) {
        if checked_pow_mod(a, phi / p, n)? == 1 {
            return Some(false);
        }
    }
    Some(true)
}

struct PrimeFactorIterator {
    n: usize,
    i: usize,
    step: usize,
    last: usize,
}

impl Iterator for PrimeFactorIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.n <= 3 {
            return None;
        }
        while self.i * self.i <= self.n {
            while self.n > 1 {
                while self.n % self.i == 0 {
                    if self.i > self.last {
                        let out = self.i;
                        self.last = self.i;
                        return Some(out);
                    }
                    self.n /= self.i;
                }
                self.i += self.step;
                self.step = 2;
            }
        }
        None
    }
}

impl PrimeFactorIterator {
    fn new(n: usize) -> Self {
        Self {
            n,
            i: 2,
            step: 1,
            last: 0,
        }
    }
}

/// Checked modular exponentiation. Computes base.pow(exponent) % modulus,
/// returning `None` if overflow occurred.
fn checked_pow_mod(mut base: usize, mut exponent: usize, modulus: usize) -> Option<usize> {
    if modulus == 1 {
        return Some(0);
    }
    // Overflow check
    (modulus - 1).checked_mul(modulus - 1)?;
    let mut result = 1;
    base %= modulus;
    while exponent > 0 {
        if exponent % 2 == 1 {
            result = (result * base) % modulus;
        }
        exponent >>= 1;
        base = (base * base) % modulus;
    }
    Some(result)
}

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn large_numbers() {
        let start = 6100000000;
        let end = 6100000005;
        let mut randrange = PrRange::try_new(start, end, 1234)
            .unwrap()
            .collect::<Vec<_>>();
        randrange.sort_unstable();
        assert_eq!(randrange, (start..end).collect::<Vec<_>>());
    }

    #[test]
    fn try_from_range() {
        let start = 13_u8;
        let end = 23_u8;
        let range = start..end;
        let mut randrange = PrRange::try_from_range(&range, 1234)
            .unwrap()
            .collect::<Vec<_>>();
        randrange.sort_unstable();
        assert_eq!(
            randrange,
            ((start as usize)..(end as usize)).collect::<Vec<_>>()
        );
    }

    #[test]
    fn itw_bug_case_1() {
        let range = 5233065207..5233065216_usize;
        let randrange = PrRange::try_from_range(&range, 400636091).unwrap();
        assert_eq!(randrange.count(), 9);
    }

    #[test]
    fn small_range() {
        let start = 1;
        let end = 2;
        let mut randrange = PrRange::try_new(start, end, 1234)
            .unwrap()
            .collect::<Vec<_>>();
        randrange.sort_unstable();
        assert_eq!(randrange, (start..end).collect::<Vec<_>>());
    }

    #[test]
    fn invalid_randrange_1() {
        let start = 1;
        let end = 0;
        assert!(PrRange::try_new(start, end, 1234).is_none());
    }

    #[test]
    fn invalid_randrange_2() {
        let start = 1;
        let end = 1;
        assert!(PrRange::try_new(start, end, 1234).is_none());
    }

    #[test]
    fn test_seeds() {
        let from = 0;
        let to = 100;
        for start in from..=to {
            for end in start + 1..=to {
                for seed in from..=to {
                    let randrange = PrRange::try_new(start, end, seed).unwrap();
                    assert_eq!(randrange.count(), end - start);
                }
            }
        }
    }
}
