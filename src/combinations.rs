/// All ways of placing `count` fermions on `sites` sites, as bit masks ranked in increasing order.
///
/// Ranking splits masks into a low and a high half (H. Q. Lin, Phys. Rev. B 42, 6561), so
/// it costs two table lookups.
pub struct Combinations {
    sites: u32,
    count: u32,
    binomials: Vec<usize>,
    half: u32,
    low_ranks: Vec<usize>,
    high_ranks: Vec<usize>,
}

impl Combinations {
    pub fn new(sites: u32, count: u32) -> Self {
        let columns = count as usize + 1;
        let mut binomials = vec![0; (sites as usize + 1) * columns];
        for n in 0..=sites as usize {
            binomials[n * columns] = 1;
            for k in 1..=(count as usize).min(n) {
                binomials[n * columns + k] =
                    binomials[(n - 1) * columns + k - 1] + binomials[(n - 1) * columns + k];
            }
        }
        let half = sites / 2;
        let mut c = Combinations {
            sites,
            count,
            binomials,
            half,
            low_ranks: Vec::new(),
            high_ranks: Vec::new(),
        };
        c.low_ranks = (0..1u32 << half).map(|low| c.partial_rank(low, 1)).collect();
        c.high_ranks = (0..1u32 << (sites - half))
            .map(|high| c.partial_rank(high << half, (count + 1).saturating_sub(high.count_ones())))
            .collect();
        c
    }

    /// `Σ C(position, ordinal)` over the set bits of `mask`, with ordinals counted from `first`.
    fn partial_rank(&self, mut mask: u32, first: u32) -> usize {
        if first == 0 || first + mask.count_ones() > self.count + 1 {
            return 0;
        }
        let mut rank = 0;
        for k in first..first + mask.count_ones() {
            rank += self.binomial(mask.trailing_zeros(), k);
            mask &= mask - 1;
        }
        rank
    }

    fn binomial(&self, n: u32, k: u32) -> usize {
        self.binomials[(n * (self.count + 1) + k) as usize]
    }

    pub fn len(&self) -> usize {
        self.binomial(self.sites, self.count)
    }

    pub fn rank(&self, mask: u32) -> usize {
        let low = mask & ((1 << self.half) - 1);
        self.low_ranks[low as usize] + self.high_ranks[(mask >> self.half) as usize]
    }

    pub fn unrank(&self, mut rank: usize) -> u32 {
        let mut mask = 0;
        let mut k = self.count;
        for n in (0..self.sites).rev() {
            if k == 0 {
                break;
            }
            let c = self.binomial(n, k);
            if rank >= c {
                mask |= 1 << n;
                rank -= c;
                k -= 1;
            }
        }
        mask
    }

    pub fn iter(&self) -> impl Iterator<Item = u32> {
        let first = ((1u64 << self.count) - 1) as u32;
        std::iter::successors(Some(first), |&mask| (mask != 0).then(|| next_combination(mask)))
            .take(self.len())
    }
}

/// Gosper's hack: the next larger integer with the same number of set bits.
fn next_combination(mask: u32) -> u32 {
    let mask = mask as u64;
    let lowest = mask.isolate_lowest_one();
    let ripple = mask + lowest;
    ((((ripple ^ mask) >> 2) / lowest) | ripple) as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_and_unrank_follow_increasing_order() {
        let c = Combinations::new(4, 2);
        let masks = [0b0011, 0b0101, 0b0110, 0b1001, 0b1010, 0b1100];
        assert_eq!(c.len(), 6);
        for (i, &m) in masks.iter().enumerate() {
            assert_eq!(c.rank(m), i);
            assert_eq!(c.unrank(i), m);
        }
        assert_eq!(c.iter().collect::<Vec<_>>(), masks);
    }

    #[test]
    fn iterator_agrees_with_unrank() {
        for (sites, count) in [(1, 0), (5, 0), (5, 5), (10, 3), (32, 2), (32, 31)] {
            let c = Combinations::new(sites, count);
            let mut n = 0;
            for (i, m) in c.iter().enumerate() {
                assert_eq!(m.count_ones(), count);
                assert_eq!(c.unrank(i), m);
                assert_eq!(c.rank(m), i);
                n += 1;
            }
            assert_eq!(n, c.len());
        }
    }
}
