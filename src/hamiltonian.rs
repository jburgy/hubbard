use rayon::prelude::*;

use crate::basis::{Basis, permute};

/// Positions of the set bits of `mask`, lowest first.
fn bits(mask: u32) -> impl Iterator<Item = usize> {
    let next = |&m: &u32| Some(m & (m - 1)).filter(|&m| m != 0);
    std::iter::successors(Some(mask).filter(|&m| m != 0), next).map(|m| m.trailing_zeros() as usize)
}

/// Sites strictly between `i` and `j`: their occupation sets the fermionic sign of a hop.
fn between(i: u8, j: u8) -> u32 {
    let (low, high) = (i.min(j), i.max(j));
    let below_high = (1u32 << high) - 1;
    let up_to_low = (1u32 << low << 1).wrapping_sub(1);
    below_high & !up_to_low
}

/// Fermions hopping from site `i` to `target[i]` along one lattice direction.
struct Direction {
    target: Vec<u8>,
    source: Vec<u8>,
    between: Vec<u32>,
}

impl Direction {
    fn new(target: &[u8]) -> Self {
        let mut source = vec![0; target.len()];
        for (i, &j) in target.iter().enumerate() {
            source[j as usize] = i as u8;
        }
        let between = target.iter().enumerate().map(|(i, &j)| between(i as u8, j)).collect();
        Direction { target: target.to_vec(), source, between }
    }

    /// Configurations reached by moving one fermion of `mask` onto an empty neighbour,
    /// and whether the move is odd.
    fn hops(&self, mask: u32) -> impl Iterator<Item = (u32, bool)> {
        let (blocked, _) = permute(&self.source, mask);
        bits(mask & !blocked).map(move |i| {
            let hopped = mask ^ (1 << i) ^ (1 << self.target[i]);
            (hopped, (mask & self.between[i]).count_ones() & 1 == 1)
        })
    }
}

/// `H = -t Σ_⟨ij⟩σ (c†_iσ c_jσ + h.c.) + U Σ_i n_i↑ n_i↓`, applied on the fly.
pub struct Hamiltonian {
    basis: Basis,
    directions: Vec<Direction>,
    hopping: f64,
    interaction: f64,
}

impl Hamiltonian {
    /// `neighbors[d][i]` is the neighbour of site `i` along direction `d`.
    pub fn new(basis: Basis, neighbors: &[Vec<u8>], hopping: f64, interaction: f64) -> Self {
        let directions = neighbors.iter().map(|target| Direction::new(target)).collect();
        Hamiltonian { basis, directions, hopping, interaction }
    }

    pub fn dimension(&self) -> usize {
        self.basis.dimension()
    }

    /// `y ← H x - β y`, row by row so that `y` can be overwritten in place.
    ///
    /// Each thread takes the contiguous states built on one spin up representative.
    pub fn apply(&self, x: &[f64], y: &mut [f64], beta: f64) {
        let mut blocks = Vec::with_capacity(self.basis.rep_count());
        let mut rest = y;
        for i in 0..self.basis.rep_count() {
            let (block, tail) = rest.split_at_mut(self.basis.rep(i).1.len());
            blocks.push((i, block));
            rest = tail;
        }
        blocks.into_par_iter().for_each(|(i, y)| {
            let (up, range) = self.basis.rep(i);
            let rows = range.zip(self.basis.rep_states(i)).zip(y);
            for ((a, (down, weight)), ya) in rows {
                *ya = self.row(a, up, down, weight, x) - beta * *ya;
            }
        });
    }

    fn row(&self, a: usize, up: u32, down: u32, weight: f64, x: &[f64]) -> f64 {
        let mut sum = self.interaction * (up & down).count_ones() as f64 * x[a];
        for direction in &self.directions {
            for (hopped, negative) in direction.hops(up) {
                sum += self.hop(hopped, down, negative, weight, x);
            }
            for (hopped, negative) in direction.hops(down) {
                sum += self.hop(up, hopped, negative, weight, x);
            }
        }
        sum
    }

    /// `⟨up, down|H|a⟩ x_b` where `b` is the basis state whose orbit contains `|up, down⟩`.
    fn hop(&self, up: u32, down: u32, negative: bool, weight: f64, x: &[f64]) -> f64 {
        let Some((b, symmetry_negative, b_weight)) = self.basis.find(up, down) else {
            return 0.0;
        };
        let t = if negative == symmetry_negative { -self.hopping } else { self.hopping };
        t * b_weight / weight * x[b]
    }

    /// `Σ_i ⟨n_i↑ n_i↓⟩` in the normalized state `psi`.
    pub fn double_occupancy(&self, psi: &[f64]) -> f64 {
        let occupancy = |(a, &p): (usize, &f64)| {
            let (up, down, _) = self.basis.state(a);
            p * p * (up & down).count_ones() as f64
        };
        psi.par_iter().enumerate().map(occupancy).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hop_sign_counts_fermions_in_between() {
        assert_eq!(between(4, 1), 0b01100);
        let direction = Direction::new(&[4, 0, 1, 2, 3]);
        let hops = |mask| direction.hops(mask).collect::<Vec<_>>();
        assert_eq!(hops(0b00001), [(0b10000, false)]);
        assert_eq!(hops(0b00011), [(0b10010, true)]);
        assert_eq!(hops(0b00111), [(0b10110, false)]);
        assert_eq!(hops(0b10100), [(0b10010, false), (0b01100, false)]);
        assert_eq!(hops(0b11111), []);
        assert_eq!(hops(0), []);
    }
}
