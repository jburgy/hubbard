use std::ops::Range;

use crate::combinations::Combinations;

/// Image of the occupation `mask` under the site permutation `image`, and whether reordering
/// the creation operators into increasing site order flips the sign.
pub fn permute(image: &[u8], mut mask: u32) -> (u32, bool) {
    let mut result = 0u32;
    let mut inversions = 0;
    while mask != 0 {
        let j = image[mask.trailing_zeros() as usize];
        inversions += (result >> j).count_ones();
        result |= 1 << j;
        mask &= mask - 1;
    }
    (result, inversions & 1 == 1)
}

/// How a spin up configuration maps onto its orbit representative.
struct Orbit {
    rep: u32,
    element: u16,
    negative: bool,
}

/// A spin up orbit representative and the spin down configurations paired with it.
///
/// When only the identity fixes `up`, every spin down configuration is allowed and
/// its rank is its position. Otherwise the allowed ones are listed explicitly, and
/// `lookup[rank(down)]` holds `±(position + 1)` of the listed configuration in the same
/// orbit, signed like the symmetry relating them, or 0 when the symmetric state vanishes.
struct Rep {
    up: u32,
    offset: usize,
    len: usize,
    stabilizer: Vec<(u16, bool)>,
    downs: Vec<u32>,
    weights: Vec<f64>,
    lookup: Vec<i32>,
}

impl Rep {
    fn is_dense(&self) -> bool {
        self.stabilizer.len() == 1
    }
}

/// Fully symmetric states `Σ_g g|up, down⟩ / √(|G| |stabilizer|)`, one per orbit.
///
/// Symmetries act on the spin up configuration first, so auxiliary storage scales
/// with the number of spin up configurations rather than with the dimension.
pub struct Basis {
    group: Vec<Vec<u8>>,
    up: Combinations,
    down: Combinations,
    orbits: Vec<Orbit>,
    reps: Vec<Rep>,
    dimension: usize,
}

impl Basis {
    pub fn new(group: Vec<Vec<u8>>, up: u32, down: u32) -> Self {
        let sites = group[0].len() as u32;
        let mut basis = Basis {
            group,
            up: Combinations::new(sites, up),
            down: Combinations::new(sites, down),
            orbits: Vec::new(),
            reps: Vec::new(),
            dimension: 0,
        };
        basis.orbits.reserve_exact(basis.up.len());
        for u in basis.up.iter() {
            let (image, element, negative) = basis.smallest_image(u);
            let orbit = if image == u {
                let rep = basis.new_rep(u);
                basis.dimension += rep.len;
                basis.reps.push(rep);
                Orbit { rep: basis.reps.len() as u32 - 1, element: 0, negative: false }
            } else {
                let rep = basis.orbits[basis.up.rank(image)].rep;
                Orbit { rep, element, negative }
            };
            basis.orbits.push(orbit);
        }
        basis
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn group_order(&self) -> usize {
        self.group.len()
    }

    pub fn rep_count(&self) -> usize {
        self.reps.len()
    }

    /// Spin up configuration and index range of the states built on representative `i`.
    pub fn rep(&self, i: usize) -> (u32, Range<usize>) {
        let rep = &self.reps[i];
        (rep.up, rep.offset..rep.offset + rep.len)
    }

    /// Spin down configurations and `√|stabilizer|` of the states built on representative `i`.
    pub fn rep_states(&self, i: usize) -> impl Iterator<Item = (u32, f64)> {
        let rep = &self.reps[i];
        // Exactly one of the two iterators is non-empty.
        let dense = rep.is_dense().then(|| self.down.iter()).into_iter().flatten();
        let listed = rep.downs.iter().copied().zip(rep.weights.iter().copied());
        dense.map(|down| (down, 1.0)).chain(listed)
    }

    fn smallest_image(&self, up: u32) -> (u32, u16, bool) {
        let images = self.group.iter().enumerate().map(|(i, g)| {
            let (image, negative) = permute(g, up);
            (image, i as u16, negative)
        });
        images.min_by_key(|&(image, _, _)| image).unwrap()
    }

    fn new_rep(&self, up: u32) -> Rep {
        let mut stabilizer = Vec::new();
        for (i, g) in self.group.iter().enumerate() {
            let (image, negative) = permute(g, up);
            if image == up {
                stabilizer.push((i as u16, negative));
            }
        }
        let mut rep = Rep {
            up,
            offset: self.dimension,
            len: 0,
            stabilizer,
            downs: vec![],
            weights: vec![],
            lookup: vec![],
        };
        if rep.is_dense() {
            rep.len = self.down.len();
            return rep;
        }
        for d in self.down.iter() {
            if let Some(size) = self.stabilizer_size(&rep, d) {
                rep.downs.push(d);
                rep.weights.push((size as f64).sqrt());
            }
        }
        rep.len = rep.downs.len();
        rep.lookup = self.down.iter().map(|d| self.lookup_entry(&rep, d)).collect();
        rep
    }

    fn lookup_entry(&self, rep: &Rep, down: u32) -> i32 {
        let images = rep.stabilizer.iter().map(|&(g, up_negative)| {
            let (image, down_negative) = permute(&self.group[g as usize], down);
            (image, up_negative ^ down_negative)
        });
        let (image, negative) = images.min_by_key(|&(image, _)| image).unwrap();
        match rep.downs.binary_search(&image) {
            Ok(local) if negative => -(local as i32 + 1),
            Ok(local) => local as i32 + 1,
            Err(_) => 0,
        }
    }

    /// `None` unless `down` is the smallest image under the stabilizer of `rep.up` and
    /// the symmetric state does not vanish.
    fn stabilizer_size(&self, rep: &Rep, down: u32) -> Option<usize> {
        let mut size = 0;
        for &(g, up_negative) in &rep.stabilizer {
            let (image, negative) = permute(&self.group[g as usize], down);
            if image < down || (image == down && negative != up_negative) {
                return None;
            }
            size += (image == down) as usize;
        }
        Some(size)
    }

    /// Spin up mask, spin down mask and `√|stabilizer|` of the representative of state `index`.
    pub fn state(&self, index: usize) -> (u32, u32, f64) {
        let rep = &self.reps[self.reps.partition_point(|rep| rep.offset <= index) - 1];
        let local = index - rep.offset;
        if rep.is_dense() {
            (rep.up, self.down.unrank(local), 1.0)
        } else {
            (rep.up, rep.downs[local], rep.weights[local])
        }
    }

    /// Index, sign and `√|stabilizer|` of the basis state whose orbit contains `|up, down⟩`.
    pub fn find(&self, up: u32, down: u32) -> Option<(usize, bool, f64)> {
        let orbit = &self.orbits[self.up.rank(up)];
        let (down, down_negative) = if orbit.element == 0 {
            (down, false)
        } else {
            permute(&self.group[orbit.element as usize], down)
        };
        let negative = orbit.negative ^ down_negative;
        let rep = &self.reps[orbit.rep as usize];
        if rep.is_dense() {
            return Some((rep.offset + self.down.rank(down), negative, 1.0));
        }
        let entry = rep.lookup[self.down.rank(down)];
        if entry == 0 {
            return None;
        }
        let local = entry.unsigned_abs() as usize - 1;
        Some((rep.offset + local, negative ^ (entry < 0), rep.weights[local]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn permutation_sign_counts_inversions() {
        let swap = [1, 0, 2];
        assert_eq!(permute(&swap, 0b011), (0b011, true));
        assert_eq!(permute(&swap, 0b101), (0b110, false));
        let cycle = [1, 2, 0];
        assert_eq!(permute(&cycle, 0b111), (0b111, false));
        assert_eq!(permute(&cycle, 0b101), (0b011, true));
    }

    #[test]
    fn identity_group_spans_every_configuration() {
        let basis = Basis::new(vec![(0..6).collect()], 2, 3);
        assert_eq!(basis.dimension(), 15 * 20);
        for index in 0..basis.dimension() {
            let (up, down, weight) = basis.state(index);
            assert_eq!(basis.find(up, down), Some((index, false, weight)));
        }
    }
}
