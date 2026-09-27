type Point = [i32; 2];
type Matrix = [[i32; 2]; 2];

const IDENTITY: Matrix = [[1, 0], [0, 1]];
/// Rotation by ½π.
const RHO: Matrix = [[0, -1], [1, 0]];
/// Reflection about the `x` axis.
const SIGMA: Matrix = [[-1, 0], [0, 1]];

fn apply(m: Matrix, x: Point) -> Point {
    [m[0][0] * x[0] + m[0][1] * x[1], m[1][0] * x[0] + m[1][1] * x[1]]
}

fn multiply(a: Matrix, b: Matrix) -> Matrix {
    let column = |j: usize| apply(a, [b[0][j], b[1][j]]);
    let (c0, c1) = (column(0), column(1));
    [[c0[0], c1[0]], [c0[1], c1[1]]]
}

/// The eight elements ρᵏσᵐ of the dihedral group D₄, with the parities of `k` and `m`.
fn point_group() -> Vec<(Matrix, bool, bool)> {
    let mut group = Vec::with_capacity(8);
    for (reflection, reflected) in [(IDENTITY, false), (SIGMA, true)] {
        let mut rotation = IDENTITY;
        for k in 0..4 {
            group.push((multiply(rotation, reflection), k % 2 == 1, reflected));
            rotation = multiply(RHO, rotation);
        }
    }
    group
}

/// Crystal momenta whose Bloch phases are all real.
#[derive(Clone, Copy, Debug, clap::ValueEnum)]
pub enum Momentum {
    Gamma,
    M,
}

/// One-dimensional irreps of C₄ᵥ, labelled as for d orbitals: B₁ is x² - y², B₂ is xy.
#[derive(Clone, Copy, Debug, clap::ValueEnum)]
pub enum Irrep {
    A1,
    A2,
    B1,
    B2,
}

impl Irrep {
    /// Whether the character of a quarter turn, and of the reflection `SIGMA`, is -1.
    fn odd(self) -> (bool, bool) {
        match self {
            Irrep::A1 => (false, false),
            Irrep::A2 => (false, true),
            Irrep::B1 => (true, false),
            Irrep::B2 => (true, true),
        }
    }
}

/// A site permutation and whether its character is -1.
pub type Symmetry = (Vec<u8>, bool);

/// Square lattice with periodic boundary conditions along `tilt` and `ρ tilt`.
///
/// Sites are sorted lexicographically, so the eight site cluster is labelled
///
/// ```text
///       5
///     1 4 7
///     0 3 6
///       2
/// ```
pub struct TiltedSquare {
    tilt: Point,
    sites: Vec<Point>,
}

impl TiltedSquare {
    /// Fails unless `n` is a positive sum of two squares.
    pub fn new(n: usize) -> Option<Self> {
        let n = i32::try_from(n).ok().filter(|&n| n > 0)?;
        let v = (0..=n.isqrt()).find(|v| (n - v * v).isqrt().pow(2) + v * v == n)?;
        let u = (n - v * v).isqrt();
        let mut lattice = TiltedSquare { tilt: [u, v], sites: Vec::new() };
        for x in -v..=u {
            for y in 0..=u + v {
                if lattice.restrict([x, y]) == [x, y] {
                    lattice.sites.push([x, y]);
                }
            }
        }
        Some(lattice)
    }

    pub fn len(&self) -> usize {
        self.sites.len()
    }

    pub fn tilt(&self) -> Point {
        self.tilt
    }

    /// Fold `x` into the cluster `0 ≤ x·tilt < n`, `0 ≤ x·(ρ tilt) < n`.
    fn restrict(&self, mut x: Point) -> Point {
        let [u, v] = self.tilt;
        let n = u * u + v * v;
        for w in [[u, v], [-v, u]] {
            let k = (x[0] * w[0] + x[1] * w[1]).div_euclid(n);
            x = [x[0] - k * w[0], x[1] - k * w[1]];
        }
        x
    }

    fn index(&self, x: Point) -> u8 {
        self.sites.binary_search(&self.restrict(x)).unwrap() as u8
    }

    /// Site permutation induced by `x ↦ r x + d`.
    pub fn permutation(&self, r: Matrix, d: Point) -> Vec<u8> {
        self.sites
            .iter()
            .map(|&s| {
                let [x, y] = apply(r, s);
                self.index([x + d[0], y + d[1]])
            })
            .collect()
    }

    /// Nearest neighbour of every site along `+x`, `-x`, `+y` and `-y`.
    pub fn neighbors(&self) -> Vec<Vec<u8>> {
        [[1, 0], [-1, 0], [0, 1], [0, -1]]
            .into_iter()
            .map(|d| self.permutation(IDENTITY, d))
            .collect()
    }

    /// A point group element is a symmetry only if it maps the superlattice onto itself.
    fn preserves_periodicity(&self, r: Matrix) -> bool {
        let [u, v] = self.tilt;
        [[u, v], [-v, u]].into_iter().all(|w| self.restrict(apply(r, w)) == [0, 0])
    }

    /// Chiral clusters have no reflections, so A₂ coincides with A₁ and B₂ with B₁.
    pub fn is_chiral(&self) -> bool {
        !self.preserves_periodicity(SIGMA)
    }

    /// Distinct space group elements `x ↦ r x + d`, identity first, with the characters of
    /// the irrep `irrep` at crystal momentum `momentum`.
    pub fn symmetries(&self, momentum: Momentum, irrep: Irrep) -> Result<Vec<Symmetry>, String> {
        let [u, v] = self.tilt;
        let staggered = matches!(momentum, Momentum::M);
        if staggered && (u + v) % 2 == 1 {
            return Err("momentum (π, π) does not fit this cluster".into());
        }
        let (odd_turn, odd_reflection) = irrep.odd();
        let mut group = Vec::new();
        for (r, turned, reflected) in point_group() {
            if !self.preserves_periodicity(r) {
                continue;
            }
            let point_negative = (turned && odd_turn) ^ (reflected && odd_reflection);
            for &d in &self.sites {
                let phase_negative = staggered && (d[0] + d[1]) % 2 != 0;
                group.push((self.permutation(r, d), point_negative ^ phase_negative));
            }
        }
        group.sort();
        group.dedup();
        if group.windows(2).any(|pair| pair[0].0 == pair[1].0) {
            return Err(format!(
                "{irrep:?} at {momentum:?} is not a representation on this cluster"
            ));
        }
        Ok(group)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eight_sites() {
        let lattice = TiltedSquare::new(8).unwrap();
        assert_eq!(lattice.tilt(), [2, 2]);
        assert_eq!(lattice.permutation(RHO, [0, 0]), [6, 3, 2, 7, 4, 1, 0, 5]);
        assert_eq!(lattice.permutation(IDENTITY, [1, 0]), [3, 4, 1, 6, 7, 0, 5, 2]);
        assert_eq!(lattice.neighbors()[1], [5, 2, 7, 0, 1, 6, 3, 4]);
    }

    #[test]
    fn only_sums_of_two_squares() {
        let sizes: Vec<usize> = (0..=20).filter(|&n| TiltedSquare::new(n).is_some()).collect();
        assert_eq!(sizes, [1, 2, 4, 5, 8, 9, 10, 13, 16, 17, 18, 20]);
    }

    #[test]
    fn chiral_clusters_lack_reflections() {
        for (n, order) in [(5, 20), (8, 64), (9, 72), (10, 40), (16, 128), (18, 144), (20, 80)] {
            let lattice = TiltedSquare::new(n).unwrap();
            let group = lattice.symmetries(Momentum::Gamma, Irrep::A1).unwrap();
            assert_eq!(group.len(), order, "{n} sites");
            assert_eq!(group[0], ((0..n as u8).collect(), false));
            assert_eq!(lattice.is_chiral(), order == 4 * n);
        }
    }

    #[test]
    fn staggered_momentum_needs_an_even_superlattice() {
        for (n, fits) in [(5, false), (8, true), (9, false), (10, true), (16, true)] {
            let lattice = TiltedSquare::new(n).unwrap();
            assert_eq!(lattice.symmetries(Momentum::M, Irrep::B1).is_ok(), fits, "{n} sites");
        }
    }

    #[test]
    fn characters_multiply() {
        let lattice = TiltedSquare::new(16).unwrap();
        let compose =
            |a: &[u8], b: &[u8]| -> Vec<u8> { b.iter().map(|&j| a[j as usize]).collect() };
        for irrep in [Irrep::A1, Irrep::A2, Irrep::B1, Irrep::B2] {
            let group = lattice.symmetries(Momentum::M, irrep).unwrap();
            for (g, g_negative) in &group {
                for (h, h_negative) in &group {
                    let product = compose(g, h);
                    let (_, negative) = group.iter().find(|(p, _)| *p == product).unwrap();
                    assert_eq!(*negative, g_negative ^ h_negative, "{irrep:?}");
                }
            }
        }
    }

    #[test]
    fn symmetries_commute_with_hopping() {
        for n in [2, 4, 5, 8, 9, 10, 13, 16] {
            let lattice = TiltedSquare::new(n).unwrap();
            let neighbors = lattice.neighbors();
            for (g, _) in lattice.symmetries(Momentum::Gamma, Irrep::A1).unwrap() {
                let compose =
                    |a: &[u8], b: &[u8]| -> Vec<u8> { b.iter().map(|&j| a[j as usize]).collect() };
                let mut before: Vec<_> = neighbors.iter().map(|d| compose(d, &g)).collect();
                let mut after: Vec<_> = neighbors.iter().map(|d| compose(&g, d)).collect();
                before.sort();
                after.sort();
                assert_eq!(before, after, "{n} sites");
            }
        }
    }
}
