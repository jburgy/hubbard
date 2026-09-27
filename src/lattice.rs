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

/// The eight elements of the dihedral group D₄: ρᵏ and ρᵏσ.
fn point_group() -> Vec<Matrix> {
    let mut group = Vec::with_capacity(8);
    for reflection in [IDENTITY, SIGMA] {
        let mut rotation = IDENTITY;
        for _ in 0..4 {
            group.push(multiply(rotation, reflection));
            rotation = multiply(RHO, rotation);
        }
    }
    group
}

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

    /// Distinct site permutations of the space group, identity first.
    pub fn symmetries(&self) -> Vec<Vec<u8>> {
        let mut group: Vec<Vec<u8>> = point_group()
            .into_iter()
            .filter(|&r| self.preserves_periodicity(r))
            .flat_map(|r| self.sites.iter().map(move |&d| self.permutation(r, d)))
            .collect();
        group.sort();
        group.dedup();
        group
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
            let group = TiltedSquare::new(n).unwrap().symmetries();
            assert_eq!(group.len(), order, "{n} sites");
            assert_eq!(group[0], (0..n as u8).collect::<Vec<_>>());
        }
    }

    #[test]
    fn symmetries_commute_with_hopping() {
        for n in [2, 4, 5, 8, 9, 10, 13, 16] {
            let lattice = TiltedSquare::new(n).unwrap();
            let neighbors = lattice.neighbors();
            for g in lattice.symmetries() {
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
