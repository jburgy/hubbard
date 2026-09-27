use std::f64::consts::TAU;

use nalgebra::{DMatrix, SymmetricEigen};

use hubbard::basis::{Basis, permute};
use hubbard::hamiltonian::Hamiltonian;
use hubbard::lanczos::{ground_state, lanczos, residual};
use hubbard::lattice::TiltedSquare;

fn identity(sites: usize) -> Vec<Vec<u8>> {
    vec![(0..sites as u8).collect()]
}

fn hamiltonian(
    lattice: &TiltedSquare,
    group: Vec<Vec<u8>>,
    up: u32,
    down: u32,
    t: f64,
    u: f64,
) -> Hamiltonian {
    Hamiltonian::new(Basis::new(group, up, down), &lattice.neighbors(), t, u)
}

fn dense(h: &Hamiltonian) -> DMatrix<f64> {
    let n = h.dimension();
    let mut matrix = DMatrix::zeros(n, n);
    let mut unit = vec![0.0; n];
    let mut column = vec![0.0; n];
    for j in 0..n {
        unit[j] = 1.0;
        h.apply(&unit, &mut column, 0.0);
        matrix.set_column(j, &column.clone().into());
        unit[j] = 0.0;
    }
    matrix
}

fn lowest(matrix: DMatrix<f64>) -> f64 {
    SymmetricEigen::new(matrix).eigenvalues.min()
}

/// Projector onto fully symmetric states, in the basis of all configurations.
fn projector(group: &[Vec<u8>], full: &Basis) -> DMatrix<f64> {
    let n = full.dimension();
    let mut p = DMatrix::zeros(n, n);
    for a in 0..n {
        let (up, down, _) = full.state(a);
        for g in group {
            let (up, up_negative) = permute(g, up);
            let (down, down_negative) = permute(g, down);
            let (b, _, _) = full.find(up, down).unwrap();
            let sign = if up_negative == down_negative { 1.0 } else { -1.0 };
            p[(b, a)] += sign / group.len() as f64;
        }
    }
    p
}

fn free_fermion_energy(lattice: &TiltedSquare, up: usize, down: usize, t: f64) -> f64 {
    let [u, v] = lattice.tilt();
    let n = lattice.len() as i32;
    let mut levels = Vec::new();
    for p in 0..n {
        for q in 0..n {
            if (p * u + q * v) % n == 0 && (q * u - p * v).rem_euclid(n) == 0 {
                let (kx, ky) = (TAU * p as f64 / n as f64, TAU * q as f64 / n as f64);
                levels.push(-2.0 * t * (kx.cos() + ky.cos()));
            }
        }
    }
    assert_eq!(levels.len(), lattice.len());
    levels.sort_by(f64::total_cmp);
    levels[..up].iter().sum::<f64>() + levels[..down].iter().sum::<f64>()
}

const CASES: [(usize, u32, u32); 8] =
    [(2, 1, 1), (4, 2, 1), (4, 2, 2), (5, 2, 2), (8, 2, 1), (8, 2, 2), (9, 2, 1), (10, 1, 2)];

#[test]
fn symmetric_sector_matches_projection_of_full_hamiltonian() {
    for (sites, up, down) in CASES {
        let lattice = TiltedSquare::new(sites).unwrap();
        let group = lattice.symmetries();
        let full = hamiltonian(&lattice, identity(sites), up, down, 1.0, 4.0);
        let symmetric = hamiltonian(&lattice, group.clone(), up, down, 1.0, 4.0);

        let p = projector(&group, &Basis::new(identity(sites), up, down));
        assert!((p.trace() - symmetric.dimension() as f64).abs() < 1e-9, "{sites} {up} {down}");

        let h = dense(&symmetric);
        assert!((&h - h.transpose()).amax() < 1e-12, "{sites} {up} {down}");

        let identity = DMatrix::identity(p.nrows(), p.ncols());
        let restricted = &p * dense(&full) * &p + (identity - &p) * 1e3;
        let expected = lowest(restricted);
        assert!((lowest(h) - expected).abs() < 1e-10, "{sites} {up} {down}");

        let energy = lanczos(&symmetric, 1e-10, 500).unwrap().energy;
        assert!((energy - expected).abs() < 1e-9, "{sites} {up} {down}");
    }
}

#[test]
fn full_basis_matches_dense_diagonalization() {
    for (sites, up, down) in CASES {
        let lattice = TiltedSquare::new(sites).unwrap();
        let h = hamiltonian(&lattice, identity(sites), up, down, 1.0, 4.0);
        let expected = lowest(dense(&h));
        let energy = lanczos(&h, 1e-10, 500).unwrap().energy;
        assert!((energy - expected).abs() < 1e-9, "{sites} {up} {down}");
    }
}

#[test]
fn noninteracting_fermions_fill_the_lowest_momenta() {
    for (sites, up, down) in
        [(4, 2, 1), (5, 2, 2), (8, 3, 2), (9, 2, 3), (10, 3, 3), (13, 2, 2), (16, 2, 1)]
    {
        let lattice = TiltedSquare::new(sites).unwrap();
        let h = hamiltonian(&lattice, identity(sites), up, down, 1.0, 0.0);
        let energy = lanczos(&h, 1e-10, 500).unwrap().energy;
        let expected = free_fermion_energy(&lattice, up as usize, down as usize, 1.0);
        assert!((energy - expected).abs() < 1e-9, "{sites} {up} {down}: {energy} vs {expected}");
    }
}

#[test]
fn atomic_limit_minimizes_double_occupancy() {
    let lattice = TiltedSquare::new(8).unwrap();
    for (up, down) in [(3, 3), (5, 4), (6, 6)] {
        let h = hamiltonian(&lattice, identity(8), up, down, 0.0, 2.0);
        let energy = lanczos(&h, 1e-12, 100).unwrap().energy;
        assert!((energy - 2.0 * (up + down).saturating_sub(8) as f64).abs() < 1e-12);
    }
}

#[test]
fn replayed_ground_state_is_an_eigenvector() {
    let lattice = TiltedSquare::new(10).unwrap();
    let h = hamiltonian(&lattice, lattice.symmetries(), 5, 5, 1.0, 4.0);
    let result = lanczos(&h, 1e-10, 500).unwrap();
    let psi = ground_state(&h, &result).unwrap();
    let (deviation, energy) = residual(&h, &psi).unwrap();
    assert!((energy - result.energy).abs() < 1e-9);
    assert!(deviation < 1e-7);
    assert!((0.0..10.0).contains(&h.double_occupancy(&psi)));
}
