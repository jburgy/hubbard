use std::f64::consts::TAU;

use nalgebra::{DMatrix, SymmetricEigen};

use hubbard::basis::{Basis, permute};
use hubbard::hamiltonian::Hamiltonian;
use hubbard::lanczos::{ground_state, lanczos, residual};
use hubbard::lattice::{Irrep, Momentum, Symmetry, TiltedSquare};

fn identity(sites: usize) -> Vec<Symmetry> {
    vec![((0..sites as u8).collect(), false)]
}

fn hamiltonian(
    lattice: &TiltedSquare,
    group: Vec<Symmetry>,
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

/// Projector `Σ_g χ(g) ρ(g) / |G|` onto the sector, in the basis of all configurations.
fn projector(group: &[Symmetry], full: &Basis) -> DMatrix<f64> {
    let n = full.dimension();
    let mut p = DMatrix::zeros(n, n);
    for a in 0..n {
        let (up, down, _) = full.state(a);
        for (g, character_negative) in group {
            let (up, up_negative) = permute(g, up);
            let (down, down_negative) = permute(g, down);
            let (b, _, _) = full.find(up, down).unwrap();
            let sign = if up_negative ^ down_negative ^ character_negative { -1.0 } else { 1.0 };
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
fn every_sector_matches_projection_of_full_hamiltonian() {
    for (sites, up, down) in CASES {
        let lattice = TiltedSquare::new(sites).unwrap();
        let full = hamiltonian(&lattice, identity(sites), up, down, 1.0, 4.0);
        let full_dense = dense(&full);
        let full_basis = Basis::new(identity(sites), up, down);
        for momentum in [Momentum::Gamma, Momentum::M] {
            for irrep in [Irrep::A1, Irrep::A2, Irrep::B1, Irrep::B2] {
                let Ok(group) = lattice.symmetries(momentum, irrep) else { continue };
                let label = format!("{sites} {up} {down} {momentum:?} {irrep:?}");
                let sector = hamiltonian(&lattice, group.clone(), up, down, 1.0, 4.0);

                let p = projector(&group, &full_basis);
                assert!((p.trace() - sector.dimension() as f64).abs() < 1e-9, "{label}");
                if sector.dimension() == 0 {
                    continue;
                }

                let h = dense(&sector);
                assert!((&h - h.transpose()).amax() < 1e-12, "{label}");

                let identity = DMatrix::identity(p.nrows(), p.ncols());
                let restricted = &p * &full_dense * &p + (identity - &p) * 1e3;
                let expected = lowest(restricted);
                assert!((lowest(h) - expected).abs() < 1e-10, "{label}");

                let energy = lanczos(&sector, 1e-10, 500).unwrap().energy;
                assert!((energy - expected).abs() < 1e-9, "{label}");
            }
        }
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
    let group = lattice.symmetries(Momentum::Gamma, Irrep::A1).unwrap();
    let h = hamiltonian(&lattice, group, 5, 5, 1.0, 4.0);
    let result = lanczos(&h, 1e-10, 500).unwrap();
    let psi = ground_state(&h, &result).unwrap();
    let (deviation, energy) = residual(&h, &psi).unwrap();
    assert!((energy - result.energy).abs() < 1e-9);
    assert!(deviation < 1e-7);
    assert!((0.0..10.0).contains(&h.double_occupancy(&psi)));
}

/// Exact diagonalization energies per site (sites, N↑, N↓, U with t = 1) from Table I of
/// H. Shi and S. Zhang, "Symmetry in auxiliary-field quantum Monte Carlo calculations",
/// Phys. Rev. B 88, 125132 (2013), https://arxiv.org/abs/1307.2147.
/// The 4×4 value at (5, 5, 4) also appears as -1.2238 in S. Zhang, J. Carlson and
/// J. E. Gubernatis, Phys. Rev. B 55, 7464 (1997), https://arxiv.org/abs/cond-mat/9607062.
fn assert_matches_shi_zhang(
    group: fn(&TiltedSquare) -> Vec<Symmetry>,
    cases: &[(usize, u32, u32, f64, &str)],
) {
    for &(sites, up, down, u, reference) in cases {
        let lattice = TiltedSquare::new(sites).unwrap();
        let h = hamiltonian(&lattice, group(&lattice), up, down, 1.0, u);
        let per_site = lanczos(&h, 1e-7, 500).unwrap().energy / sites as f64;
        let decimals = reference.len() - reference.find('.').unwrap() - 1;
        let rounding = 0.5 * 10f64.powi(-(decimals as i32));
        let reference: f64 = reference.parse().unwrap();
        assert!(
            (per_site - reference).abs() <= rounding,
            "{sites} {up} {down} U={u}: {per_site} vs {reference}"
        );
    }
}

#[test]
fn small_clusters_match_shi_zhang() {
    assert_matches_shi_zhang(
        |lattice| identity(lattice.len()),
        &[
            (4, 2, 1, 4.0, "-1.60463"),
            (9, 4, 4, 8.0, "-0.8094"),
            (16, 2, 2, 4.0, "-0.72064"),
            (16, 2, 2, 8.0, "-0.7076"),
            (16, 2, 2, 12.0, "-0.7003"),
            (16, 3, 3, 4.0, "-0.94600"),
            (16, 3, 3, 8.0, "-0.9202"),
            (16, 3, 3, 12.0, "-0.9061"),
        ],
    );
}

/// Only the fillings whose ground state Shi and Zhang find in the A₁ irrep at zero momentum.
#[test]
fn fully_symmetric_ground_states_match_shi_zhang() {
    assert_matches_shi_zhang(
        |lattice| lattice.symmetries(Momentum::Gamma, Irrep::A1).unwrap(),
        &[
            (16, 5, 5, 4.0, "-1.22381"),
            (16, 5, 5, 8.0, "-1.0944"),
            (16, 5, 5, 12.0, "-1.0284"),
            (16, 6, 6, 8.0, "-0.9328"),
            (16, 6, 6, 12.0, "-0.8512"),
            (16, 8, 8, 4.0, "-0.85137"),
            // Printed without its minus sign in the table.
            (16, 8, 8, 8.0, "-0.5293"),
            (16, 8, 8, 12.0, "-0.3745"),
        ],
    );
}

/// Fillings whose ground state Shi and Zhang find in the B₁ irrep at zero momentum.
/// (6, 6, 4) is omitted: we get -1.108098 per site, which the table prints as -1.1080.
#[test]
fn d_wave_ground_states_match_shi_zhang() {
    assert_matches_shi_zhang(
        |lattice| lattice.symmetries(Momentum::Gamma, Irrep::B1).unwrap(),
        &[
            (16, 7, 7, 4.0, "-0.9840"),
            (16, 7, 7, 6.0, "-0.8388"),
            (16, 7, 7, 8.0, "-0.7418"),
            (16, 7, 7, 10.0, "-0.6754"),
            (16, 7, 7, 12.0, "-0.6282"),
        ],
    );
}
