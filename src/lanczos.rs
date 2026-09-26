use std::collections::TryReserveError;

use rayon::prelude::*;

use crate::hamiltonian::Hamiltonian;
use crate::tridiagonal::{lowest_eigenvalue, lowest_eigenvector};

/// Lanczos coefficients `α`, `β` and the lowest Ritz pair of the tridiagonal matrix they form.
pub struct Lanczos {
    pub alpha: Vec<f64>,
    pub beta: Vec<f64>,
    pub energy: f64,
    pub ritz: Vec<f64>,
}

impl Lanczos {
    pub fn iterations(&self) -> usize {
        self.alpha.len()
    }

    /// `‖H ψ - E ψ‖` for the Ritz vector `ψ`, without constructing it.
    pub fn residual(&self) -> f64 {
        self.beta.last().unwrap() * self.ritz.last().unwrap().abs()
    }
}

fn zeros(len: usize) -> Result<Vec<f64>, TryReserveError> {
    let mut v = Vec::new();
    v.try_reserve_exact(len)?;
    v.resize(len, 0.0);
    Ok(v)
}

/// Deterministic pseudo-random unit vector (splitmix64), so both passes start identically.
fn start_vector(len: usize) -> Result<Vec<f64>, TryReserveError> {
    let mut v = zeros(len)?;
    v.par_iter_mut().enumerate().for_each(|(i, x)| {
        let mut z = (i as u64).wrapping_add(1).wrapping_mul(0x9e3779b97f4a7c15);
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        *x = (z ^ (z >> 31)) as f64 / u64::MAX as f64 - 0.5;
    });
    let norm = dot(&v, &v).sqrt();
    scale(&mut v, 1.0 / norm);
    Ok(v)
}

fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.par_iter().zip(y).map(|(a, b)| a * b).sum()
}

fn scale(x: &mut [f64], s: f64) {
    x.par_iter_mut().for_each(|a| *a *= s);
}

/// `y ← y + s x`
fn add_scaled(y: &mut [f64], s: f64, x: &[f64]) {
    y.par_iter_mut().zip(x).for_each(|(a, b)| *a += s * b);
}

/// Plain Lanczos with two vectors: `w ← H v - β w` overwrites the previous basis vector.
///
/// Stops once the Ritz residual drops below `tolerance`.
pub fn lanczos(
    h: &Hamiltonian,
    tolerance: f64,
    max_iterations: usize,
) -> Result<Lanczos, TryReserveError> {
    let mut v = start_vector(h.dimension())?;
    let mut w = zeros(h.dimension())?;
    let mut result = Lanczos { alpha: Vec::new(), beta: Vec::new(), energy: 0.0, ritz: Vec::new() };
    loop {
        h.apply(&v, &mut w, result.beta.last().copied().unwrap_or(0.0));
        let a = dot(&v, &w);
        add_scaled(&mut w, -a, &v);
        let b = dot(&w, &w).sqrt();
        result.alpha.push(a);
        result.beta.push(b);

        let off_diagonal = &result.beta[..result.beta.len() - 1];
        result.energy = lowest_eigenvalue(&result.alpha, off_diagonal);
        result.ritz = lowest_eigenvector(&result.alpha, off_diagonal, result.energy);
        if result.residual() <= tolerance || result.iterations() == max_iterations {
            return Ok(result);
        }
        scale(&mut w, 1.0 / b);
        std::mem::swap(&mut v, &mut w);
    }
}

/// Replay the recurrence with the stored coefficients and accumulate `ψ = Σ_j s_j v_j`.
pub fn ground_state(h: &Hamiltonian, lanczos: &Lanczos) -> Result<Vec<f64>, TryReserveError> {
    let mut v = start_vector(h.dimension())?;
    let mut w = zeros(h.dimension())?;
    let mut psi = zeros(h.dimension())?;
    for j in 0..lanczos.iterations() {
        add_scaled(&mut psi, lanczos.ritz[j], &v);
        if j + 1 == lanczos.iterations() {
            break;
        }
        h.apply(&v, &mut w, if j == 0 { 0.0 } else { lanczos.beta[j - 1] });
        add_scaled(&mut w, -lanczos.alpha[j], &v);
        scale(&mut w, 1.0 / lanczos.beta[j]);
        std::mem::swap(&mut v, &mut w);
    }
    drop((v, w));
    let norm = dot(&psi, &psi).sqrt();
    scale(&mut psi, 1.0 / norm);
    Ok(psi)
}

/// `‖H ψ - ⟨ψ|H|ψ⟩ ψ‖` and `⟨ψ|H|ψ⟩`, using one extra vector.
pub fn residual(h: &Hamiltonian, psi: &[f64]) -> Result<(f64, f64), TryReserveError> {
    let mut h_psi = zeros(h.dimension())?;
    h.apply(psi, &mut h_psi, 0.0);
    let energy = dot(psi, &h_psi);
    add_scaled(&mut h_psi, -energy, psi);
    Ok((dot(&h_psi, &h_psi).sqrt(), energy))
}
