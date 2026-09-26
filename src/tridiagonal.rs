//! Symmetric tridiagonal matrices with diagonal `a` and off-diagonal `b` (`b.len() + 1 == a.len()`).

/// Number of eigenvalues below `x`, from the signs of the Sturm sequence.
fn count_below(a: &[f64], b: &[f64], x: f64) -> usize {
    let mut count = 0;
    let mut q = 1.0;
    for i in 0..a.len() {
        let coupling = if i == 0 { 0.0 } else { b[i - 1] * b[i - 1] };
        q = a[i] - x - coupling / q;
        if q == 0.0 {
            q = f64::MIN_POSITIVE;
        }
        count += (q < 0.0) as usize;
    }
    count
}

/// Bisection between the Gershgorin bounds.
pub fn lowest_eigenvalue(a: &[f64], b: &[f64]) -> f64 {
    let radius = |i: usize| {
        let left = if i == 0 { 0.0 } else { b[i - 1].abs() };
        let right = if i == b.len() { 0.0 } else { b[i].abs() };
        left + right
    };
    let mut lo = (0..a.len()).map(|i| a[i] - radius(i)).fold(f64::INFINITY, f64::min);
    let mut hi = (0..a.len()).map(|i| a[i] + radius(i)).fold(f64::NEG_INFINITY, f64::max);
    loop {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            return mid;
        }
        if count_below(a, b, mid) > 0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
}

/// Normalized eigenvector for the lowest `eigenvalue`, by inverse iteration.
///
/// Shifting just below the bottom of the spectrum keeps `T - σ` positive definite,
/// so elimination without pivoting is stable.
pub fn lowest_eigenvector(a: &[f64], b: &[f64], eigenvalue: f64) -> Vec<f64> {
    let scale = a.iter().chain(b).fold(1.0f64, |m, x| m.max(x.abs()));
    let shift = eigenvalue - 1e-10 * scale;
    let mut x = vec![1.0; a.len()];
    for _ in 0..3 {
        x = solve_shifted(a, b, shift, x);
        let norm = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        x.iter_mut().for_each(|v| *v /= norm);
    }
    x
}

/// Solve `(T - σ) x = rhs` with the Thomas algorithm.
fn solve_shifted(a: &[f64], b: &[f64], shift: f64, mut x: Vec<f64>) -> Vec<f64> {
    let n = a.len();
    let mut pivots = vec![a[0] - shift; n];
    for i in 1..n {
        let l = b[i - 1] / pivots[i - 1];
        pivots[i] = a[i] - shift - l * b[i - 1];
        x[i] -= l * x[i - 1];
    }
    x[n - 1] /= pivots[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (x[i] - b[i] * x[i + 1]) / pivots[i];
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn free_particle_on_a_chain() {
        let n = 50;
        let (a, b) = (vec![2.0; n], vec![-1.0; n - 1]);
        let k = std::f64::consts::PI / (n + 1) as f64;
        let eigenvalue = lowest_eigenvalue(&a, &b);
        assert!((eigenvalue - (2.0 - 2.0 * k.cos())).abs() < 1e-14);

        let x = lowest_eigenvector(&a, &b, eigenvalue);
        let norm = (2.0 / (n + 1) as f64).sqrt();
        for (i, v) in x.iter().enumerate() {
            assert!((v.abs() - norm * (k * (i + 1) as f64).sin()).abs() < 1e-10);
        }
    }

    #[test]
    fn single_element() {
        assert_eq!(lowest_eigenvalue(&[3.0], &[]), 3.0);
        assert_eq!(lowest_eigenvector(&[3.0], &[], 3.0), [1.0]);
    }
}
