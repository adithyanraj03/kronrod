//! Gauss-Legendre rules from scratch, via two independent methods.
//!
//! * [`gauss_legendre_newton`] — Newton iteration on the three-term Legendre
//!   recurrence, seeded by `x_i = cos((i - 1/2)·π/(n + 1/2))`, `i = 1..n`.
//!   This seed never lands exactly on 0: the textbook seed
//!   `cos((2i-1)π/(2n+2))` is exactly 0 for even `n`, where `P'_n(0) = 0`
//!   and Newton divides by zero.
//! * [`gauss_legendre_golub_welsch`] — the Golub-Welsch theorem
//!   (Math. Comp. 23 (1969) 221-230): for the uniform weight on [-1, 1] the
//!   Jacobi matrix is tridiagonal with `α_k = 0` and
//!   `β_k = k²/((2k-1)(2k+1))`; the nodes are its eigenvalues and the
//!   weights are `w_i = 2·v_i[0]²`. The eigendecomposition is plain cyclic
//!   Jacobi rotation (no external BLAS).
//!
//! Both methods agree with each other and with the embedded official
//! tables ([`crate::tables`]) to full `f64` precision (nodes ~1e-16,
//! weights ~1e-15 for n ≤ 30); see `tests/gauss.rs`.

use std::f64::consts::PI;

/// Gauss-Legendre nodes and weights on [-1, 1] via Newton iteration.
///
/// Returns `(node, weight)` pairs sorted ascending. `n` must be ≥ 1.
pub fn gauss_legendre_newton(n: usize) -> Vec<(f64, f64)> {
    assert!(n >= 1, "n >= 1");
    let mut nodes: Vec<f64> = Vec::with_capacity(n);
    for i in 1..=n {
        let mut x = ((i as f64 - 0.5) * PI / (n as f64 + 0.5)).cos();
        for _ in 0..300 {
            // P_n(x) via the three-term recurrence.
            let (_p0, p1) = legendre_pair(x, n);
            let pn = p1;
            let pnm1 = if n == 1 { 1.0 } else { legendre_pair(x, n - 1).1 };
            // P'_n(x) = n (P_{n-1}(x) - x P_n(x)) / (1 - x^2)
            let d_pn = n as f64 * (pnm1 - x * pn) / (1.0 - x * x);
            let step = pn / d_pn;
            x -= step;
            if step.abs() < 1e-17 {
                break;
            }
        }
        nodes.push(x);
    }
    nodes.sort_by(f64::total_cmp);
    let weights = nodes
        .iter()
        .map(|&x| {
            let (p0, p1) = legendre_pair(x, n);
            let _ = p0;
            let pn = p1;
            let pnm1 = if n == 1 { 1.0 } else { legendre_pair(x, n - 1).1 };
            let d_pn = n as f64 * (pnm1 - x * pn) / (1.0 - x * x);
            2.0 / ((1.0 - x * x) * d_pn * d_pn)
        })
        .collect::<Vec<_>>();
    nodes.into_iter().zip(weights).collect()
}

/// (P_{n-1}(x), P_n(x)) via the three-term recurrence.
fn legendre_pair(x: f64, n: usize) -> (f64, f64) {
    if n == 0 {
        return (x, 1.0);
    }
    let (mut p0, mut p1) = (1.0, x);
    for k in 1..n {
        let p2 = ((2.0 * k as f64 + 1.0) * x * p1 - k as f64 * p0) / (k as f64 + 1.0);
        p0 = p1;
        p1 = p2;
    }
    (p0, p1)
}

/// Gauss-Legendre nodes and weights on [-1, 1] via the Golub-Welsch
/// Jacobi-matrix eigendecomposition. Returns `(node, weight)` pairs sorted
/// ascending. `n` must be ≥ 1.
pub fn gauss_legendre_golub_welsch(n: usize) -> Vec<(f64, f64)> {
    assert!(n >= 1, "n >= 1");
    // Jacobi matrix: zero diagonal, off-diagonals sqrt(beta_k).
    let mut j = vec![vec![0.0f64; n]; n];
    for k in 1..n {
        let beta = (k as f64 * k as f64) / ((2.0 * k as f64 - 1.0) * (2.0 * k as f64 + 1.0));
        let off = beta.sqrt();
        j[k - 1][k] = off;
        j[k][k - 1] = off;
    }
    let (eig, v) = jacobi_eig(j);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| eig[a].total_cmp(&eig[b]));
    let mut out = Vec::with_capacity(n);
    for &i in order.iter() {
        // Eigenvector i is COLUMN i of v (v[row][col]); its first
        // component is v[0][i]. Golub-Welsch: w_i = 2·v[0][i]².
        let v0 = v[0][i];
        out.push((eig[i], 2.0 * v0 * v0));
    }
    out
}

/// Cyclic Jacobi eigendecomposition of a small symmetric matrix.
///
/// Returns `(eigenvalues, eigenvectors)` where eigenvector `i` is column
/// `i` of the returned matrix (`eigenvectors[row][col]`), unsorted; the
/// caller sorts. Columns are orthonormal.
fn jacobi_eig(input: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = input.len();
    let mut a = input;
    let mut v: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect()).collect();

    // Givens rotation that zeroes a[i][j] when tan(2θ) = 2a[i][j]/(a[i][i]-a[j][j]):
    //   v[:,i] <- c·v[:,i] + s·v[:,j];  v[:,j] <- -s·v[:,i] + c·v[:,j]
    //   a ← G^T a G (applied to columns, then rows).
    let rot = |i: usize, j: usize, theta: f64, a: &mut Vec<Vec<f64>>, v: &mut Vec<Vec<f64>>| {
        let c = theta.cos();
        let s = theta.sin();
        for k in 0..n {
            let (vi, vj) = (v[k][i], v[k][j]);
            v[k][i] = c * vi + s * vj;
            v[k][j] = -s * vi + c * vj;
        }
        for k in 0..n {
            let (ai, aj) = (a[k][i], a[k][j]);
            a[k][i] = c * ai + s * aj;
            a[k][j] = -s * ai + c * aj;
        }
        for k in 0..n {
            let (ai, aj) = (a[i][k], a[j][k]);
            a[i][k] = c * ai + s * aj;
            a[j][k] = -s * ai + c * aj;
        }
    };

    for _sweep in 0..200 {
        let mut off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off += a[i][j] * a[i][j];
            }
        }
        if off.sqrt() < 1e-16 * (1.0 + a[0][0].abs()) {
            break;
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if a[i][j].abs() < 1e-300 {
                    continue;
                }
                let theta = 0.5 * (2.0 * a[i][j]).atan2(a[i][i] - a[j][j]);
                rot(i, j, theta, &mut a, &mut v);
            }
        }
    }

    let eig = (0..n).map(|i| a[i][i]).collect();
    (eig, v)
}
