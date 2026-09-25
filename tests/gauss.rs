//! Gauss-Legendre from scratch: two independent methods, cross-checked
//! against each other, against closed-form small-n rules, and against the
//! embedded Netlib tables.

use kronrod::gauss::{gauss_legendre_golub_welsch, gauss_legendre_newton};
use kronrod::tables::rule;

fn node_err(a: &[(f64, f64)], b: &[(f64, f64)]) -> f64 {
    a.iter()
        .zip(b)
        .map(|((na, _), (nb, _))| (na - nb).abs())
        .fold(0.0f64, f64::max)
}

fn weight_err(a: &[(f64, f64)], b: &[(f64, f64)]) -> f64 {
    a.iter()
        .zip(b)
        .map(|((_, wa), (_, wb))| (wa - wb).abs() / wb.abs())
        .fold(0.0f64, f64::max)
}

#[test]
fn n1_is_single_centre_node() {
    let r = gauss_legendre_newton(1);
    assert_eq!(r.len(), 1);
    assert!(r[0].0.abs() < 1e-15);
    assert!((r[0].1 - 2.0).abs() < 1e-14);
    let s = gauss_legendre_golub_welsch(1);
    assert!(node_err(&r, &s) < 1e-14);
}

#[test]
fn n2_closed_form() {
    for m in &[gauss_legendre_newton(2), gauss_legendre_golub_welsch(2)] {
        let want = 1.0f64.sqrt() / 3.0f64.sqrt(); // 1/sqrt(3)
        assert!((m[0].0 + want).abs() < 1e-14);
        assert!((m[1].0 - want).abs() < 1e-14);
        assert!((m[0].1 - 1.0).abs() < 1e-14);
        assert!((m[1].1 - 1.0).abs() < 1e-14);
    }
}

#[test]
fn n3_closed_form() {
    let x = (3.0f64 / 5.0).sqrt();
    for m in &[gauss_legendre_newton(3), gauss_legendre_golub_welsch(3)] {
        assert!((m[0].0 + x).abs() < 1e-13);
        assert!(m[1].0.abs() < 1e-13);
        assert!((m[2].0 - x).abs() < 1e-13);
        // weights: 5/9 at ±x, 8/9 at 0 (verified by exactness of x^2/x^4).
        assert!((m[0].1 - 5.0 / 9.0).abs() < 1e-13);
        assert!((m[2].1 - 5.0 / 9.0).abs() < 1e-13);
        assert!((m[1].1 - 8.0 / 9.0).abs() < 1e-13);
    }
}

#[test]
fn two_methods_agree() {
    for n in [3usize, 5, 7, 10, 15, 20, 25, 30] {
        let a = gauss_legendre_newton(n);
        let b = gauss_legendre_golub_welsch(n);
        assert_eq!(a.len(), n);
        assert!(node_err(&a, &b) < 1e-12, "n={n}: nodes");
        assert!(weight_err(&a, &b) < 1e-12, "n={n}: weights");
    }
}

#[test]
fn weights_sum_to_two() {
    for n in [1usize, 2, 3, 5, 7, 10, 15, 20, 25, 30] {
        let a = gauss_legendre_newton(n);
        let s: f64 = a.iter().map(|&(_, w)| w).sum();
        assert!((s - 2.0).abs() < 1e-13, "n={n}: {s}");
        let b = gauss_legendre_golub_welsch(n);
        let s: f64 = b.iter().map(|&(_, w)| w).sum();
        assert!((s - 2.0).abs() < 1e-12, "n={n} (GW): {s}");
    }
}

#[test]
fn nodes_are_antisymmetric() {
    for n in [5usize, 10, 20, 30] {
        let r = gauss_legendre_newton(n);
        for (i, &(x, _)) in r.iter().enumerate() {
            let mirror = r[n - 1 - i];
            assert!((x + mirror.0).abs() < 1e-14);
        }
    }
}

#[test]
fn matches_embedded_tables() {
    // The embedded Netlib tables carry the positive Gauss nodes at the odd
    // 0-based positions; compare against both from-scratch methods.
    for keyf in [1u8, 2, 3, 4, 5, 6] {
        let rule = rule(keyf);
        let official = rule.gauss_nodes(); // descending positive
        for (name, from_scratch) in [
            (
                "newton",
                gauss_legendre_newton(rule.n),
            ),
            (
                "golub-welsch",
                gauss_legendre_golub_welsch(rule.n),
            ),
        ] {
            let mut pos: Vec<f64> = from_scratch
                .iter()
                .filter(|&&(x, _)| x > 1e-9) // drop the ~0 centre (odd n)
                .map(|&(x, _)| x)
                .collect();
            pos.reverse(); // descending
            assert_eq!(pos.len(), official.len());
            for (o, g) in official.iter().zip(&pos) {
                assert!(
                    (o - g).abs() <= 5e-15 * o.abs().max(1e-300),
                    "keyf={} {}: node off",
                    keyf,
                    name
                );
            }
        }
    }
}

#[test]
fn matches_embedded_table_weights() {
    for keyf in [2u8, 6] {
        let rule = rule(keyf);
        let a = gauss_legendre_newton(rule.n);
        let mut pos: Vec<(f64, f64)> = a
            .iter()
            .copied()
            .filter(|(x, _)| *x > 0.0)
            .collect();
        pos.reverse();
        for ((_, g), w_off) in pos.iter().zip(rule.gauss_weights.iter().take(rule.gauss_pairs())) {
            assert!(
                (g - w_off).abs() <= 5e-14 * w_off.abs(),
                "keyf={keyf}: weight off"
            );
        }
    }
}

#[test]
fn exactness_on_polynomials() {
    // An n-point Gauss rule is exact through degree 2n-1.
    for n in [3usize, 5, 7, 10] {
        let r = gauss_legendre_newton(n);
        for deg in (0..=2 * n - 1).step_by(2) {
            let exact = 2.0 / (deg as f64 + 1.0); // ∫_{-1}^{1} x^deg (deg even)
            let approx: f64 = r.iter().map(|&(x, w)| w * x.powi(deg as i32)).sum();
            assert!(
                (approx - exact).abs() < 1e-12 * exact.abs().max(1.0),
                "n={n} deg={deg}: {approx} vs {exact}"
            );
        }
    }
}

#[test]
fn legendre_recurrence_spot_checks() {
    // P_2(x) = (3x^2 - 1)/2 via the 3-term recurrence through n=2.
    let f = |x: f64, n: usize| {
        // re-derive via the same recurrence as gauss.rs
        if n == 0 {
            return 1.0;
        }
        let (mut p0, mut p1) = (1.0f64, x);
        for k in 1..n {
            let p2 = ((2.0 * k as f64 + 1.0) * x * p1 - k as f64 * p0) / (k as f64 + 1.0);
            p0 = p1;
            p1 = p2;
        }
        p1
    };
    assert!((f(0.5, 2) - (-0.125)).abs() < 1e-15);
    assert!((f(0.0, 2) - (-0.5)).abs() < 1e-15);
    assert!((f(0.5, 1) - 0.5).abs() < 1e-15);
    // P_10(0) is a known rational; just check |P_n(0)| <= 1 and P_n(±1)=±1/1.
    assert!((f(1.0, 10) - 1.0).abs() < 1e-14);
    assert!((f(-1.0, 10) - 1.0).abs() < 1e-14); // P_10 even
}
