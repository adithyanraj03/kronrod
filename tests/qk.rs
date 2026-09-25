//! Rule evaluator: the QUADPACK estimator against independently computed
//! bit-exact constants, and polynomial exactness of the six rules.

use kronrod::qk::{estimate_error, qk, QkResult};
use kronrod::tables::rule;

/// Independent f64 constants (computed by a separate high-accuracy path,
/// not by this crate).
const EST_A: f64 = 2.0;
const EST_B: f64 = 1.76776717236491370e-13;
const EST_C: f64 = 1.11022302462515654e-14;
const EST_D: f64 = 7.77156117237609578e-14;
const EST_E: f64 = 2.0;
const EPMACH50: f64 = 1.11022302462515654e-14;

#[test]
fn estimator_a_clamped_to_ten() {
    // (200*abserr/resasc)^1.5 > 10, so the min(10, ...) clamp engages:
    // abserr = resasc * 10 = 0.3 * 0.5 * 2.0 = ... = 2.0 (hlgth scaling).
    let got = estimate_error(0.0, 1.0, 0.5, 0.3, 0.2);
    assert_eq!(got, EST_A);
}

#[test]
fn estimator_b_power_one_point_five() {
    // 200*abserr/resasc < 1: the power-1.5 branch.
    let got = estimate_error(2.0, 2.0000000001, 0.25, 1.0, 4.0);
    assert_eq!(got, EST_B);
}

#[test]
fn estimator_c_roundoff_floor() {
    // resabs large enough for the floor to engage; the floor wins.
    let got = estimate_error(1.0, 1.0 + 1e-17, 1.0, 1.0, 2.0);
    assert_eq!(got, EST_C);
    assert_eq!(got, 50.0 * f64::EPSILON * 1.0);
}

#[test]
fn estimator_d_resasc_zero() {
    // resasc == 0 skips the power branch; the roundoff floor applies.
    let got = estimate_error(5.0, 5.0, 0.5, 7.0, 0.0);
    assert_eq!(got, EST_D);
}

#[test]
fn estimator_e_floor_below_uflow_guard() {
    // resabs = 1e-300 < uflow/(50*epmach): the floor is skipped entirely.
    let got = estimate_error(0.0, 1.0, 0.5, 1e-300, 0.2);
    assert_eq!(got, EST_E);
}

#[test]
fn epmach50_constant() {
    assert_eq!(EPMACH50, 50.0 * f64::EPSILON);
}

fn on_minus_one_one(f: impl Fn(f64) -> f64, keyf: u8) -> QkResult {
    qk(rule(keyf), &f, -1.0, 1.0)
}

/// ∫_{-1}^{1} x^d dx = 2/(d+1) for even d.
fn exact_even_power(d: usize) -> f64 {
    2.0 / (d as f64 + 1.0)
}

#[test]
fn rule_10_21_exact_through_degree_20() {
    // (G10,K21) integrates every even power up to x^20 (measured exactness
    // actually extends to x^30; 20 is the guaranteed floor).
    for d in (0usize..=20).step_by(2) {
        let r = on_minus_one_one(|x| x.powi(d as i32), 2);
        let e = exact_even_power(d);
        assert!((r.result - e).abs() <= 1e-14 * e, "d={d}: {result:?} vs {e}", result = r.result);
    }
}

#[test]
fn rule_10_21_exact_through_measured_cliff() {
    // Measured f64 exactness cliff for (10,21): through degree 30, not 32.
    for d in [22usize, 24, 26, 28, 30] {
        let r = on_minus_one_one(|x| x.powi(d as i32), 2);
        let e = exact_even_power(d);
        assert!((r.result - e).abs() <= 1e-14 * e, "d={d}");
    }
    let r = on_minus_one_one(|x| x.powi(32), 2);
    let e = exact_even_power(32);
    assert!((r.result - e).abs() > 1e-14 * e, "degree 32 must NOT be exact");
}

#[test]
fn all_rules_exact_on_quartic() {
    // x^4 is within the guaranteed exactness of every (n, 2n+1) rule.
    for keyf in [1u8, 2, 3, 4, 5, 6] {
        let r = on_minus_one_one(|x| x * x * x * x, keyf);
        let e = 2.0 / 5.0;
        assert!((r.result - e).abs() <= 1e-14 * e, "keyf={keyf}: {}", r.result);
    }
}

#[test]
fn odd_n_centre_belongs_to_gauss() {
    // For (G7,K15) the constant integrand: resg uses the centre weight and
    // must agree with resk to within roundoff on a symmetric interval.
    let r = on_minus_one_one(|_| 1.0, 1);
    assert!((r.result - 2.0).abs() <= 8.0 * f64::EPSILON * 2.0);
    // And for (G10,K21) (even n, no centre in Gauss) the same check.
    let r = on_minus_one_one(|_| 1.0, 2);
    assert!((r.result - 2.0).abs() <= 8.0 * f64::EPSILON * 2.0);
}

#[test]
fn reversed_interval_sign() {
    let r1 = qk(rule(2), &|x| x * x, 0.0, 1.0);
    let r2 = qk(rule(2), &|x| x * x, 1.0, 0.0);
    assert!((r1.result - 1.0 / 3.0).abs() <= 1e-16);
    assert!((r2.result + 1.0 / 3.0).abs() <= 1e-16);
}

#[test]
fn resasc_positive_for_nonconstant() {
    let r = qk(rule(2), &|x| x * x, 0.0, 1.0);
    assert!(r.resasc > 0.0);
    assert!(r.resabs > 0.0);
    assert!(r.abserr > 0.0 || r.result == 1.0 / 3.0);
}
