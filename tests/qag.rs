//! The adaptive driver (qage.f port): KATs for exact integrals, the
//! `ier` codes, `neval` accounting, key clamping, and the faithful
//! IEEE overflow-wall behaviour on the pole integrand.

use kronrod::qag::{qag, qage_with_rule};
use std::f64::consts::PI;

fn rel(r: f64, exact: f64) -> f64 {
    (r - exact).abs() / exact.abs()
}

#[test]
fn quadratic_single_evaluation() {
    // x^2 is exact for the 21-point rule: one evaluation, early exit.
    let r = qag(&|x| x * x, 0.0, 1.0, 0.0, 1e-13, 100);
    assert_eq!(r.ier, 0);
    assert_eq!(r.neval, 21);
    assert_eq!(r.last, 1);
    assert!(rel(r.result, 1.0 / 3.0) <= 1e-14);
}

#[test]
fn exponential_on_zero_one() {
    let e = std::f64::consts::E;
    let r = qag(&f64::exp, 0.0, 1.0, 0.0, 1e-13, 500);
    assert_eq!(r.ier, 0);
    assert!(rel(r.result, e - 1.0) <= 1e-12, "{}", r.result);
    assert!(r.neval >= 21);
    assert!(r.last >= 1);
}

#[test]
fn log_singular_at_zero() {
    // ∫₀¹ ln x dx = -1; the endpoint singularity is handled (ier stays 0).
    let r = qag(&f64::ln, 0.0, 1.0, 0.0, 1e-13, 1000);
    assert_eq!(r.ier, 0);
    assert!(rel(r.result, -1.0) <= 1e-11, "{}", r.result);
}

#[test]
fn sin_5x_on_zero_pi() {
    let r = qag(&|x| (5.0 * x).sin(), 0.0, PI, 0.0, 1e-13, 500);
    assert_eq!(r.ier, 0);
    assert!(rel(r.result, 0.4) <= 1e-12);
}

#[test]
fn rational_wide_interval() {
    let exact = 10.0f64.atan(); // ∫₀^10 1/(1+x²)
    let r = qag(&|x| 1.0 / (1.0 + x * x), 0.0, 10.0, 0.0, 1e-13, 500);
    assert_eq!(r.ier, 0);
    assert!(rel(r.result, exact) <= 1e-12, "{}", r.result);
}

#[test]
fn gaussian_tail_exact_value() {
    // ∫₋₅⁵ e^{-x²} dx against the 19-digit exact constant
    // (√π·(1 − erfc(5)), from high-accuracy tail quadrature).
    const EXACT: f64 = 1.7724538509027909505;
    let r = qag(&|x| (-x * x).exp(), -5.0, 5.0, 0.0, 1e-13, 500);
    assert_eq!(r.ier, 0);
    assert!(rel(r.result, EXACT) <= 1e-12, "{}", r.result);
}

#[test]
fn quarter_disc_converges() {
    let f = |x: f64| (1.0 - x * x).sqrt();
    let exact = PI / 4.0;
    let r = qag(&f, 0.0, 1.0, 0.0, 1e-13, 500);
    assert_eq!(r.ier, 0);
    assert!(rel(r.result, exact) <= 1e-14, "{}", r.result);
    // neval must follow qage.f's scaling: (10*keyf+1)*(2n+1), keyf=2.
    let n = (r.neval - 21) / 42; // (21*(2n+1) - 21)/42
    assert_eq!(r.neval, 21 * (2 * n + 1));
    assert_eq!(r.last, n + 1);
}

#[test]
fn pole_overflow_wall_faithful_ieee() {
    // 1/sqrt(1-x²) on [0,1]: at bisection depth 46 the outermost Kronrod
    // node of the rightmost interval rounds to exactly x=1 (the pole).
    // The reference double-precision qage.f then computes errbnd = inf and
    // terminates via the IEEE comparison inf <= inf: ier=0, result=+inf,
    // neval=21*(2*46+1)=1953, last=47. The port must reproduce this.
    let f = |x: f64| (1.0 - x * x).recip().sqrt();
    let r = qag(&f, 0.0, 1.0, 0.0, 1e-13, 1000);
    assert_eq!(r.ier, 0);
    assert_eq!(r.neval, 1953);
    assert_eq!(r.last, 47);
    assert!(r.result.is_infinite() && r.result > 0.0);
    assert!(r.abserr.is_infinite() && r.abserr > 0.0);
}

#[test]
fn ier6_invalid_tolerances() {
    // epsabs <= 0 and epsrel below max(50*epmach, 0.5e-14) → ier = 6.
    let r = qag(&f64::exp, 0.0, 1.0, 0.0, 1e-30, 100);
    assert_eq!(r.ier, 6);
    assert_eq!(r.last, 0);
    assert_eq!(r.neval, 0);
    // ...and the 50*epmach boundary itself: 1e-14 < 1.1102e-14 → still 6.
    let r = qag(&f64::exp, 0.0, 1.0, 0.0, 1e-14, 100);
    assert_eq!(r.ier, 6);
    // Above the boundary the run proceeds normally.
    let r = qag(&f64::exp, 0.0, 1.0, 0.0, 1e-13, 100);
    assert_eq!(r.ier, 0);
}

#[test]
fn ier1_limit_reached() {
    let r = qag(&f64::exp, 0.0, 1.0, 0.0, 1e-8, 1);
    assert_eq!(r.ier, 1);
    assert_eq!(r.neval, 21);
    assert_eq!(r.last, 1);
}

#[test]
fn key_clamping_like_qage() {
    // key <= 0 → keyf = 1: single (7,15) evaluation → neval = 30*0+15.
    let r = qage_with_rule(&f64::exp, 0.0, 1.0, 0.0, 1e-12, 500, 0, None);
    assert_eq!(r.neval, 15);
    // key >= 7 → keyf = 6: single (30,61) evaluation → 61*(2*0+1).
    let r = qage_with_rule(&f64::exp, 0.0, 1.0, 0.0, 1e-12, 500, 99, None);
    assert_eq!(r.neval, 61);
}

#[test]
fn odd_rule_seven_fifteen_works() {
    // (G7,K15): odd-n rules include the centre in the Gauss rule.
    let e = std::f64::consts::E;
    let r = qage_with_rule(&f64::exp, 0.0, 1.0, 0.0, 1e-12, 500, 1, None);
    assert_eq!(r.ier, 0);
    assert!(rel(r.result, e - 1.0) <= 1e-12, "{}", r.result);
}

#[test]
fn rule_15_31_works() {
    let e = std::f64::consts::E;
    let r = qage_with_rule(&f64::exp, 0.0, 1.0, 0.0, 1e-12, 500, 3, None);
    assert_eq!(r.ier, 0);
    assert!(rel(r.result, e - 1.0) <= 1e-12);
    assert_eq!(r.neval, 31); // single evaluation: (10*3+1)*(2*0+1)
}

#[test]
fn reversed_limits_anti_symmetry() {
    let f = |x: f64| 1.0 / (1.0 + x * x);
    let r1 = qag(&f, 0.0, 1.0, 0.0, 1e-13, 100);
    let r2 = qag(&f, 1.0, 0.0, 0.0, 1e-13, 100);
    assert!((r1.result + r2.result).abs() <= f64::EPSILON * r1.result.abs());
}

#[test]
fn symmetric_interval_quadratic() {
    let r = qag(&|x| x * x, -1.0, 1.0, 0.0, 1e-13, 100);
    assert_eq!(r.ier, 0);
    assert!(rel(r.result, 2.0 / 3.0) <= 1e-14);
}

#[test]
fn neval_scaling_keyf1() {
    // keyf=1 scales as 30n + 15 (qage.f line 340).
    let e = std::f64::consts::E;
    let r = qage_with_rule(&|x: f64| x * x * x.exp(), 0.0, 1.0, 0.0, 1e-13, 500, 1, None);
    let n = (r.neval - 15) / 30;
    assert_eq!(r.neval, 30 * n + 15);
    assert!(rel(r.result, e - 2.0) <= 1e-11, "{}", r.result);
}

#[test]
fn trace_recapitulates_bisection() {
    let f = |x: f64| (1.0 - x * x).sqrt();
    let mut trace: Vec<kronrod::TraceStep> = Vec::new();
    let r = qage_with_rule(&f, 0.0, 1.0, 0.0, 1e-12, 500, 2, Some(&mut trace));
    assert_eq!(r.ier, 0);
    assert!(!trace.is_empty());
    // The first bisection splits [0,1] at 0.5.
    assert_eq!(trace[0].a1, 0.0);
    assert_eq!(trace[0].b1, 0.5);
    assert_eq!(trace[0].a2, 0.5);
    assert_eq!(trace[0].b2, 1.0);
    // `last` grows by one per step and matches the driver's count.
    for (i, s) in trace.iter().enumerate() {
        assert_eq!(s.last, i + 2);
    }
    assert_eq!(r.last, trace.len() + 1);
    // errsum is non-increasing in expectation... at least finite here.
    assert!(trace.iter().all(|s| s.errsum.is_finite()));
}
