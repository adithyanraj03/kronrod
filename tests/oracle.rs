//! 100-digit reference constants (π, e, √π) reconstructed from OEIS b-files.

use std::f64::consts;

use kronrod::oracle::{e, pi, sqrt_pi, E_100, PI_100, SQRT_PI_100};

#[test]
fn pi_prefix() {
    assert!(PI_100.starts_with("3.141592653589793238462643383279"));
    assert!(PI_100.starts_with("3.14159265358979323846264338327950288419716939937510"));
}

#[test]
fn e_prefix() {
    assert!(E_100.starts_with("2.718281828459045235360287471352"));
    assert!(E_100.starts_with("2.71828182845904523536028747135266249775724709369995"));
}

#[test]
fn sqrt_pi_prefix() {
    assert!(SQRT_PI_100.starts_with("1.772453850905516027298167483341"));
    assert!(SQRT_PI_100.starts_with("1.77245385090551602729816748334114518279754945612238"));
}

#[test]
fn pi_matches_std() {
    assert_eq!(pi(), consts::PI);
    assert!((pi() - 3.141592653589793).abs() < 1e-14);
}

#[test]
fn e_matches_exp() {
    assert!((e() - 2.718281828459045).abs() < 1e-14);
    // e = exp(1) to double precision.
    assert!((e() - 1.0f64.exp()).abs() <= f64::EPSILON * e());
}

#[test]
fn sqrt_pi_consistent() {
    // (sqrt(pi))^2 == pi and (sqrt(pi))^2 == pi*1 within 1 ulp-class error.
    let s = sqrt_pi() * sqrt_pi();
    assert!((s - consts::PI).abs() < 1e-14);
    // sqrt(pi) ~ 1.77245.
    assert!((sqrt_pi() - 1.772453850905516).abs() < 1e-14);
}

#[test]
fn hundred_digit_strings_are_clean() {
    for s in [PI_100, E_100, SQRT_PI_100] {
        // "d.d..." form: one dot, the rest digits.
        let body = s.split_once('.').unwrap();
        assert_eq!(body.0.len(), 1);
        assert!(body.1.chars().all(|c| c.is_ascii_digit()));
        assert!(body.1.len() >= 99);
    }
}
