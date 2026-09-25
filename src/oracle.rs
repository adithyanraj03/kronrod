//! 100-digit reference constants.
//!
//! Reconstructed from the OEIS digit b-files (`b0000796.txt`, `b001113.txt`,
//! `b002161.txt`) and cross-checked against the OEIS page displays:
//!
//! * [`PI_100`] — A000796, the decimal expansion of π.
//! * [`E_100`] — A001113, the decimal expansion of e.
//! * [`SQRT_PI_100`] — A002161, the decimal expansion of √π.
//!
//! Each string is `"<d>.<100 digits>"` (102 characters). Parsing to `f64`
//! keeps the correctly-rounded double; the extra digits are kept so the
//! provenance is inspectable and the leading-digit KATs are non-vacuous.

/// π to 100 decimal places (OEIS A000796).
pub const PI_100: &str = "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798212801335756";
/// e to 100 decimal places (OEIS A001113).
pub const E_100: &str = "2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200";
/// √π to 100 decimal places (OEIS A002161).
pub const SQRT_PI_100: &str = "1.7724538509055160272981674833411451827975494561223871282138077898529112845910321813749506567385446655359110629547";

/// π as a correctly-rounded `f64`.
pub fn pi() -> f64 {
    PI_100.parse().expect("PI_100 parses")
}

/// e as a correctly-rounded `f64`.
pub fn e() -> f64 {
    E_100.parse().expect("E_100 parses")
}

/// √π as a correctly-rounded `f64`.
pub fn sqrt_pi() -> f64 {
    SQRT_PI_100.parse().expect("SQRT_PI_100 parses")
}
