//! # Kronrod
//!
//! A QUADPACK-style adaptive Gauss-Kronrod quadrature engine, rebuilt from
//! scratch in dependency-free Rust (`std` only) and verified against the
//! official Netlib constants, including the (G10, K21) table evaluated
//! with 80-decimal-digit arithmetic by L. W. Fullerton, Bell Labs,
//! November 1981 (stored in Netlib `dqk21.f` at 33 digits and
//! cross-checked here against arbitrary-precision Laurie-method tables).
//!
//! ## What this crate contains
//!
//! * [`tables`] — the six official (G_n, K_{2n+1}) rules for the uniform
//!   weight on [-1, 1], embedded as the exact source digit strings
//!   ([`tables_gen`]) and parsed once to `f64`.
//! * [`gauss`] — Gauss-Legendre rules from scratch via two independent
//!   methods (Newton on the three-term recurrence; Golub-Welsch Jacobi
//!   eigendecomposition), cross-checked against the embedded tables to
//!   full `f64` precision.
//! * [`qk`] — the rule evaluator: a faithful port of the `dqk15`/`dqk21`/
//!   …/`dqk61` evaluation block including QUADPACK's error estimator
//!   (the `resasc·min(10, (200·abserr/resasc)^1.5)` scaling and the
//!   `50·epmach` roundoff floor).
//! * [`qpsrt`] — faithful 1-based port of `qpsrt.f` (list reordering).
//! * [`qag`] — faithful port of the `qage.f` adaptive driver: roundoff
//!   guards (`iroff1`/`iroff2`), bad-behaviour detection, `neval`
//!   accounting, and the `ier` codes (0, 1, 2, 3, 6).
//! * [`oracle`] — 100-digit reference constants (π, e, √π) reconstructed
//!   from OEIS digit b-files.
//! * [`battery`] — the estimator stress battery (16 integrands with exact
//!   values), the quarter-disc (π/4) convergence study, the
//!   rule-family × tolerance grid, and the pole overflow-wall probe; all
//!   rendered to a deterministic text report.
//! * [`pdf`] — hand-rolled deterministic PDF 1.4 writer (base-14 fonts,
//!   fixed `/ID`, fixed `CreationDate`; no wall-clock input anywhere).
//! * [`svg`] — deterministic per-frame SVG (SMIL animation + static
//!   contact sheet) of the adaptive bisection of the quarter circle.
//! * [`crypto`] — zero-dependency SHA-256 (FIPS 180-4) used for the
//!   dual-run byte-identity attestation.
//! * [`attest`] — dual-run attestation: the full battery and dossier are
//!   generated twice in-process and compared byte-for-byte (no wall
//!   clock, so identity is exact).
//!
//! ## Determinism policy
//!
//! Every artifact this crate can produce (report, dossier, attestation,
//! SVG) is a pure function of the embedded constants: no timestamps, no
//! environment input, no randomness. Two runs on the same machine produce
//! byte-identical files; this is asserted by the test suite and by
//! `kronrod attest`.
//!
//! ## Sources
//!
//! * QUADPACK, Piessens, Doncker-Kapsteijn and Uytterhoeven (1983),
//!   Netlib sources `qage.f`, `qpsrt.f`, `dqk21.f`, `qk15.f`, `qk31.f`,
//!   `qk41.f`, `qk51.f`, `qk61.f`.
//! * D. P. Laurie, "Calculation of Gauss-Kronrod Quadrature Rules",
//!   Math. Comp. 66 (1997) 1133-1145 (degree of exactness 3n+1).
//! * L. W. Fullerton, Bell Labs, "Gauss quadrature weights and Kronrod
//!   quadrature abscissae and weights as evaluated with 80 decimal digit
//!   arithmetic", November 1981 (as distributed in Netlib `dqk21.f`).
//! * P. Holoborodko, "Gauss-Kronrod Quadrature Nodes and Weights"
//!   (arbitrary-precision Laurie-method tables, 2011) — the independent
//!   oracle used to cross-check the embedded tables.
//! * G. H. Golub and J. H. Welsch, "Calculation of Gauss-Quadrature
//!   Rules", Math. Comp. 23 (1969) 221-230.

pub mod attest;
pub mod battery;
pub mod crypto;
pub mod gauss;
pub mod oracle;
pub mod pdf;
pub mod qag;
pub mod qk;
pub mod qpsrt;
pub mod svg;
pub mod tables;
pub mod tables_gen;

pub use crate::attest::{attest, render_attestation, Attestation};
pub use crate::battery::{grid, pi_study, run_battery};
pub use crate::gauss::{gauss_legendre_golub_welsch, gauss_legendre_newton};
pub use crate::oracle::{e, pi, sqrt_pi};
pub use crate::qag::{qag, qage_with_rule, QagResult, TraceStep};
pub use crate::qk::{estimate_error, qk, QkResult};
pub use crate::tables::{rule, rules, Rule, KEYF_ALL};

/// Crate version (also embedded in artifacts).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
