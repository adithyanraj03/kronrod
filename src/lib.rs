//! Kronrod — QUADPACK-style adaptive Gauss-Kronrod quadrature from scratch:
//! the exact Netlib rule tables, two-method Gauss-Legendre, and the QAG
//! engine, in std-only Rust.

/// Crate version (also embedded in artifacts).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
