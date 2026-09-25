//! Kronrod — QUADPACK-style adaptive Gauss-Kronrod quadrature from scratch:
//! the exact Netlib rule tables, two-method Gauss-Legendre, and the QAG
//! engine, in std-only Rust.

pub mod battery;
pub mod crypto;
pub mod gauss;
pub mod oracle;
pub mod pdf;
pub mod qag;
pub mod qpsrt;
pub mod qk;
pub mod svg;
pub mod tables;
pub mod tables_gen;

/// Crate version (also embedded in artifacts).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
