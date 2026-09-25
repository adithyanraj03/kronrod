//! Dual-run byte-identity attestation.

use kronrod::attest::{attest, render_attestation};
use kronrod::crypto::{sha256, to_hex};

#[test]
fn dual_run_is_byte_identical() {
    let a = attest();
    assert!(a.identical, "determinism violated");
    assert_eq!(a.report_sha256_run1, a.report_sha256_run2);
    assert_eq!(a.dossier_sha256_run1, a.dossier_sha256_run2);
}

#[test]
fn hashes_are_sha256_hex() {
    let a = attest();
    for h in [
        &a.report_sha256_run1,
        &a.report_sha256_run2,
        &a.dossier_sha256_run1,
        &a.dossier_sha256_run2,
    ] {
        assert_eq!(h.len(), 64);
        assert!(h.chars().all(|c| c.is_ascii_hexdigit()));
    }
}

#[test]
fn report_hash_matches_independent_sha256() {
    let a = attest();
    let report = kronrod::battery::render_report();
    let h = to_hex(&sha256(report.as_bytes()));
    assert_eq!(h, a.report_sha256_run1);
    let dossier = kronrod::pdf::render_dossier(&report);
    let h = to_hex(&sha256(&dossier));
    assert_eq!(h, a.dossier_sha256_run1);
}

#[test]
fn rendered_record_is_self_describing() {
    let text = render_attestation(&attest());
    assert!(text.contains("KRONROD DETERMINISM ATTESTATION"));
    assert!(text.contains("BYTE-IDENTICAL (dual-run diff empty)"));
    assert!(text.contains("kronrod 1.0.0"));
    // No wall-clock input anywhere in the record.
    assert!(!text.contains("D:2")); // PDF-style timestamps
    assert!(text.chars().all(|c| c.is_ascii()));
}

#[test]
fn attestation_version_matches_crate() {
    let a = attest();
    assert_eq!(a.version, format!("kronrod {}", kronrod::VERSION));
}
