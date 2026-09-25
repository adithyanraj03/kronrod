//! End-to-end determinism: every artifact the crate can produce is a
//! pure function of the embedded constants. Two in-process generations
//! must be byte-identical (no clock, no environment, no randomness).

use kronrod::battery::render_report;
use kronrod::crypto::{sha256, to_hex};
use kronrod::pdf::render_dossier;
use kronrod::svg::{render_anim, render_contact_sheet};

#[test]
fn report_is_byte_stable() {
    assert_eq!(render_report(), render_report());
}

#[test]
fn dossier_is_byte_stable() {
    let r = render_report();
    assert_eq!(render_dossier(&r), render_dossier(&r));
}

#[test]
fn svg_artifacts_are_byte_stable() {
    assert_eq!(render_anim(16), render_anim(16));
    assert_eq!(
        render_contact_sheet(16, 4),
        render_contact_sheet(16, 4)
    );
}

#[test]
fn full_artifact_set_hashes_are_stable() {
    let h = |bytes: &[u8]| to_hex(&sha256(bytes));
    let report = render_report();
    let dossier = render_dossier(&report);
    let anim = render_anim(16).into_bytes();
    let sheet = render_contact_sheet(16, 4).into_bytes();

    let r1 = [
        h(report.as_bytes()),
        h(&dossier),
        h(&anim),
        h(&sheet),
    ];
    let report = render_report();
    let dossier = render_dossier(&report);
    let anim = render_anim(16).into_bytes();
    let sheet = render_contact_sheet(16, 4).into_bytes();
    let r2 = [
        h(report.as_bytes()),
        h(&dossier),
        h(&anim),
        h(&sheet),
    ];
    assert_eq!(r1, r2, "artifact hashes must be reproducible");
    // And the four artifacts must be mutually distinct.
    assert_ne!(r1[0], r1[1]);
    assert_ne!(r1[1], r1[2]);
    assert_ne!(r1[2], r1[3]);
}

#[test]
fn engine_results_are_bit_stable() {
    // The engine itself is deterministic: the same integral computed
    // twice gives bit-identical results (no floating-point reassociation
    // depends on timing or ordering that could vary).
    let f = |x: f64| (1.0 - x * x).sqrt();
    let a = kronrod::qag(&f, 0.0, 1.0, 0.0, 1e-13, 500);
    let b = kronrod::qag(&f, 0.0, 1.0, 0.0, 1e-13, 500);
    assert_eq!(a, b);
}
