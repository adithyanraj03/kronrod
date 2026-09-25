//! The deterministic SVG artifacts: structure, SMIL animation, and
//! byte-stability.

use kronrod::svg::{render_anim, render_contact_sheet};

#[test]
fn animation_has_all_frames() {
    let svg = render_anim(16);
    assert!(svg.starts_with("<svg "));
    assert!(svg.ends_with("</svg>"));
    assert_eq!(svg.matches("<g opacity=\"0\">").count(), 16);
    assert_eq!(svg.matches("<animate ").count(), 16);
    assert!(svg.contains("repeatCount=\"indefinite\""));
    assert!(svg.contains("calcMode=\"discrete\""));
}

#[test]
fn animation_frames_are_deterministic() {
    assert_eq!(render_anim(16), render_anim(16));
    assert_eq!(render_anim(8), render_anim(8));
    assert_ne!(render_anim(16), render_anim(8));
}

#[test]
fn contact_sheet_layout() {
    let svg = render_contact_sheet(16, 4);
    assert!(svg.starts_with("<svg "));
    assert!(svg.ends_with("</svg>"));
    assert_eq!(svg.matches("<g transform=\"translate(").count(), 16);
    assert!(svg.contains("quarter disc"));
    assert!(svg.contains("pi/4"));
}

#[test]
fn contact_sheet_is_deterministic() {
    assert_eq!(render_contact_sheet(16, 4), render_contact_sheet(16, 4));
    assert_ne!(render_contact_sheet(16, 4), render_contact_sheet(16, 2));
}

#[test]
fn frames_are_ascii_and_finite() {
    for s in [render_anim(16), render_contact_sheet(16, 4)] {
        assert!(s.chars().all(|c| c.is_ascii()));
        assert!(!s.contains("NaN"));
        assert!(!s.contains("inf"));
    }
}

#[test]
fn frames_progress_through_bisection() {
    // The bar partitions grow: later frames must contain more rect
    // elements (each bisection replaces one interval by two).
    let svg = render_anim(16);
    let frames: Vec<&str> = svg.split("<g opacity=\"0\">").collect();
    assert_eq!(frames.len(), 17); // header + 16 frames
    let rects = |f: &str| f.matches("<rect ").count();
    let first = rects(frames[1]);
    let last = rects(frames[16]);
    // Frame 0: background rect + the single initial interval bar.
    assert_eq!(first, 2);
    assert!(last >= 11, "expected a deep bisection: {last} rects");
}
