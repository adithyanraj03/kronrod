//! Deterministic SVG artifacts: an SMIL animation (`anim.svg`) and a
//! static contact sheet (`frames.svg`, for README embedding — GitHub
//! strips SMIL).
//!
//! Scene: the adaptive bisection of ∫₀¹ √(1−x²) dx = π/4 (the unit
//! quarter disc) with the (10,21) rule. Each frame shows the quarter
//! disc with the current interval partition (bar height ∝ the
//! interval's error estimate) along the x-axis.
//! Every byte is a pure function of the engine's trace: no clock, no
//! randomness.

use std::fmt::Write;

use crate::qag::{qage_with_rule, TraceStep};
use crate::qk::qk;
use crate::tables;

const W: f64 = 320.0;
const H: f64 = 256.0;
const R: f64 = 140.0; // disc radius, px
const OX: f64 = 80.0; // disc centre x (plot origin)
const OY: f64 = 186.0; // disc centre y (x-axis)

const PAPER: &str = "#F7F6F1";
const INK: &str = "#252524";
const MUTED: &str = "#676662";
const HAIR: &str = "#DCDAD1";
const INDIGO: &str = "#5B51C7";
const GREEN: &str = "#22AC80";

fn f(x: f64) -> f64 {
    (1.0 - x * x).sqrt()
}

/// The QAG trace for the quarter-disc integral (fixed parameters).
fn trace() -> (Vec<TraceStep>, f64) {
    let mut t: Vec<TraceStep> = Vec::new();
    let _ = qage_with_rule(&f, 0.0, 1.0, 0.0, 1e-9, 60, 2, Some(&mut t));
    let rule = tables::rule(2);
    let first = qk(rule, &f, 0.0, 1.0);
    (t, first.abserr)
}

/// Reconstruct the interval partition after `i` bisections, with each
/// interval's true (recomputed) error estimate.
fn partition(steps: &[TraceStep], i: usize, e0: f64) -> Vec<(f64, f64, f64)> {
    let rule = tables::rule(2);
    let mut ivs: Vec<(f64, f64, f64)> = vec![(0.0, 1.0, e0)];
    for s in steps.iter().take(i) {
        // The bisected interval is (s.a1, s.b2); replace it by its halves.
        let pos = ivs
            .iter()
            .position(|(a, b, _)| (a - s.a1).abs() < 1e-15 && (b - s.b2).abs() < 1e-15);
        if let Some(p) = pos {
            ivs.remove(p);
        } else {
            continue;
        }
        let e1 = qk(rule, &f, s.a1, s.b1).abserr;
        let e2 = qk(rule, &f, s.a2, s.b2).abserr;
        ivs.push((s.a1, s.b1, e1));
        ivs.push((s.a2, s.b2, e2));
    }
    ivs.sort_by(|a, b| a.0.total_cmp(&b.0));
    ivs
}

/// The quarter disc: filled body + arc boundary (the integrand itself).
fn disc() -> String {
    let top = OY - R;
    let right = OX + R;
    format!(
        "<path d=\"M {OX:.6} {OY:.6} L {OX:.6} {top:.6} A {R:.6} {R:.6} 0 0 1 {right:.6} {OY:.6} Z\" fill=\"{INDIGO}\" opacity=\"0.07\"/>\
         <path d=\"M {OX:.6} {top:.6} A {R:.6} {R:.6} 0 0 1 {right:.6} {OY:.6}\" fill=\"none\" stroke=\"{INK}\" stroke-width=\"1.6\"/>"
    )
}

/// One frame: disc, axes, error bars, caption.
fn frame(i: usize, steps: &[TraceStep], e0: f64) -> String {
    let ivs = partition(steps, i, e0);
    let errsum: f64 = ivs.iter().map(|(_, _, e)| *e).sum();
    let mut s = String::new();
    let _ = write!(
        s,
        "<rect x=\"0\" y=\"0\" width=\"{W:.0}\" height=\"{H:.0}\" fill=\"{PAPER}\"/>"
    );
    let _ = write!(s, "{}", disc());
    // axes
    let _ = write!(
        s,
        "<line x1=\"{OX:.6}\" y1=\"{OY:.6}\" x2=\"{:.6}\" y2=\"{OY:.6}\" stroke=\"{HAIR}\" stroke-width=\"1\"/>",
        OX + R + 12.0
    );
    // error bars below the x-axis (height ~ error estimate)
    for (a, b, e) in &ivs {
        let h = 2.0 + 24.0 * (e / e0).min(1.0);
        let x0 = OX + a * R;
        let x1 = OX + b * R;
        let _ = write!(
            s,
            "<rect x=\"{x0:.6}\" y=\"{:.6}\" width=\"{:.6}\" height=\"{:.6}\" fill=\"{INDIGO}\" opacity=\"0.7\"/>",
            OY + 2.0,
            x1 - x0,
            h
        );
    }
    // caption
    let _ = write!(
        s,
        "<text x=\"{OX:.0}\" y=\"{:.0}\" font-family=\"Courier\" font-size=\"11\" fill=\"{INK}\">step {i:2}  intervals {n:3}  errsum {e:.4e}</text>",
        OY + 46.0,
        n = ivs.len(),
        e = errsum
    );
    let _ = write!(
        s,
        "<text x=\"{OX:.0}\" y=\"{:.0}\" font-family=\"Courier\" font-size=\"9\" fill=\"{MUTED}\">QAG (10,21) on sqrt(1-x^2), [0,1] = pi/4  (bar height ~ error estimate)</text>",
        OY + 62.0
    );
    s
}

/// The SMIL animation (loops; each frame is visible for 0.5 s).
pub fn render_anim(frames: usize) -> String {
    let (steps, e0) = trace();
    let n = frames.min(1 + steps.len());
    let total = n as f64 * 0.5;
    let mut s = String::new();
    let _ = write!(
        s,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{W:.0}\" height=\"{H:.0}\" viewBox=\"0 0 {W:.0} {H:.0}\">"
    );
    for i in 0..n {
        let ks = i as f64 / n as f64;
        let ke = (i + 1) as f64 / n as f64;
        let _ = write!(
            s,
            "<g opacity=\"0\"><animate attributeName=\"opacity\" calcMode=\"discrete\" values=\"0;1;0\" keyTimes=\"0;{ks:.9};{ke:.9}\" dur=\"{total:.6}s\" repeatCount=\"indefinite\"/>"
        );
        let _ = write!(s, "{}", frame(i, &steps, e0));
        s.push_str("</g>");
    }
    s.push_str("</svg>");
    s
}

/// The static contact sheet: `cols x rows` grid of frames.
pub fn render_contact_sheet(frames: usize, cols: usize) -> String {
    let (steps, e0) = trace();
    let n = frames.min(1 + steps.len());
    let rows = (n + cols - 1) / cols;
    let cw: f64 = W;
    let ch: f64 = H;
    let gap: f64 = 8.0;
    let title_h: f64 = 24.0;
    let width = 8.0 + cols as f64 * cw + (cols - 1) as f64 * gap;
    let height = title_h + rows as f64 * ch + (rows - 1) as f64 * gap + 8.0;
    let mut s = String::new();
    let _ = write!(
        s,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width:.0}\" height=\"{height:.0}\" viewBox=\"0 0 {width:.0} {height:.0}\">"
    );
    let _ = write!(
        s,
        "<text x=\"8\" y=\"14\" font-family=\"Courier\" font-size=\"11\" fill=\"{GREEN}\">Kronrod: adaptive bisection of the unit quarter disc (pi/4)</text>"
    );
    for i in 0..n {
        let col = i % cols;
        let row = i / cols;
        let tx = 8.0 + col as f64 * (cw + gap);
        let ty = title_h + row as f64 * (ch + gap);
        let _ = write!(s, "<g transform=\"translate({tx:.6},{ty:.6})\">");
        let _ = write!(s, "{}", frame(i, &steps, e0));
        s.push_str("</g>");
    }
    s.push_str("</svg>");
    s
}
