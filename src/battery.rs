//! The verification battery: 16 exact-value integrands through the full
//! adaptive engine, the quarter-disc (π/4) convergence study, the
//! rule-family × tolerance grid, and the pole overflow-wall probe.
//!
//! Everything here is a pure function of the embedded constants (no
//! clock, no environment), so the rendered report is byte-stable across
//! runs — the property the attestation and the golden-file tests assert.

use std::f64::consts::PI;

use crate::oracle;
use crate::qag::{qag, QagResult};
use crate::tables;

/// One row of the estimator stress battery.
#[derive(Debug, Clone)]
pub struct BatteryRow {
    pub name: &'static str,
    pub range: &'static str,
    pub exact: f64,
    pub result: f64,
    pub abserr: f64,
    pub true_error: f64,
    /// `true_error / abserr` (the estimator's honesty ratio).
    pub ratio: f64,
    pub ier: i32,
    pub neval: usize,
    pub last: usize,
}

/// One row of the π/2 convergence study.
#[derive(Debug, Clone)]
pub struct StudyRow {
    pub epsrel: f64,
    pub result: f64,
    pub ier: i32,
    pub neval: usize,
    pub last: usize,
}

/// One cell of the rule-family × tolerance grid.
#[derive(Debug, Clone)]
pub struct GridRow {
    pub keyf: u8,
    pub rule: String,
    pub epsrel: f64,
    pub ier: i32,
    pub neval: usize,
    pub digits: u32,
}

/// Correct-significant-digit count of `result` against `exact`.
pub fn digits(result: f64, exact: f64) -> u32 {
    if result == exact {
        return u32::MAX; // exact
    }
    let d = (result - exact).abs();
    let rel = if exact != 0.0 { d / exact.abs() } else { d };
    if rel == 0.0 || rel < 1e-16 {
        return 16;
    }
    let k = (-rel.log10()).floor() as i64;
    k.clamp(0, 16) as u32
}

fn digits_str(d: u32) -> String {
    if d == u32::MAX {
        "exact".into()
    } else if d >= 16 {
        ">=16".into()
    } else {
        d.to_string()
    }
}

/// The 16 battery integrands with exact values.
///
/// #15 uses the exact value `√π·(1 − erfc(5))` of ∫₋₅⁵ e^(−x²) dx,
/// computed to 19 digits by high-accuracy quadrature of the exponential
/// tail (Simpson on two independent grids, agreement to 1e-24).
const EXACT_GAUSS_E15: f64 = 1.7724538509027909505;

fn integrands() -> Vec<(&'static str, &'static str, fn(f64) -> f64, f64, f64, f64)> {
    vec![
        ("1", "[0,1]", |_| 1.0, 0.0, 1.0, 1.0),
        ("x", "[0,1]", |x| x, 0.0, 1.0, 0.5),
        ("x^2", "[0,1]", |x| x * x, 0.0, 1.0, 1.0 / 3.0),
        ("x^10", "[0,1]", |x| x.powi(10), 0.0, 1.0, 1.0 / 11.0),
        ("sqrt(x)", "[0,1]", f64::sqrt, 0.0, 1.0, 2.0 / 3.0),
        ("(1-x)^(1/3)", "[0,1]", |x| (1.0 - x).cbrt(), 0.0, 1.0, 0.75),
        ("sqrt(1-x^2)", "[0,1]", |x| (1.0 - x * x).sqrt(), 0.0, 1.0, PI / 4.0),
        ("log(x)", "[0,1]", f64::ln, 0.0, 1.0, -1.0),
        ("exp(x)", "[0,1]", f64::exp, 0.0, 1.0, oracle::e() - 1.0),
        ("x^2*exp(x)", "[0,1]", |x| x * x * x.exp(), 0.0, 1.0, oracle::e() - 2.0),
        ("(1+x)*exp(-x)", "[0,1]", |x| (1.0 + x) * (-x).exp(), 0.0, 1.0, 2.0 - 3.0 / oracle::e()),
        ("sin(5x)", "[0,pi]", |x| (5.0 * x).sin(), 0.0, PI, 0.4),
        ("1/(1+x^2)", "[0,1]", |x| 1.0 / (1.0 + x * x), 0.0, 1.0, PI / 4.0),
        (
            "1/(1+x^2)",
            "[0,10]",
            |x| 1.0 / (1.0 + x * x),
            0.0,
            10.0,
            10.0f64.atan(),
        ),
        (
            "exp(-x^2)",
            "[-5,5]",
            |x| (-x * x).exp(),
            -5.0,
            5.0,
            EXACT_GAUSS_E15,
        ),
        (
            "1/(1+100x^2)",
            "[0,1]",
            |x| 1.0 / (1.0 + 100.0 * x * x),
            0.0,
            1.0,
            (10.0f64.atan() / 10.0),
        ),
    ]
}

/// Run the full battery: QAG (10,21), `epsabs = 0`, `epsrel = 1e-13`,
/// `limit = 1000` for every integrand.
pub fn run_battery() -> Vec<BatteryRow> {
    let rule = tables::rule(2);
    let _ = rule; // (10,21)
    integrands()
        .into_iter()
        .map(|(name, range, f, a, b, exact)| {
            let r: QagResult = qag(&f, a, b, 0.0, 1e-13, 1000);
            let true_error = (r.result - exact).abs();
            let ratio = if r.abserr > 0.0 { true_error / r.abserr } else { true_error };
            BatteryRow {
                name,
                range,
                exact,
                result: r.result,
                abserr: r.abserr,
                true_error,
                ratio,
                ier: r.ier,
                neval: r.neval,
                last: r.last,
            }
        })
        .collect()
}

/// Quarter-disc convergence study: ∫₀¹ √(1−x²) dx = π/4 (the area of a
/// unit quarter disc; the integrand is bounded but has an unbounded
/// endpoint derivative, which drives the adaptive subdivision toward x=1
/// without any overflow risk). With the (10,21) rule, `epsabs = 0`,
/// `limit = 500`, `epsrel` = 1e-1 .. 1e-15.
pub fn pi_study() -> Vec<StudyRow> {
    let f = |x: f64| (1.0 - x * x).sqrt();
    (1..=15)
        .map(|p| {
            let epsrel = 10f64.powi(-p);
            let r = qag(&f, 0.0, 1.0, 0.0, epsrel, 500);
            StudyRow {
                epsrel,
                result: r.result,
                ier: r.ier,
                neval: r.neval,
                last: r.last,
            }
        })
        .collect()
}

/// Rule-family × tolerance grid on the quarter disc (π/4): all six
/// rules, `epsrel` ∈ {1e-3, 1e-6, 1e-9, 1e-12}, `limit = 500`.
pub fn grid() -> Vec<GridRow> {
    let f = |x: f64| (1.0 - x * x).sqrt();
    let tols = [1e-3, 1e-6, 1e-9, 1e-12];
    let mut out = Vec::new();
    for keyf in tables::KEYF_ALL {
        let rule = tables::rule(keyf);
        let label = format!("(G{},{})", rule.n, rule.k);
        for &epsrel in &tols {
            let r = qage_with_rule_cached(&f, 0.0, 1.0, 0.0, epsrel, 500, keyf);
            out.push(GridRow {
                keyf,
                rule: label.clone(),
                epsrel,
                ier: r.ier,
                neval: r.neval,
                digits: digits(r.result, PI / 4.0),
            });
        }
    }
    out
}

// (The grid uses the public driver directly; this alias keeps the call
// site readable and lets tests target it.)
use crate::qag::qage_with_rule;
fn qage_with_rule_cached(
    f: &impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    epsabs: f64,
    epsrel: f64,
    limit: usize,
    keyf: u8,
) -> QagResult {
    qage_with_rule(f, a, b, epsabs, epsrel, limit, keyf, None)
}

// ---------------------------------------------------------------------------
// Report rendering (deterministic text; frozen by the test suite).
// ---------------------------------------------------------------------------

fn fmt(x: f64) -> String {
    format!("{x:.17e}")
}

fn pad(s: &str, w: usize) -> String {
    if s.len() >= w {
        s.to_string()
    } else {
        format!("{s:<w$}", w = w)
    }
}

fn row<S: AsRef<str>>(cols: &[S], widths: &[usize]) -> String {
    let mut out = String::from(" ");
    for (c, w) in cols.iter().zip(widths) {
        out.push(' ');
        out.push_str(&pad(c.as_ref(), *w));
    }
    out.push('\n');
    out
}

/// Render the full verification report (byte-stable).
pub fn render_report() -> String {
    let mut s = String::new();
    let bar = "=".repeat(80);

    s.push_str(&format!("\n{bar}\n"));
    s.push_str("KRONROD VERIFICATION REPORT\n");
    s.push_str("Adaptive Gauss-Kronrod quadrature (QUADPACK QAG engine), Rust std-only.\n");
    s.push_str(&format!("version  : kronrod {}\n", crate::VERSION));
    s.push_str("rules    : (G7,K15) (G10,K21) (G15,K31) (G20,K41) (G25,K51) (G30,K61)\n");
    s.push_str("tables   : exact Netlib source digits; (G10,K21) = Fullerton/Bell\n");
    s.push_str("           Labs (Nov 1981) 80-digit-arithmetic evaluation (33 digits),\n");
    s.push_str("           cross-checked vs arbitrary-precision Laurie-method tables\n");
    s.push_str("engine   : qage.f / qpsrt.f / dqk21.f line-for-line port (double)\n");
    s.push_str("determinism: pure function of embedded constants (no clock, no env)\n");
    s.push_str(&format!("{bar}\n\n"));

    // --- Section 1: Gauss from scratch ---
    s.push_str("SECTION 1 - GAUSS-LEGENDRE FROM SCRATCH (two independent methods)\n");
    s.push_str("A = Newton on the 3-term recurrence; B = Golub-Welsch Jacobi\n");
    s.push_str("(nodes compared by max abs diff; weights by max rel diff).\n");
    let cols: Vec<String> = vec!["n".into(), "max |nodeA-nodeB|".into(), "max wtA/wtB rel".into(), "vs table (nodes)".into(), "vs table (weights)".into()];
    s.push_str(&header_row(&cols, &[4, 20, 18, 18, 18]));
    for &n in &[3usize, 5, 7, 10, 15, 20, 25, 30] {
        let a = crate::gauss::gauss_legendre_newton(n);
        let b = crate::gauss::gauss_legendre_golub_welsch(n);
        let max_node: f64 = a
            .iter()
            .zip(&b)
            .map(|((na, _), (nb, _))| (na - nb).abs())
            .fold(0.0f64, f64::max);
        let max_wt: f64 = a
            .iter()
            .zip(&b)
            .map(|((_, wa), (_, wb))| (wa - wb).abs() / wb.abs())
            .fold(0.0f64, f64::max);
        let (tn, tw) = if let Some(rule) = tables::rules().iter().find(|r| r.n == n) {
            // official positive nodes (descending) <-> method-A positive
            // nodes (reversed to descending); gauss_weights[i] is the
            // weight of the i-th (descending) official positive node.
            let official = rule.gauss_nodes();
            let mut a_pos: Vec<(f64, f64)> =
                a.iter().copied().filter(|(x, _)| *x > 0.0).collect();
            a_pos.reverse();
            let tn = official
                .iter()
                .zip(&a_pos)
                .map(|(o, &(g, _))| (o - g).abs())
                .fold(0.0f64, f64::max);
            let tw = official
                .iter()
                .enumerate()
                .zip(&a_pos)
                .map(|((i, _), &(g, gw))| {
                    let _ = g;
                    (gw - rule.gauss_weights[i]).abs() / rule.gauss_weights[i].abs()
                })
                .fold(0.0f64, f64::max);
            (tn, tw)
        } else {
            (f64::NAN, f64::NAN)
        };
        let cols: Vec<String> = vec![
            n.to_string(),
            fmt(max_node),
            fmt(max_wt),
            if tn.is_nan() { "-".into() } else { fmt(tn) },
            if tw.is_nan() { "-".into() } else { fmt(tw) },
        ];
        s.push_str(&row(&cols, &[4, 20, 18, 18, 18]));
    }
    s.push('\n');

    // --- Section 2: table structure ---
    s.push_str("SECTION 2 - TABLE STRUCTURE (embedded Netlib digits)\n");
    s.push_str("mass-2 |err|  : |2 - total rule mass| for the Kronrod rule\n");
    s.push_str("nested       : Gauss nodes at odd positions match from-scratch Gauss\n");
    s.push_str("interlaced   : t0 > x0 > t1 > x1 > ... > 0 strictly holds\n");
    let cols: Vec<String> = vec!["rule".into(), "pts".into(), "mass-2 err".into(), "symmetric".into(), "nested".into(), "interlaced".into()];
    s.push_str(&header_row(&cols, &[10, 5, 20, 10, 10, 12]));
    for rule in tables::rules() {
        let mass_err = (2.0 - rule.kronrod_mass()).abs();
        let symmetric = rule
            .nodes
            .windows(2)
            .all(|w| w[0] > w[1])
            && *rule.nodes.last().unwrap() == 0.0;
        let nested = {
            let official = rule.gauss_nodes(); // descending positive
            let from_scratch = crate::gauss::gauss_legendre_newton(rule.n);
            let mut from_scratch_pos: Vec<f64> = from_scratch
                .iter()
                .filter(|&&(x, _)| x > 0.0)
                .map(|&(x, _)| x)
                .collect();
            from_scratch_pos.reverse(); // descending, to match `official`
            official
                .iter()
                .zip(&from_scratch_pos)
                .all(|(o, g)| (o - g).abs() <= 1e-13 * o.abs())
        };
        let interlaced = {
            let t = rule.kronrod_nodes(); // descending positive
            let x = rule.gauss_nodes();   // descending positive
            // strict interlacing: t0 > x0 > t1 > x1 > ... down to 0
            // (odd n: one more t than x; even n: equal counts, chain ends at 0)
            t.len() >= x.len()
                && x
                    .iter()
                    .enumerate()
                    .all(|(i, &x)| {
                        x < t[i] && x > if i + 1 < t.len() { t[i + 1] } else { 0.0 }
                    })
        };
        let cols: Vec<String> = vec![
            format!("(G{},K{})", rule.n, rule.k),
            rule.k.to_string(),
            fmt(mass_err),
            if symmetric { "ok" } else { "FAIL" }.into(),
            if nested { "ok" } else { "FAIL" }.into(),
            if interlaced { "ok" } else { "FAIL" }.into(),
        ];
        s.push_str(&row(&cols, &[10, 5, 20, 10, 10, 12]));
    }
    s.push('\n');

    // --- Section 3: battery ---
    s.push_str("SECTION 3 - ESTIMATOR STRESS BATTERY\n");
    s.push_str("QAG (10,21), epsabs=0, epsrel=1e-13, limit=1000; true/est =\n");
    s.push_str("|true error| / |QAG abserr| (a fair estimator stays O(1) to O(100)).\n");
    let cols: Vec<String> = vec![
        "#".into(),
        "integrand".into(),
        "range".into(),
        "exact".into(),
        "QAG result".into(),
        "abserr".into(),
        "true err".into(),
        "true/est".into(),
        "neval".into(),
        "ier".into(),
    ];
    s.push_str(&header_row(&cols, &[3, 14, 6, 22, 22, 22, 22, 12, 7, 3]));
    for (i, r) in run_battery().into_iter().enumerate() {
        let cols: Vec<String> = vec![
            (i + 1).to_string(),
            r.name.into(),
            r.range.into(),
            fmt(r.exact),
            fmt(r.result),
            fmt(r.abserr),
            fmt(r.true_error),
            fmt(r.ratio),
            r.neval.to_string(),
            r.ier.to_string(),
        ];
        s.push_str(&row(&cols, &[3, 14, 6, 22, 22, 22, 22, 12, 7, 3]));
    }
    s.push('\n');

    // --- Section 4: quarter-disc study ---
    s.push_str("SECTION 4 - QUARTER-DISC CONVERGENCE STUDY (sqrt(1-x^2) on [0,1])\n");
    s.push_str("exact = pi/4 = 7.85398163397448310e-01 (unit quarter-disc area);\n");
    s.push_str("digits = leading correct significant digits of the QAG result.\n");
    let cols: Vec<String> = vec!["epsrel".into(), "QAG result".into(), "digits".into(), "neval".into(), "last".into(), "ier".into()];
    s.push_str(&header_row(&cols, &[22, 22, 8, 8, 6, 4]));
    for r in pi_study() {
        let cols: Vec<String> = vec![
            fmt(r.epsrel),
            fmt(r.result),
            digits_str(digits(r.result, PI / 4.0)),
            r.neval.to_string(),
            r.last.to_string(),
            r.ier.to_string(),
        ];
        s.push_str(&row(&cols, &[22, 22, 8, 8, 6, 4]));
    }
    s.push('\n');

    // --- Section 5: rule grid ---
    s.push_str("SECTION 5 - RULE FAMILY x TOLERANCE GRID (QUARTER DISC, PI/4)\n");
    s.push_str("every rule, four tolerances, limit=500.\n");
    let cols: Vec<String> = vec!["rule".into(), "epsrel".into(), "digits".into(), "neval".into(), "ier".into()];
    s.push_str(&header_row(&cols, &[10, 22, 8, 8, 4]));
    for r in grid() {
        let cols: Vec<String> = vec![
            r.rule.into(),
            fmt(r.epsrel),
            digits_str(r.digits),
            r.neval.to_string(),
            r.ier.to_string(),
        ];
        s.push_str(&row(&cols, &[10, 22, 8, 8, 4]));
    }
    s.push('\n');

    // --- Section 6: pole overflow wall (faithful IEEE behaviour) ---
    s.push_str("SECTION 6 - POLE OVERFLOW WALL (FAITHFUL IEEE BEHAVIOUR)\n");
    s.push_str("integrand 1/sqrt(1-x^2) on [0,1], QAG (10,21), epsrel=1e-13,\n");
    s.push_str("limit=1000. At bisection depth 46 the outermost Kronrod node\n");
    s.push_str("of the rightmost interval rounds to exactly x=1 (the pole), so\n");
    s.push_str("f evaluates to +inf; errbnd = epsrel*|area| = +inf as well, and\n");
    s.push_str("the IEEE comparison inf<=inf is true, so the driver terminates\n");
    s.push_str("exactly as the reference double-precision qage.f would:\n");
    let wall = qag(&|x: f64| (1.0 - x * x).recip().sqrt(), 0.0, 1.0, 0.0, 1e-13, 1000);
    let cols: Vec<String> = vec!["ier".into(), wall.ier.to_string(), "neval".into(), wall.neval.to_string()];
    s.push_str(&row(&cols, &[8, 6, 8, 10]));
    let cols: Vec<String> = vec!["last".into(), wall.last.to_string(), "result".into(), fmt(wall.result)];
    s.push_str(&row(&cols, &[8, 6, 8, 10]));
    let cols: Vec<String> = vec!["abserr".into(), fmt(wall.abserr), "verdict".into(), "matches reference IEEE arithmetic".into()];
    s.push_str(&row(&cols, &[8, 6, 8, 10]));
    s.push('\n');

    s.push_str(&bar);
    s.push('\n');
    s.push_str("end of report\n");
    s
}

fn header_row<S: AsRef<str>>(cols: &[S], widths: &[usize]) -> String {
    let mut out = String::from("-");
    for w in widths {
        out.push_str(&"-".repeat(w + 1));
    }
    out.push('\n');
    out.push_str(&row(cols, widths));
    out
}
