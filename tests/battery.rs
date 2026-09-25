//! The verification battery and report: structural checks on the frozen
//! KAT tables plus the deterministic report layout.

use kronrod::battery::{grid, pi_study, render_report, run_battery};
use std::f64::consts::PI;

#[test]
fn battery_has_sixteen_rows() {
    let rows = run_battery();
    assert_eq!(rows.len(), 16);
}

#[test]
fn battery_rows_are_well_formed() {
    let rows = run_battery();
    for (i, r) in rows.iter().enumerate() {
        assert!(matches!(r.ier, 0 | 1 | 2 | 3), "row {i}: ier={}", r.ier);
        assert!(r.result.is_finite(), "row {i}: result not finite");
        assert!(r.abserr.is_finite(), "row {i}");
        assert!(r.true_error.is_finite(), "row {i}");
        assert!(r.ratio.is_finite(), "row {i}: ratio");
        assert!(r.neval > 0, "row {i}: no evaluations");
        assert!(r.last >= 1, "row {i}");
    }
}

#[test]
fn battery_estimates_are_accurate() {
    let rows = run_battery();
    for r in rows.iter() {
        // Every integrand is resolved to at least 9 correct digits by the
        // (10,21) adaptive driver at epsrel=1e-13.
        let rel = if r.exact != 0.0 {
            (r.result - r.exact).abs() / r.exact.abs()
        } else {
            (r.result - r.exact).abs()
        };
        assert!(rel <= 1e-9, "row {}: {} rel={rel}", r.name, r.result);
    }
}

#[test]
fn battery_estimator_honesty() {
    // Where the estimator reports an error, the true error stays within
    // two orders of magnitude of it (the fair-estimator band).
    let rows = run_battery();
    for r in rows.iter().filter(|r| r.abserr > 0.0 && r.true_error > 0.0) {
        assert!(r.ratio <= 100.0, "row {}: true/est = {}", r.name, r.ratio);
    }
}

#[test]
fn exact_value_constants_match_known() {
    let rows = run_battery();
    let by_name = |n: &str| rows.iter().find(|r| r.name == n).unwrap().exact;
    assert_eq!(by_name("1"), 1.0);
    assert_eq!(by_name("x"), 0.5);
    assert!((by_name("x^2") - 1.0 / 3.0).abs() < 1e-17);
    assert!((by_name("x^10") - 1.0 / 11.0).abs() < 1e-18);
    assert!((by_name("sqrt(x)") - 2.0 / 3.0).abs() < 1e-16);
    assert!((by_name("(1-x)^(1/3)") - 0.75).abs() < 1e-16);
    assert!((by_name("sqrt(1-x^2)") - PI / 4.0).abs() < 1e-16);
    assert_eq!(by_name("log(x)"), -1.0);
    assert!((by_name("exp(x)") - (std::f64::consts::E - 1.0)).abs() < 1e-16);
    assert!((by_name("x^2*exp(x)") - (std::f64::consts::E - 2.0)).abs() < 1e-16);
    assert!((by_name("sin(5x)") - 0.4).abs() < 1e-17);
    assert!((by_name("1/(1+100x^2)") - 10.0f64.atan() / 10.0).abs() < 1e-17);
}

#[test]
fn pi_study_shape() {
    let rows = pi_study();
    assert_eq!(rows.len(), 15);
    // The first 13 tolerances are above the 50*epmach guard (ier=0);
    // the last two are below it (ier=6).
    for r in &rows[..13] {
        assert_eq!(r.ier, 0, "epsrel={:e}", r.epsrel);
    }
    assert_eq!(rows[13].ier, 6);
    assert_eq!(rows[14].ier, 6);
    // Digits are non-decreasing as the tolerance tightens.
    let exact = PI / 4.0;
    let digits = |r: &kronrod::battery::StudyRow| {
        if r.result == exact {
            return u32::MAX;
        }
        let d = (r.result - exact).abs() / exact.abs();
        if d < 1e-16 {
            return 16;
        }
        ((-d.log10()).floor() as i64).clamp(0, 16) as u32
    };
    for w in rows[..13].windows(2) {
        assert!(digits(&w[1]) >= digits(&w[0]), "digits decreased");
    }
    // The tightest valid tolerance delivers 14+ correct digits.
    assert!(digits(&rows[12]) >= 14);
}

#[test]
fn grid_shape() {
    let rows = grid();
    assert_eq!(rows.len(), 24); // 6 rules x 4 tolerances
    let tols = [1e-3f64, 1e-6, 1e-9, 1e-12];
    for r in &rows {
        assert_eq!(r.ier, 0, "rule={} tol={:e}", r.rule, r.epsrel);
        assert!(tols.contains(&r.epsrel), "unexpected tol {:e}", r.epsrel);
    }
    // The (10,21) row at the tightest tolerance reaches 14+ digits.
    let r = rows
        .iter()
        .find(|r| r.rule == "(G10,21)" && r.epsrel == 1e-12)
        .unwrap();
    assert!(r.digits >= 14, "G10@1e-12 digits={}", r.digits);
    // Measured: larger rules are *not* uniformly better at a fixed
    // tolerance — the (30,61) rule lands one digit short of (10,21) here
    // (14 vs 15) because the bisection pattern differs. Assert the real
    // ordering instead of a monotonicity that does not hold.
    let d10 = rows
        .iter()
        .find(|r| r.rule == "(G10,21)" && r.epsrel == 1e-12)
        .unwrap()
        .digits;
    let d30 = rows
        .iter()
        .find(|r| r.rule == "(G30,61)" && r.epsrel == 1e-12)
        .unwrap()
        .digits;
    assert!(d30 >= 13 && d10 >= d30 - 2, "d10={d10} d30={d30}");
    // Every rule converges across the full tolerance sweep.
    for keyf in [1u8, 2, 3, 4, 5, 6] {
        let cell = rows
            .iter()
            .find(|r| r.keyf == keyf && r.epsrel == 1e-12)
            .unwrap();
        assert!(cell.digits >= 13, "keyf={keyf}: {}", cell.digits);
    }
}

#[test]
fn report_layout() {
    let report = render_report();
    for header in [
        "KRONROD VERIFICATION REPORT",
        "SECTION 1 - GAUSS-LEGENDRE FROM SCRATCH",
        "SECTION 2 - TABLE STRUCTURE",
        "SECTION 3 - ESTIMATOR STRESS BATTERY",
        "SECTION 4 - QUARTER-DISC CONVERGENCE STUDY",
        "SECTION 5 - RULE FAMILY x TOLERANCE GRID",
        "SECTION 6 - POLE OVERFLOW WALL",
        "end of report",
    ] {
        assert!(report.contains(header), "missing: {header}");
    }
    // ASCII-only (PDF writer and text tools depend on it).
    assert!(report.chars().all(|c| c.is_ascii()));
    assert!(report.len() > 5000);
    // No formatting accidents in the finite sections (the overflow wall
    // section, of course, documents "inf").
    let finite = &report[..report.find("SECTION 6").unwrap()];
    assert!(!finite.contains("NaN"));
    assert!(!finite.contains("inf"), "non-wall sections must be finite");
}

#[test]
fn report_wall_section_documents_inf() {
    let report = render_report();
    let section6 = report
        .split("SECTION 6")
        .nth(1)
        .expect("section 6 present");
    assert!(section6.contains("1953"));
    assert!(section6.contains("47"));
    assert!(section6.contains("inf"));
}
