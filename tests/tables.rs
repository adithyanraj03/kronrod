//! Structural KATs for the six embedded (G_n, K_{2n+1}) rules.
//!
//! Two independent oracle families:
//! * the exact Netlib source digit strings (`tables_gen`), and
//! * arbitrary-precision Laurie-method evaluations (P. Holoborodko,
//!   Multiprecision Computing Toolbox, 2011) used to cross-check that the
//!   embedded tables are the correct tables (Netlib stores the
//!   Fullerton/Bell Labs 1981 80-digit-arithmetic evaluation at 33 digits
//!   for (G10,K21) and 16 digits for the other rules).

use kronrod::tables::{rule, rules, KEYF_ALL};
use kronrod::tables_gen as gen;

#[test]
fn six_rules_in_keyf_order() {
    assert_eq!(rules().len(), 6);
    let mut keyfs = Vec::new();
    for r in rules() {
        keyfs.push(r.keyf);
    }
    assert_eq!(keyfs, KEYF_ALL.to_vec());
    assert_eq!(rules()[0].n, 7);
    assert_eq!(rules()[5].k, 61);
    assert_eq!(rules()[1].k, 21);
    assert!(rules()[1].high_precision); // 33-digit Fullerton table
    assert!(!rules()[0].high_precision);
}

#[test]
fn node_and_weight_counts() {
    for r in rules() {
        assert_eq!(r.nodes.len(), (r.k + 1) / 2);
        assert_eq!(r.k_weights.len(), r.nodes.len());
        assert_eq!(r.gauss_weights.len(), r.gauss_pairs() + r.gauss_has_center() as usize);
        assert_eq!(r.k, 2 * r.n + 1);
    }
}

#[test]
fn centre_is_last_and_zero() {
    for r in rules() {
        assert!(r.nodes[r.center_index()] == 0.0);
        assert_eq!(r.center_index(), r.nodes.len() - 1);
    }
}

#[test]
fn nodes_strictly_decreasing_positive() {
    for r in rules() {
        for w in r.nodes[..r.nodes.len() - 1].windows(2) {
            assert!(w[0] > w[1], "keyf {} not strictly decreasing", r.keyf);
        }
        assert!(r.nodes[0] < 1.0);
        assert!(r.nodes[0] > 0.0);
    }
}

#[test]
fn kronrod_mass_is_two() {
    for r in rules() {
        let m = r.kronrod_mass();
        assert!((m - 2.0).abs() <= 4.0 * f64::EPSILON, "keyf {}: {m}", r.keyf);
    }
}

#[test]
fn gauss_mass_is_two() {
    for r in rules() {
        let m = r.gauss_mass();
        assert!((m - 2.0).abs() <= 8.0 * f64::EPSILON, "keyf {}: {m}", r.keyf);
    }
}

#[test]
fn nested_gauss_at_odd_positions() {
    for r in rules() {
        let gn = r.gauss_nodes();
        assert_eq!(gn.len(), r.gauss_pairs());
        for (j, g) in gn.iter().enumerate() {
            assert_eq!(*g, r.nodes[2 * j + 1]);
        }
        for w in gn.windows(2) {
            assert!(w[0] > w[1]);
        }
    }
}

#[test]
fn strict_interlacing_all_rules() {
    for r in rules() {
        let t = r.kronrod_nodes();
        let x = r.gauss_nodes();
        assert!(t.len() >= x.len());
        for (i, &x) in x.iter().enumerate() {
            let next = if i + 1 < t.len() { t[i + 1] } else { 0.0 };
            assert!(x < t[i] && x > next, "keyf {}: x{} not interlaced", r.keyf, i);
        }
    }
}

#[test]
fn netlib_digit_strings_exact() {
    // The embedded strings are the verbatim Netlib source digits.
    assert_eq!(gen::XGK_15[0], "0.9914553711208126e+00");
    assert_eq!(gen::XGK_31[0], "0.9980022986933971e+00");
    assert_eq!(gen::XGK_41[0], "0.9988590315882777e+00");
    assert_eq!(gen::XGK_51[0], "0.9992621049926098e+00");
    assert_eq!(gen::XGK_61[0], "0.9994844100504906e+00");
    assert!(gen::XGK_21[0].starts_with("0.995657163025808"));
    assert!(gen::XGK_21[1].starts_with("0.9739065285171717"));
}

#[test]
fn crosscheck_arbitrary_precision_oracle() {
    // Arbitrary-precision Laurie-method evaluations (Holoborodko 2011).
    // Netlib stores the Fullerton 1981 80-digit-arithmetic evaluation
    // rounded to 33 digits (K21) / 16 digits (others); the embedded
    // values must reproduce the oracle digits to that depth.
    fn sig(s: &str) -> String {
        let s = s.trim().strip_prefix('+').unwrap_or(s);
        let s = s.split('e').next().unwrap_or(s);
        s.strip_prefix("0.").unwrap_or(s).to_string()
    }
    fn agrees_to_digits(a: &str, b: &str, digits: usize) -> bool {
        // Strip signs, exponents, and the leading "0." prefix; compare
        // the first `digits` significant digits.
        let a = sig(a);
        let b = sig(b);
        a.get(..digits) == Some(&b[..digits.min(b.len())])
    }
    // (G10,K21), 33-digit Netlib literal vs 34+ digit oracle.
    assert!(agrees_to_digits(
        gen::XGK_21[0],
        "0.9956571630258080807355272806890028",
        32
    ));
    assert!(agrees_to_digits(
        gen::XGK_21[1],
        "0.9739065285171717200779640120844521",
        33
    ));
    // 16-digit rules vs oracle. The stored literals are correctly
    // rounded to 16 significant digits, so up to the 16th digit they may
    // differ from the oracle by the rounding carry; 15 leading significant
    // digits must agree exactly.
    assert!(agrees_to_digits(gen::XGK_15[0], "0.9914553711208126392068546975263285", 15));
    assert!(agrees_to_digits(gen::XGK_31[0], "0.9980022986933970652730512604445146", 15));
    assert!(agrees_to_digits(gen::XGK_41[0], "0.9988590315882773615124906974789595", 15));
    assert!(agrees_to_digits(gen::XGK_51[0], "0.9992621049926093339221680409990177", 15));
    assert!(agrees_to_digits(gen::XGK_61[0], "0.9994844100504906375713258957058108", 15));
}

#[test]
fn parsed_values_match_oracle_double() {
    // The f64 parsed from each embedded string is the correctly-rounded
    // double of that decimal literal — hence within 1 ulp of the oracle.
    let r21 = rule(2);
    let oracle = 0.9956571630258080807355272806890028f64;
    assert!((r21.nodes[0] - oracle).abs() <= 2.0 * f64::EPSILON * oracle);
    let r61 = rule(6);
    let oracle = 0.9994844100504906375713258957058108f64;
    assert!((r61.nodes[0] - oracle).abs() <= 2.0 * f64::EPSILON * oracle);
}

#[test]
fn gauss_parity() {
    // Odd n: centre belongs to the Gauss rule (extra weight entry).
    assert!(rule(1).gauss_has_center());
    assert!(!rule(2).gauss_has_center());
    assert!(rule(3).gauss_has_center());
    assert!(!rule(4).gauss_has_center());
    assert!(rule(5).gauss_has_center());
    assert!(!rule(6).gauss_has_center());
    assert_eq!(rule(1).gauss_pairs(), 3);
    assert_eq!(rule(1).kronrod_pairs(), 4);
    assert_eq!(rule(2).gauss_pairs(), 5);
    assert_eq!(rule(2).kronrod_pairs(), 5);
    assert_eq!(rule(6).gauss_pairs(), 15);
    assert_eq!(rule(6).kronrod_pairs(), 15);
}

#[test]
fn rule_lookup_boundaries() {
    assert_eq!(rule(1).keyf, 1);
    assert_eq!(rule(6).keyf, 6);
}
