//! The (G_n, K_{2n+1}) rule evaluator — a faithful port of the
//! `dqk15`/`dqk21`/`dqk31`/`dqk41`/`dqk51`/`dqk61` evaluation and
//! error-estimation block (the executable tail of `dqk21.f`).
//!
//! The evaluation walks the rule in QUADPACK's order: centre value first
//! (into `resk`; and into `resg` as well for odd `n`, where the centre
//! belongs to the Gauss rule), then the mirrored Gauss pairs
//! (`xgk(2j)`, 1-based), then the mirrored additional-Kronrod pairs
//! (`xgk(2j-1)`). The function values at the `K/2` non-centre mirrored
//! pairs are kept for the `resasc` ("residue") accumulation, exactly as
//! in `dqk21.f`.
//!
//! The error estimator is the last six lines of `dqk21.f`, exposed
//! separately as [`estimate_error`] so it can be unit-tested against
//! independently computed constants.

use crate::tables::Rule;

/// Output of one rule application on `[a, b]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QkResult {
    /// `resk * hlgth` — the (2n+1)-point estimate of the integral.
    pub result: f64,
    /// QUADPACK's absolute error estimate (see [`estimate_error`]).
    pub abserr: f64,
    /// `resabs` — scaled sum of `|weight · f|` terms; input to the
    /// `50·epmach` roundoff floor.
    pub resabs: f64,
    /// `resasc` — scaled mean-deviation "residue"; input to the
    /// `min(10, (200·abserr/resasc)^1.5)` scaling.
    pub resasc: f64,
}

/// The trailing estimator block of `dqk21.f`, verbatim:
///
/// ```fortran
/// abserr = dabs((resk-resg)*hlgth)
/// if (resasc ≠ 0 .and. abserr ≠ 0)
///    abserr = resasc * dmin1(10.0, (200.0*abserr/resasc)**1.5)
/// if (resabs > uflow/(50.0*epmach))
///    abserr = dmax1((50.0*epmach)*resabs, abserr)
/// ```
///
/// with `epmach = f64::EPSILON` and `uflow = f64::MIN` (the double
/// `d1mach(4)` / `d1mach(1)` values).
pub fn estimate_error(resg: f64, resk: f64, hlgth: f64, resabs: f64, resasc: f64) -> f64 {
    const EPMACH: f64 = f64::EPSILON;
    const UFLOW: f64 = f64::MIN;
    let mut abserr = ((resk - resg) * hlgth).abs();
    if resasc != 0.0 && abserr != 0.0 {
        abserr = resasc * ((200.0 * abserr / resasc).powf(1.5)).min(10.0);
    }
    if resabs > UFLOW / (50.0 * EPMACH) {
        abserr = (EPMACH * 50.0 * resabs).max(abserr);
    }
    abserr
}

/// Apply the (G_n, K_{2n+1}) rule of `rule` to `f` on `[a, b]`.
///
/// Mirrors `dqk21.f` line-for-line (generalized to all six rules via the
/// shared `dqk*.f` layout):
///
/// * `centr = (a+b)/2`, `hlgth = (b-a)/2`, `dhlgth = |hlgth|`
/// * centre: `resk = wgk(c)·f(centr)`, `resabs = |resk|`
/// * Gauss pairs (1-based `j = 1..n/2`, node `xgk(2j)`, weight `wg(j)`)
/// * Kronrod pairs (1-based `j = 1..n/2`, node `xgk(2j-1)`, weight `wgk(2j-1)`)
/// * `reskh = resk/2`, `resasc = wgk(c)|f(centr) - reskh| + Σ wgk(j)|fv1(j) - reskh| + |fv2(j) - reskh|`
/// * `result = resk·hlgth`; scale `resabs`, `resasc` by `dhlgth`
/// * `abserr = estimate_error(resg, resk, hlgth, resabs, resasc)`
pub fn qk(rule: &Rule, f: &impl Fn(f64) -> f64, a: f64, b: f64) -> QkResult {
    let centr = 0.5 * (a + b);
    let hlgth = 0.5 * (b - a);
    let dhlgth = hlgth.abs();

    let half_k = rule.nodes.len(); // (K+1)/2 entries incl. centre
    let gauss_pairs = rule.gauss_pairs();
    let kronrod_pairs = rule.kronrod_pairs();
    let center = rule.center_index();
    let has_center = rule.gauss_has_center();

    // For odd n the centre belongs to the Gauss rule (qk15.f: resg =
    // fc*wg(4)); for even n the Gauss rule is the mirrored pairs only.
    let fc = f(centr);
    let mut resg = if has_center {
        rule.gauss_weights[gauss_pairs] * fc
    } else {
        0.0
    };
    let mut resk = rule.k_weights[center] * fc;
    let mut resabs = resk.abs();

    // fv1/fv2 (dqk21): values at the half_k-1 mirrored non-centre pairs,
    // stored by 0-based node position.
    let mut fv1 = vec![0.0f64; half_k];
    let mut fv2 = vec![0.0f64; half_k];

    // Gauss pairs: 1-based xgk(2j) -> 0-based nodes[2j-1], weight wg(j) ->
    // gauss_weights[j-1].
    for j in 0..gauss_pairs {
        let idx = 2 * j + 1;
        let absc = hlgth * rule.nodes[idx];
        let fval1 = f(centr - absc);
        let fval2 = f(centr + absc);
        fv1[idx] = fval1;
        fv2[idx] = fval2;
        let fsum = fval1 + fval2;
        resg += rule.gauss_weights[j] * fsum;
        resk += rule.k_weights[idx] * fsum;
        resabs += rule.k_weights[idx] * (fval1.abs() + fval2.abs());
    }

    // Kronrod pairs: 1-based xgk(2j-1) -> 0-based nodes[2j-2], weight
    // wgk(2j-1) -> k_weights[2j-2].
    for j in 0..kronrod_pairs {
        let idx = 2 * j;
        let absc = hlgth * rule.nodes[idx];
        let fval1 = f(centr - absc);
        let fval2 = f(centr + absc);
        fv1[idx] = fval1;
        fv2[idx] = fval2;
        let fsum = fval1 + fval2;
        resk += rule.k_weights[idx] * fsum;
        resabs += rule.k_weights[idx] * (fval1.abs() + fval2.abs());
    }

    // reskh / resasc (dqk21 lines 168-172): the centre term uses wgk(11)
    // (= k_weights[center]); the loop j=1..10 skips the centre entry.
    let reskh = 0.5 * resk;
    let mut resasc = rule.k_weights[center] * (fc - reskh).abs();
    for idx in 0..half_k {
        if idx == center {
            continue;
        }
        resasc += rule.k_weights[idx] * ((fv1[idx] - reskh).abs() + (fv2[idx] - reskh).abs());
    }

    let result = resk * hlgth;
    let resabs = resabs * dhlgth;
    let resasc = resasc * dhlgth;
    let abserr = estimate_error(resg, resk, hlgth, resabs, resasc);

    QkResult {
        result,
        abserr,
        resabs,
        resasc,
    }
}
