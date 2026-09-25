//! `qage.f` — faithful port of the QUADPACK adaptive quadrature driver
//! (double-precision semantics: `epmach = f64::EPSILON`,
//! `uflow = f64::MIN`, `ier = 6` threshold `max(50·epmach, 0.5e-28)` as in
//! the double source `dqage.f`).
//!
//! ## `ier` codes (QUADPACK convention)
//!
//! * `0` — successful completion.
//! * `1` — `limit` subintervals reached (roundoff-free, but the requested
//!   accuracy may not have been achieved).
//! * `2` — roundoff error detected; the requested accuracy may not have
//!   been achieved.
//! * `3` — very bad integrand behaviour at a point of the integration
//!   range.
//! * `6` — invalid input (both tolerances too small).
//!
//! ## `neval` accounting
//!
//! The driver counts bisections internally; at termination it is scaled
//! exactly as `qage.f` does: `neval = (10·keyf+1)·(2n+1)` for `keyf ≥ 2`
//! and `neval = 30n + 15` for `keyf = 1`, where `n` is the bisection count
//! — i.e. the true number of integrand evaluations.

use crate::qk::{qk, QkResult};
use crate::qpsrt::qpsrt;
use crate::tables::rule as rule_for;

/// Result of an adaptive quadrature run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QagResult {
    /// The integral estimate (sum of the interval approximations).
    pub result: f64,
    /// The accumulated absolute error estimate.
    pub abserr: f64,
    /// QUADPACK `ier` status code (0, 1, 2, 3, 6).
    pub ier: i32,
    /// Total number of integrand evaluations (QUADPACK `neval`).
    pub neval: usize,
    /// Number of subintervals at termination (QUADPACK `last`).
    pub last: usize,
}

/// One bisection step, for tracing/visualization.
#[derive(Debug, Clone, Copy)]
pub struct TraceStep {
    /// `last` at this step (the index of the newly created right interval).
    pub last: usize,
    /// Left subinterval of the bisected interval.
    pub a1: f64,
    pub b1: f64,
    /// Right subinterval of the bisected interval.
    pub a2: f64,
    pub b2: f64,
    /// Combined approximation `area1 + area2`.
    pub area12: f64,
    /// Combined error estimate `error1 + error2`.
    pub erro12: f64,
    /// The error estimate `errmax` of the interval that was bisected.
    pub errmax: f64,
    /// The accumulated error `errsum` after this step.
    pub errsum: f64,
}

/// Standard QAG: the (G10, K21) rule (`keyf = 2`), as in QUADPACK `qag.f`.
pub fn qag(
    f: &impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    epsabs: f64,
    epsrel: f64,
    limit: usize,
) -> QagResult {
    qage_with_rule(f, a, b, epsabs, epsrel, limit, 2, None)
}

/// General adaptive driver for any of the six rules.
///
/// `key` is the QUADPACK selector clamped exactly as `qage.f` does
/// (`key ≤ 0 → 1`, `key ≥ 7 → 6`). `trace` optionally receives one
/// [`TraceStep`] per bisection (in order).
pub fn qage_with_rule(
    f: &impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    epsabs: f64,
    epsrel: f64,
    limit: usize,
    key: u8,
    mut trace: Option<&mut Vec<TraceStep>>,
) -> QagResult {
    const EPMACH: f64 = f64::EPSILON;
    const UFLOW: f64 = f64::MIN;
    assert!(limit >= 1, "limit >= 1");

    // qage.f lines 194-206: initialization + parameter validity test.
    let mut ier: i32 = 0;
    let mut neval: usize = 0;
    let mut last: usize = 0;
    let mut result = 0.0;
    let mut abserr = 0.0;

    let mut alist = vec![0.0f64; limit + 1];
    let mut blist = vec![0.0f64; limit + 1];
    let mut elist = vec![0.0f64; limit + 1];
    let mut rlist = vec![0.0f64; limit + 1];
    let mut iord = vec![0usize; limit + 1];
    alist[1] = a;
    blist[1] = b;
    rlist[1] = 0.0;
    elist[1] = 0.0;
    iord[1] = 0;

    if epsabs <= 0.0 && epsrel < (50.0 * EPMACH).max(0.5e-28) {
        ier = 6;
    }
    if ier == 6 {
        return QagResult {
            result,
            abserr,
            ier,
            neval,
            last,
        };
    }

    // qage.f lines 211-224: first approximation on [a, b].
    let mut keyf = key;
    if keyf == 0 {
        keyf = 1;
    }
    if keyf >= 7 {
        keyf = 6;
    }
    let rule = rule_for(keyf);
    let first: QkResult = qk(rule, f, a, b);
    let defabs_first = first.resabs;
    // qage.f passes the qk `resasc` output into its local `resabs`.
    let resasc_global = first.resasc;
    last = 1;
    rlist[1] = first.result;
    elist[1] = first.abserr;
    iord[1] = 1;

    // qage.f lines 228-233: accuracy test on the first approximation.
    let mut errbnd = epsabs.max(epsrel * first.result.abs());
    if first.abserr <= 50.0 * EPMACH * defabs_first && first.abserr > errbnd {
        ier = 2;
    }
    if limit == 1 {
        ier = 1;
    }
    if ier != 0
        || (first.abserr <= errbnd && first.abserr != resasc_global)
        || first.abserr == 0.0
    {
        result = first.result;
        abserr = first.abserr;
        return finalize(result, abserr, ier, neval, keyf, last);
    }

    // qage.f lines 239-245: initialization of the loop state.
    let mut errmax = first.abserr;
    let mut maxerr = 1usize;
    let mut area = first.result;
    let mut errsum = first.abserr;
    let mut nrmax = 1usize;
    let mut iroff1 = 0usize;
    let mut iroff2 = 0usize;

    // qage.f lines 250-329: main loop (1-based `last = 2..limit`).
    // NOTE: the loop variable is `next` so that the outer `last` (the
    // subinterval count used for the final result sum, qage.f lines
    // 334-337) tracks the Fortran do-loop index exactly.
    'outer: for next in 2..=limit {
        last = next;
        let a1 = alist[maxerr];
        let b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        let a2 = b1;
        let b2 = blist[maxerr];

        let q1: QkResult = qk(rule, f, a1, b1);
        let q2: QkResult = qk(rule, f, a2, b2);

        // qage.f lines 274-282.
        neval += 1;
        let area12 = q1.result + q2.result;
        let erro12 = q1.abserr + q2.abserr;
        errsum += erro12 - errmax;
        area += area12 - rlist[maxerr];
        // qage.f line 279: skip the roundoff counters when the estimator
        // already equals the resasc term (defab1 = qk's resasc output).
        let skip_roundoff = q1.resasc == q1.abserr || q2.resasc == q2.abserr;
        if !skip_roundoff {
            if (rlist[maxerr] - area12).abs() <= 0.1e-04 * area12.abs() && erro12 >= 0.99 * errmax {
                iroff1 += 1;
            }
            if last > 10 && erro12 > errmax {
                iroff2 += 1;
            }
        }
        rlist[maxerr] = q1.result;
        rlist[last] = q2.result;
        errbnd = epsabs.max(epsrel * area.abs());
        if errsum > errbnd {
            // qage.f lines 291-302: error flags.
            if iroff1 >= 6 || iroff2 >= 20 {
                ier = 2;
            }
            if last == limit {
                ier = 1;
            }
            if a1.abs().max(b2.abs()) <= (10.0 + 100.0 * EPMACH) * (a2.abs() + 1.0e4 * UFLOW) {
                ier = 3;
            }
        }
        // qage.f lines 306-319: append the new intervals (swap when
        // error2 > error1 so that elist(maxerr) stays the larger one).
        if q2.abserr > q1.abserr {
            alist[maxerr] = a2;
            alist[last] = a1;
            blist[last] = b1;
            rlist[maxerr] = q2.result;
            rlist[last] = q1.result;
            elist[maxerr] = q2.abserr;
            elist[last] = q1.abserr;
        } else {
            alist[last] = a2;
            blist[maxerr] = b1;
            blist[last] = b2;
            elist[maxerr] = q1.abserr;
            elist[last] = q2.abserr;
        }

        if let Some(t) = &mut trace {
            t.push(TraceStep {
                last,
                a1,
                b1,
                a2,
                b2,
                area12,
                erro12,
                errmax,
                errsum,
            });
        }

        // qage.f line 326: keep the error list ordered.
        qpsrt(limit, last, &mut maxerr, &mut errmax, &mut elist, &mut iord, &mut nrmax);

        // qage.f line 328: jump out of the do-loop.
        if ier != 0 || errsum <= errbnd {
            break 'outer;
        }
    }

    // qage.f lines 334-338: final result.
    let mut res = 0.0;
    for k in 1..=last {
        res += rlist[k];
    }
    result = res;
    abserr = errsum;

    finalize(result, abserr, ier, neval, keyf, last)
}

/// qage.f lines 339-340: scale the bisection count to the true number of
/// integrand evaluations, and pack the final result.
fn finalize(result: f64, abserr: f64, ier: i32, bisections: usize, keyf: u8, last: usize) -> QagResult {
    let neval = if keyf != 1 {
        (10 * keyf as usize + 1) * (2 * bisections + 1)
    } else {
        30 * bisections + 15
    };
    QagResult {
        result,
        abserr,
        ier,
        neval,
        last,
    }
}
