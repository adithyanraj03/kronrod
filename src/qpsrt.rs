//! `qpsrt.f` — faithful port of the QUADPACK list-reordering routine.
//!
//! Keeps the first `jupbn = min(last, limit/2 + 2)` error estimates in
//! descending order (via the `iord` permutation), inserts the freshly
//! computed `errmax` (top-down) and `errmin = elist(last)` (bottom-up),
//! and returns the index `maxerr` of the largest error estimate.
//!
//! **1-based semantics.** Like the Fortran source, `elist` and `iord` are
//! 1-indexed: the caller keeps a leading dummy element at index 0, and
//! valid entries are `1..=last`. This keeps the port line-for-line and
//! auditable against `qpsrt.f`.

/// Port of `qpsrt.f`.
///
/// * `limit` — the maximum number of subintervals (the driver's `limit`).
/// * `last` — the current number of subintervals (≥ 1).
/// * `maxerr` — in: index (1-based) of the interval bisected last; out:
///   index of the interval with the largest error estimate.
/// * `errmax` — out: the largest error estimate.
/// * `elist` — 1-indexed error estimates (length ≥ `last` + 1).
/// * `iord` — 1-indexed permutation workspace (length ≥ `last` + 1).
/// * `nrmax` — in/out: position (within the ordered prefix) of the
///   largest error estimate; the driver starts with `nrmax = 1`.
pub fn qpsrt(
    limit: usize,
    last: usize,
    maxerr: &mut usize,
    errmax: &mut f64,
    elist: &mut [f64],
    iord: &mut [usize],
    nrmax: &mut usize,
) {
    // qpsrt.f: if (last.gt.2) go to 10 / iord(1)=1; iord(2)=2; go to 90
    if last <= 2 {
        iord[1] = 1;
        iord[2] = 2;
        *maxerr = iord[*nrmax];
        *errmax = elist[*maxerr];
        return;
    }

    // label 10
    *errmax = elist[*maxerr];
    // if (nrmax.eq.1) go to 30 ; do 20 i=1,ido ...
    if *nrmax != 1 {
        let ido = *nrmax - 1;
        let mut i = 1usize;
        while i <= ido {
            let isucc = iord[*nrmax - 1];
            if *errmax <= elist[isucc] {
                break; // go to 30
            }
            iord[*nrmax] = isucc;
            *nrmax -= 1;
            i += 1;
        }
    }

    // label 30
    let mut jupbn = last;
    if last > limit / 2 + 2 {
        jupbn = limit + 3 - last;
    }
    let errmin = elist[last];
    let jbnd = jupbn - 1;
    let ibeg = *nrmax + 1;
    if ibeg > jbnd {
        // label 50
        iord[jbnd] = *maxerr;
        iord[jupbn] = last;
    } else {
        // do 40 i=ibeg,jbnd
        let mut i = ibeg;
        let mut jumped = false;
        while i <= jbnd {
            let isucc = iord[i];
            if *errmax >= elist[isucc] {
                jumped = true; // go to 60
                break;
            }
            iord[i - 1] = isucc;
            i += 1;
        }
        if !jumped {
            // label 50
            iord[jbnd] = *maxerr;
            iord[jupbn] = last;
        } else {
            // label 60
            iord[i - 1] = *maxerr;
            let mut k = jbnd;
            let mut j = i;
            let mut placed = false;
            while j <= jbnd {
                let isucc = iord[k];
                if errmin < elist[isucc] {
                    placed = true; // go to 80
                    break;
                }
                iord[k + 1] = isucc;
                k -= 1;
                j += 1;
            }
            if !placed {
                iord[i] = last;
            } else {
                // label 80
                iord[k + 1] = last;
            }
        }
    }

    // label 90
    *maxerr = iord[*nrmax];
    *errmax = elist[*maxerr];
}
