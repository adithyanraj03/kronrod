//! qpsrt: verified by replicating qage.f's exact bisection bookkeeping
//! (labels 8/10/20) with synthetic errors, and asserting the invariants
//! the driver relies on after every reorder:
//!
//!   * the maintained iord prefix is in strictly descending error order,
//!   * `errmax` equals the true maximum over `elist[1..=last]`,
//!   * `maxerr` is an index whose `elist` entry equals `errmax`.
//!
//! The driver invariant qpsrt depends on (established by labels 8/10) is
//! that after a bisection the parent's slot `maxerr` holds the LARGER
//! child error and the fresh slot `last` holds the smaller one.

use kronrod::qpsrt::qpsrt;

struct Driver {
    limit: usize,
    last: usize,
    elist: Vec<f64>,
    iord: Vec<usize>,
    maxerr: usize,
    errmax: f64,
    nrmax: usize,
}

impl Driver {
    fn new(limit: usize, first_error: f64) -> Self {
        let mut d = Driver {
            limit,
            last: 1,
            elist: vec![0.0; limit + 1],
            iord: vec![0; limit + 1],
            maxerr: 1,
            errmax: first_error,
            nrmax: 1,
        };
        d.elist[1] = first_error;
        d
    }

    /// One bisection step, replicating qage.f labels 8/10 then the qpsrt
    /// call at label 20. `e1`/`e2` are the child errors.
    fn step(&mut self, e1: f64, e2: f64) {
        self.last += 1;
        // Label 8/10 bookkeeping (error values only): the parent slot
        // `maxerr` receives the larger child error, the fresh slot `last`
        // the smaller — exactly as qage.f assigns elist(maxerr)/elist(last).
        if e2 > e1 {
            self.elist[self.maxerr] = e2;
            self.elist[self.last] = e1;
        } else {
            self.elist[self.maxerr] = e1;
            self.elist[self.last] = e2;
        }
        qpsrt(
            self.limit,
            self.last,
            &mut self.maxerr,
            &mut self.errmax,
            &mut self.elist,
            &mut self.iord,
            &mut self.nrmax,
        );
    }

    fn assert_invariants(&self) {
        // errmax is the true maximum over the live list.
        let true_max = self.elist[1..=self.last]
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            (self.errmax - true_max).abs() <= f64::EPSILON * true_max.max(1.0),
            "errmax {} != true max {}",
            self.errmax,
            true_max
        );
        // maxerr points at a slot holding that maximum.
        assert_eq!(self.elist[self.maxerr], self.errmax);
        // The maintained iord prefix is descending.
        let k = if self.last > self.limit / 2 + 2 {
            self.limit + 3 - self.last
        } else {
            self.last
        };
        for i in 1..k {
            assert!(
                self.elist[self.iord[i]] >= self.elist[self.iord[i + 1]],
                "prefix not descending at position {i} (last={})",
                self.last
            );
        }
    }
}

#[test]
fn ordered_prefix_invariant_long_run() {
    // 200 bisections; the driver caps last at limit, so limit must cover
    // the whole run (the final step exercises the capped prefix).
    let mut d = Driver::new(201, 1.0);
    let mut seed: u64 = 0x9E3779B97F4A7C15;
    let mut next = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    for _ in 0..200 {
        // Children strictly smaller than the parent (bisection-like).
        let p = d.errmax;
        let e1 = p * 0.4 * (0.5 + next());
        let e2 = p * 0.4 * (0.5 + next());
        d.step(e1, e2);
        d.assert_invariants();
    }
}

#[test]
fn last_le_2_initialisation_picks_parent_slot() {
    // Driver invariant: after the first bisection the parent slot (index
    // maxerr=1) holds the LARGER child error. qpsrt with last==2 returns
    // iord[nrmax] (nrmax=1) == index 1, so errmax must be elist[1].
    let mut d = Driver::new(50, 0.0);
    d.elist[1] = 0.0;
    d.last = 1;
    d.maxerr = 1;
    d.errmax = 0.0;
    d.nrmax = 1;
    // First bisection of the (zero) parent: children 7.0 (larger) and 3.0.
    d.step(7.0, 3.0);
    assert_eq!(d.last, 2);
    assert_eq!(d.errmax, 7.0);
    assert_eq!(d.maxerr, 1);
    assert_eq!(d.elist[1], 7.0);
    assert_eq!(d.elist[2], 3.0);
    d.assert_invariants();
}

#[test]
fn maxerr_moves_to_previous_second_place() {
    // Build a known list, then bisect the current max and confirm the
    // new max is the previous second-place interval.
    let mut d = Driver::new(50, 0.0);
    d.elist[1] = 0.0;
    d.last = 1;
    d.maxerr = 1;
    d.errmax = 0.0;
    d.nrmax = 1;
    // Step 1: children 5.0 and 2.0 -> elist[1]=5.0 (max), elist[2]=2.0.
    d.step(5.0, 2.0);
    // Step 2: bisect the max (5.0 at index 1) into small children so the
    // previous second place (2.0 at index 2) becomes the new max.
    d.step(1.0, 0.5);
    assert_eq!(d.errmax, 2.0);
    assert_eq!(d.elist[d.maxerr], 2.0);
    d.assert_invariants();
}

#[test]
fn near_limit_prefix_shrinks() {
    // As last approaches the limit, the maintained prefix shrinks to
    // limit+3-last (the qpsrt.f bound); invariants must still hold.
    let limit = 12usize;
    let mut d = Driver::new(limit, 1.0);
    let mut seed: u64 = 0x123456789ABCDEF0;
    let mut next = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    for _ in 0..(limit - 1) {
        let p = d.errmax;
        d.step(p * 0.4 * (0.5 + next()), p * 0.4 * (0.5 + next()));
        d.assert_invariants();
    }
    assert_eq!(d.last, limit);
}
