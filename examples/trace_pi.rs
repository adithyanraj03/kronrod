// Debug: trace the QAG bisection for the quarter-disc integral.
use std::f64::consts::PI;

fn main() {
    let f = |x: f64| (1.0 - x * x).sqrt();
    let mut trace: Vec<kronrod::TraceStep> = Vec::new();
    let r = kronrod::qage_with_rule(&f, 0.0, 1.0, 0.0, 1e-12, 500, 2, Some(&mut trace));
    let exact = PI / 4.0;
    println!(
        "quarter-disc: result={:.17e} abserr={:.17e} ier={} neval={} last={} rel={:.3e}",
        r.result,
        r.abserr,
        r.ier,
        r.neval,
        r.last,
        (r.result - exact).abs() / exact
    );
    for (i, s) in trace.iter().enumerate() {
        if i < 3 || i >= trace.len() - 3 {
            println!(
                "{:3} last={} [a1={:.12} b1={:.12}] a12={:.8e} e12={:.8e} emax={:.8e} esum={:.8e}",
                i, s.last, s.a1, s.b1, s.area12, s.erro12, s.errmax, s.errsum
            );
        }
    }
    println!("... ({} trace steps total)", trace.len());

    // Overflow-wall probe: pole integrand, faithful IEEE behaviour.
    let pole = |x: f64| (1.0 - x * x).recip().sqrt();
    let w = kronrod::qag(&pole, 0.0, 1.0, 0.0, 1e-13, 1000);
    println!(
        "pole wall: ier={} neval={} last={} result={:.6e} abserr={:.6e}",
        w.ier, w.neval, w.last, w.result, w.abserr
    );
}
