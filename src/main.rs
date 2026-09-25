//! `kronrod` CLI — verification commands.
//!
//! ```text
//! kronrod help              this text
//! kronrod version           version line
//! kronrod demo              short demo (pi/2 + Gauss from scratch)
//! kronrod kats              key KAT summary (structural + estimator anchors)
//! kronrod battery [PATH]    write the verification report (default: assets/report.txt)
//! kronrod dossier [PATH]    write the PDF dossier (default: assets/dossier.pdf)
//! kronrod attest [PATH]     write the attestation (default: assets/attestation.txt)
//! kronrod svg [DIR]         write anim.svg + frames.svg (default: assets)
//! ```

use std::io::Write;
use std::path::Path;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cmd = args.first().map(|s| s.as_str()).unwrap_or("help");
    let out = args.get(1).map(|s| s.as_str());
    let code = match cmd {
        "help" | "-h" | "--help" => {
            print_help();
            0
        }
        "version" | "-V" | "--version" => {
            println!("kronrod {}", kronrod::VERSION);
            0
        }
        "demo" => demo(),
        "kats" => kats(),
        "battery" => {
            let path = out.unwrap_or("assets/report.txt");
            let report = kronrod::battery::render_report();
            write_file(path, report.as_bytes());
            println!("wrote {path} ({} bytes)", report.len());
            0
        }
        "dossier" => {
            let path = out.unwrap_or("assets/dossier.pdf");
            let report = kronrod::battery::render_report();
            let pdf = kronrod::pdf::render_dossier(&report);
            write_file(path, &pdf);
            println!(
                "wrote {path} ({} bytes, sha256 {})",
                pdf.len(),
                kronrod::crypto::to_hex(&kronrod::crypto::sha256(&pdf))
            );
            0
        }
        "attest" => {
            let path = out.unwrap_or("assets/attestation.txt");
            let a = kronrod::attest::attest();
            let text = kronrod::attest::render_attestation(&a);
            write_file(path, text.as_bytes());
            print!("{text}");
            if a.identical {
                0
            } else {
                1
            }
        }
        "svg" => {
            let dir = out.unwrap_or("assets");
            let anim = kronrod::svg::render_anim(16);
            let sheet = kronrod::svg::render_contact_sheet(16, 4);
            let p1 = format!("{dir}/anim.svg");
            let p2 = format!("{dir}/frames.svg");
            write_file(&p1, anim.as_bytes());
            write_file(&p2, sheet.as_bytes());
            println!("wrote {p1} ({} bytes)", anim.len());
            println!("wrote {p2} ({} bytes)", sheet.len());
            0
        }
        other => {
            eprintln!("unknown command: {other}\n");
            print_help();
            2
        }
    };
    std::process::exit(code);
}

fn print_help() {
    println!("kronrod {} - QUADPACK-style adaptive Gauss-Kronrod quadrature, verified", kronrod::VERSION);
    println!();
    println!("commands:");
    println!("  kronrod help              this text");
    println!("  kronrod version           version line");
    println!("  kronrod demo              short demo (pi/2 + Gauss from scratch)");
    println!("  kronrod kats              key KAT summary");
    println!("  kronrod battery [PATH]    write the verification report (assets/report.txt)");
    println!("  kronrod dossier [PATH]    write the PDF dossier (assets/dossier.pdf)");
    println!("  kronrod attest [PATH]     write the attestation (assets/attestation.txt)");
    println!("  kronrod svg [DIR]         write anim.svg + frames.svg (assets)");
}

fn write_file(path: &str, bytes: &[u8]) {
    if let Some(parent) = Path::new(path).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).ok();
        }
    }
    let mut f = std::fs::File::create(path).expect("create output file");
    f.write_all(bytes).expect("write output file");
}

fn demo() -> i32 {
    let f = |x: f64| (1.0 - x * x).sqrt();
    let r = kronrod::qag(&f, 0.0, 1.0, 0.0, 1e-12, 500);
    println!("pi/4 = sqrt(1-x^2) on [0,1]  (QAG (10,21), epsrel=1e-12)");
    println!("  result = {:.17e}", r.result);
    println!("  abserr = {:.17e}", r.abserr);
    println!("  neval  = {}  ier = {}", r.neval, r.ier);
    let (d, exact) = ((r.result - std::f64::consts::PI / 4.0).abs(), std::f64::consts::PI / 4.0);
    let rel = d / exact;
    println!("  rel err = {:.3e}", rel);
    println!();
    let a = kronrod::gauss_legendre_newton(10);
    let b = kronrod::gauss_legendre_golub_welsch(10);
    let md = a
        .iter()
        .zip(&b)
        .map(|((na, _), (nb, _))| (na - nb).abs())
        .fold(0.0f64, f64::max);
    println!("Gauss-10 from scratch: Newton vs Golub-Welsch max node diff = {:.3e}", md);
    println!("  outermost node (Newton)      = {:.17}", a.last().unwrap().0);
    println!("  outermost node (Golub-Welsch)= {:.17}", b.last().unwrap().0);
    println!("  Netlib (G10,K21) table node  = 0.973906528517171720077964012084452");
    println!("  arbitrary-precision oracle   = 0.9739065285171717200779640120844521");
    0
}

fn kats() -> i32 {
    // 1. table structure
    let mut ok = true;
    for rule in kronrod::rules() {
        let mass = (2.0 - rule.kronrod_mass()).abs();
        if mass > 1e-14 {
            ok = false;
        }
        println!(
            "table (G{},K{}): mass err {:.2e}  nodes {}  weights {}  gauss {}",
            rule.n,
            rule.k,
            mass,
            rule.nodes.len(),
            rule.k_weights.len(),
            rule.gauss_weights.len()
        );
    }
    // 2. estimator anchors (independent constants)
    let checks: [(&str, f64, f64); 5] = [
        ("A clamped x10", kronrod::estimate_error(0.0, 1.0, 0.5, 0.3, 0.2), 2.0),
        ("B power 1.5", kronrod::estimate_error(2.0, 2.0000000001, 0.25, 1.0, 4.0), 1.76776717236491370e-13),
        ("C roundoff floor", kronrod::estimate_error(1.0, 1.0 + 1e-17, 1.0, 1.0, 2.0), 1.11022302462515654e-14),
        ("D resasc=0", kronrod::estimate_error(5.0, 5.0, 0.5, 7.0, 0.0), 7.77156117237609578e-14),
        ("E no floor", kronrod::estimate_error(0.0, 1.0, 0.5, 1e-300, 0.2), 2.0),
    ];
    for (name, got, want) in checks {
        if got != want {
            ok = false;
        }
        println!("estimator {name:<16}: {got:.17e}  (KAT {want:.17e})  {}", if got == want { "ok" } else { "FAIL" });
    }
    // 3. QAG polynomial exactness
    //    (10,21) integrates x^2 exactly; the early-exit test fires after the
    //    single 21-point evaluation, so neval must be exactly 21.
    let r = kronrod::qag(&|x| x * x, 0.0, 1.0, 0.0, 1e-13, 100);
    let d = (r.result - 1.0 / 3.0).abs();
    if r.ier != 0 || r.neval != 21 || d > 1e-15 {
        ok = false;
    }
    println!(
        "qag x^2 on [0,1]: result {:.17e}  ier {}  neval {}  (exact 1/3, expect ier=0 neval=21)  {}",
        r.result,
        r.ier,
        r.neval,
        if r.ier == 0 && r.neval == 21 && d <= 1e-15 { "ok" } else { "FAIL" }
    );
    println!();
    println!("summary: {}", if ok { "ALL KATS PASS" } else { "KAT FAILURES" });
    if ok { 0 } else { 1 }
}
