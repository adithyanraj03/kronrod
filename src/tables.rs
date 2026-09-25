//! Parsed quadrature tables with structural accessors.
//!
//! The digit strings are the exact source digits from the Netlib QUADPACK
//! evaluation sources ([`crate::tables_gen`]); they are parsed to `f64`
//! exactly once (process-lifetime cache) and never rounded elsewhere.
//!
//! ## Layout (mirrors `dqk*.f`)
//!
//! * `nodes` — `(K+1)/2` positive abscissae, strictly decreasing, the last
//!   entry is `0.0` (the centre).
//! * `k_weights` — the weight of each `nodes` entry; the mirrored pair
//!   `±nodes[i]` shares the weight, the centre weight is used once.
//! * `gauss_weights` — the weights of the positive Gauss nodes at the odd
//!   0-based positions `nodes[1], nodes[3], ...`; for **odd `n`** the
//!   last entry is the **centre weight** (the centre belongs to the
//!   Gauss rule, as in `qk15.f`: `resg = fc*wg((n+1)/2)`).
//!
//! Odd 0-based positions are the positive Gauss nodes, even 0-based
//! positions the additional Kronrod nodes; the centre `nodes[(K+1)/2-1]`
//! is shared by both rules.

use std::sync::OnceLock;

use crate::tables_gen as gen;

/// One (G_n, K_{2n+1}) Gauss-Kronrod pair for the uniform weight on [-1, 1].
#[derive(Debug, Clone)]
pub struct Rule {
    /// QUADPACK `keyf` selector: `1..=6`.
    pub keyf: u8,
    /// Gauss order `n` (7, 10, 15, 20, 25, 30).
    pub n: usize,
    /// Kronrod order `2n+1` (15, 21, 31, 41, 51, 61).
    pub k: usize,
    /// `(K+1)/2` positive abscissae, strictly decreasing, last == 0.0.
    pub nodes: Vec<f64>,
    /// Weight of each `nodes` entry (mirrored pair shares it).
    pub k_weights: Vec<f64>,
    /// Weights of the positive Gauss nodes (at `nodes[1], ...`); for odd
    /// `n` the final entry is the centre weight.
    pub gauss_weights: Vec<f64>,
    /// `true` for the high-precision (10,21) Fullerton/Bell Labs table
    /// (evaluated with 80-digit arithmetic, stored at 33 digits).
    pub high_precision: bool,
}

static RULES: OnceLock<Vec<Rule>> = OnceLock::new();

/// All six rules in `keyf` order (1..=6).
pub fn rules() -> &'static [Rule] {
    RULES.get_or_init(build)
}

/// The QUADPACK `keyf` selectors, in order.
pub const KEYF_ALL: [u8; 6] = [1, 2, 3, 4, 5, 6];

/// The rule for a QUADPACK `keyf` selector (panics for values outside 1..=6).
pub fn rule(keyf: u8) -> &'static Rule {
    rules()
        .iter()
        .find(|r| r.keyf == keyf)
        .unwrap_or_else(|| panic!("keyf must be in 1..=6, got {keyf}"))
}

fn parse<const M: usize>(digits: &[&str; M]) -> Vec<f64> {
    digits
        .iter()
        .map(|s| s.parse::<f64>().unwrap_or_else(|e| panic!("bad table digit string: {s}: {e}")))
        .collect()
}

fn build() -> Vec<Rule> {
    let mut out = Vec::new();
    out.push(Rule {
        keyf: 1,
        n: 7,
        k: 15,
        nodes: parse(&gen::XGK_15),
        k_weights: parse(&gen::WGK_15),
        gauss_weights: parse(&gen::WG_7),
        high_precision: false,
    });
    out.push(Rule {
        keyf: 2,
        n: 10,
        k: 21,
        nodes: parse(&gen::XGK_21),
        k_weights: parse(&gen::WGK_21),
        gauss_weights: parse(&gen::WG_10),
        high_precision: true,
    });
    out.push(Rule {
        keyf: 3,
        n: 15,
        k: 31,
        nodes: parse(&gen::XGK_31),
        k_weights: parse(&gen::WGK_31),
        gauss_weights: parse(&gen::WG_15),
        high_precision: false,
    });
    out.push(Rule {
        keyf: 4,
        n: 20,
        k: 41,
        nodes: parse(&gen::XGK_41),
        k_weights: parse(&gen::WGK_41),
        gauss_weights: parse(&gen::WG_20),
        high_precision: false,
    });
    out.push(Rule {
        keyf: 5,
        n: 25,
        k: 51,
        nodes: parse(&gen::XGK_51),
        k_weights: parse(&gen::WGK_51),
        gauss_weights: parse(&gen::WG_25),
        high_precision: false,
    });
    out.push(Rule {
        keyf: 6,
        n: 30,
        k: 61,
        nodes: parse(&gen::XGK_61),
        k_weights: parse(&gen::WGK_61),
        gauss_weights: parse(&gen::WG_30),
        high_precision: false,
    });
    out
}

impl Rule {
    /// 0-based index of the centre entry (always the last one).
    pub fn center_index(&self) -> usize {
        self.nodes.len() - 1
    }

    /// Whether the nested Gauss rule includes the centre (true for odd `n`;
    /// the centre is then the final entry of `gauss_weights`).
    pub fn gauss_has_center(&self) -> bool {
        self.n % 2 == 1
    }

    /// Number of positive (mirrored) Gauss nodes.
    pub fn gauss_pairs(&self) -> usize {
        if self.gauss_has_center() {
            (self.n - 1) / 2
        } else {
            self.n / 2
        }
    }

    /// Number of positive (mirrored) additional Kronrod nodes.
    pub fn kronrod_pairs(&self) -> usize {
        self.nodes.len() - 1 - self.gauss_pairs()
    }

    /// Total mass of the Kronrod rule: `2·Σ k_weights − k_weights[centre]`.
    pub fn kronrod_mass(&self) -> f64 {
        let c = self.center_index();
        2.0 * self.k_weights.iter().sum::<f64>() - self.k_weights[c]
    }

    /// Total mass of the nested Gauss rule. Mirrored pairs are counted
    /// twice; the centre (odd `n` only) once.
    pub fn gauss_mass(&self) -> f64 {
        let pairs: f64 = self.gauss_weights[..self.gauss_pairs()].iter().sum();
        if self.gauss_has_center() {
            let center = self.gauss_weights[self.gauss_pairs()];
            2.0 * pairs + center
        } else {
            2.0 * pairs
        }
    }

    /// The positive Gauss nodes (odd 0-based positions of `nodes`).
    pub fn gauss_nodes(&self) -> Vec<f64> {
        (0..self.gauss_pairs()).map(|j| self.nodes[2 * j + 1]).collect()
    }

    /// The positive additional Kronrod nodes (even 0-based positions).
    pub fn kronrod_nodes(&self) -> Vec<f64> {
        (0..self.kronrod_pairs()).map(|j| self.nodes[2 * j]).collect()
    }
}
