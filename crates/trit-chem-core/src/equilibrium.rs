//! Thermodynamically downhill reaction cascade.
//!
//! Given a network of candidate reactions, repeatedly drive any reaction with
//! negative actual dG to its equilibrium extent, limited by the first reactant
//! that runs out. Each move lowers total Gibbs energy, so the sweep sequence is
//! monotone and terminates.
//!
//! This is deliberately *not* a full Gibbs minimisation (White-Johnson-Dantzig
//! style element-potential solve). It is path dependent, only moves forward,
//! and therefore needs the reverse direction supplied as its own reaction. For
//! real-time games that trade-off is the point: bounded work per tick and no
//! stiff solver. See docs/DESIGN.md.

use crate::reaction::Reaction;
use crate::species::Db;
use crate::state::State;

/// One applied extent.
#[derive(Clone, Debug)]
pub struct Step {
    /// Reaction label.
    pub label: String,
    /// Extent in mol.
    pub extent: f64,
    /// Actual dG before the move, kJ/mol extent (live quotient).
    pub dg_before: f64,
    /// Standard dG at the system temperature, for ranking.
    pub dg_standard: f64,
}

/// Outcome of a relaxation run.
#[derive(Clone, Debug)]
pub struct Report {
    /// Applied moves in order.
    pub steps: Vec<Step>,
    /// Sweeps performed.
    pub sweeps: usize,
    /// True when a sweep moved less than the tolerance.
    pub converged: bool,
}

impl Report {
    /// Total extent across all moves.
    pub fn total_extent(&self) -> f64 {
        self.steps.iter().map(|s| s.extent).sum()
    }
}

/// Relax `st` over `network` until no reaction can move more than `tol` mol.
pub fn relax(db: &Db, st: &mut State, network: &[Reaction], max_sweeps: usize, tol: f64) -> Report {
    let mut steps = Vec::new();
    let mut sweeps = 0;
    let mut converged = false;
    for _ in 0..max_sweeps {
        sweeps += 1;
        let mut moved = 0.0_f64;
        for rxn in network {
            if !rxn.phase_feasible(db, st) {
                continue;
            }
            let dg_before = rxn.dg(db, st);
            if dg_before >= 0.0 {
                continue;
            }
            let xi_max = rxn.max_extent(st);
            if xi_max <= 0.0 || !xi_max.is_finite() {
                continue;
            }
            let xi = solve_extent(db, st, rxn, xi_max);
            if xi > tol {
                rxn.apply(st, xi);
                moved = moved.max(xi);
                steps.push(Step {
                    label: rxn.label.clone(),
                    extent: xi,
                    dg_before,
                    dg_standard: rxn.dg_standard(db, st.t),
                });
            }
        }
        if moved <= tol {
            converged = true;
            break;
        }
    }
    Report {
        steps,
        sweeps,
        converged,
    }
}

/// Extent that brings `rxn` to dG = 0, capped at `xi_max`.
fn solve_extent(db: &Db, st: &State, rxn: &Reaction, xi_max: f64) -> f64 {
    let f = |xi: f64| -> f64 {
        let mut probe = st.clone();
        rxn.apply(&mut probe, xi);
        rxn.dg(db, &probe)
    };
    if f(0.0) >= 0.0 {
        return 0.0;
    }
    if f(xi_max) <= 0.0 {
        return xi_max;
    }
    // Q is not guaranteed monotone once activity coefficients depend on ionic
    // strength, so find the first sign change on a coarse grid and bisect only
    // inside it. 32 probes per call; a tick with three reactions costs ~100
    // state clones, which is fine at 60 Hz but is the first thing to optimise
    // if the network grows.
    const GRID: usize = 32;
    let mut lo = 0.0_f64;
    for k in 1..=GRID {
        let hi = xi_max * (k as f64) / (GRID as f64);
        if f(hi) >= 0.0 {
            let (mut a, mut b) = (lo, hi);
            for _ in 0..64 {
                let mid = 0.5 * (a + b);
                if f(mid) > 0.0 {
                    b = mid;
                } else {
                    a = mid;
                }
            }
            return 0.5 * (a + b);
        }
        lo = hi;
    }
    xi_max
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reaction::Term;
    use crate::thermo::T0;

    #[test]
    fn complete_combustion_runs_to_limiting_reagent() {
        let db = Db::standard();
        let mut st = State::new(&db, T0, 1.0, 10.0);
        let h2 = db.require("H2(g)");
        let o2 = db.require("O2(g)");
        let h2o = db.require("H2O(l)");
        st.add(h2, 0.2);
        st.add(o2, 0.05); // limiting: 0.05 mol O2 consumes 0.1 mol H2
        let rxn = Reaction::new(
            &db,
            "H2 + 1/2 O2 -> H2O(l)",
            vec![
                Term::new(h2, -1.0),
                Term::new(o2, -0.5),
                Term::new(h2o, 1.0),
            ],
        )
        .unwrap();
        let rep = relax(&db, &mut st, &[rxn], 16, 1e-12);
        assert!(rep.converged, "{rep:?}");
        assert!((st.n[o2]).abs() < 1e-12, "O2 left {}", st.n[o2]);
        assert!((st.n[h2] - 0.1).abs() < 1e-9, "H2 left {}", st.n[h2]);
        assert!((st.n[h2o] - 0.1).abs() < 1e-9, "H2O {}", st.n[h2o]);
    }

    #[test]
    fn reaction_stops_at_equilibrium_not_at_completion() {
        // N2O4 <=> 2 NO2 has a moderate K; start from pure N2O4.
        let db = Db::standard();
        // Not in the M0 table: fall back to a weak-acid-like check with the
        // water self-ionisation reaction, whose K is 1e-14.
        let mut st = State::new(&db, T0, 1.0, 0.0);
        let h2o = db.require("H2O(l)");
        let hp = db.require("H+(aq)");
        let oh = db.require("OH-(aq)");
        st.add(h2o, 55.5);
        let rxn = Reaction::new(
            &db,
            "H2O -> H+ + OH-",
            vec![Term::new(h2o, -1.0), Term::new(hp, 1.0), Term::new(oh, 1.0)],
        )
        .unwrap();
        let _ = relax(&db, &mut st, &[rxn], 64, 1e-18);
        let kw = st.n[hp] * st.n[oh] / (1.0 * 1.0);
        assert!((kw.log10() + 14.0).abs() < 0.3, "Kw = {kw:e}");
    }
}
