//! Stoichiometry, balance checking, and Gibbs energies of reaction.

use core::fmt;

use crate::formula::{Formula, FormulaError};
use crate::species::Db;
use crate::state::State;
use crate::thermo::R;

/// One stoichiometric coefficient. Negative = reactant, positive = product.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Term {
    /// Species index in the database.
    pub species: usize,
    /// Stoichiometric coefficient.
    pub nu: f64,
}

impl Term {
    /// Convenience constructor.
    pub const fn new(species: usize, nu: f64) -> Self {
        Term { species, nu }
    }
}

/// A balanced candidate reaction.
#[derive(Clone, Debug, PartialEq)]
pub struct Reaction {
    /// Human-readable label, e.g. `"2Na(s) + 2H+(aq) -> 2Na+(aq) + H2(g)"`.
    pub label: String,
    /// Stoichiometric terms.
    pub terms: Vec<Term>,
}

/// Why a reaction could not be built.
#[derive(Clone, Debug, PartialEq)]
pub enum ReactionError {
    /// Atomic symbol could not be resolved.
    Formula(FormulaError),
    /// No terms.
    Empty,
    /// Element not conserved. `residual` = sum(nu * atoms).
    UnbalancedElement {
        /// Atomic number.
        z: u8,
        /// Non-zero residual.
        residual: f64,
    },
    /// Charge not conserved.
    UnbalancedCharge(f64),
    /// A stoichiometric coefficient is not finite or is zero.
    BadCoefficient(f64),
}

impl fmt::Display for ReactionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReactionError::Formula(e) => write!(f, "formula: {e}"),
            ReactionError::Empty => write!(f, "empty reaction"),
            ReactionError::UnbalancedElement { z, residual } => {
                write!(f, "element Z={z} residual {residual}")
            }
            ReactionError::UnbalancedCharge(q) => write!(f, "charge residual {q}"),
            ReactionError::BadCoefficient(v) => write!(f, "bad coefficient {v}"),
        }
    }
}

impl std::error::Error for ReactionError {}

impl Reaction {
    /// Build a reaction, validating coefficients, element balance and charge.
    pub fn new(db: &Db, label: impl Into<String>, terms: Vec<Term>) -> Result<Self, ReactionError> {
        if terms.is_empty() {
            return Err(ReactionError::Empty);
        }
        for t in &terms {
            if !t.nu.is_finite() || t.nu == 0.0 {
                return Err(ReactionError::BadCoefficient(t.nu));
            }
        }
        let mut elements: Vec<(u8, f64)> = Vec::new();
        let mut charge = 0.0;
        for t in &terms {
            let sp = db.get(t.species);
            let f = Formula::parse(sp.formula).map_err(ReactionError::Formula)?;
            charge += f64::from(sp.charge) * t.nu;
            for (z, n) in f.counts {
                match elements.iter_mut().find(|(zz, _)| *zz == z) {
                    Some(slot) => slot.1 += n * t.nu,
                    None => elements.push((z, n * t.nu)),
                }
            }
        }
        if let Some((z, residual)) = elements.iter().find(|(_, r)| r.abs() > 1e-9) {
            return Err(ReactionError::UnbalancedElement {
                z: *z,
                residual: *residual,
            });
        }
        if charge.abs() > 1e-9 {
            return Err(ReactionError::UnbalancedCharge(charge));
        }
        Ok(Reaction {
            label: label.into(),
            terms,
        })
    }

    /// Index of the `k`-th reactant (negative nu) or `None`.
    pub fn reactants(&self) -> impl Iterator<Item = Term> + '_ {
        self.terms.iter().copied().filter(|t| t.nu < 0.0)
    }

    /// Index of the `k`-th product (positive nu).
    pub fn products(&self) -> impl Iterator<Item = Term> + '_ {
        self.terms.iter().copied().filter(|t| t.nu > 0.0)
    }

    /// Standard Gibbs energy of reaction at `t`, kJ per extent unit.
    pub fn dg_standard(&self, db: &Db, t: f64) -> f64 {
        self.terms
            .iter()
            .map(|term| term.nu * db.get(term.species).g(t))
            .sum()
    }

    /// Standard enthalpy of reaction at `t`, kJ per extent unit.
    pub fn dh_standard(&self, db: &Db, t: f64) -> f64 {
        self.terms
            .iter()
            .map(|term| term.nu * db.get(term.species).h(t))
            .sum()
    }

    /// Standard entropy of reaction at `t`, J/(K per extent unit).
    pub fn ds_standard(&self, db: &Db, t: f64) -> f64 {
        self.terms
            .iter()
            .map(|term| term.nu * db.get(term.species).s(t))
            .sum()
    }

    /// `ln K` at `t`.
    pub fn ln_k(&self, db: &Db, t: f64) -> f64 {
        -self.dg_standard(db, t) * 1000.0 / (R * t)
    }

    /// `ln Q` for a given state.
    pub fn ln_q(&self, db: &Db, st: &State) -> f64 {
        let pairs: Vec<(usize, f64)> = self.terms.iter().map(|t| (t.species, t.nu)).collect();
        st.ln_reaction_quotient(db, &pairs)
    }

    /// Actual Gibbs energy of reaction at the state's temperature, kJ/extent.
    pub fn dg(&self, db: &Db, st: &State) -> f64 {
        self.dg_standard(db, st.t) + R * st.t * self.ln_q(db, st) / 1000.0
    }

    /// True when every reactant is available (`n > eps`, or a present solid).
    pub fn reactants_available(&self, db: &Db, st: &State) -> bool {
        self.reactants().all(|t| {
            let sp = db.get(t.species);
            if sp.is_condensed_pure() {
                st.n[t.species] > 1e-12
            } else {
                st.n[t.species] > 1e-15
            }
        })
    }

    /// True when reactants exist and every product has somewhere to live
    /// (a gas product needs headspace, a solute needs liquid volume).
    pub fn phase_feasible(&self, db: &Db, st: &State) -> bool {
        if !self.reactants_available(db, st) {
            return false;
        }
        for t in self.products() {
            match db.get(t.species).phase {
                crate::Phase::Gas if st.gas_l <= 0.0 => return false,
                crate::Phase::Aqueous if st.liquid_l <= 0.0 => return false,
                _ => {}
            }
        }
        true
    }

    /// Largest extent before the first reactant runs out, in mol.
    pub fn max_extent(&self, st: &State) -> f64 {
        self.reactants()
            .map(|t| st.n[t.species] / -t.nu)
            .fold(f64::INFINITY, f64::min)
    }

    /// Apply extent `xi` (mol) to `st`.
    pub fn apply(&self, st: &mut State, xi: f64) {
        for t in &self.terms {
            st.n[t.species] = (st.n[t.species] + t.nu * xi).max(0.0);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermo::T0;

    #[test]
    fn balance_checker_rejects_bad_oxygen() {
        let db = Db::standard();
        let h2 = db.require("H2(g)");
        let o2 = db.require("O2(g)");
        let h2o = db.require("H2O(l)");
        // 2 H2 + O2 -> H2O is not balanced (one O short).
        let err = Reaction::new(
            &db,
            "bad",
            vec![
                Term::new(h2, -2.0),
                Term::new(o2, -1.0),
                Term::new(h2o, 1.0),
            ],
        )
        .unwrap_err();
        assert!(
            matches!(err, ReactionError::UnbalancedElement { .. }),
            "{err:?}"
        );
    }

    #[test]
    fn balance_checker_rejects_bad_charge() {
        let db = Db::standard();
        let na = db.require("Na(s)");
        let nap = db.require("Na+(aq)");
        let err =
            Reaction::new(&db, "bad", vec![Term::new(na, -1.0), Term::new(nap, 1.0)]).unwrap_err();
        assert!(matches!(err, ReactionError::UnbalancedCharge(_)));
    }

    #[test]
    fn hydrogen_combustion_matches_nist() {
        let db = Db::standard();
        let h2 = db.require("H2(g)");
        let o2 = db.require("O2(g)");
        let h2o = db.require("H2O(l)");
        let r = Reaction::new(
            &db,
            "H2 + 1/2 O2 -> H2O(l)",
            vec![
                Term::new(h2, -1.0),
                Term::new(o2, -0.5),
                Term::new(h2o, 1.0),
            ],
        )
        .unwrap();
        let dg = r.dg_standard(&db, T0);
        let dh = r.dh_standard(&db, T0);
        assert!((dg - -237.13).abs() < 0.1, "dG {dg}");
        assert!((dh - -285.83).abs() < 0.1, "dH {dh}");
        // dG and dH must agree with dS.
        let ds = r.ds_standard(&db, T0);
        assert!((dh - T0 * ds / 1000.0 - dg).abs() < 0.15, "{ds}");
    }
}
