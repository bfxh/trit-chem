//! Mass-action kinetics with Arrhenius rate constants.
//!
//! Rates are in mol/(L s). Gas-phase networks are handled through the ideal
//! gas partial pressures used by [`crate::State`]. The integrator is classical
//! RK4 with a fixed step: it is adequate for the demonstration networks and
//! **not** stiffness-safe (a combustion mechanism will need an implicit method;
//! DESIGN.md lists this as an open gap, not a silent approximation).

use crate::Phase;
use crate::reaction::Reaction;
use crate::species::Db;
use crate::state::State;

/// `k = A * exp(-Ea / RT)`.
#[derive(Clone, Copy, Debug)]
pub struct Arrhenius {
    /// Pre-exponential factor in the units of the reaction's rate law.
    pub a: f64,
    /// Activation energy, kJ/mol.
    pub ea_kj: f64,
}

impl Arrhenius {
    /// Rate constant at `t`.
    pub fn k(&self, t: f64) -> f64 {
        self.a * (-self.ea_kj * 1000.0 / (crate::thermo::R * t)).exp()
    }
}

/// One reversible (or forward-only) kinetic step.
#[derive(Clone, Debug)]
pub struct KineticStep {
    /// Stoichiometry.
    pub reaction: Reaction,
    /// Forward rate parameters.
    pub forward: Arrhenius,
    /// Reverse rate parameters, if the step is reversible.
    pub reverse: Option<Arrhenius>,
    /// Volume basis used to convert mol/(L s) into mol/s.
    pub basis: Phase,
}

impl KineticStep {
    /// Net forward rate in mol/(L s) for `st`.
    pub fn rate(&self, db: &Db, st: &State) -> f64 {
        let mut rate = self.forward.k(st.t) * mass_action(db, st, &self.reaction, true);
        if let Some(rev) = self.reverse {
            rate -= rev.k(st.t) * mass_action(db, st, &self.reaction, false);
        }
        rate
    }

    /// Volume used to scale the rate into mol/s.
    fn volume(&self, st: &State) -> f64 {
        match self.basis {
            Phase::Gas => st.gas_l,
            _ => st.liquid_l,
        }
    }
}

fn mass_action(db: &Db, st: &State, rxn: &Reaction, forward: bool) -> f64 {
    let mut product = 1.0;
    for t in &rxn.terms {
        let exponent = if forward { -t.nu } else { t.nu };
        if exponent > 0.0 {
            product *= st.activity(db, t.species).max(0.0).powf(exponent);
        }
    }
    product
}

/// Time derivative of every amount, mol/s.
pub fn derivatives(db: &Db, st: &State, steps: &[KineticStep]) -> Vec<f64> {
    let mut dn = vec![0.0; st.n.len()];
    for step in steps {
        let v = step.volume(st);
        if v <= 0.0 {
            continue;
        }
        let rate = step.rate(db, st);
        for t in &step.reaction.terms {
            dn[t.species] += t.nu * rate * v;
        }
    }
    dn
}

/// Integration failure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum KineticsError {
    /// A species amount became negative beyond numerical noise.
    NegativeAmount,
    /// A non-finite rate appeared.
    NonFinite,
}

/// Integrate with fixed-step RK4 for `duration` seconds at fixed temperature.
pub fn rk4(
    db: &Db,
    st: &mut State,
    steps: &[KineticStep],
    dt: f64,
    duration: f64,
) -> Result<usize, KineticsError> {
    let mut count = 0;
    let mut time = 0.0_f64;
    while time < duration - 1e-12 {
        let h = dt.min(duration - time);
        let k1 = derivatives(db, st, steps);
        let mut s2 = st.clone();
        apply_derivative(&mut s2, &k1, h / 2.0);
        let k2 = derivatives(db, &s2, steps);
        let mut s3 = st.clone();
        apply_derivative(&mut s3, &k2, h / 2.0);
        let k3 = derivatives(db, &s3, steps);
        let mut s4 = st.clone();
        apply_derivative(&mut s4, &k3, h);
        let k4 = derivatives(db, &s4, steps);
        for i in 0..st.n.len() {
            let delta = h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
            if !delta.is_finite() {
                return Err(KineticsError::NonFinite);
            }
            st.n[i] += delta;
            if st.n[i] < -1e-9 {
                return Err(KineticsError::NegativeAmount);
            }
            st.n[i] = st.n[i].max(0.0);
        }
        time += h;
        count += 1;
    }
    Ok(count)
}

fn apply_derivative(st: &mut State, dn: &[f64], h: f64) {
    for (ni, d) in st.n.iter_mut().zip(dn) {
        *ni = (*ni + d * h).max(0.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reaction::Term;
    use crate::thermo::T0;

    #[test]
    fn arrhenius_matches_hand_value() {
        let a = Arrhenius {
            a: 1.0e10,
            ea_kj: 50.0,
        };
        let k = a.k(500.0);
        let expect = 1.0e10 * (-50_000.0_f64 / (8.314_462_618_153_24 * 500.0)).exp();
        assert!((k - expect).abs() / expect < 1e-12);
    }

    #[test]
    fn constant_activity_gives_linear_consumption() {
        // H2O(l) -> H2O(g): the liquid has unit activity, so this is
        // zeroth order in amount and consumes exactly k*t mol over time.
        let db = Db::standard();
        let mut st = State::new(&db, T0, 1.0, 1.0);
        let liq = db.require("H2O(l)");
        let gas = db.require("H2O(g)");
        st.add(liq, 1.0);
        let rxn = Reaction::new(
            &db,
            "H2O(l) -> H2O(g)",
            vec![Term::new(liq, -1.0), Term::new(gas, 1.0)],
        )
        .unwrap();
        let step = KineticStep {
            reaction: rxn,
            forward: Arrhenius { a: 1.0, ea_kj: 0.0 },
            reverse: None,
            basis: Phase::Liquid,
        };
        let n = rk4(&db, &mut st, &[step], 1e-3, 0.5).unwrap();
        assert!(n >= 499, "{n}");
        assert!((st.n[liq] - 0.5).abs() < 5e-3, "{}", st.n[liq]);
        assert!((st.n[gas] - 0.5).abs() < 5e-3, "{}", st.n[gas]);
    }
}
