//! Molar amounts, activities and reaction quotients.

use crate::Phase;
use crate::species::Db;
use crate::thermo::{P0_PA, R};

/// Debye-Huckel A parameter at 25 C (dimensionless, mol^-1/2 L^1/2).
///
/// True value: 0.5085 at 25 C, 0.4918 at 0 C, 0.6086 at 100 C. M1 keeps the
/// 25 C value for all temperatures; the resulting error in log10(gamma) is
/// under 2% between 0 and 50 C and about 20% at 100 C. Documented, not silent.
pub const DEBYE_HUCKEL_A: f64 = 0.509;

/// How solute activities are computed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ActivityModel {
    /// `a = c / 1 M`; no ionic-strength correction.
    Ideal,
    /// Davies equation, `log10(gamma) = -A z^2 (sqrt(I)/(1+sqrt(I)) - 0.3 I)`.
    ///
    /// Accurate to a few percent up to I ~ 0.1 M, approximate to ~0.5 M, and
    /// unsuitable above that (Pitzer or SIT is the next step).
    #[default]
    Davies,
}

/// A closed system: amounts in mol, one liquid volume in litres, one headspace
/// volume in litres, at a uniform temperature.
///
/// Gases are assumed ideal and share the headspace; solutes share the liquid.
/// Pure condensed phases have unit activity while present.
#[derive(Clone, Debug)]
pub struct State {
    /// Amounts in mol, indexed like the species table.
    pub n: Vec<f64>,
    /// Liquid volume, L.
    pub liquid_l: f64,
    /// Gas headspace volume, L (set to 0 for a sealed liquid-only system).
    pub gas_l: f64,
    /// Temperature, K.
    pub t: f64,
    /// Activity model for aqueous solutes.
    pub activity_model: ActivityModel,
}

impl State {
    /// Empty state with room for every species in `db`, using the Davies model.
    pub fn new(db: &Db, t: f64, liquid_l: f64, gas_l: f64) -> Self {
        State {
            n: vec![0.0; db.species().len()],
            liquid_l,
            gas_l,
            t,
            activity_model: ActivityModel::Davies,
        }
    }

    /// Empty state with the ideal (concentration) activity model.
    pub fn ideal(db: &Db, t: f64, liquid_l: f64, gas_l: f64) -> Self {
        let mut st = Self::new(db, t, liquid_l, gas_l);
        st.activity_model = ActivityModel::Ideal;
        st
    }

    /// Add `amount` mol of species `i`.
    pub fn add(&mut self, i: usize, amount: f64) {
        self.n[i] += amount;
    }

    /// Ionic strength `I = 0.5 * sum c_i z_i^2` over aqueous species, mol/L.
    pub fn ionic_strength(&self, db: &Db) -> f64 {
        if self.liquid_l <= 0.0 {
            return 0.0;
        }
        let mut sum = 0.0;
        for (i, n) in self.n.iter().enumerate() {
            let sp = db.get(i);
            if sp.phase == Phase::Aqueous && sp.charge != 0 {
                let z = f64::from(sp.charge);
                sum += (n / self.liquid_l) * z * z;
            }
        }
        0.5 * sum
    }

    /// Activity coefficient of species `i`.
    ///
    /// Ions follow Davies; everything else is 1. Neutral solutes get no
    /// salting-out term yet (M1 gap).
    pub fn activity_coefficient(&self, db: &Db, i: usize) -> f64 {
        let sp = db.get(i);
        if sp.phase != Phase::Aqueous
            || sp.charge == 0
            || self.activity_model == ActivityModel::Ideal
        {
            return 1.0;
        }
        let i_strength = self.ionic_strength(db);
        if i_strength <= 0.0 {
            return 1.0;
        }
        let z = f64::from(sp.charge);
        let root = i_strength.sqrt();
        let log10_g = -DEBYE_HUCKEL_A * z * z * (root / (1.0 + root) - 0.3 * i_strength);
        10.0_f64.powf(log10_g)
    }

    /// Activity of species `i` on the standard-state scale.
    ///
    /// Gas: partial pressure / 1 bar. Solute: molarity * gamma / 1 M.
    /// Pure condensed: 1 while present, otherwise 0.
    pub fn activity(&self, db: &Db, i: usize) -> f64 {
        let sp = db.get(i);
        match sp.phase {
            Phase::Gas => {
                if self.gas_l <= 0.0 {
                    return 0.0;
                }
                let p = self.n[i] * R * self.t / (self.gas_l * 1e-3);
                p / P0_PA
            }
            Phase::Aqueous => (self.n[i] / self.liquid_l) * self.activity_coefficient(db, i),
            Phase::Solid | Phase::Liquid => {
                if self.n[i] > 1e-12 {
                    1.0
                } else {
                    0.0
                }
            }
            Phase::Unknown => 0.0,
        }
    }

    /// `ln Q` for the given `(species, nu)` list, with activities floored at
    /// 1e-300 so the logarithm stays finite.
    pub fn ln_reaction_quotient(&self, db: &Db, terms: &[(usize, f64)]) -> f64 {
        let mut ln_q = 0.0;
        for (i, nu) in terms {
            let a = self.activity(db, *i).max(1e-300);
            ln_q += nu * a.ln();
        }
        ln_q
    }

    /// Total mass in grams, from species formulas.
    pub fn mass_g(&self, db: &Db) -> f64 {
        self.n
            .iter()
            .enumerate()
            .map(|(i, n)| {
                crate::formula::Formula::parse(db.get(i).formula)
                    .map_or(0.0, |f| f.molar_mass() * n)
            })
            .sum()
    }

    /// Net charge in mol of elementary charge.
    pub fn net_charge(&self, db: &Db) -> f64 {
        self.n
            .iter()
            .enumerate()
            .map(|(i, n)| db.get(i).charge as f64 * n)
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermo::T0;

    fn nacl(db: &Db, mol_per_l: f64) -> State {
        let mut st = State::new(db, T0, 1.0, 0.0);
        let na = db.require("Na+(aq)");
        let cl = db.require("Cl-(aq)");
        st.add(na, mol_per_l);
        st.add(cl, mol_per_l);
        st
    }

    #[test]
    fn ideal_model_never_corrects() {
        let db = Db::standard();
        let st = State::ideal(&db, T0, 1.0, 0.0);
        let i = db.require("Na+(aq)");
        let mut st = st;
        st.add(i, 1.0);
        st.add(db.require("Cl-(aq)"), 1.0);
        assert_eq!(st.activity_coefficient(&db, i), 1.0);
        assert!((st.ionic_strength(&db) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn davies_reproduces_measured_nacl_activity_coefficients() {
        let db = Db::standard();
        // Experimental mean ionic activity coefficients of NaCl at 25 C:
        // 0.001 M -> 0.965, 0.01 M -> 0.902, 0.1 M -> 0.778.
        for (c, want) in [(0.001, 0.965), (0.01, 0.902), (0.1, 0.778)] {
            let st = nacl(&db, c);
            let i_strength = st.ionic_strength(&db);
            assert!((i_strength - c).abs() < 1e-12, "I = {i_strength}");
            let g = st.activity_coefficient(&db, db.require("Na+(aq)"));
            assert!((g - want).abs() < 0.012, "{c} M: gamma {g} vs {want}");
        }
    }

    #[test]
    fn activity_includes_the_coefficient() {
        let db = Db::standard();
        let st = nacl(&db, 0.1);
        let na = db.require("Na+(aq)");
        let a = st.activity(&db, na);
        let g = st.activity_coefficient(&db, na);
        assert!((a - 0.1 * g).abs() < 1e-12);
        assert!(a < 0.1, "Davies must lower the activity of a 1:1 salt");
    }

    #[test]
    fn pure_water_and_solids_are_unaffected() {
        let db = Db::standard();
        let mut st = nacl(&db, 0.1);
        st.add(db.require("H2O(l)"), 55.5);
        st.add(db.require("NaCl(s)"), 0.01);
        assert_eq!(st.activity_coefficient(&db, db.require("H2O(l)")), 1.0);
        assert_eq!(st.activity_coefficient(&db, db.require("NaCl(s)")), 1.0);
    }
}
