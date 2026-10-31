//! Molar amounts, activities and reaction quotients.

use crate::Phase;
use crate::species::Db;
use crate::thermo::{P0_PA, R};

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
}

impl State {
    /// Empty state with room for every species in `db`.
    pub fn new(db: &Db, t: f64, liquid_l: f64, gas_l: f64) -> Self {
        State {
            n: vec![0.0; db.species().len()],
            liquid_l,
            gas_l,
            t,
        }
    }

    /// Add `amount` mol of species `i`.
    pub fn add(&mut self, i: usize, amount: f64) {
        self.n[i] += amount;
    }

    /// Activity of species `i` on the standard-state scale.
    ///
    /// Gas: partial pressure / 1 bar. Solute: molarity / 1 M.
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
            Phase::Aqueous => self.n[i] / self.liquid_l,
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
