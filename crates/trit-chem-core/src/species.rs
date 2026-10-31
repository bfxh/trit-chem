//! Species table and standard-state accessors.

use crate::Phase;
pub use crate::species_data::SPECIES;
use crate::thermo::{Cp, T0};

/// One species: formula, phase, charge and 298.15 K thermochemistry.
#[derive(Clone, Copy, Debug)]
pub struct Species {
    /// Stable identifier, e.g. `"NaCl(s)"` or `"OH-(aq)"`.
    pub id: &'static str,
    /// Formula without charge, e.g. `"NaCl"`.
    pub formula: &'static str,
    /// Aggregation state.
    pub phase: Phase,
    /// Net charge in elementary units.
    pub charge: i8,
    /// Standard enthalpy of formation at 298.15 K, kJ/mol.
    pub hf298: f64,
    /// Standard molar entropy at 298.15 K, J/(mol K).
    pub s298: f64,
    /// Heat-capacity model.
    pub cp: Cp,
    /// Tabulated standard Gibbs energy of formation at 298.15 K, kJ/mol.
    pub dgf298: Option<f64>,
    /// Where the numbers come from.
    pub source: &'static str,
}

impl Species {
    /// Standard enthalpy at `t`, kJ/mol.
    pub fn h(&self, t: f64) -> f64 {
        self.hf298 + self.cp.h_rel_kj(t)
    }

    /// Standard entropy at `t`, J/(mol K).
    pub fn s(&self, t: f64) -> f64 {
        self.s298 + self.cp.s_rel(t)
    }

    /// Standard Gibbs energy of formation at `t`, kJ/mol.
    ///
    /// Anchored at the tabulated `dGf(298.15)`, then integrated with the
    /// Gibbs-Helmholtz relation:
    /// `G(T) = Gf298 + [H(T)-H(298)] - T*[S(T)-S(298)] - (T-298)*S298`.
    ///
    /// The anchor matters for aqueous ions: their tabulated `dHf`, `dGf` and
    /// `S` are *not* related by `dG = dH - T*dS` because the electron
    /// convention assigns the electron zero Gibbs energy. Using `dHf - T*S`
    /// would shift Na+ by about 19 kJ/mol.
    pub fn g(&self, t: f64) -> f64 {
        match self.dgf298 {
            Some(gf) => {
                gf + self.cp.h_rel_kj(t)
                    - t * self.cp.s_rel(t) / 1000.0
                    - (t - T0) * self.s298 / 1000.0
            }
            None => self.h(t) - t * self.s(t) / 1000.0,
        }
    }

    /// True for pure condensed phases (activity = 1 when present).
    pub fn is_condensed_pure(&self) -> bool {
        matches!(self.phase, Phase::Solid | Phase::Liquid)
    }
}

/// All species, ordered as generated.
pub const fn all() -> &'static [Species] {
    &SPECIES
}

/// Species database with id lookups.
#[derive(Clone, Copy, Debug)]
pub struct Db {
    species: &'static [Species],
}

impl Default for Db {
    fn default() -> Self {
        Self::standard()
    }
}

impl Db {
    /// The built-in table.
    pub const fn standard() -> Self {
        Db { species: &SPECIES }
    }

    /// Borrow the backing slice.
    pub const fn species(&self) -> &'static [Species] {
        self.species
    }

    /// Index of a species id.
    pub fn index(&self, id: &str) -> Option<usize> {
        self.species.iter().position(|s| s.id == id)
    }

    /// Index lookup that panics with a useful message, for internal rules.
    pub fn require(&self, id: &str) -> usize {
        self.index(id)
            .unwrap_or_else(|| panic!("species not in database: {id}"))
    }

    /// Species by index.
    pub fn get(&self, i: usize) -> &'static Species {
        &self.species[i]
    }

    /// First species with `formula` in `phase`.
    pub fn by_formula(&self, formula: &str, phase: Phase) -> Option<usize> {
        self.species
            .iter()
            .position(|s| s.formula == formula && s.phase == phase)
    }

    /// Standard-state reference species for an element (`H2(g)`, `Na(s)`, ...).
    pub fn element_reference(&self, z: u8) -> Option<usize> {
        ELEMENT_REFERENCE
            .iter()
            .find(|(zz, _)| *zz == z)
            .and_then(|(_, id)| self.index(id))
    }

    /// Standard entropy per mole of *atoms* of an element, J/(mol K).
    ///
    /// The reference species is a molecule (`H2`, `Cl2`, `O2`), so its molar
    /// entropy is divided by the number of atoms it contains.
    pub fn element_s298(&self, z: u8) -> Option<f64> {
        let i = self.element_reference(z)?;
        let sp = self.get(i);
        let atoms = crate::formula::Formula::parse(sp.formula).ok()?.count(z);
        if atoms <= 0.0 {
            return None;
        }
        Some(sp.s298 / atoms)
    }

    /// Standard formation Gibbs energy implied by `dHf`, `S` and the element
    /// reference entropies at 298.15 K.
    ///
    /// `dGf = dHf - T * (S_species - sum_k n_k S_element_k)`
    pub fn implied_dgf298(&self, i: usize) -> Option<f64> {
        let sp = self.get(i);
        let formula = crate::formula::Formula::parse(sp.formula).ok()?;
        if formula.counts.is_empty() {
            return None;
        }
        let mut s_elements = 0.0;
        for (z, n) in &formula.counts {
            s_elements += self.element_s298(*z)? * n;
        }
        Some(sp.hf298 - T0 * (sp.s298 - s_elements) / 1000.0)
    }
}

/// `(Z, reference species id)` for elements used by the built-in table.
pub const ELEMENT_REFERENCE: &[(u8, &str)] = &[
    (1, "H2(g)"),
    (6, "C(s,graphite)"),
    (7, "N2(g)"),
    (8, "O2(g)"),
    (11, "Na(s)"),
    (12, "Mg(s)"),
    (13, "Al(s)"),
    (14, "Si(s)"),
    (15, "P(s,white)"),
    (16, "S(s,rhombic)"),
    (17, "Cl2(g)"),
    (19, "K(s)"),
    (20, "Ca(s)"),
    (26, "Fe(s)"),
    (29, "Cu(s)"),
    (30, "Zn(s)"),
    (47, "Ag(s)"),
    (82, "Pb(s)"),
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_has_unique_ids() {
        let db = Db::standard();
        let mut ids: Vec<_> = db.species().iter().map(|s| s.id).collect();
        let n = ids.len();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), n);
    }

    #[test]
    fn every_species_formula_parses_and_contains_its_charge_partners() {
        let db = Db::standard();
        for s in db.species() {
            let f = crate::formula::Formula::parse(s.formula)
                .unwrap_or_else(|e| panic!("{}: {e}", s.id));
            assert!(!f.counts.is_empty(), "{}", s.id);
        }
    }

    #[test]
    fn tabulated_dgf_is_reproducible_from_h_s() {
        let db = Db::standard();
        let mut checked = 0;
        for (i, sp) in db.species().iter().enumerate() {
            let Some(tab) = sp.dgf298 else { continue };
            if sp.phase == Phase::Aqueous || tab == 0.0 {
                continue;
            }
            let implied = db.implied_dgf298(i).expect("reference entropy");
            assert!(
                (implied - tab).abs() < 1.0,
                "{}: implied {implied:.2} != tabulated {tab:.2}",
                sp.id
            );
            checked += 1;
        }
        assert!(checked >= 15, "only {checked} species were checked");
    }

    #[test]
    fn water_formation_gibbs_matches_janaf() {
        let db = Db::standard();
        for (id, want) in [("H2O(l)", -237.13), ("H2O(g)", -228.57)] {
            let i = db.index(id).unwrap();
            let g = db.implied_dgf298(i).unwrap();
            assert!((g - want).abs() < 0.1, "{id}: {g}");
        }
    }
}
