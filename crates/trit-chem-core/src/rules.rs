//! Rule-based candidate reactions.
//!
//! These are *rules*, not a reaction table: each function composes
//! stoichiometry from element/ion identities and the engine then decides
//! whether the reaction is spontaneous. Anything the rule cannot express
//! (or whose species is missing from the table) returns an error instead of
//! guessing.

use crate::Phase;
use crate::element;
use crate::formula::Formula;
use crate::reaction::{Reaction, ReactionError, Term};
use crate::species::Db;

/// Most common aqueous cation for a metal, used when a rule needs one.
pub fn dominant_cation(symbol: &str) -> Option<i8> {
    let ox = match symbol {
        "Li" | "Na" | "K" | "Rb" | "Cs" => 1,
        "Be" | "Mg" | "Ca" | "Sr" | "Ba" | "Zn" | "Fe" | "Cu" | "Ni" | "Co" | "Mn" | "Pb" => 2,
        "Al" | "Cr" | "Sc" | "Y" | "La" => 3,
        "Ag" => 1,
        "Sn" => 2,
        _ => return element::by_symbol(symbol)?.lowest_cation(),
    };
    Some(ox)
}

fn solid_element(db: &Db, z: u8) -> Option<usize> {
    let e = element::by_z(z)?;
    db.by_formula(e.symbol, Phase::Solid)
}

fn aqueous_ion(db: &Db, z: u8, charge: i8) -> Option<usize> {
    let e = element::by_z(z)?;
    db.species()
        .iter()
        .position(|s| s.formula == e.symbol && s.phase == Phase::Aqueous && s.charge == charge)
}

fn coef(nu: f64) -> String {
    if (nu - 1.0).abs() < 1e-12 {
        String::new()
    } else if (nu.fract()).abs() < 1e-12 {
        format!("{}", nu as i64)
    } else {
        format!("{nu:.3}")
    }
}

fn describe(db: &Db, terms: &[(usize, f64)]) -> String {
    let left: Vec<String> = terms
        .iter()
        .filter(|(_, nu)| *nu < 0.0)
        .map(|(i, nu)| format!("{}{}", coef(-nu), db.get(*i).id))
        .collect();
    let right: Vec<String> = terms
        .iter()
        .filter(|(_, nu)| *nu > 0.0)
        .map(|(i, nu)| format!("{}{}", coef(*nu), db.get(*i).id))
        .collect();
    format!("{} -> {}", left.join(" + "), right.join(" + "))
}

fn build(db: &Db, terms: Vec<(usize, f64)>) -> Result<Reaction, ReactionError> {
    let label = describe(db, &terms);
    Reaction::new(
        db,
        label,
        terms.into_iter().map(|(i, nu)| Term::new(i, nu)).collect(),
    )
}

/// `2 M(s) + 2n H+(aq) -> 2 M^n+(aq) + n H2(g)`.
pub fn acid_metal(db: &Db, z: u8, n: i8) -> Result<Reaction, ReactionError> {
    let metal = solid_element(db, z).expect("metal solid in table");
    let ion = aqueous_ion(db, z, n).expect("metal ion in table");
    let h_plus = db.require("H+(aq)");
    let h2 = db.require("H2(g)");
    let nf = f64::from(n);
    build(
        db,
        vec![(metal, -2.0), (h_plus, -2.0 * nf), (ion, 2.0), (h2, nf)],
    )
}

/// `2 M(s) + 2n H2O(l) -> 2 M^n+(aq) + 2n OH-(aq) + n H2(g)`.
pub fn metal_water(db: &Db, z: u8, n: i8) -> Result<Reaction, ReactionError> {
    let metal = solid_element(db, z).expect("metal solid in table");
    let ion = aqueous_ion(db, z, n).expect("metal ion in table");
    let water = db.require("H2O(l)");
    let oh = db.require("OH-(aq)");
    let h2 = db.require("H2(g)");
    let nf = f64::from(n);
    build(
        db,
        vec![
            (metal, -2.0),
            (water, -2.0 * nf),
            (ion, 2.0),
            (oh, 2.0 * nf),
            (h2, nf),
        ],
    )
}

/// `H+(aq) + OH-(aq) -> H2O(l)`.
pub fn neutralization(db: &Db) -> Result<Reaction, ReactionError> {
    build(
        db,
        vec![
            (db.require("H+(aq)"), -1.0),
            (db.require("OH-(aq)"), -1.0),
            (db.require("H2O(l)"), 1.0),
        ],
    )
}

/// `CxHyOz + (x + y/4 - z/2) O2 -> x CO2 + y/2 H2O`.
pub fn combustion(db: &Db, fuel: usize) -> Result<Reaction, ReactionError> {
    let f = Formula::parse(db.get(fuel).formula).map_err(ReactionError::Formula)?;
    let x = f.count(6);
    let y = f.count(1);
    let z = f.count(8);
    let o2 = (x + y / 4.0 - z / 2.0).max(0.0);
    build(
        db,
        vec![
            (fuel, -1.0),
            (db.require("O2(g)"), -o2),
            (db.require("CO2(g)"), x),
            (db.require("H2O(l)"), y / 2.0),
        ],
    )
}

/// `a A^m+ + b B^n- -> A_a B_b(s)` with `a*|n| = b*m`.
pub fn precipitation(
    db: &Db,
    cation: usize,
    anion: usize,
    solid: usize,
) -> Result<Reaction, ReactionError> {
    let m = i32::from(db.get(cation).charge);
    let n = -i32::from(db.get(anion).charge);
    if m <= 0 || n <= 0 {
        return Err(ReactionError::UnbalancedCharge(f64::from(m - n)));
    }
    let g = gcd(m, n);
    let a = f64::from(n / g);
    let b = f64::from(m / g);
    build(db, vec![(cation, -a), (anion, -b), (solid, 1.0)])
}

/// `n2 M1(s) + n1 M2^n2+(aq) -> n2 M1^n1+(aq) + n1 M2(s)`.
pub fn displacement(
    db: &Db,
    solid_z: u8,
    ion_z: u8,
    ion_charge: i8,
) -> Result<Reaction, ReactionError> {
    let m1 = solid_element(db, solid_z).expect("solid");
    let ion1 = aqueous_ion(
        db,
        solid_z,
        dominant_cation(element::by_z(solid_z).unwrap().symbol).unwrap_or(1),
    )
    .expect("ion of displacing metal");
    let m2 = solid_element(db, ion_z).expect("solid");
    let ion2 = aqueous_ion(db, ion_z, ion_charge).expect("displaced ion");
    let n1 = f64::from(db.get(ion1).charge);
    let n2 = f64::from(db.get(ion2).charge);
    build(db, vec![(m1, -n2), (ion2, -n1), (ion1, n2), (m2, n1)])
}

fn gcd(a: i32, b: i32) -> i32 {
    if b == 0 { a.abs() } else { gcd(b, a % b) }
}

/// Both directions of each `(a, b)` phase pair.
///
/// The cascade only moves downhill, so melting and freezing (vaporisation and
/// condensation, sublimation and deposition) must both be present; whichever
/// direction has negative dG is the one that runs.
///
/// Pairs are given by species index: `(solid, liquid)` for melting,
/// `(liquid, gas)` for boiling and `(solid, gas)` for sublimation. The gas side
/// sits at the 1 bar standard state, so the dG of a vaporisation pair crosses
/// zero at that species' normal boiling point - which is how the transition
/// temperatures are checked without a separate table of them.
pub fn phase_transitions(db: &Db, pairs: &[(usize, usize)]) -> Vec<Reaction> {
    let mut out = Vec::new();
    for (a, b) in pairs {
        for (from, to) in [(*a, *b), (*b, *a)] {
            if let Ok(r) = build(db, vec![(from, -1.0), (to, 1.0)]) {
                out.push(r);
            }
        }
    }
    out
}

/// The six water phase transitions as one `(solid, liquid)`, `(liquid, gas)` and
/// `(solid, gas)` pair set.
pub fn water_pairs(db: &Db) -> [(usize, usize); 3] {
    let ice = db.require("H2O(s)");
    let liq = db.require("H2O(l)");
    let gas = db.require("H2O(g)");
    [(ice, liq), (liq, gas), (ice, gas)]
}

/// All rules that can be instantiated for `metal_z` with the built-in table.
pub fn candidates_for_metal(db: &Db, metal_z: u8) -> Vec<Reaction> {
    let Some(sym) = element::by_z(metal_z).map(|e| e.symbol) else {
        return Vec::new();
    };
    let Some(n) = dominant_cation(sym) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    if solid_element(db, metal_z).is_some() && aqueous_ion(db, metal_z, n).is_some() {
        if let Ok(r) = acid_metal(db, metal_z, n) {
            out.push(r);
        }
        if let Ok(r) = metal_water(db, metal_z, n) {
            out.push(r);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermo::T0;

    #[test]
    fn phase_transitions_cover_both_directions() {
        let db = Db::standard();
        let (l, g) = (db.require("NaCl(l)"), db.require("NaCl(g)"));
        let net = phase_transitions(&db, &[(l, g)]);
        assert_eq!(net.len(), 2);
        assert_eq!(net[0].label, "NaCl(l) -> NaCl(g)");
        assert_eq!(net[1].label, "NaCl(g) -> NaCl(l)");
        // The downhill direction flips at the boiling point, so both rows are
        // needed: one of them has positive dG at every temperature.
        let vap = transition_dg(&db, l, g, 1500.0);
        let back = transition_dg(&db, g, l, 1500.0);
        assert!(vap > 0.0 && back < 0.0, "{vap} / {back}");
    }

    fn transition_dg(db: &Db, from: usize, to: usize, t: f64) -> f64 {
        let r = build(db, vec![(from, -1.0), (to, 1.0)]).unwrap();
        r.dg_standard(db, t)
    }

    #[test]
    fn water_pairs_are_the_three_transitions() {
        let db = Db::standard();
        let pairs = water_pairs(&db);
        assert_eq!(pairs.len(), 3);
        assert_eq!(phase_transitions(&db, &pairs).len(), 6);
    }

    #[test]
    fn sodium_rules_are_balanced() {
        let db = Db::standard();
        let acid = acid_metal(&db, 11, 1).unwrap();
        let water = metal_water(&db, 11, 1).unwrap();
        assert_eq!(acid.label, "2Na(s) + 2H+(aq) -> 2Na+(aq) + H2(g)");
        assert!(water.dg_standard(&db, T0) < -100.0);
        assert!(acid.dg_standard(&db, T0) < water.dg_standard(&db, T0));
    }

    #[test]
    fn neutralization_and_combustion_values() {
        let db = Db::standard();
        let neut = neutralization(&db).unwrap();
        let dg = neut.dg_standard(&db, T0);
        assert!((dg - -79.9).abs() < 1.0, "{dg}");

        let methane = db.require("CH4(g)");
        let burn = combustion(&db, methane).unwrap();
        let dg = burn.dg_standard(&db, T0);
        // NIST: CH4 + 2 O2 -> CO2 + 2 H2O(l), dG = -818.0 kJ/mol.
        assert!((dg - -818.0).abs() < 3.0, "{dg}");
    }

    #[test]
    fn precipitation_balances_silver_chloride() {
        let db = Db::standard();
        let ag = db.require("Ag+(aq)");
        let cl = db.require("Cl-(aq)");
        let agcl = db.require("AgCl(s)");
        let r = precipitation(&db, ag, cl, agcl).unwrap();
        assert_eq!(r.label, "Ag+(aq) + Cl-(aq) -> AgCl(s)");
    }
}
