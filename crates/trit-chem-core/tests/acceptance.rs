//! F0 acceptance scenario: sodium in hydrochloric acid, then in water.
//!
//! Every number below is recomputed from the built-in table by the engine;
//! the comments give the NIST/CODATA reference values the engine must land on.

use trit_chem_core::equilibrium;
use trit_chem_core::rules;
use trit_chem_core::species::Db;
use trit_chem_core::state::State;
use trit_chem_core::thermo::T0;

fn setup() -> (Db, State) {
    let db = Db::standard();
    let mut st = State::new(&db, T0, 1.0, 1.0);
    st.add(db.require("Na(s)"), 0.10);
    st.add(db.require("H+(aq)"), 0.05);
    st.add(db.require("Cl-(aq)"), 0.05);
    st.add(db.require("H2O(l)"), 55.5);
    (db, st)
}

#[test]
fn standard_gibbs_of_the_two_competing_reactions() {
    let (db, _st) = setup();
    let acid = rules::acid_metal(&db, 11, 1).unwrap();
    let water = rules::metal_water(&db, 11, 1).unwrap();

    let dg_acid = acid.dg_standard(&db, T0);
    let dg_water = water.dg_standard(&db, T0);

    // 2Na + 2H+ -> 2Na+ + H2 : 2 * (-261.91) = -523.82 kJ per extent unit.
    assert!((dg_acid - -523.82).abs() < 0.5, "{dg_acid}");
    // 2Na + 2H2O -> 2Na+ + 2OH- + H2 : 2 * (-182.02) = -364.04.
    assert!((dg_water - -364.04).abs() < 0.5, "{dg_water}");
    // The design document's F0 acceptance claim is "HCl first": more negative
    // per extent, hence ranked first.
    assert!(dg_acid < dg_water);
    // Per mole of sodium, the numbers the design doc got wrong.
    assert!((dg_acid / 2.0 - -261.91).abs() < 0.3, "{}", dg_acid / 2.0);
    assert!((dg_water / 2.0 - -182.02).abs() < 0.3, "{}", dg_water / 2.0);
}

#[test]
fn cascade_consumes_acid_before_water() {
    let (db, mut st) = setup();
    let acid = rules::acid_metal(&db, 11, 1).unwrap();
    let water = rules::metal_water(&db, 11, 1).unwrap();
    let neut = rules::neutralization(&db).unwrap();
    let mass0 = st.mass_g(&db);
    let charge0 = st.net_charge(&db);

    let rep = equilibrium::relax(&db, &mut st, &[acid, water, neut], 64, 1e-12);
    assert!(rep.converged, "{rep:?}");
    assert!(rep.steps.len() >= 2, "{:?}", rep.steps);

    // The first move must be the acid reaction, and it is limited by H+:
    // extent 0.025 consumes 0.05 mol H+ and 0.05 mol Na.
    assert!(rep.steps[0].label.contains("H+(aq)"), "{:?}", rep.steps[0]);
    assert!((rep.steps[0].extent - 0.025).abs() < 1e-9);

    let h_plus = st.n[db.require("H+(aq)")];
    assert!(h_plus < 1e-9, "H+ left {h_plus}");
    assert!(st.n[db.require("Na(s)")] < 1e-9);

    // 0.05 from the acid step + 0.05 from the water step.
    assert!((st.n[db.require("Na+(aq)")] - 0.10).abs() < 1e-9);
    assert!((st.n[db.require("OH-(aq)")] - 0.05).abs() < 1e-9);
    assert!((st.n[db.require("H2(g)")] - 0.05).abs() < 1e-9);

    assert!((st.mass_g(&db) - mass0).abs() < 1e-9, "mass not conserved");
    assert!(
        (st.net_charge(&db) - charge0).abs() < 1e-12,
        "charge not conserved"
    );
}

#[test]
fn a_metal_without_a_rule_variant_reports_instead_of_guessing() {
    let db = Db::standard();
    // Titanium metal is in the element table but its ion is not in the
    // species table, so candidate generation must come back empty.
    assert!(rules::candidates_for_metal(&db, 22).is_empty());
    assert_eq!(rules::candidates_for_metal(&db, 11).len(), 2);
}
