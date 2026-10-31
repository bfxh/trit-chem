//! Activity-model wiring: Davies must change charged-unsymmetric reactions and
//! must cancel exactly where the ionic-strength correction cancels.

use trit_chem_core::equilibrium;
use trit_chem_core::rules;
use trit_chem_core::species::Db;
use trit_chem_core::state::State;
use trit_chem_core::thermo::{R, T0};

fn build(db: &Db, model: trit_chem_core::ActivityModel) -> State {
    let mut st = State::new(db, T0, 1.0, 1.0);
    st.activity_model = model;
    st.add(db.require("Na+(aq)"), 0.10);
    st.add(db.require("Cl-(aq)"), 0.05);
    st.add(db.require("OH-(aq)"), 0.05);
    st.add(db.require("H+(aq)"), 0.01);
    st.add(db.require("Na(s)"), 0.01);
    st.add(db.require("H2O(l)"), 55.5);
    st
}

#[test]
fn davies_lowers_the_water_reaction_by_four_rt_ln_gamma() {
    let db = Db::standard();
    let ideal = build(&db, trit_chem_core::ActivityModel::Ideal);
    let davies = build(&db, trit_chem_core::ActivityModel::Davies);
    let rxn = rules::metal_water(&db, 11, 1).unwrap();

    let g = davies.activity_coefficient(&db, db.require("Na+(aq)"));
    assert!((g - 0.74).abs() < 0.05, "gamma {g}");
    let expect = 4.0 * R * T0 * g.ln() / 1000.0; // two Na+ and two OH-, all |z| = 1
    let delta = rxn.dg(&db, &davies) - rxn.dg(&db, &ideal);
    assert!(
        (delta - expect).abs() < 0.15,
        "delta {delta:.3} vs expected {expect:.3}"
    );
    assert!(
        delta < -1.5,
        "Davies must make the water reaction more favourable"
    );
}

#[test]
fn davies_cancels_in_the_acid_reaction() {
    let db = Db::standard();
    let ideal = build(&db, trit_chem_core::ActivityModel::Ideal);
    let davies = build(&db, trit_chem_core::ActivityModel::Davies);
    let rxn = rules::acid_metal(&db, 11, 1).unwrap();
    // 2 H+ (reactant, |z|=1) against 2 Na+ (product, |z|=1): gamma^2 / gamma^2.
    let delta = rxn.dg(&db, &davies) - rxn.dg(&db, &ideal);
    assert!(delta.abs() < 1e-9, "delta {delta}");
}

#[test]
fn cascade_reaches_the_same_final_state_under_both_models() {
    let db = Db::standard();
    let run = |model| {
        let mut st = State::new(&db, T0, 1.0, 1.0);
        st.activity_model = model;
        st.add(db.require("Na(s)"), 0.10);
        st.add(db.require("H+(aq)"), 0.05);
        st.add(db.require("Cl-(aq)"), 0.05);
        st.add(db.require("H2O(l)"), 55.5);
        let net = [
            rules::acid_metal(&db, 11, 1).unwrap(),
            rules::metal_water(&db, 11, 1).unwrap(),
            rules::neutralization(&db).unwrap(),
        ];
        let rep = equilibrium::relax(&db, &mut st, &net, 64, 1e-12);
        assert!(rep.converged);
        (
            st.n[db.require("Na+(aq)")],
            st.n[db.require("OH-(aq)")],
            st.n[db.require("H2(g)")],
        )
    };
    let a = run(trit_chem_core::ActivityModel::Ideal);
    let b = run(trit_chem_core::ActivityModel::Davies);
    for (x, y) in [(a.0, b.0), (a.1, b.1), (a.2, b.2)] {
        assert!((x - y).abs() < 1e-9, "{x} vs {y}");
    }
}
