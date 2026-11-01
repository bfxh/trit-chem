//! Water phase transitions: melting, boiling and the triple point, computed by
//! the same downhill cascade used for chemical reactions.

use trit_chem_core::equilibrium;
use trit_chem_core::reaction::Reaction;
use trit_chem_core::rules;
use trit_chem_core::species::Db;
use trit_chem_core::state::State;
use trit_chem_core::thermo::R;

fn relax_at(db: &Db, t: f64) -> State {
    let net = rules::phase_transitions(db, &rules::water_pairs(db));
    let mut st = State::new(db, t, 0.018, 10.0);
    st.add(db.require("H2O(l)"), 1.0);
    let rep = equilibrium::relax(db, &mut st, &net, 128, 1e-14);
    assert!(rep.converged, "{t} K: {rep:?}");
    st
}

fn p_bar(db: &Db, st: &State) -> f64 {
    st.n[db.require("H2O(g)")] * R * st.t / (st.gas_l * 1e-3) / 1e5
}

fn transition(db: &Db, from: &str, to: &str) -> Reaction {
    use trit_chem_core::reaction::Term;
    Reaction::new(
        db,
        format!("{from} -> {to}"),
        vec![
            Term::new(db.require(from), -1.0),
            Term::new(db.require(to), 1.0),
        ],
    )
    .unwrap()
}

#[test]
fn ice_is_stable_below_the_melting_point() {
    let db = Db::standard();
    let cold = relax_at(&db, 263.15);
    assert!(
        cold.n[db.require("H2O(l)")] < 1e-9,
        "liquid survived below 0 C"
    );
    assert!(cold.n[db.require("H2O(s)")] > 0.99);
    // Vapour pressure over ice at -10 C is about 0.0026 bar.
    let p = p_bar(&db, &cold);
    assert!((p - 0.0026).abs() < 0.0008, "{p}");
}

#[test]
fn liquid_is_stable_above_the_melting_point() {
    let db = Db::standard();
    let warm = relax_at(&db, 283.15);
    assert!(
        warm.n[db.require("H2O(s)")] < 1e-9,
        "ice survived above 0 C"
    );
    assert!(warm.n[db.require("H2O(l)")] > 0.98);
}

#[test]
fn melting_point_is_bracketed_within_a_degree() {
    let db = Db::standard();
    assert!(relax_at(&db, 273.0).n[db.require("H2O(s)")] > 0.9);
    assert!(relax_at(&db, 273.3).n[db.require("H2O(l)")] > 0.9);
}

#[test]
fn saturation_pressure_matches_steam_tables() {
    let db = Db::standard();
    // (T/K, measured p_sat / bar)
    for (t, want) in [
        (283.15, 0.01228),
        (298.15, 0.03169),
        (323.15, 0.12352),
        (353.15, 0.47144),
        (373.15, 1.01325),
    ] {
        let st = relax_at(&db, t);
        let p = p_bar(&db, &st);
        let rel = (p - want).abs() / want;
        assert!(rel < 0.06, "{t} K: {p} vs {want} ({:.1}%)", rel * 100.0);
    }
}

#[test]
fn triple_point_pressure() {
    let db = Db::standard();
    let st = relax_at(&db, 273.15);
    let p = p_bar(&db, &st);
    // 6.11 mbar.
    assert!((p - 0.0061).abs() < 0.0015, "{p}");
}

#[test]
fn transition_enthalpies() {
    let db = Db::standard();
    let fus = transition(&db, "H2O(s)", "H2O(l)");
    let dh = fus.dh_standard(&db, 273.15);
    assert!((dh - 6.01).abs() < 1.2, "dH_fus {dh}");

    let vap = transition(&db, "H2O(l)", "H2O(g)");
    let dh = vap.dh_standard(&db, 373.15);
    assert!((dh - 40.66).abs() < 2.0, "dH_vap {dh}");

    // The transition Gibbs energies change sign at the transition temperatures.
    // Melting is uphill below the melting point and downhill above it.
    assert!(fus.dg_standard(&db, 272.0) > 0.0);
    assert!(fus.dg_standard(&db, 274.0) < 0.0);
    assert!(vap.dg_standard(&db, 372.0) > 0.0);
    assert!(vap.dg_standard(&db, 374.0) < 0.0);
}
