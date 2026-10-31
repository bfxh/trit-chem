//! Melting of a salt and a metal: the transition temperature and latent heat
//! come from measured heats of fusion, not from 298 K extrapolations.

use trit_chem_core::equilibrium;
use trit_chem_core::reaction::{Reaction, Term};
use trit_chem_core::rules;
use trit_chem_core::species::Db;
use trit_chem_core::state::State;

fn fusion(db: &Db, solid: &str, liquid: &str) -> Reaction {
    Reaction::new(
        db,
        format!("{solid} -> {liquid}"),
        vec![
            Term::new(db.require(solid), -1.0),
            Term::new(db.require(liquid), 1.0),
        ],
    )
    .unwrap()
}

#[test]
fn nacl_melts_near_1074_k_with_the_measured_latent_heat() {
    let db = Db::standard();
    let r = fusion(&db, "NaCl(s)", "NaCl(l)");
    assert!(r.dg_standard(&db, 1073.0) > 0.0);
    assert!(r.dg_standard(&db, 1075.0) < 0.0);
    let dh = r.dh_standard(&db, 1074.0);
    assert!((dh - 28.0).abs() < 1.5, "dH_fus {dh}");
}

#[test]
fn iron_melts_near_1811_k_with_the_measured_latent_heat() {
    let db = Db::standard();
    let r = fusion(&db, "Fe(s)", "Fe(l)");
    assert!(r.dg_standard(&db, 1810.0) > 0.0);
    assert!(r.dg_standard(&db, 1812.0) < 0.0);
    let dh = r.dh_standard(&db, 1811.15);
    assert!((dh - 13.8).abs() < 1.0, "dH_fus {dh}");
}

#[test]
fn cascade_keeps_the_phase_that_is_stable() {
    let db = Db::standard();
    for (name, tm) in [("NaCl", 1074.0), ("Fe", 1811.15)] {
        let solid = db.require(&format!("{name}(s)"));
        let liquid = db.require(&format!("{name}(l)"));
        let net = rules::fusion_transitions(&db, &[(solid, liquid)]);
        let run = |t: f64| {
            let mut st = State::new(&db, t, 1.0, 0.0);
            st.add(solid, 1.0);
            let rep = equilibrium::relax(&db, &mut st, &net, 16, 1e-14);
            assert!(rep.converged);
            (st.n[solid], st.n[liquid])
        };
        let (s_low, l_low) = run(tm - 50.0);
        assert!(
            s_low > 0.99 && l_low < 1e-9,
            "{name} below: {s_low}/{l_low}"
        );
        let (s_high, l_high) = run(tm + 50.0);
        assert!(
            l_high > 0.99 && s_high < 1e-9,
            "{name} above: {s_high}/{l_high}"
        );
    }
}

#[test]
fn solid_species_do_not_borrow_the_liquid_heat_capacity() {
    // Before the phase filter the solid row also carried the liquid Shomate
    // fits, so above the melting point the "solid" quietly used liquid Cp.
    // A solid's Cp at 50 K below melting must stay near the solid value.
    let db = Db::standard();
    let nacl = db.get(db.require("NaCl(s)"));
    let cp = nacl.cp.cp(1024.0);
    assert!((cp - 64.0).abs() < 5.0, "NaCl(s) Cp(1024) = {cp}");
    let fe = db.get(db.require("Fe(s)"));
    let cp = fe.cp.cp(1761.0);
    assert!((cp - 40.0).abs() < 6.0, "Fe(s) Cp(1761) = {cp}");
}
