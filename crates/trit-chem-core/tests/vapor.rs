//! Boiling and sublimation of a salt and a metal.
//!
//! The transition temperature is not a table entry: it is where the model's own
//! `(condensed -> gas)` Gibbs energy crosses zero, because the gas side of that
//! reaction sits at the 1 bar standard state. The latent heats are checked
//! against measured values.

use trit_chem_core::equilibrium;
use trit_chem_core::reaction::{Reaction, Term};
use trit_chem_core::rules;
use trit_chem_core::species::Db;
use trit_chem_core::state::State;
use trit_chem_core::thermo::R;

fn transition(db: &Db, from: &str, to: &str) -> Reaction {
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

/// Temperature where `dG` of a condensed-to-gas pair crosses zero, by bisection
/// (`dG` falls monotonically with T there).
fn boiling_point(db: &Db, r: &Reaction) -> f64 {
    let (mut lo, mut hi) = (400.0, 6000.0);
    assert!(
        r.dg_standard(db, lo) > 0.0 && r.dg_standard(db, hi) < 0.0,
        "no crossing in the bracketed range"
    );
    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        if r.dg_standard(db, mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo + hi) / 2.0
}

fn p_sat(db: &Db, r: &Reaction, t: f64) -> f64 {
    (-r.dg_standard(db, t) * 1000.0 / (R * t)).exp()
}

#[test]
fn salt_boils_within_four_percent_of_1738_k() {
    let db = Db::standard();
    let vap = transition(&db, "NaCl(l)", "NaCl(g)");
    assert!(vap.dg_standard(&db, 1700.0) > 0.0);
    assert!(vap.dg_standard(&db, 1850.0) < 0.0);
    let tb = boiling_point(&db, &vap);
    let rel = (tb - 1738.0) / 1738.0;
    assert!(
        rel.abs() < 0.04,
        "model Tb = {tb:.0} K ({:+.1}%)",
        rel * 100.0
    );
    assert!(
        (vap.dh_standard(&db, 1738.0) - 171.0).abs() < 6.0,
        "dH_vap {}",
        vap.dh_standard(&db, 1738.0)
    );
}

#[test]
fn iron_boils_within_five_percent_of_3134_k() {
    // Residual gap: the iron rows' entropy path. NIST fits iron's solid and
    // liquid only in ranges that leave the 298.15 K -> melting integral to the
    // model, and the fitted enthalpy anchor of the vapour fit does not reach
    // down to 298.15 K either; the balance lands a few J/(mol K) high, which is
    // worth more than 100 K on a 3134 K boiling point.
    let db = Db::standard();
    let vap = transition(&db, "Fe(l)", "Fe(g)");
    assert!(vap.dg_standard(&db, 2900.0) > 0.0);
    assert!(vap.dg_standard(&db, 3150.0) < 0.0);
    let tb = boiling_point(&db, &vap);
    let rel = (tb - 3134.0) / 3134.0;
    assert!(
        rel.abs() < 0.05,
        "model Tb = {tb:.0} K ({:+.1}%)",
        rel * 100.0
    );
    assert!((vap.dh_standard(&db, 3134.0) - 349.9).abs() < 10.0);
}

#[test]
fn sublimation_enthalpy_at_298_is_the_measured_one() {
    // dH_sub(298.15) is dHf(gas) - dHf(condensed), both straight from NIST:
    // NaCl 229.7 kJ/mol, Fe 416.3 kJ/mol.
    let db = Db::standard();
    for (solid, gas, want) in [("NaCl(s)", "NaCl(g)", 229.7), ("Fe(s)", "Fe(g)", 416.3)] {
        let sub = transition(&db, solid, gas);
        let dh = sub.dh_standard(&db, 298.15);
        assert!((dh - want).abs() < 2.0, "{gas}: dH_sub {dh} vs {want}");
    }
}

#[test]
fn sublimation_is_fusion_plus_vapourisation() {
    // Structural: the three pairs are the same stoichiometry added two ways, so
    // a mis-wired term shows up here rather than hiding in a tolerance.
    let db = Db::standard();
    for name in ["NaCl", "Fe"] {
        let (s, l, g) = (
            format!("{name}(s)"),
            format!("{name}(l)"),
            format!("{name}(g)"),
        );
        let (fus, vap, sub) = (
            transition(&db, &s, &l),
            transition(&db, &l, &g),
            transition(&db, &s, &g),
        );
        for t in [1200.0, 1800.0, 2600.0] {
            let want = fus.dg_standard(&db, t) + vap.dg_standard(&db, t);
            let got = sub.dg_standard(&db, t);
            assert!((want - got).abs() < 1e-9, "{name} at {t}: {want} vs {got}");
        }
    }
}

#[test]
fn vapour_pressure_rises_through_five_orders_below_the_boiling_point() {
    let db = Db::standard();
    let vap = transition(&db, "NaCl(l)", "NaCl(g)");
    let (lo, hi) = (p_sat(&db, &vap, 1100.0), p_sat(&db, &vap, 1700.0));
    assert!(lo < hi, "{lo} vs {hi}");
    assert!(
        hi > 0.5 && hi < 1.0,
        "1 bar must be reached near the boiling point: {hi}"
    );
    assert!(
        (hi / lo).ln() > 2.0 * std::f64::consts::LN_10,
        "1100->1700 K only moved {lo} -> {hi}"
    );
}

#[test]
fn cascade_follows_the_vapour_pressure_it_predicts() {
    let db = Db::standard();
    let (liq, gas, solid) = (
        db.require("NaCl(l)"),
        db.require("NaCl(g)"),
        db.require("NaCl(s)"),
    );
    let net = rules::phase_transitions(&db, &[(liq, gas), (solid, gas), (solid, liq)]);
    assert_eq!(net.len(), 6);

    // 1 mol of molten salt with a 1000 L headspace at 1900 K: all of it would
    // sit at 0.19 bar, below the vapour pressure there, so nothing holds it
    // back.
    let mut hot = State::new(&db, 1900.0, 0.05, 1000.0);
    hot.add(liq, 1.0);
    let rep = equilibrium::relax(&db, &mut hot, &net, 128, 1e-14);
    assert!(rep.converged, "hot: {rep:?}");
    assert!(hot.n[gas] > 0.99, "hot: gas {:.4}", hot.n[gas]);
    assert!(hot.n[liq] < 1e-6 && hot.n[solid] < 1e-6, "hot: {hot:?}");

    // At 1200 K the headspace saturates first, so evaporation stops at the
    // amount the same reaction's equilibrium constant implies.
    let mut cool = State::new(&db, 1200.0, 0.05, 1000.0);
    cool.add(liq, 1.0);
    let rep = equilibrium::relax(&db, &mut cool, &net, 128, 1e-14);
    assert!(rep.converged, "cool: {rep:?}");
    let vap = transition(&db, "NaCl(l)", "NaCl(g)");
    let want = p_sat(&db, &vap, 1200.0) * 1e5 * (cool.gas_l * 1e-3) / (R * cool.t);
    assert!(
        (cool.n[gas] - want).abs() / want < 0.02,
        "cascade gas {:.5} vs p_sat amount {want:.5}",
        cool.n[gas]
    );
    assert!(cool.n[liq] > 0.9, "cool: {cool:?}");
    assert!(cool.n[solid] < 1e-9, "1200 K is above the melting point");
}
