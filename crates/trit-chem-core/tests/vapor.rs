//! Boiling and sublimation of a salt and a metal.
//!
//! The transition temperature is not a table entry: it is where the model's own
//! `(condensed -> gas)` Gibbs energy crosses zero, because the gas side of that
//! reaction sits at the 1 bar standard state. The latent heats are checked
//! against measured values.

use trit_chem_core::equilibrium;
use trit_chem_core::formula::Formula;
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

/// Moles of element `z` held by a state, counted from the species formulas.
fn atoms(db: &Db, st: &State, z: u8) -> f64 {
    let mut sum = 0.0;
    for (i, n) in st.n.iter().enumerate() {
        let f = Formula::parse(db.get(i).formula).expect("formula");
        sum += f.count(z) * n;
    }
    sum
}

fn both(db: &Db, condensed: &str, gas: &str) -> (Reaction, Reaction) {
    (
        transition(db, condensed, gas),
        transition(db, gas, condensed),
    )
}

#[test]
fn vapour_pressure_increases_monotonically_up_to_the_validated_window() {
    // Clausius-Clapeyron: dS of a condensed-to-gas pair is positive, so p_sat
    // must rise at every step. Swept to 2800 K on purpose: above the last
    // fitted range the polynomial runs free and NaCl(s) (whose only fit ends at
    // 1073 K) starts *falling* - measured 77.1 bar at 2972 K down to 39.2 bar at
    // 3566 K, which means dS went negative. That defect is recorded in
    // docs/DATA.md; pretending the sweep passes would hide it.
    let db = Db::standard();
    for name in ["NaCl", "Fe"] {
        for from in ["(l)", "(s)"] {
            let r = transition(&db, &format!("{name}{from}"), &format!("{name}(g)"));
            let mut prev = f64::NEG_INFINITY;
            let mut t = 400.0;
            while t <= 2800.0 {
                let p = p_sat(&db, &r, t);
                assert!(
                    p > prev,
                    "{name}{from}: p_sat fell at {t} K ({p} after {prev})"
                );
                assert!(p.is_finite() && p > 0.0, "{name}{from} at {t} K: {p}");
                prev = p;
                t *= 1.2;
            }
        }
    }
}

#[test]
fn the_bisection_and_the_pressure_agree_at_the_boiling_point() {
    // exp(-dG/RT) must be exactly 1 bar where dG crosses zero. Units slip
    // (kJ against J) or a wrong standard pressure land here, not in the
    // calibrated tolerance of the Tb test.
    let db = Db::standard();
    for (name, want) in [("NaCl", 1738.0), ("Fe", 3134.0)] {
        let vap = transition(&db, &format!("{name}(l)"), &format!("{name}(g)"));
        let tb = boiling_point(&db, &vap);
        let p = p_sat(&db, &vap, tb);
        assert!((p - 1.0).abs() < 0.005, "{name}: p_sat(Tb) = {p} bar");
        assert!((tb - want).abs() / want < 0.05, "{name}: Tb {tb} vs {want}");
    }
}

#[test]
fn clapeyron_slope_matches_the_latent_heat() {
    // d ln p / d(1/T) = -dH_vap / R. This ties the derivative of the numeric
    // pressure curve to the enthalpy path, which the value-only checks cannot:
    // H and G can both be wrong in a way that still crosses zero near Tb.
    let db = Db::standard();
    for name in ["NaCl", "Fe"] {
        let vap = transition(&db, &format!("{name}(l)"), &format!("{name}(g)"));
        for t in [1400.0, 1800.0, 2600.0] {
            let (t1, t2) = (t - 20.0, t + 20.0);
            let slope =
                (p_sat(&db, &vap, t2).ln() - p_sat(&db, &vap, t1).ln()) / (1.0 / t2 - 1.0 / t1);
            let want = -vap.dh_standard(&db, t) * 1000.0 / R;
            let rel = (slope - want).abs() / want.abs();
            assert!(rel < 0.001, "{name} at {t} K: {slope} vs {want}");
        }
    }
}

#[test]
fn sublimation_identity_holds_across_a_dense_sweep() {
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
        for k in 0..24 {
            let t = 600.0 + k as f64 * 200.0;
            let want = fus.dg_standard(&db, t) + vap.dg_standard(&db, t);
            assert!(
                (want - sub.dg_standard(&db, t)).abs() < 1e-9,
                "{name} at {t} K: {want} vs {}",
                sub.dg_standard(&db, t)
            );
        }
    }
}

#[test]
fn only_the_downhill_direction_is_available() {
    // Both directions are in the network because the cascade refuses to climb.
    // If a pair ever has both dG <= 0, or neither, the network is wrong.
    let db = Db::standard();
    for (name, tb) in [("NaCl", 1738.0), ("Fe", 3000.0)] {
        let (gas, gas2) = (format!("{name}(l)"), format!("{name}(g)"));
        let (down, up) = both(&db, &gas, &gas2);
        for t in [1000.0, tb - 200.0, tb + 200.0, 4000.0] {
            let (a, b) = (down.dg_standard(&db, t), up.dg_standard(&db, t));
            assert!(
                (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0),
                "{name} at {t} K: {a} / {b} must straddle zero"
            );
            assert!(
                (a + b).abs() < 1e-9,
                "{name} at {t} K: reverses are {a} and {b}"
            );
        }
    }
}

#[test]
fn the_same_cascade_twice_gives_bitwise_identical_amounts() {
    let db = Db::standard();
    let (liq, gas, solid) = (
        db.require("NaCl(l)"),
        db.require("NaCl(g)"),
        db.require("NaCl(s)"),
    );
    let net = rules::phase_transitions(&db, &[(liq, gas), (solid, gas), (solid, liq)]);
    let run = || {
        let mut st = State::new(&db, 1500.0, 0.05, 1000.0);
        st.add(liq, 1.0);
        let rep = equilibrium::relax(&db, &mut st, &net, 128, 1e-14);
        assert!(rep.converged, "{rep:?}");
        st
    };
    let (a, b) = (run(), run());
    for i in 0..a.n.len() {
        assert_eq!(
            a.n[i].to_bits(),
            b.n[i].to_bits(),
            "species {i} differs between identical runs"
        );
    }
}

#[test]
fn the_liquid_amount_does_not_move_the_vapour_pressure() {
    // A pure condensed phase has unit activity while present, so a bigger
    // melt pool cannot change the equilibrium pressure - only run out.
    let db = Db::standard();
    let (liq, gas, solid) = (
        db.require("NaCl(l)"),
        db.require("NaCl(g)"),
        db.require("NaCl(s)"),
    );
    let net = rules::phase_transitions(&db, &[(liq, gas), (solid, gas), (solid, liq)]);
    let mut small = State::new(&db, 1300.0, 0.05, 1000.0);
    small.add(liq, 0.5);
    let mut big = State::new(&db, 1300.0, 5.0, 1000.0);
    big.add(liq, 5.0);
    assert!(
        equilibrium::relax(&db, &mut small, &net, 128, 1e-14).converged,
        "small pool did not converge"
    );
    assert!(
        equilibrium::relax(&db, &mut big, &net, 128, 1e-14).converged,
        "big pool did not converge"
    );
    assert!(
        (small.n[gas] - big.n[gas]).abs() < 1e-12,
        "vapour amount depends on the pool: {} vs {}",
        small.n[gas],
        big.n[gas]
    );
    assert!(small.n[liq] > 0.0 && big.n[liq] > 0.0);
}

#[test]
fn nothing_moves_without_matter() {
    let db = Db::standard();
    let (liq, gas, solid) = (
        db.require("NaCl(l)"),
        db.require("NaCl(g)"),
        db.require("NaCl(s)"),
    );
    let net = rules::phase_transitions(&db, &[(liq, gas), (solid, gas), (solid, liq)]);
    let mut st = State::new(&db, 2000.0, 0.05, 1000.0);
    let rep = equilibrium::relax(&db, &mut st, &net, 128, 1e-14);
    assert!(rep.converged, "empty state: {rep:?}");
    for (i, n) in st.n.iter().enumerate() {
        assert_eq!(*n, 0.0, "species {i} appeared from nothing");
    }
}

#[test]
fn two_substances_share_a_network_without_borrowing_atoms() {
    // 1 mol NaCl plus 1 mol Fe at 1900 K: above both melting points, and above
    // the salt's boiling point but far below the metal's, so the salt must all
    // be vapour, the iron must all be liquid, and neither may gain atoms from
    // the other.
    let db = Db::standard();
    let (nacl_l, nacl_g, fe_s, fe_l, fe_g) = (
        db.require("NaCl(l)"),
        db.require("NaCl(g)"),
        db.require("Fe(s)"),
        db.require("Fe(l)"),
        db.require("Fe(g)"),
    );
    let na = db.require("Na(s)");
    let pairs = [
        (nacl_l, nacl_g),
        (db.require("NaCl(s)"), nacl_g),
        (db.require("NaCl(s)"), nacl_l),
        (fe_s, fe_l),
        (fe_s, fe_g),
        (fe_l, fe_g),
    ];
    let net = rules::phase_transitions(&db, &pairs);
    let mut st = State::new(&db, 1900.0, 0.05, 1000.0);
    st.add(nacl_l, 1.0);
    st.add(fe_s, 1.0);
    let mass_before = st.mass_g(&db);
    let rep = equilibrium::relax(&db, &mut st, &net, 128, 1e-14);
    assert!(rep.converged, "{rep:?}");
    assert!(
        st.n[nacl_g] > 0.99,
        "salt did not evaporate: {}",
        st.n[nacl_g]
    );
    assert!(st.n[fe_l] > 0.99, "iron did not melt: {}", st.n[fe_l]);
    assert!(st.n[fe_s] < 1e-6, "solid iron survived: {}", st.n[fe_s]);
    assert!(
        st.n[fe_g] < 0.05,
        "iron evaporated too much: {}",
        st.n[fe_g]
    );
    assert!(
        st.n[na] < 1e-12,
        "NaCl split into its elements: {}",
        st.n[na]
    );
    assert!((atoms(&db, &st, 11) - 1.0).abs() < 1e-12, "sodium count");
    assert!((atoms(&db, &st, 17) - 1.0).abs() < 1e-12, "chlorine count");
    assert!((atoms(&db, &st, 26) - 1.0).abs() < 1e-12, "iron count");
    assert!(
        (st.mass_g(&db) - mass_before).abs() < 1e-9,
        "mass: {mass_before} -> {}",
        st.mass_g(&db)
    );
}
