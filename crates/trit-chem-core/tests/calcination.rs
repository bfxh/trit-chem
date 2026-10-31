//! Lime-kiln check: calcination is endothermic and becomes feasible only when
//! T*dS overtakes dH. This exercises temperature-dependent Gibbs energies built
//! from Shomate fits and, for CaCO3, a constant Cp.

use trit_chem_core::reaction::{Reaction, Term};
use trit_chem_core::species::Db;
use trit_chem_core::thermo::T0;

fn calcination(db: &Db) -> Reaction {
    Reaction::new(
        db,
        "CaCO3(s) -> CaO(s) + CO2(g)",
        vec![
            Term::new(db.require("CaCO3(s,calcite)"), -1.0),
            Term::new(db.require("CaO(s)"), 1.0),
            Term::new(db.require("CO2(g)"), 1.0),
        ],
    )
    .unwrap()
}

#[test]
fn standard_values_at_25c() {
    let db = Db::standard();
    let r = calcination(&db);
    let dg = r.dg_standard(&db, T0);
    let dh = r.dh_standard(&db, T0);
    // NBS table: -604.03 + (-394.36) - (-1128.8) = +130.41 kJ.
    assert!((dg - 130.41).abs() < 1.0, "{dg}");
    // -635.09 + (-393.51) - (-1206.9) = +178.3 kJ.
    assert!((dh - 178.3).abs() < 1.0, "{dh}");
    assert!(dg > 0.0, "calcination must not be spontaneous at 25 C");
}

#[test]
fn dg_falls_monotonically_and_crosses_zero_near_1100k() {
    let db = Db::standard();
    let r = calcination(&db);
    let temps: Vec<f64> = (0..=20).map(|k| 298.15 + f64::from(k) * 60.0).collect();
    let vals: Vec<f64> = temps.iter().map(|t| r.dg_standard(&db, *t)).collect();
    for pair in vals.windows(2) {
        assert!(pair[1] < pair[0], "dG must fall with T: {pair:?}");
    }
    let mut cross = None;
    for (i, v) in vals.iter().enumerate() {
        if *v < 0.0 {
            cross = Some((temps[i - 1], temps[i]));
            break;
        }
    }
    let (lo, hi) = cross.expect("a zero crossing below 1500 K");
    assert!(lo > 950.0 && hi < 1200.0, "crossing in [{lo}, {hi}]");
    assert!(r.dg_standard(&db, 1500.0) < -50.0);
}

#[test]
fn enthalpy_dominates_below_the_crossing() {
    let db = Db::standard();
    let r = calcination(&db);
    let t = 800.0;
    assert!(r.dh_standard(&db, t) > 0.0);
    assert!(r.ds_standard(&db, t) > 0.0);
    assert!(r.dg_standard(&db, t) > 0.0);
}

#[test]
fn entropy_integration_matches_janaf_at_1000k() {
    let db = Db::standard();
    let co2 = db.get(db.require("CO2(g)"));
    // NIST-JANAF: S(1000 K) = 269.2 J/(mol K).
    assert!((co2.s(1000.0) - 269.2).abs() < 1.5, "{}", co2.s(1000.0));
    let h2o = db.get(db.require("H2O(g)"));
    // NIST-JANAF: S(1000 K) = 232.6 J/(mol K).
    assert!((h2o.s(1000.0) - 232.6).abs() < 1.5, "{}", h2o.s(1000.0));
}
