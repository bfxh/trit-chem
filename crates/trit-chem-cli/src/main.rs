//! Command-line probes for `trit-chem-core`.

use trit_chem_core::element;
use trit_chem_core::equilibrium;
use trit_chem_core::rules;
use trit_chem_core::species::Db;
use trit_chem_core::state::State;
use trit_chem_core::thermo::T0;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cmd = args.first().map(String::as_str).unwrap_or("demo");
    match cmd {
        "demo" => demo(),
        "check" => check(),
        "elements" => {
            for a in &args[1..] {
                match element::by_symbol(a) {
                    Some(e) => println!(
                        "{:>3} {:<3} {:<12} mass={:<10.5} EN={:?} IE1={:?} eV r={:?} pm shells={} ox={:?} zone={}",
                        e.z,
                        e.symbol,
                        e.name,
                        e.mass,
                        e.electronegativity,
                        e.first_ie_ev,
                        e.atomic_radius_pm,
                        e.outer_shell_electrons,
                        e.oxidation_states,
                        e.category
                    ),
                    None => println!("no element: {a}"),
                }
            }
        }
        "thermo" => {
            let db = Db::standard();
            for id in &args[1..] {
                match db.index(id) {
                    Some(i) => {
                        let s = db.get(i);
                        println!(
                            "{:<16} H(298)={:9.2} S(298)={:8.2} G(298)={:9.2}  H(373)={:9.2} G(373)={:9.2}",
                            s.id,
                            s.h(T0),
                            s.s(T0),
                            s.g(T0),
                            s.h(373.15),
                            s.g(373.15)
                        );
                    }
                    None => println!("no species: {id}"),
                }
            }
        }
        other => {
            eprintln!("unknown command {other:?}; use demo | check | elements | thermo");
            std::process::exit(2);
        }
    }
}

fn demo() {
    let db = Db::standard();
    let mut st = State::new(&db, T0, 1.0, 1.0);
    let na = db.require("Na(s)");
    let hp = db.require("H+(aq)");
    let cl = db.require("Cl-(aq)");
    let water = db.require("H2O(l)");
    st.add(na, 0.10);
    st.add(hp, 0.05);
    st.add(cl, 0.05);
    st.add(water, 55.5);

    let acid = rules::acid_metal(&db, 11, 1).expect("Na + acid rule");
    let wet = rules::metal_water(&db, 11, 1).expect("Na + water rule");
    let neut = rules::neutralization(&db).expect("neutralisation rule");
    let network = [acid, wet, neut];

    println!("== ranked by standard dG at 298.15 K (kJ per extent unit) ==");
    println!("(standard, not the live quotient: products start at zero, so the live");
    println!(" quotient is -inf until the first move -- it is reported per move instead)");
    let mut ranked: Vec<_> = network
        .iter()
        .map(|r| (r.dg_standard(&db, st.t), r))
        .collect();
    ranked.sort_by(|a, b| a.0.total_cmp(&b.0));
    for (dg, r) in &ranked {
        println!("{dg:9.3}  {}", r.label);
    }

    let before_mass = st.mass_g(&db);
    let before_charge = st.net_charge(&db);
    let rep = equilibrium::relax(&db, &mut st, &network, 64, 1e-12);
    println!(
        "\n== cascade: {} moves, {} sweeps, converged={} ==",
        rep.steps.len(),
        rep.sweeps,
        rep.converged
    );
    for s in &rep.steps {
        println!(
            "{:+10.5} mol  (dG0 {:9.3}, live at move {:9.3})  {}",
            s.extent, s.dg_standard, s.dg_before, s.label
        );
    }

    println!("\n== final composition ==");
    for (i, n) in st.n.iter().enumerate() {
        if *n > 1e-12 {
            println!("{:>10.6} mol  {}", n, db.get(i).id);
        }
    }
    let h = st.n[hp] / st.liquid_l;
    let oh = st.n[db.require("OH-(aq)")] / st.liquid_l;
    // The M0 network has no water autoionisation pair, so H+ and OH- are not
    // tied by Kw. Report pH from whichever is present and say which one.
    let ph = if h > 1e-12 {
        format!("{:.3} (from H+)", -h.log10())
    } else if oh > 1e-12 {
        format!("{:.3} (from OH-, assuming Kw=1e-14)", 14.0 + oh.log10())
    } else {
        "undefined (no H+ or OH-)".to_string()
    };
    let ph2 = st.n[db.require("H2(g)")] * 8.314_462_618_153_24 * st.t / (st.gas_l * 1e-3) / 1e5;
    println!("\npH = {ph}   p(H2) = {ph2:.4} bar");
    println!(
        "mass {:.6} -> {:.6} g (delta {:.2e}); charge {:.3e} -> {:.3e}",
        before_mass,
        st.mass_g(&db),
        st.mass_g(&db) - before_mass,
        before_charge,
        st.net_charge(&db)
    );
}

fn check() {
    let db = Db::standard();
    let mut worst = (0.0_f64, "");
    let mut n = 0;
    for (i, sp) in db.species().iter().enumerate() {
        let Some(tab) = sp.dgf298 else { continue };
        if sp.phase == trit_chem_core::Phase::Aqueous {
            continue;
        }
        let Some(implied) = db.implied_dgf298(i) else {
            println!("{:>16}: no element reference entropy", sp.id);
            continue;
        };
        let err = (implied - tab).abs();
        if err > worst.0 {
            worst = (err, sp.id);
        }
        n += 1;
    }
    println!(
        "checked {n} condensed/gas species; worst |dGcalc - dGtab| = {:.3} kJ/mol ({})",
        worst.0, worst.1
    );
}
