#!/usr/bin/env python3
"""Generate crates/trit-chem-core/src/species_data.rs.

Merges:
  * curated 25 C thermochemistry (dHf, dGf, S) from the verified table in
    data/external/thermo_table.txt (sources: NIST-JANAF / CODATA compilations);
  * Shomate blocks scraped from NIST WebBook (data/external/shomate.json).

Aqueous ions carry Cp = 0 (standard partial molar heat capacities are small and
the M0 model is documented as near-298 K only).  Every row keeps `dgf` so a
test can recompute dG = dH - T*dS and compare.
"""
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
SHOM = ROOT / "data" / "external" / "shomate.json"
OUT = ROOT / "crates" / "trit-chem-core" / "src" / "species_data.rs"

SRC_JANAF = "NIST-JANAF 4th ed. (Chase 1998) via NIST WebBook"
SRC_NBS = "NBS/CODATA 25 C table (data/external/thermo_table.txt)"

# id, formula, phase, charge, dHf298 kJ/mol, dGf298 kJ/mol, S298 J/mol/K, source
SPECIES = [
    ("H2(g)", "H2", "gas", 0, 0.0, 0.0, 130.68, SRC_JANAF),
    ("O2(g)", "O2", "gas", 0, 0.0, 0.0, 205.14, SRC_JANAF),
    ("N2(g)", "N2", "gas", 0, 0.0, 0.0, 191.61, SRC_JANAF),
    ("Cl2(g)", "Cl2", "gas", 0, 0.0, 0.0, 223.07, SRC_JANAF),
    ("H2O(g)", "H2O", "gas", 0, -241.82, -228.57, 188.83, SRC_JANAF),
    ("CO(g)", "CO", "gas", 0, -110.53, -137.17, 197.67, SRC_JANAF),
    ("CO2(g)", "CO2", "gas", 0, -393.51, -394.36, 213.74, SRC_JANAF),
    ("CH4(g)", "CH4", "gas", 0, -74.81, -50.72, 186.26, SRC_JANAF),
    ("NH3(g)", "NH3", "gas", 0, -46.11, -16.45, 192.45, SRC_JANAF),
    ("SO2(g)", "SO2", "gas", 0, -296.83, -300.19, 248.22, SRC_JANAF),
    ("HCl(g)", "HCl", "gas", 0, -92.31, -95.30, 186.91, SRC_JANAF),
    ("H2O(l)", "H2O", "liquid", 0, -285.83, -237.13, 69.91, SRC_JANAF),
    ("H2SO4(l)", "H2SO4", "liquid", 0, -813.99, -690.00, 156.90, SRC_NBS),
    ("Na(s)", "Na", "solid", 0, 0.0, 0.0, 51.21, SRC_JANAF),
    ("NaCl(s)", "NaCl", "solid", 0, -411.15, -384.12, 72.13, SRC_JANAF),
    ("NaOH(s)", "NaOH", "solid", 0, -425.61, -379.49, 64.46, SRC_JANAF),
    ("Fe(s)", "Fe", "solid", 0, 0.0, 0.0, 27.28, SRC_JANAF),
    ("Cu(s)", "Cu", "solid", 0, 0.0, 0.0, 33.15, SRC_JANAF),
    ("Zn(s)", "Zn", "solid", 0, 0.0, 0.0, 41.63, SRC_JANAF),
    ("Al2O3(s)", "Al2O3", "solid", 0, -1675.7, -1582.3, 50.92, SRC_JANAF),
    ("S(s,rhombic)", "S", "solid", 0, 0.0, 0.0, 31.80, SRC_JANAF),
    ("C(s,graphite)", "C", "solid", 0, 0.0, 0.0, 5.740, SRC_JANAF),
    ("CaCO3(s,calcite)", "CaCO3", "solid", 0, -1206.9, -1128.8, 92.9, SRC_NBS),
    ("SiO2(s,quartz)", "SiO2", "solid", 0, -910.94, -856.64, 41.84, SRC_NBS),
    ("Al(s)", "Al", "solid", 0, 0.0, 0.0, 28.33, SRC_NBS),
    ("Ca(s)", "Ca", "solid", 0, 0.0, 0.0, 41.42, SRC_NBS),
    ("Si(s)", "Si", "solid", 0, 0.0, 0.0, 18.83, SRC_NBS),
    ("Mg(s)", "Mg", "solid", 0, 0.0, 0.0, 32.68, SRC_NBS),
    ("Pb(s)", "Pb", "solid", 0, 0.0, 0.0, 64.81, SRC_NBS),
    ("Ag(s)", "Ag", "solid", 0, 0.0, 0.0, 42.55, SRC_NBS),
    ("AgCl(s)", "AgCl", "solid", 0, -127.7, -109.79, 96.2, SRC_NBS),
    ("K(s)", "K", "solid", 0, 0.0, 0.0, 64.18, SRC_NBS),
    ("P(s,white)", "P", "solid", 0, 0.0, 0.0, 41.09, SRC_NBS),
    # aqueous ions (Cp = 0; see module doc)
    ("H+(aq)", "H", "aqueous", 1, 0.0, 0.0, 0.0, SRC_NBS),
    ("OH-(aq)", "OH", "aqueous", -1, -229.99, -157.24, -10.75, SRC_NBS),
    ("Na+(aq)", "Na", "aqueous", 1, -240.12, -261.91, 59.0, SRC_NBS),
    ("Cl-(aq)", "Cl", "aqueous", -1, -167.16, -131.23, 56.5, SRC_NBS),
    ("K+(aq)", "K", "aqueous", 1, -252.38, -283.27, 102.5, SRC_NBS),
    ("Ca2+(aq)", "Ca", "aqueous", 2, -542.83, -553.58, -53.1, SRC_NBS),
    ("Mg2+(aq)", "Mg", "aqueous", 2, -466.85, -454.8, -138.1, SRC_NBS),
    ("Fe2+(aq)", "Fe", "aqueous", 2, -89.1, -78.90, -137.7, SRC_NBS),
    ("Fe3+(aq)", "Fe", "aqueous", 3, -48.5, -4.7, -315.9, SRC_NBS),
    ("Cu2+(aq)", "Cu", "aqueous", 2, 64.77, 65.49, -99.6, SRC_NBS),
    ("Cu+(aq)", "Cu", "aqueous", 1, 71.67, 49.98, 40.6, SRC_NBS),
    ("Zn2+(aq)", "Zn", "aqueous", 2, -153.89, -147.06, -112.1, SRC_NBS),
    ("Ag+(aq)", "Ag", "aqueous", 1, 105.58, 77.11, 72.68, SRC_NBS),
    ("Pb2+(aq)", "Pb", "aqueous", 2, -1.7, -24.43, 10.5, SRC_NBS),
    ("NH4+(aq)", "NH4", "aqueous", 1, -132.51, -79.31, 113.4, SRC_NBS),
    ("NO3-(aq)", "NO3", "aqueous", -1, -205.0, -108.74, 146.4, SRC_NBS),
    ("SO4 2-(aq)", "SO4", "aqueous", -2, -909.27, -744.53, 20.1, SRC_NBS),
]

# Fallback constant Cp (J/mol/K) where NIST has no Shomate block for M0.
CP_CONST = {
    "C(s,graphite)": 8.527,
    "CaCO3(s,calcite)": 81.88,
    "SiO2(s,quartz)": 44.59,
    "H2SO4(l)": 138.9,
}


def rust_ident(s):
    out = []
    for ch in s:
        if ch.isalnum():
            out.append(ch)
        else:
            out.append("_")
    return "".join(out).upper()


def dedupe(blocks):
    seen = set()
    keep = []
    for b in blocks:
        key = (b["label"], round(b["lo"], 3), round(b["hi"], 3), round(b["a"], 6))
        if key in seen:
            continue
        seen.add(key)
        keep.append(b)
    keep.sort(key=lambda b: b["lo"])
    return keep


def main():
    shom = json.loads(SHOM.read_text(encoding="utf-8")) if SHOM.exists() else {}
    lines = [
        "// @generated by scripts/gen_thermo.py -- do not edit by hand.",
        "use crate::{Cp, Phase, ShomateBlock, Species};",
        "",
    ]
    statics = []
    rows = []
    for sid, formula, phase, charge, hf, gf, s, src in SPECIES:
        blocks = dedupe(shom.get(sid, []))
        if blocks:
            name = "SHOMATE_" + rust_ident(sid)
            statics.append(
                f"static {name}: [ShomateBlock; {len(blocks)}] = [\n"
                + "".join(
                    "    ShomateBlock { lo: %r, hi: %r, a: [%s] },\n"
                    % (float(b["lo"]), float(b["hi"]),
                       ", ".join(repr(float(b[k])) for k in "abcdefgh"))
                    for b in blocks
                )
                + "];\n"
            )
            cp = f"Cp::Shomate(&{name})"
        else:
            cp = f"Cp::Const({CP_CONST.get(sid, 0.0)})"
        rows.append(
            f'    Species {{ id: "{sid}", formula: "{formula}", phase: Phase::{phase.capitalize()}, '
            f"charge: {charge}, hf298: {hf}, s298: {s}, cp: {cp}, "
            f"dgf298: Some({gf}), source: {json.dumps(src, ensure_ascii=False)} }},"
        )
    lines.extend(statics)
    lines.append("/// Built-in species table, ordered as listed in scripts/gen_thermo.py.")
    lines.append(f"pub static SPECIES: [Species; {len(rows)}] = [")
    lines.extend(rows)
    lines.append("];")
    lines.append("")
    OUT.write_text("\n".join(lines), encoding="utf-8", newline="\n")
    print(f"wrote {OUT} ({len(rows)} species, {len(statics)} shomate tables)")


if __name__ == "__main__":
    main()
