#!/usr/bin/env python3
"""Rebuild data/external/elements.json from mendeleev (MIT).

Why this source: the previous table (Bowserinator/Periodic-Table-JSON) is
CC BY-SA 3.0, whose ShareAlike clause would propagate into the generated Rust
table. mendeleev is MIT, so the derived facts can live in an MIT repository
with attribution. The package ships a SQLite database; this script downloads a
pinned wheel, verifies its SHA-256, and reads the database without installing
anything.

Outputs:
  data/external/elements.json          derived facts (committed)
  data/external/MENDELEEV-LICENSE      MIT notice (committed)
"""
import hashlib
import json
import pathlib
import re
import sqlite3
import urllib.request
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "external" / ".cache"
OUT = ROOT / "data" / "external" / "elements.json"
LIC = ROOT / "data" / "external" / "MENDELEEV-LICENSE"

VERSION = "1.3.0"
WHEEL = f"mendeleev-{VERSION}-py3-none-any.whl"
URL = ("https://files.pythonhosted.org/packages/77/c0/"
       "c6fe26d8eda3b6eb1c4e3265ced9879dd95f0ca51d3569b042781f98ea76/" + WHEEL)
SHA256 = "76f415767b654634e51dead755513ef431abf200e892eab8468aa659ff218d1c"

SERIES = {
    "Nonmetals": "nonmetal",
    "Noble gases": "noble gas",
    "Alkali metals": "alkali metal",
    "Alkaline earth metals": "alkaline earth metal",
    "Metalloids": "metalloid",
    "Halogens": "halogen",
    "Poor metals": "post-transition metal",
    "Transition metals": "transition metal",
    "Lanthanides": "lanthanide metal",
    "Actinides": "actinide metal",
}


def wheel_path() -> pathlib.Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / WHEEL
    if not p.exists():
        print(f"downloading {URL}")
        urllib.request.urlretrieve(URL, p)
    got = hashlib.sha256(p.read_bytes()).hexdigest()
    if got != SHA256:
        raise SystemExit(f"wheel hash mismatch: {got} != {SHA256}")
    return p


def outer_shell(config: str) -> int:
    total, max_n = 0, 0
    for n, _l, k in re.findall(r"(\d)([spdf])(\d*)", config):
        n = int(n)
        k = int(k) if k else 1
        if n > max_n:
            max_n, total = n, k
        elif n == max_n:
            total += k
    return total


def phase_at_298(melt, boil) -> str:
    if melt is None and boil is None:
        return "unknown"
    if melt is None:
        return "solid" if boil is not None and boil > 298.15 else "unknown"
    if melt > 298.15:
        return "solid"
    if boil is None or boil > 298.15:
        return "liquid"
    return "gas"


def main():
    wp = wheel_path()
    with zipfile.ZipFile(wp) as z:
        db_bytes = z.read("mendeleev/elements.db")
        LIC.write_bytes(z.read(f"mendeleev-{VERSION}.dist-info/licenses/LICENSE"))
    db = CACHE / "elements.db"
    db.write_bytes(db_bytes)
    c = sqlite3.connect(db)
    series = dict(c.execute("select id, name from series"))

    rows = c.execute(
        "select atomic_number, symbol, name, atomic_weight, en_pauling, atomic_radius,"
        " density, electronic_configuration, series_id from elements order by atomic_number"
    ).fetchall()

    out = []
    for z, sym, name, mass, en, radius, density, config, sid in rows:
        ie = c.execute(
            "select ionization_energy from ionizationenergies"
            " where atomic_number=? and ion_charge=0",
            (z,),
        ).fetchone()
        pts = c.execute(
            "select melting_point, boiling_point from phasetransitions"
            " where atomic_number=? order by id",
            (z,),
        ).fetchall()
        melt = next((p[0] for p in pts if p[0] is not None), None)
        boil = next((p[1] for p in pts if p[1] is not None), None)
        ox = sorted(
            r[0]
            for r in c.execute(
                "select oxidation_state from oxidationstates"
                " where atomic_number=? and category='main'",
                (z,),
            )
        )
        out.append(
            {
                "z": z,
                "symbol": sym,
                "name": name,
                "mass": mass,
                "electronegativity": en,
                "first_ie_ev": None if ie is None else ie[0],
                "atomic_radius_pm": radius,
                "outer_shell_electrons": outer_shell(config),
                "oxidation_states_mendeleev": ox,
                "phase_298": phase_at_298(melt, boil),
                "melt_k": melt,
                "boil_k": boil,
                "density_g_cm3": density,
                "category": SERIES.get(series.get(sid), "unknown"),
                "electron_configuration": config,
            }
        )

    OUT.write_text(
        json.dumps(
            {
                "provenance": {
                    "package": "mendeleev",
                    "version": VERSION,
                    "license": "MIT",
                    "wheel": WHEEL,
                    "wheel_sha256": SHA256,
                    "url": URL,
                },
                "elements": out,
            },
            indent=1,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    missing = {
        k: sum(1 for e in out if e[k] is None)
        for k in ("electronegativity", "first_ie_ev", "atomic_radius_pm", "melt_k", "boil_k")
    }
    print(f"wrote {OUT} ({len(out)} elements); null counts: {missing}")


if __name__ == "__main__":
    main()
