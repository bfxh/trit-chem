# Data provenance

Every generated Rust table has a generator; nothing was typed by hand into
`src/`.

| Data | Source | Fetched by | Licence status |
| --- | --- | --- | --- |
| 118 elements (mass, EN, IE1, radius, shells, phase, melt/boil, density, category) | mendeleev 1.3.0 SQLite database (MIT, Lukasz Mentel) | `scripts/fetch_elements.py` (pinned wheel, SHA-256 checked) -> `data/external/elements.json` | **MIT** — compatible with this repo; MIT notice kept in `data/external/MENDELEEV-LICENSE` |
| Common oxidation states | CRC Handbook 97th ed., table "Oxidation states of the elements" | curated dict in `gen_elements.py` | facts |
| Shomate coefficients (20 species, up to 3 ranges each) | NIST WebBook / NIST-JANAF (Chase 1998) | `scripts/fetch_shomate.py` | US Government compilation; attribution required |
| 25 C `dHf`, `dGf`, `S` for 50 species | NIST-JANAF and NBS/CODATA, transcribed in `data/external/thermo_table.txt` | curated list in `gen_thermo.py` | facts with source column per row |
| Standard reduction potentials (38 couples) | same 25 C compilation, Table 3 | curated list in `redox.rs` | facts |

Each `Species` row carries a `source` field pointing at `NIST-JANAF 4th ed.
(Chase 1998)` or `NBS/CODATA 25 C table`. Keep it that way; an unsourced row is
a bug.

## Licensing

**Element data: cleared.** The first source (Bowserinator/Periodic-Table-JSON)
was CC BY-SA 3.0; ShareAlike would have propagated into the generated table and
into downstream game code. It has been replaced with mendeleev 1.3.0, which is
MIT (Copyright (c) 2015 Lukasz Mentel). The wheel is downloaded by
`fetch_elements.py`, its SHA-256 is verified, and only derived facts are
committed as `data/external/elements.json`. The MIT notice is kept verbatim in
`data/external/MENDELEEV-LICENSE`. No ShareAlike, no per-file copyleft.

Oxidation states are a small curated overlay (common states, CRC Handbook)
because mendeleev's "main" set is sparse for some metals. `gen_elements.py`
asserts that mendeleev's main set is a subset of the overlay; a typo in either
fails generation. These are facts, not a copied compilation.

The NIST WebBook pages state: "Data compilation copyright by the U.S. Secretary
of Commerce on behalf of the U.S.A." NIST data is distributed for public use
with attribution; do not strip the source strings.

## Regeneration

From the repository root:

```text
python scripts/fetch_elements.py   # network; pinned mendeleev wheel -> elements.json
python scripts/gen_elements.py     # elements.json -> elements_data.rs
python scripts/fetch_shomate.py    # network; writes data/external/shomate.json
python scripts/gen_thermo.py       # writes species_data.rs
cargo test --workspace             # data consistency gates run here
```

`fetch_shomate.py` is the only step that needs the network. It keeps every
fitted range, not just the first, because NIST lists ranges as table columns.

## Consistency gates

These run in `cargo test`:

1. **Element table**: 118 rows, strictly ordered by Z, spot values (Na mass and
   EN, Cl first ionization energy).
2. **Every formula parses** and yields at least one element.
3. **`dGf` recomputation**: for every non-ionic species, `dGf` recomputed from
   `dHf`, `S` and the element reference entropies must match the tabulated value
   within 1.0 kJ/mol. Current worst: 0.65 kJ/mol (AgCl(s), a rounding
   difference between the two compilations).
4. **Shomate sanity**: the H2 fit must reproduce S(298.15) = 130.68 J/mol/K and
   H(298.15) = 0.
5. **Balance checker** rejects deliberately bad oxygen and charge balances.
6. **Acceptance**: the sodium/acid/water scenario reproduces -523.82 and
   -364.04 kJ per extent, consumes acid first, and conserves mass and charge.

Run `cargo run -p trit-chem-cli -- check` to print the worst recomputation error
without running the test suite.

## Known holes

* **Water vapour below 500 K.** NIST's fit starts at 500 K. The engine
  extrapolates; steam at 25 C is not used in the demo for this reason.
* **Aqueous ion Cp** is zero by fiat. Ion `Cp` values are small and often
  negative; leaving them at zero is safer than inventing them.
* **No F2, Br2, I2 reference species**, so those elements have no reference
  entropy and are skipped by the `dGf` gate. Add them before using the
  corresponding ions.
* **`AgCl(s)` and other solids** carry only 25 C data; `shomate.json` has no
  entry, so they use constant Cp. `SiO2(s)` and `CaCO3(s)` fall back to manually
  curated constant Cp values (`CP_CONST` in `gen_thermo.py`); their provenance is
  weaker than the Shomate rows and should be upgraded.
* **Reduction potentials** cover 38 couples, not the full series. They are data
  for rules and for `dG = -nFE`, not yet wired into the reaction generator.
