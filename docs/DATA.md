# Data provenance

Every generated Rust table has a generator; nothing was typed by hand into
`src/`.

| Data | Source | Fetched by | Licence status |
| --- | --- | --- | --- |
| 118 elements (mass, EN, IE1, radius, shells, phase, melt/boil, density, category) | mendeleev 1.3.0 SQLite database (MIT, Lukasz Mentel) | `scripts/fetch_elements.py` (pinned wheel, SHA-256 checked) -> `data/external/elements.json` | **MIT** — compatible with this repo; MIT notice kept in `data/external/MENDELEEV-LICENSE` |
| Common oxidation states | CRC Handbook 97th ed., table "Oxidation states of the elements" | curated dict in `gen_elements.py` | facts |
| Shomate coefficients (32 WebBook pages fetched; 25 of the 56 species carry fits, up to 5 ranges each) | NIST WebBook / NIST-JANAF (Chase 1998) | `scripts/fetch_shomate.py` | US Government compilation; attribution required |
| 25 C `dHf`, `dGf`, `S` for 56 species | NIST-JANAF and NBS/CODATA, transcribed in `data/external/thermo_table.txt` | curated list in `gen_thermo.py` | facts with source column per row |
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
python scripts/derive_liquid_anchors.py  # prints the NaCl(l)/Fe(l) dHf/S/dGf anchors
python scripts/gen_thermo.py       # writes species_data.rs
cargo test --workspace             # data consistency gates run here
```

`fetch_shomate.py` is the only step that needs the network. It keeps every
fitted range, not just the first, because NIST lists ranges as table columns,
and it reads NIST's `x10-n` cells as scientific notation: dropping the exponent
is not a rounding error (Fe's liquid B coefficient is -1.884667e-8, and the two
stray integers left behind also shift every column after them). `gen_thermo.py`
prepends a constant-Cp bridge for rows whose fits start above 298.15 K - see
`docs/DESIGN.md` §2.3 - and derives the molten-salt and molten-iron anchors
from the measured heat of fusion at the melting point.

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
7. **Cp range sanity**: a Shomate row whose first fit starts above 298.15 K must
   reproduce its own boundary heat capacity within 30% there, or be bridged;
   bridged rows may step at most 20% in Cp and 0.1 kJ/mol in H at the seam.
8. **Clausius-Clapeyron consistency**: `d ln p_sat / d(1/T)` must equal
   `-dH_vap/R` within 0.1% at 1400/1800/2600 K for both substances; this is what
   caught the moving-node Simpson integral.
9. **Boiling invariants**: `p_sat` strictly increasing up to 2800 K, 1 bar at the
   bisected boiling point, `dG_sub = dG_fus + dG_vap` over a 24-point sweep,
   both directions of a pair straddling zero, identical cascade results run twice
   (bitwise), vapour amount independent of pool size, empty state stays empty,
   and a shared NaCl+Fe network conserving each element to 1e-12.

Run `cargo run -p trit-chem-cli -- check` to print the worst recomputation error
without running the test suite.

## Known holes

* **Water vapour below 500 K.** NIST's fit starts at 500 K. The engine
  extrapolates; steam at 25 C is not used in the demo for this reason.
* **Above the highest fitted range.** Same extrapolation going up, and there it
  breaks monotonicity: `NaCl(s)` (single fit ending at 1073 K) gives a *falling*
  vapour pressure above about 3000 K - 77.1 bar at 2972 K down to 39.2 bar at
  3566 K, which is `dS < 0`. No phase transition is used up there, so nothing
  shipped depends on it, but do not ask the model for a 4000 K salt vapour until
  the plateau above `hi` is generated too.
* **Iron's solid entropy path.** NIST publishes five overlapping solid ranges
  whose seam values disagree, so integrating `Cp/T` from 298.15 K lands
  8.0 J/(mol K) below the absolute entropy of the fit that covers 1811 K (engine
  84.2 against 92.2). The molten-iron row inherits the offset through its
  fusion anchor, which is why iron boils 4.2% low (3003 K against 3134 K). The
  vapour row is accurate to 0.2 J/(mol K), and the NaCl solid and liquid rows
  agree with their own fits to within 0.15 J/(mol K) at the melting point.
* **Bridge heat capacities are effective, not measured.** Each bridge constant is
  chosen so `H` meets the published fit at the seam; the phase often does not
  exist at those temperatures (molten NaCl below 1074 K), so only differences
  above the seam are load-bearing.
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
