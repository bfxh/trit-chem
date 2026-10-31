# Changelog

## 0.1.1 (2026-10-31) — element source and licence

### Changed

- Element table rebuilt from **mendeleev 1.3.0** (MIT, (c) Lukasz Mentel) instead
  of Bowserinator/Periodic-Table-JSON (**CC BY-SA 3.0**). ShareAlike would have
  propagated into this MIT repository and into downstream game code.
- `scripts/fetch_elements.py` downloads a pinned wheel, verifies its SHA-256
  (`76f41576...218d1c`) and reads the SQLite database without installing the
  package. Only derived facts are committed (`data/external/elements.json`);
  the MIT notice is kept in `data/external/MENDELEEV-LICENSE`.
- Ionization energies are now taken directly in eV (no kJ/mol conversion).
- `phase_298` is derived from melting/boiling points; subliming elements with no
  melting point (carbon) are classified solid.
- Oxidation states remain a curated CRC overlay because mendeleev's main set is
  sparse for some metals (Cu: [2] only). Generation asserts that mendeleev's
  main set is a subset of the overlay; Se/Te gained +2 from that check.

### Added

- `tests/data_license.rs`: the MIT notice, the pinned hash and the absence of
  the CC BY-SA file are now test-enforced.

### Data deltas worth knowing

- Null coverage after the swap: electronegativity 33, first IE 10, radius 28,
  melting point 17, boiling point 30 (mostly superheavy elements).
- Densities for gases are g/cm3 (e.g. H = 8.2e-5), same unit as before.

## 0.1.0 (2026-10-31) — M0

First slice of the chemistry foundation for AE-ENGINE F0.

### Added

- `trit-chem-core`: 118-element table, formula parser, Shomate/constant-Cp
  thermochemistry, balance-checked reactions, dG/dH/dS/K, downhill reaction
  cascade, Arrhenius + RK4 kinetics, rule-based candidate reactions, 38
  reduction half-reactions with Nernst.
- `trit-chem-cli`: `demo` (F0 acceptance scenario), `check` (data
  self-consistency), `elements`, `thermo`.
- Data generators: `gen_elements.py`, `fetch_shomate.py`, `gen_thermo.py`,
  `module_budget.py`.

### Verified

- 2Na + 2H+ -> 2Na+ + H2: **-523.82 kJ per extent** (-261.91 per mol Na).
- 2Na + 2H2O -> 2Na+ + 2OH- + H2: **-364.04 kJ per extent** (-182.02 per mol Na).
- H+ + OH- -> H2O: -79.89 kJ.
- Cascade consumes HCl first, then water; mass and charge conserved.
- `dGf` recomputed from `dHf`, `S` and element references for 33 species:
  worst deviation 0.65 kJ/mol (AgCl(s), source rounding).

### Corrections to the design draft

- The draft's `-402 / -368 kJ` for the two sodium reactions do not reproduce;
  see `docs/DESIGN.md` §4.
- "Standard electrode potential" is a half-reaction property, not an element
  property (`redox::HALF_REACTIONS`).
- "Valence electrons" for d-block elements is not the outer-shell population;
  both are stored separately.

### Retracted during development

- An earlier draft of `docs/DESIGN.md` claimed that using `dG = dH - T*S`
  instead of anchoring at `dGf` would break electrolyte reactions by
  ~19 kJ/mol. Mutation testing showed the per-species offset cancels in
  charge-balanced reactions (reaction-level difference 0.03 kJ/mol). The anchor
  is still required for per-species correctness; the claim was narrowed.

### Known gaps

Listed in `README.md` and `docs/DESIGN.md` §6: no activity coefficients, zero
ion Cp, no phase transitions, forward-only cascade, non-stiff RK4, data licence
gate before publishing.
