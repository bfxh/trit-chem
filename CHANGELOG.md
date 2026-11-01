# Changelog

## 0.5.0 (2026-11-01) — boiling and sublimation of the salt and the metal

### Fixed

- **NIST scientific notation is parsed.** `fetch_shomate.py` carried an exponent
  helper that no call site used, so cells rendered as `-1.884667×10-8` were read
  as `-1.884667` plus two stray integers: the value truncated by up to eight
  orders of magnitude *and* every later column shifted. `shomate.json` held
  truncated liquid fits for Al2O3, CaO, Fe, SiO2 and Zn; nothing shipped was
  wrong only because phase filtering discarded those blocks until boiling needed
  a liquid fit. The parser now requires exactly one value per fitted range and
  drops the block otherwise. Falsified: re-fetching with the old parser fails
  `iron_melts_near_1811_k_with_the_measured_latent_heat` and
  `cascade_keeps_the_phase_that_is_stable`.
- **Heat capacity is no longer extrapolated as a polynomial below a fit's
  range.** `Cp(Fe gas, 298.15 K)` came out 79 J/(mol K) against the real 25.15,
  `Cp(NaCl liquid, 298.15 K)` 428 against 70.4 at the fit's own lower edge.
  `gen_thermo.py` prepends a constant-Cp bridge (298.15 K to the first fit)
  whenever a row's backward extrapolation misses that boundary value by more
  than 25%; the bridge constant is what makes `H` meet the published fit exactly
  at the seam. Falsified: with the bridge disabled the new data-shape gate goes
  red.

### Added

- `NaCl(g)` and `Fe(g)` species rows. `NaCl(l)` and `Fe(l)` now take their Cp
  shape from NIST's own liquid fits (1074-2500 K and 1809-3133 K) instead of a
  constant taken at the melting point, with the fusion anchors re-derived.
- `rules::phase_transitions` replaces `fusion_transitions` and
  `water_transitions`: one helper emitting both directions of any `(a, b)` phase
  pair, with `rules::water_pairs` for the three water pairs.
- `trit-chem-cli boil`: model boiling point, dG of vaporisation and sublimation,
  and the vapour pressure they imply.
- `scripts/derive_liquid_anchors.py` recomputes the `NaCl(l)` / `Fe(l)`
  `dHf`/`S`/`dGf` anchors from the fetched fits plus the measured latent heat, so
  the two hand-embedded rows stay recomputable instead of folklore.
- `tests/vapor.rs` (6 tests) plus `thermo::tests::no_row_extrapolates_a_runaway_
  heat_capacity_down_to_298_k` as a data-shape gate.
- Rejected experiment, recorded: splitting the entropy integral at every range
  edge. Against a 4000-slice reference (molten NaCl, 1738 K) the single 48-point
  pass is out by +0.05 J/(mol K) and the split by -0.11, because each segment's
  endpoint still samples the block claiming the boundary. The simpler code stays.

### Verified

| quantity | model | reference |
| --- | ---: | ---: |
| NaCl Tb (dG_vap = 0) | 1791 K (+3.1%) | 1738 K |
| NaCl dH_vap(1738 K) | 166.5 kJ/mol | 171 |
| NaCl dH_sub(298.15 K) | 229.7 kJ/mol | 229.7 |
| Fe Tb (dG_vap = 0) | 3001 K (-4.3%) | 3134 K |
| Fe dH_vap(3134 K) | 358.7 kJ/mol | 349.9 |
| Fe dH_sub(298.15 K) | 415.5 kJ/mol | 416.3 |
| dG_sub - (dG_fus + dG_vap) | < 1e-9 kJ | identity |
| 1 mol NaCl(l) in 1000 L | all vapour at 1900 K; stops at p_sat within 2% at 1200 K | - |

NaCl and Fe melting readings are unchanged (dG_fus still brackets 1074 K and
1811.15 K, dH_fus 28.0 and 13.8 kJ/mol) and the water validation table in
`docs/DESIGN.md` §4 still passes.

### Known gap

Iron's solid entropy path lands 8.0 J/(mol K) below the absolute entropy of the
fit that covers it (engine 84.2 against 92.2 at 1811 K), because NIST publishes
five overlapping solid ranges; the liquid row inherits the offset through its
fusion anchor, which is most of the -4.3% on iron's boiling point. See
`docs/DATA.md` and DESIGN §4.

## 0.4.0 (2026-10-31) — salt and metal melting; Shomate phase filter

### Fixed

- **Shomate fits are filtered by phase.** A NIST page for a solid also lists its
  liquid and gas fits; the generator kept them all, so above the melting point
  the "solid" species used liquid Cp and the phase change silently vanished.
  `gen_thermo.py` now keeps only the fits matching the species phase, and
  `tests/fusion.rs` pins the solid Cp just below melting.

### Added

- `NaCl(l)` and `Fe(l)`, derived from the measured heats of fusion
  (28.0 kJ/mol at 1074 K; 13.8 kJ/mol at 1811 K). Published 298 K liquid values
  are extrapolations that disagree with those enthalpies, so they are replaced
  by values anchored at the melting point.
- `rules::fusion_transitions` and `trit-chem-cli melt`.

### Verified

- NaCl: dG_fus > 0 at 1073 K, < 0 at 1075 K, dH_fus = 28.0 kJ/mol.
- Fe: dG_fus > 0 at 1810 K, < 0 at 1812 K, dH_fus = 13.8 kJ/mol.
- Cascade keeps the stable phase 50 K on either side of both melting points.
- Mutation: changing the NaCl liquid Cp from 70.4 to 100.0 fails two tests.

## 0.3.0 (2026-10-31) — water phase transitions

### Added

- `H2O(s)` (ice Ih) and the six water transitions as reactions
  (`rules::water_transitions`): melting/freezing, evaporation/condensation,
  sublimation/deposition. Both directions are present because the cascade only
  moves downhill.
- `trit-chem-cli phase`: equilibrium composition of 1 mol water in a 10 L
  headspace from 250 to 450 K.
- `tests/phase.rs`: ice stable below 0 C, liquid above, melting bracketed within
  0.3 K, saturation pressure at five temperatures, triple-point pressure,
  transition enthalpies.

### Verified

- p_sat: 0.0317 bar at 25 C (ref 0.03169), 0.9988 at 100 C (ref 1.01325),
  triple point 0.0061 bar (ref 0.00611).
- dH_fus = 6.01 kJ/mol (ref 6.008), dH_vap ~ 40.7 (ref 40.66).

### Notes

- Ice is not stable at 298 K, so its dGf(298.15) is anchored to the measured
  melting point rather than extrapolated. The dGf recomputation check still
  passes.
- Limits: ice Ih only, no volume change, no supercooling.

## 0.2.1 (2026-10-31) — Shomate boundary fix, lime-kiln demo

### Fixed

- `H(T)-H(298.15)` no longer subtracts a same-block value at 298.15 K. Each NIST
  block is already referenced to 298.15; the old code extrapolated the block
  covering T back to 298.15 and produced a step at every range boundary
  (CO2 at 1200 K: 3.1 kJ/mol).
- `S(T)-S(298.15)` is now integrated as `integral Cp/T` (48-point Simpson)
  instead of subtracting block entropies, for the same reason. Matches
  NIST-JANAF S(1000 K) for CO2 and H2O(g) within 0.1 J/(mol K).

### Added

- `CaO(s)` (NBS 25 C values + NIST Shomate fit).
- `trit-chem-cli kiln`: CaCO3 -> CaO + CO2 over 298-1600 K. dG(298) = +130.4
  kJ/mol, zero crossing at 1073 K, endothermic throughout.
- `tests/calcination.rs`: standard values, monotone dG, crossing bracket,
  entropy check against JANAF.

## 0.2.0 (2026-10-31) — ionic-strength activity coefficients

### Added

- `ActivityModel::{Ideal, Davies}`; `State` defaults to Davies.
- `State::ionic_strength()` and `State::activity_coefficient()`; the CLI prints
  both the ionic strength and the coefficients.
- `tests/activity.rs`: Davies reproduces measured NaCl mean activity
  coefficients (0.965 / 0.902 / 0.778 at 0.001 / 0.01 / 0.1 M); the correction
  cancels exactly in the acid reaction and shifts the water reaction by
  `4 RT ln(gamma)`; both models reach the same final state in the F0 scenario.

### Changed

- `solve_extent` no longer assumes Q is monotone: it scans 32 points for the
  first sign change and bisects inside that bracket (activity coefficients make
  Q non-monotone for high-charge ions).

### Known limits

- Davies is good to I <= 0.1 M, indicative to 0.5 M, unusable above; Pitzer/SIT
  is the next step.
- A is fixed at the 25 C value (0.509); at 100 C this costs ~20% in
  log10(gamma).
- Neutral solutes have no salting-out term.

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
