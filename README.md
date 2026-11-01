# trit-chem

Deterministic, dependency-free thermochemistry and reaction engine in Rust.
Built as the chemistry foundation (AE-ENGINE `F0`) for the simulation game: real
elements, real formation energies, reactions ranked by Gibbs energy instead of a
hand-written recipe table.

Status: **M0, local only, not published.** 28 tests pass (`cargo test --workspace`).

## What works today

```text
$ cargo run -p trit-chem-cli -- demo
cargo run -p trit-chem-cli -- kiln      # CaCO3 -> CaO + CO2 over temperature
cargo run -p trit-chem-cli -- phase     # ice / liquid / vapour equilibrium vs T
cargo run -p trit-chem-cli -- melt      # NaCl and Fe melting, phase at each T
cargo run -p trit-chem-cli -- boil      # NaCl and Fe boiling / sublimation vs T

== ranked by standard dG at 298.15 K (kJ per extent unit) ==
 -523.820  2Na(s) + 2H+(aq) -> 2Na+(aq) + H2(g)
 -364.040  2Na(s) + 2H2O(l) -> 2Na+(aq) + 2OH-(aq) + H2(g)
  -79.890  H+(aq) + OH-(aq) -> H2O(l)

== cascade: 2 moves, 2 sweeps, converged=true ==
  +0.02500 mol  2Na(s) + 2H+(aq) -> 2Na+(aq) + H2(g)
  +0.02500 mol  2Na(s) + 2H2O(l) -> 2Na+(aq) + 2OH-(aq) + H2(g)
```

That is the F0 acceptance scenario from the game design: hydrochloric acid is
consumed first, then the remaining sodium reacts with water. Mass and charge are
conserved to machine precision.

Delivered in M0:

| Piece | Content |
| --- | --- |
| Elements | all 118, with Pauling electronegativity, first ionization energy, radius, outer-shell electrons, oxidation states, phase/melt/boil/density |
| Thermodynamics | standard-state H, S, Cp, G; NIST Shomate fits where available, constant Cp otherwise |
| Species | 56 species (inorganic gases, metals, oxides, aqueous ions, molten and vapour salt/iron) with dHf/dGf/S at 25 C |
| Activities | Davies ionic-strength correction (switchable to ideal); ionic strength reported by the CLI |
| Reactions | element + charge balance checking, dG/dH/dS, K, reaction quotient |
| Kinetics | Arrhenius rate constants, mass action, fixed-step RK4 |
| Equilibrium | forward-only downhill cascade limited by the first exhausted reagent |
| Rules | acid + metal, metal + water, neutralization, combustion, precipitation, displacement, water/salt/metal phase transitions (melt, boil, sublime and the reverse) |
| Redox | 38 half-reactions with E0, Nernst equation, cell potential |
| CLI | `demo`, `check`, `elements`, `thermo`, `kiln`, `phase`, `melt`, `boil` |

## Quick start

```text
cargo test --workspace
cargo run -p trit-chem-cli -- demo
cargo run -p trit-chem-cli -- kiln      # CaCO3 -> CaO + CO2 over temperature
cargo run -p trit-chem-cli -- phase     # ice / liquid / vapour equilibrium vs T
cargo run -p trit-chem-cli -- melt      # NaCl and Fe melting, phase at each T
cargo run -p trit-chem-cli -- boil      # NaCl and Fe boiling / sublimation vs T
cargo run -p trit-chem-cli -- thermo "H2O(l)" "NaCl(s)"
cargo run -p trit-chem-cli -- check     # data self-consistency
```

## Honest limits

- **Activity coefficients.** Ions use the Davies equation (a few percent at
  I <= 0.1 M, approximate to ~0.5 M); neutral solutes and gases are ideal.
  Above 0.5 M, Pitzer or SIT is needed and is not implemented.
- **Constant Cp for ions.** Aqueous ions assume Cp = 0, so their temperature
  dependence is flat. Fine near 25 C, wrong near boiling.
- **Phase transitions are partial.** Water (ice Ih / liquid / vapour) and the
  fusion, boiling and sublimation of NaCl and Fe are modelled as separate
  species with both directions in the network; the transition temperature is
  where the model's own condensed-to-gas dG crosses zero rather than a table
  entry. Boiling: NaCl 1738 -> 1788 K (+2.9%), Fe 3134 -> 3003 K (-4.2%).
  Latent heats of vaporisation are within 2.6% (salt) and 2.5% (iron); the
  sublimation enthalpies at 25 C are within 0.2%. No nucleation barrier (no
  supercooled water), no volume change, no polymorphs, and the iron gap is an
  entropy-path error explained in `docs/DESIGN.md` §4.
- **The cascade is not a Gibbs minimiser.** It only moves downhill, one reaction
  at a time, and is path dependent. A White-Johnson-Dantzig element-potential
  solve is the M2 target; see `docs/PAPERS.md`.
- **RK4 is not stiffness-safe.** Combustion mechanisms will need an implicit
  integrator.
- **Element data** comes from mendeleev 1.3.0 (MIT). The build script
  downloads a pinned wheel, verifies its SHA-256 and commits only derived
  facts. The CC BY-SA source used in the first draft has been removed; see
  `docs/DATA.md`.

## Layout

```
crates/trit-chem-core/    engine (std only, no dependencies)
crates/trit-chem-cli/     probe binary
docs/DESIGN.md            model choices, error bands, audit of the v7.1 draft
docs/PAPERS.md            annotated bibliography
docs/DATA.md              data provenance and regeneration
scripts/                  data generators (run from the repo root)
```

## Regenerating data

```text
python scripts/fetch_elements.py   # pinned mendeleev wheel -> elements.json
python scripts/gen_elements.py     # elements.json -> elements_data.rs
python scripts/fetch_shomate.py    # NIST WebBook -> shomate.json
python scripts/gen_thermo.py       # curated table + shomate.json -> species_data.rs
cargo test --workspace
```

Licence: MIT. See `LICENSE`.
