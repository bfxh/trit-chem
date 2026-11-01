# Design and model boundaries (M0)

This document explains what M0 computes, how each number is obtained, what the
error bands are, and where the game-design draft (`AE-ENGINE v7.1`) does not
survive contact with the data.

## 1. Pipeline

```
formula string ──► Formula { (Z, count) } ──► balance check
                                                │
species table ──► H(T), S(T), G(T) ─────────────┴─► reaction dG(T), K(T)
                                                     │
state {mol, V, T} ──► activities ──► Q ──► live dG
                                                     │
                                    ┌────────────────┴──────────────┐
                                    │                               │
                            downhill cascade                  Arrhenius + RK4
                            (equilibrium-ish,                 (explicit kinetics,
                             deterministic, O(network))        stiff-hostile)
```

Two independent solvers are provided on purpose:

* **Cascade** (`equilibrium::relax`) is what the game uses per tick. It picks
  every reaction with negative live dG, drives it to `Q = K`, capped by the first
  reactant that runs out, and repeats. Each move lowers total Gibbs energy, so it
  terminates. It cannot run a reaction backwards; reverse directions must be
  supplied as their own reaction. This is a *coordinate descent on G*, not a
  global minimiser.
* **Mass-action RK4** (`kinetics::rk4`) is for cases where the rates are known.
  It is explicit and therefore unsuitable for stiff combustion chemistry.

## 2. Thermochemistry

### 2.1 Anchor: dGf, not dHf - T*S

For every species, `G(298.15)` is anchored at the tabulated standard Gibbs
energy of formation, then integrated with the Gibbs-Helmholtz relation:

```
G(T) = Gf(298) + [H(T) - H(298)] - T*[S(T) - S(298)] - (T - 298)*S(298)
```

Why the anchor matters: for aqueous ions the tabulated `dHf`, `dGf` and `S` are
**not** related by `dG = dH - T*dS`. The electron convention assigns the
solvated electron zero Gibbs energy, so `dHf - T*S` shifts Na+ by about
**19 kJ/mol**.

That offset is proportional to charge and **cancels in any charge-balanced
reaction**: `dHf - T*S` and the anchored formula give the same reaction dG (they
differ by 0.03 kJ/mol on the sodium test, which is data rounding). So this is not
a reaction-level trap; it is a per-species trap. `dHf - T*S` produces a G for
Na+ that is not the Na+ formation Gibbs energy, which breaks partial queries,
single-species lookups and any future unbalanced or constrained calculation.

**Retraction:** an earlier draft of this section claimed the wrong formula
"silently breaks every electrolyte reaction" and quoted a 19 kJ/mol reaction
error. A mutation test disproved that: the reaction energy is unchanged. The
anchor is still the correct implementation, but for the reason above.

For neutral species the two definitions agree to rounding; the test suite checks
that agreement for all 33 non-ionic species (worst deviation 0.65 kJ/mol,
AgCl(s), a source-rounding discrepancy).

### 2.2 Element reference entropies are per mole of atoms

`implied_dgf298()` recomputes `dGf` from `dHf` and `S` using
`dGf = dHf - T*(S_species - sum_k n_k S_element)`. The reference species for an
element is a molecule (`H2`, `O2`, `Cl2`), so its molar entropy must be divided
by the number of atoms it contains before multiplying by the atom count. Getting
this backwards shifts H2O(g) by ~11 kJ/mol and is exactly the kind of error that
looks plausible in spot checks; the self-consistency test catches it.

### 2.3 Temperature dependence

* NIST Shomate blocks (8 coefficients, up to three ranges per species) from
  NIST-JANAF via the WebBook, fetched by `scripts/fetch_shomate.py`.
* If a temperature falls outside every fitted range, the nearest block is used
  and the result is extrapolated. This is a documented approximation, not an
  error state. Water vapour, for example, has no NIST fit below 500 K.
* **Except where the polynomial is not extrapolable.** NIST fits some phases
  only from their melting or boiling point upward, and the Shomate `E/t^2` term
  blows up below the range: measured `Cp(Fe gas, 298.15 K) = 79 J/(mol K)`
  against the real 25.15, and `Cp(NaCl liquid, 298.15 K) = 428` against 70.4 at
  the fit's own lower edge. `gen_thermo.py` therefore prepends a constant-Cp
  bridge segment (298.15 K to the first fit) whenever a row's backward
  extrapolation misses its boundary value by more than 25%. The bridge's Cp is
  fixed so that `H` meets the NIST fit exactly at the seam, which keeps H, S and
  the entropy integral continuous; it is an *effective* heat capacity for a
  stretch of temperature the phase may not even occupy, not a measurement.
  Rows within tolerance (H2O gas from 500 K, NaCl gas from 2500 K, Fe liquid,
  whose fit is essentially constant) keep the published polynomial.
* Each NIST block is already referenced to H(298.15), so `H(T)-H(298.15)` is
  taken straight from the block covering T. Subtracting the same block's value
  at 298.15 would extrapolate it outside its range and put a step at every
  range boundary (measured 3.1 kJ/mol for CO2 at 1200 K before the fix).
* `S(T)-S(298.15)` is integrated numerically as `integral Cp/T` (48-point
  Simpson). Using the blocks' absolute entropies directly also steps at range
  boundaries, because only the first block is anchored at the tabulated
  entropy. The integration reproduces NIST-JANAF S(1000 K) to within 0.1 J/K
  for CO2 and H2O(g). Splitting the integral at every range edge was tried and
  rejected: against a 4000-slice reference (molten NaCl, 1738 K) the single pass
  is out by +0.05 J/(mol K) and the split by -0.11, because each segment's
  endpoint still samples the block that claims the boundary.
* Aqueous ions use `Cp = 0`. Standard partial molar heat capacities of ions are
  small and frequently negative; pretending they are zero is better than
  inventing values, but it means ion thermodynamics is a 25 C tool.

### 2.4 Standard states

Gas: ideal, 1 bar. Solute: 1 mol/L. Pure condensed phase: unit activity while
present. Gas partial pressure uses the ideal gas law with the headspace volume.

## 3. Activity model

`a = p/p0` for gases, `a = 1` for present pure condensed phases, and for solutes
`a = (c / 1 M) * gamma`.

`gamma` is 1 under `ActivityModel::Ideal`. Under `ActivityModel::Davies` (the
default) ions follow

```
log10(gamma) = -A z^2 ( sqrt(I)/(1+sqrt(I)) - 0.3 I ),   A = 0.509
I = 0.5 * sum c_i z_i^2
```

Checked against measured NaCl mean ionic activity coefficients at 25 C:
0.001 M -> 0.965 (measured 0.965), 0.01 M -> 0.902 (0.902), 0.1 M -> 0.782
(0.778). Accuracy degrades above 0.1 M; at 0.5 M Davies gives ~0.73 against a
measured 0.681, so treat 0.1-0.5 M as indicative and refuse to trust > 0.5 M.
Neutral solutes get no salting-out term. `A` is held at its 25 C value, which
costs about 20% in log10(gamma) at 100 C.

Two consequences worth knowing:

* The correction cancels exactly when the ionic charges on both sides match
  (2 H+ vs 2 Na+), so acid-metal reactions are unchanged; it does *not* cancel
  for the metal-water reaction (2 ions of each sign produced), which becomes
  about `4 RT ln(gamma)` more favourable.
* Because `gamma` depends on `I`, the reaction quotient is no longer guaranteed
  monotone in extent. The extent solver scans 32 points for the first sign
  change and bisects inside it instead of assuming monotonicity.

## 4. Phase transitions

Phase transitions are ordinary reactions between phase-specific species, so
they run through the same cascade:

| reaction | meaning |
| --- | --- |
| H2O(s) -> H2O(l) | melting |
| H2O(l) -> H2O(g) | evaporation |
| H2O(s) -> H2O(g) | sublimation |
| the three reverse directions | freezing, condensation, deposition |

Both directions must be in the network because the cascade only moves downhill.

Ice data is *derived*, not tabulated at 298.15 K, because ice is not stable
there. dHf is anchored to the measured heat of fusion (6.008 kJ/mol at
273.15 K) with a constant Cp of 38.0 J/(mol K); dGf(298.15) is then chosen so
that `G_ice = G_liquid` at 273.15 K exactly. That keeps the species row
self-consistent (the dGf recomputation check still passes) and puts the
transition temperatures where they belong instead of ~0.5 K off.

Validation (`tests/phase.rs`):

| quantity | model | reference |
| --- | ---: | ---: |
| dG_fus sign | + at 272 K, - at 274 K | melts at 273.15 K |
| p_sat(283.15 K) | 0.0123 bar | 0.01228 |
| p_sat(298.15 K) | 0.0317 | 0.03169 |
| p_sat(323.15 K) | 0.1231 | 0.12352 |
| p_sat(353.15 K) | 0.4699 | 0.47144 |
| p_sat(373.15 K) | 0.9988 | 1.01325 (1.4% low) |
| triple-point p(273.15 K) | 0.0061 | 0.00611 |
| dH_fus(273.15 K) | 6.01 kJ | 6.008 |
| dH_vap(373.15 K) | ~40.7 kJ | 40.66 |

`trit-chem-cli phase` prints the equilibrium composition of 1 mol of water in a
10 L headspace from 250 to 450 K.

Melting is also modelled for NaCl (1074 K, dH_fus 28.0 kJ/mol) and Fe
(1811 K, dH_fus 13.8 kJ/mol). Their liquid rows are *derived*: the measured
heat of fusion fixes H and S at the melting point, while the Cp shape comes from
NIST's own liquid fit (1074-2500 K and 1809-3133 K) bridged down to 298.15 K.
The published 298 K values for these liquids are extrapolations that disagree
with the measured fusion enthalpy by 5-10 kJ, so they are not used.

Boiling and sublimation are the same mechanism with a gas partner
(`NaCl(l)/NaCl(g)`, `NaCl(s)/NaCl(g)`, `Fe(l)/Fe(g)`, `Fe(s)/Fe(g)` through
`rules::phase_transitions`). No boiling point is stored anywhere: the gas side
of a condensed-to-gas reaction sits at the 1 bar standard state, so the
temperature where its dG crosses zero *is* the normal boiling point, and
`exp(-dG/RT)` is the vapour pressure in bar below it.

Validation (`tests/vapor.rs`):

| quantity | model | reference |
| --- | ---: | ---: |
| NaCl Tb (dG_vap = 0) | 1791 K (+3.1%) | 1738 K |
| NaCl dH_vap(1738 K) | 166.5 kJ/mol | 171 |
| NaCl dH_sub(298.15 K) | 229.7 kJ/mol | 229.7 |
| Fe Tb (dG_vap = 0) | 3001 K (-4.3%) | 3134 K |
| Fe dH_vap(3134 K) | 358.7 kJ/mol | 349.9 |
| Fe dH_sub(298.15 K) | 415.5 kJ/mol | 416.3 |
| dG_sub - (dG_fus + dG_vap) | < 1e-9 kJ | identity |
| 1 mol NaCl(l), 1000 L headspace | all vapour at 1900 K; stops at p_sat within 2% at 1200 K | - |

The 4% iron gap is traced rather than mysterious. The entropy path for solid
iron lands 8.0 J/(mol K) below the absolute entropy of the fit that covers it
(engine 84.2 against fit 92.2 at 1811 K) because NIST publishes five
overlapping solid ranges whose seam values disagree, and the liquid row
inherits the same offset through the fusion anchor. The vapour side is right to
0.2 J/(mol K) (engine 236.5 against fit 236.7 at 3134 K). Cancelling the offset
would move iron's boiling point to about 3220 K - the same size error the other
way - so the fix is to choose one solid fit per row instead of letting the range
list decide, which is the next data task.

**Bug fixed on the way:** a NIST page for a solid also carries its liquid and
gas fits, and the generator used to keep all of them. Above the melting point
the "solid" species therefore used the liquid heat capacity, which silently
erased the phase change. `gen_thermo.py` now keeps only fits whose phase label
matches the species; `tests/fusion.rs` pins the solid Cp near the melting
point so the filter cannot be removed silently.

**Second bug, same family:** `fetch_shomate.py` dropped NIST's scientific
notation. Those cells arrive as `-1.884667×10-8` once the tags are stripped, and
the parser read `-1.884667` plus two stray integers - truncating the value by
eight orders of magnitude *and* shifting every column after it. The exponent
handling existed in the script as a helper that no call site used, so
`shomate.json` carried truncated liquid fits for Al2O3, CaO, Fe, SiO2 and Zn.
Nothing shipped was wrong only because phase filtering discarded those blocks;
the boiling milestone is the first thing that needed a liquid fit.

Limits: ice Ih only (no other polymorphs); fusion, boiling and sublimation for
two substances only; no volume change on melting, boiling or freezing (the
liquid volume is held fixed); no nucleation or bubble barrier, so the cascade
jumps straight to the stable phase; the gas is ideal, which is why the
saturation pressure comes out 1-3% low near the boiling point.

## 5. Audit of the AE-ENGINE v7.1 F0 numbers

The draft claims the F0 acceptance scenario is
`钠+盐酸→H₂（ΔG=-402 优先），HCl 耗尽后钠+水（ΔG=-368）`.

Recomputed from NIST-JANAF / CODATA standard values at 298.15 K:

| Reaction | Draft | Recomputed |
| --- | ---: | ---: |
| 2Na(s) + 2H+(aq) → 2Na+(aq) + H2(g) | -402 | **-523.82** per extent; **-261.91 per mol Na** |
| 2Na(s) + 2H2O(l) → 2Na+(aq) + 2OH-(aq) + H2(g) | -368 | **-364.04** per extent; **-182.02 per mol Na** |
| H+(aq) + OH-(aq) → H2O(l) | — | -79.89 |

The ordering claim (acid before water) holds: the acid reaction is more negative.
The magnitudes do not. `-402` is close to the formation Gibbs energy of NaCl(aq)
itself (-393.1), and `-368` is close to twice the water reaction (-364.0), which
suggests the draft mixed a per-mole and a per-two-mole convention somewhere. The
engine reports both conventions explicitly so this cannot recur silently.

A second data-model error in the draft: it lists "standard electrode potential"
as a property of each of the 118 elements. A potential belongs to a half
reaction (`Cu2+/Cu` = +0.34 V, `Cu2+/Cu+` = +0.15 V); there is no single
per-element value. That data lives in `redox::HALF_REACTIONS`.

Third: "118 elements × valence electrons" conflates two different quantities.
The table stores `outer_shell_electrons` (the outermost shell) and
`oxidation_states` separately; d-block elements have more valence electrons than
their outer shell population.

## 6. What "first principles" can and cannot mean here

The draft says "reaction/function is computed, not looked up" and rejects
preset reaction tables. The engine honours the spirit, with one unavoidable
qualification:

* **Stoichiometry and feasibility are computed.** Candidate reactions come from
  rules; products are never keyed off a recipe. `dG < 0` decides.
* **The species themselves must exist.** A species' `dHf`, `dGf` and `S` are
  measured or computed data. No engine can derive the formation energy of an
  arbitrary compound at runtime from element properties alone. The honest options
  are: a curated database (chosen), group additivity estimates (Benson/Joback,
  ±5–20 kJ/mol), or offline quantum chemistry (DFT), which is the only way to
  cover arbitrary new molecules.

So "no preset reaction table" is achievable; "no preset species data" is not,
and the draft should say so.

## 7. Roadmap

M1 (correctness within the current model):

1. ~~Activity coefficients: Debye-Huckel / Davies first~~ (done in 0.2.0); next: Pitzer/SIT for brines.
2. Ion Cp data, or a documented "298 K only" mode that refuses to extrapolate.
3. ~~Phase-change species and transition reactions (ice/water/steam, molten
   salts)~~ (done through 0.5.0: water plus NaCl and Fe melting, boiling and
   sublimation). Next: choose one NIST solid fit per row so iron's entropy stops
   losing 8 J/(mol K) at its seam, then cover more substances (Cu, Al, Si).
4. Data expansion: the ~200 species a survival/crafting game actually touches
   (fuels, ores, alloys, acids, gases), with a per-species source field.
5. ~~Replace the CC BY-SA element source~~ (done in M0.1: mendeleev 1.3.0, MIT).

M2 (model upgrades):

1. White-Johnson-Dantzig element-potential Gibbs minimisation with condensed
   phase complementarity, replacing the forward-only cascade.
2. Implicit/adaptive ODE integration for stiff kinetics.
3. Group-additivity estimator for unknown organics, with uncertainty carried into
   the ranking.
4. Optional CALPHAD-style solution phases for alloys (smelting is a game system,
   not a lab).

M3 (game integration):

1. `no_std` core and a stable serialised state.
2. Fixed-point / deterministic mode if the game's netcode needs bit-identical
   replay.
3. Spatial coupling: the cascade is per-cell; heat and mass transport between
   cells belongs to the physics/environment layer, not here.
