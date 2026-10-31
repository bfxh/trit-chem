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
* Aqueous ions use `Cp = 0`. Standard partial molar heat capacities of ions are
  small and frequently negative; pretending they are zero is better than
  inventing values, but it means ion thermodynamics is a 25 C tool.

### 2.4 Standard states

Gas: ideal, 1 bar. Solute: 1 mol/L. Pure condensed phase: unit activity while
present. Gas partial pressure uses the ideal gas law with the headspace volume.

## 3. Activity model

`a = p/p0` for gases, `a = c/(1 M)` for solutes, `a = 1` for present pure
condensed phases. No activity coefficients.

Consequence: ionic-strength effects are missing. At 0.1 M the mean ionic
activity coefficient of a 1:1 salt is roughly 0.75, worth ~0.7 kJ/mol in `RT ln
a` per ion; at 1 M it can be several kJ/mol and changes qualitative behaviour
(undersaturated vs supersaturated). The demo stays at 0.05–0.1 M and should be
read with that in mind.

## 4. Audit of the AE-ENGINE v7.1 F0 numbers

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

## 5. What "first principles" can and cannot mean here

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

## 6. Roadmap

M1 (correctness within the current model):

1. Activity coefficients: Debye-Huckel / Davies first, Pitzer for brines.
2. Ion Cp data, or a documented "298 K only" mode that refuses to extrapolate.
3. Phase-change species and transition reactions (ice/water/steam, molten salts).
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
