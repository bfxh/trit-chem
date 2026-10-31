# Annotated bibliography

Links below were opened while writing M0. Entries marked *citation only* were
not fetched; they are standard references, listed without a URL rather than
guessing one.

## A. Thermodynamic data and temperature fits

**M. W. Chase, Jr., NIST-JANAF Thermochemical Tables, 4th ed.**
J. Phys. Chem. Ref. Data, Monograph 9 (1998).
<https://webbook.nist.gov/chemistry/> — the source behind every Shomate block
and the 25 C formation values used here. WebBook pages expose `dHf`, `S` and
Shomate coefficients per phase; `scripts/fetch_shomate.py` scrapes them.

**B. J. McBride, M. J. Zehe, S. Gordon, NASA Glenn Coefficients for Calculating
Thermodynamic Properties of Individual Species.** NASA/TP-2002-211556 (2002).
<https://ntrs.nasa.gov/citations/20020085330> — the 7- and 9-coefficient NASA
polynomial format (the standard for combustion codes) and coefficients for about
2000 species. Not used yet; it is the M1 upgrade path for species NIST lacks.

**S. Gordon, B. J. McBride, Computer Program for Calculation of Complex Chemical
Equilibrium Compositions and Applications.** NASA RP-1311 (1994).
<https://ntrs.nasa.gov/api/citations/19950013764/downloads/19950013764.pdf> —
CEA: element-potential equilibrium with condensed phases. The reference
implementation for M2.

**D. D. Wagman et al., The NBS Tables of Chemical Thermodynamic Properties.**
J. Phys. Chem. Ref. Data 11, Suppl. 2 (1982). *Citation only.* Source of the
aqueous-ion values (and of the electron-convention caveat in DESIGN.md §2.1).

## B. Chemical equilibrium algorithms

**W. B. White, S. M. Johnson, G. B. Dantzig, Chemical Equilibrium in Complex
Mixtures.** J. Chem. Phys. 28, 751 (1958).
<https://www.rand.org/pubs/papers/P1059.html> — the element-potential /
Gibbs-minimisation formulation that every modern equilibrium code descends from.
It is also the reference for why the problem is a convex program.

**W. R. Smith, R. W. Missen, Chemical Reaction Equilibrium Analysis: Theory and
Algorithms.** Wiley (1982). *Citation only.* Stoichiometric vs non-stoichiometric
formulations, and the stability questions the M0 cascade does not answer.

**W. C. Reynolds, The element potential method for chemical equilibrium
analysis: implementation in the interactive program STANJAN.** Stanford (1986).
*Citation only.* Practical Newton implementation notes.

**P. Koukkari, R. Pajarre, Introducing mechanistic kinetics to the
Gibbsian equilibrium calculation.** Comput. Chem. Eng. 30, 1189 (2006).
*Citation only.* How to constrain a Gibbs minimiser with kinetic barriers — the
principled version of what the M0 cascade approximates with a fixed order.

## C. Kinetics

**S. Arrhenius, Über die Reaktionsgeschwindigkeit bei der Inversion von
Rohrzucker durch Säuren.** Z. Phys. Chem. 4, 226 (1889). *Citation only.* The
`k = A exp(-Ea/RT)` form implemented in `kinetics.rs`.

**D. G. Goodwin et al., Cantera: an open-source suite for chemical kinetics,
thermodynamics and transport.**
<https://cantera.org/> — the reference open implementation. Useful as an oracle
for M1/M2: build the same mechanism in Cantera and compare.

**R. J. Kee et al., CHEMKIN-III** (Sandia report SAND96-8216, 1996).
*Citation only.* Mechanism files and the stiff-integrator conventions combustion
workflows expect.

## D. Solutions and electrolytes

**P. Debye, E. Hückel, Zur Theorie der Elektrolyte.** Phys. Z. 24, 185 (1923).
*Citation only.* Limiting law; the first correction to add for dilute solutions.

**K. S. Pitzer, Thermodynamics of Electrolytes. I. Theoretical Basis and General
Equations.** J. Phys. Chem. 77, 268 (1973).
<https://en.wikipedia.org/wiki/Pitzer_equations> (overview) — the virial
treatment that stays accurate into concentrated brines. Needed before seawater,
acid baths or molten-salt gameplay.

**C. W. Davies, Ion Association.** Butterworths (1962). *Citation only.* The
Davies equation, implemented in `state::activity_coefficient` and checked
against measured NaCl activity coefficients at 0.001/0.01/0.1 M.

**A. Fredenslund, R. L. Jones, J. M. Prausnitz, Group-contribution estimation of
activity coefficients in nonideal liquid mixtures.** AIChE J. 21, 1086 (1975).
*Citation only.* UNIFAC; the organic-mixture counterpart to Pitzer.

## E. Estimating unknown species

**S. W. Benson, Thermochemical Kinetics, 2nd ed.** Wiley (1976).
*Citation only.* Group additivity for `dHf`, `S` and `Cp`. The only tractable way
to invent chemically plausible organics at runtime, with ±5–20 kJ/mol error that
must be carried into the ranking.

**K. G. Joback, R. C. Reid, Estimation of pure-component properties from
group-contributions.** Chem. Eng. Commun. 57, 233 (1987). *Citation only.*
A coarser, wider-coverage alternative to Benson groups.

**RMG-Py thermochemistry estimation chapter.**
<https://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/thermo.html> —
worked description of the group-additivity database format, useful if the game
needs to grow its own species.

## F. Alloys and phase diagrams

**N. Saunders, A. P. Miodownik, CALPHAD (Calculation of Phase Diagrams): A
Comprehensive Guide.** Pergamon (1998). *Citation only.* The framework behind
real smelting/alloy gameplay; M2/M3 material.

## G. Games

**P. Purho, Exploring the Tech and Design of Noita.** GDC 2019.
<https://www.gdcvault.com/play/1025695/Exploring-the-Tech-and-Design> — the
falling-sand cellular automaton with per-pixel reactions. The design draft
explicitly rejects its preset reaction table; this talk is the best description
of what that alternative buys and costs (speed, determinism, authoring effort).

## What M0 actually takes from each

| Source | Taken |
| --- | --- |
| JANAF via WebBook | species `dHf`, `S`, Shomate blocks |
| NBS/CODATA table | aqueous ions, reduction potentials |
| White-Johnson-Dantzig | the equilibrium concept and the M2 target |
| Arrhenius | rate law |
| Cantera | intended cross-check, not yet wired |
| Noita GDC | the rejected design, kept as a performance yardstick |
