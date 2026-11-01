#!/usr/bin/env python3
"""Recompute the liquid-row anchors (dHf298, S298, dGf298) for NaCl(l) and Fe(l).

Mirrors the Rust engine and scripts/gen_thermo.py: H(T) = dHf298 + h_rel(T),
S(T) = S298 + integral of Cp/T from 298.15 (composite Simpson, 48 intervals),
with the same generated bridge segment in front of fits that start too high.

Anchors are chosen so that at the measured melting point
    dH_fus(Tm) = the measured latent heat  and  dS_fus(Tm) = dH_fus / Tm,
which is what puts dG_fus(Tm) at zero.
"""
import json
import math
import pathlib
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace", write_through=True)
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gen_thermo as G  # noqa: E402  reuse the generator's own model

T0 = G.T0
SHOM = json.loads((pathlib.Path(__file__).resolve().parents[1] / "data" / "external" / "shomate.json").read_text(encoding="utf-8"))


def row_blocks(sid, phase, s298):
    blocks = [b for b in G.dedupe(SHOM.get(sid, [])) if b["label"] == phase]
    return G.with_bridge(blocks, s298)


def pick(blocks, t):
    def dist(b):
        return max(b["lo"] - t, t - b["hi"], 0.0)
    for b in blocks:
        if b["lo"] <= t <= b["hi"]:
            return b
    return min(blocks, key=dist)


def cp_at(blocks, t):
    return G.cp_at(pick(blocks, t), t)


def h_rel(blocks, t):
    return G.h_rel(pick(blocks, t), t)


def s_rel(blocks, t):
    """Mirror of Cp::s_rel: one composite Simpson pass, 48 intervals."""
    a, b, n = T0, t, 48
    h = (b - a) / n

    def f(x):
        return cp_at(blocks, x) / x
    total = f(a) + f(b)
    for k in range(1, n):
        total += (4 if k % 2 else 2) * f(a + h * k)
    return total * h / 3.0


# solid row: (id, dHf298, S298, element entropy sum for the formula)
SOLIDS = {
    "NaCl": ("NaCl(s)", -411.15, 72.13, 51.21 + 0.5 * 223.07),
    "Fe": ("Fe(s)", 0.0, 27.28, 27.28),
}
FUSION = {"NaCl": (1074.0, 28.0), "Fe": (1811.15, 13.8)}

for name, (sid, hf_s, s_s, s_el) in SOLIDS.items():
    tm, dh_fus = FUSION[name]
    sb = row_blocks(sid, "solid", s_s)
    lb = row_blocks(f"{name}(l)", "liquid", 80.0)  # s298 only feeds the bridge's G
    bridged = len(lb) > len([b for b in SHOM.get(f"{name}(l)", []) if b["label"] == "liquid"])
    h_s, s_s_at = hf_s + h_rel(sb, tm), s_s + s_rel(sb, tm)
    dh298 = h_s + dh_fus - h_rel(lb, tm)
    ds298 = s_s_at + dh_fus * 1000.0 / tm - s_rel(lb, tm)
    # the bridge's G constant is cosmetic, but regenerate it consistently
    lb = row_blocks(f"{name}(l)", "liquid", ds298)
    dgf = dh298 - T0 * (ds298 - s_el) / 1000.0
    print(f"{name}(l): bridged={bridged} Tm={tm} Cp(l,Tm)={cp_at(lb, tm):.3f} Cp(l,298)={cp_at(lb, T0):.3f}")
    print(f"   H_s(Tm)={h_s:.3f} S_s(Tm)={s_s_at:.3f} h_rel_l(Tm)={h_rel(lb, tm):.3f} s_rel_l(Tm)={s_rel(lb, tm):.3f}")
    print(f"   hf298={dh298:.3f}  s298={ds298:.3f}  dgf298={dgf:.3f}")

# what the vapour rows imply, with the same model
for sid, tb in [("NaCl(g)", 1738.0), ("Fe(g)", 3134.0)]:
    vb = row_blocks(sid, "gas", 0.0)
    print(f"{sid}: bridged={len(vb) > 1} Cp(298)={cp_at(vb, T0):.2f} Cp(Tb)={cp_at(vb, tb):.2f} "
          f"H(Tb)-H298={h_rel(vb, tb):.2f} kJ")
