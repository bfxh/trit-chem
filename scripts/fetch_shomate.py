#!/usr/bin/env python3
"""Fetch Shomate coefficients from NIST WebBook (NIST-JANAF data).

Writes data/external/shomate.json: species id -> list of blocks
{label, lo, hi, a..h}.  Every block is kept; the engine picks the one that
covers the requested temperature.
"""
import json
import pathlib
import re
import sys
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "external" / "shomate.json"
UA = {"User-Agent": "Mozilla/5.0 (research; trit-chem data build)"}

# (output id, NIST cbook ID, mask)  mask: 1 = gas, 2 = condensed
TARGETS = [
    ("H2(g)", "C1333740", 1), ("O2(g)", "C7782447", 1), ("N2(g)", "C7727379", 1),
    ("Cl2(g)", "C7782505", 1), ("H2O(g)", "C7732185", 1), ("CO2(g)", "C124389", 1),
    ("CO(g)", "C630080", 1), ("CH4(g)", "C74828", 1), ("NH3(g)", "C7664417", 1),
    ("SO2(g)", "C7446095", 1), ("HCl(g)", "C7647010", 1),
    ("H2O(l)", "C7732185", 2), ("Na(s)", "C7440235", 2), ("NaCl(s)", "C7647145", 2),
    ("NaOH(s)", "C1310732", 2), ("Fe(s)", "C7439896", 2), ("Cu(s)", "C7440508", 2),
    ("Zn(s)", "C7440666", 2), ("C(s,graphite)", "C7782425", 2),
    ("SiO2(s)", "C7631869", 2), ("Al2O3(s)", "C1344281", 2), ("CaCO3(s)", "C471341", 2),
    ("H2SO4(l)", "C7664939", 2), ("S(s,rhombic)", "C7704349", 2),
    ("Fe2O3(s)", "C1309371", 2), ("Fe3O4(s)", "C1317619", 2),
    ("CaO(s)", "C1305788", 2), ("Ca(OH)2(s)", "C1305620", 2),
]

SECTION = re.compile(r"(Gas|Liquid|Solid) Phase Heat Capacity \(Shomate Equation\)")


def fetch(cid, mask):
    url = f"https://webbook.nist.gov/cgi/cbook.cgi?ID={cid}&Units=SI&Mask={mask}"
    req = urllib.request.Request(url, headers=UA)
    return urllib.request.urlopen(req, timeout=45).read().decode("utf-8", "ignore")


def plain(html):
    txt = re.sub(r"<[^>]+>", " ", html)
    txt = txt.replace("&nbsp;", " ")
    return re.sub(r"[ \t]+", " ", txt)


def blocks(html):
    """Parse NIST Shomate tables.

    The table lists every fitted range as one column:
        Temperature (K) | 100. to 700. | 700. to 2000. ...
        A | v1 | v2 | v3
    """
    out = []
    marks = list(SECTION.finditer(html))
    for i, m in enumerate(marks):
        lab = m.group(1).lower()
        end = marks[i + 1].start() if i + 1 < len(marks) else len(html)
        sec = html[m.end():end]
        for chunk in re.split(r"Temperature \(K\)", sec)[1:]:
            txt = plain(chunk)
            ranges = re.findall(r"([\d.]+)\.?\s*to\s*([\d.]+)\.?", txt)
            if not ranges:
                continue
            cut = txt.find("Reference")
            body = txt[:cut] if cut > 0 else txt
            cols = {}
            letters = "ABCDEFGH"
            for j, k in enumerate(letters):
                start = re.search(r"(?<![A-Za-z])" + k + r"(?![A-Za-z])", body)
                if not start:
                    cols = {}
                    break
                nxt = body[start.end():]
                stop = re.search(r"(?<![A-Za-z])" + letters[j + 1] + r"(?![A-Za-z])", nxt) if j + 1 < 8 else None
                seg = nxt[:stop.start()] if stop else nxt
                nums = [float(x) for x in re.findall(r"-?\d+\.?\d*", seg)]
                if len(nums) < len(ranges):
                    cols = {}
                    break
                cols[k.lower()] = nums[:len(ranges)]
            if len(cols) != 8:
                continue
            for c, (lo, hi) in enumerate(ranges):
                out.append({"label": lab, "lo": float(lo), "hi": float(hi),
                            **{k: cols[k][c] for k in cols}})
    return out


def main():
    data = {}
    for sid, cid, mask in TARGETS:
        try:
            data[sid] = blocks(fetch(cid, mask))
            print(f"{sid:14s} {[(b['label'], b['lo'], b['hi']) for b in data[sid]]}")
        except Exception as e:  # noqa: BLE001 - report and continue
            print(f"{sid:14s} FAIL {e}", file=sys.stderr)
            data[sid] = []
    OUT.write_text(json.dumps(data, indent=1), encoding="utf-8")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
