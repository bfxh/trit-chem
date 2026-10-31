#!/usr/bin/env python3
"""Module budget: no hand-written Rust file above 800 lines.

Generated data files are exempt -- they are data, not modules.
"""
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
LIMIT = 800
GENERATED = {"elements_data.rs", "species_data.rs"}


def main():
    rows = []
    for p in sorted(ROOT.glob("crates/**/*.rs")):
        if p.name in GENERATED:
            continue
        n = sum(1 for _ in p.open(encoding="utf-8", errors="ignore"))
        rows.append((n, p.relative_to(ROOT).as_posix()))
    rows.sort(reverse=True)
    for n, name in rows[:10]:
        print(f"{n:5d}  {name}")
    over = [(n, name) for n, name in rows if n > LIMIT]
    if over:
        for n, name in over:
            print(f"OVER BUDGET: {name} has {n} lines (limit {LIMIT})", file=sys.stderr)
        return 1
    print(f"module budget ok: {len(rows)} files, max {rows[0][0]} lines")
    return 0


if __name__ == "__main__":
    sys.exit(main())
