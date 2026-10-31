# trit-chem repository rules

## Gates (run both before every commit)

```text
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
python scripts/module_budget.py
```

## Rules

1. **std only.** Core has no dependencies. Data lives in generated tables, not
   in ad-hoc constants inside algorithms.
2. **Every number has a source.** A `Species` row without a `source` string, or
   a constant without a comment naming its compilation, is a bug.
3. **No silent physics.** If a model is an approximation (ideal solution, zero
   ion Cp, forward-only cascade, extrapolated Shomate range), the doc comment
   says so and `docs/DESIGN.md` lists it.
4. **Tests must be decisive.** A test that would still pass if the sign, the
   reference-state basis or the element-entropy divisor were flipped is not a
   test. The acceptance test must fail if `dGf` anchoring or per-atom entropy is
   reverted.
5. **Hand-written files stay under 800 lines** (`scripts/module_budget.py`).
   Generated `elements_data.rs` / `species_data.rs` are exempt.
6. **Data licences must stay compatible with MIT.** Element facts come from
   mendeleev (MIT); NIST data requires attribution. Do not vendor a CC BY-SA
   table. See `docs/DATA.md`.
