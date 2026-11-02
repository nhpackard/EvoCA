# 11 — `python/evoca_nb.py`: shared lab-notebook helpers (2026-10-18)

**Request (N).** "Tools used by the example runs" should live in a Python module; the six
lab notebooks each defined private copies of the same helpers (load a recipe reproducibly,
build a sim from a scan row, ARGB→RGB frames, export a recipe).

## What was built

`python/evoca_nb.py`:

| Function | Replaces |
|---|---|
| `rgb(sim, colormode)`, `argb_to_rgb(a)` | four private `rgb`/`render` copies (RD, Resource, Lineage, Criticality) |
| `load_recipe(name, N, seed, recipe, seed_v, shadow, **overrides)` | Activity `load_recipe`, RD `load` (RD keeps `seed_v=False`, its old behavior) |
| `sim_from_scan(scan, config_idx, N, seed, **overrides)` | Resource `r1_sim` (recipe written to a temp dir that is removed; no more `mkdtemp` litter) |
| `params_from_row(row)` | Lineage `params_from_row` |
| `make_sim(N, seed, **params)` | Ring `make_sim` |
| `export_recipe_once(sim, descriptor, workstream, **kw)` | Ring's glob guard; fixes Resource and Criticality writing a new dated `.evoca` on every Run All |
| `slope(t, y)` | Activity/Lineage `slope` |

`evoca_explore.build_sim` gained keyword-only `lineage=False`, `egenome='uniform'` (defaults
unchanged, pinned by a test); it now reproduces the three Lineage hand setups exactly.

**Single source for scan-row corrections.** `evoca_explore.correct_scan_row_params()` holds
the load-time corrections for params read from a scan row (today: pre-A9 rows without
`param_food_cap` get `food_cap=0`). Both `evoca_from_scan` and `evoca_nb.params_from_row` call
it, so a scan row means the same run whichever path loads it (test shown to fail without it).
Consequence: Lineage §1's example runs (pre-A9 scan rows) run uncapped, matching their stored
outputs.

`Research/Resource-driving/Runs/2026-10-15_R1cfg208_intermediate_fi0p016.evoca` now carries
`food_cap: 0`: it was saved from an uncapped run.

## Verification

- `tests/test_evoca_nb.py` (26 tests): every replaced helper is copied verbatim and compared
  with its replacement by a state fingerprint after 50 steps; `rgb` against all four old
  variants in all six colormodes; `export_recipe_once` idempotent (fails against the naive
  export). Suite: 143 passed.
- All six lab notebooks executed in scratch copies (outputs discarded; committed outputs
  untouched); all 25 interactive snippets run with a stubbed window; no files left in `Runs/`.
- Cell counts unchanged, so the cell numbers in guides 07/08 still hold.

## Open

- Reaction-diffusion runs keep the unseeded `v` field (`seed_v=False`), so they are not
  reproducible run to run; switching to `seed_v=True` would change its examples.
- `build_sim`'s signature change is a shared-harness change; a candidate for review if
  `/code-review ultra` is run on the drift-shadow merge.
