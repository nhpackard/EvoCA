# food-cap branch: A9 mouthful cap (2026-10-15)

N approved: "cap the mouthful at available food on its own branch". This
branch, `food-cap`, comes off main `a0798d1` and is not merged; it waits on
N's review.

## Plan
1. Cap the Phase-3 mouthful at the food present, behind a `food_cap` toggle
   (default 1 = capped; 0 = legacy) so past runs can still be reproduced.
2. Plumb the toggle through C setter/getter → `EvoCA.init/params/_DEFAULTS/
   update_food_cap/get_food_cap` → `_PARAM_KEYS` → recipe export/import,
   with a legacy-recipe notice.
3. Tests, parity, a small impact study, a `Docs/model.md` update, and a
   pre-merge notebook.

## What was done
- **C** (`C/evoca.c` Phase 3): with `gfood_cap`, `mouthful = min(mouthful,
  max(F,0), max(1−f,0))`. Otherwise the legacy headroom-only clamp applies.
  All downstream accounting (chem store via `chem_cell_add`, lineage
  intake, `eg_food`, `g_sum_mouthful` / dissipation `mouthful_step`) reads
  the same realized `mouthful`, so no accounting needed changing. No
  other path withdraws food: regeneration only adds and clamps at 1;
  diffusion is a convex average; tax death clamps f to 0; the reproduction
  split halves a store ≥ 1.
- **Python**: `food_cap` in `init()`, `params()`, `_DEFAULTS`,
  `_init_metaparams`, recipe `metaparams_final`, `update_food_cap()`,
  `get_food_cap()`, and `_PARAM_KEYS`. `run_sim` always echoes the
  effective `param_food_cap`. `import_run` of a recipe with no `food_cap`
  key runs capped and prints a one-line notice when `m_scale > 1`.
- **Tests**: `tests/test_food_cap.py` (20 tests). The regression test
  fails without the fix (F<0 at step 2), and the run_sim reachability test
  fails without the `_PARAM_KEYS` entry. `tests/chem_parity.py` pins
  `food_cap=0`, because its golden predates the fix and runs at m_scale=3.
- **Parity**: `tools/parity_probe.py run <root> <out> 0` gives IDENTICAL /
  BAD 0 against main. At m_scale ≤ 1, capped == legacy bit-for-bit. The
  probe fixture's m_scale=1.5 means its default path changes (BAD 120).
- **Impact**: `Research/Resource-driving/Scans/2026-10-15_food_cap_impact/`.
  The R1 productive `m_scale` interior does not survive the cap (see
  notes.md there).
- **Notebook**: `tests/notebooks/food-cap-test.ipynb`.

## Open for N
- Default capped (as implemented) or legacy?
- Re-run list in the impact notes (R1, pure-evo #1/#2 low-food_inc cells,
  #3c/#8, m_scale>1 recipes, criticality oscillation).
- Pre-existing on main: `EvoCA.grouped_params()` already raises, because
  `food_per` and `n_chem` are in neither `_ENV_KEYS` nor `_AGENT_KEYS`.
  `food_cap` (a model-rule switch) adds a third unclassified key.
