# The capped rule as the standard model (2026-10-23)

Branch `food-cap-standard`. N's decisions after the refractory-gate experiment (`Research/Reaction-diffusion/plans/2026-10-22_food-thresh.md`):

- **Conclusion about the model class.** The standard model makes no reaction-diffusion spiral (scroll) waves, because it has no mechanism that creates a refractory period.
- **Ruled out.** Tweaks that force a refractory period but violate the model's simplicity, e.g. a starvation tax on cells whose site is low in food, which has no modeling justification. The only acceptable path is the waste/inhibitor field, i.e. the chem model.
- **The default.** `food_cap=1` (no negative food) is the standard model and the default everywhere.
- **The notebook.** Restructure it to eliminate or sequester the `food_cap=0` results.

N also asked whether the boom-bust and patchy waves depend on `food_cap=0`.

## Answer: boom-bust and patchy waves exist under the standard rule

Example runs, N=256, seed 0, 3000 steps; alive fraction in the second half of each run:

| config | capped (standard) | legacy |
|---|---|---|
| boom-bust (test.ipynb corner, m_scale 1.2) | 0.05-0.71, std 0.195, period ~150 | 0.01-0.73, std 0.244 |
| foodinc013_osc (m_scale 3.5) | 0.03-0.75, std 0.248 | 0.02-0.69, std 0.224 |
| patchy cfg102 (m_scale 0.8) | 0.14-0.74, std 0.186 (identical: the cap cannot bind at m_scale <= 1) | same |
| patchy cfg246 (m_scale 1.2) | 0.36-0.52, std 0.026 (no waves) | 0.11-0.72, std 0.202 |

Over 10000 steps the capped boom-bust oscillation persists with mutation off. With mutation on, evolution damps it to a fixed point by t~6000 (window std 0.005). Under the legacy rule the oscillation grows instead. These are one-seed leads.

## Plan

1. Make every loader default to the standard rule. The C library, `init` and recipe import already did; pre-A9 scan rows still loaded `food_cap=0`.
2. Restructure the Reaction-diffusion notebook: main thread standard-rule only; one final section for legacy-rule results.
3. Replace legacy-rule scan metrics used by the main thread.
4. ISSUES, board, docs; a pre-merge test notebook.

## Design conflict found, and its resolution

Rescanning only the 04-27 rows where the cap can bind (`m_scale > 1`: 65 push_edges rows) turned out not to be enough. **The April rows describe an older model, not only the legacy rule.** Since then:

- max-match eating across egene slots and the egene-aware genome hash arrived (2026-04-28);
- the eating score became ternary egenes with ±1 scoring (2026-05-10);
- the neutral shadow got its own RNG stream (2026-10-15).

On 40 push_edges rows where the cap cannot bind, today's code gave a lower correlation length in 93% of rows (median 30.3 -> 14.3). Replacing only the 65 rows would have put them on a different scale from the rest of the pool. Resolution: rescan all 600 rows under today's standard model (`Scans/2026-10-23_rescan_0427_standard`, ~10 s per run). The 65-row both-rules rescan (`Scans/2026-10-23_push_edges_m1p2_rescan`) is kept for the legacy section, where it separates the rule from the model change.

## Summary of what was done

- **Default** (`python/evoca_explore.py` `correct_scan_row_params`):
  - A pre-A9 scan row (no `param_food_cap` column) now loads `food_cap=1`, written explicitly. When `m_scale > 1`, a one-line notice says the row's metrics were measured under the legacy rule.
  - A row that recorded `food_cap` keeps it, and `food_cap=0` can always be passed explicitly.
  - One source serves `params_from_row`, `evoca_from_scan`, `sim_from_scan` and gallery entries, so [RXD-D4](../Research/Reaction-diffusion/ISSUES.md) is settled: both doors agree.
  - Docstrings in `evoca_nb.py` and `evoca_py.py` were updated, along with `Docs/model.md`.
- **Tests:** 5 updated in `test_food_cap.py`, `test_evoca_gallery.py` and `test_evoca_nb.py`. The `evoca_from_scan` test now also checks the notice, its absence for recorded rows, and that an explicit 0 is kept. With the change reverted, those 5 fail. Full suite: 162 passed.
- **Notebook** (`Research/Reaction-diffusion/Reaction-diffusion.ipynb`, executed):
  - *Main thread, standard rule only:*
    - waves: the boom-bust recipe `2026-10-23_rd_waves_boom_bust_cap1.evoca`, cfg102 and `foodinc013_osc` capped;
    - a new section, "Spiral waves: not in the standard model", with N's conclusion, the capped 5 x 5 grid around the legacy spiral recipe (all stationary texture) and the gate result;
    - the 04-27 section on the rescan, with a recorded-vs-rescan comparison;
    - the gallery table's scan columns from the rescan;
    - the two-objective pool, picks, and observations rewritten from the new numbers.
  - *Final section, "Legacy rule (`food_cap=0`)":* (1) the 10-21 wave configs under both rules, including cfg246 and the legacy boom-bust variants; (2) the spirals (the 10-21 section, demoted, with the gate diagnosis added to its mechanism note); (3) the 65-row comparison: recorded / legacy-now / standard.
- **Galleries:**
  - new: `rd_spiral_grid_cap1`, `rd_waves_legacy`, `rd_boom_variants_legacy`;
  - redefined: `rd_waves` and `rd_boom_variants` (standard rule), and `rd_spiral_sweep` (legacy only, without the capped grid);
  - rebuilt: the pareto galleries.
- **Recipes:** `rd_waves_boom_bust_cap1`, and 12 `2026-10-23_pareto_*` recipes for the new box/wedge-only picks (`food_cap=1` recorded). The `2026-10-21_pareto_*` and `rd_waves_*` recipes stay; they record `food_cap=0`.
- **ISSUES:**
  - [RXD-D4](../Research/Reaction-diffusion/ISSUES.md) moved to Fixed (decided);
  - [RXD-Q5](../Research/Reaction-diffusion/ISSUES.md) rewritten for the standard model;
  - [RXD-Q6](../Research/Reaction-diffusion/ISSUES.md) answered (no spirals in the standard model);
  - [RXD-N6](../Research/Reaction-diffusion/ISSUES.md) re-aimed at the chem inhibitor field.

  The board was regenerated.
- **`Research/Reaction-diffusion/plans/2026-10-22_food-thresh.md`** is on this branch with its status. Branch `food-thresh` is not merged; it is kept as the record of that experiment.
- **Pre-merge test:** `tests/notebooks/food-cap-standard-test.ipynb`, 25 PASS checks.
  - The default rule at every door.
  - Rescan coverage and parameters.
  - Determinism: two independent rescans agree; a single-row re-run reproduces the rescan.
  - The notebook: no errors, and no main-thread code sets `food_cap=0` (an AST check, itself tested).

## What the rescan changed in the Pareto section

- Only 4 of the April composite top-20 remain (cfg102, 106, 148, 246). The April #1, cfg236, is no longer picked.
- Nine `initial`-scan configs enter the picks. They fill the N=128 lattice (alive ~0.97). Their correlation length (19-45) belongs to an almost full alive mask, a second artifact of the spatial score next to the halfplane transient.
- [RXD-N1](../Research/Reaction-diffusion/ISSUES.md)/N2 are noted in the notebook: re-score on post-burn-in metrics and exclude filled lattices before comparing selection methods.

## Structural risk and open items

- **The default flip reaches other workstreams.** Pre-A9 scans with `m_scale > 1` rows:
  - Activity-metrics: `evo_focus` 122/250, `pure_evo_joint` 237/300, `pure_evo_lut_only` 151/180;
  - Lineage-co-evolution: `pure_evo_egene_on_winners` 288/288;
  - Resource-driving: `R1_mscale_bracket` 300/300.

  Their notebooks load scan rows. Stored outputs are unchanged, but the next re-execution runs those rows capped against metrics recorded under the legacy rule (and the April/May model); the notice flags it. The same model-version caveat as for RD applies to any scan recorded before 2026-05-10.
- **Moved layer boundary.** The load-time correction changed meaning (from "reproduce the run" to "standard rule"). This is one of the moments where an independent review (`/code-review ultra`, user-triggered and billed) is worth considering.
- **Open:**
  - a 3-seed run of the capped boom-bust recipe (evolution damping the waves);
  - [RXD-N2](../Research/Reaction-diffusion/ISSUES.md) re-scoring with a filled-lattice exclusion;
  - whether to rescan the other workstreams' pre-A9 scans.

---

# Part 2: every research workstream on the standard model (2026-10-23, later)

N asked for every research notebook to be rescanned where necessary, with results that use `food_cap=0` sequestered as in Reaction-diffusion. Scans ran on torque.

## Why every notebook had to change

- **Stored outputs predate the fix.** All research notebooks were last executed on 2026-10-15 around 16:30 UTC, before the A9 fix (19:28 UTC). So every stored example run with `m_scale > 1` was a legacy-rule result.
- **Scans with `m_scale > 1` rows and no `food_cap` column:**
  - Activity-metrics: `evo_focus`, `pure_evo_joint`, `pure_evo_lut_only`, `causal_control` v1/v2.
  - Lineage: `pure_evo_egene_on_winners`, `3c_coevolution_substrate`.
  - Resource-driving: `R1_mscale_bracket`.

## Rescans (`Research/*/Scans/2026-10-23_*_standard`, torque, about 36 minutes in all)

- **How they were made.** Each is a copy of the source driver plus a shim forcing `food_cap=1`. The April `evo_focus` driver additionally maps its pre-rename `mu_egenome` to `mu_egene`. The launcher is `tools/rescan_standard_2026-10-23.sh`, and each directory has a `notes.md`.
- **Checks.**
  - The same configs as the source, column by column.
  - 0 errors.
  - 3c `determinism_ok` on every row.
- **Not rescanned:**
  - `rd_under_pure_evo` (`m_scale` ≤ 0.6);
  - `food_cap_impact` (a deliberate both-rules study);
  - `8_lineage_coevolution` (never run).

## Notebooks (one agent per workstream; I reviewed each report)

| Workstream | What changed under the standard model |
|---|---|
| Resource-driving | **R1 reverses:** extinction 36% → 10%, no U; productive `m_scale` moves from interior (2.5–3.5) to the 5.0 grid edge; `mu_lut` optimum high → low. `food_cap_impact` shows A9 built R1's headline at both ends. New [RES-Q8](../Research/Resource-driving/ISSUES.md). |
| Criticality-bifurcation | **Extinction edge** `food_inc` ~0.012 → ~0.005 (a legacy artifact). New candidate transition: amplitude collapse between 0.015 and 0.017. Still no period-doubling. Its figures were never saved before (an `Agg` backend line); fixed. |
| Lineage-co-evolution | **#3 "3a ≫ 3b" → 3a ≈ 3b** (3b extinction 46% → 1%). **#3c:** GoL-substrate death gone (arm B 10/10 → 0/10 extinct); `dyn_excess_pc` artifact gone (21/30 → 0/30 runs > 50). C > A survives (27/30). A legacy-rule rerun on today's code reproduces the recorded outcomes, so the change is the rule, not code drift. Reducer tax gap 20–23% → 0.5% ([LIN-D3](../Research/Lineage-co-evolution/ISSUES.md) no longer blocks). |
| Activity-metrics | **pure-evo headlines hold**; "egene freeze reduces dyn" is unsupported under either rule. **Causal controls:** nulls exact; eg synergy 2.1x → 1.4x. **`evo_focus`:** winner moves from `mu_lut` 0.06 to 0.003–0.01, 16 → 0 extinct (model change plus rule). The `dyn_report` verdict on the winner flips silent-churn → functional. Section-2 example runs passed `mu_egenome` meaning Negene growth; fixed. |
| Ring-complexity | **No finding depends on the rule** (no legacy section); counts move a few percent. Found that `state()` drew `v` unseeded (below). |

## Library fix found on the way: `EvoCA.state()` initial `v`

- **The defects.**
  - `state()` drew `v(x)` from an unseeded numpy generator, so `make_sim` + `set_seed` + `state()` runs were never reproducible (the C PRNG does not touch `v`).
  - `v_density` was only recorded; the draw was always 50/50.
- **Who was affected.** Scans (`build_sim`/`run_sim`) and galleries (`load_recipe(seed_v=True)`) were not: they overwrite `v` from a seeded generator.
- **The fix.**
  - New `state(..., v_seed=None)`, recorded in recipes.
  - `v_density` now sets the density.
  - Default behavior is otherwise unchanged, so the C RNG stream and the golden fingerprints are untouched.
- **Tests.** `tests/test_state_v_seed.py` (5). They fail without the change; the density cases fail with the old 50/50 draw.

## Tests

- `tests/notebooks/food-cap-standard-test.ipynb`, 47 PASS. It now also checks every rescan against its source, and every research notebook and digest: no errors, no main-thread `food_cap=0` (AST), and no stale intro.
- Full suite: 167 passed.

## Open

- The project memory note on the `dyn_excess_pc` artifact is updated: it was a legacy-rule result in #3c, and the guard is kept.
- Historical devlog entries (`Docs/devlog/2026-05-17_research-campaigns.md`, `research_board_archival.md`) still state the legacy headlines as recorded. They are left as history; add a pointer if N wants one.
- Lineage added `Scans/2026-10-23_legacy_now_subset` (legacy rule, today's code) as its control.
