# Refractory eating gates `food_min` and `mouthful_min` (2026-10-22)

Branch `food-thresh`. Prompted by N after the spiral-wave finding (`Research/Reaction-diffusion/plans/2026-10-21_spirals.md`): the only spirals in EvoCA need the legacy uncapped mouthful (`food_cap=0`, defect A9), whose negative-food trough behind each front acts as a refractory zone. N proposed a capped rule with a refractory gate, M'(x) = 0 for F(x) < θ, otherwise M(x), and chose to implement both readings of θ as separate parameters.

## Plan

1. C: two metaparameters, both default 0 (off): `food_min` (no eating while F(x) < food_min) and `mouthful_min` (no mouthful, after the caps, below mouthful_min). Each acts only when its threshold is > 0, under either `food_cap` rule.
2. Python: `init` / `update_*` / `get_*`, `params()`, `_DEFAULTS`, recipe `metaparams_init` and `metaparams_final`, `evoca_explore._PARAM_KEYS` and `run_sim` echo (`param_food_min`, `param_mouthful_min`), gallery table columns. `Docs/model.md`.
3. Tests through the public doors; each shown to fail without the change.
4. Experiment: can either gate give refractory waves or spirals under the corrected rule?

## Design decisions (flag for N)

- **No error for `food_cap=0` with a gate.** I had proposed refusing that combination. Overrides are applied key by key (`evoca_nb`, `Entry.load`), so such a check would fire or not depending on key order. Instead each gate is guarded by `θ > 0` and applies under both rules (under the legacy rule it also blocks eating at negative-F sites). θ = 0 is bit-for-bit today's code in both rules.
- **`mouthful_min` relighting.** The mouthful is `min(k·F, F, 1 − f)`, `k = m_scale·score/25`, so a site relights at `F ≥ mouthful_min` for strong eaters (`k ≥ 1`) and at `F ≥ mouthful_min/k` for weaker ones; a cell with headroom below `mouthful_min` pauses.

## Summary of what was done

- `C/evoca.c`, `C/evoca.h`: globals `gfood_min`, `gmouthful_min`, setters clamp at 0, one gate after the cap in Phase 3. Downstream (chem store, lineage intake, eg_food, accounting) sees a gated cell as having eaten nothing.
- `python/evoca_py.py`, `evoca_explore.py`, `evoca_gallery.py`; `Docs/model.md` (parameter table, eating rule, Phase 3, setter list).
- `tests/test_food_thresh.py` (16): gates off == main b374612 (SHA1 fingerprints, capped and legacy); each gate changes the run; a constructed one-step fixture (no regrowth, diffusion, tax or mutation, so F_before − F_after is each site's mouthful) checks `food_min` and `mouthful_min` site by site against the ungated step; F ≥ 0 and the energy identity with gates on; defaults and re-init reset; recipe round trip (init/final); `run_sim` door (warnings as errors, so a key missing from `_PARAM_KEYS` fails); gallery door (`Entry.load` override and table); `food_min` under the legacy rule. With the C gate disabled, 6 behavior tests fail; with the keys removed from `_PARAM_KEYS`, the `run_sim` test fails. Full suite: 178 passed.

## Results: no spirals or refractory waves from either gate (example runs, N=512, 2500 steps, seed 0, `food_cap=1`)

| Setting | Gate values | Outcome |
|---|---|---|
| Spiral recipe (R1 row 208: `m_scale` 3.5, `food_inc` 0.016, `tax` 0.035, `gdiff` 0.08) | `food_min` 0-0.7; `mouthful_min` 0.05-0.5 | uniform grazing texture, alive 0.32-0.39, std < 0.015 |
| Same, slower refill and higher tax | `food_min` 0.3/0.5/0.7 × `gdiff` 0/0.02/0.08 × `tax` 0.035/0.06 | viable only at `gdiff` 0.08, `tax` 0.035 (uniform); all 15 others extinct |
| `cfg102` (`m_scale` 0.8, patchy waves) and boom-bust recipe | `food_min` 0-0.5; `mouthful_min` 0.02-0.1 | waves at θ ≤ 0.1, damped as θ rises (std 0.15 → 0.005 by `food_min` 0.5); no rotation |
| Wider grid | `m_scale` 1.5/2/2.5/3.5 × `food_inc` 0.010-0.020 × `food_min` 0.1/0.2 | all 32 uniform |

**Why (diagnosis, partly measured).** Under the capped rule live cells sit on nearly empty sites (median F under live cells 0.013 at the spiral recipe) and live on regrowth plus food diffusing in from rich neighbors (about 0.03 per step at `gdiff` 0.08, twice `food_inc`). A gate makes each cell wait for its own site and then gulp, so every cell cycles on its own and nothing propagates. Slowing refill enough to make the wait outlast a cell's reserve (`f/tax`) starves the whole population instead. The legacy debt did something the gates do not: it **killed** cells that entered a trough (negative mouthfuls drain `f`), so troughs were emptied and later re-invaded by a coherent front. The missing ingredient is a refractory zone that is cleared of life, not merely one where food is withheld.

## Status (decided by N, 2026-10-23)

- **Conclusion.** The standard model has no mechanism that creates a refractory period, so it makes no spiral (scroll) waves; this is a conclusion about the model class. Refractory tweaks that violate the simplicity of the model are unacceptable. That rules out a starvation cost (extra tax where the environment is low in food), which has no modeling justification. The only acceptable path is the waste/inhibitor field of the chem model (chem milestone 2).
- **Branch `food-thresh` (06a181f) is not merged.** It is kept as the record of this experiment; the gates do not belong in the model.
- The notebook restructure (capped rule as standard; legacy results in their own section) is `Docs/2026-10-23_food-cap-standard.md`.
