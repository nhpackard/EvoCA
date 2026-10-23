# 2026-10-23: legacy rule with today's code, a subset of #3 and #3c

**Why.** The standard-model rescans of #3 and #3c (`2026-10-23_*_standard`) differ from the recorded scans in two ways: the eating rule (`food_cap=1` vs the legacy `food_cap=0`) and the code since May (mainly the neutral shadow's own RNG stream, 2026-10-15). This subset reruns the legacy rule with today's code, so recorded vs legacy-now is the code change and legacy-now vs standard is the rule.

**What.** `scan.py` imports the two rescan drivers unchanged and forces `food_cap=0` instead of 1. Rows (31, local, 14 workers, ~6 min):
- #3 (`src=3`): `config_idx` 120-127, 168-175, 264-271 (base regimes 5, 7, 11, the first four egene combos, both modes; every 3b row here was extinct in the recorded scan). N=256, 8000 steps, seed 0, as the source.
- #3c (`src=3c`): arm B seeds 0-3, and the A/C pair at phase-1 `mu_lut=0.01` seeds 0-2. Same protocol as the source.

**Result.** The legacy rule reproduces the recorded outcomes with today's code: 12/12 3b rows extinct (0/12 standard), 4/4 B extinct (0/4 standard), and the A-arm `dyn_excess_pc_slope` artifact in 2 of 3 A runs (+413, +7, +484; recorded +675, +425, +723; standard -0.6, -1.3, -2.1). So the changes between the recorded scans and the rescans are due to the rule, not the code drift.

**Used by:** the "Legacy rule" section of `Lineage-co-evolution.ipynb`.
