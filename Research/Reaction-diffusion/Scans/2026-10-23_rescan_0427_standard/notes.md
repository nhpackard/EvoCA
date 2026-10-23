# 2026-10-23: rescan of the three 2026-04-27 RD scans under the standard model

**Why.** The 04-27 scans (`initial`, `refined`, `push_edges`; 600 rows) ran the legacy uncapped mouthful (`food_cap=0`) and an older model. Since then: max-match eating across egene slots and the egene-aware genome hash (2026-04-28), ternary egenes with +/-1 scoring (2026-05-10), and the neutral shadow's own RNG stream (2026-10-15). On 40 push_edges rows where the cap cannot bind (`m_scale <= 1`), today's code gave a lower `correlation_length_mean` in 93% of rows (median 30.3 -> 14.3). So the recorded metrics do not describe the current model, whatever the rule. N made the capped rule the standard model on 2026-10-23.

**What.** Every row re-run with its own parameters (`evoca_nb.params_from_row`: the pre-rename `param_mu_egenome` column read as `mu_egene`) and `food_cap=1`, using the protocol of its source scan: N=256, 5000 steps, seed 0, shadow on, halfplane start, and the scan's own `sample_every` (100 for initial/refined, 25 for push_edges). Columns `scan` and `config_idx` identify the source row. 600 rows, no errors, ~10 s per run.

**Findings (Reaction-diffusion notebook, section "2026-04-27: the three RD scans").**
- Candidate yields are 23 / 50 / 100% (recorded 17 / 50 / 98%).
- Push_edges' median correlation length is 9.7 (recorded 17.3). The rescan is lower in 71% of rows for correlation length and in 68% for `unique_top_genomes`.
- The candidate pool is now 383 (369 recorded) on 42 Pareto layers. Only 4 of the April composite top-20 remain in it.
- Deterministic: the 65 `m_scale` 1.2 rows equal those of the independent `2026-10-23_push_edges_m1p2_rescan`, and a single-row re-run reproduces its row (`tests/notebooks/food-cap-standard-test.ipynb`).

**Used by:** `load_0427()` in the notebook (default `rule='standard'`), the two-objective pool, and the gallery table's scan columns.
