# food_cap impact (A9 fix), 2026-10-15

Branch `food-cap`. This study pairs the legacy uncapped mouthful
(`food_cap=0`, the rule behind every pre-fix run) with the capped mouthful
(`food_cap=1`) on the same seeds. It is deliberately small: N=128, 8000
ticks, 4 seeds, 24 configs, 192 runs, 0 errors, 442 s locally on 14
workers. See `scan.py` for the config groups. The `c_*` columns come from
a companion run with the shadow off, so they follow a different
trajectory with the same params and seed.

**Finding (1 line):** the R1 "productive `m_scale` interior 2.5–3.5" does
**not** survive the cap. Under the cap, viability is 100% for `m_scale`
1–5 at all three `food_inc` values. Density, food stock and birth rate
change only smoothly and weakly with `m_scale`, `dyn_excess_pc_slope`
is flat (13–18) and `eg_excess_pc_slope` rises monotonically up to the
grid edge at 5.0. The interior was built by A9 at both ends:

- **Low/intermediate end.** Under food scarcity (`food_inc` 0.010/0.013),
  legacy runs die at `m_scale` 2.0–3.0 (viability 0–0.75). This is R1's
  "U-shaped viability, worst at intermediate m_scale". Capped runs at the
  same configs are 4/4 viable, so that U is an A9 artifact.
- **High end.** At legacy `m_scale` 5.0, 17–21% of sites hold F<0 and
  births run 25–35× the capped rate (debt-fed churn). Specificity is
  pushed to ~17.5 and `eg_excess` collapses to 9–14, against 76–95
  capped. That collapse is what made 5.0 look unproductive in R1.

**Named configs.** The **pure-evo joint winner** (`m_scale` 2.5, `food_inc`
0.018) is essentially unchanged: legacy F<0 on 1.2% of sites, density
0.392→0.413, both excess slopes within seed noise. **R1 rows 225 and 174**
(`m_scale` 3.5) change little: excess is slightly higher under the cap and
births fall about 20%. In the **#3c/#8 regime** (`m_scale` 3.0, `food_inc`
0.015, `tax` 0.035) mean density is unchanged (0.378 vs 0.370), but births
per cell fall 40% (0.037→0.022) and `dyn_excess_pc_slope` rises from 3.8 to
6.2. The final-density snapshot differs a lot (0.65 vs 0.20) because this
regime oscillates. Use the mean, not the final value.

**Caveats.** N=128 rather than the scans' 256, so slope magnitudes are not
comparable to R1 in absolute terms. There are 4 seeds and a single base
genome (the joint-winner axes), and the extinctions use run_sim's
trajectory, which has the shadow on. Under the cap, eating saturates once
`m_scale·score ≥ 25`: the cell then eats all of F(x). So `m_scale` above
~25/score is not a stronger eater; it only lowers the score threshold for
saturation. That is why `cog_specificity` *falls* with `m_scale` under
the cap (11.5 at m_scale 1 → 4.2 at 5.0, `food_inc` 0.010). Whether the rising `eg_excess`
there is productive novelty or relaxed selection still needs
interpretation.

**Needs re-running under `food_cap=1` (not launched):**
- R1 `2026-05-18_R1_mscale_bracket`: all of it. Its headline conclusions
  (interior `m_scale` optimum, U-shaped viability, mu_lut × m_scale
  interaction) are suspect.
- Pure-evo #1/#2 (`2026-05-17_pure_evo_joint`, `_lut_only`,
  `_egene_on_winners`; `m_scale` up to 2.5): the low-`food_inc` cells at
  `m_scale` ≥ 2 changed most; the winner region at `food_inc` 0.018 hardly
  changed.
- #3c and #8 (`m_scale` 3.0): demography (births) changed, so the #8
  lineage statistics and the #3c paired C−A deltas need re-checking.
- Recipes with `m_scale` > 1: 28 RD recipes at 1.2 (mild; F<0 is rare at
  1.2–1.5), 9 at 1.5, Activity-metrics 2.5 (×2), Ring 2.5, Resource-
  driving 3.5, Criticality `2026-10-15_foodinc013_osc` (3.5, so the
  boom-bust oscillation should be re-checked under the cap).
- A5's nominal-vs-actual tax gap should shrink toward the 0.4% seen at
  `m_scale`=1.
