# rd_under_pure_evo (#4)

Take the 2026-04-27 evo+spatial (RD) balanced-composite winners and
run them long (12k ticks) under pure evolution; ask whether their
RD-scale spatial structure *persists* (final≈mean) or *decays*
(final≪mean). 10 distinct winner recipes × 3 seeds = 30 runs,
0 extinct, 89 s on torque.

**Finding (1-line):** RD spatial structure is **not robust under
pure-evo** — `correlation_length` final/mean median **0.23** (all 30
runs decay), so the RD structure was an artifact of the 2026-04-27
*spatial* objective, not a pure-evo attractor. Two regimes: 8/10 keep
dynamic turnover (`dyn_excess` 12–27) while spatial correlation washes
out; 2/10 (`cfg49`,`cfg30`: m_scale 0.6 / mu_lut 0.01 / food_inc 0.025)
collapse to a **static colony-filling blob** (largest_patch f/m≈1.00,
~70 % of grid frozen) and have the **lowest** `dyn_excess` (~5.8) —
exactly the "suspect the metrics" degenerate cases NP predicted.

Full analysis → `../2026-05-17_pure_evo_analysis.md`
Reasoning trail → `../../Docs/devlog/2026-05-17_research-campaigns.md`
