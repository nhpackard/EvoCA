# 2026-10-23: push_edges rows at m_scale 1.2, both eating rules

The 65 rows of `2026-04-27_push_edges` with `m_scale` 1.2 are the only 04-27 rows where the mouthful cap can bind (`m_scale > 1`). Each row was re-run with today's code and the push_edges protocol (N=256, 5000 steps, sample_every 25, seed 0, shadow on, halfplane) under `food_cap=1` and `food_cap=0`, 130 runs.

**Purpose.** To separate the effect of the rule from the change in the model since April. Recorded vs. legacy-now is the model change; legacy-now vs. standard is the rule.

**Findings (Reaction-diffusion notebook, Legacy rule section, part 3).**
- The median correlation length goes 20.9 (recorded) -> 10.3 (legacy now) -> 9.9 (standard). The model change is far larger than the rule.
- The rule matters for the wave configs. Five rows wave under the legacy rule (alive std > 0.1); three of them drop below 0.1 under the cap, e.g. cfg246 0.18 -> 0.09.

The `food_cap=1` rows are identical to the same rows of `2026-10-23_rescan_0427_standard`.
