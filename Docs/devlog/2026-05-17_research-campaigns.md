# 2026-05-17 — research campaigns (pure-evolution program + others)

Compute: torque.protolife.com (32 cores). Local: causal-control v2.
Conventions: this file = findings/why; `commit-log.md` = ledger;
`research_board.md` = status/decisions.

---

## Campaign #1 — pure-evo, LUT + egene  (`Scans/2026-05-17_pure_evo_joint`)

**Setup.** 300 random configs, productive box with `m_scale` and
`mu_lut` brackets pushed upward, N=256, 8000 ticks, neutral shadows
on. Ranked post-hoc by the *corrected* turnover-invariant composite:
`excess_pc_slope` (component-normalised whole-LUT-genome, §2),
`eg_excess_pc_slope` (fixed-space egene shadow, S2d),
`dyn_excess_pc_slope` (fixed-space LUT-entry shadow, S2d),
`n_distinct_genomes_temporal_std`, `unique_top_genomes`. 284/300
candidates (16 extinct, 0 saturated).

**Headline finding (significant).** Across **every** top config:
- `excess_pc_slope` ≈ **0.000–0.001** — the component-normalised
  *whole-genome-hash* excess shows essentially **no signal** anywhere.
- `dyn_excess_pc_slope` strongly **positive** (~23–60) and
  `eg_excess_pc_slope` strongly **positive** (~40–135).

Interpretation: at these parameters selection is real and measurable
in the **effective dynamics** (which LUT-entry transitions the alive
population exercises) and in **egene cognition**, but is *not*
detectable as whole-genome-hash excess. This directly vindicates the
§2 thesis: whole-genome FNV-hash churn ≠ effective selection;
`dyn_activity`/`eg_activity` are the sharper probes. Corollary: the
2026-04-27 scans, which optimised raw ΣG−ΣN whole-genome excess, were
chasing a quantity that under correct normalisation is ~flat — their
"`mu_lut` wants 0.06" optimum was an artdefact of the broken metric.

**New optimum (contradicts the old scans).** Under corrected metrics
the top configs favour **low `mu_lut` (0.01)** + **high `m_scale`
(2.5)** — the opposite of the old high-`mu_lut` regime. Low mutation
preserves coherent rule structure while strong eating drives the
effective-dynamics/cognition selection that the corrected metrics
actually see.

**Open design question (→ board R1).** `m_scale = 2.5` is again at
the **grid ceiling** — still unbracketed after pushing it from 1.5 to
2.5. A follow-up needs `m_scale ∈ {2.5, 3.5, 5.0}` to bracket it.

**Status.** Done (491 s on torque). `results.csv` versioned at
`e6525b5`. Feeds #2/#3 (LUT-only → egene-on-winners) and #4
(does spatial RD collapse under this objective?).
