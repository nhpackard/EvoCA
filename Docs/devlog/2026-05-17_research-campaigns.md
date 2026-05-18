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

---

## Campaign #2 — pure-evo, LUT-only (egene frozen)  (`Scans/2026-05-17_pure_evo_lut_only`)

**Setup.** 180 configs, same grid minus egene axes; `mu_egene =
mu_egenome = tax_per_egene = 0`. 258 s on torque.

**Findings.**
- **0/180 extinct** — with the fully-specified egene frozen, the
  population is extremely robust across the whole box.
- `dyn_excess_pc_slope` is **positive but markedly smaller** than #1
  joint (top configs ~3–20 vs #1's ~23–60). Freezing egene
  *reduces* effective-rule-transition excess. Cross-campaign causal
  inference (the causal-control logic at scan scale): **egene
  co-evolution amplifies LUT effective-dynamics selection** — the two
  subsystems are not independent; cognition evolving alongside the
  rule makes the rule's effective selection stronger.
- `excess_pc_slope` again ≈ 0 (occasionally tiny +ve) — whole-genome
  excess remains a non-signal, consistent with #1.
- **Optimum (contrasts #1):** LUT-only favours **high `mu_lut`
  (0.15, ceiling)** + **high `m_scale` (2.0–2.5, ceiling)** + **low
  `food_inc` (0.010, low edge)**. The mu_lut optimum is therefore
  *conditional on egene*: with egene frozen, high LUT mutation is
  favoured (raw rule exploration must substitute for the adaptive
  cognition that's switched off); with egene co-evolving (#1), low
  mu_lut wins and cognition does the adaptive work. This
  mu_lut↔egene interaction is exactly the kind of metaparameter
  relationship the program is after.
- Spatial: top LUT-only configs keep `correlation_length` ~55–75 —
  RD-scale structure persists under LUT-only pure-evo (relevant to
  #4).
- Multiple axes unbracketed at edges (`m_scale`, `mu_lut`,
  `food_inc`-low) — feeds R1 and suggests a low-`food_inc` extension
  too.

**Status.** Done. `results.csv` versioned. Winners feed #3.

---

## Campaign #8 (QUEUED, fires after #3) — direct lineage co-evolution metric

Goal: replace the *indirect* excess-comparison evidence for LUT↔egene
co-evolution with a *direct* measurement on the S2a lineage field.

Design (ready to execute):
- Run lineage ON at a #3a winner regime (LUT-pure-optimal substrate
  with egene co-evolving) + the JOINT productive corner, multi-seed,
  post-transient window.
- Per parent→child edge, decompose change into ΔLUT (LUT-only hash
  changed?) and Δegene (active egene key changed?) — snapshot per-cell
  (birth_id, parent_id, lut-hash, egene-key) periodically; rebuild
  chains via parent_id.
- Statistic: correlation / excess joint-retention of (ΔLUT, Δegene)
  across *surviving* edges **vs the neutral product of the two
  marginal mutation frequencies**. Mutations are independent by
  construction at the mutation step, so any correlation in surviving
  edges beyond the mutation-rate product is *selection jointly
  retaining co-adapted LUT+egene* — the clean, direct co-evolution
  signal.
- Decisive cross-check with #3: if 3a (good substrate) shows positive
  joint-retention while 3b (frozen LUT) cannot, synergy is real
  co-evolution, not the GoL-frozen-substrate handicap confound.

Status: queued; auto-launches once #3 results are in.

---

## Causal-control v2 (local)  (`Scans/2026-05-17_causal_control_v2`)

Bug-corrected rerun of the 2026-05-16 v1 control (which had exposed
the shadow-scope bug, EGENE_only excess ≈ −26). 3 arms × 4 seeds,
N=256, 15000 ticks. Means over seeds (cog_spec, intake, n_distinct,
**excess_pc**, **eg_excess_pc**, **dyn_excess_pc**):

| arm | cog_spec | intake | n_distinct | excess_pc | eg_excess_pc | dyn_excess_pc |
|-----|---:|---:|---:|---:|---:|---:|
| JOINT      | 11.1 | 0.044 | 16182 | 0.0002 | **64.8** | **33.0** |
| LUT_only   | 25.0 | 0.055 |  9377 | 0.0001 | **0.000** | 17.3 |
| EGENE_only | 20.4 | 0.062 |   583 | −46.4  | 31.0 | **0.000** |

**Methodological validation (the headline).** The fixed-space shadows
pass their null controls *exactly*:
- LUT_only (egene frozen) → `eg_excess_pc` = **0.000** (no egene
  evolution ⇒ zero egene excess).
- EGENE_only (LUT frozen) → `dyn_excess_pc` = **0.000** (no rule
  evolution ⇒ zero effective-dynamics excess).
The neutral-null invariant S2d was built for now holds on the *actual*
causal controls at full scale. The metric is sound.

**Bug closed.** v1 EGENE_only excess ≈ −26 (artifact). v2: the
*fixed-space* `eg_excess_pc` = **+31** (egene evolution genuinely
beats neutral), while the *whole-genome* `excess_pc` = **−46** for the
same arm — i.e. whole-genome-hash excess remains actively pathological
for egene-driven runs, exactly the §2 prediction and precisely why
the fixed-space shadow was necessary. Use `eg_/dyn_excess_pc`; never
`excess_pc` for egene-driven evolution.

**Science.** Cognition is individually selectable (EGENE_only
`eg_excess_pc` = +31, LUT frozen). Effective rule dynamics are
individually selectable (LUT_only `dyn_excess_pc` = +17, egene
frozen). In JOINT *both* are larger than their isolated values
(eg 65 > 31; dyn 33 > 17) ⇒ **LUT and egene co-evolution mutually
amplify selection in both subsystems** — synergistic, not
independent. This corroborates the scan-scale #1-vs-#2 contrast with
a clean controlled experiment (two independent methods, same
conclusion).
