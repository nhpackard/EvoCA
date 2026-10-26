# Lineage-co-evolution: guide to decisions [LIN-D1](../ISSUES.md), [LIN-D2](../ISSUES.md), [LIN-D3](../ISSUES.md) (2026-10-15)

For N. This is the homework behind the three Lineage decisions on the desk
(`06-seeding-findings.md` §D item 4). It unpacks each one, adds short
verification runs, and gives an ordered reading guide to the lab notebook and
the related records. Nothing in the repo was changed except this file. Scratch
runs live in `/private/tmp/claude-505/lin/`.

**Evidence labels used throughout.**
- **[V]**: verified here, by a headless run or by reading the code.
- **[R]**: reported by the seeding agent, the notebook, or an earlier record,
  and not re-checked here.
- **[A]**: arithmetic or an estimate from the code, not measured.

---

## 0. One-screen summary

**The three decisions, as I would now answer them**

| Decision | What `ISSUES.md` recommends | My recommendation | What changed |
|---|---|---|---|
| **[LIN-D1](../ISSUES.md)**: redesign of #8 (does selection jointly retain LUT and egene changes?) | Option A (Hamming magnitudes) with `restricted_mu=1`, plus a synthetic-coupling power check | **Keep #8, but use Option A′.** (i) Measure the LUT magnitude **post hoc on expressed bits**: Hamming distance restricted to the lineage's own used-mask (`get_eff_used_mask`). (ii) Leave `restricted_mu=0`. (iii) Fit an outcome model on **all** edges with a ΔLUT×Δegene interaction, not a correlation among survivors. (iv) Use the child's **offspring count** as the outcome. (v) Run the power check post hoc on real edge tables. | `restricted_mu` masks with the **population union** of queried bits, which is 180–250 of 250 bits whenever the population isn't crashing [V]. It is nearly a no-op, and it makes the LUT mutation rate depend on demography. Correlation among survivors is also **scale-dependent**: additive harm with no co-adaptation gives ρ = −0.055 among survivors in a toy model [V]. |
| **[LIN-D2](../ISSUES.md)**: the lineage unit for H1 | (a) restrict to `mu_lut` ≤ 0.003, plus (b) a Python clade roll-up along `parent_id` | **Use barcoded clades, tracked exactly in Python every step.** Every alive cell gets a barcode at a marking time, and children inherit it. Re-mark every few hundred steps. Fitness is a clade's log-frequency growth against a neutral null driven by the same demography. Keep the genotype rows for the D3 punctuation work. Low `mu_lut` becomes a cross-check, not the main regime. | A vectorized per-step tracker is **exact**: 0–3 unresolved parent lookups in 114k–327k births. It costs a few ms per step at N=256 [V]. Clades coalesce fast: 9216 marks fall to 2–5 survivors in 1200 steps at N=96 [V], so re-marking is required. Genotype churn is driven by the **egene** rates as well as by `mu_lut` [A], so restricting `mu_lut` alone does not fix the unit. |
| **[LIN-D3](../ISSUES.md)**: actual-tax column in the reducer | Add the column (small C change) before H1 relies on `net` | **Defer the column. Decide a new upstream defect first: food overdraw at `m_scale` > 1.** If overdraw is fixed, nominal and actual tax differ by about 0.4% [V] and the column is optional. If overdraw is intended, define an energy-consistent net per clade instead. | **New finding [V]:** eating can take more than the food present, `mouthful = (m/25)·score·F` > F, whenever `m_scale·score/25` > 1. F_food goes negative in 22–44% of cells, live cells carry negative `f_priv`, and "eviction dissipation" goes negative (−349k in the 4b control). At `m_scale=3` (the #3c/#8/notebook regime), the 20–23% gap between nominal and actual tax is almost entirely this effect. At `m_scale=1` the gap is 0.4%. |

**New cross-cutting defect for the A-list (proposed A9) [V].**

The eating step does not clamp the mouthful to the food available. The C
comment at `C/evoca.c` ~1403 says "mouthful never goes negative"; that holds
only while F_food ≥ 0, and overdraw breaks that.

The energy identity still closes to 0.04 food units, because the debt shows up
as negative sinks. But the population is subsidized by debt that is erased when
indebted cells die. In the 4b control, the population removes 764k of tax from
417k of intake, and the difference is −349k of "eviction loss".

Affected records and plans:
- campaign #3 (`m_scale` 1.5–2.5);
- #3c and the #8 spec (3.0);
- the pure-evo winners (up to 2.5);
- Resource-driving R1 (2.5–5) and the planned [RES-N1](../../Resource-driving/ISSUES.md) (`m_scale` ∈ {2.5, 3.5}).

This is N's call. A fix belongs on its own branch with the parity probe at
`m_scale` ≤ 1, where it should be a no-op.

**Read first** (about 40 minutes; §7 has the full guide):
1. Notebook §4b, **cells 30–32**: the nominal/actual tax table and the
   "−0.020 per-capita net" observation. These are the D3 numbers. Read them
   alongside §2 of this guide, which explains them as overdraw.
2. Notebook §4a, **cells 26–28**: the genotype-singleton table. Note that the
   `mu_lut` 0.01/0.03 columns are crashed populations of 320–410 cells (D2).
3. Notebook §3, **cells 20–23**: the binary-indicator saturation and Hamming
   variance (D1).
4. `Research/Lineage-co-evolution/plans/plan-8-lineage-coevolution.md`, "Possible fixes", **Option F
   caveat**. It already warned that the active footprint is most of 250 bits;
   `ISSUES.md` dropped that warning (D1).
5. `Docs/checkpoint 3 July 2026.md`, lines ~1376–1435. The original D2 intent
   was to "reduce by lineage root (walk `parent_id` chains)", that is, clades
   rather than genotypes. The H1 program, step 2, asks for a per-lineage shadow
   that does not exist (D2, H1).

---

## 1. What the workstream is trying to establish, and with which instruments

### 1.1 The two scientific questions

- **Part 1 (coupling):** do LUT and egene changes co-adapt? Selection would
  retain *combinations* beyond what retention of each separately predicts.
  Only #8 addresses this directly.
- **Part 2 (H1, adaptation vs drift):** can lineage-level quantities separate
  adaptive lineages from drifting ones? This comes from checkpoint §H,
  expanded at `Docs/checkpoint 3 July 2026.md` ~1407. Its four steps are:
  1. instrument: lineage field plus the D2 reducer;
  2. per-lineage neutral control: "does a lineage's realised growth/net-return
     exceed what its *shadow* (same birth/death counts, non-selective
     inheritance) produces?";
  3. two cross-checked readouts, A_new and lineage net-return staircases;
  4. the decisive §B `food_inc/tax` sweep, with the falsifier that raw activity
     should *not* be inverted-U.

H1 is not a single hypothesis. It is a program whose central claim is:
*adaptive fraction and shadow-beating lineage count are inverted-U in
`food_inc/tax`, while raw activity is not.*

### 1.2 What each instrument measures (from the code) [V]

| Quantity | Where | What it is exactly |
|---|---|---|
| `lin_birth_id[i]` | `C/evoca.c` ~1835 | Monotonic id assigned at each Phase-4 birth into slot i. It is 0 for founders, meaning cells alive when lineage was enabled. |
| `lin_parent_id[i]` | same | The parent's `lin_birth_id` at reproduction (0 if the parent was a founder). Only the current occupant's link is retained; **dead cells' links are lost** unless sampled. |
| `lin_parent_hash[i]` | same | The parent's **full-genome** FNV-1a hash (`lut_hash_cache`, which is LUT bytes ‖ sorted *active* egene value/mask bytes). |
| `lin_intake_cell`, `lin_tax_cell` | ~1588, ~1636 | Per-cell accumulators over the **current life**, reset at birth. Intake is the realized mouthful after the headroom clamp, **including negative mouthfuls** (§2.1). Tax is the **nominal** `t` charged each step, including the step on which the cell dies. |
| Reducer row (`get_lineage_return`) | `lret_*` ~260–310, ~1027–1090 | Keyed by **full-genome hash**, so a row is a genotype, not a clade. Fields: `intake`, `tax` (= dead + live members), `net = intake − tax`, `births` (births *into* this genotype, counted at birth against the post-mutation hash), `pop`, `first_step`, and `last_step` (only as current as the last `lineage_return_update()`). Never compacted. |
| What `net` excludes | — | (i) **The birth endowment**: the parent's half of `f_priv` passes to the child, so a mutant genotype starts with 0.5 food that never appears as intake. (ii) Food destroyed on eviction. For genotype rows, the endowment crosses between rows (parent genotype → mutant child genotype), so short-lived genotypes that live off their endowment show `net` < 0 by construction. In 4b, 70% of evolving lineages with ≥ 3 births have `net` < 0 [R]. For a clade, endowments are internal except at the founding mark. |
| `diss_tax_cum` (`get_dissipation`) | ~1590–1607 | Tax **actually removed**: the full `t` for survivors, and the pre-tax `f_before` for a tax death. `f_before` can be **negative** (§2.1). |
| `diss_death_cum` | ~1698 | `f_priv` of an evicted cell. It is negative when the evicted cell is in debt (§2.1). |
| `restricted_mu` mask | ~1480–1511 | `lut_active` = OR over **all alive cells** of the LUT bit each cell queried **this step**. Mutation draws from the `n_active` set bits. It is **not** per-cell. |
| `eff_used_mask` (per cell) | ~1490, ~1716 | The bits a cell has queried, **inherited from the parent at birth** and accumulated after that. It is a lineage-accumulated expressed-bit set. This is what Option A′ uses. |

### 1.3 How the mutation operator bears on #8 [V]

LUT flips, egene active-bit flips, and egene value/mask flips are separate
Poisson draws (`C/evoca.c` ~1727–1800). So ΔLUT and Δegene **are**
independent at birth, and the premise of #8 holds.

The egene value/mask pass walks all 8 slots × 12 bits, including inactive
slots. Only flips that land in **active** slots change the genome hash.

---

## 2. Verified findings that change the decision menu

### 2.1 Food overdraw at `m_scale` > 1 (bears on D3, H1, and the regime choice for D1)

**Mechanism [V].** In Phase 3:

```
mouthful = (m_scale/25) * eat_score * F_food[idx]   // capped only by headroom 1 - f_priv
F_food[idx] -= mouthful
```

There is no `min(mouthful, F_food)`. When `m_scale·eat_score > 25`, a cell
takes more than the food present and F_food goes negative. The next eater on
that site receives a **negative** mouthful, so its `f_priv` falls and can go
below 0. Live cells with `f_priv` down to −1.29 were observed. A cell in debt
then dies, by tax (actual tax removed = negative `f_before`) or by eviction
(negative death sink), and the debt disappears. That is net energy creation.

**Magnitude (N=96, seed 2, #3c regime except `m_scale`, `mu_lut`=0.0003,
`mu_egene`=`mu_egenome`=0.003, 2000 steps) [V]:**

| `m_scale` | min F_food | mean fraction of cells with F < 0 | nominal / actual tax | `death_cum` | final pop |
|---|---|---|---|---|---|
| 1.0 | 0.000 | 0.000 | 1.004 | +8,239 | 355 |
| 1.25 | −0.107 | 0.224 | 1.006 | +41,942 | 411 |
| 1.5 | −0.163 | 0.249 | 1.005 | +40,655 | 506 |
| 2.0 | −0.243 | 0.264 | 1.006 | +47,126 | 5,100 |
| 2.5 | −0.598 | 0.364 | 1.165 | −76,021 | 6,702 |
| 3.0 | −0.663 | 0.443 | 1.200 | −209,861 | 8,982 |

**The 4b runs reproduced exactly (N=96, 3000 steps) [V].**
- Evolving run: intake 414,067; nominal tax 887,729; actual tax 738,159;
  `death_cum` **−325,345**.
- Control: intake 416,977; nominal tax 939,192; actual tax 763,713;
  `death_cum` **−348,943**.
- The whole-system identity closes to |residual| < 0.05 in both, so the
  accountant is self-consistent and is faithfully booking the debt.
- The same control at `m_scale=1.0` gives nominal 349,367 vs actual 347,794
  (0.45%), positive `death_cum`, and pop 331.

**Consequences.**
- **[LIN-D3](../ISSUES.md)'s 20–23% gap is overdraw, not "tax a dying cell could not pay".**
  An actual-tax column would faithfully record negative tax for indebted
  cells; it would encode the defect rather than correct for it.
- **[LIN-Q3](../ISSUES.md)'s "constant per-capita net −0.020"** (control, cell 32) is the
  leak's signature. A population of 9,136 paying 0.035 each per step on an
  intake of about 139 per step can only persist on erased debt. The structural
  point of [LIN-Q3](../ISSUES.md) still holds without the leak: intake ≈ source, so raw net is
  a pool share. But the negative sign and the high density are artifacts.
- **Density at `m_scale=3` (0.97 at N=96) far exceeds `food_inc/tax` = 0.43.**
  The Resource-driving pilot reports density ≈ `food_inc/tax` below 1 [R].
  The #3c and #8 regime is therefore partly leak-supported. Every #3c result
  stands as a record of *this* model, but read "viable dense regime" there
  with this in mind.
- The flux merge test (`tests/notebooks/metric-flux-test.ipynb`) ran at
  `m_scale` 2.0, where `death_cum` is still positive, so it could not have
  caught this [V: grep of the notebook].

**Fix options (N's decision, not done here).**
- (a) **Clamp:** `mouthful = min(mouthful, max(F_food, 0))`. This is
  dynamics-altering for every run in which `m_scale·score > 25` ever occurs,
  and a no-op otherwise, so the parity probe at `m_scale` ≤ 1 should be
  bit-identical. It needs its own branch, a pre-merge notebook, and a review
  in the `/tmp` worktree.
- (b) **Keep overdraw as intended semantics** ("food debt") and document it.
  Then every energy-based lineage quantity must include eviction loss and use
  actual tax.

I recommend (a). Option (b) makes `net` uninterpretable as "metabolic
success", because a lineage can profit by dying in debt.

### 2.2 `restricted_mu=1` is close to a no-op and couples mutation to demography (D1) [V]

`n_active` is the per-step **population union** of queried bits (§1.2).
Measured from `v`/`alive` every 100 steps:

| Run | n_active (of 250) |
|---|---|
| N=256, `mu_lut`=0.03 | 225–250 |
| N=96, `mu_lut`=0.003 | 110–236 |
| N=96, `mu_lut`=0.03 | 9–242. It tracks a boom-bust cycle (pop 207 ↔ 8,066, period ~700 steps): about 240 at booms, 9–42 at crashes. |

So with `restricted_mu=1`:
- The Poisson rate falls by only 0–30% in normal phases.
- Flips are **not** made phenotypic for the child. Each cell queries one bit
  per step, and the union is everyone's bits.
- During crashes, all LUT mutation is concentrated on about 10 bits.

The LUT mutation spectrum then depends on demographic phase, which is a new
confound for any survival-based statistic. It also makes #8 study a different
evolutionary process from #3c and H1.

`Research/Lineage-co-evolution/plans/plan-8-lineage-coevolution.md` already noted under Option F that "active
footprint at N=256 is often most of 250 → may not fully de-saturate";
`ISSUES.md` lost that caveat. The binary indicator would still saturate:
`1 − exp(−0.03·240)` = 0.9993.

**Better route (post hoc, no dynamics change).** Each cell's lineage-inherited
`eff_used_mask` has a median of 60–130 set bits by step ~1500 in these runs
[V]. "Expressed Hamming" is

```
popcount((L_parent XOR L_child) & M_child_at_first_observation)
```

It counts flips on entries the lineage actually exercises. It costs one more
N²×32-byte read per snapshot, the same as the LUT bytes. It does not change
the model, so #8 stays comparable with #3c and H1.

### 2.3 Correlation among survivors is a scale-dependent test, not a co-adaptation test (D1) [V]

Given independent ΔL and ΔE at birth, the joint distribution among survivors
is ∝ p(ΔL)·p(ΔE)·P(survive | ΔL, ΔE). The two stay independent **iff**
survival factorizes multiplicatively.

So a nonzero survivor correlation means "non-multiplicative survival", and
that includes plain truncation or probability-scale additivity with **no**
co-adaptation (Berkson's paradox). Toy model, n=20,000 edges,
ΔL ~ Poisson(7.5), ΔE ~ Poisson(1):

| True survival model (no co-adaptation in any row) | Spearman among survivors |
|---|---|
| multiplicative, log-additive | −0.010 (null) |
| additive on probability scale | **−0.055** (about 4.5σ at n≈6,600) |
| truncation on a linear index | **−0.030** |

**Implications.**
- (i) Negative association is uninterpretable as biology.
- (ii) *Positive* (synergistic or compensatory) interaction, consistent on both
  the log and logit scales, is the defensible co-adaptation signature.
- (iii) Fit outcome ~ ΔL + ΔE + ΔL·ΔE on **all** captured edges, not only on
  survivors. That also removes the arbitrary `T_surv` cut.

The original plan's statement that the within-survivor permutation null
"cancels survivor-capture bias" is true for *capture*, but not for this
selection-induced dependence.

### 2.4 Genotype-keyed "lineages": A6 is certain, and the reason is broader than `mu_lut` (D2) [V]/[A]

**[V] A6 holds across seeds.** N=96, 2000 steps, `mu_egene`=`mu_egenome`=0.01,
three seeds each:

| `mu_lut` | lineages with ≥ 10 births | max births in one genotype | fraction of cells sharing a genotype | min / final pop |
|---|---|---|---|---|
| 0.003 | 158–181 | 2,223–2,313 | 0.29–0.34 | ~2,400–2,900 / 3,000–4,700 |
| 0.01 | 8–10 | 323–376 | 0.01–0.06 | **326–410** / 410–7,900 |
| 0.03 | **0** | 2–6 | 0.000 | **207–222** / 212–514 |
| 0.03, egene rates 0 | **0** | 2–4 | 0.000–0.001 | 367–464 / 564–3,131 |

**[A] The singleton result is arithmetic, not a finding about the runs.** A
genotype's expected number of faithful copies per birth is set by the total
per-birth genome mutation rate:

```
U = 250·mu_lut + 96·mu_egene·(Negene/8) + 8·mu_egenome
```

At `mu_lut`=0.03 the LUT term alone gives P(child = parent genotype) ≈ e^−7.5
= 0.0006. The egene term matters too: at `mu_egene`=`mu_egenome`=0.01 and
Negene 1–3 it is 0.2–0.44, so 18–36% of births mint a new genotype even with
the LUT frozen. The 4b "slow" evolving run (`mu_lut`=0.0003, egene rates
0.003) still produced 1.14 M genotypes in 3000 steps [R]. Restricting
`mu_lut` (option (a) in `ISSUES.md`) does not make genotypes persist unless
the egene rates come down as well.

**[V] The `mu_lut` ≥ 0.01 rows in notebook cell 26 are crashed or bottlenecked
populations** (320–410 alive of 9,216 at t=2000 in seed 2). "Productive
corner" is a stretch for those exact runs. At N=256 the same `mu_lut`=0.03
start keeps 10k–43k cells over 600 steps [V, one seed], so the bottleneck is
at least partly an N=96 effect.

### 2.5 Exact clade (barcode) tracking is cheap in pure Python (D2) [V]

**Method.**
1. At a marking step, give every alive cell a barcode, its lattice index.
2. Every step, newborns are cells whose `birth_id` changed.
3. Look up each newborn's `parent_id` in the *previous* step's
   `birth_id → barcode` table (a sorted array with `searchsorted`, fully
   vectorized). The parent cannot be born in the same step, because a newborn
   holds at most 0.5 food. An evicted-after-reproducing parent is still in
   the previous step's table.

**Results.**

| Run | Unresolved lookups | Step + tracking time | Surviving clades |
|---|---|---|---|
| N=96, `mu_lut`=0.03, 1,200 steps | 3 in 114,410 births | 0.8 ms per step | 5 of 9,216 marks |
| N=96, `mu_lut`=0.003 | 0 in 197,709 births | 1.0 ms per step | 2 |
| N=256, `mu_lut`=0.03, 500 steps | 2 in 327,211 births | 5.7 ms per step | 69 of 65,536 |

The unresolved lookups are probably children of founders (`parent_id`=0).
Script: `/private/tmp/claude-505/lin/barcode.py`.

**Implications.**
- (i) [LIN-D2](../ISSUES.md) option (b) is **exact and cheap if sampled every step**. With
  periodic snapshots it is lossy, because dead intermediates' links are lost.
  The #8 spec's S=50 has the same loss.
- (ii) Coalescence is fast (boom-bust bottlenecks), so H1 needs **re-marking**
  every ~100–300 steps and a neutral coalescence reference under the realized
  demography.
- (iii) The same per-step tracker yields each birth's parent exactly. That
  gives #8 an **offspring-count outcome** per edge and an exact edge table with
  no `T_surv`/S capture bias. This is one shared tool for D1 and D2.

### 2.6 The #3c "C beats A" result is phase-1-`mu_lut`-specific and demographically confounded (context for D1) [V]

This is from `Scans/2026-05-19_3c_coevolution_substrate/results.csv`,
re-analyzed. Paired C−A `eg_excess_pc_slope`:

| Phase-1 `mu_lut` | Mean C−A | C > A |
|---|---|---|
| 0.01 | +0.03 | 4/10 |
| 0.03 | +1.02 | 9/10 |
| 0.06 | +0.80 | 9/10 |

Overall: 22/30, sign test p=0.016, Wilcoxon p=0.0009.

The arms also differ hugely in demography. Mean final pop is A 48.9k vs C
27.1k; mean min pop is A 21.7k vs C 5.6k. C cycles boom-bust, and the
cumulative excess moves with the cycles (notebook cell 18). The modest
"co-evolving beats freezing" result is therefore real as a record. It is
not a coupling measurement, and that motivates keeping a direct #8 test.

### 2.7 Scaling hazards in the genotype reducer (D2, H1 campaign) [A]

- **Table size.** An `lret` entry is about 56 B plus a 4 B key, at load
  factor ≤ 0.7, and the table is never compacted. The 4b evolving run's 1.14 M
  genotypes (N=96, 3000 steps) is about 2^21 slots, roughly 126 MB. An H1 run
  at N=256 for 8000 steps scales about ×19: roughly 22 M genotypes and about
  2 GB per run.
- **Hash collisions.** Keys are 32-bit, so expected colliding pairs ≈ n²/2³³:
  about 150 at 1.1 M genotypes and about 56k at 22 M. That is harmless for
  totals but corrupts parent-hash chains.

A clade-keyed reducer is bounded by N² barcodes per mark and has neither
problem.

---

## 3. [LIN-D1](../ISSUES.md): the redesign of #8

### 3.1 The hypothesis #8 tests

LUT and egene changes arise independently at birth, so if selection is blind
to their *combination*, any association among retained changes is fully
explained by the separate effects. A direct test asks whether
the fate of an edge depends on the **interaction** of ΔLUT and Δegene beyond
separate effects. It needs no shadow and no frozen arm.

### 3.2 Why the current design fails

1. **Saturated binary indicator.** P(ΔLUT>0) = 1−exp(−250·`mu_lut`): 0.918 at
   0.01 and 0.9995 at 0.03. This was measured to three digits in notebook cell
   18 [R: notebook output; arithmetic V]. The 2×2 table has no variance.
2. **Silent flips dilute any magnitude.** Most flipped bits are never queried
   by the lineage (§2.2).
3. **Survival over 40–200 steps is a weak proxy, and the survivor-only
   statistic is scale-dependent** (§2.3). Notebook cell 23 already notes that
   the survival filter "barely moves P(ΔLUT>0)" [R].
4. **Snapshot capture (S=50 in the spec) drops short-lived edges.** The
   permutation null does not repair this, because it permutes within what was
   captured.
5. **The regime is the overdraw regime** (`m_scale=3`, §2.1).

### 3.3 Options

| Option | What it is | Pros | Cons | Build cost |
|---|---|---|---|---|
| Retire #8 | Rely on H1 | Free | Part 1 (coupling) stays untested; H1 cannot answer it | 0 |
| **A (ISSUES)** | Hamming magnitudes, Spearman among survivors, `restricted_mu=1` | Has variance; harness exists | `restricted_mu` is nearly a no-op and demography-coupled (§2.2); scale-dependent statistic (§2.3); silent-bit dilution persists | ~1 day Python |
| **A′ (recommended)** | Expressed Hamming (`eff_used_mask`, post hoc) for ΔLUT; active-slot Hamming for Δegene; exact per-step edge table (§2.5); outcome = child's **offspring count** (plus lifespan as secondary); a GLM with interaction on **all** edges, fitted on log (Poisson) and logit scales; a pre-registered **positive** interaction stable on both scales; `restricted_mu=0` | No dynamics change; comparable with #3c and H1; no `T_surv`/S artifacts; correct estimand | More analysis code; expressed-bit masks are lineage-accumulated, so they are an approximation of "expressed for this child" | ~2–3 days Python (the tracker is shared with D2) |
| B | Change in the `lut_complexity` level | Phenotypic; links to the Ring workstream | Low variance; needs a per-cell accessor (small C) | small C plus Python |
| C | Mutual information on binned magnitudes | Catches non-monotone coupling | Weaker power; binning choices | Python |
| D | Joint persistence across generations along chains | Temporal retention | Needs exact genealogy (now cheap, §2.5) and longer runs | Python, after A′ |
| E | Conditional per-bit retention rate | Fitness-cost reading | Baseline subtleties | Python |
| F via `restricted_mu` | — | — | Superseded by post-hoc expressed bits (§2.2) | — |

**The power check, made concrete (free, post hoc).** Take a real edge table
(ΔL_expr, ΔE, offspring). Replace the outcome with a synthetic one drawn from
a model with known interaction β_LE, on the log and logit scales in turn.
Recompute the statistic and plot power against β_LE. Also run the
*additive-probability, no-interaction* synthetic outcome, to measure the
false-positive rate of the scale artifact.

This calibrates detectability. It cannot say whether a biological effect that
large exists.

### 3.4 Unresolved for D1

- **Regime.** Run #8 at `m_scale=3` for continuity with #3c, or in a
  non-overdraw regime? This depends on the A9 decision. Suggestion: run a
  pilot in both.
- **Δegene measure.** Byte Hamming over active slots (current) versus a
  functional measure, such as the change in best-match score against the
  local fiducial pattern. Hamming is the first pass.
- **Is the offspring count sufficient at N=256?** The power check answers
  this.

---

## 4. [LIN-D2](../ISSUES.md): the lineage unit for H1

### 4.1 The problem

D2 keys on the full-genome hash, so a "lineage" is a genotype. At any
`mu_lut` where the LUT evolves, and at moderate egene rates, genotypes are
near-singletons (§2.4).

Per-genotype `net` is then a lifetime quantity of individuals, biased negative
by the uncounted endowment (§1.2). The original D2 intent was a clade unit:
the checkpoint says "reduce by lineage root (walk `parent_id` chains)"
(`Docs/checkpoint 3 July 2026.md` ~1380).

### 4.2 Options

| Option | Pros | Cons | Cost |
|---|---|---|---|
| (a) Low `mu_lut` (≤ 0.003) | No dev; genotypes persist | Excludes the regime of #3c/#8; **also needs low egene rates** (§2.4); bottleneck dynamics differ | 0 |
| (b) Python roll-up from periodic snapshots | No C | Lossy: chains break at unsampled dead intermediates; hash↔clade mapping is ambiguous | moderate |
| **(b+) Per-step barcode tracker (recommended)** | **Exact** (§2.5); cheap; any mutation rate; a clade's log-frequency growth is the relative-fitness quantity [LIN-Q3](../ISSUES.md) asks for, immune to pool-share pinning; endowments internal; table bounded by N² per mark | Clades coalesce fast, so re-marking is needed; per-clade intake/tax come from per-cell accumulators (`get_lineage_intake`/`_tax`) and lose about one step per death | ~1–2 days Python |
| (c) C clade reducer, keyed by an inherited per-cell barcode | Fast at scale; exact intake/tax | C change on merged code (additive, parity-safe); only needed if (b+) is too slow for N=256 campaigns | small C |
| (d) Full birth-event log (child, parent, step, hash) | Any post-hoc unit | About 24 B per birth: ~2.4 GB per N=256/8000-step run | moderate C/IO |
| Add `parent_hash` to the reducer row | Exact genotype phylogeny; solves the [LIN-N1](../ISSUES.md) gap (i) for D3 punctuation | Does not fix the H1 unit on its own | ~10 lines of additive C |

**The neutral null for H1, which is the real open item.** H1 step 2 assumes a
per-lineage shadow. The neutral shadow is population-level and LUT-only (see
the memory note *shadow-models-lut-only-mismatch*), and it currently perturbs
the RNG stream (A2). Candidate nulls for clade growth:

- (i) A Python neutral resampler driven by the recorded per-step birth and
  death counts: non-selective parent choice, the same demography. This is
  cheap and is literally §H step 2.
- (ii) Barcodes carried in the neutral shadow (C, after A2 merges).
- (iii) A coalescence-time comparison: real versus neutral resampler.
  Selective sweeps coalesce faster than neutral.

### 4.3 Recommendation

Use (b+) as the H1 unit, with null (i). Keep the genotype reducer for D3
punctuation work, and add `parent_hash` to its rows in the same small additive
C change as anything else that touches the reducer. Use (a) only as a sanity
cross-check, at low `mu_lut` **and** low egene rates. Decide on (c) only if
(b+) proves too slow at N=256.

---

## 5. [LIN-D3](../ISSUES.md): tax basis in the reducer

### 5.1 The problem as originally framed

The reducer books nominal `t`, while the accountant books tax actually
removed. Reported gap: 20–23% (4b), with `net` biased low and more so for
short-lived lineages [R; reproduced V].

### 5.2 What the gap actually is

The gap is overdraw (§2.1). It is 0.4–0.6% at `m_scale` ≤ 2 and 16–20% at
2.5–3.0 [V].

The original explanation, "a dying cell could not pay its full `t`", accounts
for the small residual only. So the decision is downstream of A9:

| A9 outcome | What to do about D3 |
|---|---|
| Overdraw fixed (clamp) | The gap is ~0.4%. The column is optional exactness; add it only if some other change is touching the reducer anyway, together with `parent_hash`. |
| Overdraw kept as debt semantics | Do not add a bare actual-tax column. Define the per-clade energy budget `intake − actual tax − eviction loss = Δstore`, booked per clade at the barcode level (§4). Report net only alongside that identity. |
| Undecided | Treat `net` as relative within a run only, not absolute, and keep H1 fitness on births and clade growth rather than `net` ([LIN-Q3](../ISSUES.md)). |

### 5.3 Unresolved for D3

Is overdraw intended by N? The C comment (~1403: "mouthful never goes
negative") suggests not.

---

## 6. Interaction with A2 (shadow RNG) and the `shadow-rng` branch

- **Example runs that enable the shadow will change** when A2 merges:
  notebook §1 (cells 8–11) and §2 (cells 15–18). Their conclusions
  (frozen-GoL viability is a coin flip; the A/C fork is bit-identical) should
  survive, but the specific outcomes, such as config 253 surviving or dying,
  will move.
- **Example runs without the shadow should stay bit-identical**, if the fix
  only re-routes the shadow's draws: §3 (cell 20 `edge_table`) and §4
  (cells 26 and 30 `run_lineage`) never call `neutral_enable()` [V: code
  read]. That makes them a convenient parity check on the `shadow-rng` branch.
- **[LIN-Q4](../ISSUES.md) is answered** by A2 [R/V in `06-seeding-findings.md`]. Enabling
  the shadow does perturb the real run's stream. Close it once A2 merges.
- **The #3c fork design stays valid**, because both arms ran shadow-on
  identically.
- **#8 is shadow-free**, so A2 does not affect it.
- **H1 benefits most**: after the fix, shadow-on and shadow-off runs are
  paired, so neutral-null comparisons can run in the same trajectory.

---

## 7. Reading guide (ordered)

Cell indices are 0-based positions in
`Research/Lineage-co-evolution/Lineage-co-evolution.ipynb`. Headings are given
so you can find them without counting.

| # | Where | What to look at | Bears on |
|---|---|---|---|
| 1 | §4b "Per-lineage net-return trajectories…", **cell 30** (output table) and **cell 32** (observations) | `total intake` ≈ source (414,067 vs 414,720); nominal vs actual tax; the "−0.020 per-capita net" bullet. Read it knowing the dense, negative-net steady state is the overdraw leak (§2.1). The 70% `net` < 0 also reflects the uncounted endowment (§1.2). | D3, [LIN-Q3](../ISSUES.md), A9 |
| 2 | §4a "What does a hash-keyed 'lineage' resolve…", **cell 26** (table) and **cell 28** | The "alive cells" row: 6,385 → 3,160 → **410 → 320**. The `mu_lut` ≥ 0.01 columns are bottlenecked populations. 582k genotypes even at `mu_lut`=0.001 shows the egene-driven churn (§2.4). | D2 |
| 3 | §3 "Campaign #8 … binary-indicator degeneracy", **cell 20** (table), **cell 22** (figure), **cell 23** | Analytic vs measured P(ΔLUT>0); the Hamming sd growing 0.5→2.7; Spearman ≈ 0 with ~5–11k surviving edges (noise floor about ±0.01). Notice that "surviving" means seen across ≥ 40 steps; that is the cut A′ removes. | D1 |
| 4 | §2 "Campaign #3c…", **cell 15** (table), **cell 17** (figure), **cell 18** | C's boom-bust cycles in seed 1; the dyn-excess artifact in A; paired C−A of −3.9/−1.0/−0.3/+26 across seeds. Then see §2.6 here for the by-`mu_lut` split of the recorded campaign. | D1 motivation, [LIN-Q2](../ISSUES.md) |
| 5 | §1 "Campaign #3…", **cell 8** (table) and **cell 11** | Config 253 is extinct here but survived in the record; the shadow flip; 3b at 46% extinct. A quick read. | A2, [LIN-Q4](../ISSUES.md) |
| 6 | §4c "D3 preview…", **cells 34–35** | Negene saturates by ~500 steps; no population-level staircase. The per-lineage linkage need is what `parent_hash` in the reducer would solve. | D2 (genotype rows retained), [LIN-N1](../ISSUES.md) |
| 7 | §5 "Current state and directions", **cell 36** | Note that it lists only [LIN-D1](../ISSUES.md) and [LIN-D2](../ISSUES.md); **[LIN-D3](../ISSUES.md) is missing** there (a minor inconsistency with `ISSUES.md`). | — |
| 8 | `Research/Lineage-co-evolution/plans/plan-8-lineage-coevolution.md`: "Possible fixes" (Options A–F, especially the **F caveat**), then "Original spec" D3/D5 and "Caveats" | The spec's claim that the permutation null "cancels survivor-capture bias" (true for capture, not for the §2.3 selection effect); the regime `m_scale=3.0`, S=50, `T_surv`=200. | D1 |
| 9 | `Docs/checkpoint 3 July 2026.md` ~1376–1381 (original D2 intent: reduce by lineage root) and ~1407–1435 (H1 four-step program) | What H1 promises; step 2's per-lineage shadow, which does not exist. | D2, H1 |
| 10 | `Docs/Archive/Oct2026/01-branch-review.md`, "Semantic note for N" and "The `lineage-return` crash" | The nominal/actual definitions. The `enable_lineage()` trap: `init(lineage=True)` is the safe path. | D3 |
| 11 | `Research/Lineage-co-evolution/Scans/2026-05-19_3c_coevolution_substrate/notes.md` | The recorded #3c finding as written in May. | D1 context |
| 12 | `Research/Lineage-co-evolution/digests/coevolution-substrate.ipynb`, cell 0 | Draft status line is stale ("Research/_drafts/"); its outputs were run from `/private/tmp/evoca_restructure`. It needs the §2.6 split before it is finalized ([LIN-N6](../ISSUES.md)). | housekeeping |

To reproduce the §2.6 split (5 lines):

```python
c3 = pd.read_csv('Research/Lineage-co-evolution/Scans/2026-05-19_3c_coevolution_substrate/results.csv')
a = c3[c3.arm=='A'].set_index(['param_mu_lut_p1','seed']); c = c3[c3.arm=='C'].set_index(['param_mu_lut_p1','seed'])
d = c.eg_excess_pc_slope - a.eg_excess_pc_slope
print(d.groupby(level=0).agg(['mean', lambda s: (s>0).sum()]))
print(a.final_pop.mean(), c.final_pop.mean(), a.min_pop.mean(), c.min_pop.mean())
```

---

## 8. Other results worth examining

- `tests/notebooks/lineage-return-test.ipynb`: the reducer's merge test,
  run at `m_scale=0.5`, so overdraw cannot appear there.
- `tests/notebooks/metric-flux-test.ipynb`: the energy identity test at
  `m_scale` 2.0 (`death_cum` +121,640). This is the edge of the overdraw
  regime, so it is not a counterexample.
- `Research/Lineage-co-evolution/Scans/2026-05-17_pure_evo_egene_on_winners/results.csv`:
  campaign #3, `m_scale` 1.5–2.5 (median 2.5), so partly in the overdraw
  regime. It has no `notes.md` or `scan_config.json` ([LIN-Q6](../ISSUES.md)).
- `Research/Lineage-co-evolution/Scans/2026-05-19_8_lineage_coevolution/scan.py`:
  the harness. Its reusable machinery is the incremental `birth_id` map and
  the permutation null. `_egene_keys` is the D1 egene key.
- `Docs/devlog/2026-05-17_research-campaigns.md` ~185: the #8 rationale as
  first written, including the "vs mutation-rate product" framing.
- `Research/Resource-driving/ISSUES.md` ([RES-Q5](../../Resource-driving/ISSUES.md), [RES-N1](../../Resource-driving/ISSUES.md)) and
  `Research/Resource-driving/Scans/2026-05-18_R1_mscale_bracket/`: the
  bracket covers `m_scale` 2.5–5, and the planned [RES-N1](../../Resource-driving/ISSUES.md) uses {2.5, 3.5}.
  Both are in the overdraw regime. The "U-shaped viability in `m_scale`" may
  partly be the leak. This is cross-workstream; flagging it only.
- Memory notes: `dyn-shadow-frozen-rich-substrate` (do not read
  `dyn_excess_pc` in `mu_lut`=0 arms) and `shadow-models-lut-only-mismatch`
  (excess is undefined for egene-driven runs). Both constrain any H1 use of
  the shadow.

---

## 9. Concrete next steps

**A9 (prerequisite; N decides).**
1. Decide whether overdraw is intended.
2. If not, open a branch with the clamp. Run `tools/parity_probe.py` at
   `m_scale` ≤ 1; it should be IDENTICAL. Add a regression test that enters
   through `run_sim`: at `m_scale=3`, assert `min(F_food) ≥ 0` and
   `death_cum ≥ 0`, and show that it fails without the fix.
3. Add a pre-merge notebook re-running notebook §4b at `m_scale=3`.
4. Add an A9 row to `06-seeding-findings.md` and notify the Resource-driving
   workstream.

**[LIN-D1](../ISSUES.md).**
1. Approve A′ (§3.3), with `restricted_mu=0`.
2. Build the per-step tracker (shared with D2), the expressed-Hamming
   computation, and the GLM with interaction on log and logit scales.
3. Run the post-hoc power check on N=96 edge tables at two regimes (one with
   `m_scale`=3, one without overdraw).
4. Only then run the 2 arms × 8 seeds of the original spec. The frozen arm
   becomes optional, since ΔL≡0 is a trivial null.

**[LIN-D2](../ISSUES.md).**
1. Approve (b+) barcodes as the H1 unit, with re-marking every ~200 steps and
   the neutral-resampler null.
2. Pilot at N=96: the clade log-frequency growth distribution, real vs
   resampler, plus coalescence times.
3. Decide on the (c) C reducer only if N=256 runtime demands it.
4. Separately, approve the ~10-line additive `parent_hash` field in `lret`
   rows (for D3 punctuation).

**[LIN-D3](../ISSUES.md).**
1. Park until A9 is decided.
2. If A9 is fixed, fold the actual-tax column into the same additive reducer
   change as `parent_hash`, or drop it.
3. If A9 is kept, implement the per-clade energy identity in the barcode
   analysis instead.

**Housekeeping (no decision needed).**
- Close [LIN-Q4](../ISSUES.md) when A2 merges.
- Add [LIN-D3](../ISSUES.md) to notebook §5.
- Update the digest's status line.
- Record in `ISSUES.md` that `restricted_mu` is population-union (§2.2), so
  that Option F is not re-proposed.

---

## 10. Verification ledger

| Claim | Status | How |
|---|---|---|
| Overdraw: F_food < 0, `f_priv` < 0 in live cells, negative `death_cum` at `m_scale` ≥ 2.5 | **V** | `energy.py`, `fneg.py`, `Fneg.py`, `msweep.py` in scratch; code read at `C/evoca.c` ~1627–1631 |
| 4b numbers (intake, nominal and actual tax) | **V** | Exact reproduction |
| Nominal/actual gap ≈ 0.4% at `m_scale=1` | **V** | Seed 2 only |
| `restricted_mu` uses the per-step population union; n_active 180–250 normally, 9–42 at crashes | **V** | Code read plus `nactive.py` (single seed per config) |
| Lineage `eff_used_mask` median 60–130 bits | **V** | `nactive.py`, one seed |
| Berkson / scale dependence | **V** | Toy simulation (`berkson.py`), not EvoCA data |
| A6 (no genotype reaches 10 births at 0.03) | **V** | 3 seeds plus an egene-off variant (`a6.py`) |
| P(ΔLUT>0) = 1 − exp(−250·`mu_lut`) | **V** (arithmetic) + **R** (notebook cell 20 measured) | — |
| Egene-driven genotype churn (18–36% per birth) | **A** | From the mutation code; not measured directly |
| Exact per-step barcode tracking and its cost | **V** | `barcode.py`, three configs |
| #3c C−A split by phase-1 `mu_lut`, demographic gap | **V** | Re-analysis of the recorded `results.csv` |
| `lret` memory and collision estimates | **A** | Struct size and birthday bound |
| Shadow perturbs the RNG stream (A2) | **R** here (V in `06`) | — |
| Notebook §3/§4 runs don't enable the shadow | **V** | Code read of cells 20, 26 |
| Resource-driving density ≈ `food_inc/tax` | **R** | `06-seeding-findings.md` §B |

---

## 11. Structural and scaling notes

- **One tool, two decisions.** The per-step tracker (§2.5) should be a single
  module, e.g. `python/lineage_track.py`, owned by this workstream and used by
  both #8 (edges) and H1 (clades). This avoids a second edge-reconstruction
  path next to `scan.py`'s snapshot map.
- **The reducer is the wrong owner for clade fitness.** Keep it as the
  genotype ledger for D3. Clade metrics belong in the tracker, which avoids
  two sources of "lineage net".
- **Scale.** The genotype reducer grows without bound (§2.7). Any N=256 H1
  campaign should either leave genotype tracking off or accept about 2 GB per
  run.
- **Structural risk flag.** The A9 clamp changes a dynamics boundary that many
  recorded campaigns sit on. That is the kind of change where an independent
  review (`/code-review ultra`, N-triggered) is worth considering before
  merge.
