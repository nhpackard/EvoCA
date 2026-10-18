# 07 - ACT-D1: choosing the primary activity objective (guide and reading list)

For N. Written 2026-10-15 after reading the committed Activity-metrics lab notebook
(`git show HEAD:Research/Activity-metrics/Activity-metrics.ipynb`; your working copy has
uncommitted edits and was not touched), its two digests, the Lineage and Resource digests
that bear on activity, the scan analyses, and the metric code in `C/evoca.c`,
`python/evoca_explore.py` and `python/dyn_analysis.py`. I also ran a few short headless
checks (scratch scripts, outside the repo; see section 10). Every claim is tagged
**[verified]** (I checked it this pass, by code reading or by a run), **[recorded]** (it is
in a committed CSV, notebook output or analysis document and I re-read it), or
**[reported]** (someone asserted it and I did not check it).

---

## 0. Summary (one screen)

**The decision.** ACT-D1 asks which readout becomes the primary "adaptive activity"
objective for re-ranking the pure-evo scans and for selection in future campaigns. It
replaces the 5-key pure-evo composite. The seeding pass (`06-seeding-findings.md` D2,
`ISSUES.md` ACT-D1) recommended a **two-channel (`dyn`, `eg`) vector of fixed-space excess,
selected by Pareto box**, with `A_new` and eff-species kept as diagnostics.

**What this pass adds, and why it changes the recommendation.** Three verified facts change
the picture.

1. **The fixed-space shadows have no memory.** `fs_space_step` rebuilds the shadow's
   composition every step from the *observed* composition of that same step, plus a
   mutation smear toward uniform. Two things follow.
   - Whenever the flux is nonzero, the reported excess reduces exactly to
     `excess_pc = (cum. activity of occupied buckets)/D_G − (total mass)/K`. That is a
     **concentration index** (mean cumulative activity per occupied bucket, minus a
     uniform baseline). The shadow contributes nothing except the constant `M/K`.
   - The celebrated "nulls pass exactly" is true **by construction**: zero flux sets the
     shadow equal to the observed composition. It validates the code path. It does not
     validate the statistic.
2. **The per-bucket z-scores are not a selection test.** The count of buckets with z>3
   equals, within 0–2 buckets, the number of buckets whose cumulative share exceeds `1/K`
   (7 EvoCA configs). In a stylized, purely **neutral** population run through the same C
   construction:
   - 113–173 of 500 buckets reach z>3.
   - With skewed bucket accessibility (as dyn inputs are), the neutral model reproduces the
     EvoCA pattern: about 225 buckets occupied, about 110 at z>3, about 385 at z<−3, and a
     cumulative excess of about 2.1e4, against about 2.3–2.9e4 for the pure-evo winner at
     the same N and length.
3. **Genotype-keyed measures are degenerate at productive `mu_lut`.** At `mu_lut=0.01`, a
   birth carries no LUT flip with probability exp(−2.5) ≈ 8%, so genotypes are
   near-singletons. Full-genome species, `A_new`, and every eff-species key I tried
   (population union, per-cell mask, mask on the inputs carrying 99% of queries) stay at
   87–100% of the full-genome species count.

Two further points come from re-reading the recorded data.

- **The eg/dyn trade-off is mostly a `tax_per_egene` factor.** Among the 145 joint
  candidates with egene cost on, `eg` and `dyn` excess are *positively* correlated
  (Spearman 0.68). The "front" comes from pooling these with the cost-free half.
- **"Egene co-evolution amplifies `dyn` excess" is not in the full #1/#2 data.** Mean
  `dyn_excess_pc_slope` is 24.6 for LUT-only, 24.0 for joint with no egene cost, and 24.6
  for joint with cost.

**Options** (section 5):

- **A.** Current fixed-space excess, two channels, Pareto.
- **B.** The same two-channel structure on **repaired, stateful fixed-space shadows**
  (dyn: one two-allele drift process per input; eg: a neighbor mutation kernel), selected by
  quantile-box Pareto **within egene-cost strata**.
- **C.** A repaired `A_new` scalar.
- **D.** An eff-species key.
- **E.** Interim: descriptive readouts only.

**Recommendation.** Decide ACT-D1 as **Option B**.

- Keep the *structure* of the seeding recommendation: two channels, Pareto, no scalar
  composite. Drop its *measurement* claim: the current excess is not adaptive activity.
- Until B lands and passes a neutral calibration, use **Option E**. Record everything, call
  the current `eg`/`dyn` numbers "concentration", report `D_G` and population next to them,
  and select nothing as "adaptive".
- Do not invest in C or D (section 5 explains why).
- The re-ranking campaign (ACT-N1) does not need to wait. Have `run_sim` dump the per-bucket
  occupancy and flux time series, so the runs can be re-scored post hoc under any shadow.

**Read first** (committed notebook indices; in your working copy, cells after 4b are
shifted by +2, so use the headings as anchors):

1. Notebook **4a, cell 27** (objective comparison and eg/dyn scatter), read with section 2.4
   below.
2. **Section 3, cells 20–23** (causal controls; `A_new` = 0).
3. **4b/4c, cells 30–37** (winner vs. frozen control; z histograms; `dyn_report` says
   `silent-churn`).
4. `C/evoca.c` **`fs_space_step` (about lines 3376–3445) and `fs_excess_pc` (about
   3496–3511)**. They are about 70 lines, and they are the whole of finding 1.
5. `Research/Activity-metrics/digests/pure-evo-regime.ipynb` **cells 7–8**, whose headline
   the full data contradict.

Full reading guide in section 7.

---

## 1. What ACT-D1 decides, and what "adaptive activity" means here

**Bedau–Packard framing.**
- A *component* is a heritable unit.
- Its *activity* is its cumulative presence (summed population count over time).
- Evolutionary activity is *adaptive* when components persist and spread more than they
  would under a **neutral shadow**: the same population, births and mutation supply, but
  with fitness-blind reproduction.
- **Excess activity** is real minus shadow.
- The classical statistics are diversity `D`, mean cumulative activity `Ā = ΣA/D`, and
  **new activity** `A_new`. In the Bedau–Snyder–Packard 1998 usage (from memory, please
  check), `A_new` counts the activity of components that have recently crossed a
  neutral-derived threshold. It does not count every component above the threshold.

**EvoCA has three candidate component spaces.**

| Space | Component | Size | Mutation moves a unit to... |
|---|---|---|---|
| Whole genome | FNV hash of LUT ‖ active egenes | open, unbounded | a new hash (essentially always new) |
| `dyn` | (LUT input index, output bit) transitions actually taken by alive cells | 500 fixed | the *other* output of the *same* input; input frequencies are set by the v-field, not by mutation |
| `eg` | ternary (value, mask) key of each active egene | 729 fixed | a one-digit neighbor key (or none, for presence-bit flips) |

**What ACT-D1 must pin down.**
- (i) The component space or spaces.
- (ii) The neutral reference.
- (iii) The aggregate per run.
- (iv) How runs are ranked: scalar, vector, or Pareto.

The review's expectation (`00-review-of-prompt1.md` §7) was "`A_new` on eff-species,
z-scored against the fixed-space shadow". That fixed (i) as the eff-genome, (ii) as the
fixed-space shadow and (iii) as `A_new`. Each of these turns out to be unavailable as
built: sections 2 and 3.

**What the objective is for.** It ranks configurations for "adaptive evolutionary activity"
in ACT-N1 and in the selection step of later campaigns (Resource-driving L1/§B, Lineage H1,
the RD cross-check ACT-N7). It must be:
- valid in both subsystems (LUT and egene);
- a null under no selection, and not merely under no mutation;
- comparable across configurations with different population sizes and dynamics.

The third requirement has not been stated explicitly before. Section 2 shows it bites.

---

## 2. New findings from this pass (the ones that move the decision)

### 2.1 The fixed-space shadow is memoryless; the excess is a concentration index [verified]

**Code** (`C/evoca.c`, `fs_space_step` about 3376–3445, called every step from
`fs_shadow_step` about 3461):

- Each step it sets `q[k] = (1−m)·obs_f[k] + m/K`, plus a negligible `m/Np` drift nudge,
  where `obs_f` is **this step's observed composition**. Nothing carries over from the
  previous step's `q`.
- It then adds `q·T` to the shadow's cumulative activity `nact`.
- With any `m>0`, every bucket gets mass. Once each bucket has accumulated one unit (a
  handful of steps), `D_N = K` and `ΣN = M`, the total mass. (`ΣN = M` was verified:
  `dyn_shadow_total/M = 1.000`.)
- So `fs_excess_pc` (about 3496) returns `A_occ/D_G − M/K`. Here `A_occ` is the cumulative
  observed activity of currently occupied buckets, `D_G` their number, and `K` the space
  size (500 or 729).
- The result does not depend on the flux magnitude, on the drift term, or on any neutral
  dynamics. With `m=0`, `q` equals the observed composition and the excess is 0 exactly.

**Check** (N=96, 3000 steps, 8 configs: the pure-evo winner and variants, the three
causal-control arms, the `dyn` top config 293, and the evo_focus top):

- The C value and the formula agree at the final sample to within 0.6 counts in both
  channels, in every config with nonzero flux.
  - Examples: winner `dyn` 28598.12 vs 28597.63; `eg` 148540.09 vs 148539.71.
  - At a few intermediate samples the two differ by up to a few hundred counts (two configs
    up to about 4,600 or 29,000). I attribute this, without checking, to samples taken
    right after a step with zero births, when the shadow collapses onto the observed
    support.
- The design intent (`Docs/EvoCA_research_directions.md` §2: "Wright–Fisher on 729 alleles
  at the realised mutation flux") was a *stateful* drift process. The implementation is not
  one.

**Consequences.**

- **The reciprocal-freeze "validation" is tautological** (v2: `Scans/2026-05-17_causal_control_v2`;
  digest `neutral-model-methodology` cells 0–3; notebook cells 20–23 and 32).
  - Freezing a subsystem sets its flux to 0, and flux 0 sets shadow = observed by
    construction.
  - The null tests whether the code wires flux correctly. It says nothing about whether a
    *neutral* process with nonzero flux reads near 0.
- **The excess is extensive in population.**
  - At steady state the slope is roughly `T·(1/D_G − 1/K)`. Here `T` is the per-step mass:
    alive cells for `dyn`; active egenes for `eg`, which is about Negene × population.
  - So ranking across configs mixes "concentrated" with "big". This connects to the
    Resource-driving finding that density ≈ `food_inc/tax` in the food-limited regime.
  - The rough steady-state estimate tracked the actual `dyn` slope within a factor of about
    1.4 in all 7 configs with LUT flux. It usually overestimates, because mass in buckets
    that went extinct drops out of `A_occ`; config 293 is the one underestimate.
- **The evo_focus artifact (ACT-Q3) is this identity at `D_G=1`.**
  - In that base `mu_egenome` flips presence bits (flux > 0) while only one key is occupied.
  - The slope is then about `T(1 − 1/729)` ≈ active egenes per step ≈ 1.6e4, matching
    16,914 (notebook cell 11/18) and 16,322 in my run.
- **The `dyn` space has the wrong baseline even in principle.**
  - A LUT flip at input `i` moves a cell between buckets `(i,0)` and `(i,1)`. Input
    frequencies come from the v-field statistics.
  - A uniform `1/500` baseline therefore treats the rarity of input configurations, such as
    four live cells in every ring, as "concentration". So "high `dyn` excess" partly
    rewards dynamically simple v-fields, which use few inputs.
  - Extreme case: a quiescent field (`D_G=1`, as in every frozen-GoL arm here) would
    score the maximum possible as soon as any LUT flux appears.

### 2.2 Per-bucket z-scores: share-above-uniform counts with an uncalibrated false-positive rate [verified, plus a stylized neutral model]

**In EvoCA** (the same 8 runs):
- `n(z>3)` for `dyn` = 98, 110, 72, 84, 85, 72, 121.
- The number of buckets with cumulative share > 1/500 = 98, 111, 72, 86, 86, 74, 122.
- The rest are at z<−3 (for example 402 for the winner).
- So z-tail counts and z-tail mass are functions of the occupancy distribution, like 2.1.
- `eg` is similar, with z>3 slightly below the share count (35 vs 43, 60 vs 71).

**Neutral false-positive check** (stylized model, Python):
- Setup: a Moran-type population, Np=3000, 5% replaced per step by random parents, a
  mutation per birth sending the child to a random bucket, and the C shadow and variance
  formulas applied verbatim.
- Uniform targets: `n(z>3)` = 147–173 of 500, against roughly 0.7 expected for a calibrated
  test. The excess slope is about 0, because a neutral population with uniform accessibility
  spreads to `D_G ≈ K`.
- Skewed targets (geometric over bucket index; a stand-in for skewed input frequencies):
  - `D_G` ≈ 225, `n(z>3)` ≈ 113, `n(z<−3)` ≈ 384.
  - Final cumulative excess ≈ 2.1e4, against the EvoCA winner's 2.35e4 (`run_sim`, cell
    28) and 2.86e4 (`dyn_report`, cell 36) at N=96 / 3000 steps.
- **Caveats.** The model is stylized and its parameters are mine, not fitted. It is an
  existence proof that **a selection-free process produces EvoCA-sized `dyn` readouts**.
  It does not show that EvoCA has no selection.
- This answers ACT-Q9 for the current construction: the z-tail count is not usable as a
  selection test.

### 2.3 Genotype-keyed measures are near-singleton dominated [verified]

**`A_new` blindness (A4), now verified independently.**
- EGENE_only arm, seeds 2 and 3 (the notebook used 0 and 1), through `run_sim`:
  `a_new_final = 0` both times, with `eg` slope +15.5 and +17.0.
- Mechanism: the Channon shadow (LUT-only) holds **one** bucket, the frozen GoL LUT, so its
  p90 equals that single bucket's total activity (thresholds 40,773 and 56,525). No real
  whole-genome component, split over 85–156 egene variants, can exceed it.

**`A_new` threshold collapse (new pathology).**
- In evolving runs the p90 shadow threshold is **2–5 activity units**: winner 5.0 (seeds 2
  and 3), LUT_only 2.0.
- The Channon shadow reproduces with the full `mu_lut` operator, so it is a no-heredity
  model (`Docs/EvoCA_research_directions.md` §2), and its components are ephemeral.
- `A_new` therefore counts real genotypes that have lasted more than a few samples (about
  18% of present genotypes for the winner: 668 of 3,720). That is a persistence-weighted
  diversity count, not new adaptive activity.

**Eff-species keys do not coarsen** (winner, winner at `mu_lut=0.10`, LUT-only; N=96, 3000
steps).

| Key | Winner | `mu_lut`=0.10 | LUT-only |
|---|---|---|---|
| Full LUT species | 3,592 | 3,826 | 4,653 |
| Population-union key (as merged) | union mask 250/250 bits, so equal to full | same | same |
| Per-cell key (LUT & own mask, plus mask) | 3,462 | 3,826 | 4,601 |
| Mask on inputs carrying 99% of queries (150 / 166 / 104 inputs) | 3,347 | 3,826 | 4,011 |

- The root cause is not silent churn. It is that about 92% of births carry at least one LUT
  flip at `mu_lut=0.01`, and most land on frequently used inputs. This is the same
  degeneracy as Lineage A6, P(ΔLUT) = 1−exp(−250·`mu_lut`).
- The 99%-mask row assumes little-endian bit packing of `get_lut`, which I did not verify.
  The conclusion holds either way, because no row drops below 86%.

**A contradiction to fix.**
- `06-seeding-findings.md` A3 says the per-cell used-mask "is the variant that can actually
  coarsen".
- The C design note (`C/evoca.c` about lines 219–233) says the per-cell key *splits*
  genomes, with 3–10x more species. Here it does neither: it is within 4% of full.
- A3's sentence should be corrected.

### 2.4 Re-reading the recorded pure-evo data [verified from committed CSVs]

All numbers are on the 284 joint candidates and 180 LUT-only candidates.

**The composite** (5 keys, positive-clipped, min–max, mean):
- Spearman with `eg` 0.689 and with `dyn` 0.409.
- Its top-5 {28, 199, 26, 80, 69} and the `dyn` top-5 {293, 240, 146, 256, 67} are disjoint.
- These match notebook cell 27.

**The trade-off is a stratification effect.**

| Group | n | median `eg` slope | median `dyn` slope | Spearman(`eg`, `dyn`) |
|---|---|---|---|---|
| joint, `tax_per_egene`=0 | 139 | 0.1 | 16.9 | 0.06 |
| joint, `tax_per_egene`=8e-4 | 145 | 34.4 | 24.8 | **0.68** |
| pooled | 284 | | | 0.34 |

- With egene cost off, the `eg` channel is essentially void. By 2.1, egene key diversity
  `D_G` is high: 253 occupied keys in my run of config 293.
- With cost on, the two channels **agree**.
- The first Pareto layer on (`eg`, `dyn`) has 8 configs: two cost-free (293, 146) at the
  `dyn` end and six costed. Three of the composite's top-5 (28, 199, 80) are on it. So the
  composite did not pick dominated points so much as **one end** of a front whose other end
  is a different egene-cost regime.

**No `dyn` amplification by egene co-evolution in the full data.**
- Mean `dyn` slope: LUT-only 24.6; joint, no cost 24.0; joint, cost 24.6. The medians are
  18.0 / 16.9 / 24.8.
- Within the joint scan the median difference tracks egene *cost*, not egene *mutation*:
  joint with no cost ≈ LUT-only.
- So the claims "freezing egene reduces `dyn` excess" (`Scans/2026-05-17_pure_evo_analysis.md`
  §Campaign #2; digest `pure-evo-regime` headline and cells 7–8) and "the `mu_lut` optimum
  is conditional on egene" rest on top-ranked configs under the composite, not on the pool.
  ACT-Q4 suspected this; the stratification settles it for the recorded data.

### 2.5 `dyn_report` and `dyn` excess disagree structurally, not by calibration (A8) [verified by code reading]

- `dyn_report` (`python/dyn_analysis.py`) reads **low** `dyn_eff` (participation ratio,
  meaning a concentrated transition distribution) as **silent churn**.
- `dyn_excess_pc` is, by 2.1, roughly `T(1/D_G − 1/K)`, which is **high** when the
  distribution is concentrated.
- The two are monotone in the same underlying quantity, the concentration of the `dyn`
  distribution, with opposite interpretations. Their disagreement on the winner (cell 37:
  `silent-churn`, `dyn_eff` 38/500, `dyn` excess strongly positive) is therefore expected
  and will recur. Recalibrating the gates (ACT-Q7) cannot reconcile them.
- The scientific question beneath both is:

  > **Is concentration of the transition distribution evidence of selection or of
  > silence?**

  Only a neutral reference with memory can answer it. This is the core unresolved issue
  behind ACT-D1.

---

## 3. The candidate metrics, one by one

| Candidate | Measures (precisely) | Known pathologies (evidence) | Defensible as primary objective? |
|---|---|---|---|
| `excess_activity_slope` (legacy ΣG−ΣN) | Slope of total whole-genome activity minus total Channon-shadow activity | Inverts at high `mu_lut`, because the shadow is no-heredity (research directions §2). Scope bug: real side is whole genome, shadow is LUT-only (v1 EGENE_only −9.1; notebook cells 20, 23). Sign unstable across N, seed and length (cell 18: +2.7 vs −1.8 recorded) | **No** |
| `excess_pc_slope` (whole-genome, per component) | Slope of ΣG/D_G − ΣN/D_N | Same scope bug (EGENE_only −16.8 at N=96, −46.4 in v2). About 0 everywhere in #1 (median 2e-4). Not a null under freeze (+4,000 cumulative on one component, cell 32) | **No** |
| `dyn_excess_pc_slope` (500 transitions) | `A_occ/D_G − M/K` whenever LUT flux > 0, i.e. mean cumulative activity per used transition minus a uniform baseline, **extensive in population** (2.1). Exactly 0 when flux = 0 | Memoryless shadow, null tautological (2.1). Uniform baseline wrong for input-frequency-structured space (2.1). Neutral model reproduces its magnitude (2.2). Frozen-rich-substrate artifact (#3c arm A +335; memory note; digest `neutral-model-methodology` cell 6–7). Population-confounded (recorded Spearman with density 0.32) | **Not as is.** A valid *descriptor* of transition concentration; becomes a candidate only on a stateful, input-conditioned shadow (Option B) |
| `eg_excess_pc_slope` (729 keys) | Same identity over active-egene keys; T = active egenes | All of the above, plus a `D_G=1` blow-up (ACT-Q3: 1.7e4 in the evo_focus base; cell 18). Void when `tax_per_egene=0` (median 0.1). Presence-bit flips count as flux although they never change a key (`EG_BITS_PER_COMPONENT`, about line 3326) | **Not as is**; same route as `dyn` |
| `A_new` (G1; `a_new_final/slope`) | Number of present whole-genome components whose cumulative activity exceeds the p90 of present Channon-shadow components (`run_sim` about lines 262–282) | Blind in egene-evolving / LUT-frozen runs: 0 in 4 of 4 seeds, verified (2.3). Not a null under freeze (= 1, cell 32). Threshold collapses to 2–5 units, so it becomes persistence-weighted diversity (2.3). p99 never recorded (A7). Differs from the BSP "new activity" semantics | **No**; repair is possible but buys little (Option C) |
| Eff-species / `eff_activity` (F1) | Species keyed by LUT masked with the population-union used mask, plus egenes | Union saturates at 250/250 bits, so eff = full species in every run (cells 7, 30, 39). Alternative keys also fail to coarsen at productive `mu_lut` (2.3) | **No** as an objective; at most a churn diagnostic |
| Per-bucket z-scores (G2) | `(obs_act − nact)/sqrt(Σ m T² q(1−q)/Np)` per bucket | z>3 count ≈ count of buckets above uniform share (2.2). Neutral false-positive rate far above nominal (147–173 of 500). Bidirectional with hundreds of large negative values (cell 33) | **No** as built; meaningful only against a stateful shadow with a calibrated variance |
| `dyn_report` verdict (D1) | Crosses genome-churn rate with `dyn_eff` level (gate 60/500) | Gates tuned on `mu_lut=0.03` demo regimes. Structurally anti-aligned with `dyn` excess (2.5) | **No** as a selection criterion; keep as descriptor |
| Pure-evo 5-key composite | Mean of min–max-normalized, positive-clipped {`excess_pc`, `eg`, `dyn`, `n_distinct_temporal_std`, `unique_top_genomes`} | `eg`-dominated (0.69 vs 0.41). One term is void (`excess_pc` normalized SD 0.06). `unique_top_genomes` anti-correlates with `dyn` (−0.55) and is likely sample-cap bound. Pool-relative. Compensatory (review §6) | **No** (section 4) |
| Legacy `EVO_METRICS` mean | 4 keys including `excess_activity_slope` and `unique_top_genomes` | Internally contradictory (Spearman −0.68 between two of its keys; cell 25/28). Picked `mu_lut=0.06` as an artifact | **No**; keep under an explicit legacy name (ACT-D3) |
| Lineage net-return (D2), Φ (L1) | Adaptive return per lineage; per-capita dissipation | Not activity; hash lineages are near-singletons (A6); nominal vs actual tax (A5) | Second axes, not ACT-D1 candidates |

---

## 4. Why the composite misleads (exact pointers)

1. **It is a weighted `eg` objective in disguise.**
   - Spearman with `eg` 0.69, with `dyn` 0.41 [recorded; notebook 4a, cell 27 output; also
     reproduced].
   - Mechanism: `excess_pc_slope` has almost no spread after normalization (SD 0.06 vs
     0.21–0.25 for the others), so it contributes a near-constant.
   - `unique_top_genomes` has the largest spread (0.245) but anti-correlates with `dyn`
     (−0.55, cell 25).
   - So the effective weights favor `eg` and penalize `dyn`.
2. **Its top-5 shares no config with the `dyn` top-5** (cell 27).
   - The `eg` / composite leaders all have `tax_per_egene=8e-4`. The `dyn` leaders all have
     `tax_per_egene=0`, `food_inc=0.010` and `m_scale` 1.0–1.5.
   - Section 2.4 adds that this is two egene-cost regimes, not a smooth trade-off: within
     the costed half the channels agree (0.68).
3. **The headline levers are channel-specific.**
   - `m_scale=2.5`, the composite winner, is the **worst** marginal `dyn` value (14.5 vs
     26.7–31.5; cell 26). ACT-Q5.
   - The "conditional `mu_lut` optimum" is invisible in marginal means (cell 26) and in
     cost-stratified means (2.4).
4. **The #1 vs #2 contrast confounds egene mutation with egene cost.**
   - #2 zeroed `tax_per_egene` as well as the mutation rates (`Scans/2026-05-17_pure_evo_lut_only/notes.md`;
     notebook cell 28, last bullet).
   - Inside #1, the cost-free half already looks like #2 in `dyn` (2.4).
   - The v2 causal controls keep the cost in all arms, but their comparison is the
     tautological-null construction of 2.1, and their `eg` synergy does not reproduce at
     N=96 (cell 23).
   - The N=96 neighborhood re-ranking reverses the contrast (cell 41: LUT-only `dyn` 8.9–11.3
     vs joint 6.2–6.3). That reversal is consistent with the population term in 2.1:
     LUT-only runs carry no egene tax and hold more cells (4,743 vs 3,791).
5. **Each of the six neighborhood configs wins under some readout** (cell 39/41). A
   composite of disagreeing readouts mostly records which readouts were included.

The composite is not *wrong* about everything. Three of its top-5 lie on the (`eg`, `dyn`)
first Pareto layer. But it hides which channel it is optimizing, and both channels measure
concentration, not adaptation (section 2).

---

## 5. Options for ACT-D1

### Option A: current fixed-space excess, two channels, Pareto (the seeding recommendation)

- **What.** (`dyn_excess_pc_slope`, `eg_excess_pc_slope`) with guards: `mu_lut>0`,
  `D_G≥k`, not extinct. Select by quantile-box via `conepareto` (`qbox`).
- **Pros.**
  - Exists today and is scope-valid in both subsystems.
  - Pareto removes the compensation and silent-weighting problems.
  - `conepareto` is parity-tested (45 tests).
- **Cons.**
  - Both axes are concentration indices with a uniform baseline and no neutral content
    (2.1–2.2), extensive in population, and anti-aligned with `dyn_report` (2.5).
  - Selecting on them steers toward small, concentrated, dynamically simple populations,
    with no evidence that this is adaptation.
- **Build.** Nothing in C. Python: guards plus a `conepareto` call.
- **Verdict.** Not as the *adaptive* objective. Acceptable as an explicitly labeled
  *concentration* descriptor (Option E).

### Option B (recommended): two channels on repaired, stateful fixed-space shadows; Pareto within egene-cost strata

**What.** Replace the memoryless shadow with neutral processes that carry state.

- **`dyn`: a per-input two-allele neutral process.**
  - For each input `i`, track the shadow's output share `p_i(t)`, started from the observed
    share at `t0`.
  - Each step, apply the realized per-birth flip rate at input `i` (flips land on specific
    entries, so this is countable) and Wright–Fisher / Moran drift with the realized
    births.
  - Weight by the **observed** input frequency `n_i(t)`, because input frequency is set by
    the v-field and is not a genetic choice.
  - Excess then measures whether outputs, given inputs, are maintained more consistently
    than drift plus mutation would maintain them. That is selection on the effective rule.
- **`eg`: a stateful 729-key process.**
  - Use the actual mutation kernel: one-digit neighbors for value and mask flips.
    Presence-bit flips change egene count, not key.
  - Births drive the drift, not the whole population every step.
- **Readout per channel.** Either a component-normalized excess or a calibrated z-tail mass,
  chosen by the calibration in next step 3.
- **Selection.** Two-objective `qbox` selection (`conepareto`), run **separately within
  `tax_per_egene` strata**, or with egene cost as an explicit design factor, never pooled
  across it.

**Pros.**
- Keeps the parts of the current program that are right: fixed spaces, two channels, no
  scalar, Pareto.
- Gives a real neutral reference, so "concentration = selection?" (2.5) becomes answerable.
- Calibratable: a stateful neutral process run on the stylized model of 2.2 must give about
  nominal false positives.
- Cheap at run time (K ≤ 729).

**Cons and unresolved.**
- (a) **Input-frequency feedback.**
  - Evolved rules change the v-field and so change input frequencies.
  - Conditioning on observed `n_i(t)` puts that route of selection into the null, which
    makes the test conservative: it detects selection on outputs given inputs, not
    selection acting through which inputs occur.
  - I do not see a clean fix. It should be stated as the estimand.
- (b) **The drift model.** EvoCA's reproduction is local (nearest neighbor with lowest
  food), so effective population size is not the census `Np`. Spatial structure reduces Ne,
  which widens neutral variance. Ne needs estimating (for example from the
  identity-by-descent decay of a marked neutral bit), or the variance must be inflated
  conservatively.
- (c) A **neutral calibration arm** inside EvoCA needs a selection-free mode. The Channon
  shadow population's random-parent reproduction is one template. Until one exists,
  calibration rests on the stylized model.

**Build.**
- C: an additive, opt-in stateful shadow per channel with new getters; the old getters stay
  for comparability. Per-birth flip positions are already counted (`g_lut_bitflips_last`);
  per-input attribution needs one extra counter array (250).
- Python: `run_sim` keys `dyn2_*`, `eg2_*`, plus `D_G`, `T` and flux columns.
- Tests at the public boundary (`run_sim`):
  - frozen null;
  - stylized-neutral calibration (Python reference implementation of the same shadow);
  - a constructed selection fixture (an output forced at one input) that must read
    positive;
  - each test shown to fail without its fix.
- Estimate: a focused branch with a test notebook (pre-merge protocol). This makes the
  shadow **pluggable** and adds a second measurement path for the same quantity, which is a
  structural-risk moment. `/code-review ultra` is worth offering when it is ready.

### Option C: one scalar, repaired `A_new`

- **What.** Key both sides on the LUT only, or build an eg-space `A_new` as well. Fix the
  threshold, for example with a heredity-preserving shadow or BSP's recent-crossing window.
- **Pros.** It is the canonical literature statistic; one number.
- **Cons.**
  - Genotypes are near-singletons at productive `mu_lut` (2.3), so in any open genotype
    space `A_new` measures persistence-weighted diversity.
  - The Channon shadow itself needs the research-directions §2 "make the shadow actually
    neutral" rewrite: a no-mutation copy plus injection at the realized novelty rate.
  - Egene coverage forces a second channel anyway, so the "scalar" advantage disappears.
- **Build.** The C rewrite of `neutral_apply_demography`, plus keying.
- **Verdict.** Not worth doing for ACT-D1. Possibly worth it for H1 at low `mu_lut`, where
  genotypes persist.

### Option D: an eff-species key (per-cell, windowed, or frequency-masked)

- **What.** A coarser genotype unit, then `A_new` or excess on it.
- **Evidence against** (2.3): no key tried gets below 86% of full-genome species at
  productive `mu_lut`.
- **Better reading.** The right coarsening of "effective rule" *is* the `dyn` space
  conditioned on input, which is Option B.
- **Verdict.** Drop as an objective. Keep `eff_activity` as a churn diagnostic. Close
  ACT-Q2 with this evidence.

### Option E (interim, recommended until B is validated): descriptive readouts, no "adaptive" selection

- **What.** Record and report, per channel:
  - `D_G` (occupied buckets);
  - `T` (mass per step);
  - effective number (participation ratio);
  - the current excess, renamed "concentration";
  - flux;
  - validity flags (ACT-Q6).
- Rank, if at all, for coverage. For example, pick configs spread across the (`dyn`
  concentration, `eg` concentration) plane with `conepareto` `wedge` or `band` regions to
  trace the frontier, rather than maximizing.
- **Pros.** Honest, unblocks ACT-N1, costs only Python.
- **Cons.** Does not answer "which config maximizes adaptive activity".

### Summary table

| | Valid in both subsystems | Null under no *selection* (not just no mutation) | Population-comparable | Build cost | Answers "adaptive?" |
|---|---|---|---|---|---|
| A | yes | **no** (2.2) | no | none | no |
| B | yes | yes, after calibration | intensive form chosen in calibration | C + Python branch | yes, for outputs given inputs |
| C | only with second channel | needs shadow rewrite | partly | larger C rewrite | weakly (diversity-like at high `mu_lut`) |
| D | n/a | n/a | n/a | small | no (does not coarsen) |
| E | yes | n/a (makes no claim) | reports `T`, `D_G` | Python only | no, by design |

---

## 6. Recommendation, and what it changes

**Decide ACT-D1:**

> The primary activity objective is a two-channel (`dyn`, `eg`) vector computed against
> **stateful, input-conditioned fixed-space neutral shadows** (Option B). Configurations are
> selected by quantile-box Pareto within egene-cost strata (ACT-D2: yes, with
> stratification). Until those shadows exist and pass a neutral calibration, the current
> `eg`/`dyn` excess is reported as **concentration**, alongside `D_G` and `T`, and is not
> used to claim or select "adaptive" activity (Option E). `A_new`, eff-species, z-tail counts
> and `dyn_report` verdicts are diagnostics only.

**What this changes elsewhere** (suggested edits; none made):

| Document | Current text | Change |
|---|---|---|
| `Research/Activity-metrics/ISSUES.md`, ACT-D1 | "Recommendation: (c) now" (current excess) | Option B, with E as interim; cite sections 2.1–2.2 |
| `ISSUES.md`, ACT-Q9 | "untested" | Answered for the current construction: `n(z>3)` ≈ share-above-uniform count; stylized neutral false positives 113–173 of 500 |
| `ISSUES.md`, ACT-Q2 / Q7 / Q4 | open | Q2: close (no key coarsens). Q7: structural, not calibration. Q4: no pooled `dyn` difference (24.0–24.6) |
| `06-seeding-findings.md` A3 | "per-cell mask can actually coarsen" | Contradicted by the C note and by 2.3 |
| `06-seeding-findings.md` D2 | "two-channel ... on fixed-space shadows" | Two-channel on *repaired* shadows |
| `00-review-of-prompt1.md` §7 expectation | "A_new on eff-species, z-scored against the fixed-space shadow" | Superseded; each of the three parts is unavailable as built (2.1–2.3) |
| `digests/neutral-model-methodology.ipynb` headline | "pass the null controls exactly ... trustworthy instrument" | The null is by construction; the instrument measures concentration |
| `digests/pure-evo-regime.ipynb`, `Scans/2026-05-17_pure_evo_analysis.md` §Campaign #2 | "egene-freeze reduces `dyn` excess" | Not supported by pooled or stratified means |
| Memory note `dyn-shadow-frozen-rich-substrate` | "neutral baseline collapses toward 0" | Mechanism is the memoryless construction (shadow = observed when flux = 0); worth updating |

**What would change my mind.**
- If the stateful shadow of Option B, applied to EvoCA runs, gives excess and z-tails close
  to the current values, the current numbers were approximately right, for a different
  reason than believed, and Option A becomes acceptable after the fact.
- If B finds no selection on outputs given inputs at the pure-evo winners, the headline
  "selection is real and measurable in effective dynamics" (pure-evo analysis) has to be
  withdrawn. Activity work would then shift to the eg channel and to lineage net-return.

---

## 7. Reading guide (ordered)

Notebook = `Research/Activity-metrics/Activity-metrics.ipynb` (44 cells since 2026-10-18, when the
commented "run interactively" cells were added; all cell numbers in this guide are 0-based and were
remapped then). Headings are the stable anchors.

| # | Where | What to look for | Bearing on ACT-D1 |
|---|---|---|---|
| 1 | Notebook **§4a, cell 27** (top-5 tables, Spearman line, eg-vs-dyn scatter colored by `tax_per_egene`) and **cell 28** (observations) | The two clouds in the scatter: blue (`tax_per_egene=0`) along the `dyn` axis, red (8e-4) on a rising diagonal | The "trade-off" is a cost-regime split; within red the channels agree (0.68, section 2.4). Decides whether Pareto is needed *within* strata or only across them |
| 2 | Notebook **§4a, cells 25–26** (Spearman matrix; marginal `dyn` means by `mu_lut` and `m_scale`) | `unique_top_genomes` vs `dyn` −0.55; flat `mu_lut` marginals; `m_scale=2.5` the worst `dyn` | Why the composite and legacy `EVO_METRICS` mislead; ACT-Q4/Q5 |
| 3 | Notebook **§3, cells 20–23** (v1/v2 recorded tables, N=96 reproduction, null asserts) | `a_new_final` = 0.0 for EGENE_only; nulls hold to 1e-3 | `A_new` blindness (now verified on 2 more seeds). Read the null assertions knowing they are true by construction (2.1) |
| 4 | Notebook **§4b, cells 30–33** (winner vs frozen control panels, summary table, z histograms) | Control: every layer zero, except whole-genome `excess_pc` (+4,000) and `A_new` (1). Winner: z histogram with a large negative mass | The frozen "null holds at every layer" is flux = 0. The bidirectional z pattern is what a share-above-uniform statistic looks like (2.2). Note the control also drops `tax_per_egene` (cell 28 text), so it is not cost-matched |
| 5 | Notebook **§4c, cells 36–37** (`dyn_report`) | `silent-churn`, `dyn_eff` 38/500, while `dyn_excess_pc_final` ≈ 2.9e4 | The A8 disagreement; structural, not calibration (2.5) |
| 6 | Notebook **§4d, cells 39–41** (six-config re-ranking) | Every config wins under some readout; LUT-only `dyn` > joint at N=96; eff = full species | Objective choice drives the answer; the population term in 2.1 explains the LUT-only reversal |
| 7 | Notebook **§2, cells 11–18** (evo_focus winner, `mu_lut` sweep, cfg 98) | `eg_excess_pc_slope` ≈ 17,000 with `eg_keys_present = 1`; legacy excess rising with `mu_lut` while `dyn` is flat | The `D_G=1` case of the 2.1 identity (ACT-Q3); legacy metric retraction |
| 8 | Notebook **§1, cells 7–9** | eff species = full species (1,614 = 1,614); union 250/250 | ACT-Q2 / A3; with 2.3, drop eff keys as an objective |
| 9 | Notebook **§5, cell 42** | The "settled core" bullet | The bullet "fixed-space shadows pass their null controls" is the one this guide revises |
| 10 | `C/evoca.c` **`fs_space_step` (about 3376–3445)**, **`fs_shadow_step` (about 3461)**, **`fs_excess_pc` (about 3496)**, **`fs_zscores` (about 3573)** | `q[k]` built from this step's `obs_f`, with no use of the previous `q` | The mechanism of 2.1–2.2; about 120 lines in total |
| 11 | `python/evoca_explore.py` **`run_sim`, about lines 262–282 (A_new) and 386–447 (excess keys)** | p90 of present shadow components; slopes of cumulative series | `A_new` construction; threshold collapse (2.3) |
| 12 | `python/dyn_analysis.py` **module docstring and `dyn_report` docstring** | "Concentration = LOW & FLAT `dyn_eff`" means silent | The opposite reading of concentration from `dyn` excess (2.5) |

## 8. Other results worth examining

| Path | What to look for | Bearing |
|---|---|---|
| `Research/Activity-metrics/digests/neutral-model-methodology.ipynb`, cells 0, 3, 6–7 | Headline "pass the null controls exactly"; Bug 2 (#3c arm A `dyn` +335) | The null claim needs the 2.1 qualification. Bug 2 is another face of the same memoryless construction |
| `Research/Activity-metrics/digests/pure-evo-regime.ipynb`, cells 0, 6, 7–8 | "Egene-freeze reduces `dyn`-excess"; the `mu_lut`-by-regime plot | Contradicted by pooled means (2.4); the plot's own lines cross |
| `Research/Activity-metrics/Scans/2026-05-17_pure_evo_analysis.md`, §Campaign #2, §Causal-control v2, §Caveats | "Markedly smaller (3–20 vs 23–60)" is a top-config comparison; the frozen-substrate confound | Source of the amplification and conditional-`mu_lut` claims |
| `Research/Activity-metrics/Scans/2026-05-17_causal_control_v2/notes.md` and `results.csv` | 4 seeds, N=256, 15k steps | The only multi-seed evidence; same tautological-null caveat |
| `Research/Activity-metrics/Scans/2026-04-27_evo_focus/notes.md`, §"Excess-activity-slope dropped" | The "Channon quirk" explanation | Historical; it is the no-heredity shadow (research directions §2) |
| `Research/Lineage-co-evolution/digests/coevolution-substrate.ipynb`, cells 0, 2, 6 | #3c: C−A `eg` +0.6 (22/30); `eg` slopes ≲ 1; arm A `dyn` artifact | `eg` signal under controlled substrate is tiny next to #1's 20–135. Consistent with the excess tracking occupancy rather than selection strength |
| `Research/Resource-driving/digests/viability-brackets.ipynb`, cell 6 | Top-quartile `dyn` configs by `m_scale` and `food_inc` | R1 ranked on `dyn` alone; re-read per channel and with the population term in mind (ACT-Q5) |
| `Docs/EvoCA_research_directions.md` §2 (about lines 230–312) and §3 | Design intent: WF on 729 alleles at realized flux; "make the shadow actually neutral" for the Channon shadow | The repair target for Option B; Option C's prerequisite |
| `Docs/Oct2026/05-cone-pareto-package.md` and `~/Projects/conepareto/README.md` | `qbox`, `wedge`, `band`; `on_shortfall`; validity guard on the first objective only | Use `qbox` for "good at both", `wedge`/`band` to trace the front. Guards on `D_G` must be a pre-filter, since the validity guard covers only the first objective |
| `Research/Activity-metrics/ISSUES.md` ACT-Q1–Q9, `06-seeding-findings.md` A2–A8 | The open items this guide resolves or reframes (section 6 table) | Bookkeeping |
| `tests/notebooks/dyn-promote-demo.ipynb` (not re-read this pass) | The `silent-churn` vs `functional` demo regimes the D1 gates were tuned on | Useful as a first test bed for Option B: does the stateful shadow separate them? |

---

## 9. Next steps if you accept the recommendation

| # | Step | Kind | Notes |
|---|---|---|---|
| 1 | Record ACT-D1 = B (interim E) and ACT-D2 = yes-with-strata in `ISSUES.md`; apply the section 6 corrections to `06-seeding-findings.md` and the two digests (marked as revisions) | docs | Regenerate the board with `tools/build_board.py` |
| 2 | Add notebook **section 6, "Auditing the fixed-space shadow"**: the identity check, the z vs share count, the stylized neutral model, the cost-stratified re-read, the `A_new` threshold values | analysis | Makes this pass's evidence durable; it currently lives in scratch scripts (section 10). About 1 s per run |
| 3 | **Option B branch**: stateful `dyn` (per-input two-allele) and `eg` (neighbor-kernel) shadows, additive and opt-in. Python reference implementation of the same shadows for the stylized neutral model. Tests through `run_sim`: (a) frozen null; (b) stylized-neutral false-positive rate ≤ nominal; (c) constructed selection fixture reads positive; (d) old keys unchanged (bit parity). Each test shown to fail without its fix | needs dev (C + Python) | Structural-risk moment (new pluggable shadow, second path for the same quantity): offer `/code-review ultra`. Pre-merge test notebook. Ne estimation (5-B caveat b) is part of the calibration |
| 4 | **ACT-N2 extended**: `run_sim` emits per channel `D_G`, `T`, participation ratio, flux, validity flags. Add an option to **dump per-sample bucket occupancy (500 + 729 counts) and per-interval births and flips** | needs dev (Python) | Lets ACT-N1 run now and be re-scored post hoc under any shadow. At 80 samples it is about 100k ints per run |
| 5 | **ACT-N1 campaign**: pure-evo grid with `tax_per_egene` as an explicit factor (matched across joint and LUT-only), 3+ seeds, N=128, 5000 steps, with dumps | campaign | Laptop scale. Rank under E now and under B when it lands; `conepareto` `qbox` within strata, `wedge` for the frontier figure |
| 6 | **ACT-D3**: retire the composite and `EVO_METRICS` default after step 5; keep the legacy set as `EVO_METRICS_LEGACY` | dev (small) | Prevents silent mixing of old and new rankings |
| 7 | Downstream re-reads under the new objective: R1 per channel (ACT-Q5), RD winners (ACT-N7), #3c / H1 inputs | analysis | Coordinate with Resource-driving, Reaction-diffusion, Lineage |

Anti-patterns to avoid while B is pending:
- Pooling joint configs across `tax_per_egene` in one Pareto pool.
- Treating a frozen-arm null as validation of a shadow.
- Ranking on any excess without `T` and `D_G` beside it.
- Investing in eff-species keys or the Channon shadow for ACT-D1.

---

## 10. Verification ledger

**Verified this pass.** Scratch scripts are in `/private/tmp/claude-505/act/` (outside the
repo, not durable; next step 2 moves them into the notebook). All runs were N=96 headless
through `build_sim` / `run_sim`.

- `verify_fs.py`:
  - the identity `excess_pc = A_occ/D_G − M/K` at the final sample (within 0.6 counts,
    both channels, every config with nonzero flux);
  - shadow mass = observed mass (ratio 1.000);
  - `n(z>3)` ≈ count of buckets above uniform share (7 configs, within 0–2 for `dyn`).
- `neutral_z.py` and `neutral_wf.py`: stylized neutral false positives and magnitudes
  (2.2). These are *stylized*, with parameters chosen by me.
- `anew.py`: `A_new = 0` in EGENE_only, seeds 2 and 3; p90 thresholds 2–5 in evolving runs;
  single-bucket Channon shadow in EGENE_only.
- `effkey.py`: no eff key coarsens (2.3). The bit-order assumption for the 99% mask is
  unverified.
- Recorded-CSV re-reads (2.4): Spearman values, stratified means, first Pareto layer,
  composite overlap.
- Code reading: `fs_space_step`, `fs_excess_pc`, `fs_zscores`, `run_sim` `A_new` block,
  `dyn_report` docstrings, the eff union-mask design note.

**Recorded, re-read, not re-run:** v1/v2 causal-control tables; #3c arm table; R1
top-quartile counts; evo_focus notes.

**Reported, not checked:**
- the #3c arm A mechanism (I offer the memoryless-shadow explanation as a hypothesis only);
- A5 and A6 (Lineage);
- the "3–10x more species" figure in the C note (my one run found about 0.96x);
- the BSP 1998 definition of `A_new` (from memory);
- intermediate-sample deviations from the identity (cause not investigated).
