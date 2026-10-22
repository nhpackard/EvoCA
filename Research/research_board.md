# EvoCA research board

**Living coordination surface for the EvoCA research program** — in-flight
status, decisions-needed, and the running log. This is the *process*
layer; it is distinct from and points to the *results* layer
(`Research/README.md` = curated results index / reading order;
`Research/<campaign>/` digests = authoritative findings). The board
carries campaign **status** (headline + link), not the full result, so
nothing lives in two places that can drift. Big-picture design + roadmap:
[`../Docs/checkpoint 3 July 2026.md`](../Docs/checkpoint%203%20July%202026.md) (§K′/§K″).

Conventions: Tier-2 / dynamics-altering branches never auto-merge (human
merge gate); Tier-1 additive branches merge under the standing approval
below. Each code branch runs in its own `git worktree`; analysis-only work
runs on `main`. Read "Decisions needed" first; everything else is status.

Pre-checkpoint (2026-05) board state is archived, not authoritative:
see [`research_board_archival.md`](research_board_archival.md).

> ### How to read the codes
>
> **Action codes** such as `D2`, `G1`, `H1`, `L1` are `<letter><n>`: the
> **letter** is a section **§A–§L** of
> [`../Docs/checkpoint 3 July 2026.md`](../Docs/checkpoint%203%20July%202026.md)
> (the checkpoint) and the **number** is one action inside that section.
> The full list is the **§I′ consolidated action table (~line 1470)**; the
> prioritized roadmap is **§K′ (~line 1499)** and the parallel execution
> plan is **§K″ (~line 1594)**. Example: `G1` = section G (per-gene neutral
> model), action 1 = Bedau-Packard A_new. Letters are historical labels,
> not meaningful in themselves; use the table below for the plain meaning.
>
> **`§X` references** (e.g. "§B figure", "§H program") point to the
> discussion section with that letter in the same checkpoint doc
> (§A env/agent split, §B resource-driving inverted-U, §C food_per, §D
> innovation/punctuation, §E chem, §F relevant activity, §G per-gene neutral
> model, §H lineage backbone / adaptive-vs-drift, §I dev logistics, §K
> roadmap, §L energy gradient). Section headings start around lines
> 213 (§A), 278 (§B), 400 (§D), 616 (§H), 738 (§L).
>
> **`W1-a` ... `W1-d`** are *decision IDs* for the Wave-1 review (see
> "Decisions needed"), unrelated to the action codes.
>
> **Old IDs in the archive** (S1, S2a-e, #1-#8, R1, M1, D1-D4) are an
> **earlier, independent numbering** from the 2026-05 board. In particular
> **D1-D4 there are NOT checkpoint §D actions** (checkpoint D1 = promote
> dyn_activity; archive D1 = a worktree-base quirk). Never mix them.

---

## Decisions needed (updated 2026-10-15) — read this first

| # | Workstream | Decision / status |
|---|---|---|
| W1-d | Integration | **RESOLVED by N 2026-10-15.** Tier-1 merges after a green combined suite + bit-for-bit default-param parity, one combined review notebook. Tier-2 (`resource-driving`, `mu-genome`, `chem`) **merged 2026-10-15** on N's authorization after the same checks (see Dev branches). |
| W1-a | Activity-metrics | **Proceeding on default** (review §3e): add A_new / z-scores as ranking *options* without changing `EVO_METRICS` defaults (the new-activity choice is the §7 candidates' job); keep legacy excess for comparability. Override any time. |
| W1-b | Activity-metrics | **Proceeding on default:** keep the p90 neutral threshold for A_new and record p99 alongside. Override any time. |
| W1-c | Activity-metrics (eff_activity) | **Proceeding on default:** keep the per-cell used-mask diagnostic (+8.4 MB @N=512) as **opt-in**. Override any time. |

Branch deletions and Tier-1 integration: **done 2026-10-15** (see Log).

---

## Wave 1 (2026-07-22) — status (the July batch of 11 branches)

11 agents spawned from checkpoint §K″, each in its own worktree on a
named branch off `main@f6ccdab`. **All complete, all suites green,
nothing merged/pushed** (human gate). **Update 2026-10-15:** all 7 Tier-1
branches merged (`5de867d`, `6548106`); Tier 2 still at N's gate.

**Tier 1 — main-merge candidates (additive / no-op):**

| Branch | Commit | Action | Tests |
|---|---|---|---|
| `env-agent-split` | a333434 | A1 env/agent split + A3 tax_ring doc | 42 |
| `dyn-promote` | 73ba66e | D1 dyn_activity → primary signal (extracts `build_sim`) | 42 |
| `metric-anew` | 92b3987 | G1 Bedau–Packard A_new | 42 |
| `metric-flux` | 1beaaf0 | L1 dissipation Φ accountant | 42 |
| `metric-zscore` | 90056c1 | G2 per-bucket eg/dyn z-scores | 45 |
| `eff-activity` | fd57b0b | F1 effective-rule species (union-mask key) | 47 |
| `lineage-return` | 4cff481 | D2 per-lineage net-return | 47 |

**Tier 2 — feature branches (own merge gates; review in worktree; keep):**

| Branch | Commit | Action | Tests |
|---|---|---|---|
| `resource-driving` | 72adb59 | C1 food_per (temporal regen) | 50 |
| `mu-genome` | 6c0a360 | A2 per-cell rate-gene storage (mutation off) | 42 |
| `chem` | 74cb1e5 | E1/E3 N_chem=1 parity + sparse store | 45 |

**Drafts:** `research-drafts` (84fc29b) → 6 candidate digests in
`Research/_drafts/` (5 real-figure, 1 skeleton).

**Integration hotspots (hand-resolution, not auto-merge):**
`python/evoca_explore.py` (dyn-promote `build_sim` refactor + metric-anew
+ metric-flux → merge dyn-promote first, rebase the metric branches on
it); `C/evoca.c` phase loops (metric-flux 2c/4, lineage-return 2c/3,
eff-activity 1, resource-driving 2, mu-genome 4, chem 2c/3/4 — additive
but overlapping).

---
---

## Active workstreams (post-checkpoint)

The checkpoint §K′/§K″ roadmap is the current program. N approved six
workstream names on 2026-10-15 (review `../Docs/Oct2026/00-review-of-prompt1.md`
§3a, §11): **Activity-metrics** (was "Wave"), **Resource-driving**,
**Lineage-co-evolution**, **Reaction-diffusion**, **Criticality-bifurcation**,
**Ring-complexity**; **Chem** gets a workstream dir later, once it has
results. **"Wave 1" now names only the July 2026 batch of 11 branches**
(status section above), not a scientific workstream; "Wave 2" is the set of
follow-on actions gated on Wave-1 merges.

**Kind legend:** *built* = already implemented on a branch, lands when that
branch merges; *needs dev* = not yet written (analysis or code); *campaign* =
a scan/run program using already-built tools; *infra* = no workstream.

### Activity-metrics (was "Wave")

| Code | Plain meaning | Kind | Status / blocker |
|---|---|---|---|
| G1 | Bedau-Packard A_new: new-activity, per component, shared threshold | built | **merged** 2026-10-15 (`5de867d`) |
| G2 | per-bucket z-scores in the fixed eg(729)/dyn(500) spaces | built | **merged** 2026-10-15 (`5de867d`) |
| F1 | `eff_activity`: species keyed by the effective (actually-used) rule | built | **merged** 2026-10-15 (`5de867d`) |
| D1 | promote `dyn_activity`/`dyn_excess` to the primary phenotypic signal | built | branch `dyn-promote` (extracts `build_sim`); merges first, metric branches rebase on it |
| H1 | adaptive-vs-drift program: lineage + A_new + shadow-beating test on the §B sweep | needs dev (analysis) | blocked on D2 + G1 merging |

### Resource-driving

| Code | Plain meaning | Kind | Status / blocker |
|---|---|---|---|
| C1 | `food_per`: temporal regeneration period for food | built | **merged** 2026-10-15 (`c9ecf32`) |
| C2/L2 | float `env_mask` + relaxation toward mask; spatial gradient niches | needs dev | new branch from `main` (C1 merged) |
| L1 | dissipation accountant: sum of tax + death + overflow; per-capita throughput Φ (the §B x-axis) | built | **merged** 2026-10-15 (`5de867d`) |
| §B figure | scan summary of outcome vs `food_inc` and `food_inc/tax` (the inverted-U test); needs B1/B2 notebook | campaign | blocked on L1 + G1 merging |
| C3 | "temporal lumpiness gives evolutionary opportunity?" as a falsifiable campaign | campaign | unblocked (C1 merged); RES-N5 |

### Lineage-co-evolution

| Code | Plain meaning | Kind | Status / blocker |
|---|---|---|---|
| D2 | per-lineage net metabolic return reducer on the lineage field | built | **merged** 2026-10-15 (`6548106`) |
| D3 | single-run punctuation detector: lineage net-return dip-then-staircase at `Negene` increments | needs dev | blocked on D2 merging |
| H1 | (see Activity-metrics) adaptive-vs-drift test | needs dev | blocked on D2 + G1 |

Old #8 (binary-indicator lineage metric, on hold) is superseded by D2 + H1.

### Infrastructure and Chem (no workstream yet)

| Code | Plain meaning | Kind | Status / blocker |
|---|---|---|---|
| A1 | env/agent two-dict metaparam split (derived from one recipe) | built | **merged** 2026-10-15 (`5de867d`) |
| A3 | `model.md` Phase 2c: add the `tax_ring` term (doc fix) | built | **merged** 2026-10-15 (`5de867d`) |
| A2 | `mu-genome`: per-cell rate-gene storage (mutation off), later evolvable rates | built (storage only); rest needs dev | storage **merged** 2026-10-15 (`28e28be`); evolvable rates need dev |
| E1/E3 | chem branch: N_chem machinery (sparse store, N_chem=1 parity) | built (milestone 1) | **merged** 2026-10-15 (`702143a`), N_chem=1 only |
| E2 | chem milestone 2: food + waste (N_chem=2), needs a reaction-execution engine | needs dev | new branch from `main`; high risk |

### Criticality-bifurcation (carried 2026-07-22)

Parameter-sweep bifurcation diagrams over global time series (pop,
activity); harness `python/bifurcation.py` already merged. "Bifurcation"
generally: Feigenbaum period-doubling has **not** been observed in EvoCA,
so we don't presume that route. Subsumes the control-parameter /
finite-size exponent questions from `../Docs/EvoCA_research_directions.md`
§5. Kind: campaign. Status: no sweeps launched yet.

### Ring-complexity (carried 2026-07-22)

Genelife ring-ladder A/B (old #5) + ring-complexity tax + the ring-stability
two-track plan (`../Docs/EvoCA_research_stability.md`: Track-1 ring-class
dominance analysis-on-main, Track-2 bounded ring-mutation in its own
worktree). See memory `ring-stability-two-track-plan`. Kind: campaign +
needs dev (Track 2). Status: not started.

### Reaction-diffusion (NEW topic 2026-07-22)

Curate and showcase EvoCA's photogenic RD regimes with example `.evoca`
recipes in a digest notebook. **Distinct from** campaign #4 (RD-robustness,
which found RD structure is *not* robust under pure-evo optimization):
this documents the RD dynamics' existence and visual character as a
first-class result. Source recipes in `Reaction-diffusion/Runs/*R-D*` and
`Reaction-diffusion/Runs/2026-04-27_*`. Kind: campaign. Status: digest to be drafted.

---

## Workstream issues and directions (generated)

Pulled from each `Research/<Workstream>/ISSUES.md` — the single source; edit there, then run
`python3 tools/build_board.py` (`--check` exits 1 if this section is stale). Do not edit below by hand.

<!-- BEGIN workstream-issues (generated by tools/build_board.py; edit Research/<WS>/ISSUES.md) -->

### Activity-metrics

Source: [`Activity-metrics/ISSUES.md`](Activity-metrics/ISSUES.md)

_Updated 2026-10-15_

Workstream question: what is the right measure of adaptive (vs. drift/turnover) evolutionary activity in EvoCA, and which configurations maximize it? Evidence labeled "example run" below comes from `Activity-metrics.ipynb` (N=96, 3000 steps, one or two seeds; not campaign results).

#### Decisions needed (N)
- **ACT-D1**: Which readout is the primary activity objective for re-ranking the pure-evo scans (and for any new campaign's selection step)? Candidates are listed in `Docs/Oct2026/00-review-of-prompt1.md` section 7; the review's expectation was "A_new on eff-species, z-scored against the fixed-space shadow". First-look evidence from the example runs weakens that single-scalar plan:
  - `A_new` as merged compares whole-genome-keyed real activity with the LUT-only shadow, so it is blind in egene-evolving runs: in the EGENE_only arm it is 0 while `eg` excess is +9 to +21 (notebook section 3).
  - `eff_activity` (F1) with the population-union mask collapses nothing: union mask is 250/250 bits and eff species equal full-genome species in every evolving example run (sections 1, 2, 4), so "A_new on eff-species" currently equals "A_new on full genome".
  - Per-bucket z-scores (G2) are valid in both subsystems, exactly 0 under a freeze, but are bidirectional (hundreds of buckets at z below -60 plus a positive tail) and need an aggregation rule (count above z=3? tail mass?) with a calibrated false-positive rate.
  - **Options:** (a) one scalar, `A_new` after repairing it (ACT-Q1, ACT-Q2); (b) one scalar built from z-tail counts; (c) a **two-channel vector** (`dyn`, `eg`), each as per-component excess or z-tail count on the fixed-space shadows, with `mu_lut>0` and diversity guards, selected by ACT-D2. **Recommendation: (c) now, as the recorded #1 data already show the two channels trade off** (below); keep `A_new` and eff species as diagnostics until ACT-Q1/Q2 are resolved. Blocks: the re-ranking campaign (ACT-N1), H1.
- **ACT-D2**: Adopt Pareto selection (`~/Projects/conepareto`, report `Docs/Oct2026/05-cone-pareto-package.md`) for two-objective questions in this workstream? The recorded `pure_evo_joint` candidates (284) give a concrete case: top-5 by the 5-key composite and top-5 by `dyn_excess_pc_slope` share no config; the composite correlates 0.69 with `eg` and 0.41 with `dyn` excess; `eg` winners have `tax_per_egene = 8e-4` and `food_inc` 0.013-0.018, `dyn` winners have `tax_per_egene = 0`, `food_inc = 0.010`, `m_scale` 1.0-1.5 (notebook 4a). **Recommendation: yes, for (`dyn`, `eg`) selection, using a box or quantile-box region so a config strong in one channel and dead in the other is not admitted (the compensation weakness of the mean, review section 6); keep scalarization only for single-objective composites.** The other two-objective questions (activity vs. Φ from L1, activity vs. lineage net-return from D2, RD structure vs. activity) belong to Resource-driving, Lineage-co-evolution, and Reaction-diffusion respectively. Blocks: how ACT-N1 selects follow-up configs.
- **ACT-D3**: Change the `EVO_METRICS` default once D1 is settled? It currently contains legacy `excess_activity_slope` (invalid for egene-driven runs) and `unique_top_genomes` (Spearman -0.55 with `dyn` excess and -0.68 with the legacy excess among the joint candidates), so the mean mixes turnover with its opposite. W1-a (2026-10-15) keeps defaults unchanged. **Recommendation: leave unchanged until D1, then replace with the D1 readouts and keep the legacy set under an explicit name for comparability.** Blocks: nothing urgent; it keeps old and new rankings from silently mixing.

#### Open questions
- **ACT-Q1**: Can `A_new` be made valid for egene-driven runs? The real side is keyed by the whole genome while `N` (Channon shadow) carries LUT bytes only (established finding, `Docs/devlog/2026-05-17_research-campaigns.md`). Observed: EGENE_only arm `A_new` = 0 (both seeds) with positive `eg` excess; a frozen-mutation control gives `A_new` = 1, not 0. Possible repairs: key both sides on the LUT only, or build a shadow in the eff key space, or define `A_new` on the fixed spaces (count of buckets above a shared z or activity threshold). Also unresolved: p90 vs p99 threshold (W1-b says keep p90 and record p99 alongside; I found no p99 key in `run_sim` on main, so that recording is not yet implemented).
- **ACT-Q2**: Does `eff_activity` (F1) ever collapse species? The union-mask key saturates (250/250 bits exercised) in every evolving example run, so eff species equal full-genome species. The F1 design note says the per-cell mask splits identical genomes by neighborhood history (3-10x more species). Is a windowed or recent-use union, or a per-cell mask with a coarser key, the intended coarsening? Test on runs where silent churn is known to dominate (the D1 `silent-churn` regime in `tests/notebooks/dyn-promote-demo.ipynb`).
- **ACT-Q3**: `eg_excess_pc` is an artifact when the real run has one egene key present: in the evo_focus base (`mu_egenome` set, `mu_egene` unset) it reads about +1.7e4 per tick with one key populated, versus 12-65 in egene-evolving runs. Per-component normalization (`sum G / D_G - sum N / D_N`) explodes when `D_G` is 1. Needs a guard or a reported `D_G`; same family as the `dyn_excess_pc` frozen-rich-substrate artifact.
- **ACT-Q4**: Is the "mu_lut optimum is conditional on egene" claim robust? (a) Scan #2 freezes the egene by zeroing `tax_per_egene` as well as the mutation rates, so #1 vs #2 confounds egene mutation with egene cost (the v2 controls keep `tax_per_egene = 8e-4` in all arms and do not). (b) Marginal mean `dyn_excess_pc_slope` by `mu_lut` is flat and non-monotone in both scans (joint 21.3, 26.0, 28.9, 23.9, 22.1; LUT-only 27.5, 21.2, 26.2, 25.0, 23.6), so the claim lives in the composite's top ranks. (c) Example runs at N=96 reverse the #1-vs-#2 `dyn` contrast (joint 6.2-6.3 vs LUT-only 8.9-11.3). Needs a matched-cost re-run with replicates (ACT-N1).
- **ACT-Q5**: Is `m_scale = 2.5` (the composite winner, at the grid ceiling) a property of the `eg` channel only? Marginal mean `dyn_excess_pc_slope` by `m_scale` (joint): 26.7, 31.5, 26.3, 20.6, 14.5 for 1.0-2.5, so 2.5 is the worst `dyn` value on average. R1 (the `m_scale` bracket, owned by Resource-driving) should be re-read per channel.
- **ACT-Q6**: Whole-genome excess is not interpretable under freezes or egene evolution: negative (-17 to -46) in EGENE_only and not a null in a single-genome frozen control (cumulative `excess_pc` about +4000 from a 5% real/shadow difference on one component). Should `run_sim` emit a validity flag (egene or `mu_lut=0` conditions, extinct, `D<k`) next to every excess key so scans cannot rank on an invalid readout? Extinct runs print large meaningless constants (about -1.9e5 for `eg`, -1.2e5 for `dyn` in extinct N=96 runs).
- **ACT-Q7**: `dyn_report` (D1) labels the pure-evo winner `silent-churn` (`dyn_eff` 38/500, under the 60 gate) while its `dyn` excess is strongly positive and 106/500 buckets have z>3. The gates were tuned on `mu_lut=0.03` demo regimes. Recalibrate on pure-evo winners before using the verdict for selection, or drop the verdict from ranking and keep it as a diagnostic.
- **ACT-Q8**: Seed and scale sensitivity. Scans used one seed (`seed=0`); the example runs show large seed spread (EGENE_only `eg` slope 9.0 vs 21.5 across two seeds at N=96), and the `eg` synergy (JOINT above EGENE_only) recorded at N=256 / 15000 steps does not appear at N=96 / 3000. Any re-ranking needs at least 3 seeds and a scale check.
- **ACT-Q9**: False-positive rate of the z-score test. Frozen arms give z=0 trivially; the rate of z>3 under neutral drift with mutation flux >0 (the case that matters) is untested. A selection-free run with mutation on (for example the shadow compared with a second independent shadow) would calibrate it.

#### Next directions
- **ACT-N1**: Re-run the pure-evo grid (joint and LUT-only, matched `tax_per_egene`) with the readouts now on main (`dyn`/`eg` excess, `A_new`, z-tail counts, eff species) and 3+ seeds, then rank under each candidate objective and under the ACT-D2 selection. Kind: campaign. Blockers: ACT-D1 (which readouts to record as primary), ACT-D2, and ACT-N2 for z-tail and diversity columns. Cost note: an N=96, 3000-step run with every probe takes about 1 s here, so 300 configs x 3 seeds at N=128 and 5000 steps is laptop-scale. Cross-ref checkpoint D1, G1, G2.
- **ACT-N2**: Promote the notebook-local `trace()` quantities into `run_sim`: z-tail counts and tail mass, number of `eg`/`dyn` buckets present, eff species count, a per-run validity flag (ACT-Q6), and an extinct mask. Kind: needs dev (touches `python/evoca_explore.py`, so outside this workstream's write scope; N to authorize). Cross-ref G2, F1.
- **ACT-N3**: H1, the adaptive-vs-drift program: lineage net-return (D2) plus `A_new`/z-tail plus a shadow-beating test on the section-B sweep, including the adaptive fraction `excess/raw` (B2). Kind: analysis. Unblocked by the 2026-10-15 merges of D2 and G1; coordinate with Lineage-co-evolution (D2/D3) and Resource-driving (B1/B2, L1). Cross-ref checkpoint H1.
- **ACT-N4**: Repair or replace `A_new` and `eff_activity` per ACT-Q1/ACT-Q2 (keying both sides consistently; a coarsening that actually collapses silent churn). Kind: needs dev (C and python). Blocker: ACT-D1 (only worth doing if the scalar route is chosen). Cross-ref G1, F1.
- **ACT-N5**: Long-run (>= 50k ticks) at a pure-evo winner to characterize absolute (not just differential) behavior of the excess metrics; listed as step 4 of `Scans/2026-05-17_pure_evo_analysis.md`. Kind: campaign; N=128 is enough to start.
- **ACT-N6**: Calibration study for the z-score aggregate (ACT-Q9) and the diversity guards (ACT-Q3). Kind: analysis, then dev. Needed before z-tail counts are used as an objective.
- **ACT-N7**: Reaction-diffusion cross-check: the 04-27 balanced and evo-focused winners were picked on the legacy composite; re-rank them under the D1 objective (the RD-robustness #4 inputs came from those). Kind: analysis; coordinate with Reaction-diffusion.

#### Parked
- Per-gene open-space (LUT-hash) neutral shadow with zero-fill (checkpoint G3): decided against (changes the estimand; pooled shadow is the correct estimator there).
- evo_focus "Scan 5" (push `mu_lut` to 0.10-0.20, `m_scale` to 2.0+) and N=512 validation of the 04-27 top configs: obsolete under the corrected-metric finding that the 04-27 optimum was a metric artifact.
- Per-config reproducibility fix for pooled workers: done (`set_seed`, S1); the 2026-04-27 scans remain non-reproducible per config, only in distribution.
- Recipe `Runs/2026-04-27_high_excess_activity_slope.evoca` is a copy of `push_edges` cfg 98 (descriptor `top_9_cfg98_score0.78`), not a high-excess-slope winner; rename or re-describe when convenient (low priority).

### Resource-driving

Source: [`Resource-driving/ISSUES.md`](Resource-driving/ISSUES.md)

_Updated 2026-10-15_

#### Decisions needed (N)
- **RES-D1** (resolved 2026-10-15): N authorized the Tier-2 merge; `resource-driving` (C1 `food_per`) is on `main` after bit-for-bit parity at `food_per=1`. Record: `Docs/Oct2026/01-branch-review.md`, Tier-2 section.
- **RES-D2**: Evolution arms for the RES-N1 campaign — options: (a) LUT-only (`mu_egene = mu_egenome = 0`), where `A_new` is well posed; (b) joint LUT + egene as in R1 and the checkpoint's two "exception" hypotheses, where `A_new` is not well posed (RES-Q3) and only the fixed-space `eg_excess` / `dyn_excess` and eff-activity are interpretable; (c) both as separate arms. Recommendation: (c), with (a) as the primary inverted-U figure. Blocks: RES-N1 design freeze.

#### Open questions
- **RES-Q1**: Is the high-resource flank a drift regime (raw diversity high, excess near zero), as checkpoint section B predicts? R1 never sampled `food_inc` above 0.018, so the high half of the inverted U is untested at campaign scale. The 2026-10-15 pilot (1 seed, `tax` 0.035, LUT-only, N=128, 3000 steps; notebook final section) shows `A_new` rising, peaking in the interior and collapsing at `food_inc/tax` about 2.9 while raw distinct-genome count saturates — suggestive only, and confounded by occupancy (below ratio 1 more cells means more genomes; judge excess per capita or against a population-matched baseline). Source: `Docs/checkpoint 3 July 2026.md` section B; notebook pilot.
- **RES-Q2**: Does per-capita throughput Phi carry information beyond `food_inc/tax`? In the pilot, Phi sits at about `tax` whenever the population is food-limited and nothing overflows, and rises toward `food_inc` only on the overflow-dominated rich flank, so as an x-axis it is nearly degenerate across the viable interior. This follows from the L1 energy identity: at quasi-steady state Phi = max(`tax`, `food_inc`), and below ratio 1 the final alive density equals `food_inc/tax` (0.403, 0.600, 0.897 vs 0.4, 0.6, 0.9 in the pilot), so the ratio is also the carrying occupancy. Test: Phi vs the ratio across `tax` levels in RES-N1; also try total throughput per cell and Phi restricted to cells with positive mouthful.
- **RES-Q3**: `A_new` (G1) is built on the LUT-hash Channon shadow, which models LUT bytes only, while the activity hash covers the full genome. Per the recorded `shadow-models-lut-only-mismatch` finding (and `Docs/EvoCA_research_directions.md` section 2), excess metrics are not well posed when egenes evolve. Is `A_new` interpretable in the joint LUT + egene arm? Needs either an egene-aware shadow for `A_new`, or restriction to LUT-only arms. Do not interpret `A_new` from R1-style joint runs.
- **RES-Q4**: Does the viability threshold and the productive peak collapse onto `food_inc/tax`? Weak hint in existing R1 rows: ordering extinction by the ratio has AUC 0.863 vs 0.839 by `food_inc` alone, and at `food_inc` 0.013 extinction is 0% at `tax` 0.030/0.035 but 30% at 0.040 (`m_scale` and mutation axes pooled, so confounded). R1's `tax` range (0.030 to 0.040, a 33% span) is far too narrow for a real collapse test. Source: `Research/Resource-driving/Scans/2026-05-18_R1_mscale_bracket/results.csv`, notebook R1 section.
- **RES-Q5**: R1 viability was U-shaped in `m_scale` (worst at intermediate under scarcity). Since intake is `(m/25) * matches * F`, is the real control axis closer to `food_inc * m_scale / tax` than `food_inc / tax`? This is a hypothesis from the intake formula, not a finding; RES-N1 should cross `m_scale` with the ratio to test it.
- **RES-Q6**: Are N's two "exceptions" real: does egenome (cognitive load) or LUT eating-speed evolution persist at high resource (checkpoint prompt; section B)? The falsifying outcome (cognitive load still rising at high `food_inc`) would be the more interesting result. Readouts: `cog_load`, per-lineage net metabolic return (D2, merged).
- **RES-Q7**: R1's "low `mu_lut` wins" is conditional on `m_scale` <= 2.5 (optimum shifts high when `m_scale` is bracketed up). Does the mutation optimum also depend on the ratio? Fold into RES-N1 as a `mu_lut` axis at two levels.

#### Next directions
- **RES-N1**: Full inverted-U campaign (the section-B figure, checkpoint B1/LB2), unblocked by the 2026-10-15 merges of L1 and G1 — campaign, no blocker (uses `run_sim` as is). Spec: x-axis `food_inc/tax`, log-spaced about 10 points over 0.25 to 4 (denser at 0.25 to 0.8 where the pilot shows a sharp viability threshold, and extending well past R1's 0.018); `tax` at three levels spanning at least a factor of 5 (for example 0.01, 0.035, 0.07) so the collapse is a real test; `m_scale` at {2.5, 3.5}; `mu_lut` at {0.03, 0.15}; three seeds; N=256, 8000 steps (pilot at N=128/3000 steps shows the activity metrics still growing); arms per RES-D2. Record: `extinct`, `final_pop`, `a_new_final/slope`, `dyn_excess_pc_slope`, `eg_excess_pc_slope`, `n_distinct_genomes_mean`, `cog_load_mean`, `diss_phi_mean/_final`, `diss_tax/death/overflow_step_mean`. Size: about 10 x 3 x 2 x 2 x 3 = 360 runs per arm (R1: 300 runs in 743 s on 32 cores). Falsification criteria: (i) curves at different `tax` superimpose against the ratio but not against `food_inc`; (ii) raw diversity flat/rising on the high flank while `A_new` and `dyn_excess` fall; (iii) if `A_new` or `cog_load` keeps rising at high ratio, the flat-landscape claim is falsified. Cross-ref: B1, B2, L1, G1, C-none.
- **RES-N2**: Promote `digests/resource-driving-inverted-U.ipynb` from skeleton to a result digest after RES-N1: raw vs excess overlay against the ratio and against Phi, 0.013 window marked, markers encode other scan params (LB2); add the adaptive-fraction readout `excess/raw` (B2). Kind: analysis; blocked on RES-N1.
- **RES-N3**: Adaptive-vs-drift test on the RES-N1 sweep (H1): lineage net return (D2, merged) plus `A_new` plus shadow-beating test at starved/intermediate/rich ratios. Kind: analysis; cross-workstream with Lineage-co-evolution and Activity-metrics; blocked on RES-N1.
- **RES-N4**: Add Phi and the dissipation breakdown to the interactive `ts` probe (it currently logs `t,pop,F_env,f_priv,lut_div,eg_ent,activity_flux`; no dissipation column) so a run can be watched crossing from the starved to the rich flank in the SDL/notebook UI. Kind: needs dev (small, python/controls.py plus C getter already exists); no blocker.
- **RES-N5**: C3 campaign: does temporal lumpiness of the food field create evolutionary opportunity? With `food_per` > 1 (C1), hold mean `F_env` fixed while raising its temporal/spatial variance and test whether excess activity rises at the *rich* ratio, where uniform regeneration gives drift. Prediction and negative-result reading in checkpoint section C. Analysis caveat: with `food_per` > 1 the per-step dissipation scalars `run_sim` samples every `sample_every` steps (`diss_*_step_*`, `diss_phi_*`) alias with the regen period, so use the `*_cum_final` totals or pick `sample_every` deliberately relative to `food_per`. Kind: campaign; unblocked (C1 merged 2026-10-15).
- **RES-N6**: C2/L2: float `env_mask` ceiling plus relaxation toward the mask (LB3 option (a)), and a spatial source/sink gradient; test whether a more efficient lineage extends the viable distance from the source (checkpoint section L item 2). Kind: needs dev (new branch from `main`; the `resource-driving` branch was merged and deleted); no blocker.
- **RES-N7**: Localize the viability threshold `ratio_crit(tax)` by bisection at several `tax` levels (the pilot shows extinction-to-viable between ratio 0.3 and 0.4 at `tax` 0.035; near-threshold runs are highly sensitive to the random stream: in the pilot the same seed gave 21 cells inside `run_sim` and 255 when the identical configuration was stepped by hand, because the per-sample activity/shadow updates consume the C PRNG; use replicate seeds and report the extinction fraction). It is the cheapest direct test of the ratio hypothesis (RES-Q4). Kind: campaign (small); no blocker.

#### Parked
- Food + waste model (checkpoint LB6 / E2: waste with its own `W_env` field as a mechanism that keeps selection alive at abundance): belongs to chem milestone 2 (E2, needs a reaction engine; chem milestone 1 merged 2026-10-15 as N_chem=1 scaffolding only), not this workstream.
- Motility / navigating a gradient (M1): deferred by N on 2026-07-19; characterize sessile EvoCA first. Not re-pitched here.
- Per-agent evolvable `tax` (A2, mu-genome): would let harsh environments evolve low-tax agents (N's section A response); infrastructure (A2 storage scaffold merged 2026-10-15; the gene is not yet read or mutated), tracked on the board, not this workstream.

### Lineage-co-evolution

Source: [`Lineage-co-evolution/ISSUES.md`](Lineage-co-evolution/ISSUES.md)

_Updated 2026-10-15_

The workstream question has two parts: (1) is there adaptive LUT<->egene *coupling* (co-evolution) or just joint drift, and (2) can lineage-level net metabolic return separate adaptation from drift. H1/D3 address part 2; the #8 redesign is the only direct handle on part 1. Example-run evidence cited below comes from `Lineage-co-evolution.ipynb` (single seed, N<=128, not campaign results).

#### Decisions needed (N)
- **LIN-D1**: #8 redesign: which fix option, now that D2 exists? (`Docs/plan-8-lineage-coevolution.md` lists A-F; the board lists old #8 as "superseded by D2 + H1", but D2/H1 answer part 2, not the part-1 coupling question.) Options:
  - *Retire #8*, rely on H1 only. Cheapest; leaves "do LUT and egene changes get retained *jointly*?" untested.
  - **Recommended: Option A (Hamming magnitude of LUT bytes vs egene bytes per surviving edge) with `restricted_mu=1` (the existing toggle, i.e. Option F for free), reusing the harness.** Notebook section 3 shows the binary indicator saturating exactly as 1-exp(-250 mu) (0.918 at 0.01, 0.9995 at 0.03) while Hamming keeps mean 250*mu and nonzero spread, and bulk LUT bytes come from one pointer read, so no per-cell loop. Pure Python, no C change. Add the D2 join as phase 2 (below).
  - Option B (add lut_complexity level change) as a secondary read-out; needs a tiny per-cell accessor (dev).
  - Options C, D, E: heavier; keep in reserve if A shows signal but non-monotone coupling (C) or needs temporal persistence (D).
  - Caveats found in the example run: (a) Spearman(LUT Hamming, egene Hamming) was ~0 at every `mu_lut` there, which is the expected null *and* what unpowered silent-bit noise looks like, so run a power check (inject a synthetic coupling) before the campaign; (b) `restricted_mu=1` changes the Poisson rate to `mu_lut * n_active`, so `mu_lut` must be re-chosen; not tested here.
  - Blocks: any #8 campaign. Does not block H1/D3.
- **LIN-D2**: Lineage unit for H1 when `mu_lut` is productive. D2 keys on the full-genome hash, i.e. a *genotype*. Example run: at `mu_lut` >= 0.01 almost every live cell is a unique genotype (fraction sharing a genotype 0.01 at 0.01, 0.0 at 0.03; no hash reaches 10 births at 0.03), so per-hash net return is a per-individual quantity, not a clade property. Options: (a) restrict H1 to `mu_lut` <= ~0.003 where genotypes persist (no dev, but excludes the productive corner used by #3c/#8); (b) clade roll-up in Python along `parent_id` chains from the existing getters (analysis, moderate cost); (c) a C-side clade reducer (dev, only if (b) is too slow). Recommendation: (a) for a first H1 pass plus (b); decide (c) only if needed. Blocks: H1 at the productive corner of the Sec. B sweep.
- **LIN-D3**: Tax basis in the net-return reducer. It counts nominal tax; the dissipation accountant counts tax actually removed. In the example runs the reducer tax was 20-23% above the accountant's (887,729 vs 738,159 evolving; 939,192 vs 763,713 control, N=96, 3000 steps), so per-lineage `net` is biased low, more for short-lived lineages. Options: document and correct in analysis (cheap, loses per-lineage resolution), or add an actual-tax column to the reducer (small C change on merged code, needs N's gate). Recommendation: add the column before H1 conclusions rely on `net`. Blocks: trusting absolute `net` in H1/D3 (relative comparisons within one run are less affected).

#### Open questions
- **LIN-Q1**: Is there any adaptive LUT<->egene co-evolution at all? Recorded evidence: #3 3a >> 3b but confounded; #3c says the C-vs-A residual is small (22/30 positive, mean +0.6) and the egene signal weak and plateaued by 7-8k ticks (`Scans/2026-05-19_3c_coevolution_substrate/notes.md`). The direct test (#8) has produced no data (`Scans/2026-05-19_8_lineage_coevolution/` has `scan.py` only). Unanswered.
- **LIN-Q2**: Is the #3c plateau real or a saturation of the cumulative-excess metric? In the 4-seed N=96 replica, C's cumulative `eg_excess_pc` rises and falls with population boom-bust cycles (seed 1: ~250 to ~7000 cells, period ~800 steps). The shadow is not re-baselined at the phase switch (follow-up F1 in `research_board_archival.md`), so absolute excess is biased. Needs F1 or a lineage-based read-out to settle.
- **LIN-Q3**: Raw net return is pinned in a saturated, resource-limited lattice. Example run: total intake equals the food source (~414,000 vs `food_inc*N^2*T` = 414,720) in both arms, and the control's dominant genotype has constant per-capita net -0.020 = `food_inc*N^2/pop - tax`. So what is the right *relative* fitness quantity for H1: net return relative to the same-time population mean, per-lineage growth rate vs shadow (checkpoint Sec. H step 2), or births per capita? D2 output alone does not answer adaptive-vs-drift.
- **LIN-Q4**: Frozen-GoL (3b/B) viability is not reproducible across N, seed, or the RNG stream. Example: config 253 survives at N=128 seed 0 without the neutral shadow and goes extinct with it (seeds 1, 2 extinct either way); the recorded N=256 run survived. Does enabling the shadow perturb the real run's RNG stream (the shadow draws random kills, `C/evoca.c` ~line 2418; not verified)? If so, shadow-on and shadow-off runs are different experiments, which matters for comparing shadow-based campaigns with the shadow-free #8/H1.
- **LIN-Q5**: How often does the `dyn_excess_pc` artifact fire? Recorded: 21 of 30 A runs have slope > 50 (median +249; larger at lower phase-1 `mu_lut`). The N=96 replica showed it in 1 of 4. Understand the N/length dependence before using any dyn-based metric with `mu_lut=0` (rule: only interpret `dyn_excess_pc` where `mu_lut>0`).
- **LIN-Q6**: Provenance gap in campaign #3: `Scans/2026-05-17_pure_evo_egene_on_winners/` has `scan.py` and `results.csv` but no `scan_config.json` and no `notes.md`, so `evoca_from_scan` cannot rebuild its configs. Needs a notes file and a `save_scan_config`-style record (N=256, 8000 steps, halfplane init, GoL LUT; values read from `results.csv`).

#### Next directions
- **LIN-N1**: D3 punctuation detector: lineage net-return dip-then-staircase at Negene increments, single unforced run. Kind: needs dev (analysis). Unblocked by the D2 merge. Gaps found: (i) a Negene change creates a *new* hash row, so the detector must link child hash to parent hash; `lineage_parent_hash` gives the parent's hash per newborn but no per-cell full-genome hash accessor exists, and no parent-hash link table is exposed; (ii) `last_step` is only as current as the last `lineage_return_update()`, so trajectories need periodic sampling (longevity resolution = sampling interval); (iii) per LIN-Q3, detect the step in a *relative* quantity. Checkpoint action D3; N's 2026-07-19 note on Sec. D proposes punctuation as a random-slot-addition event (slot activation alone may have negligible effect).
- **LIN-N2**: H1 adaptive-vs-drift program on the Sec. B `food_inc/tax` sweep: lineage net return + A_new + shadow-beating lineage count, with the falsifier (raw activity inverted-U too means the drift-collapse story is wrong). Kind: campaign + analysis. Blockers: G1 (A_new) merged; LIN-D2 (lineage unit), LIN-D3 (tax basis), LIN-Q3 (relative fitness). Checkpoint action H1 (shared with Activity-metrics).
- **LIN-N3**: #8 relaunch per the LIN-D1 outcome: synthetic-coupling power check first, then 2 arms x 8 seeds as in the original spec with `mu_lut` retuned if `restricted_mu=1`. Kind: campaign (pure Python harness exists, flawed observable only). Blocked on LIN-D1.
- **LIN-N4**: Join D2 to the #8 edges: use the child genotype's reducer `births` (reproductive success, not `net`, given LIN-Q3) as the outcome and test the dLUT x dEgene interaction beyond additive effects, with a permutation null. Kind: needs dev (analysis). Blocked on LIN-D1, LIN-D2.
- **LIN-N5**: Re-baseline the shadow at a phase switch (F1) and rerun a #3c-style A/C comparison with lineage read-outs instead of `dyn_excess_pc`. Kind: needs dev (small C change) then campaign. Not started; N to decide priority.
- **LIN-N6**: Finish the `coevolution-substrate` digest (currently a draft in `digests/`) with the section 2 findings, and write `notes.md`/`scan_config.json` for campaign #3 (LIN-Q6). Kind: doc/analysis, no blocker.

#### Parked
- Option C/D/E variants of the #8 observable: revisit only if LIN-D1's Option A shows signal that is non-monotone or temporal.
- Evolving-rate genes (`mu-genome`, A2) as a co-evolution lever: A2 storage scaffold merged 2026-10-15 (`get/set_mu_lut_gene`; inherited, not yet mutated or read); the evolvable-rate follow-on is infrastructure work, outside this workstream until it lands.
- Recipes: none exported; the notebook's example runs are scripted (seeded) rather than saved as `.evoca`, because several need mid-run switches (#3c) or post-init `enable_lineage()` that the recipe format does not appear to capture (not checked in detail).

### Reaction-diffusion

Source: [`Reaction-diffusion/ISSUES.md`](Reaction-diffusion/ISSUES.md)

_Updated 2026-10-21_

#### Decisions needed (N)
- **RXD-D1**: What is the second objective for RD selection? — Options: (a) evolutionary activity as before (`n_distinct_genomes_temporal_std` + `unique_top_genomes`, which saturates at the sample cap, see `Scans/2026-04-27_scan_analysis.md` obs. 8); (b) dyn_excess / A_new on eff_activity species (Activity-metrics, merged 2026-10-15); (c) lineage net-return (D2) or per-capita dissipation Phi (L1) as a more adaptive axis (`Docs/Oct2026/00-review-of-prompt1.md` §7). Recommendation: run RXD-N1 first with (a), then repeat with (b); defer (c). Blocks the final form of RXD-N1.
- **RXD-D2**: Is the halfplane start the RD protocol, or should RD be scored from a uniform random start (`init='fraction'`) and/or after a burn-in? — The halfplane start makes the first ~N steps an invasion front, so every 04-27 spatial metric partly measures that transient (observed in this notebook's N=128 runs, not a campaign result). Recommendation: keep halfplane for continuity but add a burn-in variant (RXD-N2). Blocks any re-run of the scans.
- **RXD-D4**: Should legacy recipes (written before 2026-10-15, no `food_cap` key) load with `food_cap=0` by default? — They ran with the uncapped mouthful, but `import_run` / `load_recipe` give them the new default `food_cap=1` and print a notice only. Scan rows already get this correction at load time (`evoca_explore.correct_scan_row_params`); recipes do not. This matters for RD: cfg246's waves (alive 0.11-0.72) vanish under the cap (0.41-0.48), and the boom-bust waves are damped out by evolution by ~t=5000 under the cap but persist to t=10000 without it (notebook section "Where the big RD waves are", N=256 example runs). The gallery and decay figures in the notebook currently run the 04-27 recipes with `food_cap=1`, so their "ex" vs "scan" comparisons mix two eating rules. Recommendation: apply the same load-time correction in `import_run` (a missing `food_cap` key means 0), then re-run the gallery. Quick alternative: pass `food_cap=0` in the notebook's `run_track`. Blocks a faithful gallery.
- **RXD-D3**: Keep the 85 recipes in `Runs/` as is, or dedupe? — About 55 distinct configs; `top_N`, `top_evo_N`, `auto_*` and three `best_corr_length` recipes repeat the same metaparameters under different ranking labels. Recommendation: keep (names record which ranking picked a config; disk cost trivial), add a note in the README only. Low urgency.

#### Open questions
- **RXD-Q1**: Is the correlation-length plateau (~75-80, about N/3) a grid-size ceiling? — Recorded as a hypothesis in `Scans/2026-04-27_scan_analysis.md` (recommended step 2); never tested at N=384/512.
- **RXD-Q2**: Does the #4 decay verdict hold for the higher-scoring push_edges winners? — The campaign's ten bases are all `initial`-scan recipes (`top_N_cfg*_score0.55-0.75`); the push_edges and evo_focus winners (corr_len 60-80, the actual "best" spatial configs) were not in `Scans/2026-05-18_rd_under_pure_evo/scan.py::load_bases`'s output. Cheap to extend with the existing harness.
- **RXD-Q3**: Which of the visual regimes are attractors versus transients? — In this notebook's N=128 example runs the dot-lattice and salt-and-pepper states are end states, the invasion-front and sparse-colony regimes are transient; not established at N=256/long horizon. Needs a classification metric beyond correlation length (the labyrinth, cfg53, has corr_len ~1 yet clear texture).
- **RXD-Q4**: Is there a persistent regime at `gdiff` far above the scanned range? — `fill_dead_region_with_evolving_patches` uses gdiff=20 (frozen domain mosaic); every scan grid stops at gdiff 0.20. Unscanned.
- **RXD-Q5** (answered 2026-10-21): Where are the big sustained RD waves? — Not in `waves_in_blank_half`. They are in the push_edges corner (m_scale >= 1, food_inc ~0.013, gdiff ~0.06, tax ~0.035) with `food_cap=0`. Two recipes now hold them: `2026-10-21_rd_waves_boom_bust.evoca` (the "Scan-3 productive corner" config from `test.ipynb`; synchronized lattice-wide bloom and crash, period ~100-130 steps) and `2026-10-21_rd_waves_patchy_cfg246.evoca` (push_edges cfg246; out-of-phase patches). Also wave-like: `top_16_cfg102` (either cap) and `Criticality-bifurcation/Runs/2026-10-15_foodinc013_osc.evoca` (either cap). History of the old question: `waves_in_blank_half`'s saved *initial* metaparams go extinct by ~t=200 (N=128/256/512, seeds 0-1, example runs); with the saved *final* metaparams it is a one-shot invasion front, not a sustained wave. Two-parameter descriptor ambiguity in recipes (`metaparams_init` vs `metaparams_final`) is worth a note on any recipe promoted to a gallery.
- **RXD-Q6** (partly answered 2026-10-21): Does EvoCA make spiral waves? — Yes, under the legacy uncapped mouthful. One saved recipe does: `Resource-driving/Runs/2026-10-15_R1cfg208_intermediate_fi0p016.evoca` from its final metaparams (`m_scale` 3.5, `food_inc` 0.016, `tax` 0.035, `gdiff` 0.08, `food_cap=0` recorded), best seen at N=512 (cramped into arcs at N=256). Found by running all 103 saved recipes (both metaparam sets) and inspecting env-food frames. Notebook section "Spiral waves" (N=512 example runs, mostly seed 0): spirals in a narrow window, `m_scale` 2.5-3.5 with scarce food just above extinction (none at `m_scale` 1.2 or 5.0, none at `food_inc` >= 0.020); they need `gdiff` ~0.08 (extinct at <= 0.02, uniform at >= 0.3); 3 seeds all show them, lattice-filling only for seed 0. With `food_cap=1` there are none anywhere on the same 5 x 5 grid. Under `food_cap=0`, 7.2% of the lattice has `F` < 0 (min -0.61) at t=1600, so a post-front food deficit acting as the refractory zone is the working hypothesis, not demonstrated. The RD scans never reached them because they stop at `m_scale` 1.2. Open: whether the corrected (capped) model makes spirals anywhere (RXD-N6).

#### Next directions
- **RXD-N1**: Re-select the 04-27 winners with conepareto box and wedge and compare with the composite top-20 — analysis (no new dev; package exists at `~/Projects/conepareto`). Concrete first experiment: pool the candidates of `initial`/`refined`/`push_edges`; spatial score = mean of normalized `correlation_length_mean` and `largest_patch_temporal_std`, evolution score = mean of normalized `n_distinct_genomes_temporal_std` and `unique_top_genomes` (a preview with these choices is in the notebook, section "Two objectives"); form top-20 by (i) composite mean, (ii) box with 0.5/0.5 exclusion floors, (iii) wedge with `max_imbalance=0.15`; report overlap and Pareto layers; then run each distinct set (3 seeds, N=256, 12000 steps, the #4 harness) and score by persistence (`correlation_length_final/mean`, `largest_patch_final/mean`) plus final `dyn_excess_pc`: the selection method that keeps structure with turnover wins. Cross-ref: H1 (adaptive vs drift) for the evolution axis.
- **RXD-N2**: Burn-in / persistence-based spatial metrics — analysis, needs a small harness change (`python/evoca_explore.py` is outside this workstream; request via dev). Report spatial metrics after a burn-in (e.g. last half of the run) and final/mean as first-class scan columns.
- **RXD-N3**: Extend campaign #4 to the push_edges/evo_focus winners (RXD-Q2) and to N=384/512 (RXD-Q1) — campaign, existing harness, torque-scale compute (30 workers).
- **RXD-N4**: Scan gdiff in [0.2, 20] together with m_scale (RXD-Q4) — campaign.
- **RXD-N5**: Add a visual-regime classifier or texture metric (e.g., structure-factor peak, anisotropy) so regime names in the notebook gallery become measured labels — needs dev (analysis helper), cross-ref Criticality-bifurcation for order-parameter ideas.
- **RXD-N6**: Spirals in the corrected model (RXD-Q6) — campaign (small) plus analysis. The spirals found so far need the legacy overdraw (`food_cap=0`). Under `food_cap=1`, search for a refractory zone that does not depend on overdraw: lower `food_inc` with higher `m_scale` near extinction, `gdiff` 0.04-0.2, and periodic regrowth (`food_per` > 1, on `main` since 2026-10-15), at N=512, with a spiral detector (phase singularities of the local food oscillation, or the free-end count of thresholded fronts) instead of inspection by eye. Also test the deficit hypothesis directly under `food_cap=0` by clamping `F` at 0 after eating while keeping the uncapped intake (needs a C option; dev). Cross-ref: Criticality-bifurcation (`foodinc013_osc` makes target waves, not spirals).

#### Fixed
- **2026-10-21, `mu_egenome` read with its post-rename meaning.** The 2026-04-27 scans (initial, refined, push_edges, and Activity-metrics' evo_focus) record the egene mutation rate as `param_mu_egenome`, its name before the 2026-04-28 rename to `mu_egene`; the name `mu_egenome` now means the Negene slot-growth rate. `evoca_from_scan` / `params_from_row` / `import_run` passed it through under the new meaning: mu_egene=0 plus nonzero slot growth. This affected the 20 `2026-04-29_top_evo_*` recipes and both `best_corr_length` recipes (the files carry the old key). No campaign consumed them. Fix: one load-time correction, `evoca_py.correct_legacy_mu_egenome` (a dict with `mu_egenome` but no `mu_egene` is pre-rename), used by `import_run` and `evoca_explore.correct_scan_row_params`; tests in `tests/test_legacy_mu_egenome.py`, which fail without the fix. Found because push_edges cfg246's waves vanished when re-exported from its scan row.

#### Parked
- Larger-N interactive visual validation of scan winners in the SDL viewer (recommended in the 04-27 analysis, step 1); superseded by the headless N=256 gallery in the notebook.
- Fixing `evoca_set_seed` per task in pool workers is done (`build_sim` reseeds); the 04-27 note about shared RNG state predates it.

### Criticality-bifurcation

Source: [`Criticality-bifurcation/ISSUES.md`](Criticality-bifurcation/ISSUES.md)

_Updated 2026-10-15_

#### Decisions needed (N)
- **CRIT-D1**: First control parameter, range and N for the torque campaign — options: (a) `food_inc` 0.010-0.06 fine grid (the free-energy source; Phi is its natural partner; the first example diagram already shows extinction, oscillation, amplitude collapse and saturation inside this range), (b) `m_scale` 1.5-6 (R1: interior productive optimum ~2.5-3.5, U-shaped viability), (c) both as a 2-D map. Recommend (a) first at fixed R1-informed values (m_scale 3.5, tax 0.035, mu_lut 0.06, gdiff 0.06), then (b) at food_inc 0.018. Blocks CRIT-N1.
- **CRIT-D2**: What counts as a transition (a pre-registered criterion, so a "flat diagram" is a valid outcome) — recommend all of: extinction boundary (fraction of replicate seeds extinct crosses 0.5), onset of sustained oscillation (spectral peak above a noise floor, not the raw extrema count), amplitude collapse (tail std/mean falling through a threshold), and, for any candidate, an N-scaling test. Blocks CRIT-N1/N4.

#### Open questions
- **CRIT-Q1**: Is the observed boom-bust oscillation a finite-size synchronized-lattice cycle (relative amplitude decays with N) or intrinsic? The single-seed example check (notebook) suggests the large-amplitude cycle at food_inc 0.016 persists with N (std/N^2 0.335, 0.318, 0.297 at N=32/64/128) while the small-amplitude regime at 0.022 shrinks (relative amplitude 0.033 -> 0.009 from N=64 to 128); neither is replicated. Decides whether diagram columns mean anything in the thermodynamic limit; couples to the Docs/EvoCA_research_directions.md §5 H_critical / H_ceiling question (xi ~ N/3 plateau).
- **CRIT-Q2**: Is `recurrent_extrema` the right reducer for noisy, irregular-amplitude series? In the example runs the extrema count is a filled band (tens to hundreds of points), so it cannot distinguish limit cycle from irregular dynamics. Candidates: autocorrelation/spectral peak, return map of successive maxima, Poincare-section sampling, largest-Lyapunov-style divergence of nearby seeds. See `python/bifurcation.py` `recurrent_extrema`.
- **CRIT-Q3**: Which side of any boundary is a bifurcation in the dynamical sense versus a stochastic extinction (absorbing-state) transition? EvoCA has an absorbing state (pop=0), so the extinction edge may belong to the directed-percolation-like class rather than a deterministic bifurcation; the R1 extinction table (m_scale x food_inc) is endpoint data only.
- **CRIT-Q4**: Hysteresis: does ramping `food_inc` up versus down give different columns (the directions doc, R3 (c), suggests a slow ramp as an alternative to a fixed grid)? Untested; needs a ramp driver (CRIT-N3).
- **CRIT-Q5**: Does evolution matter? Runs with `mu_lut=0` freeze the substrate (see memory notes on frozen rich substrates); compare mu_lut=0 against mu_lut>0 columns to separate ecological (food-limited) oscillation from evolutionary change.

#### Next directions
- **CRIT-N1**: Real `food_inc` bifurcation campaign — campaign (torque), N in {64, 128, 256}, 8000+ ticks, 41+ values x >=5 seeds, observable `pop` first; blocked on CRIT-D1/D2. Uses `python/bifurcation.py` `sweep` unchanged.
- **CRIT-N2**: Multi-observable harness — needs dev (small, Python only): extend `scalar_series` to record Phi (`get_dissipation`), births/deaths, `dyn_activity` total, `A_new` (needs neutral shadow, computed as in `evoca_explore.run_sim`), `lut_div`, lineage return; the example notebook has a working per-tick `record()` prototype. Cross-ref checkpoint actions G1 (A_new), L1 (Phi), D2/D3 (lineage return).
- **CRIT-N3**: Slow-ramp / hysteresis driver — needs dev (Python): change `food_inc` (or `m_scale`) during a run through the live sim attribute instead of re-initializing, up and down.
- **CRIT-N4**: Control-parameter and finite-size follow-up on whichever transition passes CRIT-D2 — analysis/campaign: order parameter = alive density or oscillation amplitude, susceptibility, Binder cumulant across N; fixes the Q1 H_critical / H_ceiling question (directions doc §5, R3 (a)/(b)).
- **CRIT-N5**: Second axes — campaign: `m_scale` (R1 optimum), `tax`, `mu_lut`, then a 2-D map (m_scale x food_inc) where the R1 non-monotone viability suggests structure; start from R1 `results.csv` to pre-select regions.
- **CRIT-N6**: Replace the reducer's raw extrema count by a classifier (fixed point / periodic / irregular / extinct / saturated) for sweep summaries — needs dev; feeds CRIT-Q2. Keep `feigenbaum_delta` but do not report it unless a clean cascade appears.

#### Parked
- Feigenbaum delta extraction: period-doubling has not been observed in EvoCA (N's position); helper retained, not used.
- Motility-related parameters: motility branch on hold (memory, 2026-07-19).

### Ring-complexity

Source: [`Ring-complexity/ISSUES.md`](Ring-complexity/ISSUES.md)

_Updated 2026-10-15_

#### Decisions needed (N)
- **RING-D1**: Commission the Track-1 ring-class dominance campaign now? — Options: (a) yes, as an analysis-only campaign on main (no C change; Python-side level-k LUT constructors, LUT hash as the A/B tag at `mu_lut=0`); (b) first build the harness pieces (RING-N2) and the viability screen as a short pilot, then commission. Recommendation: (b) then (a). The notebook pilot (`Ring-complexity.ipynb`, 2026-05-19 section) shows the plan is runnable with no C change; the null spread is modest (F sd about 0.09) but untaxed class effects are of similar size, so the design needs many LUT draws per class. Blocks RING-N1. Source: `Docs/EvoCA_research_stability.md` sections C.5, C.10, F.
- **RING-D2**: Tax schedule for the ring ladder — keep `tax_ring * level` (linear, current, exactly tested in `tests/test_ring_tax.py`), or add a nonlinear / effective-dependence variant? — Linear makes the cost differential between classes a built-in selection pressure, so a "level 1 wins under tax" result is partly by construction. Recommendation: keep linear for the first A/B, and use `eff_activity`-based used-entry counts (merged F1) as the sharper alternative later (research directions section 8). Blocks nothing immediately; shapes RING-N3.
- **RING-D3**: Which "climb" metric defines success for the ring-ladder A/B (#5)? — Options: population fraction at level 3 at the end; time to level-3 majority; level-weighted mean; or net fraction relative to a neutral-mutation baseline (RING-N4). No spec exists (board lists #5 as queued with no headline). Recommendation: report the fraction relative to the neutral baseline. Blocks RING-N3.

#### Open questions
- **RING-Q1**: How much of "ring level climbs under mutation" is mutational entropy rather than selection? — In the example runs (notebook, ring-ladder section; N=96, one seed, not a campaign result) a GoL or level-1 start with `mu_lut=0.002` and no tax goes to >90% level 3 within ~500 steps. Hypothesis (untested): level-1/2 LUTs are a tiny tied subset of LUT space and one bit flip breaks the tie. Needs a neutral-mutation control (RING-N4) before any "evolution climbs" claim.
- **RING-Q2**: Is a ring class stable against competition at `mu_lut=0`, controlling for LUT identity? — Pilot (notebook, 2026-05-19 section; N=64, 3 seeds, 8 LUTs per class): untaxed classes coexist with a weak seed-dependent level-1 gain; null control (6 seeds, two random level-3 pools) gives final F sd about 0.09. Class effects of that size are not distinguishable from LUT-draw noise with so few draws; needs many draws per class (`EvoCA_research_stability.md` C.6, G).
- **RING-Q3**: What is the genelife-like "no resource flow" arm in EvoCA? — `food_inc=0` is not frozen: the initial food stock lets the population fill the lattice, then turnover largely stops (example run: 67% level 3 at t=2000 vs 100% at `food_inc=0.01`; one seed, confounded with number of births). With any substantial `tax_ring` and `food_inc<=0.003` the population goes extinct. Need an explicit design of the null arm (matched birth count, or a low-but-viable food_inc) before the #5 A/B (research directions section 8) is interpretable.
- **RING-Q4**: Does the selection survive in a regime where the lattice is not saturated? — In all example runs the lattice fills (alive = N*N) within a few hundred steps, so selection acts through food and private-food tax rather than space. Not yet explored whether a lower-food regime changes the class outcome.
- **RING-Q5**: Declared level versus effective level — does a level-3 LUT whose n3 entries are rarely queried behave as level 2? — The tax and probe use declared level (exact bit-structure); `restricted_mu` / `eff_activity` give used entries. Relevant to whether the tax penalizes reach that is actually used.

#### Next directions
- **RING-N1**: Track-1 ring-class dominance campaign — kind: campaign (analysis-only, main). Preparations: ring 1 vs 2, 2 vs 3, 1 vs 3 (`EvoCA_research_stability.md` "different rings"); geometries band / mosaic; `mu_lut=0`; many LUT draws per class; viability screen (step 0); null control; early vs late windows; frequency `F` plus front position. Blocker: RING-D1 (and the RING-N2 helpers). Related checkpoint actions: none directly (Lineage H1/D3 are independent).
- **RING-N2**: ~~Add `tax_ring` to `_PARAM_KEYS`~~ **done 2026-10-15 (`d15d0ad`)**: run_sim/build_sim/bifurcation now honor `tax_ring` and warn on any dropped param; regression test via run_sim. Remaining (needs dev, small):  a level-k LUT constructor helper (`make_random_lut(n_init=k)` exists) and an optional mosaic/band initializer.
- **RING-N3**: Ring-ladder A/B (old #5) — kind: campaign, blocked on RING-D3, RING-Q3, RING-N2. Axes: `tax_ring` x `food_inc` x initial population (GoL, random n_init=1/2/3), track level fractions via `get_lut_complexity()`. Status: spec is one paragraph (research directions section 8 plus devlog S2b); no scan config exists.
- **RING-N4**: Neutral-mutation baseline for ring level — kind: analysis (no C change possible with the neutral shadow, which is LUT-only per the memory note; note the known LUT-only mismatch before reusing excess-activity). Run `tax_ring=0` with random selection removed (for example `food_inc` high enough that all survive, or compare against the entropy expectation computed analytically: probability that a mutated level-k LUT stays level k). Informs RING-Q1.
- **RING-N5**: Track-2 ring-bounded mutation branch ("bounded <=": a level-2 LUT may collapse to <=2 but never promote to 3) — kind: needs dev, dynamics-altering, own worktree, merge gate, pre-merge test notebook, N reviews in a /tmp worktree. Not started; no open dev branch (board "Dev branches by workstream"). Suggested toggle pattern: parallel to `restricted_mu`; default off must be bit-for-bit no-op. Consider offering a `/code-review ultra` pass at the pluggable-operator moment.
- **RING-N6**: Optional additive `a_b_tag` C field (per-cell uint8, parent-to-child copy, no mutation; parallels lineage fields) for evolutionary A/B preparations at `mu_lut>0` — kind: needs dev, lightweight but merge-gated (`EvoCA_research_stability.md` C.3). Not needed for Track 1.

#### Parked
- Effective-dependence (mutual-information) tax — deferred from S2b (more code; needs an online MI estimate); revisit if RING-D2 moves away from linear.
- Patch-transfer 2x2 (old #7) — archived, superseded by ring-stability.
- Motility and EvoCA-vs-Bugs framing — on hold (see memory "motility deferred"); not part of this workstream.

<!-- END workstream-issues -->

---

## Dev branches by workstream

**Tier-1 merged 2026-10-15** (N-approved: green combined suite + bit-for-bit
default-param parity). Six branches in `5de867d`; `lineage-return` in
`6548106` after fixing a pre-existing `enable_lineage()` call-order
segfault (`fadb242`). Branch refs deleted. Record:
[`../Docs/Oct2026/01-branch-review.md`](../Docs/Oct2026/01-branch-review.md).
**Tier-2 merged 2026-10-15** (N-authorized "if bug free"): `resource-driving`
`c9ecf32`, `mu-genome` `28e28be`, `chem` `702143a` (two chem fixes; see the
record). Branch refs deleted. Sketches below are plans, not commitments.

### Activity-metrics

- `dyn-promote` (D1): **merged** `5de867d` 2026-10-15.
- `metric-anew` (G1): **merged** `5de867d` 2026-10-15.
- `metric-zscore` (G2): **merged** `5de867d` 2026-10-15.
- `eff-activity` (F1): **merged** `5de867d` 2026-10-15.
- `research-drafts` (84fc29b): merge, then the restructure moves the drafts
  from `Research/_drafts/` into `Research/<Workstream>/digests/` (still
  unapproved drafts, flagged as such), then delete the branch.

### Resource-driving

- `metric-flux` (L1): **merged** `5de867d` 2026-10-15.
- `resource-driving` (C1): **merged** `c9ecf32` 2026-10-15. C2/L2 (RES-N6)
  start on a new branch from `main`; C3 (RES-N5) is unblocked.

### Lineage-co-evolution

- `lineage-return` (D2): **merged** `6548106` 2026-10-15; D3 and H1 can now start on `main`.

### Reaction-diffusion, Criticality-bifurcation, Ring-complexity

- No open dev branches. (Ring Track-2 bounded mutation will get its own
  worktree branch when started.)

### Infrastructure

- `env-agent-split` (A1/A3): **merged** `5de867d` 2026-10-15.
- `mu-genome` (A2 storage): **merged** `28e28be` 2026-10-15. The gene is
  inherited but not mutated or read; the distribution-seeding assay and
  evolvable rates are follow-on branches.

### Chem (no workstream dir until results exist)

- `chem` (E1/E3 milestone 1): **merged** `702143a` 2026-10-15 (N_chem=1 only;
  fixes: `stamp_patch` chem-store resync, `chem_energy` clamped to 1).
  Milestone 2 (E2) starts on a new branch from `main`.

### Other branches

- `bigLUT`: archival 5-ring (1407 B/cell) version. Keep, do not merge.
- Deleted 2026-10-15 (N-approved): `softmouth` (merged PR #1/#2; local +
  origin), `pre6todo` (one Mar-16 test-cell commit `d413c18`; local +
  origin), `lut-hashes-accessor` (merged `f8d6d66`), all 11
  `worktree-agent-*` stubs, and the 7 Tier-1 branches after merge.
- `research-drafts`: see Activity-metrics above.

---

## Campaign results → Research/

Authoritative results live in `Research/<campaign>/` digests, indexed by
`Research/README.md`. Current drafts under `Research/_drafts/` awaiting
approval (Notebook link gated on approval per the README drafting
workflow): `neutral-model-methodology` (flagship), `pure-evo-regime`
(#1/#2), `coevolution-substrate` (#3/#3c), `RD-robustness` (#4),
`viability-brackets` (R1) — real figures from `Research/*/Scans/*/results.csv`;
`resource-driving-inverted-U` — skeleton (sweep gated on
metric-flux/metric-anew); `reaction-diffusion-gallery` — NEW, digest to
be drafted. Pre-checkpoint campaign results (#1–#4, 3c, R1,
causal-control) are in `research_board_archival.md`, mirrored by these drafts.

---

## Reconciliation with checkpoint 3 July 2026 (2026-07-22)

The checkpoint reframing (metaparam split, resource-driving, chem,
energy-gradient, activity metrics) is largely orthogonal to the earlier
program and had silently dropped several tracked threads. Triage (N):

| Pre-checkpoint thread | Disposition |
|---|---|
| Bifurcation / criticality (old #6 + FSS / control-param exponents) | **CARRIED** active — drop the "Feigenbaum" label (period-doubling unobserved) |
| Ring program (old #5 + ring-stability two-track) | **CARRIED** active |
| Lineage co-evolution metric (old #8) | **FOLDED** into lineage-return (D2) + §H |
| Patch-transfer 2×2 (old #7) | **ARCHIVED / superseded** by ring-stability (per `EvoCA_research_stability.md`) |
| Reaction-diffusion patterns | **NEW Research topic** (gallery) |

---


---

## Log (newest first)

- **2026-10-15** — **Tier-2 merged.** `resource-driving` (C1 `food_per`),
  `mu-genome` (A2 rate-gene storage) and `chem` (E1/E3 milestone 1) rebased
  onto main and merged (`c9ecf32`, `28e28be`, `702143a`): 88 tests (from 68),
  clean `-Wall`, parity probe BAD 0 / IDENTICAL after each (viable fixture,
  alive 0.34). Two chem bugs fixed with regression tests shown to fail
  without the fix: `stamp_patch` wrote `f_priv` behind the chem store (host
  food leaked into stamped organisms), and `chem_energy != 1` broke the
  energy identity. New params `food_per`, `n_chem` added to `_PARAM_KEYS`
  with run_sim reachability tests. Merge-test notebooks in
  `tests/notebooks/`. Record: `Docs/Oct2026/01-branch-review.md`.
- **2026-10-15** — **Research restructure + workstream seeding.** Top-level
  Runs/Scans eliminated -> `Research/<WS>/{Runs,Scans,digests}` (names kept;
  `python/evoca_paths.py`; drivers find root by marker; drafts in digests/).
  Six lab notebooks + ISSUES.md seeded (Sonnet agents; all execute clean);
  this board's issues section is generated from them (`tools/build_board.py`).
  Seeding surfaced merged-code problems: `tax_ring` silently dropped by
  run_sim (fixed `d15d0ad`); neutral shadow shares the dynamics RNG
  (observer effect, open); eff_activity union mask saturates (open); A_new
  blind in egene-only runs (open). See `Docs/Oct2026/06-seeding-findings.md`.
- **2026-10-15** — **Tier-1 merged + branch cleanup.** Six additive branches
  merged in `5de867d` (dyn-promote, metric-anew, metric-flux, metric-zscore,
  eff-activity, env-agent-split): 50/50 tests, clean `-Wall`, bit-for-bit
  default-param parity with pre-merge main. `lineage-return` held back for a
  segfault, root-caused to a pre-existing `enable_lineage()` call-order trap
  from `lineage-field` (no existing results affected: every caller used
  `init(lineage=True)`), fixed in `fadb242` with regression tests shown to
  fail without the fix, then merged in `6548106` (57/57, parity identical).
  Pushed. Stale/merged branches deleted (see "Other branches"). Record:
  `Docs/Oct2026/01-branch-review.md`.
- **2026-10-15** — **Board restructure.** Per N's Oct-2026 prompt and review
  `../Docs/Oct2026/00-review-of-prompt1.md` §4: moved the pre-checkpoint
  archive (old board, verbatim) to `research_board_archival.md`; added the
  "How to read the codes" box; rewrote "Active workstreams" as per-workstream
  tables annotating every code with plain meaning, kind (built / needs dev /
  campaign) and blocker; adopted the six approved workstream names (Wave 1
  stays the name of the July branch batch); added "Dev branches by workstream";
  updated "Decisions needed" (W1-d resolved by N; W1-a/b/c proceeding on
  review §3e defaults). Tier-1 integration since completed (see the next entry)
  (green combined suite + bit-for-bit default-param parity); Tier-2 stays at
  N's gate.
