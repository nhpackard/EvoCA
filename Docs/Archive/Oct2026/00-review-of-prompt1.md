# Review of `Docs/prompt1-Oct2026.md` (2026-10-15)

Claude's response to N's October prompt: the state of the repo it found, what it
thinks of each request, the conflicts that need N's decision before agents start,
direct answers to the "remind me" questions, and a proposed parallel execution
plan.

**Where responses live (proposal).** Answers go in `Docs/Oct2026/`, one file per
query, numbered in the order of the prompt:

| File | Contents | Status |
|---|---|---|
| `00-review-of-prompt1.md` | this file: review, decisions needed, plan | written |
| `01-branch-review.md` | per-branch disposition + merge log | after decisions |
| `02-runs-vs-scans.md` | Runs vs Scans (answered in §5 below; split out if it grows) | answered here |
| `03-dual-objective-ranking.md` | how multi-metric ranking works now, every place it is used | answered here in §6; split out once N confirms |
| `04-activity-candidates.md` | candidate activity metrics to optimize | answered here in §7 |
| `05-cone-pareto-package.md` | design of the pdt-server → generic Pareto package port | design sketch in §8 |

The prompt itself stays as written; this directory is the reply.

---

## 1. State of the repo (what I found on refreshing)

- **Wave 1 has been sitting unmerged since 2026-07-22 (12 weeks).** 11 branches
  were built from checkpoint §K″, each one commit ahead of `main@f6ccdab`, with
  green suites and nothing merged. Their worktrees under `.claude/worktrees/`
  still exist.
- **The board on `main` is stale. The current version is uncommitted.** HEAD
  `660b9ce` carries the pre-reconciliation board. The 2026-07-22 reconciliation
  (Decisions W1-a…d, Wave 1 tables, Active workstreams, Reconciliation) exists
  only in the working copy, along with uncommitted edits to `Research/README.md`,
  the checkpoint doc, `Docs/devlog/README.md` and `evoca_test.ipynb`. **This must
  be committed before any agent branches off `main`**, or every agent will work
  from the stale board.
- **Most Runs are untracked.** 14 of 87 `.evoca` files are in git. Several of the
  untracked ones are load-bearing: `Runs/2026-04-27_top_*` are the bases of
  campaign #4 (RD robustness). The disk cost is trivial (Runs 348 KB, Scans
  1.4 MB).
- **Branches outside Wave 1:**

| Branch | vs main | Assessment |
|---|---|---|
| `bigLUT` | 0 ahead, 157 behind | archival reference (5-ring LUT); CLAUDE.md cites it → **keep** |
| `softmouth` | 0 ahead | merged as PR #1/#2 → **delete** (local + origin) |
| `lut-hashes-accessor` | 0 ahead | merged M1 → **delete**; prune `/private/tmp/evoca_lut_hashes` |
| `pre6todo` | 1 ahead (2026-03-16, "test cell for 6 item To Do debugging") | almost certainly superseded → inspect, then delete or archive-tag |
| `worktree-agent-*` ×11 | 0 ahead, all at `f6ccdab` | leftover agent scaffolding → **delete** |

## 2. Overall take

The prompt is the right move. The program has more built than integrated: 11
finished branches, 6 drafted digests and about 85 runs, with none of it wired
into a structure that makes the next step obvious. In the prompt's order,
"clean up branches → per-workstream structure → per-workstream notebooks" is
also the correct dependency order. Three places need sharpening before agents
fire. Each is a real design conflict, not a wording issue:

1. **"Wave" is a delivery batch, not a scientific workstream** (§3a).
2. **A growing notebook and the authoritative campaign digests do different jobs** (§3b).
3. **Moving Runs/Scans per workstream conflicts with cross-workstream reuse and
   with path coupling in the harness** (§3c).

## 3. Design conflicts to resolve

### 3a. The workstream list

N's list is Wave, Criticality-bifurcation, Ring-complexity, Lineage-co-evolution,
Reaction-diffusion. "Wave 1" was the name of the 11-branch batch spawned from
§K″, and it mixes unrelated things: metrics (A_new, Φ, z-scores, eff_activity,
dyn-promote), infrastructure (env-agent-split, mu-genome) and new physics
(resource-driving, chem). Several existing results also have no home in the
5-name list: neutral-model methodology, pure-evo #1/#2, R1 viability brackets
and the resource-driving inverted-U.

**Proposed workstream set** (owner of each existing campaign shown):

| Workstream dir | Question it owns | Absorbs |
|---|---|---|
| `Activity-metrics` (= N's "Wave") | What is the right measure of adaptive evolutionary activity? | neutral-model-methodology, pure-evo #1/#2, causal-control v1/v2, metric-anew/zscore/eff-activity/dyn-promote, the 2026-04-27 evo_focus scan |
| `Resource-driving` | Does evolvability peak at intermediate driving, and is the axis `food_inc/tax`? | R1 viability brackets, `resource-driving` branch (C1/C3), metric-flux Φ (§L), §B figure |
| `Lineage-co-evolution` | Do LUT and egene co-evolve adaptively, or drift together? | #3/#3c coevolution-substrate, #8 (on hold), `lineage-return` (D2), H1, D3 |
| `Reaction-diffusion` | What RD regimes exist, and are they robust? | 2026-04-27 initial/refined/push_edges, #4 RD-robustness, the gallery |
| `Criticality-bifurcation` | Where are the transitions; what are their exponents? | `bifurcation.py` harness, future sweeps |
| `Ring-complexity` | How does the effective rule's ring reach evolve, and what does it cost? | ring-tax, ring-ladder A/B (#5), ring-stability two-track |
| *(later)* `Chem` | Vector resources | `chem` branch; no results yet, so no dir until there are some |

`env-agent-split` and `mu-genome` are infrastructure and belong to no workstream.
They are tracked on the board only.

### 3b. Growing lab notebook vs. digest

N wants a growing notebook per workstream, with the newest work at the bottom
and example runs in every section. The current README sets a different rule: one
digest notebook per campaign, committed with outputs, as the authoritative
finding. These do different jobs and fail differently if merged into one file. A
single growing notebook with committed figures from every section will reach
tens of MB, be slow to render on GitHub, and need a long re-execution if anyone
re-runs it top to bottom.

**Proposal: keep both, in the same directory.**

```
Research/<Workstream>/
  <Workstream>.ipynb        lab notebook: chronological, newest section at bottom,
                            each section = question → example runs (tailored probes)
                            → observations → directions/decisions
  ISSUES.md                 the single source of open questions + next directions
  digests/<campaign>.ipynb  finished campaign results (what Research/_drafts/ graduates into)
  Runs/   Scans/            this workstream's raw artifacts
```

A lab-notebook section that concludes gets promoted to a digest. The index in
`Research/README.md` links digests, not lab notebooks.

### 3c. Per-workstream Runs/Scans

I agree with the goal: the top-level `Runs/` is already a morass (62 files dated
04-27), and per-workstream dirs let file names shorten. Three issues need a
policy first:

- **Cross-workstream reuse.** `Runs/2026-04-27_top_*` were selected by the
  balanced RD+evo composite. They are the RD gallery's material, campaign #4's
  input, and comparators for pure-evo. Copying them into several dirs is a
  "fact stated twice" hazard. **Policy:** each artifact has exactly one owning
  workstream, and other workstreams refer to it by path. A small
  `Research/artifacts.csv` (old path → new path → owner) is the single map.
- **Path coupling in code.**
  - `evoca_explore._resolve_scan_dir()` defaults to the newest subdir of
    top-level `Scans/`.
  - `.evoca` recipes record `source_scan` by bare name.
  - Scan drivers and two root notebooks (`evoca_scans.ipynb`,
    `evoca_explore.ipynb`) hardcode `Scans/...`/`Runs/...` paths.
  - Torque launch/pull scripts assume the top-level layout.

  The harness needs a `workstream=` (or `root=`) argument plus a resolver that
  maps legacy `source_scan` names through `artifacts.csv`. Moves must be
  `git mv` so history follows. Commit the untracked Runs first, or they are
  lost from the record.
- **Timing.** Restructure **after** the branch merges. Several Wave-1 branches
  edit `python/evoca_explore.py`, and a path refactor landing in the middle of
  that integration would multiply the conflicts.

### 3d. Coordinating ISSUES with notebook and board

Markdown has no import, but both consumers can pull from one file:

- **Notebook:** the first cell of the newest section runs
  `display(Markdown(Path('ISSUES.md').read_text()))`. The committed output shows
  on GitHub but is a snapshot as of the last run, so the cell carries the
  file's git date.
- **Board:** a ~40-line `tools/build_board.py` that regenerates a marker-
  delimited section of `research_board.md`
  (`<!-- BEGIN workstreams -->…<!-- END -->`) by concatenating each
  `Research/*/ISSUES.md` under its heading. It runs by hand or as a pre-commit
  hook.

`ISSUES.md` is then the single source and the other two are derived, which is
the rule from the global lessons. The alternative of having the board only link
to each `ISSUES.md` is simpler but loses the "read one page" property that makes
the board useful.

### 3e. Merging authority

"Merge everything that seems mergeable" conflicts with three standing
conventions: the human merge gate, the pre-merge test-notebook protocol, and
review in a `/tmp` worktree. I read the prompt as N **granting** merge authority
for the additive/no-op Tier-1 branches, provided the combined suite is green and
each merge is verified. If so, I propose:

- **Tier 1** (7 additive branches): rebase in `/tmp/evoca_integ` in order
  `dyn-promote` → `metric-anew`, `metric-flux` (rebased on `build_sim`) →
  `metric-zscore`, `eff-activity`, `lineage-return`, `env-agent-split`. Run the
  combined suite and a **bit-for-bit parity check** (default-param run, same
  seed, identical trajectories before/after; this is the external check that
  "additive" is true). Then merge. One combined review notebook replaces seven
  per-branch ones.
- **Tier 2** (`resource-driving`, `mu-genome`, `chem`): stay at the human gate.
  Each claims a no-op default (food_per off, mutation off, N_chem=1 parity), so
  each gets the same parity check plus a test notebook. They merge only on N's
  explicit go. They could be dynamics-altering, which is the case the gate
  exists for.
- **`research-drafts`:** merge, then the restructure moves the drafts into
  `Research/<Workstream>/digests/` (still unapproved drafts, flagged as such).
- **W1-a/b/c** must be answered first or defaulted. W1-d is answered by this
  prompt. My defaults:
  - W1-a: add A_new/z-scores to ranking options **without changing**
    `EVO_METRICS` defaults; the new-activity choice is §7's job.
  - W1-b: keep p90 and record p99 alongside.
  - W1-c: keep the per-cell used-mask as **opt-in**.

## 4. Board restructure (`research_board.md`)

All agreed. Concretely:

- Move everything under "Archive: pre-checkpoint board" (lines ~138–357) to
  `research_board_archival.md`, leaving a one-line pointer.
- Add a **"How to read the codes"** box. The labels are `<letter><n>`, where the
  letter is a section of `Docs/checkpoint 3 July 2026.md` and the number is an
  action within it. The full table is **§I′ (consolidated action table, ~line
  1470)** and the roadmap is **§K′/§K″ (lines 1499/1594)**. For example, G1 = §G
  per-gene neutral model, action 1 = A_new. Each code is also annotated in place
  with its status, which answers N's question directly:

| Code | What | Kind | Status |
|---|---|---|---|
| D2 | per-lineage net-return reducer | **built** | on branch `lineage-return` → lands with merge |
| G1 | Bedau–Packard A_new | **built** | branch `metric-anew` |
| L1 | dissipation Φ accountant | **built** | branch `metric-flux` |
| H1 | adaptive-vs-drift program | **needs dev** (analysis) | blocked on D2+G1 merge |
| D3 | punctuation detector | **needs dev** | blocked on D2 |
| C3 | resource-driving campaign | **campaign** | blocked on C1 (`resource-driving`) |
| §B fig | `food_inc/tax` inverted-U figure | **campaign** | blocked on L1+G1 |

- The new section N asked for, **"Dev branches by workstream"**, lists each
  branch under its workstream with a 2–3 line conclude-and-merge sketch. Once
  the merges are done it will be short: Tier 2 plus any future branch.

## 5. Runs vs Scans (reminder)

- **Run** = one `.evoca` **recipe**: a JSON file with the full metaparams
  (init/final), the initialization (LUT/egenome/layout), N, seed, and provenance
  (`source_scan`, `source_config_idx`, `descriptor`). It is replayable
  interactively (SDL window, probes) or headless. Most were *emitted from* scans
  by `evoca_from_scan(_top)` / `nearest_*`, which is why names look like
  `top_7_cfg107_score0.80`.
- **Scan** = one **campaign directory**: `scan.py` (grid + `multiprocessing`
  driver, torque-capable), `scan_config.json` (fixed init), `results.csv` (one
  row of summary metrics per config, from `run_sim`), and `notes.md`. It is
  headless only. Docs: `Scans/README.md`.
- The relationship is many-to-one: a scan produces metrics, and selected rows
  become runs. That provenance is why per-workstream co-location makes sense.

## 6. How "dual-param" optimization is handled now

**It is not Pareto, and it is not really dual.** It is **scalarization**: an
unweighted mean of min-max-normalized metrics. There is no optimizer anywhere.
Scans are grids or LHS-free factorial grids, and "optimization" means
**post-hoc ranking** of grid rows.

```
candidates = rows with extinct==False and alive_density_mean <= 0.95
m̂_i = (m_i − min m) / (max m − min m)        per metric, within the scan's pool
score_i = mean over metrics in score_keys of m̂_i
```

- Implemented in `python/evoca_explore.py:722` (`evoca_from_scan_top`).
- Explained in `Scans/2026-04-27_scan_analysis.md` §"Composite score"
  (~line 130), including a **global** normalization variant across the 4 scans.
- The pure-evo variant (positive-clipped, 5 keys) is defined in
  `Scans/2026-05-17_pure_evo_analysis.md` ~line 61.
- A second, related tool is `nearest_evo` / `nearest_spatial`, a
  mean-|rank-difference| distance in a metric set, used to find neighbors of a
  chosen config rather than to optimize.

The "structure AND activity" case for RD was handled by **pooling** spatial
keys (`correlation_length_mean`, `largest_patch_temporal_std`) with evo keys
(`n_distinct_genomes_temporal_std`, `unique_top_genomes`) into one 4-key mean,
called the "balanced composite". Known weaknesses, which the cone/box approach
addresses:

- **(a) Compensation.** A config with strong structure and zero activity can
  outrank a balanced one. Scalarization rewards lopsidedness.
- **(b) Concave regions are unreachable.** A weighted mean only finds points on
  the convex hull of the front.
- **(c) Pool-relative normalization.** The same config scores differently in
  different scans. The 04-27 analysis had to bolt on a global normalization to
  compare scans.
- **(d) Silent weighting.** The weights are implicit (2 spatial + 2 evo keys
  happens to give 50/50 between axes; the 5-key pure-evo composite does not).

**Complete list of uses:**

| Where | Keys | Purpose | Produced |
|---|---|---|---|
| `Scans/2026-04-27_initial`, `_refined`, `_push_edges` | balanced 4-key (spatial+evo) | rank RD/evo regimes | `Runs/2026-04-27_top_1…20` (+auto_*) |
| `Scans/2026-04-27_evo_focus` (`scan.py:63`) | 4-key composite to pick base configs from source scan; then `EVO_METRICS` | evo-only ranking | `Runs/2026-04-27_top_evo_*`, `Runs/2026-04-29_top_evo_*` |
| `Scans/2026-04-27_scan_analysis.md` | balanced + EVO-only, global-normalized | cross-scan comparison | tables only |
| `Scans/2026-05-17_pure_evo_joint` / `_lut_only` | pure-evo 5-key / 4-key | rank pure-evo regimes | top tables in analysis md |
| `Scans/2026-05-17_pure_evo_egene_on_winners` (#3) | LUT-only pure-evo composite | choose top-12 base regimes | #3 inputs |
| `Scans/2026-05-18_rd_under_pure_evo` (#4) | inherits balanced composite (bases = 04-27 top_*) | input selection | #4 inputs |
| `evoca_scans.ipynb`, `evoca_explore.ipynb` (root) | ad hoc | interactive ranking | — |

The only **genuinely two-objective** question in that list is RD structure vs.
evolutionary activity (rows 1 and 6). The others are single objectives assembled
from several correlated activity proxies. That distinction matters for §8.

## 7. Candidate "activity" metrics to optimize after the metrics work concludes

| Candidate | Source | Status | Caveat |
|---|---|---|---|
| `excess_activity_slope` (legacy, whole-genome hash) | main | used in `EVO_METRICS` | **invalid for egene-driven runs**: shadow is LUT-only, activity hashes full genome (memory: shadow-models-lut-only-mismatch) |
| `excess_pc_slope` (component-normalized, Channon-2006 correction) | main | used in pure-evo composite | ≈0 everywhere in #1; weak discriminator |
| `dyn_excess_pc_slope` (500 fixed LUT-transition buckets) | main | pure-evo composite | **artifact on frozen-rich substrates** when `mu_lut=0` (memory: dyn-shadow-frozen-rich-substrate) |
| `eg_excess_pc_slope` (729 egene keys) | main | pure-evo composite | egene-only; 0 by construction in LUT-only runs |
| `dyn_activity` / `dyn_excess` as primary readout | `dyn-promote` (D1) | branch | same frozen-rich caveat |
| **A_new** (Bedau–Packard new activity, shared neutral threshold) | `metric-anew` (G1) | branch | one free parameter: threshold p90 vs p99 (W1-b) |
| per-bucket **z-scores** of eg/dyn vs shadow | `metric-zscore` (G2) | branch | scale-free, so comparable across scans (fixes §6c) |
| **eff_activity** (effective-rule species, union-mask key) | `eff-activity` (F1) | branch | counts only phenotypically realized rule content |
| lineage **net-return** | `lineage-return` (D2) | branch | adaptive signal rather than activity; candidate *second axis* |
| **Φ** per-capita dissipation | `metric-flux` (L1) | branch | not activity; the natural **second objective** / x-axis for §B |

My expectation, to be tested rather than assumed: **A_new on eff_activity
species, z-scored against the fixed-space shadow**, becomes the primary activity
objective. Net-return or Φ is the principled second axis, and then the cone/box
machinery of §8 has a real two-objective problem to work on.

## 8. Abstracting the restricted Pareto approach from pdt-server

What exists (`~/Projects/pdt-server/backend`, **R**):

- `ParetoFunctions.R`: `ParetoLayer2D` (non-dominated layer peeling),
  `ResolveParetoObjectiveRanges`, `ParetoHypervolume2D`.
- `ParetoAcceptanceFunctions.R`: regions over normalized `[0,1]²`. These are
  `SymmetricWedgeAcceptance` (the "cone", one parameter: half-angle about the
  diagonal), `BoxAcceptance`/`BoxQuantileAcceptance` (floors),
  `InverseBox*`, `Band*` and `ValidityMask`, plus `RelaxAcceptance`, which is
  the only thing allowed to widen a region.
- `SecondObjectiveFunctions.R`: pdt-specific (dose norms). Not ported.
- The header documents hard-won conclusions:
  - a wedge bounds **ratio**, not magnitude, so it gives no floors;
  - a box gives floors and admits lopsided candidates;
  - the asymmetric wedge was removed (it aims at the origin rather than
    pushing away from it, and it adds a second hand-set parameter);
  - normalization is pool-relative;
  - the number of Pareto-optimal points grows like log n.

  These must travel with the code, or someone will reimplement the asymmetric
  wedge.

**Proposal:** a standalone Python package (working name `conepareto`, own repo
under `~/Projects/`, `pip install -e` into EvoCA's `.venv`) containing:

- `pareto_layers(s1, s2)`
- `normalize(…, ranges=None | 'pool' | fixed)`
- the acceptance regions (wedge, box, box-quantile, band, inverse box)
- `relax(region, …, k)`
- `select(scores, region, k)`, which peels layers inside the region until k
  are taken, mirroring the pdt sampler
- `hypervolume_2d`
- the design header, ported as module docs

Its correctness check is **external**: a parity test that loads a saved pdt pool
(e.g. `DualObjectiveParetoSweepDiagnosticsGen1.RData`, exported to CSV once via
Rscript) and asserts that the Python layers and selections reproduce the R
output exactly.

Decisions this raises:

- **Language and home.** Python is needed for EvoCA, but pdt-server stays R.
  The two implementations of one algorithm are then held together only by the
  parity test, so the test is mandatory, not optional. The alternative (an R
  package called via rpy2) is worse for EvoCA's torque workflow.
- **Axis composition in EvoCA.** Pareto is 2-D. EvoCA metric sets have 4–5
  keys, so each axis must be defined first. For RD that is axis 1 = structure
  (a z-scored spatial composite, or `correlation_length` alone) and axis 2 =
  activity (§7's winner). Without that decision the package has nothing
  principled to rank.
- **Which region by default.** Per the pdt header: use a **box** when the
  question is "good at both" (RD gallery selection) and a **wedge** to trace
  the trade-off curve (characterizing the structure/activity frontier as a
  result in itself).
- **First validation in EvoCA:** re-select from the 04-27 balanced scans with
  box and wedge, and report the overlap with the composite's `top_1…20`. If the
  overlap is high, the composite was fine for that pool; if not, the difference
  set is the interesting material.

## 9. Proposed execution plan (parallel, with model per task)

**Phase 0: me, now, sequential.**
- Get N's answers to §10.
- Commit the uncommitted board/README/docs on `main`.
- `git add` the untracked Runs.

**Phase 1: parallel.**

| Agent | Model | Task | Isolation |
|---|---|---|---|
| Integrator | Opus (me) | Tier-1 rebase, parity check, combined suite, combined review notebook, merges; Tier-2 parity + notebooks; branch cleanup; `01-branch-review.md` | `/tmp/evoca_integ` |
| Pareto port | **Sonnet** | `conepareto` package + R parity test + docs; `05-cone-pareto-package.md` | new repo |
| Artifact census | **Sonnet** | read-only: assign every Run/Scan/digest/notebook to a workstream, draft `artifacts.csv`, flag ambiguous ones | read-only on main |
| Board split | **Sonnet** | archival split, codes glossary, status annotations | main (docs only, after Phase-0 commit) |

**Phase 2: after the merges.**
- Harness path refactor + resolver (Opus, small but load-bearing).
- `git mv` per `artifacts.csv`.
- `tools/build_board.py`.

**Phase 3: parallel, one Sonnet agent per workstream (6).**
- Seed `Research/<WS>/` with `ISSUES.md` from board + checkpoint.
- A lab notebook whose sections replay curated example runs headless, with the
  probes relevant to that workstream (activity strip, eg_pop, spatial snapshots,
  repro-age histogram, etc.).
- Move the `_drafts` digest into `digests/`.

Agents write only inside their own workstream dir. I review each for figure
sanity against the original scan results before committing.

**Phase 4.** First EvoCA use of the Pareto package (§8 validation). Activity
re-selection waits until the metrics work settles W1-a/b.

Total agent count is about 10 across phases, so separate `Agent` spawns suffice;
no workflow orchestration is needed.

## 10. Decisions needed from N

1. **Workstream set:** accept §3a's six (renaming "Wave" → `Activity-metrics`,
   adding `Resource-driving`)?
2. **Notebook model:** lab notebook + `digests/` per workstream (§3b), or truly
   one notebook?
3. **Merge authority:** am I authorized to merge Tier-1 after a green suite and
   bit-for-bit parity, with a single combined review notebook in place of the
   per-branch protocol? Tier-2 stays gated?
4. **W1-a/b/c:** accept my defaults (§3e)?
5. **Pareto package:** standalone Python repo + R parity test (§8)? Name/location?
6. **Branch deletions:** `softmouth`, `lut-hashes-accessor`, `worktree-agent-*`
   (and `pre6todo` after inspection), including the origin copies?
7. **Response location:** this `Docs/Oct2026/` layout OK?

## 11. N's decisions (2026-10-15)

- Workstreams: **accept the six** of §3a (Wave → `Activity-metrics`; add `Resource-driving`; `Chem` later).
- Notebook model: **lab notebook + `digests/`** per workstream (§3b).
- Merge authority: **Tier-1 merges after a green combined suite + bit-for-bit default-param parity**, with one combined review notebook. Tier-2 stays at N's gate.
- Pareto: **standalone Python repo** + R parity test (§8).
- Not asked, proceeding on defaults (override any time):
  - W1-a/b/c per §3e.
  - Responses in `Docs/Oct2026/`.
  - Delete fully-merged *local* branches and the stale `worktree-agent-*` branches.
  - **Origin** branch deletion (`softmouth`) and `pre6todo` wait for an explicit OK.
