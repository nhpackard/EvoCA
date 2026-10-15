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
| W1-d | Integration | **RESOLVED by N 2026-10-15.** Tier-1 merges after a green combined suite + bit-for-bit default-param parity, one combined review notebook. Tier-2 (`resource-driving`, `mu-genome`, `chem`) stays at N's gate. |
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
| C1 | `food_per`: temporal regeneration period for food | built | branch `resource-driving`; Tier-2, N's gate |
| C2/L2 | float `env_mask` + relaxation toward mask; spatial gradient niches | needs dev | follows C1, same branch |
| L1 | dissipation accountant: sum of tax + death + overflow; per-capita throughput Φ (the §B x-axis) | built | **merged** 2026-10-15 (`5de867d`) |
| §B figure | scan summary of outcome vs `food_inc` and `food_inc/tax` (the inverted-U test); needs B1/B2 notebook | campaign | blocked on L1 + G1 merging |
| C3 | "temporal lumpiness gives evolutionary opportunity?" as a falsifiable campaign | campaign | blocked on C1 merging |

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
| A2 | `mu-genome`: per-cell rate-gene storage (mutation off), later evolvable rates | built (storage only); rest needs dev | branch `mu-genome`; Tier-2, N's gate |
| E1/E3 | chem branch: N_chem machinery (sparse store, N_chem=1 parity) | built (milestone 1) | branch `chem`; Tier-2 |
| E2 | chem milestone 2: food + waste (N_chem=2), needs a reaction-execution engine | needs dev | blocked on milestone-1 merge; high risk |

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

## Dev branches by workstream

**Tier-1 merged 2026-10-15** (N-approved: green combined suite + bit-for-bit
default-param parity). Six branches in `5de867d`; `lineage-return` in
`6548106` after fixing a pre-existing `enable_lineage()` call-order
segfault (`fadb242`). Branch refs deleted. Record:
[`../Docs/Oct2026/01-branch-review.md`](../Docs/Oct2026/01-branch-review.md).
Tier-2 (`resource-driving`, `mu-genome`, `chem`) stays at N's gate.
Sketches below are plans, not commitments.

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
- `resource-driving` (72adb59, C1): Tier-2. Claims `food_per` off = no-op:
  verify bit-for-bit parity, add a pre-merge test notebook, N reviews in a
  `/tmp` worktree and decides. C2/L2 continue on the branch after C1 merges.

### Lineage-co-evolution

- `lineage-return` (D2): **merged** `6548106` 2026-10-15; D3 and H1 can now start on `main`.

### Reaction-diffusion, Criticality-bifurcation, Ring-complexity

- No open dev branches. (Ring Track-2 bounded mutation will get its own
  worktree branch when started.)

### Infrastructure

- `env-agent-split` (A1/A3): **merged** `5de867d` 2026-10-15.
- `mu-genome` (6c0a360, A2 storage): Tier-2. Storage-first with mutation off;
  parity check + test notebook; N's gate. Distribution-seeding assay and
  evolvable rates are follow-on branches/commits.

### Chem (no workstream dir until results exist)

- `chem` (74cb1e5, E1/E3 milestone 1): Tier-2. N_chem=1 parity + sparse
  store; parity check + test notebook; N's gate. Milestone 2 (E2) continues
  on the branch afterward.

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
