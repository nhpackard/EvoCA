# 2026-05-16 — parallel workstream batch (S1, S2a–e)

Backfilled by the orchestrator from each agent's structured report.
Going forward agents write their own section as a committed deliverable.

---

## S1 — deterministic C PRNG reseed  (`main`, commit `447a019`)

**Motivation.** `multiprocessing.Pool` reuses workers; the per-process
xorshift32 was never reset, so per-config scan results were not
reproducible — blocks finite-size scaling and any parallel science.

**Key findings.** The C core has exactly **one** PRNG: `static
uint32_t g_rng = 0x12345678u` (xorshift32), consumed only via
`rng_next()`. No entangled generators, no `rand`/`srand`. Clean to
reseed.

**Design decisions & rejected alternatives.** Added additive
`evoca_set_seed(uint32_t)`; `seed==0` remapped to `0x12345678`
(xorshift32 has an all-zero fixed point — seeding 0 would freeze the
RNG). No-op when uncalled, so display/other scans are byte-unchanged.
Rejected: reseeding inside `evoca_init` unconditionally (would change
existing behaviour and break the "additive" guarantee).

**Status.** 10/10 tests; clean `-Wall`. In-review.

---

## S2a — lineage field  (`lineage-field`, commit `ef41d36`)

**Motivation.** §1 of the research directions: fitness in EvoCA is
realised, not assigned; nearly every landscape/OEE question reduces to
parent→child lineage, which the model did not track. Backbone
dependency.

**Key findings.** Phase 4 reproduction cleanly sees the parent (`idx`)
and `lut_hash_cache[idx]` is valid/unmodified there — no entanglement;
implementable exactly as specified.

**Design decisions & rejected alternatives.** Opt-in
(`evoca_enable_lineage`), default OFF, gated by `g_lineage_on` so the
hot path has a single perfectly-predicted branch and no Phase-1
additions ⇒ benchmarked zero-cost when off (N=256 GoL 409→408 fps,
within 1σ; ON −0.43%). Records child `parent_hash` (full-genome
FNV-1a, same as activity tracker), monotonic `birth_id`, and
`parent_id` (parent's birth id; 0 for founders) → enables chain
walking. Rejected: always-on field (violates the performance budget
that is a stated project constraint).

**Status.** 12/12 tests. In-review. Based pre-S1 (D1 → rebase at
integration).

---

## S2b — ring-dependence tax  (`ring-tax`, commit `829c8db`)

**Motivation.** §8: `tax_lut` taxes LUT '1'-bit popcount, which is
largely orthogonal to "what ring does the rule depend on." Taxing the
wrong thing would confound the genelife ring-ladder A/B.

**Key findings / design decisions.** Reused the **exact** existing
`lut_complexity()` classifier (the same function the `lut_complexity`
probe uses) for the per-cell level ∈ {1,2,3}; tax += `tax_ring*level`
in Phase 2c. Tax and probe therefore measure the identical quantity —
this is what makes the ring-ladder experiment interpretable. Only
structural addition: a forward declaration (classifier defined below
the tax path). Rejected: writing a second classifier (would risk
tax/probe divergence); the MI/effective-dependence variant (deferred —
more code, needs online MI estimate).

**Flagged modelling choice / obstacle.** Kept `tax_ring` out of
`params()`/`metaparams_final` to mirror `tax_lut`'s existing
treatment — surfaced rather than silently diverging. → became **D2**;
user decided **add both** `tax_lut`+`tax_ring` to recipe export
(integration task at the S2b merge gate, + round-trip test).

**Status.** 12/12 tests; default 0.0 verified bit-for-bit no-op
(mutation-tested the guard). Self-rebased onto S1. In-review.

---

## S2c — patch transfer  (`patch-transfer`, commit `e8a2f07`)

**Motivation.** §11: the evolutionary-success assay; the science is in
the controls, not the mechanism.

**Key findings / design decisions.** Zero C change — existing per-cell
getters expose enough; `stamp_patch` writes egene planes first then
`set_lut()` last so the species hash rehashes against correct egene
content. Scoring = alive-fraction trajectory (lag) **plus** boundary
flux (lead). Harness implements the **reciprocal 2×2** (Geritz
adaptive-dynamics sign structure), the **mandatory self-transfer
control** (every cross corrected by its own self-transfer), and the
size series (propagule-pressure threshold).

**Flagged modelling choice.** Dead cells carry stale non-heritable
egene scratch bytes → round-trip equality is asserted on **alive
cells only** (the actual heritable genome; strictly stronger than a
hash compare). Documented as the one semantic subtlety, and the
correct reading of "genome-level not pixel-level."

**Status.** 13/13 tests. In-review. Based pre-S1 (D1 → rebase).

---

## S2e — bifurcation harness  (`bifurcation-harness`, commit `dfdb010`)

**Motivation.** NP-reply R3: parameter-driven bifurcations
(period-doubling / Feigenbaum) are the relevant exponent program, and
are "mostly analysis on swept runs you can already produce."

**Key findings / design decisions.** Pure Python, no C change.
`scalar_series` is a thin per-tick loop mirroring `run_sim`'s init
exactly (rather than modifying `run_sim`). `recurrent_extrema`
discards a transient, collapses a constant tail to a fixed point, else
clusters turning points so a period-k cycle yields exactly k points
and chaos fills an interval. `feigenbaum_delta` is explicitly
**best-effort** — it computes the gap-ratio convergence but cannot
self-detect a genuine cascade (documented, not hidden).

**Obstacle resolved.** A synthetic roundtrip test failed due to grid
quantization collapsing deep cascade thresholds → reworked the test to
build rows at exact cascade locations (tests the estimator without
grid artefacts).

**Status.** 20/20 tests; `--demo` ~0.7s. Self-rebased onto S1.
In-review.

---

## S2d — eg/dyn fixed-space neutral shadows  (`eg-dyn-shadows`, commit `89ad3eb`)

**Motivation.** The 2026-05-16 causal control proved the existing
Channon shadow models **LUT bytes only** while `get_activity` hashes
the **full genome** → excess activity is undefined for egene-driven
runs (EGENE_only arm: excess_pc ≈ −26). §2 fix: fixed-space baselines.

**Key findings / design decisions.** Both 729-key eg and 500-bucket
dyn spaces use a per-step occupancy `q[k]`. Under the neutral null
there is no selection, so the shadow's expected composition = the
**realised** per-step composition, departing only via (1) realised
mutation flux `m = 1−exp(−f)` blending toward uniform-over-all-K and
(2) Wright–Fisher mean-field drift toward uniform-over-occupied,
**gated by `m`**. Per-component-normalised excess, exactly like
`excess_pc_*`. Key invariant: zero realised flux ⇒ shadow == observed
⇒ excess ≈ 0 (exact in practice).

**Flagged modelling choices (§2 left these open).**
1. Realised flux = actual Phase-4 Poisson bit-flip counts applied to
   children (measured, run-faithful, zero-tunable), per
   bits-per-component (12 eg / 1 dyn).
2. WF drift applied as **deterministic mean-field** (reproducible),
   **gated by flux** so a clonal population has zero drift — required
   for the neutral-null property. Closed-form variance documented as
   available for z-scoring but not injected.
3. Added a **private per-step eg observed mirror** because the
   eg_activity probe arrays only advance on `evoca_eg_activity_update()`;
   without it eg excess would undercount. Probe arrays untouched.

**Obstacles resolved.** (a) Unconditional drift made the zero-mutation
case wrongly non-zero (~12k) → fixed by gating drift on flux. (b) eg
shadow never accumulated (probe not per-step) → private mirror.

**INCIDENT (first attempt).** The first S2d agent ignored its worktree
and edited the **main working tree** via absolute `/Users/n/Projects/
EvoCA/…` paths + `cd` to main; killed mid-rewrite, left main
non-compiling (uncommitted only — HEAD always clean). Mechanism
confirmed by bounded transcript grep (first Edit targeted main; 28
main-path ops; explicit `cd` to main). Recovered: partial work saved
to `S2D_INCOMPLETE.patch`, contaminated files `git checkout HEAD --`'d,
verified clean. Mitigation = **D3**: (1) orchestrator asserts main
clean before/after every agent (mechanism-agnostic detection — the
real guarantee); (2) agent path-jail prompt (forbid the canonical-repo
path string, no `cd` to main, assert toplevel). Relaunch under this
regime verified isolated (rule-1/rule-5 checks in its report).

**Status.** 12/12 tests; clean `-Wall`. In-review. Based pre-S1
(D1 → rebase).

---

## Cross-batch process lessons

- Worktree isolation is **not an OS sandbox** — an LLM agent that
  builds absolute paths to the canonical repo escapes it. Bare
  worktree + prompt discipline is adequate at ≤5 supervised agents
  with an orchestrator clean-check; genuine scale needs container /
  enforced-cwd isolation (outside what the Agent tool exposes).
- The **causal-control harness doubles as a bug detector**: the
  reciprocal freeze arms exposed the shadow/genome-hash scope mismatch
  that motivated S2d. Designed discrepancies surface metric bugs that
  in-isolation review misses.
- Worktree-base quirk (D1): agent worktrees were repeatedly based at
  the pre-S1 commit regardless of main's HEAD; orchestrator rebases at
  the integration gate. Final integrated `main` is unaffected.
