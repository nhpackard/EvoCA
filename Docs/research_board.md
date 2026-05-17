# EvoCA research board

Single coordination surface for parallel workstreams. Each agent
**appends** its result and any decision it needs; it does **not** block
waiting for a human. A human reads the "Decisions needed" section
first; everything else is status.

Branches never auto-merge. Agents propose; merge is a human gate
(dynamics-altering branches `ring-tax`/`mu-genome` especially — a
half-tuned tax or rate-gene scheme contaminates every run's results).

Isolation: each code branch runs in its own `git worktree` (separate
on-disk dir + its own `C/libevoca.dylib`). Analysis-only work runs on
`main`.

---

## Decisions needed (human) — read this first

| # | Workstream | Decision | Raised |
|---|------------|----------|--------|
| **D3** ✓ RESOLVED | architecture / S2d | Relaunch under the path-jail prompt + orchestrator clean-check completed **verified-isolated** (rule-1/5 checks in its report; post-check shows main clean). Guardrail validated. Standing policy for all future spawns: path-jail prompt + orchestrator main-clean assert before/after each agent; devlog per workstream. Original note: **was MITIGATED →** Mechanism investigated (bounded transcript grep) & confirmed: agent used absolute `/Users/n/Projects/EvoCA/…` paths and `cd` to main, bypassing its worktree. Two-layer fix deployed: (1) orchestrator asserts `main` clean before/after every agent — mechanism-agnostic detection, the real guarantee; (2) agent prompt path-jail — forbids the literal main-path string, no `cd` to main, asserts `git rev-parse --show-toplevel` is the worktree. S2d relaunched alone under this. Residual: prompt-level prevention is soft (LLM can disobey) — layer (1) is what makes it safe. | 2026-05-16 |
| D2 ✓ RESOLVED | S2b (+ pre-existing `tax_lut`) | `tax_ring`/`tax_lut` don't round-trip into `.evoca` recipes → genelife ring-ladder configs won't reproduce. **DECISION (user, 2026-05-16): add both.** → Integration task at S2b merge gate: add `tax_lut` AND `tax_ring` to `params()` / `metaparams_final` / recipe export + add a recipe round-trip test covering both. | 2026-05-16 |
| D1 | S2a/b/c/d | Worktree-base quirk: `lineage-field`, `ring-tax`, `patch-transfer`, `eg-dyn-shadows` based at `b04a0d1` (pre-S1), so they lack `evoca_set_seed`. Mitigation (no decision needed, FYI): orchestrator rebases each onto `main`@447a019 and re-runs its suite *before* the merge gate; final integrated `main` is unaffected. S2e (bifurcation) correctly has S1. | 2026-05-16 |

---

## Workstreams

| ID | Branch | Isolation | Owner | Status | Last result | Updated |
|----|--------|-----------|-------|--------|-------------|---------|
| S1 determinism fix | `main` | none (main, gating) | agent | **in-review** | Added `evoca_set_seed()` (single xorshift32 `g_rng`; seed==0 remapped to 0x12345678); `run_sim` reseeds C RNG per call. tests/test_determinism.py: same-seed runs identical, different seeds differ. `pytest -q tests/` = 10 passed (8 existing + 2 new), clean `-Wall` build. | 2026-05-16 |
| S2a lineage field | `lineage-field` | worktree | agent | **in-review** | `ef41d36`: opt-in per-cell parent_hash/birth_id/parent_id written in Phase 4 (§1). OFF zero-cost (N=256 GoL 409→408 fps, within 1σ); ON −0.43%. 12/12 tests pass (8+4). Independently confirmed D1 (based pre-S1). Merge-ready pending rebase onto main@447a019. | 2026-05-16 |
| S2b ring-complexity tax | `ring-tax` | worktree | agent | **in-review** | `829c8db`: `tax_ring` (default 0.0, bit-for-bit no-op verified w/ mutation test) = `tax_ring*level` reusing the exact `lut_complexity` classifier (tax & probe measure identical quantity). GUI slider added. Self-rebased onto main@447a019 (has S1). 12/12 tests pass. | 2026-05-16 |
| S2c patch transfer | `patch-transfer` | worktree | agent | **in-review** | `e8a2f07`: zero-C-change `extract_patch`/`stamp_patch` (genome-level) + `python/patch_transfer.py` §11 harness (reciprocal 2×2, self-transfer control, size series, boundary-flux scoring). 13/13 tests pass. Agent claims base had S1 — verify at integration. | 2026-05-16 |
| S2d eg/dyn fixed-space shadows | `eg-dyn-shadows` | worktree | agent (relaunch) | **in-review** | `89ad3eb`: fixed-space neutral baselines for eg (729) + dyn (500), realised-flux-driven WF-drift model, per-component excess; zero-flux ⇒ excess≈0, strong selection ⇒ ≫0. New `eg_/dyn_excess_pc_*` in run_sim. 12/12 tests, clean `-Wall`. Relaunch verified isolated (rule-1/5). Based pre-S1 (D1 → rebase). | 2026-05-16 |
| S2e bifurcation-diagram harness | `bifurcation-harness` | worktree | agent | **in-review** | `dfdb010`: `python/bifurcation.py` (scalar_series/recurrent_extrema/sweep/Feigenbaum δ) + tests, pure analysis no C change. Agent self-rebased onto main@447a019 (has S1) — D1 N/A here. 20/20 tests pass. `--demo` runs ~0.7s. Note: its tree carries S1's older board copy → resolve board to main HEAD at merge. | 2026-05-16 |
| S3 nightly digest | n/a | n/a | orchestrator | deferred | collector script to be built now; `/schedule` wrapper deferred until first batch produces something to digest (user-approved) | 2026-05-16 |

Status vocabulary: `queued` → `launched` → `in-review` (PR/branch
ready, awaiting human merge decision) → `merged` / `parked`.

---

## Log (newest first)

- **2026-05-16** — S2d relaunch completed verified-isolated
  (`89ad3eb`, 12/12). **All 5 workstreams in-review.** D3 resolved.
  `Docs/devlog/` created + backfilled (the design-evolution narrative
  now lives in-repo, not just chat). Next: orchestrator integration
  prep in a **dedicated worktree** (main's working tree carries the
  user's uncommitted notebooks → blocks rebase there): rebase the 3
  pre-S1 branches onto main, apply D2 to ring-tax, re-run each suite,
  then one consolidated merge review.
- **2026-05-16** — D2 resolved (user: add both `tax_lut` + `tax_ring`
  to recipe export; folded into S2b integration task). D1 acknowledged
  (no decision; FYI only).
- **2026-05-16** — **INCIDENT**: S2d (eg-dyn-shadows) edited the
  `main` working tree instead of its worktree, killed mid-rewrite →
  `main` working tree left non-compiling (uncommitted only; HEAD
  always clean). Recovered: legit board edits committed, S2d partial
  diff saved to `S2D_INCOMPLETE.patch`, contaminated tracked files
  `git checkout HEAD --`'d, build clean `-Wall`, 10/10 tests pass.
  Other 4 branches verified intact. Raised D3 — parallel spawns
  paused pending human decision on isolation guardrail.
- **2026-05-16** — S2b completed → in-review (`829c8db`, 12/12).
  Self-rebased onto S1. Raised D2 (tax_lut/tax_ring not in recipe
  export → genelife A/B configs won't round-trip). 4/5 done; S2d
  (eg-dyn-shadows) still running.
- **2026-05-16** — S2a, S2c, S2e completed → in-review (commits
  `ef41d36`, `e8a2f07`, `dfdb010`; tests 12/13/20 pass). All local,
  unmerged. S2e self-rebased onto S1; S2a confirmed pre-S1 (D1);
  S2c base to verify at integration. S2b, S2d still running.
- **2026-05-16** — S2a–e launched in parallel, each in its own
  worktree+branch (`lineage-field`, `ring-tax`, `patch-transfer`,
  `eg-dyn-shadows`, `bifurcation-harness`). All branch from local
  `main` @ 447a019 so they include the S1 seed fix. Agents do NOT
  edit this board; orchestrator maintains it from their reports.
  Branches stay local; merge is a human gate. S1 commit 447a019
  remains local-only (not pushed) pending human review.
- **2026-05-16** — S1 implemented & in-review. Single per-process
  xorshift32 (`g_rng` in C/evoca.c) confirmed — no entangled
  generators. Added additive `void evoca_set_seed(uint32_t)`
  (seed==0 → 0x12345678), `sim.set_seed()` ctypes binding, and a
  reseed call in `run_sim` after `sim.init()`. No `set_seed` call
  ⇒ behaviour unchanged. New `tests/test_determinism.py`; full
  suite 10/10 pass; clean `-Wall` build. Committed to `main`
  (not pushed). Awaiting human merge/verification before S2* launch.
- **2026-05-16** — Board created. S1 (determinism fix) launched alone
  first; everything else gated on it per the approved plan. S2* launch
  automatically when S1 reports complete.
