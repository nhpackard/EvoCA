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

Convention: **unresolved decisions at the top** (act on these); then
resolved/done; **FYI items (no decision) last**. None open right now.

| # | Workstream | Decision | Raised |
|---|------------|----------|--------|
| R1 | research / next pure-evo scan | Campaign #1 found the corrected-metric optimum at **`m_scale=2.5` (grid ceiling) + low `mu_lut`**, still unbracketed. Decision: a follow-up scan with `m_scale ∈ {2.5, 3.5, 5.0}` to bracket it. Non-blocking (the #1–#4 program proceeds first); recommend running it after #4. | 2026-05-17 |
| D4 ✓ DONE | recipe export | Extend `metaparams_final`/`params()`/`_DEFAULTS` to the full param set (`mu_egenome`, `p_dup_egene`, `tax_per_egene` were also missing). **User: definitely extend.** Done `a1f70e2` (on `ring-tax`, merged) + round-trip test extended. | 2026-05-17 |
| D2 ✓ DONE | S2b / `tax_lut` | Add `tax_lut`+`tax_ring` to recipe export so genelife ring-ladder configs reproduce. **User: add both.** Done `dca4afa` + `test_recipe_roundtrip.py`. | 2026-05-16 |
| D3 ✓ RESOLVED | architecture / S2d | Worktree isolation breach (agent used absolute main paths + `cd` to main). Two-layer fix: (1) orchestrator asserts `main` clean before/after every agent — mechanism-agnostic, the real guarantee; (2) agent path-jail prompt. Relaunch verified-isolated. Standing policy for all future spawns. | 2026-05-16 |
| D1 (FYI) | S2a/b/c/d | Worktree-base quirk: agent worktrees based pre-S1; orchestrator rebased each onto `main` at the integration gate (done — all 5 S1✓, merged). No decision was needed. | 2026-05-16 |

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

## Research campaigns (2026-05-17, R-track) — IN PROGRESS

Compute: torque (32 cores) + local. Findings in
`Docs/devlog/2026-05-17_research-campaigns.md`.

| # | Campaign | Status | Headline |
|---|----------|--------|----------|
| 1 | pure-evo LUT+egene | ✅ done | whole-genome excess≈0 everywhere; dyn/eg excess strongly +ve; low mu_lut + high m_scale wins (contradicts old broken-metric optimum) |
| 2 | pure-evo LUT-only | ✅ done | 0/180 extinct; egene-freeze *reduces* dyn-excess (egene co-evo amplifies rule selection); LUT-only wants high mu_lut+m_scale, low food_inc — mu_lut optimum is conditional on egene |
| 3 | #2 winners → egene scan (3a/3b) | ✅ done | 3a≫3b (eg_excess 20 vs 6; 3b 46% extinct); strong co-dependence BUT 3b froze GoL not evolved-LUT → confound not cleanly resolved; #8 now decisive, #3c queued |
| 4 | past evo+spatial winners → pure-evo: does RD vanish? | queued | — |
| 5 | genelife ring-ladder A/B | queued | — |
| 6 | bifurcation sweep | queued | — |
| 7 | patch-transfer 2×2 | queued | — |
| 8 | direct lineage co-evolution metric (ΔLUT×Δegene joint-retention vs neutral) | ⏳ next (DECISIVE — promoted) | needs no freezing; resolves the confound directly |
| 3c | evolve→extract evolved LUT→freeze→evolve egene (clean substrate test, uses S2c) | queued | the test 3b was meant to be |
| — | causal-control v2 (local) | ✅ done | shadows pass null controls EXACTLY (eg/dyn excess=0 when frozen); bug closed (eg_excess_pc +31 vs v1's −26); LUT+egene co-evo mutually amplify selection — corroborates #1/#2 |

## Integration status (2026-05-17) — ✅ MERGED & PUSHED

Done in a dedicated `/tmp/evoca_integ` worktree (main's tree carries
the user's uncommitted notebooks). Hygiene fix first: untracked
`tests/__pycache__/*.pyc` (was blocking rebases) → `289191b` on main.

| Branch | Commits over main | Delivers | Tests (own+base) | S1 | Rebased on main |
|--------|------|----------|-------|----|----|
| `bifurcation-harness` | 1 (`dfdb010`) | R3 bifurcation/Feigenbaum analysis (pure Py) | 20 | ✓ | ✓ |
| `ring-tax` | 2 (`829c8db`+`dca4afa` D2) | §8 ring-dependence tax + D2 recipe round-trip | 15 | ✓ | ✓ |
| `patch-transfer` | 1 (`e8a2f07`) | §11 extract/stamp + reciprocal-2×2 assay | 15 | ✓ | ✓ |
| `lineage-field` | 1 (`ef41d36`) | §1 opt-in lineage (zero-cost off) | 14 | ✓ | ✓ |
| `eg-dyn-shadows` | 1 (`89ad3eb`) | §2 fixed-space eg/dyn neutral shadows | 14 | ✓ | ✓ |

- Individual rebases onto main: **clean, no conflicts.**
- Trial sequential merge of all 5: **zero cross-branch conflicts**,
  combined build clean, **full combined suite 38/38 passed**.
- Branch refs are the rebased versions; merge-ready.
- **MERGED** into `main` (order lineage→eg-dyn→ring-tax→patch→bifn):
  merge commits `4900145 1164ecb 8783059 2bb9cd5 9fdb3e6`, build
  clean, **40/40 combined suite**. D4 also executed (`a1f70e2`).
- **Pushed**: `b04a0d1..9fdb3e6 → origin/main`.
- `S2D_INCOMPLETE.patch` deleted (S2d redone cleanly). See
  `Docs/devlog/commit-log.md` for the running ledger.

## Log (newest first)

- **2026-05-17** — **All 5 merged to `main` + pushed**
  (`b04a0d1..9fdb3e6`). D4 executed (user: extend). 40/40 suite.
  `S2D_INCOMPLETE.patch` deleted. Devlog merge-review appended;
  `commit-log.md` ledger started (standing convention). Decisions
  table reordered (unresolved-first; none open). Batch complete.
- **2026-05-17** — Integration prep complete. Repo-hygiene fix
  (`289191b`: untrack `tests/__pycache__`). All 5 branches rebased
  onto main individually clean; D2 executed on `ring-tax`
  (`dca4afa`, +D4 flagged); trial all-5 merge = zero cross-branch
  conflicts, 38/38 combined. Verified in `/tmp/evoca_integ`; nothing
  merged/pushed. Awaiting merge gate.
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
