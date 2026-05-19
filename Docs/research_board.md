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
resolved/done; **FYI items (no decision) last**. None open right now
(M1 merged 2026-05-19).

| # | Workstream | Decision | Raised |
|---|------------|----------|--------|
| M1 ✓ MERGED | #8 enabler | **User-approved & merged 2026-05-19** (`f8d6d66`, `--no-ff`; reviewed in worktree via `lut-hashes-accessor-test.ipynb`). dylib rebuilt on `main`, 42/42, clean `-Wall`. #8 campaign now buildable. | 2026-05-18 |
| R1 ✓ DONE | research / next pure-evo scan | **Resolved 2026-05-18.** Both levers bracketed: m_scale productive optimum is **interior ≈2.5–3.5** (NOT runaway; 5.0 maximally viable but unproductive); food_inc optimum **HIGH** 0.013–0.018 (sub-0.010 closed). New finding: #1's "low mu_lut wins" is **conditional on m_scale ≤ 2.5**; U-shaped viability (36 % extinct, worst at intermediate m_scale under scarcity). Detached run survived a connectivity blackout. `Scans/2026-05-18_R1_mscale_bracket`. | 2026-05-17 |
| D4 ✓ DONE | recipe export | Extend `metaparams_final`/`params()`/`_DEFAULTS` to the full param set (`mu_egenome`, `p_dup_egene`, `tax_per_egene` were also missing). **User: definitely extend.** Done `a1f70e2` (on `ring-tax`, merged) + round-trip test extended. | 2026-05-17 |
| D2 ✓ DONE | S2b / `tax_lut` | Add `tax_lut`+`tax_ring` to recipe export so genelife ring-ladder configs reproduce. **User: add both.** Done `dca4afa` + `test_recipe_roundtrip.py`. | 2026-05-16 |
| D3 ✓ RESOLVED | architecture / S2d | Worktree isolation breach (agent used absolute main paths + `cd` to main). Two-layer fix: (1) orchestrator asserts `main` clean before/after every agent — mechanism-agnostic, the real guarantee; (2) agent path-jail prompt. Relaunch verified-isolated. Standing policy for all future spawns. | 2026-05-16 |
| F1 (FYI / follow-up) | metrics | #3c surfaced: S2d eg/dyn fixed-space shadows accumulate across a mid-run phase switch; slope-windowing mitigates but absolute excess is biased and `dyn_excess_pc` is an artifact on frozen-rich substrates. Small additive C change (re-baseline shadows at a switch point) recommended before any future phase-switch campaign leans on absolute excess. Not blocking #8 (lineage metric needs no shadow). | 2026-05-19 |
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
| 4 | past evo+spatial winners → pure-evo: does RD vanish? | ✅ done | RD spatial structure NOT robust under pure-evo (corrL final/mean med 0.23, 0/30 extinct); 8/10 keep turnover while spatial washes out, 2/10 collapse to static colony-blob (lowest dyn_excess) — NP's "suspect-the-metric" cases. Recovered after DNS-drop (scan finished autonomously on torque) |
| 5 | genelife ring-ladder A/B | queued | — |
| 6 | bifurcation sweep | queued | — |
| 7 | patch-transfer 2×2 | queued | — |
| 8 | direct lineage co-evolution metric | **enabler MERGED (`f8d6d66`, M1).** `evoca_get_lut_hashes()` on `main`. Campaign harness now buildable: periodic `get_lut_hashes()` + active-egene-key + lineage (birth_id/parent_id) → joint-retention vs neutral marginal product. | 2026-05-17 |
| 3c | evolve LUT-only → in-situ freeze → evolve egene | ✅ **done** (mixed/methodological) | 30/30 A/C determinism held. **3b's gap was mostly GoL-substrate death** (B 10/10 extinct; B′ viable but eg_excess≈0). Residual genuine freezing penalty C>A modest+consistent (22/30, mean +0.6, scales with substrate quality) but signal weak (≲1 vs #1's 20–135) & plateaued by 8k. **dyn_excess_pc artifact** on frozen-rich substrate (A≈+335, not biology). Decisive clean test now → #8. Spec `Docs/plan-3c-coevolution-substrate.md` |
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

- **2026-05-19** — **#3c done — mixed/methodological result.**
  30/30 A/C determinism held (shared-fork valid). #3 3a≫3b gap is
  **mostly GoL-substrate death** (B 10/10 extinct; B′ viable but
  egene-excess≈0), with a **small consistent residual** freezing
  penalty C>A (22/30, mean +0.6, scales w/ substrate quality) on a
  **weak, plateaued** egene signal. Surfaced a **dyn_excess_pc
  artifact** for frozen-rich substrates (memory:
  dyn-shadow-frozen-rich-substrate). Clean decisive co-evolution
  test now firmly = **#8**. Notes/devlog written; new follow-up F1.
- **2026-05-19** — **#3c spec finalized & archived**
  (`Docs/plan-3c-coevolution-substrate.md`). User approved all 3
  decisions: R1-informed single regime (not #2-winner substrates —
  R1 confound), phase-1 3000 / phase-2 8000, both B & B′ GoL
  anchors. Added a first-class long-timescale ongoing-evolution
  diagnostic (user: wants evidence evolution still proceeds robustly
  at 8000 ticks). Crux = A-vs-C shared-fork via S1 determinism.
  Next: build harness on `main`.
- **2026-05-19** — **M1 merged.** `lut-hashes-accessor` →
  `main` (`f8d6d66`, `--no-ff`), user-approved after worktree review.
  dylib rebuilt on `main`, 42/42, clean `-Wall`. #8 campaign
  unblocked. #3c parameters under discussion (user leanings +
  agent recommendation integrating the R1 finding).
- **2026-05-19** — **Pre-merge test-notebook protocol established.**
  New standing convention (user): a brief runnable test notebook
  accompanies a branch at its merge gate for hands-on human review
  (exceptions allowed). First instance:
  `lut-hashes-accessor-test.ipynb` committed on the branch
  (`41550c3`), 5 cells mirroring `test_lut_hashes.py`, validated
  end-to-end (all PASS). Branch `lut-hashes-accessor` is checked out
  in worktree `/private/tmp/evoca_lut_hashes` — cannot also be
  `git checkout`'d in the main folder (git one-branch-per-worktree
  rule); inspect there or request the worktree be released.
- **2026-05-18** — **#8 enabler built; #3c rescoped (parallel).**
  Per user, #8 accessor and #3c proceed in parallel (independent:
  C-change-in-worktree vs Python-only-on-main). #8: `lut-hashes-
  accessor` branch — `evoca_get_lut_hashes()` (per-cell LUT-only
  FNV-1a, reuses existing `lut_hash_fn`, additive, no approx;
  C+ctypes+wrapper+2 tests; 42/42, clean `-Wall`, dylib byte-for-byte
  unchanged). Local, **in-review → M1 merge gate**. #3c: user dropped
  patch/S2c — now a single continuous sim, evolve LUT-only then flip
  `mu_lut=0`+enable egene **in place** (runtime setters already
  exist → no C change, runs on `main`). Design settled; awaiting
  user's parameter writeup (phase lengths, phase-1 params, configs,
  control arms).
- **2026-05-18** — **R1 done, recovered & analysed.** Launched
  detached on torque (`nohup`+`disown`, `.venv` python); a
  multi-minute VPN/torque connectivity blackout (torque drops ICMP
  under load — SSH is the true health check; `charge` same-subnet
  ping proved the tunnel healthy) left the detached run untouched —
  finished in 743 s, recovered exactly like #4. **Findings:**
  m_scale productive optimum interior ≈2.5–3.5 (runaway disproved);
  food_inc optimum high 0.013–0.018 (low extension closed); #1's
  low-mu_lut optimum is conditional on m_scale ≤ 2.5 (new
  mu_lut×m_scale interaction); U-shaped viability (36 % extinct).
  notes/devlog/board written; board R1 resolved.
- **2026-05-18** — **Campaign #4 recovered & analysed.** Local DNS
  failure → reboot dropped the launching SSH session; the #4 scan had
  already **finished autonomously on torque** (`run.log`: "Done in
  89s", 30 rows, 0 errors). Results pulled to local; analysis,
  `notes.md`, devlog section written. **Finding:** RD spatial
  structure not robust under pure-evo (corrL final/mean median 0.23,
  0/30 extinct); the 2 static colony-blob configs are the
  lowest-`dyn_excess` ones — corroborates NP's intuition. R1
  unblocked (run detached this time to survive disconnects).
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
