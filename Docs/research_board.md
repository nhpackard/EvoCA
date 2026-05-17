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
| — | (none yet) | | |

---

## Workstreams

| ID | Branch | Isolation | Owner | Status | Last result | Updated |
|----|--------|-----------|-------|--------|-------------|---------|
| S1 determinism fix | `main` | none (main, gating) | agent | **in-review** | Added `evoca_set_seed()` (single xorshift32 `g_rng`; seed==0 remapped to 0x12345678); `run_sim` reseeds C RNG per call. tests/test_determinism.py: same-seed runs identical, different seeds differ. `pytest -q tests/` = 10 passed (8 existing + 2 new), clean `-Wall` build. | 2026-05-16 |
| S2a lineage field | `lineage-field` | worktree | (pending S1) | queued | — | 2026-05-16 |
| S2b ring-complexity tax | `ring-tax` | worktree | (pending S1) | queued | — | 2026-05-16 |
| S2c patch transfer | `patch-transfer` | worktree | (pending S1) | queued | — | 2026-05-16 |
| S2d eg/dyn fixed-space shadows | `main` | none (analysis) | (pending S1) | queued | — | 2026-05-16 |
| S2e bifurcation-diagram harness | `main` | none (analysis) | (pending S1) | queued | — | 2026-05-16 |
| S3 nightly digest | n/a | n/a | orchestrator | proposed | — | 2026-05-16 |

Status vocabulary: `queued` → `launched` → `in-review` (PR/branch
ready, awaiting human merge decision) → `merged` / `parked`.

---

## Log (newest first)

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
