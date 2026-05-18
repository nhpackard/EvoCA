# Commit / merge log

Running ledger of every commit and merge that lands on `main`
(newest first). One line each: `<short-hash> <type> — <summary>`.
Convention: the orchestrator appends here whenever it commits to or
merges into `main` (and after every push, note the push range). This
is the at-a-glance "what has changed in the repo" record; the devlog
batch files carry the *why*.

**Hash rule (avoids the self-reference trap):** a commit cannot
contain its own final hash. So a substantive commit's hash is
**backfilled on the next ledger touch**, not written self-
referentially. Pure ledger-maintenance commits are not themselves
logged. If an entry ever shows a hash not in `git log main`, it is a
stale pre-amend hash — backfill it from the real history.

## 2026-05-17

- `pending` commit — #3 finding + honest design-limitation; #8 promoted, #3c queued
- `e5ed170` commit — consolidated scan analysis doc + per-dir notes.md
- `37187ee` commit — queue campaign #8 (lineage co-evolution metric)
- `b8c63a5` commit — causal-control v2 finding (metric validated, bug closed)
- `9bda882` commit — #3 scan (egene-on-winners, 3a+3b paired)
- `86f2b3f` commit — board: #3 running
- `a6dbcb3` commit — #2 finding (LUT-only; mu_lut optimum conditional on egene)
- `1eaafcc` commit — campaign devlog + board R-track + #1 finding
- `0abb72d` commit — pure-evo #2 (LUT-only) scan + #1 results.csv
- `bfee1b1` commit — hygiene: gitignore Scans run.log + .venv
- `e6525b5` commit — Scans/2026-05-17_pure_evo_joint/scan.py (priority #1 pure-evo LUT+egene scan)
- `5c74caf` commit — `tests/README.md`: test-suite index (40 tests / 10 files, one-line descriptions)
- `6ba31ad` commit — docs: merge review appended; commit-log ledger; board finalised
- `9fdb3e6` merge — `bifurcation-harness` (R3 bifurcation/Feigenbaum)
- `2bb9cd5` merge — `patch-transfer` (§11 extract/stamp + 2×2 assay)
- `8783059` merge — `ring-tax` (§8 ring tax + D2 + D4)
- `1164ecb` merge — `eg-dyn-shadows` (§2 fixed-space neutral shadows)
- `4900145` merge — `lineage-field` (§1 opt-in lineage)
- `a1f70e2` commit — D4: full metaparam set in recipe export *(on ring-tax, merged via 8783059)*
- `dca4afa` commit — D2: tax_lut/tax_ring recipe round-trip *(on ring-tax)*
- `97777e1` commit — board: integration verified
- `289191b` commit — hygiene: gitignore __pycache__; untrack tests pyc
- `0017dc7` commit — devlog created (design-evolution narrative)
- *(pushed: `b04a0d1..9fdb3e6 → origin/main`)*

## 2026-05-16 (pre-integration, on main)

- board / D-tracking commits (see research_board.md log for detail)
- `447a019` commit — S1: deterministic C PRNG reseed (`evoca_set_seed`)
- `b04a0d1` commit — exp: reciprocal causal-control campaign (§7.3)
- `febd57e` commit — test: pytest suite (GoL/tax/metrics)
- `35e249f` commit — evoca_explore: component-normalised excess (§2 fix)
- `6ade765` commit — docs: deep-analysis appendix + NP reply + refs
