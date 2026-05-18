# Commit / merge log

Running ledger of every commit and merge that lands on `main`
(newest first). One line each: `<short-hash> <type> — <summary>`.
Convention: the orchestrator appends here whenever it commits to or
merges into `main` (and after every push, note the push range). This
is the at-a-glance "what has changed in the repo" record; the devlog
batch files carry the *why*.

## 2026-05-17

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
