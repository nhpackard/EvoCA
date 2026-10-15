# Branch review and merge record (2026-10-15)

This answers the "cleanup of dev branches" part of `Docs/prompt1-Oct2026.md`.
N's decisions are recorded in `00-review-of-prompt1.md` §11 and in this
session.

## Outcome

| Branch | Disposition | Commit on `main` |
|---|---|---|
| `dyn-promote` (D1) | **merged** | `5de867d` |
| `metric-anew` (G1) | **merged** | `5de867d` |
| `metric-flux` (L1) | **merged** | `5de867d` |
| `metric-zscore` (G2) | **merged** | `5de867d` |
| `eff-activity` (F1) | **merged** | `5de867d` |
| `env-agent-split` (A1/A3) | **merged** | `5de867d` |
| `lineage-return` (D2) | **merged after fix** `fadb242` | `6548106` |
| `resource-driving` (C1) | Tier 2: waiting for N's gate | — |
| `mu-genome` (A2) | Tier 2: waiting for N's gate | — |
| `chem` (E1/E3) | Tier 2: waiting for N's gate | — |
| `research-drafts` | merges with the Research restructure (drafts move into `Research/<WS>/digests/`) | — |
| `bigLUT` | archival, kept | — |
| `softmouth`, `pre6todo` | deleted locally and on origin (N's OK). `pre6todo` was one test-cell commit from Mar 16, recoverable as `d413c18` | — |
| `lut-hashes-accessor`, 11 × `worktree-agent-*`, the 7 merged Tier-1 branches | deleted locally | — |

`main` is pushed to origin at `6548106`.

## Method

1. Created the `/tmp/evoca_integ` worktree on `main`. Each Tier-1 branch was
   one commit, so each was cherry-picked in order (`dyn-promote` first,
   because it extracts `build_sim`). After every step:
   - `make` with `-Wall`, which gave no warnings;
   - the full `pytest` suite;
   - a parity probe against pre-merge `main`.
2. **Parity probe**, now in the repo as `tools/parity_probe.py`. It uses a
   viable pure-evo config (N=128, mutation and egene on, final alive density
   0.34) and checks:
   - `run_sim` scalars for 2 seeds;
   - SHA1 of v/alive/f/F/egenome/LUT-hashes every 100 steps of a
     1200-step trajectory.

   An additive branch must give no differing scalar and an identical
   trajectory; new keys are allowed.

   **Caveat found along the way:** the first fixture went extinct and would
   have passed trivially. The tool now prints the extinction status so a pass
   can be trusted.
3. Results:
   - All six branches gave BAD 0 and IDENTICAL.
   - The test count rose from 42 to 50.
   - New `run_sim` keys: `a_new_*` and the `diss_*` family.
   - `metric-zscore` and `eff-activity` add accessors only.

## The conflict

`lineage-return` and `metric-flux` both edit the Phase-2c tax block and the
Phase-4 eviction block of `C/evoca.c`. Both edits are additive, so the
resolution keeps both:
- `float f_before` (flux) plus the lineage tax accumulation;
- the evicted-food sink (flux) plus `lret_flush_cell(child)` (lineage).

The duplicated `g_deaths_last++` was removed.

**Semantic note for N (no action needed):**
- The lineage reducer counts **nominal** tax `t`.
- The flux accountant counts tax **actually removed**, which is the pre-tax
  `f_priv` when the tax kills the cell.

These agree except on tax-death steps. Each is right for its own purpose:
lineage cost charged vs. energy that left the system. They should not be
cross-checked against each other without that adjustment.

## The `lineage-return` crash

- **Symptom:** segfault at `C/evoca.c` `lin_tax_cell[i] += t` (the address
  was NULL) whenever lineage was turned on with `sim.enable_lineage(True)`.
- **Root cause:** a pre-existing API trap from `lineage-field` (`c9c6a0a`,
  merged May 2026). `lineage-return` only added another array the trap
  reaches.
  - The C code allocates the lineage arrays only inside `evoca_init`, when the
    flag is already set.
  - Python's `init()` always sets the flag from its own `lineage=False`
    default.
  - So `enable_lineage()` **before** `init()` was silently undone: zero
    lineages recorded. **After** `init()`, it turned on hot-path writes into
    NULL arrays and crashed.
  - The bug was also present on pre-merge `main`.
- **Impact on past results: none.** Every caller (the #8 lineage scan and all
  tests) uses `init(..., lineage=True)`, which was always correct.
- **Fix (`fadb242`):**
  - `evoca_enable_lineage` now allocates the arrays and reducer table when
    lineage is enabled on an initialized lattice. Cells alive at that moment
    become founders.
  - Python's `enable_lineage()` before `init()` raises a `RuntimeError` that
    points to `init(lineage=True)`.
- **Regression tests:**
  - The crash case runs in a subprocess, so a regression reports as a test
    failure instead of crashing the whole run.
  - It also checks reducer consistency: `sum(pop) == alive`, and lineage
    intake/tax sums equal the run totals.
  - Both tests were **shown to fail with the fix reverted** (rc=-11, no
    raise).
- **Verification after the fix:**
  - 57/57 tests;
  - parity identical with `main`;
  - lineage ON vs OFF gives an identical state and identical per-step
    dissipation.

## Not done as planned

- **No combined review notebook was built.** The review proposed one in place
  of per-branch notebooks. The branches' own pre-merge test notebooks came in
  with the merges and now sit at the **repo root**:
  - `dyn-promote-demo.ipynb` (+ `_fig.png`);
  - `{eff-activity, env-agent-split, lineage-return, metric-anew, metric-flux, metric-zscore}-test.ipynb`;
  - the older `lut-hashes-accessor-test.ipynb`.

  Proposal: during the Research restructure, move them to
  `Research/Activity-metrics/digests/merge-tests/` (or the owning workstream)
  rather than build a new combined notebook. They already cover each branch.

## Tier 2 next steps

For each Tier-2 branch:
1. Rebase onto `6548106`.
2. Run `tools/parity_probe.py` with the new feature at its no-op default
   (`food_per` off, mutation off, N_chem=1).
3. Run the suite.
4. Add a pre-merge test notebook.
5. N reviews in a `/tmp` worktree and decides.
