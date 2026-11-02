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
| `resource-driving` (C1) | **merged** (Tier 2) | `c9ecf32` |
| `mu-genome` (A2) | **merged** (Tier 2) | `28e28be` |
| `chem` (E1/E3) | **merged after two fixes** (Tier 2) | `702143a` |
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

## Tier 2 merge record (2026-10-15)

N authorized merging the three Tier-2 branches "if they seem bug free". All
three were merged, `chem` after two fixes.

### Summary

| Branch | Verdict | Merge | Notes |
|---|---|---|---|
| `resource-driving` (C1) | merged | `c9ecf32` | conflict with metric-flux resolved; `food_per` in `_PARAM_KEYS` |
| `mu-genome` (A2) | merged | `28e28be` | storage only; added pytest coverage (branch had none) |
| `chem` (E1/E3) | merged after fixes | `702143a` | two bugs fixed; `n_chem` in `_PARAM_KEYS` and recipes |

- Tests: 68 → **88**, all green after each step.
- `make` with `-Wall`: no warnings after each step.
- Parity probe vs pre-merge `main` (`86d80e1`) after each step: **BAD 0,
  trajectory IDENTICAL**, no new or lost keys. The fixture was viable
  (not extinct, both seeds; final alive density 0.338).
- Step time at N=256 unchanged within noise (about 2.4 ms/step before and
  after; the chem store adds 72 B/cell, about 19 MB at N=512).

### Method

1. Created the `/tmp/evoca_tier2` worktree on `main`.
2. Cherry-picked each one-commit branch in order (`resource-driving`,
   `mu-genome`, `chem`), resolving conflicts and adding tests, then reran the
   build, the suite and the parity probe.
3. Rebuilt the history as one `--no-ff` merge commit per branch and
   fast-forwarded `main`.

### `resource-driving` (C1, `food_per`)

- **Conflict with metric-flux (Phase 2).** The Phase-2 body is now gated on
  `g_step % food_per == 0` and adds `inc = food_inc*food_per`. The flux
  accountant's `g_source_step` and the overflow sink count `inc` on firing
  steps. At `food_per=1`, `inc == food_inc` exactly.
- **Energy identity.** It holds every step at `food_per=6`; a test was added.
- **Reachability.** `food_per` was added to `_PARAM_KEYS`. A `run_sim` test
  checks that the cumulative source equals `food_inc·K·cells·fires`. It
  **fails without the `_PARAM_KEYS` entry**.
- **Feature-on smoke test.** With `run_sim` at `food_per=8`, every dynamic
  scalar moves. At 300 steps the cumulative source is 37 fires × 8, not 300
  steps' worth.
- **Caveat for N (analysis, not a bug).** The per-step `diss_*_step_*` and
  `diss_phi_*` scalars that `run_sim` samples every `sample_every` steps
  alias with the regen period. Use the `*_cum_final` totals. This is recorded
  in RES-N5.

### `mu-genome` (A2, rate-gene storage)

- **Review.** The array is always allocated and freed. It is seeded to
  `mu_lut` with no RNG use, copied at reproduction, and never read by the
  dynamics.
- **Tests.** The branch had none, so `tests/test_mu_lut_gene.py` was added.
  It checks seeding, inertness (a wild gene field leaves the full state
  bit-identical) and inheritance without mutation.
- **Caveat.** `extract_patch`/`stamp_patch` do not carry the gene, so stamped
  cells keep the host slot's value. This is harmless while the gene is
  unread, but it must be added when the gene starts to matter. The
  dead-cell setters also leave the gene in place.

### `chem` (E1/E3 milestone 1)

**Conflicts with metric-flux and lineage-return (tax and eating).** Tax and
intake now go through `chem_cell_add`. The flux `f_before`, the dissipation
sinks and the lineage accumulators are kept.

The branch's own golden-snapshot test, captured before the Tier-1 merges,
still passes bit-for-bit.

**Bug 1: `stamp_patch` bypassed the chem store.**
- *Cause.* `python/evoca_py.py: stamp_patch` writes `f_priv` straight
  through the `evoca_get_f()` pointer. With chem, the sparse store is the
  storage of record, so the next `chem_cell_add` (the tax) restored the
  host cell's old private food into the stamped organism.
- *Effect.* An overwrite stamp onto a host fed at f=0.8 left 64 of 64
  stamped cells alive after one step. On `main` it leaves 0, since stamped
  cells start at f=0 and are taxed to death.
- *Scope.* Invade mode onto dead ground was unaffected, because the store
  was already 0 there.
- *Fix.* New `evoca_chem_sync_from_f()`, called at the end of
  `stamp_patch`.
- *Test.* The regression test **fails with the call removed**.
- *Rule going forward.* Any future Python write through `evoca_get_f()`
  must also call `evoca_chem_sync_from_f()`.

**Bug 2: `set_chem_energy` accepted weights other than 1.**
- *Cause.* Intake, tax, the split and the eating headroom move amounts. A
  weight of 2 on compound 0 therefore let f exceed 1.
- *Effect.* It broke the dissipation energy identity: the residual was about
  2e3 over 50 steps, against 3e-5 at weight 1.
- *Fix.* Milestone 1 now clamps the weight to 1.0 with a stderr notice,
  mirroring `set_n_chem`. Milestones 2–3 must define energy-weighted intake
  before lifting this.
- *Test.* The regression test **fails without the clamp**.

**Parameters.** `n_chem` was added to:
- `_PARAM_KEYS`;
- `params()`, `_DEFAULTS` and the recipe `metaparams_init`/`_final`.

The `run_sim` reachability test (`n_chem=4` reaches the C setter, which
prints its clamp notice) **fails without the `_PARAM_KEYS` entry**.

**Feature-on smoke test.** Only `N_chem=1` exists at milestone 1, so there is
no "on" state beyond the clamp checks above.

### Files relocated

- Merge-test notebooks went from the repo root to `tests/notebooks/`:
  `resource-driving-test.ipynb`, `mu-genome-test.ipynb`,
  `chem-test.ipynb`.
- Each now finds the repo root via the `C/evoca.h` marker and was executed
  with outputs.
- The `resource-driving` notebook no longer writes a PNG into the repo root.
- `chem-test.ipynb` gained a cell for the `stamp_patch` fix.
- `Docs/plan-chem.md` stays in `Docs/`. There is no Chem workstream dir until
  there are results.
- None of the branches wrote to `Runs/` or `Scans/`.

### Deferred / for N

- **Lab notebook snapshot.** `Research/Resource-driving/Resource-driving.ipynb`
  still shows the pre-merge RES-D1/N5/N6 text as a snapshot. `ISSUES.md` and
  the board are updated.
- **C2/L2 (RES-N6) and chem milestone 2 (E2)** were planned to continue on
  their branches. Those branches are now merged and deleted, so each needs a
  new branch from `main`.
