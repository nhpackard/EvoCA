# Findings from the workstream seeding pass (2026-10-15)

Six Sonnet agents seeded `Research/<Workstream>/` (lab notebook +
`ISSUES.md`); all six notebooks execute clean. While building their example
runs, several of them hit problems in **merged code or past methodology**.
This file collects the cross-cutting ones and marks what Claude verified
independently.

Per-workstream detail lives in each `ISSUES.md` (the source of truth); the
board's generated section mirrors them.

## A. Defects and degeneracies in merged code

| # | Finding | Found by | Verified | Status |
|---|---|---|---|---|
| A1 | `run_sim`/`build_sim`/`bifurcation` **silently dropped `tax_ring`**, because it was missing from `_PARAM_KEYS`. | Ring | yes | **fixed** in `d15d0ad`. The new `init_kwargs()` helper warns on any key it would drop. Regression test goes through `run_sim` and was shown to fail without the fix. No past scan set `tax_ring`, so no recorded result is affected. |
| A2 | **The neutral shadow perturbs the dynamics.** Its neutral mutations draw from the same C RNG stream (`rng_next`, `C/evoca.c` ~2556), so `shadow=True` and `shadow=False` with the same seed give different trajectories. Sampling activity does not perturb anything. | Lineage, Resource (independently) | yes. Same seed: no shadow `d62cffd337`, shadow `1c4457f6dd`; sampling has no effect. | **Open, N's decision.** The distributions are presumably unchanged (different draws of the same process), so past statistical results stand, and all arms of every campaign ran shadow-on. But a "same seed, shadow on vs off" comparison is *not paired*, and frozen-GoL viability flips with it (LIN). The fix is a separate RNG stream for the shadow, making it a pure observer. That changes every shadow-on trajectory relative to the past, so it is not a no-op merge. **Update 2026-10-15:** N approved; implemented on branch `shadow-rng` (dedicated shadow stream, tests in `tests/test_shadow_rng.py`, notebook `tests/notebooks/shadow-rng-test.ipynb`), awaiting N's review. |
| A3 | **`eff_activity` (F1) barely coarsens.** Its key masks each LUT by the *population-union* used-mask, which saturates: 236/250 bits at step 1 and 250/250 by step 500. Effective-rule species therefore come out at 95–99% of alive cells, essentially full-genome species. | Activity | yes (N=96 pure-evo config) | **Open.** The union key is degenerate in evolving populations, the same pattern as #8's binary indicator. The per-cell used-mask (opt-in, W1-c) is the variant that can actually coarsen; it needs promoting to a species key. This weakens the review §7 expectation of "A_new on eff-species". |
| A4 | **`A_new` (G1) is blind in egene-only runs.** The real side is keyed by the whole genome, while the shadow is LUT-only, so A_new=0 in EGENE_only even though eg-excess > 0. The frozen control gives A_new=1. | Activity | not independently. It is consistent with the known LUT-only shadow mismatch. | **Open** → ACT-D1. |
| A5 | **The lineage reducer counts nominal tax.** The dissipation accountant counts tax actually removed; the two differ by 20–23% in runs with tax deaths. | Lineage | consistent with the merge-time note in `01-branch-review.md` | **Open** → LIN-D3 (add an actual-tax column before H1 relies on `net`). |
| A6 | **Hash-keyed "lineages" are genotypes.** At productive `mu_lut` they are near-singletons (no hash reaches 10 births at 0.03). Raw net return is also pinned by pool share. | Lineage | not independently | **Open** → LIN-D2 (lineage unit: clade roll-up, or start H1 at low `mu_lut`). |
| A7 | **W1-b's "record p99 alongside p90" was not implemented** in `metric-anew`. | Activity | — | Open (ACT ISSUES). |
| A8 | `dyn_report` labels a config `silent-churn` while its dyn excess is strongly positive (106/500 buckets with z>3). | Activity | not independently | Open (ACT ISSUES). |

## B. Methodology findings that change plans

- **Activity:** the 5-key pure-evo composite tracks eg-excess (r=0.69), not dyn-excess (0.41). Its top-5 shares no config with the dyn-excess top-5, and eg and dyn winners sit in different parameter corners. #1 vs #2 also confounds egene mutation with egene cost, because #2 zeroed `tax_per_egene`. → ACT-D1 recommends a **two-channel (dyn, eg) objective**, with Pareto selection (ACT-D2) instead of a scalar composite.
- **Resource-driving:** below `food_inc/tax`=1, final density ≈ `food_inc/tax` (pilot). The rising flank of any activity curve is therefore confounded with occupancy. Φ ≈ max(tax, food_inc) is a weak x-axis inside the food-limited regime. Both findings are built into the RES-N1 campaign spec.
- **Reaction-diffusion:** much of the 04-27 spatial structure is the halfplane-start transient (correlation length peaks within ~600 steps). `waves_in_blank_half` reproduces only from its `metaparams_final`. #4's ten bases all came from the `initial` scan. → RXD-D2 asks whether halfplane is the RD protocol; RXD-Q2 asks to extend #4. In the cone/box preview, the box shares 18/20 with the composite top-20 while the wedge shares 12/20, so the wedge is what changes the answer (RXD-N1).
- **Criticality:** the first `food_inc` diagram shows a stochastic extinction edge, a large irregular boom-bust cycle (it persists to N=128 in one seed), amplitude collapse and saturation, with no period doubling. The extrema reducer is unusable on noisy series (CRIT-N6). Pop=0 is absorbing, so the extinction edge may be an absorbing-state transition (CRIT-Q3).
- **Ring:** without tax, a GoL start climbs to more than 90% ring-level 3 under mutation within ~500 steps. That may be mutational entropy rather than selection, so it needs a neutral baseline (RING-Q1/N4).
- **Lineage:** #8's degeneracy reproduces exactly (P(ΔLUT)=1−exp(−250·mu_lut)). LIN-D1 recommends Option A (Hamming magnitudes) with `restricted_mu=1`, plus a synthetic-coupling power check first.

## C. Housekeeping notes

- `Research/Activity-metrics/Runs/2026-04-27_high_excess_activity_slope.evoca` is a copy of push_edges cfg 98 (`top_9`), not a high-slope winner. It is left unrenamed and flagged in ACT ISSUES.
- Campaign #3 has no `scan_config.json`/`notes.md`, and #8 has no `results.csv`.
- The RD lab notebook is 7.9 MB with its gallery. If the repo grows uncomfortable, a future option is to keep gallery images as separate PNGs.
- Notebook "example runs" are labeled as such throughout. None of them is a campaign result.

## D. Decisions now on N's desk (in suggested order)

1. **A2:** a separate RNG stream for the shadow? (Recommended yes, as its own branch with a test notebook, since it breaks bit-parity with past shadow-on runs.)
2. **ACT-D1 / A3 / A4:** the primary activity objective: a two-channel (dyn, eg) vector on fixed-space shadows, with A_new and eff-species kept as diagnostics until fixed. This replaces the review §7 expectation.
3. **Tier-2 gates:** `resource-driving` (RES-D1; unblocks C3, C2/L2), `mu-genome`, `chem`.
4. **LIN-D1/D2/D3:** #8 redesign, the H1 lineage unit, the actual-tax column.
5. **RXD-D1/D2, CRIT-D1/D2, RING-D1/D3, RES-D2:** campaign scoping (see each ISSUES.md).
