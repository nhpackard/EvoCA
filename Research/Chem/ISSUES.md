# Chem: issues and directions
_Updated 2026-10-23_

## Decisions needed (N)
- **CHEM-D1**: Merge chem milestone 2 (branch `chem`) into `main`? — Additive: at `n_chem=1` none of the new code runs and every golden (milestone-1 parity, food-cap, food-thresh fingerprints) is unchanged. Pre-merge test notebook `tests/notebooks/chem-m2-test.ipynb`; `/code-review ultra` is worth running first (new C mechanism). Recommendation: merge after review.
- **CHEM-D2**: Should environmental waste saturate? — Currently `W` is unbounded and mass-conserving: excretion stops at `W >= 1`, but deaths still return waste, so where the sink is slower than production (`waste_sink` < tax x density) `W` piles up to 5-12, cells can never excrete, and waste acts as uniform aging (lifespan ~ 1/tax). Option: cap `W` at 1 with the excess discarded and booked as an overflow sink, the counterpart of the food `F <= 1` regeneration cap (N's "parallel to food" principle). Recommendation: keep unbounded for now (conserving, and the balance regime is where structure appears); revisit if the aging regime matters.
- **CHEM-D3**: wgenes — evolvable excretion (N's open question in the checkpoint)? — Without them excretion is uniform (`waste_disp`), so selection can act on waste only through the evolvable tax terms (`tax_per_egene`, `tax_lut`) and through where cells live. Recommendation: first characterize the no-wgene model (CHEM-N1/N2).

## Open questions
- **CHEM-Q1**: Does food/waste give the standard (capped) model a refractory medium with propagating waves or spirals? — First look (notebook, 2026-10-23, N=256, one seed): at the legacy-spiral parameter point under the cap, waste produces a sustained oscillating patchy regime near the balance `waste_sink ~ tax x density` (`waste_sink` 0.014-0.015), absent without waste and on either side of the balance. Whether the patches propagate as waves, and whether any rotate, is not yet established.
- **CHEM-Q2**: Does waste keep selection alive at abundance (checkpoint section B / LB6 hypothesis)? — Untested. Needs evolution metrics (dyn_excess, A_new) at high `food_inc` with waste on and off.
- **CHEM-Q3**: In the sink-limited (aging) regime, is there selection for cheaper genomes (lower evolvable tax terms), since lifespan is ~ 1/tax?

## Next directions
- **CHEM-N1**: Wave survey around the balance — campaign (torque): N=512, `waste_sink` x `waste_disp` x `wdiff` x `food_inc`, films, and the spiral detector shared with RXD-N6 (phase singularities of the local food/waste oscillation).
- **CHEM-N2**: Selection at abundance (section B high flank) — campaign: high `food_inc`, waste on/off, dyn_excess / A_new, 3+ seeds.
- **CHEM-N3**: wgenes (evolvable excretion) by symmetry with eating — dev, after CHEM-D3.
- **CHEM-N4**: Milestone 3, the reaction network (FKP-style autocatalysis) and egene -> compound-profile eating — dev, `Docs/plan-chem.md`.

## Fixed

## Parked
