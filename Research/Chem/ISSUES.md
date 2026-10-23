# Chem: issues and directions
_Updated 2026-10-23_

## Decisions needed (N)
- **CHEM-D2**: Should environmental waste saturate? — Currently `W` is unbounded and mass-conserving: excretion stops at `W >= 1`, but deaths still return waste, so where the sink is slower than production (`waste_sink` < tax x density) `W` piles up to 5-12, cells can never excrete, and waste acts as uniform aging (lifespan ~ 1/tax). Option: cap `W` at 1 with the excess discarded and booked as an overflow sink, the counterpart of the food `F <= 1` regeneration cap (N's "parallel to food" principle). Recommendation: keep unbounded for now (conserving, and the balance regime is where structure appears); revisit if the aging regime matters.

## Decided
- **CHEM-D1** (2026-10-23): merged to `main` (`1a41c5d`, fast-forward) after a local review and the cloud `/code-review ultra`; findings fixed (CHEM-F1), shared `diffuse()` parked (CHEM-P1).
- **CHEM-D3** (N, 2026-10-23): yes to wgenes, after a short look at global chem (CHEM-N5). Hereditary local chemistries are the longer-term goal; wgenes are the first step: the chemistry stays global and one exchange becomes heritable. Design (to confirm in the dev plan): wgenes use the egene representation (ternary D4 patterns, max signed match `wscore` clamped at 0); excretion `d = disp * (b + (1-b) * wscore/25) * a * (1-E)`, with a floor `b` so that `b = 1` reproduces today's uniform excretion bit-for-bit (the first parity test); `tax_per_wgene` and `mu_wgene` by symmetry with eating; wgenes enter the genome hash and the neutral shadow mutates them (see the egene shadow mismatch). Built as per-compound exchange genes (egenes = uptake of compound 0, wgenes = excretion of compound 1), separate sets from egenes at first. Motivation: eating and excreting then demand different local CA patterns, which can select for temporal alternation, a visible dynamics change driven by heritable function (the activity-to-visual link, `Docs/activity-thoughts.md`).

## Open questions
- **CHEM-Q1**: Does food/waste give the standard (capped) model a refractory medium with propagating waves or spirals? — First look (notebook, 2026-10-23, N=256, one seed): at the legacy-spiral parameter point under the cap, waste produces a sustained oscillating patchy regime near the balance `waste_sink ~ tax x density` (`waste_sink` 0.014-0.015), absent without waste and on either side of the balance. Whether the patches propagate as waves, and whether any rotate, is not yet established.
- **CHEM-Q2**: Does waste keep selection alive at abundance (checkpoint section B / LB6 hypothesis)? — Untested. Needs evolution metrics (dyn_excess, A_new) at high `food_inc` with waste on and off.
- **CHEM-Q3**: In the sink-limited (aging) regime, is there selection for cheaper genomes (lower evolvable tax terms), since lifespan is ~ 1/tax?

## Next directions
- **CHEM-N1**: Wave survey around the balance — campaign (torque): N=512, `waste_sink` x `waste_disp` x `wdiff` x `food_inc`, films, and the spiral detector shared with RXD-N6 (phase singularities of the local food/waste oscillation).
- **CHEM-N2**: Selection at abundance (section B high flank) — campaign: high `food_inc`, waste on/off, dyn_excess / A_new, 3+ seeds.
- **CHEM-N3**: wgenes (evolvable excretion) by symmetry with eating — dev, decided (CHEM-D3); starts after CHEM-N5.
- **CHEM-N5**: Robustness of the balance regime (`waste_sink ~ tax x density`): seeds, N=256/512, a `waste_sink` x `food_inc` grid; the global-chem baseline wgenes will be compared against — short campaign, before CHEM-N3.
- **CHEM-N4**: Milestone 3, the reaction network (FKP-style autocatalysis) and egene -> compound-profile eating — dev, `Docs/plan-chem.md`.

## Fixed
- **CHEM-F1** (2026-10-23): review fixes, `Docs/plan-chem.md` §8. Setters and stamps that kill a cell, replace its occupant or rewrite its store (`stamp_patch`, `set_f_all`, `set_alive*`, `set_w`) now keep the waste-mass balance. Waste can no longer carry energy (an energy-1 product made the tax free). `set_n_chem`, a refused `init` and the controls Restart keep `params()`/recipes in agreement with C (Restart had silently dropped the food/waste chemistry). Colormode 6 shows waste with metabolism off. NaN scan-row cells are dropped. Per-compound `get_chem_env`/`get_chem_priv` were added. Ordinary food/waste runs are unchanged: the notebook outputs are identical.

## Parked
- **CHEM-P1**: One shared `diffuse(field, tmp, d)` for the three copies of the 3x3 periodic box blur (food full and partial passes, compound fields). A separate refactor; it must keep the food goldens bit-exact.
