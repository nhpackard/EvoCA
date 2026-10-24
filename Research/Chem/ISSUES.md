# Chem: issues and directions
_Updated 2026-10-23_

## Decisions needed (N)
- **CHEM-D2**: Should environmental waste saturate? — Currently `W` is unbounded and mass-conserving: excretion stops at `W >= 1`, but deaths still return waste, so where the sink is slower than production (`waste_sink` < tax x density) `W` piles up to 5-12, cells can never excrete, and waste acts as uniform aging (lifespan ~ 1/tax). Option: cap `W` at 1 with the excess discarded and booked as an overflow sink, the counterpart of the food `F <= 1` regeneration cap (N's "parallel to food" principle). Recommendation: keep unbounded for now (conserving, and the balance regime is where structure appears); revisit if the aging regime matters.

- **CHEM-D5**: Chem model v2, a metabolism-based agent as a separate model (`Docs/2026-10-23_chem-model-v2.md`, §13 lists 9 decisions). Core points: under pinned decision B a genetic food->waste rate evolves to 0, so v2 needs an energy currency produced by reactions; closed mass (F + W) with W -> F regeneration by a light field; autocatalysis F + W -> 2W offered as a gene choice to make waste needed; mass-action rates on a discrete log grid, reaction genes (A, B, C, k) in the egene slot layout; separate `C/chemca.c` with a shared `core.c` (retires CHEM-P1), recipe `model` key, `controls.py` moved onto a sim interface. Conflicts with CHEM-D3/N3 (wgenes are CA-pattern genes; recommendation: pause N3 pending D5). Also proposes recalibrating the wave statistics against N's eye before any ranking or Pareto search (the one-way tox 0.4 / sink 0.019 point was never measured by them).

## Decided
- **CHEM-D4** (N, 2026-10-23): two-way (passive membrane) waste exchange is adopted **as an option**: metaparameter `waste_exchange` (0 = one-way excretion, the default, so every recipe is unchanged; 1 = two-way `d = waste_disp * (w - W)`), a checkbox in the controls. Bit-identical to the experiment patch; `waste_exchange=0` is bit-identical to the previous `main`. Tests `tests/test_waste_exchange.py` and the two-way cases of the waste-mass and energy identities. The wave/spiral survey under it is CHEM-N1.
- **CHEM-D1** (2026-10-23): merged to `main` (`1a41c5d`, fast-forward) after a local review and the cloud `/code-review ultra`; findings fixed (CHEM-F1), shared `diffuse()` parked (CHEM-P1).
- **CHEM-D3** (N, 2026-10-23): yes to wgenes, after a short look at global chem (CHEM-N5, done). Hereditary local chemistries are the longer-term goal; wgenes are the first step: the chemistry stays global and one exchange becomes heritable. Design (to confirm in the dev plan): wgenes use the egene representation (ternary D4 patterns, max signed match `wscore` clamped at 0); excretion `d = disp * (b + (1-b) * wscore/25) * a * (1-E)`, with a floor `b` so that `b = 1` reproduces today's uniform excretion bit-for-bit (the first parity test); `tax_per_wgene` and `mu_wgene` by symmetry with eating; wgenes enter the genome hash and the neutral shadow mutates them (see the egene shadow mismatch). Built as per-compound exchange genes (egenes = uptake of compound 0, wgenes = excretion of compound 1), separate sets from egenes at first. Motivation: eating and excreting then demand different local CA patterns, which can select for temporal alternation, a visible dynamics change driven by heritable function (the activity-to-visual link, `Docs/activity-thoughts.md`).

## Open questions
- **CHEM-Q1**: Does food/waste give the standard (capped) model a refractory medium with propagating waves or spirals?
  - **First look** (2026-10-23, N=256, one seed): a sustained oscillating patchy regime near the balance.
  - **CHEM-N5** (2026-10-23, campaign): under the one-way exchange rule, **no**.
    - The balance regime is reproducible but only 0.001-0.002 wide in `waste_sink`.
    - Sites cycle locally without propagating, because the refractory zone is soft.
    - Lowering `waste_tox` hardens the zone everywhere at once, so it gives no contrast (toxicity sweep in the notebook).
  - **Two-way exchange** (CHEM-D4; CHEM-N1 campaign, 2026-10-23) gives **traveling waves with stationary sources** across a broad region, 3 seeds each.
    - Wavelength about 200 cells; R about 0.1, against 0.4 for one-way.
    - Defect survival over 4 periods reaches about twice the one-way level.
    - **No spirals resolved in 2D**: the phase field is noisy at 32-64 cells.
- **CHEM-Q4**: What would make the two-way medium's phase coherent below the ~200-cell wavelength, so spirals could be resolved?
  - The colonies are a scattered ~33% speckle, so the local cycle phase is noisy.
  - Candidates:
    - denser colonies behind a sharper front, which the egene/LUT, `m_scale` and `food_inc` set;
    - evolution itself (wgenes, CHEM-D3, may select tighter colony cycles);
    - a coarser readout.
  - Whether spirals matter enough to pursue is N's call.
- **CHEM-Q2**: Does waste keep selection alive at abundance (checkpoint section B / LB6 hypothesis)? — Untested. Needs evolution metrics (dyn_excess, A_new) at high `food_inc` with waste on and off.
- **CHEM-Q3**: In the sink-limited (aging) regime, is there selection for cheaper genomes (lower evolvable tax terms), since lifespan is ~ 1/tax?

## Next directions

- **CHEM-N2**: Selection at abundance (section B high flank) — campaign: high `food_inc`, waste on/off, dyn_excess / A_new, 3+ seeds.
- **CHEM-N3**: wgenes (evolvable excretion) by symmetry with eating — dev, decided (CHEM-D3); starts after CHEM-N5.
- **CHEM-N4**: Milestone 3, the reaction network (FKP-style autocatalysis) and egene -> compound-profile eating — dev, `Docs/plan-chem.md`.

## Fixed
- **CHEM-F2** (2026-10-23, N):
  - Colormode 6 (food-waste) is built like env-food, and colormode 7 (env-waste) is new.
  - The Run/Pause label fix.
  - `waste_tox` is a metaparameter (default 1.0, the old fixed value).
  - Sliders for `waste_sink`, `wdiff`, `waste_disp` and `waste_tox`.
- **CHEM-N1 done** (2026-10-23): `Scans/2026-10-23_N1_twoway_waves`; results in CHEM-Q1. Detector `python/evoca_waves.py`.
- **CHEM-N5 done** (2026-10-23): `Scans/2026-10-23_N5_global_chem`; results in CHEM-Q1.
- **CHEM-F1** (2026-10-23): review fixes, `Docs/plan-chem.md` §8. Setters and stamps that kill a cell, replace its occupant or rewrite its store (`stamp_patch`, `set_f_all`, `set_alive*`, `set_w`) now keep the waste-mass balance. Waste can no longer carry energy (an energy-1 product made the tax free). `set_n_chem`, a refused `init` and the controls Restart keep `params()`/recipes in agreement with C (Restart had silently dropped the food/waste chemistry). Colormode 6 shows waste with metabolism off. NaN scan-row cells are dropped. Per-compound `get_chem_env`/`get_chem_priv` were added. Ordinary food/waste runs are unchanged: the notebook outputs are identical.

## Parked
- **CHEM-P1**: One shared `diffuse(field, tmp, d)` for the three copies of the 3x3 periodic box blur (food full and partial passes, compound fields). A separate refactor; it must keep the food goldens bit-exact.
