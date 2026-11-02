# Chem: issues and directions
_Updated 2026-10-23_

## Plans and design docs
- [`plans/plan-chem.md`](plans/plan-chem.md): The `chem` branch plan: N_chem machinery (milestone 1, done), food/waste (milestone 2, done, on `main`). Milestone 3, the reaction network, is open as CHEM-N4. Review fixes in §8 (CHEM-F1).
- [`plans/2026-10-23_chem-model-v2.md`](plans/2026-10-23_chem-model-v2.md): Chem model v2 (`chemca`), a metabolism-based agent: design only, Track B (CHEM-D5). Its implementation is on the unmerged `chemca` branch.

## Decisions needed (N)
- **CHEM-D2**: Should environmental waste saturate? — Currently `W` is unbounded and mass-conserving: excretion stops at `W >= 1`, but deaths still return waste, so where the sink is slower than production (`waste_sink` < tax x density) `W` piles up to 5-12, cells can never excrete, and waste acts as uniform aging (lifespan ~ 1/tax). Option: cap `W` at 1 with the excess discarded and booked as an overflow sink, the counterpart of the food `F <= 1` regeneration cap (N's "parallel to food" principle). Recommendation: keep unbounded for now (conserving, and the balance regime is where structure appears); revisit if the aging regime matters.

- **CHEM-D5**: Chem model v2 (`chemca`), a metabolism-based agent as a separate model (`Research/Chem/plans/2026-10-23_chem-model-v2.md`, §13 lists 9 decisions). **Status (N, 2026-10-24): design only; Track B, later.** Current work stays in vanilla EvoCA food/waste (Track A): its dynamics and their evolution. An evolutionary Gray-Scott (three species: F + 2A -> 3A and A -> P inside agents, P -> F by light) is agreed as the reference starting point for chemca. CHEM-N3 (wgenes) is not paused: wgenes belong to vanilla.

## Decided
- **CHEM-D4** (N, 2026-10-23): two-way (passive membrane) waste exchange is adopted **as an option**: metaparameter `waste_exchange` (0 = one-way excretion, the default, so every recipe is unchanged; 1 = two-way `d = waste_disp * (w - W)`), a checkbox in the controls. Bit-identical to the experiment patch; `waste_exchange=0` is bit-identical to the previous `main`. Tests `tests/test_waste_exchange.py` and the two-way cases of the waste-mass and energy identities. The wave/spiral survey under it is CHEM-N1.
- **CHEM-D1** (2026-10-23): merged to `main` (`1a41c5d`, fast-forward) after a local review and the cloud `/code-review ultra`; findings fixed (CHEM-F1), shared `diffuse()` parked (CHEM-P1).
- **CHEM-D3** (N, 2026-10-23): yes to wgenes, after a short look at global chem (CHEM-N5, done). Hereditary local chemistries are the longer-term goal; wgenes are the first step: the chemistry stays global and one exchange becomes heritable. Design (to confirm in the dev plan): wgenes use the egene representation (ternary D4 patterns, max signed match `wscore` clamped at 0); excretion `d = disp * (b + (1-b) * wscore/25) * a * (1-E)`, with a floor `b` so that `b = 1` reproduces today's uniform excretion bit-for-bit (the first parity test); `tax_per_wgene` and `mu_wgene` by symmetry with eating; wgenes enter the genome hash and the neutral shadow mutates them (see the egene shadow mismatch). Built as per-compound exchange genes (egenes = uptake of compound 0, wgenes = excretion of compound 1), separate sets from egenes at first. Motivation: eating and excreting then demand different local CA patterns, which can select for temporal alternation, a visible dynamics change driven by heritable function (the activity-to-visual link, `Research/Activity-metrics/plans/activity-thoughts.md`).

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
  - **CHEM-N6** (2026-10-24): **excitation waves**, also under **one-way** exchange at `waste_tox` 0.4 / `waste_sink` 0.019 (N's interactive point), at N = 256-1024.
    - Space-time plots show fronts leaving sources, traveling at about 1-2 cells/step and annihilating on collision.
    - Sites fire with one period, 115-120 steps.
    - The fronts are irregular in space (many sources), so the spatial spectrum has no ring, only a knee at the patch scale.
    - The N5 "no propagation" statement holds only at `waste_tox` 1.0; toxicity 0.4 makes the one-way medium excitable.
    - Spiral versus target is still unresolved.
- **CHEM-Q4**: What would make the two-way medium's phase coherent below the ~200-cell wavelength, so spirals could be resolved?
  - The colonies are a scattered ~33% speckle, so the local cycle phase is noisy.
  - Candidates:
    - denser colonies behind a sharper front, which the egene/LUT, `m_scale` and `food_inc` set;
    - evolution itself (wgenes, CHEM-D3, may select tighter colony cycles);
    - a coarser readout.
  - Whether spirals matter enough to pursue is N's call.
- **CHEM-Q5**: How does the wave regime shape evolution? CHEM-N7 (2026-10-24, N=512, 20000 steps, 4 seeds):
  - **Genotype persistence above drift collapses in waves.** Activity excess and A_new are 5-400 times lower than vanilla, per step and per generation. Mortality there is massive and largely genotype-blind: half of all deaths are toxic, and the generation time is 42 steps against 217.
  - **The phenotype moves directionally.** Egenes per cell rise from 1.0 to 1.38 (one-way waves) and 1.28 (two-way), against 1.02 in vanilla, and intake is about 15% higher. So: strong selection on a few traits, short-lived genotypes.
  - Vanilla's high activity is persistence (stasis), not innovation; these statistics cannot tell the two apart.
  - One-way at tox 1.0 (T1) has waves early and loses them by 6000-15000 steps in every seed: evolution or a transient? Needs a mutation-off control.
  - wdiff 0.08 at tox 0.4 (S1) makes waste harmless, so it behaves like vanilla.
  - Next: mutation-off controls, transplant or competition assays (`patch_transfer.py`), phenotype trajectories.
- **CHEM-Q2**: Does waste keep selection alive at abundance (checkpoint section B / LB6 hypothesis)? — Untested. Needs evolution metrics (dyn_excess, A_new) at high `food_inc` with waste on and off.
- **CHEM-Q3**: In the sink-limited (aging) regime, is there selection for cheaper genomes (lower evolvable tax terms), since lifespan is ~ 1/tax?

## Next directions

- **CHEM-N2**: Selection at abundance (section B high flank) — campaign: high `food_inc`, waste on/off, dyn_excess / A_new, 3+ seeds.
- **CHEM-N3**: wgenes (evolvable excretion) by symmetry with eating — dev, decided (CHEM-D3); starts after CHEM-N5.
- **CHEM-N4**: Milestone 3, the reaction network (FKP-style autocatalysis) and egene -> compound-profile eating — dev, `Research/Chem/plans/plan-chem.md`.

## Fixed
- **CHEM-F2** (2026-10-23, N):
  - Colormode 6 (food-waste) is built like env-food, and colormode 7 (env-waste) is new.
  - The Run/Pause label fix.
  - `waste_tox` is a metaparameter (default 1.0, the old fixed value).
  - Sliders for `waste_sink`, `wdiff`, `waste_disp` and `waste_tox`.
- **CHEM-N7 done** (2026-10-24): `Scans/2026-10-24_N7_evolution_waves`, run locally; results in CHEM-Q5. `run_sim` gained an optional `step_hook`.
- **CHEM-N6 done** (2026-10-24): `Scans/2026-10-24_N6_layered_fronts`; results in CHEM-Q1. New statistics in `python/evoca_waves.py`: `ring_spectrum` (with ring contrast), `kw_spectrum`, `site_period`. Measures tried and dropped are listed in the notes.
- **CHEM-N1 done** (2026-10-23): `Scans/2026-10-23_N1_twoway_waves`; results in CHEM-Q1. Detector `python/evoca_waves.py`.
- **CHEM-N5 done** (2026-10-23): `Scans/2026-10-23_N5_global_chem`; results in CHEM-Q1.
- **CHEM-F1** (2026-10-23): review fixes, `Research/Chem/plans/plan-chem.md` §8. Setters and stamps that kill a cell, replace its occupant or rewrite its store (`stamp_patch`, `set_f_all`, `set_alive*`, `set_w`) now keep the waste-mass balance. Waste can no longer carry energy (an energy-1 product made the tax free). `set_n_chem`, a refused `init` and the controls Restart keep `params()`/recipes in agreement with C (Restart had silently dropped the food/waste chemistry). Colormode 6 shows waste with metabolism off. NaN scan-row cells are dropped. Per-compound `get_chem_env`/`get_chem_priv` were added. Ordinary food/waste runs are unchanged: the notebook outputs are identical.

## Parked
- **CHEM-P1**: One shared `diffuse(field, tmp, d)` for the three copies of the 3x3 periodic box blur (food full and partial passes, compound fields). A separate refactor; it must keep the food goldens bit-exact.
