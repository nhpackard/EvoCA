# [CHEM-N5](../../ISSUES.md): global-chem baseline (2026-10-23)

**Why.** Before wgenes ([CHEM-D3](../../ISSUES.md)) there are two questions. Is the three-state oscillating regime at the waste balance (`waste_sink ~ tax x density`) robust? And does the food/waste model give the refractory medium that spirals need (N's question, 2026-10-23; [RXD-Q6](../../../Reaction-diffusion/ISSUES.md))? So far there is one seed, at N=256.

**Design** (driver `scan.py`; recipe `Runs/2026-10-23_chem_fw_spiralpoint.evoca`, standard rule, mutation on, halfplane start, 6000 steps, statistics over steps 3000-6000):
- **A, robustness (N=256, 99 runs):** `waste_sink` 0.010-0.020 in steps of 0.001, x `food_inc` {0.013, 0.016, 0.020}, x seeds 0-2.
- **B, wave probe (N=512, 60 runs):** `wdiff` {0.01, 0.03, 0.08, 0.2, 0.5}, x `waste_sink` {0.013-0.016}, x seeds 0-2, at `food_inc` 0.016. Seed-0 frames of F, W and alive at t = 4000, 5000 and 6000 go to `frames/` (gitignored).

**Hypothesis for B.**
- A slowly diffusing inhibitor (`wdiff` < `gdiff` = 0.08) leaves a sharp refractory tail behind a front, which gives waves.
- A fast one spreads ahead of the front (lateral inhibition), which gives stationary patches.

**Caveat: the refractory zone is soft.** A lineage can live in waste-saturated ground (E >= 1, so no excretion) if it reproduces fast enough. Each split halves private waste, and metabolism adds `tax` per step, so the steady load is about 2 x tax x T_rep. Survival needs T_rep < 1/(2 tax), about 14 steps at tax 0.035. Behind a front the food is grazed, so reproduction is slow and the zone acts as refractory. A fast lineage could still creep in.

**Detector.**
- Phase-singularity count on 16x16-block means of (F, W), measured as absolute deviations.
- Plaquettes with low local cycle amplitude are ignored.
- Cores must persist (same sign, within 1 coarse cell) over 60 steps, about half the ~120-step period.
- The chance level is the earlier frame shifted by half the lattice.
- `n_sing_rotating` additionally requires the phase to advance around the core. Under a bulk oscillation, noise defects also pass that test, so the reliable signal is `n_sing_persist` well above `n_sing_persist_null`.

Calibration on built test fields:

| Field | Result |
|---|---|
| One planted rotating spiral | found; it persists and rotates |
| The same spiral frozen | 0 rotating |
| A static mosaic of random patches | 118 persistent, 0 rotating |
| Small-amplitude noise | 0 singularities |

Smoke test at t <= 1500 (balance, N=512, `wdiff` 0.03): persistence 108 against a chance level of 82. That is no cores, matching the frame by eye: patches of about 50 cells, no fronts.

**Results** (159 runs, 0 errors, 0 extinctions; analysis and figure in `Chem.ipynb`, "No spirals yet").
- **The balance rule holds, as a knife edge.**
  - The structured regime (largest alive std, period 100-140) sits at the last `waste_sink` step below the measured waste production:

    | `food_inc` | Structured at `waste_sink` | Measured production |
    |---|---|---|
    | 0.013 | 0.011-0.012 | 0.0125 |
    | 0.016 | 0.015 | 0.0153 |
    | 0.020 | 0.019 | 0.0192 |

  - All 3 seeds agree.
  - Further below, waste piles up (mean W up to 5, 20 and 39 at the three food supplies) and nearly every death is toxic.
  - One step above, waste never matters.
- **Part B: no spiral cores.**
  - Persistent defects run 1.04-1.69x chance, with the highest ratios in the pile-up regime, where nothing propagates.
  - The frames show patches, with no fronts or rotation.
  - `wdiff` 0.01-0.03 keeps the regime one step further (`waste_sink` 0.016).
- **Follow-ups.**
  - The toxicity sweep is in `Chem.ipynb`.
  - The two-way exchange experiment is in `../../Experiments/2026-10-23_twoway_exchange/`.
