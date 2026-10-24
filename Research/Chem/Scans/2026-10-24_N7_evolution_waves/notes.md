# CHEM-N7: evolution in vs out of the wave regime (2026-10-24)

**Question (N).** Does evolution differ in the wave regime found by CHEM-N6?

**Runs.** `scan.py`, run locally (torque was overheating). Recipe `2026-10-23_chem_fw_spiralpoint`
with mutation on, N=512, 20000 steps, 4 seeds per condition, 0 errors, about 12 minutes on 10
workers.

| | condition | regime measured (site period, windows ending 2k/6k/10k/15k/20k) |
|---|---|---|
| W1 | one-way, tox 0.4, sink 0.019, wdiff 0.01 | waves in every window, every seed, period 116-117 |
| W2 | two-way, tox 1.0, sink 0.015, wdiff 0.08 | waves in every window, every seed, period 120-121 |
| S1 | as W1 with wdiff 0.08 | no waves, and **waste is harmless**: no toxic deaths, private waste ~0.1 against toxicity 0.4, deaths as in C0 |
| T1 | one-way, tox 1.0, sink 0.019, wdiff 0.01 | **waves early, then lost**: period 128-133 at 2k in all seeds, gone by 6k-15k |
| C0 | waste off (`n_chem=1`) | no waves |

S1 turned out to be effectively vanilla, the same as C0 on every measure, so it is not the
intended "waste-active, no-waves" control. T1 after its waves vanish is the closest thing to that
control.

## Results (seed means; seed spreads are small, see `summary.csv` and `evolution.png`)

| | W1 | W2 | T1 | S1 | C0 |
|---|---|---|---|---|---|
| alive density | 0.325 | 0.347 | 0.373 | 0.383 | 0.383 |
| generation time (1 / births per alive cell per step) | 42 | 52 | 64 | 216 | 217 |
| dyn excess slope, per step | 3.5 | 4.2 | 11 | 83 | 84 |
| dyn excess per generation | 145 | 218 | 709 | 17 900 | 18 100 |
| eg excess slope, per step | 79 | 85 | 100 | 257 | 252 |
| A_new (final) | 47 | 1086 | 60 | 18 600 | 18 900 |
| distinct genomes (mean) | 83 800 | 89 700 | 96 300 | 99 000 | 99 000 |
| **egenes per cell** (starts at 1.0) | **1.38** | **1.28** | 1.20 | 1.02 | 1.02 |
| egene specificity | 7.8 | 7.5 | 5.5 | 7.6 | 7.6 |
| mean intake | 0.048 | 0.045 | 0.044 | 0.042 | 0.042 |

## Reading

1. **Waves change evolution a lot, and in two directions that need separating.**
   - **Genotype persistence above drift collapses in the wave regime.** The activity-excess
     measures and A_new are 5-400 times lower in W1 and W2 than in vanilla, per step and also per
     generation.
     - These measures accumulate the presence of genotypes over time against neutral baselines,
       so they reward genotypes that persist.
     - In the wave regime, mortality is massive and mostly genotype-blind: in W1 about half of
       all deaths are toxic, deaths run at about 3.7 times the vanilla rate, and the generation
       time is 42 steps against 217.
     - A front of waste kills whatever stands in the refractory zone. That weakens selection
       relative to drift and shortens every lineage.
   - **Yet the phenotype moves directionally.** Egenes per cell rise from 1.0 to 1.38 (W1) and
     1.28 (W2), despite the per-egene tax, while vanilla stays at 1.02. Intake is up about 15%.
     Waves select for broader eating (more egenes): plausibly a race to refill after each
     refractory passage, when food is abundant and briefly uncontested.
   - So, in the wave regime: **strong selection on a few traits, short-lived genotypes**.
2. **"Vanilla has high activity" is not "vanilla is innovating".** C0's large A_new and excess
   come with the stasis N accepted on 2026-10-23: long-lived genotypes accumulate cumulative
   activity. These statistics reward persistence, so they must not be read as "waves suppress
   evolution".
   - The measure built for this is the **activity flux** (N, 2026-10-24): the population of
     genomes crossing a low activity band (p20-p30) per tick, the light-blue/teal trace of the
     controls time-series probe. See "Activity flux" below for why it could not be used here.
3. **Evolution may break waves (T1).** T1's waves disappear after 2000-15000 steps in every seed.
   Whether that is evolution (for example egene or LUT change damping the excitability) or a slow
   transient of the medium needs a mutation-off control.
4. **The dominant-genome changes measure is saturated** (about 100 of 101 samples in every
   condition): at N=512 with about 90 000 genomes the top genome is a tiny fraction. It is
   uninformative here.

## Activity flux (`flux.py`, 2026-10-24): two defects, no result

The run (five conditions x 4 seeds, activity updated every tick) gave flux of about 8190 in
**every** condition and seed, although births differ 5-fold. Two problems:

1. **Clipping.** `EvoCA.activity_crossings` collects crossing samples into a buffer of
   `max_n = 65536`. When it fills, the sum stops, and `activity_flux = sum / window` pins at
   65536 / 8 = 8192. Checked with a larger buffer at t = 3000:
   - W1: 72136 samples, flux 9019 (reported 8194);
   - C0: 82624 samples, flux 10333 (reported 8197).

   Nothing signals the truncation. This affects the GUI probe trace too, whenever the band holds
   more than 65536 genome-ticks (large N, high diversity).
2. **Degenerate when every genome is a singleton.** In this recipe D/pop = 1.00: every cell
   carries its own LUT genome. A singleton's activity is its age, and each genome inside the
   band adds 1 per tick, so the unclipped flux is about the number of genomes in one decile,
   D/10 (W1: 9019 / 78045 = 0.116; C0: 10333 / 100673 = 0.103). It measures diversity, not
   clones rising.

The runs were stopped. To use the flux as intended:

- fix the clipping (a quick fix in Python: pass `max_n` >= live buckets x window; or the proper
  fix in C: compute the sum without a buffer);
- measure where genomes have populations, at a lower `mu_lut` (D << pop), or on a coarser
  identity (egene keys, effective-rule species).

## Caveats

- **Confounds.** Population (0.33 against 0.38), food flux and turnover all differ between
  conditions. The per-generation numbers remove the turnover difference only to first order.
- The activity measures' known caveats apply (the shadow models; `Docs/EvoCA_research_directions.md` §2).
- One recipe. Everything here is one point of the food/waste plane.

## Next (proposed)

- **Mutation-off controls** (`mu_lut = mu_egene = mu_egenome = 0`) for W1 and T1: does T1 lose its waves without evolution, and does W1 keep its period?
- **Transplant or competition assays** (`python/patch_transfer.py`): evolved W1 genomes against ancestral ones in the W1 medium and in C0, and evolved C0 genomes in the W1 medium. Adaptation shows up directly as fitness, without depending on the activity statistics.
- **Phenotype trajectories over time** (egene count, specificity, LUT ring level) to see when and how fast waves select.
- A **waste-active no-wave control** chosen on purpose, for example T1 late or a tox and sink point with toxic deaths but no period.

## Infrastructure

`run_sim` has a new optional `step_hook` (`.start(sim)`, `.step(sim, t)`, `.summary()`) for
per-scan extras without copying its step loop. It may not overwrite `run_sim` keys.
`tests/test_run_sim_step_hook.py` checks that the hook sees every step and reaches the result,
that clashes are refused, and that a passive hook leaves the run bit-identical. Dropping the
per-step call or the clash check each fails a test.
