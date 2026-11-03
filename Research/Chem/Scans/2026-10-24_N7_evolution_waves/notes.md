# [CHEM-N7](../../ISSUES.md): evolution in vs out of the wave regime (2026-10-24)

**Question (N).** Does evolution differ in the wave regime found by [CHEM-N6](../../ISSUES.md)?

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
     controls time-series probe. See "Activity flux" below.
3. **Evolution may break waves (T1).** T1's waves disappear after 2000-15000 steps in every seed.
   Whether that is evolution (for example egene or LUT change damping the excitability) or a slow
   transient of the medium needs a mutation-off control.
4. **The dominant-genome changes measure is saturated** (about 100 of 101 samples in every
   condition): at N=512 with about 90 000 genomes the top genome is a tiny fraction. It is
   uninformative here.

## Activity flux (`flux.py`)

The first flux run (2026-10-24) was discarded: the flux summed a 65536-sample buffer and read
8192 in every condition. The cap is fixed (2026-11-03: `evoca_activity_flux_sum` sums without a
buffer, and the activity-table readers no longer truncate; `tests/test_activity_flux.py`). The
N7 evolution measures above were not affected: their table read (200 000) holds the whole
table at N = 512 (about 100 000 entries, checked).

A second limit is not a bug: at the recipe's `mu_lut` = 0.03 nearly every birth is a new
genome (D/pop = 1), and a singleton's activity is its age, so the flux then counts genomes
aging through the band (about D/10). The rerun therefore records the clonal part of the flux
(buckets with pop >= 2) and adds a `mu_lut` sweep.

### Activity flux, rerun (2026-11-03, `flux.csv`, `flux.log`)

20000 steps, second half sampled every 8 ticks, run locally (12 workers, 11-38 min per run).
Means over seeds (A: 4 seeds, B: 2); D = distinct genomes.

| part | cond | mu_lut | alive | births/tick | flux | flux / D | clonal frac | clonal flux | D/pop | top genome |
|---|---|---|---|---|---|---|---|---|---|---|
| A | C0 | 0.03 | 0.39 | 469 | 10509 | 0.10 | 0.002 | 25 | 1.00 | 2 |
| A | S1 | 0.03 | 0.39 | 470 | 10509 | 0.10 | 0.002 | 25 | 1.00 | 2 |
| A | T1 | 0.03 | 0.38 | 1561 | 11676 | 0.12 | 0.001 | 9 | 1.00 | 2 |
| A | W1 | 0.03 | 0.33 | 2054 | 10318 | 0.12 | 0.001 | 8 | 1.00 | 2 |
| A | W2 | 0.03 | 0.35 | 1756 | 10889 | 0.12 | 0.001 | 8 | 1.00 | 2 |
| B | C0 | 0.01 | 0.39 | 464 | 10294 | 0.11 | 0.09 | 976 | 0.93 | 5 |
| B | C0 | 0.003 | 0.39 | 490 | 8840 | 0.12 | 0.30 | 2641 | 0.70 | 11 |
| B | C0 | 0.001 | 0.40 | 464 | 6771 | 0.15 | 0.49 | 3348 | 0.45 | 25 |
| B | W1 | 0.01 | 0.33 | 2058 | 10289 | 0.13 | 0.09 | 946 | 0.93 | 5 |
| B | W1 | 0.003 | 0.34 | 2074 | 8956 | 0.17 | 0.47 | 4174 | 0.61 | 19 |
| B | W1 | 0.001 | 0.37 | 2174 | 7046 | 0.21 | 0.70 | 4901 | 0.35 | 50 |

Seed spread of the flux is small (sd across seeds <= 64).

1. **At the recipe's mutation rate the uncapped flux is a diversity count.** It is 0.10-0.12 D in
   every condition, clones carry 0.1-0.2% of it, and the five conditions differ by at most 13%,
   though turnover differs 4.4-fold. It cannot compare the conditions.
2. **Lowering `mu_lut` makes it a clone measure.** At 0.001 clones carry half (C0) to 70% (W1)
   of the flux, and the top genome reaches 25-50 cells. 0.01 is not low enough (9%).
3. **In the wave regime more of the flux is clonal**, at every `mu_lut` below 0.01: W1 0.47 and
   0.70 vs C0 0.30 and 0.49, with twice the top-genome size. Per tick the clonal flux is 1.5x
   C0's; per birth it is a third of C0's (2.3 vs 7.2 at 0.001), because W1 turns over 4.7x
   faster. Which normalization is right is the same per-step vs per-generation question as
   above; neither says whether the rise is selection.
4. **Missing: a neutral reference.** Clones also rise through the band by drift. The test of
   selection is the same clonal flux on the neutral shadow (identical births and deaths, random
   parents), whose table the fix now reads in full: real clonal flux above the shadow's is
   genomes rising faster than drift.
5. Lowering `mu_lut` barely changes the ecology (alive and births within 10%), so the low-mu
   runs are a fair comparison point.

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
