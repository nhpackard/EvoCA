# Scan 2026-04-27 — push gdiff and m_scale beyond refined-scan edges

## Question

The 2026-04-27_refined scan top configs sat with gdiff = 0.08 and
m_scale = 0.8 — both the *upper limits* of the refined search ranges,
suggesting the optimum is at or beyond. This scan extends both axes:
gdiff up to 0.20, m_scale up to 1.2. We also drop sample_every from 100
to 25 so `unique_top_genomes` (which maxed out at 51/51 in the previous
two scans) becomes a useful discriminator again — should now max at
~200/200.

## Search axes

```
food_inc      ∈ {0.008, 0.010, 0.013, 0.018}        (4)
tax           ∈ {0.030, 0.035, 0.040}                (3)
gdiff         ∈ {0.06, 0.08, 0.12, 0.15, 0.20}       (5)  ← extended up
mu_lut        ∈ {0.001, 0.003, 0.010, 0.030}         (4)
mu_egenome    ∈ {0.001, 0.010}                       (2)
m_scale       ∈ {0.7, 0.8, 1.0, 1.2}                 (4)  ← extended up
restricted_mu ∈ {True, False}                        (2)
```

Full grid = 7,680 combos. Random sample of 250.

Fixed: N=256, n_steps=5000, **sample_every=25** (was 100), init='halfplane',
seed=0, shadow=True.

## Diminishing-returns question

Compare against `2026-04-27_refined`:
- Are gdiff > 0.08 and m_scale > 0.8 actually better, or did refined
  hit the corner?
- Does `unique_top_genomes` (now able to range up to 200) reveal new
  discrimination among the top configs?

## Results — 179 s on 12 cores

- 250 configs run, **245 candidates** after extinct/saturated filtering
- Yield is **98 %**, up from 50 % (refined) and 17 % (initial). The
  extended `gdiff/m_scale` region is overwhelmingly survivable — we
  have firmly left the death-edge regime.

**Per-metric: refined → push_edges**

| metric | refn max | push max | ratio | refn top-5 mean | push top-5 mean | ratio |
|--------|---------:|---------:|------:|----------------:|----------------:|------:|
| correlation_length_mean        |  82.7  |   75.8  | **0.92** |  75.2 |  74.4 | **0.99** |
| largest_patch_temporal_std     | 15824  |  17277  | 1.09 | 11358 | 15747 | 1.39 |
| n_distinct_genomes_temporal_std| 12984  |  13602  | 1.05 |  9519 | 12636 | 1.33 |
| unique_top_genomes             |    51  |    194  | **3.80** |    51 |   192 | 3.77 |
| F_std_mean                     |   0.5  |    0.4  | 0.84 |   0.4 |   0.4 | 0.85 |

Top-genomes 3.8× rise is the un-saturation effect (sample bins went
from 51 → 200), not a real regime change. The honest comparison is
the *fraction* of bins that saw a different dominant: refined hit
51/51 = 100 % at sample_every=100; push_edges hit 192/200 ≈ 96 % —
turnover is *qualitatively the same*, both runs effectively never
have the same top genome at consecutive samples.

**The genuine signals:**
- `correlation_length_mean` has **plateaued** (slight decline). The
  refined-scan max of 82.7 looks like a real ceiling at N=256 — pushing
  gdiff past 0.08 didn't help, and may slightly hurt.
- `largest_patch_temporal_std` and `n_distinct_genomes_temporal_std`
  top-5 means each gained ~30–40 %, so temporal turbulence is still
  improving.
- Composite max 0.853 → 0.898, top-5 mean 0.714 → 0.850 — modest gains
  because correlation length stalled but temporal metrics rose.

## Top 5 push_edges configs

```
 score  f_inc    tax   gdiff  mu_lut  mu_eg  m_sc  r_mu | corr_L  patch_std  div_std top_g  alive
 0.898  0.013  0.035  0.080   0.003  0.001   1.2  False|   73.9      15484  10706    188  0.322
 0.884  0.018  0.040  0.060   0.003  0.010   1.2  False|   74.6      15242  10467    182  0.360
 0.849  0.013  0.030  0.060   0.010  0.001   1.0  False|   69.8      13781  12525    158  0.403
 0.818  0.008  0.030  0.060   0.030  0.001   1.2  True |   61.7      12466  13602    158  0.217
 0.801  0.010  0.035  0.060   0.030  0.001   1.2  True |   66.5      11215  12290    168  0.231
```

## Distribution of top-25 in parameter space

```
   gdiff: 0.06: 16  0.08:  9             ← only 0.06 & 0.08 represented;
                                            none of {0.12, 0.15, 0.20}
                                            survived to top-25
 m_scale: 0.7:  2  0.8:  6  1.0: 6  1.2: 11   ← skewed strongly toward 1.2
food_inc: 0.008: 4  0.01: 8  0.013: 10  0.018: 3
     tax: 0.03: 11  0.035: 10  0.04: 4
  mu_lut: 0.001: 2  0.003: 9  0.01: 5  0.03: 9
```

**Clear signals from the distribution:**
- **gdiff is bracketed**: the optimum is at 0.06–0.08. Pushing higher
  was wrong-direction — values above 0.08 are *all* suboptimal. The
  refined-scan boundary at 0.08 was right at the peak.
- **m_scale wants to go higher still**: 11 of 25 top configs are at the
  new upper limit 1.2. The optimum may be at 1.5–2.0. *Not bracketed.*
- food_inc / tax / mu_lut distributions are spread fairly broadly across
  their values — no single sharp optimum, soft optima around food_inc=
  0.013, tax=0.030, mu_lut=0.003 or 0.030 (bimodal — two regimes).

## Diminishing returns?  Mixed.

- **Spatial scale: yes, plateaued.** correlation_length_mean has hit a
  ceiling at ~75–80 cells (≈ 30 % of N=256 grid). Probably grid-size
  bound — bigger grids may show bigger structures.
- **Temporal metrics: still gaining** modestly — largest_patch_temporal_std
  and div_temporal_std top-5 means rose 30–40 %.
- **Yield: hugely better** (98 % vs 50 %), so the population of
  *interesting* runs is denser.

## Next steps

- [x] Run push_edges scan
- [x] Compare to refined
- [ ] **Validate top configs at N=512 in SDL.** Visual confirmation
      they actually do "big patches + ongoing genomic turnover"
      is the only way to know if the metrics are pointing where we
      want.
- [ ] (Optional) **Scan 4: push m_scale further.** m_scale ∈
      {1.0, 1.2, 1.5, 1.8, 2.0}, gdiff fixed at {0.05, 0.06, 0.08}
      to bracket the spatial peak, food_inc / tax / mu_lut at the
      best-distribution values. ~150 combos, ~2 min.
- [ ] (Optional) **Larger N** (384 or 512) to break the
      correlation_length_mean ceiling and see whether the dynamics
      retain their character at bigger scale.
