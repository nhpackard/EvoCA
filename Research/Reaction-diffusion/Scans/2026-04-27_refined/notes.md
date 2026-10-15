# Scan 2026-04-27 — refined search around the productive region

## Question

The 2026-04-27_initial scan found 25/150 candidates surviving the
extinct/saturated filters, concentrated around `food_inc=0.015,
tax=0.035, gdiff∈{0.01, 0.05}, m_scale=0.6`. Two qualitative branches:
- **gdiff=0.05** branch (cfg49, cfg34, cfg30): corr_L≈21–49, more
  saturated alive (~0.4–0.7), strong genomic turnover
- **gdiff=0.01** branch (cfg70): corr_L≈33, lower density (~0.5),
  large patches with strong temporal variation

This refinement pushes the grid:
- **denser around the productive food_inc/tax/gdiff range**
- **wider** in gdiff (testing 0.005–0.08) to find boundaries of both branches
- **finer** in mu_lut and m_scale (3 values instead of 2)
- adds **food_inc=0.012** to probe the death edge

## Search axes

```
food_inc      ∈ {0.010, 0.013, 0.015, 0.018, 0.022}    (5 values)
tax           ∈ {0.025, 0.030, 0.035, 0.040, 0.045}    (5)
gdiff         ∈ {0.005, 0.010, 0.015, 0.020, 0.030, 0.050, 0.080}  (7)
mu_lut        ∈ {0.001, 0.003, 0.010, 0.030}           (4)
mu_egenome    ∈ {0.001, 0.010}                         (2)
m_scale       ∈ {0.4, 0.6, 0.8}                        (3)
restricted_mu ∈ {True, False}                          (2)
```

Full grid = 16 800 combos. Random sample of 200.

Fixed: N=256, n_steps=5000, sample_every=100, init='halfplane',
seed=0, shadow=True. Same as initial scan for direct comparability.

## Diminishing-returns question

Compare against the initial scan: are the new top configs *meaningfully
better* than the previous top, or have we hit a plateau? Specifically:

- For each metric M (correlation_length_mean, largest_patch_temporal_std,
  n_distinct_genomes_temporal_std, unique_top_genomes), what is
  `max(refined[M]) / max(initial[M])`? A ratio close to 1.0 means
  diminishing returns; substantially > 1.0 means the refinement
  uncovered new territory.
- Where do the new top configs sit in parameter space — are they
  *different* from the initial winners, or just nearby variants?
- Distribution of composite scores: is the median candidate score
  rising, or is only the top moving?

## Results — 103 s on 12 cores

- 200 configs run, 99 candidates after extinct/saturated filtering
- Compare to initial scan: 25/150 candidates → 99/200 (17 % → 50 % yield)

**Per-metric improvement (refined max / initial max):**

| metric | init max | refn max | ratio | init top-5 mean | refn top-5 mean | ratio |
|--------|---------:|---------:|------:|----------------:|----------------:|------:|
| correlation_length_mean        |  50.6 |    82.7 | **1.63** |   40.7 |  75.2 | **1.85** |
| largest_patch_temporal_std     |  8478 |   15825 | **1.87** |   8131 | 11358 | **1.40** |
| n_distinct_genomes_temporal_std|  8085 |   12984 | **1.61** |   7753 |  9519 | 1.23 |
| unique_top_genomes             |    51 |      51 | 1.00     |     50 |    51 | 1.02 |
| F_std_mean                     |   0.5 |     0.4 | 0.98     |    0.4 |   0.4 | 0.97 |

`unique_top_genomes` is already at the sample-count ceiling (51 of 51
samples = every sample saw a different dominant genome) — saturated as a
discriminator at this `sample_every`.

**Composite score (default 4-metric ranker):** initial max 0.750 → refined
max **0.853** (1.14×). Modest because the top-genomes axis is saturated;
the *individual* spatial metrics improved much more.

## Top 5 refined configs (composite-ranked)

```
 score  f_inc    tax   gdiff  mu_lut  mu_eg  m_sc  r_mu | corr_L  patch_std  div_std top_g  alive
 0.853  0.010  0.035  0.080   0.010  0.001   0.8  True |   50.8      13612  12984.4    51  0.254
 0.748  0.013  0.035  0.080   0.001  0.010   0.8  False|   59.3      15825   4110.3    50  0.334
 0.688  0.022  0.040  0.050   0.030  0.010   0.8  True |   67.6       8291   7701.3    50  0.490
 0.644  0.018  0.035  0.050   0.010  0.010   0.6  True |   59.2       8164   6891.4    50  0.456
 0.635  0.010  0.030  0.050   0.030  0.010   0.8  True |   59.7       5938   8469.8    51  0.273
```

Notice:
- **gdiff = 0.08 wins twice** at the top — the highest value in our new
  range, suggesting the optimum is at or above 0.08 (not yet bracketed).
- **m_scale = 0.8 dominates 4/5** of the top — also the highest value in
  the new range; not yet bracketed.
- Top configs sit at low `alive_density` (0.25–0.49) — sparser populations
  than the initial winners (0.37–0.70). System is operating closer to the
  death edge with bigger spatial gradients.

## Diminishing returns?  No.

Parameter-space novelty: each top-5 refined config's normalised L2 distance
to its nearest initial-top-5 neighbour:

```
refn rank  refn cfg  min L2 dist
        1       149       1.187
        2       198       1.071
        3        74       0.713
        4        71       0.718
        5       161       0.633
```

The top refined config sits 1.19 normalised units from any initial winner
(normalised across the union grid; ~1.2 is "different parameter regime"
not "nearby variant"). The refinement has discovered **new territory** —
in particular the `gdiff=0.08, m_scale=0.8` corner — where spatial scale
and patch turbulence are both substantially larger than anything the
initial scan reached.

**We are not at a plateau yet.** Both `gdiff` and `m_scale` winners are
at the *upper edge* of our refined ranges, suggesting another iteration
that pushes those further (gdiff ∈ [0.06, 0.20], m_scale ∈ [0.7, 1.2])
should keep gaining.

## Next steps

- [x] Run refined scan
- [x] Diminishing-returns analysis vs initial scan
- [ ] Visualise top refined configs in SDL at N=512 to confirm the
      "bigger spatial structure with active evolution" story matches
      the metrics
- [ ] Push further if visuals look promising: a third scan with
      `gdiff ∈ {0.06, 0.08, 0.12, 0.15, 0.20}` and
      `m_scale ∈ {0.7, 0.8, 1.0, 1.2}`, food_inc / tax fixed near
      `(0.010–0.018, 0.030–0.040)`, mu_lut spread.
- [ ] Increase `sample_every` from 100 → 25 (4× more samples per run)
      so `unique_top_genomes` un-saturates and becomes a useful
      discriminator again.
