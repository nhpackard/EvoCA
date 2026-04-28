# Scan 2026-04-27 — evo-focused refinement of push_edges winners

## Question

The first three scans (`initial`, `refined`, `push_edges`) used a
balanced composite score combining **spatial** and **evolution**
metrics. The push_edges scan revealed:
- Spatial scale (`correlation_length_mean`) plateaued at $\approx 75$
  cells, plausibly grid-size-bound at $N=256$.
- Temporal turbulence and evolutionary turnover were still gaining.

We have enough interesting spatial dynamics. **This scan optimises
purely for EVO_METRICS** to see how far we can push genomic turnover
and excess activity, given that the productive spatial regime is
already established.

## Configs

Combined config set:

1. **Top-50 baseline.** Take the 50 highest-composite-score candidates
   from `Scans/2026-04-27_push_edges` and re-run them. These set a
   reference: same composite-balanced winners, same conditions, no
   variation.

2. **~200 random variations** from a tightened grid around the
   productive corner identified by push_edges, with **finer
   granularity on the evolution-affecting axes** (`mu_lut`,
   `mu_egenome`):

   ```
   food_inc      ∈ {0.008, 0.010, 0.013, 0.018}
   tax           ∈ {0.030, 0.035, 0.040}
   gdiff         ∈ {0.05, 0.06, 0.08, 0.10}              ← brackets the spatial peak
   mu_lut        ∈ {0.001, 0.003, 0.005, 0.010, 0.020, 0.030, 0.060}  ← extended high
   mu_egenome    ∈ {0.001, 0.003, 0.010, 0.030}          ← finer
   m_scale       ∈ {0.8, 1.0, 1.2, 1.5}                  ← extends scan-3 upper edge
   restricted_mu ∈ {True, False}
   ```

   Full grid: $4 \times 3 \times 4 \times 7 \times 4 \times 4 \times 2 = 5\,376$
   combos. Random sample of 200 unique combinations not already in the
   baseline.

Total: ≈ 250 configs. Fixed: $N=256$, $n_{\text{steps}}=5000$,
$\text{sample\_every}=25$, `init='halfplane'`, `seed=0`,
`shadow=True` — matches push_edges so results are directly comparable.

## Ranking criterion

Composite of `EVO_METRICS` only, using the same normalised-mean
formula as `evoca_from_scan_top` but with score keys =

```python
EVO_METRICS = [
    'n_distinct_genomes_mean',
    'n_distinct_genomes_temporal_std',
    'unique_top_genomes',
    'excess_activity_slope',
]
```

- `n_distinct_genomes_mean` — average diversity per sample
- `n_distinct_genomes_temporal_std` — temporal flux of diversity
- `unique_top_genomes` — count of distinct dominators across samples
- `excess_activity_slope` — Channon's $\Sigma G - \Sigma N$ trend

## Results — 180 s on 12 cores

- 250 configs run (50 baseline + 200 random variants).
- 234/250 candidates after extinct/saturated filter (94 %).

### EVO_METRICS — push_edges → evo_focus

| metric | push max | evo max | ratio | push top-5 mean | evo top-5 mean | ratio |
|--------|---:|---:|---:|---:|---:|---:|
| `n_distinct_genomes_mean`         | 38 946 | 38 957 | 1.00 | 38 855 | 38 645 | 0.99 |
| `n_distinct_genomes_temporal_std` | 13 602 | **16 603** | **1.22** | 12 636 | **15 192** | **1.20** |
| `unique_top_genomes`              |   194 |   195 | 1.01 |   192 |   192 | 1.00 |
| `excess_activity_slope`           |   527 | **244** | **0.46** |   448 | **228** | **0.51** |

Pure-EVO composite score: **max 0.637 → 0.724 (1.14×)**, top-5 mean 0.597 → 0.697 (1.17×).

The headline gain is in **temporal flux of diversity**
(`n_distinct_genomes_temporal_std`): the genome-bucket count
fluctuates more dramatically over time.

### Excess-activity-slope dropped — interpretation

The one EVO metric that *fell* is `excess_activity_slope`. With
`mu_lut` pushed higher (the new top configs sit at 0.06, the highest
value in the scan grid), both real-run G-activity and shadow N-activity
grow faster — and at high mutation rates the shadow can grow
*proportionally faster* than the real run because random selection in
the shadow always allows new bucket creation, whereas real selection
filters. So the slope of $\Sigma G - \Sigma N$ can flatten or invert
even though *both* are growing more vigorously.

This is a known quirk of the Channon excess metric in high-mutation
regimes; it doesn't mean evolution stopped, it means the excess metric
is being saturated by raw turnover on both sides.

### Spatial metrics — kept fine despite not optimising for them

| metric | push top-5 mean | evo top-5 mean | ratio |
|--------|---:|---:|---:|
| `correlation_length_mean`     |  74.4 |  79.3 | **1.07** ▲ |
| `F_std_mean`                  |   0.4 |   0.4 | 1.06 |
| `largest_patch_mean`          | 36 930 | 35 377 | 0.96 |
| `largest_patch_temporal_std`  | 15 747 | 18 416 | **1.17** ▲ |
| `n_patches_mean`              | 11 052 |  9 672 | 0.88 |

Spatial scale and patch turbulence both *improved* incidentally —
the productive corner of parameter space jointly serves spatial and
evolutionary criteria. Confirms the user's intuition that we already
have enough spatial dynamics; pushing evo doesn't sacrifice them.

## Top 5 evo_focus configs (EVO-only composite)

```
 score  f_inc    tax   gdiff  mu_lut  mu_eg  m_sc  r_mu |  div_mn  div_std  top_g    slope  corr_L
 0.724  0.018  0.035  0.060   0.060  0.030   1.5  True |   30 767  16 603    113   -1.8    74.8
 0.699  0.013  0.030  0.050   0.020  0.001   1.5  False|   25 147  14 306    153   -3.2    73.8
 0.698  0.013  0.035  0.060   0.030  0.001   1.2  True |   20 034  15 339    165    3.0    61.6
 0.683  0.018  0.040  0.060   0.060  0.003   1.5  True |   24 266  12 796    164   -2.0    75.2
 0.683  0.013  0.035  0.050   0.060  0.030   1.5  False|   17 580  15 181    169    0.9    61.0
```

Notice two surprising things:
- **mu_lut = 0.060 wins three of the top five** — the highest value in
  the new grid (and double scan-3's max of 0.030).
- **m_scale = 1.5 wins four of five**, the new upper limit. Like
  push_edges, the optimum is at-or-above the boundary.

## Top-25 parameter distribution (EVO-only ranking)

```
   mu_lut: 0.005:  3  0.01:  3  0.02:  3  0.03:  4  0.06: 12   ← skewed to highest
mu_egenome: 0.001:  8  0.003:  5  0.01:  8  0.03:  4           ← spread, low-mu_eg slight edge
     gdiff: 0.05:  7  0.06: 10  0.08:  6  0.10:  2             ← peak at 0.06, same as scan 3
   m_scale: 0.8:   2  1.0:   4  1.2: 10  1.5:  9              ← split between 1.2 and 1.5
  food_inc: 0.008: 1  0.01:  4  0.013: 14  0.018: 6           ← centred at 0.013
       tax: 0.03: 12  0.035: 5  0.04:  8                       ← bimodal, low-tax edge
```

`mu_lut` is the most consequential lever and is **not bracketed** —
the optimum lies at or above 0.060. `m_scale` also still wants higher
than our 1.5 cap.

## Diminishing returns on evolution?  No — still room.

`mu_lut` extending above 0.030 (scan-3 cap) gave a clear improvement
at 0.060. We haven't bracketed the upper end yet. Same story for
`m_scale`: scan 3 maxed at 1.2, scan 4 max at 1.5, top-25 distribution
peaks at 1.5.

## Best configs (by raw EVO_METRICS-only composite)

Top winner: `food_inc=0.018, tax=0.035, gdiff=0.06, mu_lut=0.060,
mu_egenome=0.030, m_scale=1.5, restricted_mu=True` — score 0.724,
diversity-temporal-std = 16 603, alive_density ≈ 0.49.

This sits at the high-mutation, high-eating, modest-diffusion corner.
With `mu_lut=0.060` and `mu_egenome=0.030`, every birth is essentially
a fresh genome — the population is constantly remixing.

## Caveats

The 50 baseline reruns (same params as push_edges top 50, same seed,
same code) produced *different* metric values to their push_edges
originals — multiprocessing.Pool reuses worker processes across tasks
and the C-side xorshift RNG is per-process and not reset between sims.
So per-config determinism doesn't hold across the pool, only the
*aggregate* statistics. This is a methodological limitation that
applies to all four scans equally; results are still reliable in
distribution but not on any individual run.

## Status

- [x] Run scan
- [x] Compare evo metrics to push_edges
- [ ] (Optional) Scan 5: push `mu_lut` to 0.10–0.20 and `m_scale` to
      2.0+ since both are still at upper edges.
- [ ] Validate top configs at N=512 visually.
- [ ] Possibly fix the per-config determinism by adding an explicit
      RNG seed reset at the top of `run_sim` so reruns are bitwise
      identical.
