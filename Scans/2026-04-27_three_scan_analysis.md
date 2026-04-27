# Three-Scan Analysis: Searching for Spatial Structure + Open-Ended Evolution

**Date**: 2026-04-27
**Scans**: `2026-04-27_initial/`, `2026-04-27_refined/`, `2026-04-27_push_edges/`

## Question

We are searching for a parameter regime in EvoCA where the simulation
exhibits **both**:

1. **Large-scale spatial structure** — patches of alive cells that span
   tens of grid cells, persist long enough to be visible, and continue
   changing rather than locking into a static "garden plot" attractor.
2. **Ongoing genomic turnover** — the dominant genome keeps being
   displaced by mutants, rather than selection finding one stable
   attractor and freezing there.

The earlier interactive exploration showed that low-mutation regimes
quickly lock in: selection finds a static spatial attractor and
dynamics stop. We want to know whether there is a parameter pocket
where the two timescales (spatial reaction-diffusion and evolutionary
drift) are well-matched, so that one keeps disrupting the other.

## Methodology

Each scan: random sample from a discrete parameter grid, run headless
via `multiprocessing.Pool` (12 workers), 5000 ticks at $N = 256$, with
the neutral (Channon) shadow active for excess-activity calibration,
left-halfplane initialisation, fixed seed 0.

| scan dir | # configs | grid combos | `sample_every` | runtime |
|----------|-----:|---------:|---------:|--------:|
| `2026-04-27_initial` | 150 | 1 152 | 100 | 132 s |
| `2026-04-27_refined` | 200 | 16 800 | 100 | 103 s |
| `2026-04-27_push_edges` | 250 | 7 680 | 25 | 179 s |

Fixed across all three scans: $N = 256$, $n_{\text{steps}} = 5000$,
`init='halfplane'` (left half alive), `seed=0`, neutral shadow on,
`lut='gol'`, `egenome='uniform'` (value $0\text{b}000011$ = plus-shape
fiducial), `v_density=0.5`, `f_init=0.5`, `F_init=1.0`.

## Metrics

### Filtering

A scan row is a *candidate* if the run neither died nor saturated:

$$\text{candidate}(r) \iff \neg\,r.\text{extinct} \;\wedge\; \overline{\text{alive\_density}}_r \le 0.95$$

i.e. the alive count never reached zero AND the time-mean alive
fraction stayed below 95 %. The remainder is the "interesting middle":
runs with non-trivial spatial structure that didn't go extinct.

### Spatial-structure metrics

**`correlation_length_mean`** — first-zero-crossing radius of the 2-D
autocovariance of the alive mask, averaged over angle then over the
sample times.

For each snapshot $A(\mathbf{x}, t)$ of the binary alive mask, define
the normalised autocovariance:

$$\tilde{A}(\mathbf{r}, t) = \frac{\mathcal{F}^{-1}\!\left[\;|\mathcal{F}[A - \bar{A}]|^2\;\right](\mathbf{r})}{\mathcal{F}^{-1}\!\left[\;|\mathcal{F}[A - \bar{A}]|^2\;\right](\mathbf{0})}$$

For each radius $r$, average over 16 evenly-spaced angles
$\theta_k = 2\pi k/16$:

$$C(r, t) = \frac{1}{16}\sum_{k=0}^{15} \tilde{A}\!\big(\lfloor r\cos\theta_k\rceil,\ \lfloor r\sin\theta_k\rceil,\ t\big)$$

The correlation length at time $t$ is

$$L(t) \;=\; \min\big\{r \in \mathbb{Z}_{>0} \;:\; C(r, t) \le 0\big\}$$

and the metric reports $\langle L(t) \rangle_{t \in T_{\text{sample}}}$.

**`largest_patch_mean`** — time-averaged size of the largest connected
alive component (4-neighbour connectivity, computed with
`scipy.ndimage.label`).

**`largest_patch_temporal_std`** — temporal standard deviation of the
same. High value indicates patches that grow and shrink dramatically
over the run.

**`F_std_mean`** — time average of $\sigma_{\mathbf{x}}[F(\mathbf{x}, t)]$,
the spatial standard deviation of the food field. High value indicates
spatial heterogeneity in food.

**`n_patches_mean`** — time-averaged count of distinct alive
components.

### Evolution metrics

**`n_distinct_genomes_mean`** — time average of $|\{h : \text{pop}_h(t) > 0\}|$,
where the index $h$ runs over LUT-content hashes (FNV-1a over the 32-byte
LUT). One bucket per distinct LUT.

**`n_distinct_genomes_temporal_std`** — temporal standard deviation of
the same. High value indicates the bucket count fluctuates over time,
i.e. ongoing genomic flux.

**`unique_top_genomes`** — count of distinct genomes ever holding the
#1-population spot across sample times:

$$U \;=\; \left|\bigcup_{t \in T_{\text{sample}}} \big\{\arg\max_h \text{pop}_h(t)\big\}\right|$$

High value = dominance changes hands often. Capped at $|T_{\text{sample}}| = n_{\text{steps}}/\text{sample\_every}$
(50 in scans 1 and 2, 200 in scan 3). Saturation at the cap means
"every sample saw a different dominator" and the metric stops
discriminating among the strongest configs.

**`excess_activity_slope`** — slope of the linear regression of
Channon's excess-activity time series $\Sigma G(t) - \Sigma N(t)$ vs $t$,
where $G$ and $N$ are real-run and shadow-run total activity over all
buckets. Positive slope = real run is pulling away from the neutral
shadow over time.

### Composite score (used by `evoca_from_scan_top` and ranking tables)

For each candidate $i$ and each score-key metric $m \in K$, normalise
to $[0, 1]$ across the candidate set:

$$\hat{m}_i \;=\; \frac{m_i - \min_{j \in \mathcal{C}} m_j}{\max_{j \in \mathcal{C}} m_j - \min_{j \in \mathcal{C}} m_j}$$

then take the unweighted mean:

$$\text{score}_i \;=\; \frac{1}{|K|} \sum_{m \in K} \hat{m}_i \;\;\in\; [0, 1].$$

Default $K = \{$ `correlation_length_mean`, `largest_patch_temporal_std`, `n_distinct_genomes_temporal_std`, `unique_top_genomes` $\}$ — one
spatial-scale, one spatial-temporal, one evolution-temporal, one
evolution-cumulative.

## Scan Coverage

| axis | initial | refined | push_edges |
|------|---------|---------|------------|
| `food_inc`      | $\{0.015,\,0.025,\,0.040\}$ | $\{0.010,\,0.013,\,0.015,\,0.018,\,0.022\}$ | $\{0.008,\,0.010,\,0.013,\,0.018\}$ |
| `tax`           | $\{0.020,\,0.025,\,0.035\}$ | $\{0.025,\,0.030,\,0.035,\,0.040,\,0.045\}$ | $\{0.030,\,0.035,\,0.040\}$ |
| `gdiff`         | $\{0.005,\,0.01,\,0.02,\,0.05\}$ | $\{0.005,\,0.01,\,0.015,\,0.02,\,0.03,\,0.05,\,0.08\}$ | $\{0.06,\,0.08,\,0.12,\,0.15,\,0.20\}$ |
| `mu_lut`        | $\{0.001,\,0.003,\,0.01,\,0.03\}$ | $\{0.001,\,0.003,\,0.01,\,0.03\}$ | $\{0.001,\,0.003,\,0.01,\,0.03\}$ |
| `mu_egenome`    | $\{0.001,\,0.01\}$ | $\{0.001,\,0.01\}$ | $\{0.001,\,0.01\}$ |
| `m_scale`       | $\{0.4,\,0.6\}$ | $\{0.4,\,0.6,\,0.8\}$ | $\{0.7,\,0.8,\,1.0,\,1.2\}$ |
| `restricted_mu` | $\{T,\,F\}$ | $\{T,\,F\}$ | $\{T,\,F\}$ |

The grid was tightened around the productive region between scan 1
and scan 2, then translated upward (in `gdiff` and `m_scale`) for
scan 3.

## Yield (candidate fraction)

| scan | extinct | saturated | candidates | yield |
|------|---:|---:|---:|---:|
| initial    | 31 | 94 | 25 | 17 % |
| refined    | ~25 | ~76 | 99 | 50 % |
| push_edges | ~3 | ~2 | 245 | 98 % |

The "interesting middle" went from 17 % of runs to almost everything.
We have moved firmly out of the death-edge / saturation regimes the
initial scan stumbled around in.

## Per-Metric Maxima Across Scans

| metric | initial max | refined max | push_edges max | overall ratio |
|--------|---:|---:|---:|---:|
| `correlation_length_mean`        |   50.6 |   82.7 |   75.8 | 1.50× |
| `largest_patch_temporal_std`     |   8478 |  15825 |  17277 | 2.04× |
| `n_distinct_genomes_temporal_std`|   8085 |  12984 |  13602 | 1.68× |
| `unique_top_genomes`             |     51 |     51 |    194 | 3.80× ($\dagger$) |
| `F_std_mean`                     |   0.49 |   0.43 |   0.40 | 0.82× |
| **composite max**                | 0.750 | 0.853 | 0.898 | 1.20× |

$\dagger$ The 3.8× rise in `unique_top_genomes` is mostly the
sample-cap rising from 51 to 200 (sample_every dropped from 100 to 25).
The fraction of sample bins that saw a unique dominator is essentially
the same in scans 2 and 3 ($\approx 0.96$–$1.00$): turnover is
qualitatively saturated in both cases — every sample sees a different
top genome.

## Top-5 Composite Configs Per Scan

### Initial (1 152-combo grid, 25 candidates)

| score | food_inc | tax   | gdiff | mu_lut | mu_eg | m_sc | r_mu  | corr_L | patch_std | div_std | top_g | alive |
|-----:|---------:|------:|------:|-------:|------:|----:|:-----:|------:|----------:|--------:|-----:|------:|
| 0.750 |   0.025 | 0.035 | 0.050 |  0.010 | 0.001 | 0.6 | T     |  21.1 |     8 314 |   7 546 |   43 | 0.696 |
| 0.721 |   0.025 | 0.035 | 0.050 |  0.010 | 0.010 | 0.6 | F     |  21.2 |     7 689 |   7 647 |   43 | 0.698 |
| 0.681 |   0.015 | 0.020 | 0.010 |  0.010 | 0.010 | 0.6 | F     |  32.6 |     7 717 |   7 850 |   19 | 0.511 |
| 0.630 |   0.015 | 0.035 | 0.050 |  0.030 | 0.010 | 0.6 | T     |  49.2 |     4 462 |   5 226 |   50 | 0.373 |
| 0.611 |   0.015 | 0.020 | 0.010 |  0.003 | 0.001 | 0.6 | F     |  32.4 |     7 900 |   5 846 |   23 | 0.518 |

### Refined (16 800-combo grid, 99 candidates)

| score | food_inc | tax   | gdiff | mu_lut | mu_eg | m_sc | r_mu  | corr_L | patch_std | div_std | top_g | alive |
|-----:|---------:|------:|------:|-------:|------:|----:|:-----:|------:|----------:|--------:|-----:|------:|
| 0.853 |   0.010 | 0.035 | 0.080 |  0.010 | 0.001 | 0.8 | T     |  50.8 |    13 612 |  12 984 |   51 | 0.254 |
| 0.748 |   0.013 | 0.035 | 0.080 |  0.001 | 0.010 | 0.8 | F     |  59.3 |    15 825 |   4 110 |   50 | 0.334 |
| 0.688 |   0.022 | 0.040 | 0.050 |  0.030 | 0.010 | 0.8 | T     |  67.6 |     8 291 |   7 701 |   50 | 0.490 |
| 0.644 |   0.018 | 0.035 | 0.050 |  0.010 | 0.010 | 0.6 | T     |  59.2 |     8 164 |   6 891 |   50 | 0.456 |
| 0.635 |   0.010 | 0.030 | 0.050 |  0.030 | 0.010 | 0.8 | T     |  59.7 |     5 938 |   8 470 |   51 | 0.273 |

### Push_edges (7 680-combo grid, 245 candidates)

| score | food_inc | tax   | gdiff | mu_lut | mu_eg | m_sc | r_mu  | corr_L | patch_std | div_std | top_g | alive |
|-----:|---------:|------:|------:|-------:|------:|----:|:-----:|------:|----------:|--------:|-----:|------:|
| 0.898 |   0.013 | 0.035 | 0.080 |  0.003 | 0.001 | 1.2 | F     |  73.9 |    15 484 |  10 706 |  188 | 0.322 |
| 0.884 |   0.018 | 0.040 | 0.060 |  0.003 | 0.010 | 1.2 | F     |  74.6 |    15 242 |  10 467 |  182 | 0.360 |
| 0.849 |   0.013 | 0.030 | 0.060 |  0.010 | 0.001 | 1.0 | F     |  69.8 |    13 781 |  12 525 |  158 | 0.403 |
| 0.818 |   0.008 | 0.030 | 0.060 |  0.030 | 0.001 | 1.2 | T     |  61.7 |    12 466 |  13 602 |  158 | 0.217 |
| 0.801 |   0.010 | 0.035 | 0.060 |  0.030 | 0.001 | 1.2 | T     |  66.5 |    11 215 |  12 290 |  168 | 0.231 |

## Top-25 Parameter Distribution (Scan 3)

The distribution among the 25 top push_edges configs reveals which
axes are bracketed and which still want to push:

```
   gdiff: 0.06: 16   0.08:  9                         (none of {0.12, 0.15, 0.20})
 m_scale: 0.7:  2   0.8:  6   1.0:  6   1.2: 11      (skewed toward 1.2)
food_inc: 0.008: 4  0.010: 8  0.013: 10  0.018: 3    (centred at 0.013)
     tax: 0.030: 11 0.035: 10 0.040:  4              (low-tax preferred)
  mu_lut: 0.001: 2  0.003: 9  0.010:  5  0.030: 9    (bimodal: low or high)
```

`gdiff` and `tax` are bracketed. `m_scale` is *not* bracketed — its
distribution piles up at the upper search limit, suggesting the optimum
is at or above $m_{\text{scale}} = 1.2$. `mu_lut` is bimodal: two
distinct evolutionary regimes (low-rate "preserve a working LUT" vs
high-rate "constantly remix") both produce competitive configs.

## Diminishing Returns?

Per metric, comparing scan 3 to scan 2:

| metric | refn → push (top-5 means) | trend |
|--------|---:|---|
| `correlation_length_mean`         | 75.2 → 74.4 | **plateau** |
| `largest_patch_temporal_std`      | 11 358 → 15 747 | still gaining (1.39×) |
| `n_distinct_genomes_temporal_std` | 9 519 → 12 636 | still gaining (1.33×) |
| `unique_top_genomes`              | 51 / 51 → 192 / 200 | metric saturates each time |
| candidate yield                   | 50 % → 98 % | huge gain |
| composite score (top-5)           | 0.714 → 0.850 | gaining |

**Spatial scale has plateaued.** `correlation_length_mean` topped out at
$\approx 75$–$83$ across both refined and push_edges scans, which is
$\approx 30\,\%$ of the $N = 256$ grid. This looks like a grid-size
ceiling, not a parameter ceiling: with periodic boundaries, an
autocovariance can decorrelate within a half-period. The way to test
that hypothesis is to scan at $N = 384$ or $512$.

**Temporal turbulence is still gaining.** Both `largest_patch_temporal_std`
and `n_distinct_genomes_temporal_std` continue to improve. Patches
keep growing and shrinking more dramatically; the genome population
keeps fluctuating more.

**Survival yield is essentially solved.** From 17 % at the initial
scan to 98 % at push_edges — the $(\text{food\_inc},\,\text{tax},\,\text{gdiff},\,m_{\text{scale}})$
joint distribution we now sample in is overwhelmingly viable.

## Parameter-Space Novelty Across Scans

Normalised L2 distance from each top-5 config to its nearest top-5
neighbour in the *previous* scan (axes min-max scaled to $[0, 1]$
over the union of all rows):

| comparison | rank-1 dist | rank-2 dist | rank-3 dist | rank-4 dist | rank-5 dist |
|-----------|---:|---:|---:|---:|---:|
| refined vs initial    | 1.187 | 1.071 | 0.713 | 0.718 | 0.633 |
| push_edges vs refined | 1.120 | 0.675 | 1.087 | 0.905 | 0.858 |

A top config moving by $\sim 1$ normalised unit means it sits in a
different parameter region than any winner of the previous scan. Each
iteration genuinely opened new territory rather than just refining
existing winners. By scan 3, however, the dominant `gdiff` value
shifted only from 0.05 to 0.06–0.08 (one grid step), so the rate of
parameter migration is slowing.

## Observations

1. **The "interesting middle" regime exists and is large.** With 98 %
   yield at scan 3, almost any reasonable parameter combination in the
   $(\text{food\_inc} \in [0.008, 0.018],\,\text{tax} \in [0.030, 0.040],\,\text{gdiff} \in [0.06, 0.08],\,m_{\text{scale}} \in [0.7, 1.2])$
   box gives a population that survives but does not saturate. The
   regime is robust, not fragile.

2. **Correlation length is grid-size-bound.** The plateau at $\approx 75$
   cells $\approx N/3$ suggests we have hit the natural limit imposed
   by the autocovariance reaching a half-period before zero crossing.
   Scaling $N$ should test this directly.

3. **`m_scale` is genuinely the surprise lever.** Initial scans capped
   it at 0.6; refined extended to 0.8 and the optimum moved up; push
   extended to 1.2 and the distribution still piles at the top. Eating
   aggressiveness matters more than the initial ranges suggested.

4. **`gdiff` peaks at 0.06–0.08.** Fully bracketed by scan 3. Higher
   diffusion smears patches; lower diffusion fragments the population
   into death.

5. **`mu_lut` is bimodal.** Both 0.003 and 0.030 produce competitive
   configs. The 0.003 regime preserves coherent CA dynamics with
   occasional invasions; the 0.030 regime constantly remixes. Two
   different routes to the same composite outcome.

6. **`unique_top_genomes` is informative only as a fraction of
   sample bins.** All three scans hit 96–100 % of the sample-cap — the
   raw count is just measuring the cap. Future analyses should
   probably compute and report
   $U / |T_{\text{sample}}|$ instead.

## Recommended Next Steps

1. **Visual validation at $N = 512$.** Pull the top configs from scan 3
   into the SDL viewer. Confirm the metrics ("big patches with ongoing
   genomic turnover") match what one actually sees. Easy first move.

2. **Larger-grid scan ($N \in \{384, 512\}$).** Test whether the
   correlation-length plateau is a grid-size artefact. Smaller config
   set, fixed `food_inc=0.013, tax=0.035, gdiff∈{0.06, 0.08}, mu_lut∈{0.003, 0.030}`,
   `m_scale ∈ {0.8, 1.0, 1.2}`. ~24 combos, would take 10–30 min depending
   on $N$ (compute scales as $N^2$).

3. **`m_scale` extension scan.** Bracket the upper bound. `m_scale ∈
   {1.0, 1.2, 1.5, 1.8, 2.0}`, gdiff fixed at 0.06–0.08, others at
   the top of the scan-3 distribution. ~150 combos, ~2 min.

4. **Multi-seed runs** of the top configs to confirm dynamics are
   robust to initial conditions, not artefacts of `seed=0`.

## Reproducing this analysis

```python
# Assumes you're in the EvoCA repo root.
import sys, csv
sys.path.insert(0, 'python')

from evoca_explore import (run_sim, save_scan_config,
                           evoca_from_scan_top,
                           nearest_params, nearest_evo, nearest_spatial)

# Re-run any scan:
#   cd Scans/2026-04-27_<name> && python3 scan.py

# Inspect a results.csv:
with open('Scans/2026-04-27_push_edges/results.csv') as f:
    rows = list(csv.DictReader(f))

# Top configs (composite score):
top = evoca_from_scan_top('Scans/2026-04-27_push_edges', top_k=5)
# top[0][1] is a .evoca path; load with evoca_py.import_run

# Top by a single metric:
top_corr = evoca_from_scan_top(
    'Scans/2026-04-27_push_edges', top_k=5,
    score_keys=['correlation_length_mean'])

# Configs near a target:
near = nearest_params(target_config_idx=236,  # cfg236 = scan-3 top
                      n=5, scan_dir='Scans/2026-04-27_push_edges')
```
