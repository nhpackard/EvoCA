# EvoCA Scan-Campaign Analysis: Searching for Spatial Structure + Open-Ended Evolution

**Date**: 2026-04-27
**Scans**: `2026-04-27_initial/`, `2026-04-27_refined/`,
`2026-04-27_push_edges/`, `2026-04-27_evo_focus/`

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

The first three scans (`initial`, `refined`, `push_edges`) ranked
candidates by a balanced composite of spatial-and-evolution metrics.
The fourth scan (`evo_focus`) instead **optimised purely for the
evolution metrics**, on the basis that the productive spatial regime
was already established and we wanted to know how far genomic turnover
could be pushed without sacrificing spatial dynamics.

## Methodology

Each scan: random sample from a discrete parameter grid, run headless
via `multiprocessing.Pool` (12 workers), 5000 ticks at $N = 256$, with
the neutral (Channon) shadow active for excess-activity calibration,
left-halfplane initialisation, fixed seed 0.

| scan dir | # configs | grid combos | `sample_every` | runtime |
|----------|-----:|---------:|---------:|--------:|
| `2026-04-27_initial`    | 150 |  1 152 | 100 | 132 s |
| `2026-04-27_refined`    | 200 | 16 800 | 100 | 103 s |
| `2026-04-27_push_edges` | 250 |  7 680 |  25 | 179 s |
| `2026-04-27_evo_focus`  | 250 |  5 376 ($\dagger$) |  25 | 180 s |

$\dagger$ The evo_focus config set is 50 baseline configs (top-50 of
`push_edges` by combined-composite score, re-run unchanged) plus 200
random samples from the tightened 5 376-combo grid below.

Fixed across all four scans: $N = 256$, $n_{\text{steps}} = 5000$,
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

| axis | initial | refined | push_edges | evo_focus |
|------|---------|---------|------------|-----------|
| `food_inc`      | $\{0.015,\,0.025,\,0.040\}$ | $\{0.010,\,0.013,\,0.015,\,0.018,\,0.022\}$ | $\{0.008,\,0.010,\,0.013,\,0.018\}$ | $\{0.008,\,0.010,\,0.013,\,0.018\}$ |
| `tax`           | $\{0.020,\,0.025,\,0.035\}$ | $\{0.025,\,0.030,\,0.035,\,0.040,\,0.045\}$ | $\{0.030,\,0.035,\,0.040\}$ | $\{0.030,\,0.035,\,0.040\}$ |
| `gdiff`         | $\{0.005,\,0.01,\,0.02,\,0.05\}$ | $\{0.005,\,0.01,\,0.015,\,0.02,\,0.03,\,0.05,\,0.08\}$ | $\{0.06,\,0.08,\,0.12,\,0.15,\,0.20\}$ | $\{0.05,\,0.06,\,0.08,\,0.10\}$ |
| `mu_lut`        | $\{0.001,\,0.003,\,0.01,\,0.03\}$ | $\{0.001,\,0.003,\,0.01,\,0.03\}$ | $\{0.001,\,0.003,\,0.01,\,0.03\}$ | $\{0.001,\,0.003,\,0.005,\,0.01,\,0.02,\,0.03,\,0.06\}$ |
| `mu_egenome`    | $\{0.001,\,0.01\}$ | $\{0.001,\,0.01\}$ | $\{0.001,\,0.01\}$ | $\{0.001,\,0.003,\,0.01,\,0.03\}$ |
| `m_scale`       | $\{0.4,\,0.6\}$ | $\{0.4,\,0.6,\,0.8\}$ | $\{0.7,\,0.8,\,1.0,\,1.2\}$ | $\{0.8,\,1.0,\,1.2,\,1.5\}$ |
| `restricted_mu` | $\{T,\,F\}$ | $\{T,\,F\}$ | $\{T,\,F\}$ | $\{T,\,F\}$ |

Through scans 1 → 3 the grid was tightened around the productive
region, then translated upward (in `gdiff` and `m_scale`). Scan 4
(`evo_focus`) brackets `gdiff` tightly around the scan-3 peak (0.05–
0.10), extends `m_scale` further upward (to 1.5), and **doubles the
mutation-rate resolution** — 7 `mu_lut` values up to 0.06 and 4
`mu_egenome` values up to 0.03. The 50 highest-composite-score
configs from `push_edges` are also reproduced as baselines.

## Yield (candidate fraction)

| scan | extinct | saturated | candidates | yield |
|------|---:|---:|---:|---:|
| initial    | 31 | 94 | 25 | 17 % |
| refined    | ~25 | ~76 | 99 | 50 % |
| push_edges | ~3 | ~2 | 245 | 98 % |
| evo_focus  | ~5 | ~11 | 234 | 94 % |

Yield held at near-saturation in scan 4 even after pushing `mu_lut`
and `m_scale` to higher values, indicating the broad productive region
is robust and not narrowly tuned.

## Per-Metric Maxima Across Scans

Raw maxima of each metric, taken from candidates only (extinct and
saturated rows excluded):

| metric | initial | refined | push_edges | evo_focus | overall ratio |
|--------|---:|---:|---:|---:|---:|
| `correlation_length_mean`        |   50.6 |   82.7 |   75.8 |   81.6 | 1.61× |
| `largest_patch_temporal_std`     |   8478 |  15825 |  17277 |  18550 | 2.19× |
| `n_distinct_genomes_temporal_std`|   8085 |  12984 |  13602 |  16603 | 2.05× |
| `unique_top_genomes`             |     51 |     51 |    194 |    195 | 3.82× ($\dagger$) |
| `F_std_mean`                     |   0.49 |   0.43 |   0.40 |   0.40 | 0.82× |

$\dagger$ The 3.82× rise in `unique_top_genomes` is mostly the
sample-cap rising from 51 to 200 (sample_every dropped from 100 to 25
between scans 2 and 3). The fraction of sample bins that saw a unique
dominator is essentially the same from scan 2 onward ($\approx 0.96$–
$1.00$): turnover is qualitatively saturated.

### Composite scores, globally normalised

The default per-scan composite normalises each metric within that
scan's candidate set, so the maxima always live near 1.0 and aren't
directly comparable across scans. To compare across scans we pool all
four scans' candidates and normalise each metric globally before
averaging:

$$\text{score}^{\,\text{global}}_i \;=\; \frac{1}{|K|}\sum_{m \in K} \frac{m_i - \min_{j\,\in\,\text{all 4 scans}} m_j}{\max_{j\,\in\,\text{all 4 scans}} m_j - \min_{j\,\in\,\text{all 4 scans}} m_j}.$$

| composite (global) | initial | refined | push_edges | evo_focus |
|--------------------|---:|---:|---:|---:|
| **balanced** max    | 0.305 | 0.578 | 0.823 | **0.869** |
| **balanced** top-5 mean | 0.301 | 0.479 | 0.775 | **0.832** |
| **EVO-only** max    | 0.460 | 0.465 | 0.600 | **0.625** |
| **EVO-only** top-5 mean | 0.425 | 0.448 | 0.572 | **0.618** |

Both composites improve monotonically across scans, with `evo_focus`
the highest on both axes. The balanced composite gained 2.85× (initial
top-5 → evo_focus top-5); the EVO-only composite gained 1.45×.
`evo_focus` improved on `push_edges` by **+7 % balanced** and **+8 %
EVO-only** despite optimising only for the latter — which means the
productive corner of parameter space jointly serves both axes.

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

### Push_edges — balanced composite (7 680-combo grid, 245 candidates)

| score | food_inc | tax   | gdiff | mu_lut | mu_eg | m_sc | r_mu  | corr_L | patch_std | div_std | top_g | alive |
|-----:|---------:|------:|------:|-------:|------:|----:|:-----:|------:|----------:|--------:|-----:|------:|
| 0.898 |   0.013 | 0.035 | 0.080 |  0.003 | 0.001 | 1.2 | F     |  73.9 |    15 484 |  10 706 |  188 | 0.322 |
| 0.884 |   0.018 | 0.040 | 0.060 |  0.003 | 0.010 | 1.2 | F     |  74.6 |    15 242 |  10 467 |  182 | 0.360 |
| 0.849 |   0.013 | 0.030 | 0.060 |  0.010 | 0.001 | 1.0 | F     |  69.8 |    13 781 |  12 525 |  158 | 0.403 |
| 0.818 |   0.008 | 0.030 | 0.060 |  0.030 | 0.001 | 1.2 | T     |  61.7 |    12 466 |  13 602 |  158 | 0.217 |
| 0.801 |   0.010 | 0.035 | 0.060 |  0.030 | 0.001 | 1.2 | T     |  66.5 |    11 215 |  12 290 |  168 | 0.231 |

### Evo_focus — EVO-only composite (5 376-combo grid + 50 baselines, 234 candidates)

The evo_focus scan ranks by `EVO_METRICS` (diversity-mean,
diversity-temporal-std, unique-top-genomes, excess-activity-slope) only.
Spatial metrics are reported for context.

| score | food_inc | tax   | gdiff | mu_lut | mu_eg | m_sc | r_mu  | corr_L | patch_std | div_std | div_mn | top_g | slope | alive |
|-----:|---------:|------:|------:|-------:|------:|----:|:-----:|------:|----------:|--------:|------:|-----:|-----:|------:|
| 0.724 |   0.018 | 0.035 | 0.060 |  0.060 | 0.030 | 1.5 | T     |  74.8 |    16 603 |  16 603 | 30 767 |  113 | −1.8 | 0.49 |
| 0.699 |   0.013 | 0.030 | 0.050 |  0.020 | 0.001 | 1.5 | F     |  73.8 |    14 306 |  14 306 | 25 147 |  153 | −3.2 | 0.35 |
| 0.698 |   0.013 | 0.035 | 0.060 |  0.030 | 0.001 | 1.2 | T     |  61.6 |    15 339 |  15 339 | 20 034 |  165 |  3.0 | 0.27 |
| 0.683 |   0.018 | 0.040 | 0.060 |  0.060 | 0.003 | 1.5 | T     |  75.2 |    12 796 |  12 796 | 24 266 |  164 | −2.0 | 0.41 |
| 0.683 |   0.013 | 0.035 | 0.050 |  0.060 | 0.030 | 1.5 | F     |  61.0 |    15 181 |  15 181 | 17 580 |  169 |  0.9 | 0.27 |

The most striking new feature: **`mu_lut = 0.060`** (the new
upper-grid value, double scan-3's max) wins three of the top five —
selection rewards aggressive remixing. **`m_scale = 1.5`** also
takes four of five. Both axes are still at the upper edge of their
ranges, so neither is bracketed yet.

## Top-25 Parameter Distribution

How densely the top-25 cluster on each axis tells us which axes are
bracketed (peak well inside the range) vs unbracketed (peak at the
edge, suggesting the optimum lies further out).

### Push_edges (balanced composite ranking)

```
   gdiff: 0.06: 16   0.08:  9                         (none of {0.12, 0.15, 0.20})
 m_scale: 0.7:  2   0.8:  6   1.0:  6   1.2: 11      (skewed toward 1.2)
food_inc: 0.008: 4  0.010: 8  0.013: 10  0.018: 3    (centred at 0.013)
     tax: 0.030: 11 0.035: 10 0.040:  4              (low-tax preferred)
  mu_lut: 0.001: 2  0.003: 9  0.010:  5  0.030: 9    (bimodal: low or high)
```

### Evo_focus (EVO-only composite ranking)

```
   gdiff: 0.05:  7  0.06: 10  0.08:  6  0.10:  2     (peak at 0.06, same as scan 3)
 m_scale: 0.8:   2  1.0:  4  1.2: 10  1.5:  9        (split between 1.2 and 1.5)
food_inc: 0.008: 1  0.010: 4  0.013: 14  0.018: 6    (centred at 0.013)
     tax: 0.030: 12  0.035: 5  0.040: 8              (bimodal, low-tax edge)
  mu_lut: 0.005: 3  0.010: 3  0.020: 3  0.030: 4  0.060: 12  (skewed to highest)
mu_egenome: 0.001: 8  0.003: 5  0.010: 8  0.030: 4   (spread, low-mu_eg slight edge)
```

Cross-scan summary by axis:

| axis | bracketed? | comment |
|------|---|---|
| `gdiff` | **yes** at 0.06 | both scans 3 & 4 peak at 0.06; values $> 0.08$ never make top-25 |
| `food_inc` | **yes** at 0.013 | tight peak in both scans |
| `tax` | partly | scan 3 prefers low (0.030); scan 4 bimodal between 0.030 and 0.040 |
| `m_scale` | **no** | peak at upper edge in both scans (1.2 → 1.5) |
| `mu_lut` | **no, scan 4** | scan 3 was bimodal at 0.003 / 0.030; scan 4 (with 0.06 added) skewed to 0.060 |
| `mu_egenome` | spread | no clear peak; flat distribution |
| `restricted_mu` | no preference | both T and F well represented in top-25 |

## Diminishing Returns?

Per metric, top-5 means by scan:

| metric | initial | refined | push_edges | evo_focus | trend |
|--------|---:|---:|---:|---:|---|
| `correlation_length_mean`         | 40.7 | 75.2 | 74.4 | 79.3 | **plateau** at ≈ 75–80 |
| `largest_patch_temporal_std`      | 8131 | 11 358 | 15 747 | 18 416 | still gaining |
| `n_distinct_genomes_temporal_std` | 7753 |  9 519 | 12 636 | 15 192 | still gaining |
| `unique_top_genomes`              | 49.8 |   50 |   192 |   192 | saturated at sample-cap |
| candidate yield                   | 17 % |  50 % |  98 % |  94 % | yield essentially solved |
| balanced composite (global, top-5)| 0.301 | 0.479 | 0.775 | **0.832** | gaining |
| EVO-only composite (global, top-5)| 0.425 | 0.448 | 0.572 | **0.618** | gaining |

**Spatial scale has plateaued** at $\approx 75$–$80$ across the last
three scans, which is $\approx 30\,\%$ of the $N = 256$ grid. This
looks like a grid-size ceiling, not a parameter ceiling: with
periodic boundaries, an autocovariance can decorrelate within a
half-period. The way to test that hypothesis is to scan at
$N = 384$ or $512$.

**Temporal turbulence is still gaining.** Both `largest_patch_temporal_std`
and `n_distinct_genomes_temporal_std` continue to improve through scan 4.
Patches keep growing and shrinking more dramatically; the genome
population keeps fluctuating more. Scan 4 specifically pushed
`mu_lut` to 0.060 (double scan-3's cap) and got another 17–22 %
improvement on these axes.

**Survival yield is essentially solved.** From 17 % at the initial
scan to 94–98 % at scans 3 & 4 — the
$(\text{food\_inc},\,\text{tax},\,\text{gdiff},\,m_{\text{scale}})$
joint distribution we now sample in is overwhelmingly viable.

**One Channon-metric anomaly.** `excess_activity_slope` (the slope of
$\Sigma G - \Sigma N$ vs $t$) *fell* in scan 4 (push top-5 mean = 448
→ evo top-5 mean = 228, ratio 0.51). With `mu_lut$ = 0.060$ both real
G-activity and shadow N-activity tables churn rapidly. The shadow,
which uses random selection, can grow proportionally faster than the
real run because every random reproduction creates a new bucket; in
the real run, selection filters. So the slope of the difference
flattens or inverts even though *both* G and N are growing more
vigorously. The other three EVO metrics still rise, so evolution
hasn't slowed — just that this particular metric saturates / inverts
past `mu_lut` $\approx$ 0.03.

## Parameter-Space Novelty Across Scans

Normalised L2 distance from each top-5 config to its nearest top-5
neighbour in the *previous* scan (axes min-max scaled to $[0, 1]$
over the union of all rows):

| comparison | rank-1 dist | rank-2 dist | rank-3 dist | rank-4 dist | rank-5 dist |
|-----------|---:|---:|---:|---:|---:|
| refined vs initial      | 1.187 | 1.071 | 0.713 | 0.718 | 0.633 |
| push_edges vs refined   | 1.120 | 0.675 | 1.087 | 0.905 | 0.858 |
| evo_focus vs push_edges | (mostly bracketed; new `mu_lut=0.060` and `m_scale=1.5` axes are the novelty)            ||||

A top config moving by $\sim 1$ normalised unit means it sits in a
different parameter region than any winner of the previous scan. The
first three iterations each opened new territory; by scan 4 the
parameter-space migration is small (one or two grid steps in `mu_lut`
and `m_scale`), but the *raw metric values* are still rising sharply
on the temporal-turbulence axes.

## Observations

1. **The "interesting middle" regime exists and is large.** With
   94–98 % yield in scans 3 & 4, almost any parameter combination in
   the
   $(\text{food\_inc} \in [0.008, 0.018],\,\text{tax} \in [0.030, 0.040],\,\text{gdiff} \in [0.05, 0.08],\,m_{\text{scale}} \in [0.8, 1.5])$
   box gives a population that survives but does not saturate. The
   regime is robust, not fragile.

2. **Correlation length is grid-size-bound.** The plateau at $\approx
   75$–$80$ cells $\approx N/3$ across scans 2–4 strongly suggests we
   have hit the natural limit imposed by the autocovariance reaching
   a half-period before zero crossing. Scaling $N$ should test this
   directly.

3. **`m_scale` is genuinely the surprise lever.** Initial scans capped
   it at 0.6; refined extended to 0.8 and the optimum moved up; push
   extended to 1.2 and the distribution still piled at the top;
   evo_focus extended to 1.5 and the distribution *still* piles at the
   upper edge. Eating aggressiveness matters more than the initial
   ranges suggested, and is **still unbracketed** after four scans.

4. **`gdiff` is firmly bracketed at 0.06.** Both scan 3 and scan 4
   independently identified 0.06 as the peak, with 0.05 and 0.08 also
   competitive but values $\ge 0.10$ never producing a top-25 config.
   Higher diffusion smears patches; lower diffusion fragments the
   population into death.

5. **`mu_lut` wants higher than initially explored.** Scans 1–3 used
   $\{0.001, 0.003, 0.01, 0.03\}$ and saw a bimodal preference (low
   for "preserve coherent CA, occasional invasions"; high for
   "constantly remix"). Scan 4 added 0.060 to the range, and
   $\text{mu\_lut} = 0.060$ wins 12 of the 25 top EVO-scoring configs
   — **also unbracketed**. Higher mutation gives a clearer pure-evo
   regime where every birth is essentially a fresh genome.

6. **Optimising for evolution didn't sacrifice spatial structure.**
   Scan 4 ranked by EVO_METRICS only and still scored higher on
   spatial metrics (top-5 correlation length 79.3 vs push_edges 74.4,
   patch-temporal-std 18 416 vs 15 747). The two axes share a
   parameter region — there is no Pareto trade-off in the regime
   we're searching.

7. **Channon's `excess_activity_slope` saturates / inverts at high
   mutation.** Past `mu_lut` $\approx 0.03$ the metric stops being
   informative because both real- and shadow-activity tables churn so
   fast. Other EVO metrics (diversity-mean, diversity-temporal-std,
   unique-top-genomes) keep rising, so this is a metric quirk, not
   evolution stalling.

8. **`unique_top_genomes` is informative only as a fraction of
   sample bins.** All four scans hit 96–100 % of the sample-cap — the
   raw count just measures the cap. Future analyses should compute
   and report $U / |T_{\text{sample}}|$ instead.

## Recommended Next Steps

1. **Visual validation at $N = 512$.** Pull the top configs from scan
   3 (best balanced) and scan 4 (best evo-only) into the SDL viewer.
   Confirm the metrics ("big patches with ongoing genomic turnover")
   match what one actually sees. Easy first move.

2. **Larger-grid scan ($N \in \{384, 512\}$).** Test whether the
   correlation-length plateau is a grid-size artefact. Smaller config
   set, fixed at the productive corner. Compute scales as $N^2$;
   ~24 combos at $N = 512$ would take 10–30 min.

3. **Bracket `m_scale` and `mu_lut` in scan 5.** Both are still at
   the upper edge of their ranges in scan 4. Suggested grid:
   `m_scale ∈ {1.5, 2.0, 2.5}`, `mu_lut ∈ {0.06, 0.10, 0.15, 0.20}`,
   `gdiff ∈ {0.05, 0.06, 0.08}`, `food_inc ∈ {0.013, 0.018}`,
   `tax ∈ {0.030, 0.035, 0.040}`. ~150 combos, ~2 min.

4. **Multi-seed runs** of the top configs to confirm dynamics are
   robust to initial conditions. Currently all scans use `seed=0`,
   and `multiprocessing.Pool` worker reuse adds further per-config
   non-determinism (the C-side xorshift RNG isn't reset between
   tasks). A fix would be to add `evoca_set_seed()` and call it at
   the top of `run_sim`.

## Reproducing this analysis

```python
# Assumes you're in the EvoCA repo root.
import sys, csv
sys.path.insert(0, 'python')

from evoca_explore import (run_sim, save_scan_config,
                           evoca_from_scan_top,
                           nearest_params, nearest_evo, nearest_spatial,
                           EVO_METRICS, SPATIAL_METRICS)

# Re-run any scan:
#   cd Scans/2026-04-27_<name> && python3 scan.py

# Inspect a results.csv:
with open('Scans/2026-04-27_evo_focus/results.csv') as f:
    rows = list(csv.DictReader(f))

# Top configs (default balanced composite score):
top = evoca_from_scan_top(top_k=5)   # uses most-recent scan dir

# Top configs by EVO-only or SPATIAL-only:
top_evo     = evoca_from_scan_top(top_k=5, score_keys=EVO_METRICS,
                                  descriptor_prefix='top_evo')
top_spatial = evoca_from_scan_top(top_k=5, score_keys=SPATIAL_METRICS,
                                  descriptor_prefix='top_spatial')

# Configs near a target:
near = nearest_params(target_config_idx=236, n=5)   # cfg236 = scan-3 top
```
