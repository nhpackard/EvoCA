# EvoCA Parameter Scans

Each subdirectory `<YYYY-MM-DD>_<name>/` is one scan campaign:

```
<date>_<name>/
  scan.py            — driver: builds param grid, fans out to multiprocessing.Pool
  scan_config.json   — fixed init parameters used by every config in this scan
                       (lut, egenome, v_density, alive layout, N, n_steps, …)
  results.csv        — output: one row per config, columns are summary metrics
  notes.md           — hypothesis, axes, post-hoc observations
```

The harness module is `python/evoca_explore.py`. `run_sim(params, ...)` runs
one configuration headless and returns a dict of summary metrics. Drivers
import that and feed configurations through it; results are gathered to
CSV.

## Running a scan from the shell

```bash
cd Scans/<date>_<name>
python3 scan.py
```

## Running a scan from a notebook

```python
import sys, os
SCAN_DIR = 'Scans/2026-04-27_initial'
sys.path.insert(0, 'python')                    # make harness importable
sys.path.insert(0, SCAN_DIR)                    # make scan.py importable
import scan                                     # defines main()
results = scan.main()                           # runs to completion (~few minutes)
```

Or if you want to run a single configuration interactively without the
full multiprocessing apparatus:

```python
import sys
sys.path.insert(0, 'python')
from evoca_explore import run_sim

metrics = run_sim(
    params={'food_inc': 0.015, 'm_scale': 0.6, 'gdiff': 0.05,
            'mu_lut': 0.003, 'mu_egenome': 0.01,
            'tax': 0.035, 'restricted_mu': True},
    n_steps=5000, sample_every=100, N=256, seed=0,
    shadow=True, init='halfplane',
)
print(metrics)
```

## Inspecting scan results from a notebook

```python
import csv
with open('Scans/2026-04-27_initial/results.csv') as f:
    rows = list(csv.DictReader(f))
# rows[i] is one config's metrics; param_food_inc, correlation_length_mean, etc.
```

For ranking and visualisation, the helpers `evoca_from_scan` and
`evoca_from_scan_top` pull rows out and emit `.evoca` recipe files
(see below).

## Promoting a scan row to an interactive sim

`evoca_from_scan(scan_dir, config_idx, ...)` reads one CSV row plus the
scan's `scan_config.json`, builds a recipe dict, and writes a `.evoca`
file to `Runs/`. That file can then be loaded by `import_run()` and
launched in the SDL window:

```python
import sys
sys.path.insert(0, 'python')
from evoca_explore import evoca_from_scan, evoca_from_scan_top
from evoca_py import import_run
from controls import run_with_controls

# Pick by config_idx (matches the 'config_idx' column in results.csv)
path = evoca_from_scan('Scans/2026-04-27_initial', config_idx=70,
                       descriptor='large_patch_dynamics',
                       probes={'activity': True, 'q_activity': True,
                               'n_activity': True, 'env_food': True,
                               'priv_food': True, 'ts': True},
                       colormode=4)
sim, kw = import_run(path)
run_with_controls(sim, **kw)
```

Or grab the top-K by composite score automatically:

```python
top = evoca_from_scan_top('Scans/2026-04-27_initial', top_k=5,
                          probes={...}, colormode=4)
# top is a list of (config_idx, filepath); index 0 = highest score
sim, kw = import_run(top[0][1])
run_with_controls(sim, **kw)
```

## Composite score used by `evoca_from_scan_top`

After filtering out runs that went extinct or saturated above 95%
alive density, four metrics are normalised to `[0, 1]` across the
remaining candidates and **unweighted-averaged**:

| metric | what it captures |
|--------|------------------|
| `correlation_length_mean` | mean spatial scale of the alive field over the run (large = big patches) |
| `largest_patch_temporal_std` | temporal variation of the biggest connected alive component (large = patches keep changing size) |
| `n_distinct_genomes_temporal_std` | temporal fluctuation in G-bucket count (large = genomic turnover, not locked in) |
| `unique_top_genomes` | how many distinct genomes ever held the top-pop spot across samples (large = dominance changes hands) |

```
score(row) = mean(normalised(correlation_length_mean),
                  normalised(largest_patch_temporal_std),
                  normalised(n_distinct_genomes_temporal_std),
                  normalised(unique_top_genomes))
```

Each component contributes equally. Pass `score_keys=[...]` to
`evoca_from_scan_top` to use a different combination of metrics, e.g.
`score_keys=['correlation_length_mean']` for pure spatial-scale ranking.

The score is **purely a ranking helper, not a fitness function**. Two
configs at score 0.75 may be in very different qualitative regimes:
in the 2026-04-27 initial scan, the top-1 (cfg49) lived at
`gdiff=0.05, food_inc=0.025, tax=0.035` while top-3 (cfg70) lived at
`gdiff=0.01, food_inc=0.015, tax=0.020` — different physics, both
producing large-scale dynamics. Run several top-K configs, look at
each, and decide which dynamics you want.
