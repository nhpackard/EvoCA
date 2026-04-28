"""
Fourth-iteration scan: optimise for evolution metrics only.

Config set:
- Top 50 from Scans/2026-04-27_push_edges (by combined evo+spatial
  composite) re-run as baselines.
- ~200 random variations from a tightened grid around the productive
  corner, with finer granularity on the evolution-affecting axes.

Run conditions match push_edges (5000 ticks, N=256, sample_every=25,
neutral shadow on, halfplane init), so results.csv is directly
comparable to push_edges/results.csv.
"""

import csv
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PYTHON_DIR = os.path.normpath(os.path.join(HERE, '..', '..', 'python'))
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)

from evoca_explore import (run_sim, save_scan_config,   # noqa: E402
                            _load_scan)


SOURCE_SCAN = os.path.normpath(
    os.path.join(HERE, '..', '2026-04-27_push_edges'))

# Tightened grid around the productive corner identified by push_edges,
# with finer mu_lut and mu_egenome since those are the main evolution
# levers.
PARAM_AXES = {
    'food_inc':      [0.008, 0.010, 0.013, 0.018],
    'tax':           [0.030, 0.035, 0.040],
    'gdiff':         [0.05, 0.06, 0.08, 0.10],
    'mu_lut':        [0.001, 0.003, 0.005, 0.010, 0.020, 0.030, 0.060],
    'mu_egenome':    [0.001, 0.003, 0.010, 0.030],
    'm_scale':       [0.8, 1.0, 1.2, 1.5],
    'restricted_mu': [True, False],
}

N_TOP_BASELINE      = 50
N_RANDOM_VARIATIONS = 200
N_WORKERS = 12
N_STEPS   = 5000
N_GRID    = 256
SAMPLE_EVERY = 25
SEED = 0


def _num(s, d=0.0):
    try: return float(s)
    except (TypeError, ValueError): return d


def get_top_baseline_from_source(n=N_TOP_BASELINE):
    """Top-n configs from SOURCE_SCAN by composite (default 4-key) score,
    after standard candidate filtering (not extinct, not saturated).
    Returns a list of param dicts with the right key set/types."""
    cfg, rows = _load_scan(SOURCE_SCAN)
    cands = [r for r in rows
             if r.get('extinct', '').lower() == 'false'
             and _num(r.get('alive_density_mean', 0)) <= 0.95]
    keys = ['correlation_length_mean', 'largest_patch_temporal_std',
            'n_distinct_genomes_temporal_std', 'unique_top_genomes']
    arrs = {k: np.array([_num(r.get(k, 0)) for r in cands]) for k in keys}

    def _norm(a):
        if a.max() == a.min(): return np.zeros_like(a)
        return (a - a.min()) / (a.max() - a.min())

    score = sum(_norm(arrs[k]) for k in keys) / len(keys)
    order = np.argsort(-score)[:n]

    out = []
    for i in order:
        r = cands[i]
        c = {}
        for k in PARAM_AXES.keys():
            v = r[f'param_{k}']
            if k == 'restricted_mu':
                c[k] = (v in ('True', 'true', True))
            else:
                c[k] = float(v)
        out.append(c)
    return out


def random_grid_sample(axes, n, seed=42, exclude=None):
    """Random unique samples from the discrete grid, skipping any tuples
    already in `exclude`."""
    rng = np.random.default_rng(seed)
    keys = list(axes.keys())
    sizes = [len(axes[k]) for k in keys]
    total = int(np.prod(sizes))

    seen = set()
    if exclude:
        for e in exclude:
            seen.add(tuple(e[k] for k in keys))

    out = []
    attempts = 0
    while len(out) < n and attempts < total * 5:
        idx = int(rng.integers(0, total))
        c = {}
        i = idx
        for k, sz in zip(keys, sizes):
            c[k] = axes[k][i % sz]
            i //= sz
        t = tuple(c[k] for k in keys)
        if t in seen:
            attempts += 1
            continue
        seen.add(t)
        out.append(c)
    return out


def run_one(arg):
    cfg_idx, params = arg
    try:
        r = run_sim(params, n_steps=N_STEPS, sample_every=SAMPLE_EVERY,
                    N=N_GRID, seed=SEED, shadow=True, init='halfplane')
        r['config_idx'] = cfg_idx
        r['error'] = ''
    except Exception as e:
        r = {'config_idx': cfg_idx, 'error': repr(e),
             **{f'param_{k}': v for k, v in params.items()}}
    return r


def main():
    save_scan_config(HERE, N=N_GRID, n_steps=N_STEPS,
                     sample_every=SAMPLE_EVERY, seed=SEED, shadow=True)

    print(f"Pulling top {N_TOP_BASELINE} from {SOURCE_SCAN}…")
    baseline = get_top_baseline_from_source(N_TOP_BASELINE)
    print(f"  {len(baseline)} baseline configs")

    print(f"Drawing {N_RANDOM_VARIATIONS} random variations from the "
          f"tightened grid (excluding the baseline)…")
    variations = random_grid_sample(PARAM_AXES, N_RANDOM_VARIATIONS,
                                     exclude=baseline)
    print(f"  {len(variations)} unique variations")

    configs = baseline + variations
    print(f"Total: {len(configs)} configs across {N_WORKERS} workers "
          f"({N_STEPS} ticks at N={N_GRID}, sample_every={SAMPLE_EVERY})")

    t0 = time.perf_counter()
    args = list(enumerate(configs))
    results = []
    with Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(run_one, args):
            results.append(r)
            done = len(results)
            if done % 25 == 0 or done == len(args):
                dt = time.perf_counter() - t0
                rate = done / dt
                eta = (len(args) - done) / rate if rate > 0 else 0
                print(f"  {done}/{len(args)}  {dt:5.0f}s elapsed  "
                      f"~{eta:.0f}s remaining")
    dt = time.perf_counter() - t0
    print(f"Done in {dt:.0f}s ({dt / max(len(results), 1):.1f}s avg).")

    out_path = os.path.join(HERE, 'results.csv')
    all_keys = sorted({k for r in results for k in r})
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=all_keys)
        w.writeheader()
        w.writerows(results)
    print(f"Wrote {len(results)} rows to {out_path}")
    return results


if __name__ == '__main__':
    main()
