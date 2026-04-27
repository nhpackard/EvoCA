"""
Third-iteration scan: push gdiff and m_scale beyond the refined-scan
upper edges (0.08 / 0.8), where the refined-scan top configs sat. Also
quadruples the sampling rate (sample_every: 100 → 25) so the
unique_top_genomes metric un-saturates and discriminates again.
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

from evoca_explore import run_sim, save_scan_config   # noqa: E402


PARAM_AXES = {
    'food_inc':      [0.008, 0.010, 0.013, 0.018],
    'tax':           [0.030, 0.035, 0.040],
    'gdiff':         [0.06, 0.08, 0.12, 0.15, 0.20],
    'mu_lut':        [0.001, 0.003, 0.010, 0.030],
    'mu_egenome':    [0.001, 0.010],
    'm_scale':       [0.7, 0.8, 1.0, 1.2],
    'restricted_mu': [True, False],
}

N_SAMPLES = 250
N_WORKERS = 12
N_STEPS   = 5000
N_GRID    = 256
SAMPLE_EVERY = 25
SEED = 0


def random_grid_sample(axes, n, seed=42):
    rng = np.random.default_rng(seed)
    keys = list(axes.keys())
    sizes = [len(axes[k]) for k in keys]
    total = int(np.prod(sizes))
    n = min(n, total)
    idxs = rng.choice(total, size=n, replace=False)
    out = []
    for raw in idxs:
        cfg = {}
        i = int(raw)
        for k, sz in zip(keys, sizes):
            cfg[k] = axes[k][i % sz]
            i //= sz
        out.append(cfg)
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

    configs = random_grid_sample(PARAM_AXES, N_SAMPLES)
    print(f"Running {len(configs)} configs across {N_WORKERS} workers "
          f"({N_STEPS} ticks at N={N_GRID}, sample_every={SAMPLE_EVERY})...")

    t0 = time.perf_counter()
    args = list(enumerate(configs))
    results = []
    with Pool(N_WORKERS) as pool:
        for i, r in enumerate(pool.imap_unordered(run_one, args)):
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
