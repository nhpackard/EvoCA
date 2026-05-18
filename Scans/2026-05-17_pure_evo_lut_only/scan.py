"""Priority #2 — pure-evolution scan, LUT-only (egene frozen).

Identical methodology to 2026-05-17_pure_evo_joint but with the egene
subsystem frozen: mu_egene = mu_egenome = tax_per_egene = 0, so the
seeded uniform egene never changes. This isolates whether LUT-rule
evolution alone produces excess over the neutral baselines —
specifically dyn_excess_pc_slope (effective rule transitions) and the
component-normalised whole-genome excess_pc_slope.

Grid axes match #1 (minus the egene knobs) so the two scans are
directly comparable and #3 can take #2's LUT winners and re-open the
egene knobs on top of them.

Designed for torque (32 cores). Post-hoc ranking by LUT-evo keys
(eg_* excluded — egene is frozen here).
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

MODE = 'lut_only'       # egene frozen (mu_egene/mu_egenome/tax_per_egene = 0)

PARAM_AXES = {
    'food_inc':      [0.010, 0.013, 0.018],
    'tax':           [0.030, 0.035, 0.040],
    'gdiff':         [0.05, 0.06, 0.08],
    'm_scale':       [1.0, 1.2, 1.5, 2.0, 2.5],
    'mu_lut':        [0.01, 0.03, 0.06, 0.10, 0.15],
    'restricted_mu': [True, False],
}
FROZEN = dict(mu_egene=0.0, mu_egenome=0.0, tax_per_egene=0.0)

N_RANDOM  = 180        # smaller grid (no egene axes) → fewer combos
N_WORKERS = 30
N_STEPS   = 8000
N_GRID    = 256
SAMPLE_EVERY = 100
SEED = 0


def random_grid_sample(axes, n, seed=42):
    rng = np.random.default_rng(seed)
    keys = list(axes.keys())
    sizes = [len(axes[k]) for k in keys]
    total = int(np.prod(sizes))
    seen, out, attempts = set(), [], 0
    while len(out) < min(n, total) and attempts < total * 5:
        idx = int(rng.integers(0, total))
        c, i = {}, idx
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
    p = dict(params, **FROZEN)
    try:
        r = run_sim(p, n_steps=N_STEPS, sample_every=SAMPLE_EVERY,
                    N=N_GRID, seed=SEED, shadow=True, init='halfplane')
        r['config_idx'] = cfg_idx
        r['error'] = ''
    except Exception as e:
        r = {'config_idx': cfg_idx, 'error': repr(e),
             **{f'param_{k}': v for k, v in p.items()}}
    return r


def main():
    save_scan_config(HERE, N=N_GRID, n_steps=N_STEPS,
                     sample_every=SAMPLE_EVERY, seed=SEED, shadow=True)
    configs = random_grid_sample(PARAM_AXES, N_RANDOM)
    print(f"[pure_evo {MODE}] {len(configs)} configs / {N_WORKERS} workers "
          f"/ {N_STEPS} ticks @ N={N_GRID}", flush=True)
    t0 = time.perf_counter()
    args = list(enumerate(configs))
    results = []
    with Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(run_one, args):
            results.append(r)
            d = len(results)
            if d % 25 == 0 or d == len(args):
                dt = time.perf_counter() - t0
                rate = d / dt if dt else 0
                eta = (len(args) - d) / rate if rate else 0
                print(f"  {d}/{len(args)}  {dt:5.0f}s  ~{eta:.0f}s left",
                      flush=True)
    dt = time.perf_counter() - t0
    print(f"Done in {dt:.0f}s ({dt/max(len(results),1):.1f}s avg).")
    out_path = os.path.join(HERE, 'results.csv')
    all_keys = sorted({k for r in results for k in r})
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=all_keys)
        w.writeheader()
        w.writerows(results)
    print(f"Wrote {len(results)} rows -> {out_path}")


if __name__ == '__main__':
    main()
