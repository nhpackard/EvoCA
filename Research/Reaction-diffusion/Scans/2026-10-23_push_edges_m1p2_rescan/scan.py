"""
Rescan of the 65 rows of 2026-04-27_push_edges with m_scale = 1.2, the
only 04-27 rows where the capped mouthful can bind (m_scale > 1). The
04-27 scans ran the legacy uncapped rule; the capped rule is the standard
model (N, 2026-10-23). Same protocol as push_edges (N=256, 5000 steps,
sample_every 25, seed 0, shadow on, halfplane start), today's code, each
row under both rules: food_cap=1 replaces the row in the Reaction-diffusion
notebook's candidate pool; food_cap=0 is the control that separates the
rule's effect from the RNG-stream change since April (today's code does not
reproduce the April metrics bit for bit even where the cap cannot bind).
config_idx is the push_edges config_idx.
"""

import csv
import os
import sys
import time
from multiprocessing import Pool

import pandas as pd

def _repo_root(p):
    """Walk up to the repo root (marker C/evoca.h) — depth-independent."""
    while not os.path.exists(os.path.join(p, 'C', 'evoca.h')):
        if os.path.dirname(p) == p:
            raise FileNotFoundError('EvoCA repo root not found')
        p = os.path.dirname(p)
    return p


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = _repo_root(HERE)
PYTHON_DIR = os.path.join(ROOT, 'python')
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)

from evoca_explore import run_sim, save_scan_config   # noqa: E402
from evoca_nb import params_from_row                   # noqa: E402

SOURCE = os.path.join(HERE, '..', '2026-04-27_push_edges', 'results.csv')
M_SCALE = 1.2
N_WORKERS = 10
N_STEPS   = 5000
N_GRID    = 256
SAMPLE_EVERY = 25
SEED = 0


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
    src = pd.read_csv(SOURCE).set_index('config_idx')
    src = src[src.param_m_scale == M_SCALE]
    args = []
    for cap in (1, 0):
        for cfg in src.index:
            p = params_from_row(src.loc[cfg])        # mu_egenome column read as mu_egene
            p['food_cap'] = cap
            args.append((int(cfg), p))
    print(f"Running {len(args)} runs ({len(src)} rows x 2 rules) across {N_WORKERS} workers...")
    t0 = time.perf_counter()
    results = []
    with Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(run_one, args):
            results.append(r)
            if len(results) % 25 == 0 or len(results) == len(args):
                print(f"  {len(results)}/{len(args)}  {time.perf_counter() - t0:5.0f}s", flush=True)
    results.sort(key=lambda r: (-int(r.get('param_food_cap', 0)), r['config_idx']))
    out_path = os.path.join(HERE, 'results.csv')
    all_keys = sorted({k for r in results for k in r})
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=all_keys)
        w.writeheader()
        w.writerows(results)
    print(f"Wrote {len(results)} rows to {out_path}")


if __name__ == '__main__':
    main()
