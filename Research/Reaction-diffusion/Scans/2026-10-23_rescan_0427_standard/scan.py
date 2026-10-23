"""
Rescan of all three 2026-04-27 RD scans (initial 150, refined 200,
push_edges 250 rows) under today's standard model. The 04-27 rows were
recorded under the legacy uncapped mouthful (food_cap=0) and, more
importantly, under the April model: ternary egenes with +/-1 scoring
(2026-05-10), max-match eating across egene slots and the egene-aware
genome hash (2026-04-28), and the shadow's own RNG stream (2026-10-15)
all came later. On 40 push_edges rows where the cap cannot bind
(m_scale <= 1), today's code gives a lower correlation_length_mean in 93%
of rows (median 30.3 -> 14.3) and fewer unique_top_genomes in 90%, so
the recorded metrics do not describe the current model.

Same protocol as each source scan (N=256, 5000 steps, seed 0, shadow
on, halfplane start; sample_every from that scan's scan_config.json),
params from each row via evoca_nb.params_from_row (mu_egenome column read
as mu_egene) with food_cap=1. Columns `scan` (initial/refined/push_edges)
and config_idx identify the source row.
"""

import csv
import json
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

from evoca_explore import run_sim   # noqa: E402
from evoca_nb import params_from_row  # noqa: E402

SCANS = ('initial', 'refined', 'push_edges')
N_WORKERS = int(os.environ.get('N_WORKERS', 8))


def run_one(arg):
    scan, cfg_idx, params, spec = arg
    try:
        r = run_sim(params, n_steps=spec['n_steps'], sample_every=spec['sample_every'],
                    N=spec['N'], seed=spec['seed'], shadow=spec['shadow'], init='halfplane')
        r['error'] = ''
    except Exception as e:
        r = {'error': repr(e), **{f'param_{k}': v for k, v in params.items()}}
    r['scan'], r['config_idx'] = scan, cfg_idx
    return r


def main():
    args, specs = [], {}
    for s in SCANS:
        d = os.path.join(HERE, '..', f'2026-04-27_{s}')
        spec = json.load(open(os.path.join(d, 'scan_config.json')))
        assert spec['initialization']['alive'] == 'halfplane'
        specs[s] = {k: spec[k] for k in ('N', 'n_steps', 'sample_every', 'seed', 'shadow')}
        src = pd.read_csv(os.path.join(d, 'results.csv')).set_index('config_idx')
        for cfg in src.index:
            p = params_from_row(src.loc[cfg])
            p['food_cap'] = 1
            args.append((s, int(cfg), p, specs[s]))
    json.dump({'source_scans': [f'2026-04-27_{s}' for s in SCANS], 'specs': specs, 'food_cap': 1,
               'created': time.strftime('%Y-%m-%dT%H:%M:%S')},
              open(os.path.join(HERE, 'scan_config.json'), 'w'), indent=2)
    print(f"Running {len(args)} rows across {N_WORKERS} workers...", flush=True)
    t0 = time.perf_counter()
    results = []
    with Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(run_one, args):
            results.append(r)
            if len(results) % 50 == 0 or len(results) == len(args):
                print(f"  {len(results)}/{len(args)}  {time.perf_counter() - t0:5.0f}s", flush=True)
    results.sort(key=lambda r: (SCANS.index(r['scan']), r['config_idx']))
    all_keys = sorted({k for r in results for k in r})
    with open(os.path.join(HERE, 'results.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=all_keys)
        w.writeheader()
        w.writerows(results)
    print(f"Wrote {len(results)} rows")


if __name__ == '__main__':
    main()
