"""Priority #3 — egene scan on the LUT-only (#2) winners. BOTH modes.

Reads ../2026-05-17_pure_evo_lut_only/results.csv, ranks by the
LUT-only pure-evo composite, takes the top-12 distinct base regimes,
and for each runs an egene-knob sweep in two modes:

  3a seeded-joint : keep the winner's full params (incl. its winning
                    mu_lut/m_scale), additionally vary the egene knobs.
                    → does cognition on top of LUT-pure-optimal regimes
                      beat the co-optimised #1?
  3b frozen-LUT   : same base but mu_lut=0 (LUT frozen at the winning
                    regime's other params); vary egene knobs only.
                    → clean isolation: does cognition evolve / produce
                      eg-excess on a fixed LUT-pure substrate?

The 3a/3b pair is the scientifically strongest readout (interaction
answer + clean control). Each result row carries 'mode' and
'base_cfg' so the pair can be compared against #1/#2.

torque (32 cores), same run length as #1/#2 for comparability.
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

LUT_ONLY_CSV = os.path.normpath(os.path.join(
    HERE, '..', '2026-05-17_pure_evo_lut_only', 'results.csv'))

N_BASE = 12
EGENE_AXES = {
    'mu_egene':      [0.003, 0.01],
    'mu_egenome':    [0.003, 0.01, 0.03],
    'tax_per_egene': [0.0, 8e-4],
}
BASE_KEYS = ['food_inc', 'tax', 'gdiff', 'm_scale', 'mu_lut',
             'restricted_mu']

N_WORKERS = 30
N_STEPS = 8000
N_GRID = 256
SAMPLE_EVERY = 100
SEED = 0


def _num(s, d=0.0):
    try:
        return float(s)
    except (TypeError, ValueError):
        return d


def top_lut_only_bases(n):
    rows = list(csv.DictReader(open(LUT_ONLY_CSV)))
    cand = [r for r in rows
            if r.get('extinct', '').lower() == 'false'
            and _num(r.get('alive_density_mean', 0)) <= 0.95
            and not r.get('error')]
    keys = ['excess_pc_slope', 'dyn_excess_pc_slope',
            'n_distinct_genomes_temporal_std', 'unique_top_genomes']
    A = {k: np.array([_num(r.get(k, 0)) for r in cand]) for k in keys}

    def nrm(a):
        a = np.clip(a, 0, None)
        return ((a - a.min()) / (a.max() - a.min())
                if a.max() > a.min() else np.zeros_like(a))
    score = sum(nrm(A[k]) for k in keys) / len(keys)
    order = np.argsort(-score)
    bases, seen = [], set()
    for i in order:
        r = cand[i]
        b = {}
        for k in BASE_KEYS:
            v = r.get(f'param_{k}', '')
            b[k] = (v in ('True', 'true', True)) if k == 'restricted_mu' \
                else float(v)
        t = tuple(b[k] for k in BASE_KEYS)
        if t in seen:
            continue
        seen.add(t)
        bases.append(b)
        if len(bases) >= n:
            break
    return bases


def egene_combos():
    ks = list(EGENE_AXES)
    sizes = [len(EGENE_AXES[k]) for k in ks]
    out = []
    for idx in range(int(np.prod(sizes))):
        c, i = {}, idx
        for k, s in zip(ks, sizes):
            c[k] = EGENE_AXES[k][i % s]
            i //= s
        out.append(c)
    return out


def build_configs():
    bases = top_lut_only_bases(N_BASE)
    combos = egene_combos()
    jobs = []
    for bi, base in enumerate(bases):
        for eg in combos:
            for mode in ('3a', '3b'):
                p = dict(base, **eg)
                if mode == '3b':
                    p['mu_lut'] = 0.0          # freeze LUT
                jobs.append((mode, bi, p))
    return jobs


def run_one(arg):
    cfg_idx, (mode, base_i, params) = arg
    try:
        r = run_sim(params, n_steps=N_STEPS, sample_every=SAMPLE_EVERY,
                    N=N_GRID, seed=SEED, shadow=True, init='halfplane')
        r['error'] = ''
    except Exception as e:
        r = {'error': repr(e),
             **{f'param_{k}': v for k, v in params.items()}}
    r['config_idx'] = cfg_idx
    r['mode'] = mode
    r['base_cfg'] = base_i
    return r


def main():
    save_scan_config(HERE, N=N_GRID, n_steps=N_STEPS,
                     sample_every=SAMPLE_EVERY, seed=SEED, shadow=True)
    jobs = build_configs()
    print(f"[#3 egene-on-winners] {len(jobs)} configs "
          f"({N_BASE} bases x {len(egene_combos())} egene x 2 modes) "
          f"/ {N_WORKERS} workers / {N_STEPS} ticks", flush=True)
    t0 = time.perf_counter()
    args = list(enumerate(jobs))
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
    print(f"Done in {dt:.0f}s.")
    out_path = os.path.join(HERE, 'results.csv')
    all_keys = sorted({k for r in results for k in r})
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=all_keys)
        w.writeheader()
        w.writerows(results)
    print(f"Wrote {len(results)} rows -> {out_path}")


if __name__ == '__main__':
    main()
