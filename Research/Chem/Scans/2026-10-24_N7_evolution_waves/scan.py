"""CHEM-N7: evolution in vs out of the wave regime.

Recipe 2026-10-23_chem_fw_spiralpoint (mutation on), N=512, 20000 steps,
seeds 0-3, five conditions:
  W1  one-way, waste_tox 0.4, waste_sink 0.019, wdiff 0.01   waves (CHEM-N6)
  S1  as W1 with wdiff 0.08                                  no waves (N6 speckle)
  T1  one-way, waste_tox 1.0, waste_sink 0.019, wdiff 0.01   expected no waves
  W2  two-way, waste_tox 1.0, waste_sink 0.015, wdiff 0.08   waves (N1/N6)
  C0  waste off (n_chem = 1)                                 vanilla baseline
Every run records run_sim's evolution measures (activity excess against the
dyn/eg drift shadows and the LUT-hash shadow, A_new, genome diversity,
dominant-genome turnover, egene count and specificity, intake), plus, via a
read-only step_hook:
  - the site period (evoca_waves.site_period, 4-cell coarse-graining) of the
    middle row over 600-step windows ending at 2000, 6000, 10000, 15000 and
    20000: whether the run is in the wave regime, and whether evolution
    makes or breaks it;
  - births per alive cell per step over the second half (turnover);
  - the mean reproduction interval over the second half, from the repro-age
    histogram with t0 = 10000. It sees only reproducing parents, which
    reproduce in bursts when fed, so it is the burst interval; the
    generation time is 1 / births_per_alive_step.
The middle row of the last window is saved (frames/, not committed)."""

import csv
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np


def _repo_root(p):
    while not os.path.exists(os.path.join(p, 'C', 'evoca.h')):
        if os.path.dirname(p) == p:
            raise FileNotFoundError('EvoCA repo root not found')
        p = os.path.dirname(p)
    return p


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = _repo_root(HERE)
sys.path.insert(0, os.path.join(ROOT, 'python'))

import evoca_waves as ew                                     # noqa: E402
from evoca_explore import run_sim, save_scan_config          # noqa: E402

RECIPE = os.path.join(ROOT, 'Research', 'Chem', 'Runs', '2026-10-23_chem_fw_spiralpoint.evoca')
N_WORKERS = 10
N, N_STEPS, SAMPLE_EVERY = 512, 20000, 200
WINDOW_ENDS = (2000, 6000, 10000, 15000, 20000)
WINDOW = 600
HALF = N_STEPS // 2
SEEDS = [0, 1, 2, 3]
WASTE_KEYS = ('wdiff', 'waste_sink', 'waste_disp', 'waste_tox', 'waste_exchange')
CONDITIONS = {
    'W1': dict(waste_exchange=0, waste_tox=0.4, waste_sink=0.019, wdiff=0.01),
    'S1': dict(waste_exchange=0, waste_tox=0.4, waste_sink=0.019, wdiff=0.08),
    'T1': dict(waste_exchange=0, waste_tox=1.0, waste_sink=0.019, wdiff=0.01),
    'W2': dict(waste_exchange=1, waste_tox=1.0, waste_sink=0.015, wdiff=0.08),
    'C0': dict(n_chem=1),
}


def params_for(cond):
    p = dict(json.load(open(RECIPE))['metaparams_init'])
    p.update(CONDITIONS[cond])
    if p.get('n_chem') == 1:
        for k in WASTE_KEYS:
            p.pop(k, None)
    return p


class WaveEvoHook:
    """Read-only: get_alive, get_births_last and the repro-age histogram."""

    def __init__(self):
        self.rows = {e: [] for e in WINDOW_ENDS}
        self.births = 0.0
        self.alive_sum = 0.0

    def start(self, sim):
        sim.set_repro_age_t0(HALF)

    def step(self, sim, t):
        for e in WINDOW_ENDS:
            if e - WINDOW < t <= e:
                self.rows[e].append(sim.get_alive().reshape(N, N)[N // 2].astype(np.uint8))
        if t > HALF:
            self.births += sim.get_births_last()
            self.alive_sum += float(sim.get_alive().sum())
        if t == N_STEPS:
            h = np.asarray(sim.get_repro_age_hist(), float)
            self.repro_interval = float((h * np.arange(len(h))).sum() / h.sum()) if h.sum() else float('nan')

    def summary(self):
        out = {}
        for e in WINDOW_ENDS:
            if len(self.rows[e]) == WINDOW:
                per, pk = ew.site_period(np.array(self.rows[e], float))
            else:
                per, pk = float('nan'), 0.0
            out[f'period_t{e}'], out[f'period_peak_t{e}'] = per, pk
        out['births_per_alive_step'] = self.births / self.alive_sum if self.alive_sum else float('nan')
        out['repro_interval_mean'] = getattr(self, 'repro_interval', float('nan'))
        return out


def run_one(task):
    cond, seed = task
    t0 = time.perf_counter()
    hook = WaveEvoHook()
    try:
        row = run_sim(params_for(cond), n_steps=N_STEPS, sample_every=SAMPLE_EVERY, N=N, seed=seed,
                      shadow=True, init='halfplane', step_hook=hook)
        row['error'] = ''
    except Exception as e:                       # noqa: BLE001
        row = {'error': repr(e)}
    row.update(condition=cond, seed=seed, wall_s=time.perf_counter() - t0)
    last = hook.rows[WINDOW_ENDS[-1]]
    if len(last) == WINDOW:
        os.makedirs(os.path.join(HERE, 'frames'), exist_ok=True)
        np.savez_compressed(os.path.join(HERE, 'frames', f'{cond}_s{seed}.npz'), kymo_alive=np.array(last))
    return row


def main():
    cfg = save_scan_config(HERE, N=N, n_steps=N_STEPS, sample_every=SAMPLE_EVERY, seed=0, shadow=True)
    c = json.load(open(cfg))
    c.update({'seed': '0-3', 'recipe': os.path.relpath(RECIPE, HERE), 'conditions': CONDITIONS})
    json.dump(c, open(cfg, 'w'), indent=2)
    tasks = [(c, s) for c in CONDITIONS for s in SEEDS]
    print(f"Running {len(tasks)} runs across {N_WORKERS} workers...", flush=True)
    t0 = time.perf_counter()
    rows = []
    with Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(run_one, tasks):
            rows.append(r)
            print(f"  {r['condition']} s{r['seed']}: {r['wall_s']:5.0f}s  period_t20000 "
                  f"{r.get('period_t20000')}  {r['error']}  [{time.perf_counter() - t0:5.0f}s]", flush=True)
    rows.sort(key=lambda r: (r['condition'], r['seed']))
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(HERE, 'results.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows; errors: {sum(bool(r['error']) for r in rows)}")


if __name__ == '__main__':
    main()
