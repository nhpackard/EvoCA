"""
CHEM-N1: wave and spiral survey under two-way waste exchange (CHEM-D4).

Recipe 2026-10-23_chem_fw_spiralpoint (food/waste, standard rule, mutation
on, halfplane start) with waste_exchange=1, 6000 steps:

  grid (N=512): waste_sink {0.013, 0.015, 0.017, 0.019} x wdiff {0.01, 0.03,
        0.08} x waste_tox {1.0, 0.6, 0.4} x seeds {0, 1, 2}
  large (N=1024, seed 0): wdiff {0.01, 0.08} at waste_sink 0.017, tox 1.0
  control (N=512, one-way, seeds 0-2): waste_sink 0.015, wdiff 0.01, tox 1.0

Diagnostics (python/evoca_waves.py), on 8 x 8 blocks:
  - phase relative to each block's time-mean state over steps 2000-3000;
  - Kuramoto R (synchronized bulk pulsing near 1; waves/spirals low);
  - phase defects every 60 steps over 3000-6000, and their survival: the
    fraction of defects at t + lag with a same-sign defect within 1 block
    at t, for lag 60, 120, 240, 480 (half a period to four periods), against
    the chance level (the earlier frame shifted by half the lattice), and the
    excess (observed - chance) / (1 - chance).
    Noise defects decay to chance; spiral cores stay above it at long lags.
Seed-0 runs save a space-time plot (middle row, colormode 6, steps
5400-6000) and the colormode-6 frame at t=6000 (kymo/, not committed).
"""

import csv
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np


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

import evoca_waves as ew                                     # noqa: E402
from evoca_explore import build_sim, save_scan_config        # noqa: E402

RECIPE = os.path.join(ROOT, 'Research', 'Chem', 'Runs', '2026-10-23_chem_fw_spiralpoint.evoca')
N_WORKERS = 28
N_STEPS = 6000
CENTER_FROM, STATS_FROM = 2000, 3000
SAMPLE_EVERY = 10
DEFECT_EVERY = 60
LAGS = (60, 120, 240, 480)
MATCH_TOL = 1                  # blocks (Chebyshev); defects are dense, so keep it tight
BLOCK = 8
KYMO_STEPS = 600

SINKS = [0.013, 0.015, 0.017, 0.019]
WDIFFS = [0.01, 0.03, 0.08]
TOXES = [1.0, 0.6, 0.4]
SEEDS = [0, 1, 2]


def rgb6(sim, N):
    c = sim.colorize(np.zeros(N * N, np.int32), 6).view(np.uint32).reshape(N, N)
    return np.dstack([(c >> 16) & 255, (c >> 8) & 255, c & 255]).astype(np.uint8)


def run_one(task):
    group, N, seed, over = task
    params = dict(json.load(open(RECIPE))['metaparams_init'])
    params.update(over)
    row = {'group': group, 'N': N, 'seed': seed, **{f'param_{k}': v for k, v in params.items()}}
    try:
        sim = build_sim(params, N=N, seed=seed, init='halfplane')
        Fsum = Wsum = 0.0; ncen = 0
        center = None
        alive_s, R_s, amp_s = [], [], []
        frames = {}                      # t -> (positions, signs)
        n_def, tox, deaths = [], 0, 0
        kymo = []
        for t in range(1, N_STEPS + 1):
            sim.step()
            if CENTER_FROM < t <= STATS_FROM:
                Fsum = Fsum + ew.block_means(sim.get_F(), BLOCK)
                Wsum = Wsum + ew.block_means(sim.get_W(), BLOCK)
                ncen += 1
                if t == STATS_FROM:
                    center = (Fsum / ncen, Wsum / ncen)
            if t > STATS_FROM:
                tox += sim.get_tox_deaths_last(); deaths += sim._lib.evoca_get_deaths_last()
                if t % SAMPLE_EVERY == 0:
                    alive_s.append(sim.get_alive().mean())
                    th, amp = ew.phase_field(sim.get_F(), sim.get_W(), BLOCK, center=center)
                    R_s.append(ew.order_parameter(th, amp)); amp_s.append(float(amp.mean()))
                    if t % DEFECT_EVERY == 0:
                        frames[t] = ew.cores(th, amp)
                        n_def.append(len(frames[t][0]))
            if seed == 0 and t > N_STEPS - KYMO_STEPS:
                kymo.append(rgb6(sim, N)[N // 2])
            if t % 100 == 0 and sim.get_alive().sum() == 0:
                break
        last = rgb6(sim, N) if seed == 0 else None
        sim.free()
        n = N // BLOCK
        surv = {}
        for lag in LAGS:
            obs = ch = tot = 0
            for t0, (pa, sa) in frames.items():
                if t0 + lag not in frames:
                    continue
                pb, sb = frames[t0 + lag]
                obs += ew.matched(pb, sb, pa, sa, n, tol=MATCH_TOL)
                ch += ew.matched(pb, sb, (pa + n // 2) % n, sa, n, tol=MATCH_TOL)
                tot += len(pb)
            surv[lag] = (obs / tot if tot else float('nan'), ch / tot if tot else float('nan'))
        a = np.asarray(alive_s)
        per, peak = ew.period_of(a, SAMPLE_EVERY) if len(a) else (float('nan'), 0.0)
        row.update({
            'extinct': int(len(a) == 0 or a[-1] == 0),
            'alive_mean': float(a.mean()) if len(a) else 0.0,
            'alive_std': float(a.std()) if len(a) else 0.0,
            'period': per, 'acf_peak': peak,
            'R_mean': float(np.nanmean(R_s)) if R_s else float('nan'),
            'amp_mean': float(np.mean(amp_s)) if amp_s else float('nan'),
            'n_defects': float(np.mean(n_def)) if n_def else float('nan'),
            'toxic_share': tox / deaths if deaths else float('nan'),
            **{f'surv_{lag}': surv[lag][0] for lag in LAGS},
            **{f'surv_{lag}_chance': surv[lag][1] for lag in LAGS},
            # survival above chance, normalized: 0 = chance, 1 = every defect survives
            **{f'surv_{lag}_excess': (surv[lag][0] - surv[lag][1]) / (1 - surv[lag][1])
               if surv[lag][1] < 1 else float('nan') for lag in LAGS},
            'error': '',
        })
        if seed == 0 and kymo:
            os.makedirs(os.path.join(HERE, 'kymo'), exist_ok=True)
            tag = (f"{group}_N{N}_ws{over.get('waste_sink')}_wd{over.get('wdiff')}"
                   f"_tox{over.get('waste_tox')}_x{over.get('waste_exchange')}")
            np.savez_compressed(os.path.join(HERE, 'kymo', tag + '.npz'),
                                kymo=np.array(kymo), last=last)
    except Exception as e:
        row['error'] = repr(e)
    return row


def tasks():
    out = []
    for ws in SINKS:
        for wd in WDIFFS:
            for tox in TOXES:
                for s in SEEDS:
                    out.append(('grid', 512, s, dict(waste_exchange=1, waste_sink=ws, wdiff=wd, waste_tox=tox)))
    for wd in (0.01, 0.08):
        out.append(('large', 1024, 0, dict(waste_exchange=1, waste_sink=0.017, wdiff=wd, waste_tox=1.0)))
    for s in SEEDS:
        out.append(('control', 512, s, dict(waste_exchange=0, waste_sink=0.015, wdiff=0.01, waste_tox=1.0)))
    out.sort(key=lambda t: -t[1])        # longest first
    return out


def main():
    cfg_path = save_scan_config(HERE, N=512, n_steps=N_STEPS, sample_every=SAMPLE_EVERY,
                                seed=0, shadow=False)
    cfg = json.load(open(cfg_path))
    cfg.update({'N': 'per row: 512 (grid, control), 1024 (large)', 'seed': 'per row: 0, 1, 2',
                'recipe': os.path.relpath(RECIPE, HERE)})
    json.dump(cfg, open(cfg_path, 'w'), indent=2)
    ts = tasks()
    print(f"Running {len(ts)} runs across {N_WORKERS} workers...", flush=True)
    t0 = time.perf_counter()
    results = []
    with Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(run_one, ts):
            results.append(r)
            if len(results) % 10 == 0 or len(results) == len(ts):
                print(f"  {len(results)}/{len(ts)}  {time.perf_counter() - t0:5.0f}s", flush=True)
    results.sort(key=lambda r: (r['group'], r.get('param_waste_sink', 0), r.get('param_wdiff', 0),
                                -r.get('param_waste_tox', 0), r['seed']))
    keys = sorted({k for r in results for k in r})
    with open(os.path.join(HERE, 'results.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(results)
    print(f"Wrote {len(results)} rows; errors: {sum(bool(r['error']) for r in results)}")


if __name__ == '__main__':
    main()
