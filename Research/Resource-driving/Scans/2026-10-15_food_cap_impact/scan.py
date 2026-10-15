"""food_cap impact study (A9 fix, branch `food-cap`, 2026-10-15).

Re-runs a few representative configs from the "productive" m_scale regime
under the legacy uncapped mouthful (food_cap=0, every pre-fix run) and the
corrected capped mouthful (food_cap=1), paired by seed. Deliberately small:
N=128, 8000 ticks, 4 seeds. It answers "does the R1 productive m_scale
interior 2.5-3.5 survive the cap?" and flags what needs a full re-run; it
is not itself that re-run.

Config groups:
  sweep  : pure-evo joint winner base (Research/Activity-metrics/Runs/
           2026-10-15_pure_evo_joint_winner.evoca) with
           m_scale in {1.0,1.5,2.0,2.5,3.0,3.5,5.0} x food_inc in
           {0.010,0.013,0.018}. m_scale=2.5, food_inc=0.018 IS the winner.
  r1_225 : R1 bracket row 225 (top dyn_excess_pc_slope, m_scale 3.5)
  r1_174 : R1 bracket row 174 (top eg_excess, m_scale 3.5)
  lin_3c8: #3c / #8 regime (m_scale 3.0, food_inc 0.015, tax 0.035,
           gdiff 0.06; #8 coevo arm mu_lut 0.03, mu_egene=mu_egenome=0.01)

Per (config, food_cap, seed): run_sim metrics (shadow on, the scan door)
plus a companion headless run (build_sim, same params/seed, shadow off, so
a different RNG trajectory) for food stocks and demography over the second
half: mean F, mean alive f, births per alive cell per step, fraction of
sites with F<0, min F, min alive f.

    python3 scan.py            # writes results.csv next to this file
"""
import csv
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
PYTHON_DIR = os.path.normpath(os.path.join(_repo_root(HERE), 'python'))
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)

from evoca_explore import run_sim, build_sim, save_scan_config   # noqa: E402

N_GRID = 128
N_STEPS = 8000
SAMPLE_EVERY = 100
SEEDS = (0, 1, 2, 3)
N_WORKERS = int(os.environ.get('N_WORKERS', '14'))

WINNER = dict(food_inc=0.018, m_scale=2.5, gdiff=0.08, mu_lut=0.01,
              mu_egene=0.003, mu_egenome=0.003, p_dup_egene=1.0, tax=0.035,
              tax_per_egene=0.0008, tax_lut=0.0, tax_ring=0.0,
              restricted_mu=0)
R1_225 = dict(food_inc=0.013, gdiff=0.08, m_scale=3.5, mu_egene=0.01,
              mu_egenome=0.03, mu_lut=0.01, restricted_mu=1, tax=0.03,
              tax_per_egene=0.0008)
R1_174 = dict(food_inc=0.018, gdiff=0.08, m_scale=3.5, mu_egene=0.003,
              mu_egenome=0.003, mu_lut=0.06, restricted_mu=0, tax=0.03,
              tax_per_egene=0.0008)
LIN_3C8 = dict(m_scale=3.0, food_inc=0.015, tax=0.035, gdiff=0.06,
               mu_lut=0.03, mu_egene=0.01, mu_egenome=0.01, tax_per_egene=0.0)

M_SWEEP = (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 5.0)
FI_SWEEP = (0.010, 0.013, 0.018)


def configs():
    out = []
    for fi in FI_SWEEP:
        for m in M_SWEEP:
            out.append(('sweep', dict(WINNER, m_scale=m, food_inc=fi)))
    out += [('r1_225', R1_225), ('r1_174', R1_174), ('lin_3c8', LIN_3C8)]
    return out


def _companion(params, seed):
    sim = build_sim(params, N=N_GRID, seed=seed)
    half = N_STEPS // 2
    Fm, fm, br, negF = [], [], [], []
    minF, minf = np.inf, np.inf
    for t in range(1, N_STEPS + 1):
        sim.step()
        F = sim.get_F()
        alive = sim.get_alive().astype(bool)
        minF = min(minF, float(F.min()))
        n_alive = int(alive.sum())
        if n_alive:
            minf = min(minf, float(sim.get_f()[alive].min()))
        if t > half:
            Fm.append(float(F.mean()))
            negF.append(float((F < 0).mean()))
            if n_alive:
                fm.append(float(sim.get_f()[alive].mean()))
                br.append(sim.get_births_last() / n_alive)
    final_density = float(sim.get_alive().mean())
    sim.free()
    nan = float('nan')
    return dict(c_F_mean=float(np.mean(Fm)),
                c_f_alive_mean=float(np.mean(fm)) if fm else nan,
                c_birth_rate=float(np.mean(br)) if br else nan,
                c_negF_frac=float(np.mean(negF)),
                c_min_F=minF, c_min_f_alive=minf if np.isfinite(minf) else nan,
                c_density_final=final_density)


def run_one(arg):
    group, params, food_cap, seed = arg
    p = dict(params, food_cap=food_cap)
    try:
        r = run_sim(p, n_steps=N_STEPS, sample_every=SAMPLE_EVERY,
                    N=N_GRID, seed=seed, shadow=True, init='halfplane')
        r.update(_companion(p, seed))
        r['error'] = ''
    except Exception as e:     # keep the row, surface the error
        r = {'error': repr(e), **{f'param_{k}': v for k, v in p.items()},
             'seed': seed}
    r['group'] = group
    return r


def main():
    save_scan_config(HERE, N=N_GRID, n_steps=N_STEPS,
                     sample_every=SAMPLE_EVERY, seed=SEEDS[0], shadow=True)
    # (seed varies per row: see the `seed` column; SEEDS above.)
    args = [(g, p, fc, s) for g, p in configs() for fc in (0, 1)
            for s in SEEDS]
    print(f"[food_cap impact] {len(args)} runs / {N_WORKERS} workers / "
          f"{N_STEPS} ticks @ N={N_GRID}", flush=True)
    t0 = time.perf_counter()
    rows = []
    with Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(run_one, args):
            rows.append(r)
            if len(rows) % 24 == 0 or len(rows) == len(args):
                print(f"  {len(rows)}/{len(args)}  "
                      f"{time.perf_counter() - t0:5.0f}s", flush=True)
    keys = sorted({k for r in rows for k in r})
    with open(os.path.join(HERE, 'results.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"done in {time.perf_counter() - t0:.0f}s; "
          f"errors: {sum(1 for r in rows if r.get('error'))}", flush=True)


if __name__ == '__main__':
    main()
