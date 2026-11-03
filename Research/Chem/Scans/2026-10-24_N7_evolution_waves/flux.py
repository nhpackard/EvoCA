"""CHEM-N7b: activity flux (N, 2026-10-24: the measure that does not reward
long-lived genomes), rerun 2026-11-03 with the uncapped flux.

The first run (2026-10-24) was discarded: the flux summed a 65536-sample
buffer and pinned at 8192 everywhere. Fixed in C/evoca.c
(evoca_activity_flux_sum, no buffer); tests/test_activity_flux.py.

The flux (EvoCA.activity_flux, the light-blue/teal trace of the controls
time-series probe) is the population of LUT-genome buckets whose cumulative
activity crosses the p20-p30 band of the live distribution, per tick over the
last 8 updates. It needs the activity table updated every tick, as the GUI
does, so this is a separate minimal runner: no shadows, activity updated
every tick, flux sampled every 8 ticks over the second half.

In the recipe (mu_lut 0.03, restricted) nearly every birth is a new genome
(D/pop = 1), and a singleton's activity is just its age, so the flux then
counts genomes aging through the band, about D/10. Two things separate that
from clones rising:
  - flux_clonal: the part of the flux carried by buckets with pop >= 2 (each
    crossing sample is the bucket's population at that tick);
  - part B: a mu_lut sweep, to find where genomes have populations.

Part A: the five N7 conditions x 4 seeds, recipe mu_lut, 20000 steps.
Part B: W1 and C0 x mu_lut {0.01, 0.003, 0.001} x 2 seeds, 20000 steps.

    python3 flux.py            # both parts -> flux.csv, progress to stdout
"""
import csv
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scan                                                  # noqa: E402  (params_for, CONDITIONS, SEEDS, N)
from evoca_explore import build_sim                          # noqa: E402

N_STEPS = 20000
HALF = N_STEPS // 2
EVERY = 8
N_WORKERS = 12
MU_SWEEP = (0.01, 0.003, 0.001)

TASKS = ([('A', c, None, s) for c in scan.CONDITIONS for s in scan.SEEDS]
         + [('B', c, mu, s) for c in ('W1', 'C0') for mu in MU_SWEEP for s in (0, 1)])


def run_one(task):
    part, cond, mu, seed = task
    t0 = time.perf_counter()
    p = scan.params_for(cond)
    if mu is not None:
        p['mu_lut'] = mu
    sim = build_sim(p, N=scan.N, seed=seed, init='halfplane')
    fl, fl_cl, fl_pc, births, alive_n, dpop, pop_max = [], [], [], 0, 0.0, [], []
    for t in range(1, N_STEPS + 1):
        sim.step()
        sim._lib.evoca_activity_update()
        if t > HALF:
            births += sim.get_births_last()
            if t % EVERY == 0:
                a = float(sim.get_alive().sum())
                c = sim.activity_crossings(window=EVERY)
                f = sim.activity_flux(window=EVERY)
                assert abs(f - c.sum() / EVERY) <= 1e-6 * max(f, 1.0)
                fl.append(f)
                fl_cl.append(float(c[c >= 2].sum()) / EVERY)
                fl_pc.append(f / a if a else np.nan)
                alive_n += a
                if t % 400 == 0:
                    dpop.append(sim.n_distinct_genomes() / a if a else np.nan)
                    G = sim.get_activity()
                    pop_max.append(int(G['pop_count'].max()) if len(G['pop_count']) else 0)
    sim.free()
    n_samp = len(fl)
    bpt = births / HALF
    return {'part': part, 'condition': cond, 'mu_lut': p['mu_lut'], 'seed': seed,
            'flux_mean': float(np.mean(fl)), 'flux_sd': float(np.std(fl)),
            'flux_clonal': float(np.mean(fl_cl)),
            'clonal_frac': float(np.mean(fl_cl) / np.mean(fl)) if np.mean(fl) else np.nan,
            'flux_per_alive': float(np.nanmean(fl_pc)),
            'births_per_tick': bpt,
            'flux_per_birth': float(np.mean(fl)) / bpt if births else np.nan,
            'alive_mean': alive_n / n_samp / scan.N ** 2,
            'distinct_per_alive': float(np.nanmean(dpop)),
            'top_genome_pop': float(np.mean(pop_max)),
            'wall_s': time.perf_counter() - t0}


def main():
    print(f"Running {len(TASKS)} flux runs across {N_WORKERS} workers...", flush=True)
    rows = []
    with Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(run_one, TASKS):
            rows.append(r)
            print(f"  {r['part']} {r['condition']} mu {r['mu_lut']} s{r['seed']}: flux {r['flux_mean']:.0f}"
                  f"  clonal {r['clonal_frac']:.2f}  per birth {r['flux_per_birth']:.2f}"
                  f"  D/pop {r['distinct_per_alive']:.2f}  top {r['top_genome_pop']:.0f}"
                  f"  {r['wall_s']:.0f}s", flush=True)
    rows.sort(key=lambda r: (r['part'], r['condition'], r['mu_lut'], r['seed']))
    with open(os.path.join(HERE, 'flux.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print('done', flush=True)


if __name__ == '__main__':
    main()
