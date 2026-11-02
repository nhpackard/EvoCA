"""SUPERSEDED (no result): the flux clips at the crossings buffer (65536) and is ~D/10 when
every genome is a singleton, as here; see notes.md "Activity flux".

CHEM-N7b: activity flux (N, 2026-10-24: the measure that does not reward
long-lived genomes) for the five N7 conditions, same seeds, N=512, 20000 steps.

The flux (EvoCA.activity_flux, the light-blue/teal trace of the controls
time-series probe) is the population of LUT-genome buckets whose cumulative
activity crosses the p20-p30 band of the live distribution, per tick over the
last 8 updates. It needs the activity table updated every tick, as the GUI
does; run_sim updates it only at its samples (and its A_new/excess would break
if it were updated per tick, since the neutral shadow is not). So this is a
separate, minimal runner: no shadows, activity updated every tick, flux sampled
every 8 ticks over the second half. Also recorded: births per tick and distinct
genomes per alive cell (D / pop), because when nearly every cell carries its
own genome the flux mostly counts newborns aging through the band."""
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
N_WORKERS = 10


def run_one(task):
    cond, seed = task
    t0 = time.perf_counter()
    sim = build_sim(scan.params_for(cond), N=scan.N, seed=seed, init='halfplane')
    fl, fl_pc, births, alive_n, dpop = [], [], 0, 0.0, []
    for t in range(1, N_STEPS + 1):
        sim.step()
        sim._lib.evoca_activity_update()
        if t > HALF:
            births += sim.get_births_last()
            if t % EVERY == 0:
                a = float(sim.get_alive().sum())
                f = sim.activity_flux(window=EVERY)
                fl.append(f)
                fl_pc.append(f / a if a else np.nan)
                alive_n += a
                if t % 400 == 0:
                    dpop.append(sim.n_distinct_genomes() / a if a else np.nan)
    sim.free()
    n_samp = len(fl)
    return {'condition': cond, 'seed': seed,
            'flux_mean': float(np.mean(fl)), 'flux_sd': float(np.std(fl)),
            'flux_per_alive': float(np.nanmean(fl_pc)),
            'births_per_tick': births / HALF,
            'flux_per_birth': float(np.mean(fl)) / (births / HALF) if births else np.nan,
            'alive_mean': alive_n / n_samp / scan.N ** 2,
            'distinct_per_alive': float(np.nanmean(dpop)),
            'wall_s': time.perf_counter() - t0}


def main():
    tasks = [(c, s) for c in scan.CONDITIONS for s in scan.SEEDS]
    print(f"Running {len(tasks)} flux runs across {N_WORKERS} workers...", flush=True)
    rows = []
    with Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(run_one, tasks):
            rows.append(r)
            print(f"  {r['condition']} s{r['seed']}: flux {r['flux_mean']:.0f}  per birth {r['flux_per_birth']:.2f}"
                  f"  D/pop {r['distinct_per_alive']:.2f}  {r['wall_s']:.0f}s", flush=True)
    rows.sort(key=lambda r: (r['condition'], r['seed']))
    with open(os.path.join(HERE, 'flux.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print('done')


if __name__ == '__main__':
    main()
