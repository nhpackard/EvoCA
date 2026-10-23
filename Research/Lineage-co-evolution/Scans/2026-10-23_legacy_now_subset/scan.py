"""2026-10-23. Attribution check: rerun subsets of #3 and #3c with TODAY's code under the
legacy rule (food_cap=0), to separate the rule change from other model drift
(shadow RNG stream, 2026-10-15) between the recorded scans and the rescans."""
import os, sys, importlib.util, json, time
import multiprocessing as mp
import pandas as pd, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
WT = HERE
while not os.path.exists(os.path.join(WT, 'C', 'evoca.h')):
    WT = os.path.dirname(WT)
SC = WT + '/Research/Lineage-co-evolution/Scans'
sys.path.insert(0, WT + '/python')

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

m3c = load('d3c', SC + '/2026-10-23_3c_coevolution_substrate_standard/scan.py')
import evoca_py
def _init_legacy(self, *a, **k):
    k['food_cap'] = 0
    return m3c._evoca_init_orig(self, *a, **k)
evoca_py.EvoCA.init = _init_legacy

m3 = load('d3', SC + '/2026-10-23_pure_evo_egene_on_winners_standard/scan.py')
import evoca_explore as ee
CFG3 = m3.build_configs()

def job(j):
    t0 = time.time()
    if j[0] == '3c':
        _, kind, seed, p1 = j
        rows = m3c.run_one((-1, {'kind': kind, 'seed': seed, 'p1mu': p1}))
        for r in rows: r['src'] = '3c'
        return rows
    _, idx = j
    mode, bi, p = CFG3[idx]
    p = dict(p); p['food_cap'] = 0
    r = ee.run_sim(p, n_steps=8000, sample_every=100, N=256, seed=0, shadow=True, init='halfplane')
    r.update(src='3', config_idx=idx, mode=mode, base_cfg=bi, sec=time.time()-t0)
    return [r]

if __name__ == '__main__':
    old3 = pd.read_csv(SC + '/2026-05-17_pure_evo_egene_on_winners/results.csv')
    # 3b configs from bases 5, 7, 11 (recorded 100% extinct) + their 3a twins, first 4 egene combos each
    sel = []
    for b in (5, 7, 11):
        d = old3[(old3.base_cfg == b)].sort_values('config_idx')
        sel += list(d.config_idx[:8])
    jobs = [('3', int(i)) for i in sel]
    jobs += [('3c', 'B', s, 0.0) for s in range(4)]
    jobs += [('3c', 'AC', s, 0.01) for s in range(3)]
    print(len(jobs), 'jobs', flush=True)
    out = []
    with mp.get_context('fork').Pool(14) as pool:
        for rows in pool.imap_unordered(job, jobs):
            out.extend(rows); print(len(out), flush=True)
    pd.DataFrame(out).to_csv(os.path.join(HERE, 'results.csv'), index=False)
    print('done')
