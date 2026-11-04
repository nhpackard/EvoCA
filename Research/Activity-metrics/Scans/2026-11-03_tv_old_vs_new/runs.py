"""Old vs new fixed-space readouts (N, 2026-11-03): the per-component excess
sum G / D_G - sum N / D_N (eg_excess_pc, dyn_excess_pc) against the paired
total-variation distance 1/2 sum_k |G_k - N_k| / sum G (eg_tv, dyn_tv).

Two sets, N = 128, 6000 steps, sampled every 100:
  decomp  the median-eg-slope config of each egene-tax stratum of the
          pure-evo joint rescan (cfg 57 untaxed, cfg 275 taxed), seeds 0-2;
  causal  the causal-control arms of Activity-metrics.ipynb section 3
          (JOINT, LUT_only, EGENE_only) plus FROZEN (no mutation), seeds 0-3.

The dyn TV is also split into its input-marginal and within-input parts
(dyn_tv_split), because the dyn drift shadow spreads mutation mass across
LUT inputs, which mutation cannot do (found here, 2026-11-03).

    python3 runs.py      # -> runs.npz (time series, final per-key arrays), runs.csv
"""
import math
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', '..', 'python'))
from evoca_explore import build_sim  # noqa: E402

N, STEPS, EVERY, N_WORKERS = 128, 6000, 100, 12
SCAN = os.path.join(HERE, '..', '2026-10-23_pure_evo_joint_standard', 'results.csv')

BASE = dict(food_inc=0.013, m_scale=1.2, gdiff=0.06, tax=0.035, tax_per_egene=8e-4)
ARMS = {'JOINT': dict(BASE, mu_lut=0.01, mu_egene=0.003, mu_egenome=0.005),
        'LUT_only': dict(BASE, mu_lut=0.01, mu_egene=0.0, mu_egenome=0.0),
        'EGENE_only': dict(BASE, mu_lut=0.0, mu_egene=0.003, mu_egenome=0.005),
        'FROZEN': dict(BASE, mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0)}
DECOMP = (57, 275)
SERIES = ('t', 'alive', 'eg_excess_pc', 'eg_tv', 'eg_D', 'eg_top', 'dyn_excess_pc', 'dyn_tv', 'dyn_D',
          'dyn_tv_input', 'dyn_tv_within')


def dyn_tv_split(G, Nn):
    """Split the dyn TV (bucket = input * 2 + output) into the input-marginal
    part (how often each of the 250 LUT inputs is used: set by the lattice
    dynamics, which LUT mutation cannot move) and the within-input part
    (which output each input gets, weighted by the real run's input use)."""
    G2, N2 = G.reshape(-1, 2), Nn.reshape(-1, 2)
    Gi, Ni = G2.sum(1), N2.sum(1)
    if Gi.sum() <= 0:
        return 0.0, 0.0
    marg = 0.5 * np.abs(Gi / Gi.sum() - Ni / Ni.sum()).sum()
    cG = G2 / np.where(Gi > 0, Gi, 1)[:, None]
    cN = N2 / np.where(Ni > 0, Ni, 1)[:, None]
    within = 0.5 * ((Gi / Gi.sum())[:, None] * np.abs(cG - cN)).sum()
    return float(marg), float(within)


def scan_params(cfg):
    d = pd.read_csv(SCAN)
    r = d[d.config_idx == cfg].iloc[0]
    return {c[6:]: r[c] for c in d.columns if c.startswith('param_')}


def tasks():
    out = [('decomp', f'cfg{c}', scan_params(c), s) for c in DECOMP for s in (0, 1, 2)]
    out += [('causal', a, p, s) for a, p in ARMS.items() for s in (0, 1, 2, 3)]
    return out


def run_one(task):
    group, name, p, seed = task
    sim = build_sim(p, N=N, seed=seed)
    ts = {k: [] for k in SERIES}
    m_sum = {'eg': 0.0, 'dyn': 0.0}
    m_n = 0
    for t in range(1, STEPS + 1):
        sim.step()
        # the shadows' per-tick mixing weight m = 1 - exp(-rate), rate as in
        # fs_shadow_step (C/evoca.c): flips / (births x bits per component)
        r = sim._lib.evoca_get_repro_last()
        if r > 0:
            m_sum['eg'] += 1 - math.exp(-sim._lib.evoca_get_eg_bitflips_last() / (r * 12.0))
            m_sum['dyn'] += 1 - math.exp(-sim._lib.evoca_get_lut_bitflips_last() / r)
            m_n += 1
        if t % EVERY:
            continue
        sim._lib.evoca_eg_activity_update()
        eg_pop = sim.get_eg_activity()['pop_count']
        dyn_pop = sim.get_dyn_activity()['pop_count']
        dp = sim.get_dyn_shadow()
        tv_in, tv_wi = dyn_tv_split(dp['obs'], dp['shadow'])
        for k, v in (('t', t), ('alive', sim.get_alive().mean()),
                     ('eg_excess_pc', sim.eg_excess_pc()), ('eg_tv', sim.eg_tv()),
                     ('eg_D', int((eg_pop > 0).sum())),
                     ('eg_top', eg_pop.max() / max(eg_pop.sum(), 1)),
                     ('dyn_excess_pc', sim.dyn_excess_pc()), ('dyn_tv', sim.dyn_tv()),
                     ('dyn_D', int((dyn_pop > 0).sum())),
                     ('dyn_tv_input', tv_in), ('dyn_tv_within', tv_wi)):
            ts[k].append(v)
    ts = {k: np.asarray(v, dtype=float) for k, v in ts.items()}
    for ch, pair in (('eg', sim.get_eg_shadow()), ('dyn', sim.get_dyn_shadow())):
        ts[f'{ch}_G_final'], ts[f'{ch}_N_final'] = pair['obs'], pair['shadow']
        ts[f'{ch}_m_mean'] = np.array(m_sum[ch] / m_n if m_n else 0.0)
    sim._lib.evoca_eg_activity_update()
    ts['eg_pop_final'] = sim.get_eg_activity()['pop_count'].astype(float)
    ts['dyn_pop_final'] = sim.get_dyn_activity()['pop_count'].astype(float)
    sim.free()
    return group, name, seed, ts


def main():
    T = tasks()
    print(f'{len(T)} runs on {N_WORKERS} workers', flush=True)
    res = {}
    with Pool(N_WORKERS) as pool:
        for group, name, seed, ts in pool.imap_unordered(run_one, T):
            res[(group, name, seed)] = ts
            print(f'  {group} {name} s{seed}: eg_tv {ts["eg_tv"][-1]:.3f} dyn_tv {ts["dyn_tv"][-1]:.3f}'
                  f'  eg_excess {ts["eg_excess_pc"][-1]:.0f}', flush=True)
    np.savez_compressed(os.path.join(HERE, 'runs.npz'),
                        **{f'{g}|{n}|{s}|{k}': v for (g, n, s), ts in res.items() for k, v in ts.items()})
    rows = []
    for (g, n, s), ts in sorted(res.items()):
        h = len(ts['t']) // 2
        row = dict(group=g, name=n, seed=s, alive=ts['alive'][h:].mean())
        for ch in ('eg', 'dyn'):
            row[f'{ch}_excess_pc_final'] = ts[f'{ch}_excess_pc'][-1]
            row[f'{ch}_excess_pc_slope'] = np.polyfit(ts['t'], ts[f'{ch}_excess_pc'], 1)[0]
            row[f'{ch}_tv_final'] = ts[f'{ch}_tv'][-1]
            row[f'{ch}_tv_mean'] = ts[f'{ch}_tv'][h:].mean()
            row[f'{ch}_D_final'] = ts[f'{ch}_D'][-1]
            G, Nn, pop = ts[f'{ch}_G_final'], ts[f'{ch}_N_final'], ts[f'{ch}_pop_final']
            row[f'{ch}_m_mean'] = float(ts[f'{ch}_m_mean'])
            row[f'{ch}_tv_vs_uniform'] = 0.5 * np.abs(G / G.sum() - 1 / len(G)).sum() if G.sum() else 0.0
            pres = pop > 0
            # the old excess recomputed from the real run alone. With m > 0 the shadow
            # holds mass on every key it has touched (D_N = keys with N > 0); with
            # m = 0 it equals the real composition and the C excess is exactly 0.
            row[f'{ch}_excess_real_only'] = (G[pres].sum() / pres.sum() - G.sum() / (Nn > 0).sum()
                                             if row[f'{ch}_m_mean'] > 0 and pres.any() else 0.0)
        for k in ('dyn_tv_input', 'dyn_tv_within'):
            row[f'{k}_final'] = ts[k][-1]
            row[f'{k}_mean'] = ts[k][h:].mean()
        rows.append(row)
    pd.DataFrame(rows).to_csv(os.path.join(HERE, 'runs.csv'), index=False)
    print('done', flush=True)


if __name__ == '__main__':
    main()
