"""Drift-shadow validation campaign (Docs/shadows.md §7).

Usage (from the repo root, after `make`):
  python3 tools/drift_validation.py a       [n_runs=400]   # synthetic neutral
  python3 tools/drift_validation.py b       [n_runs=100]   # synthetic power
  python3 tools/drift_validation.py c                      # EvoCA runs
  python3 tools/drift_validation.py timing
  python3 tools/drift_validation.py all

Results go to tests/notebooks/drift-shadow-validation.json (one key per
part, merged into the existing file), which the pre-merge notebook
tests/notebooks/drift-shadow-test.ipynb reads and plots.

(a) Synthetic neutral: "real" = one extra replicate of the shadow process
    driven by the same event stream (python/drift_synth.py). Cells
    theta in {0.5, 2, 10, 50} x accessibility {uniform, geo0.7, geo0.9};
    each run carries dyn-like (250 loci, weights = accessibility) AND
    eg-like (729-key) and marker variation, scored separately. R = 19
    (the rank p-value at R = 8 cannot reach 0.05); R = 8 statistics are
    recomputed from the first 8 shadows of the same runs. Plus the
    demographic variant (real parents without replacement).
(b) Synthetic power at Np = 500, theta = 10: purifying selection on the
    top input (s = 0.05, k = 1) and on 5 top inputs (s = 0.02); s and
    generation sweeps for monotonicity; directional selection (sweeps);
    eg family selection.
(c) EvoCA runs (run_sim): frozen-clonal marker null, draft floor (JOINT,
    LUT_only), standing-variation freeze, pure-evo joint winner and frozen
    control recipes, dyn_expression_bias.
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'python'))
OUT = os.path.join(ROOT, 'tests', 'notebooks', 'drift-shadow-validation.json')
NPROC = max(1, (os.cpu_count() or 2) - 2)

NP, B = 500, 25                     # synthetic population, births = deaths / step
THETAS = (0.5, 2.0, 10.0, 50.0)
ACCESS = ('uniform', 'geo0.7', 'geo0.9')
STATS_KEEP = ('dyn_drift_consistency', 'dyn_drift_z', 'dyn_drift_p',
              'dyn_drift_concentration', 'dyn_drift_concentration_z',
              'dyn_drift_concentration_p', 'dyn_drift_tv_z', 'dyn_drift_tv_p',
              'dyn_drift_frac_z3', 'dyn_drift_nhi',
              'eg_drift_concentration', 'eg_drift_z', 'eg_drift_p',
              'eg_drift_tv_z', 'eg_drift_tv_p', 'eg_drift_frac_z3', 'eg_drift_nhi',
              'marker_drift_concentration', 'marker_drift_z', 'marker_drift_p',
              'marker_drift_tv_z', 'marker_drift_tv_p', 'marker_drift_frac_z3',
              'dyn_drift_consistency_runmean', 'dyn_drift_z_runmean')


def _sub(acc, R):
    """The same run scored with only the first R shadows."""
    a = dict(acc)
    for k in ('eg_shadow', 'dyn_shadow', 'marker_shadow', 'psum_shadow', 'negene_shadow'):
        if a.get(k) is not None:
            a[k] = np.asarray(a[k])[:R]
    a['R'] = R
    return a


def _scores(acc):
    from drift_stats import drift_statistics
    st = drift_statistics(acc)
    return {k: float(st[k]) for k in STATS_KEEP if k in st}


def _synth_job(job):
    from drift_synth import synth_run
    seed, kw, Rs = job
    acc = synth_run(seed, **kw)
    return {R: _scores(_sub(acc, R)) for R in Rs}


def _pool_map(fn, jobs):
    with Pool(NPROC) as p:
        return p.map(fn, jobs, chunksize=4)


def _summ(vals, R):
    """Calibration summary of one statistic over runs (p and z lists)."""
    p = np.asarray([v[0] for v in vals], float)
    z = np.asarray([v[1] for v in vals], float)
    zf = z[np.isfinite(z)]
    # Chi-square goodness of fit of the discrete rank p-values to uniform
    # on {1/(R+1), ..., 1} (KS is not valid for a discrete law).
    from scipy.stats import chisquare
    grid = np.arange(1, R + 2) / (R + 1)
    idx = np.clip(np.rint(p * (R + 1)).astype(int) - 1, 0, R)
    counts = np.bincount(idx, minlength=R + 1)
    unif_p = float(chisquare(counts).pvalue)
    return {'n': int(p.size), 'reject_005': float(np.mean(p <= 0.05 + 1e-12)),
            'reject_floor': float(np.mean(p <= 1.0 / (R + 1) + 1e-12)),
            'z_mean': float(zf.mean()) if zf.size else float('nan'),
            'z_sd': float(zf.std(ddof=1)) if zf.size > 1 else float('nan'),
            'n_inf_z': int((~np.isfinite(z)).sum()),
            'unif_chi2_p': unif_p, 'p_hist': counts.tolist(),
            'grid': grid.tolist()}


def part_a(n_runs=400):
    t0 = time.time()
    res = {'cells': [], 'variant': None, 'settings': dict(
        Np=NP, b=B, steps=1500, generations=75, R_primary=19, R_default=8,
        n_runs=n_runs)}
    jobs, keys = [], []
    for th in THETAS:
        mu = th / (2 * NP)
        for acc in ACCESS:
            kw = dict(R=19, Np=NP, b=B, steps=1500, mu_lut=mu, mu_eg=mu,
                      mu_act=mu / 4, mu_mk=mu, w=acc)
            for i in range(n_runs):
                jobs.append((10_000 + 1000 * len(keys) + i, kw, (19, 8)))
            keys.append((th, acc))
    # demographic variant: real = without-replacement parents
    kw_v = dict(R=19, Np=NP, b=B, steps=1500, mu_lut=0.01, mu_eg=0.01,
                mu_act=0.0025, mu_mk=0.01, w='geo0.7', sel_mode=3)
    for i in range(n_runs):
        jobs.append((90_000 + i, kw_v, (19, 8)))
    out = _pool_map(_synth_job, jobs)
    for ci, (th, acc) in enumerate(keys):
        runs = out[ci * n_runs:(ci + 1) * n_runs]
        res['cells'].append({'theta': th, 'access': acc,
                             **_cell_summary(runs)})
    res['variant'] = {'theta': 10.0, 'access': 'geo0.7', 'real_parents': 'without replacement',
                      **_cell_summary(out[len(keys) * n_runs:])}
    res['seconds'] = time.time() - t0
    return res


def _cell_summary(runs):
    cs = {}
    for R in (19, 8):
        rr = [r[R] for r in runs]
        cs[f'R{R}'] = {
            'dyn_consistency': _summ([(r['dyn_drift_p'], r['dyn_drift_z']) for r in rr], R),
            'dyn_concentration': _summ([(r['dyn_drift_concentration_p'],
                                         r['dyn_drift_concentration_z']) for r in rr], R),
            'dyn_tv': _summ([(r['dyn_drift_tv_p'], r['dyn_drift_tv_z']) for r in rr], R),
            'eg_concentration': _summ([(r['eg_drift_p'], r['eg_drift_z']) for r in rr], R),
            'eg_tv': _summ([(r['eg_drift_tv_p'], r['eg_drift_tv_z']) for r in rr], R),
            'marker_concentration': _summ([(r['marker_drift_p'], r['marker_drift_z'])
                                           for r in rr], R),
            'dyn_frac_z3': float(np.mean([r['dyn_drift_frac_z3'] for r in rr])),
            'eg_frac_z3': float(np.mean([r['eg_drift_frac_z3'] for r in rr])),
            'marker_frac_z3': float(np.mean([r['marker_drift_frac_z3'] for r in rr])),
            'dyn_nhi_mean': float(np.mean([r['dyn_drift_nhi'] for r in rr])),
            'eg_nhi_mean': float(np.mean([r['eg_drift_nhi'] for r in rr])),
        }
    return cs


def part_b(n_runs=100):
    t0 = time.time()
    mu = 0.01                     # theta = 2 * 500 * 0.01 = 10
    base = dict(R=19, Np=NP, b=B, mu_lut=mu, mu_eg=mu, mu_act=mu / 4, mu_mk=mu,
                w='geo0.7')
    cells = [
        # (name, kwargs)
        ('purify s=0.05 k=1 gen150', dict(steps=3000, sel_mode=1, s=0.05, loci=[0])),
        ('purify s=0.02 k=5 gen150', dict(steps=3000, sel_mode=1, s=0.02, loci=[0, 1, 2, 3, 4])),
        ('purify s=0.01 k=1 gen75', dict(steps=1500, sel_mode=1, s=0.01, loci=[0])),
        ('purify s=0.02 k=1 gen75', dict(steps=1500, sel_mode=1, s=0.02, loci=[0])),
        ('purify s=0.05 k=1 gen75', dict(steps=1500, sel_mode=1, s=0.05, loci=[0])),
        ('purify s=0.01 k=1 gen150', dict(steps=3000, sel_mode=1, s=0.01, loci=[0])),
        ('purify s=0.02 k=1 gen150', dict(steps=3000, sel_mode=1, s=0.02, loci=[0])),
        ('directional s=0.05 k=1 gen150', dict(steps=3000, sel_mode=1, s=-0.05, loci=[0])),
        ('directional s=0.02 k=5 gen150', dict(steps=3000, sel_mode=1, s=-0.02,
                                               loci=[0, 1, 2, 3, 4])),
        ('eg family s=0.05 gen150', dict(steps=3000, sel_mode=2, s=0.05)),
    ]
    jobs = []
    for ci, (_, kw) in enumerate(cells):
        for i in range(n_runs):
            jobs.append((200_000 + 1000 * ci + i, dict(base, **kw), (19, 8)))
    out = _pool_map(_synth_job, jobs)
    res = {'cells': [], 'settings': dict(Np=NP, b=B, theta=10, weights='geo0.7',
                                         n_runs=n_runs)}
    for ci, (name, kw) in enumerate(cells):
        runs = out[ci * n_runs:(ci + 1) * n_runs]
        cell = {'name': name}
        for R in (19, 8):
            alpha = 0.05 if R == 19 else 1.0 / 9 + 1e-12
            rr = [r[R] for r in runs]
            cell[f'R{R}'] = {
                'alpha': alpha,
                'power_dyn_consistency': float(np.mean([r['dyn_drift_p'] <= alpha for r in rr])),
                'power_dyn_concentration': float(np.mean([r['dyn_drift_concentration_p'] <= alpha
                                                          for r in rr])),
                'power_dyn_tv': float(np.mean([r['dyn_drift_tv_p'] <= alpha for r in rr])),
                'power_eg_concentration': float(np.mean([r['eg_drift_p'] <= alpha for r in rr])),
                'power_eg_tv': float(np.mean([r['eg_drift_tv_p'] <= alpha for r in rr])),
                'power_marker_concentration': float(np.mean([r['marker_drift_p'] <= alpha
                                                             for r in rr])),
                'mean_dyn_z': float(np.nanmean([r['dyn_drift_z'] for r in rr])),
                'mean_dyn_consistency': float(np.nanmean([r['dyn_drift_consistency'] for r in rr])),
                'mean_dyn_tv_z': float(np.nanmean(np.clip([r['dyn_drift_tv_z'] for r in rr],
                                                          -1e3, 1e3))),
                'mean_eg_z': float(np.nanmean([r['eg_drift_z'] for r in rr])),
            }
        res['cells'].append(cell)
    res['seconds'] = time.time() - t0
    return res


# ── (c) EvoCA ──────────────────────────────────────────────────────────

WINNER = dict(food_inc=0.018, m_scale=2.5, gdiff=0.08, mu_lut=0.01, mu_egene=0.003,
              mu_egenome=0.003, p_dup_egene=1.0, tax=0.035, tax_per_egene=0.0008,
              tax_lut=0.0, tax_ring=0.0, restricted_mu=0)
FROZEN = dict(WINNER, mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0, tax_per_egene=0.0)
CC_BASE = dict(food_inc=0.013, m_scale=1.2, gdiff=0.06, tax=0.035, tax_per_egene=8e-4)
JOINT = dict(CC_BASE, mu_lut=0.01, mu_egene=0.003, mu_egenome=0.005)
LUT_ONLY = dict(CC_BASE, mu_lut=0.01, mu_egene=0.0, mu_egenome=0.0)

C_KEEP = STATS_KEEP + ('eg_drift_tv', 'dyn_drift_tv', 'marker_drift_tv', 'negene_drift',
                       'dyn_expression_bias', 'dyn_expression_bias_signed',
                       'eg_drift_z_last', 'dyn_drift_z_last', 'dyn_drift_tv_z_last',
                       'eg_drift_tv_z_last', 'marker_drift_tv_z_last',
                       'eg_drift_concentration_last', 'dyn_drift_consistency_last',
                       'marker_drift_concentration_last', 'eg_drift_tv_last',
                       'dyn_drift_tv_last', 'marker_drift_tv_last',
                       'marker_drift_tv_p', 'final_pop', 'min_pop', 'extinct',
                       'eg_excess_pc_final', 'dyn_excess_pc_final')


def _recipe_params(name):
    path = os.path.join(ROOT, 'Research', 'Activity-metrics', 'Runs', name)
    with open(path) as f:
        return dict(json.load(f)['metaparams_init'])


def _c_job(job):
    """One EvoCA run. job = (label, params, N, seed, n_steps, reps, marker_mu,
    lut_init) with lut_init None (GoL founder) or 'random'."""
    from evoca_explore import run_sim
    label, params, N, seed, n_steps, reps, marker_mu, lut_init = job
    if lut_init == 'random':
        out = _run_random_founder(params, N, seed, n_steps, reps, marker_mu)
    else:
        out = run_sim(params, n_steps=n_steps, sample_every=max(n_steps // 10, 1),
                      N=N, seed=seed, shadow=False, fs_shadow='both', fs_reps=reps,
                      marker_mu=marker_mu)
    row = {k: (float(out[k]) if k in out else float('nan')) for k in C_KEEP}
    row.update(label=label, N=N, seed=seed, marker_mu=marker_mu, reps=reps)
    return row


def _run_random_founder(params, N, seed, n_steps, reps, marker_mu):
    """Standing-variation freeze: each cell a random LUT (n_init=3), all
    mutation off. Mirrors run_sim's drift block (statistics from t=0)."""
    from evoca_explore import build_sim
    from drift_stats import drift_statistics, window
    sim = build_sim(params, N=N, seed=seed)
    sim.set_lut_random(n_init=3, density=0.5, seed=seed)
    sim.set_marker_mu(marker_mu)
    sim.fs_drift_enable(reps)
    snap = None
    for t in range(n_steps):
        if t == n_steps // 2:
            snap = sim.fs_drift_get()
        sim.step()
    acc = sim.fs_drift_get()
    out = drift_statistics(acc)
    for k, v in drift_statistics(window(acc, snap)).items():
        out[k + '_last'] = v
    pop = int(sim.get_alive().sum())
    out.update(final_pop=pop, min_pop=pop, extinct=pop == 0)
    sim.free()
    return out


def part_c(seeds_null=20, steps=5000):
    t0 = time.time()
    jobs = []
    for N in (64, 128):
        for mm in (0.001, 0.003, 0.01):
            for sd in range(seeds_null):
                jobs.append(('frozen_clonal_marker', FROZEN, N, 100 + sd, steps, 19, mm, None))
    for name, P in (('draft_JOINT', JOINT), ('draft_LUT_only', LUT_ONLY)):
        for sd in range(4):
            jobs.append((name, P, 128, 300 + sd, steps, 19, 0.003, None))
    for sd in range(20):
        jobs.append(('standing_variation_freeze', FROZEN, 64, 500 + sd, 3000, 19, 0.003,
                     'random'))
    win = _recipe_params('2026-10-15_pure_evo_joint_winner.evoca')
    frz = _recipe_params('2026-10-15_pure_evo_frozen_control.evoca')
    for sd in range(4):
        jobs.append(('recipe_joint_winner', win, 128, 700 + sd, steps, 19, 0.003, None))
        jobs.append(('recipe_frozen_control', frz, 128, 700 + sd, steps, 19, 0.003, None))
    rows = _pool_map(_c_job, jobs)
    res = {'rows': rows, 'settings': dict(steps=steps, seeds_null=seeds_null, reps=19)}
    # Marker null summary
    summ = {}
    for N in (64, 128):
        for mm in (0.001, 0.003, 0.01):
            rr = [r for r in rows if r['label'] == 'frozen_clonal_marker'
                  and r['N'] == N and r['marker_mu'] == mm and not r['extinct']]
            p = np.array([r['marker_drift_p'] for r in rr])
            ptv = np.array([r['marker_drift_tv_p'] for r in rr])
            summ[f'N{N}_mu{mm}'] = {
                'n': len(rr),
                'reject_conc_005': float(np.mean(p <= 0.05)) if len(rr) else float('nan'),
                'reject_tv_005': float(np.mean(ptv <= 0.05)) if len(rr) else float('nan'),
                'conc_q95': float(np.quantile([r['marker_drift_concentration'] for r in rr], 0.95))
                if rr else float('nan'),
                'abs_conc_q95': float(np.quantile([abs(r['marker_drift_concentration'])
                                                   for r in rr], 0.95)) if rr else float('nan'),
                'tv_q95': float(np.quantile([r['marker_drift_tv'] for r in rr], 0.95))
                if rr else float('nan'),
                'z_mean': float(np.mean([r['marker_drift_z'] for r in rr])) if rr else float('nan'),
                'z_sd': float(np.std([r['marker_drift_z'] for r in rr], ddof=1))
                if len(rr) > 1 else float('nan'),
            }
    res['marker_null'] = summ
    sv = [r for r in rows if r['label'] == 'standing_variation_freeze' and not r['extinct']]
    res['standing_variation'] = {
        'n': len(sv),
        'n_dyn_z_pos': int(sum(r['dyn_drift_z'] > 0 for r in sv)),
        'n_dyn_consistency_pos': int(sum(r['dyn_drift_consistency'] > 0 for r in sv)),
        'n_dyn_tv_p_floor': int(sum(r['dyn_drift_tv_p'] <= 0.05 for r in sv)),
    }
    res['seconds'] = time.time() - t0
    return res


def timing():
    from evoca_explore import build_sim
    P = JOINT
    res = {}
    for label, R in (('off', 0), ('drift_R1', 1), ('drift_R8', 8), ('drift_R16', 16)):
        sim = build_sim(P, N=256, seed=1)
        for _ in range(1500):             # past the founder transient
            sim.step()
        if R:
            sim.fs_drift_enable(R)
        for _ in range(50):
            sim.step()
        t = time.perf_counter()
        n = 400
        for _ in range(n):
            sim.step()
        ms = (time.perf_counter() - t) / n * 1e3
        res[label] = {'ms_per_step': ms, 'alive': int(sim.get_alive().sum())}
        sim.free()
    base = res['off']['ms_per_step']
    for k in res:
        res[k]['overhead_pct'] = 100.0 * (res[k]['ms_per_step'] / base - 1.0)
    # run_sim wall time per mode (includes per-sample metric work)
    from evoca_explore import run_sim
    for mode in ('tracking', 'drift', 'both'):
        t = time.perf_counter()
        run_sim(P, n_steps=600, sample_every=100, N=256, seed=1, shadow=False,
                fs_shadow=mode)
        res[f'run_sim_{mode}_ms_per_step'] = (time.perf_counter() - t) / 600 * 1e3
    return res


def _save(part, data):
    allres = {}
    if os.path.exists(OUT):
        with open(OUT) as f:
            allres = json.load(f)
    allres[part] = data
    with open(OUT, 'w') as f:
        json.dump(allres, f, indent=1)


if __name__ == '__main__':
    what = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else None
    if what in ('a', 'all'):
        _save('a', part_a(n or 400))
        print('a done')
    if what in ('b', 'all'):
        _save('b', part_b(n or 100))
        print('b done')
    if what in ('c', 'all'):
        _save('c', part_c())
        print('c done')
    if what in ('timing', 'all'):
        _save('timing', timing())
        print(json.dumps(json.load(open(OUT))['timing'], indent=1))
