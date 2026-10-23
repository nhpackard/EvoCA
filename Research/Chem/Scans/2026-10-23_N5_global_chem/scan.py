"""
CHEM-N5: the global-chem baseline before wgenes (CHEM-D3).

Two parts, both from the recipe 2026-10-23_chem_fw_spiralpoint (food/waste,
standard rule, mutation on, halfplane start), 6000 steps, statistics over
the last 3000:

  A  robustness of the balance regime (N=256):
     waste_sink 0.010..0.020 (step 0.001) x food_inc {0.013, 0.016, 0.020}
     x seeds {0, 1, 2}. Tests whether the three-state oscillating regime
     sits where waste_sink ~ tax x density, and whether it survives seeds
     and food supply (the first look was one seed at food_inc 0.016).
  B  first wave/spiral probe (N=512, where the legacy spirals fit):
     wdiff {0.01, 0.03, 0.08, 0.2, 0.5} x waste_sink {0.013, 0.014, 0.015,
     0.016} x seeds {0, 1, 2} at food_inc 0.016. Hypothesis: a slowly
     diffusing inhibitor (wdiff < gdiff = 0.08) leaves a sharp refractory
     tail behind a front (waves); a fast one spreads ahead of it (lateral
     inhibition, stationary patches). Seed-0 runs save quantized F/W/alive
     frames (frames/, not committed) for inspection by eye.

Wave/spiral evidence is a phase-singularity count: each b x b block gets a
cycle phase theta = atan2(W~, F~) from its standardized mean food and
waste (the excitable cycle food -> colony -> waste -> food is a loop in the
(F, W) plane); a plaquette whose phase winds by +-2 pi is a singularity, a
spiral core if it persists between frames PERSIST_LAG steps apart while
the phase turns around it. Block phases use
absolute (not standardized) deviations, and plaquettes with low local cycle
amplitude are ignored, so a uniform texture's noise does not count.
Persistence (same sign, within 1 coarse cell, PERSIST_LAG steps later) is reported
with a chance level (the earlier frame shifted by half the lattice), and
counts are also read against the no-structure runs above the balance.
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

from evoca_explore import build_sim, save_scan_config, _correlation_length   # noqa: E402

RECIPE = os.path.join(ROOT, 'Research', 'Chem', 'Runs',
                      '2026-10-23_chem_fw_spiralpoint.evoca')
N_WORKERS = 28
N_STEPS = 6000
STATS_FROM = 3000          # statistics over steps (STATS_FROM, N_STEPS]
SAMPLE_EVERY = 10          # alive / W / F means
SPATIAL_EVERY = 250        # correlation lengths and singularity counts
BLOCK = 16                 # coarse-graining block for the phase field (patch scale ~50)
AMP_MIN = 0.1              # minimum local (F, W) cycle amplitude for a singularity
PERSIST_LAG = 60           # steps between the frames compared (~half the ~120-step period)
ROT_MIN = 0.5              # rad over PERSIST_LAG: phase advance around a rotating core
FRAME_STEPS = (4000, 5000, 6000)

SINKS_A = [round(0.010 + 0.001 * i, 3) for i in range(11)]
FOOD_INC_A = [0.013, 0.016, 0.020]
WDIFF_B = [0.01, 0.03, 0.08, 0.2, 0.5]
SINKS_B = [0.013, 0.014, 0.015, 0.016]
SEEDS = [0, 1, 2]


def phase_singularities(F, W, b=BLOCK, amp_min=AMP_MIN):
    """Positions (k x 2 int array, coarse grid) and signs of the phase
    singularities of theta = atan2(W~, F~) on the b x b coarse grid, where
    F~, W~ are block means minus their lattice means (absolute units, not
    standardized, so a uniform texture's small fluctuations do not count).
    Plaquettes whose 3 x 3 coarse neighborhood has mean amplitude
    hypot(F~, W~) below amp_min are ignored."""
    n = F.shape[0] // b
    Fc = F[:n * b, :n * b].reshape(n, b, n, b).mean(axis=(1, 3))
    Wc = W[:n * b, :n * b].reshape(n, b, n, b).mean(axis=(1, 3))
    x = Fc - Fc.mean()
    y = Wc - Wc.mean()
    th = np.arctan2(y, x)
    t01 = np.roll(th, -1, axis=1)
    t11 = np.roll(t01, -1, axis=0)
    t10 = np.roll(th, -1, axis=0)

    def wrap(d):
        return (d + np.pi) % (2 * np.pi) - np.pi

    wind = wrap(t01 - th) + wrap(t11 - t01) + wrap(t10 - t11) + wrap(th - t10)
    q = np.rint(wind / (2 * np.pi)).astype(int)
    r = np.hypot(x, y)
    r_loc = sum(np.roll(np.roll(r, i, 0), j, 1)
                for i in (-1, 0, 1) for j in (-1, 0, 1)) / 9.0
    q[r_loc < amp_min] = 0
    pos = np.argwhere(q != 0)
    return pos, q[q != 0], n, th


def rotating(pos_a, q_a, th_a, pos_b, q_b, th_b, n, tol=1, dth_min=ROT_MIN):
    """Persistent singularities (same sign within tol between frame b and
    the later frame a) around which the phase advanced: the mean wrapped
    phase change over the 3 x 3 coarse neighborhood exceeds dth_min in
    magnitude. A spiral core stays put while the phase turns around it; a
    defect of a stationary patch mosaic stays put with the phase frozen."""
    if len(pos_a) == 0 or len(pos_b) == 0:
        return 0
    d = np.abs(pos_a[:, None, :] - pos_b[None, :, :])
    d = np.minimum(d, n - d).max(axis=2)
    keep = ((d <= tol) & (q_a[:, None] == q_b[None, :])).any(axis=1)
    dth = (th_a - th_b + np.pi) % (2 * np.pi) - np.pi
    count = 0
    for (i, j) in pos_a[keep]:
        nb = dth[np.ix_([(i + k) % n for k in (-1, 0, 1, 2)],
                        [(j + k) % n for k in (-1, 0, 1, 2)])]
        if abs(nb.mean()) > dth_min:
            count += 1
    return count


def persistent(pos_a, q_a, pos_b, q_b, n, tol=1):
    """Number of singularities in a with a same-sign one in b within tol
    coarse cells (periodic Chebyshev distance)."""
    if len(pos_a) == 0 or len(pos_b) == 0:
        return 0
    d = np.abs(pos_a[:, None, :] - pos_b[None, :, :])
    d = np.minimum(d, n - d).max(axis=2)
    same = q_a[:, None] == q_b[None, :]
    return int(((d <= tol) & same).any(axis=1).sum())


def period_of(series, dt):
    """Dominant period (steps): the first local maximum of the
    autocorrelation after its first zero crossing; nan if that peak is
    below 0.2."""
    x = np.asarray(series, float) - np.mean(series)
    if len(x) < 4 or x.std() == 0:
        return float('nan'), 0.0
    acf = np.correlate(x, x, 'full')[len(x) - 1:]
    acf /= acf[0]
    half = len(acf) // 2
    neg = np.where(acf[:half] < 0)[0]
    if len(neg) == 0:
        return float('nan'), 0.0
    for k in range(neg[0] + 1, half - 1):
        if acf[k] >= acf[k - 1] and acf[k] >= acf[k + 1]:
            peak = float(acf[k])
            return (float(k * dt) if peak >= 0.2 else float('nan')), peak
    return float('nan'), 0.0


def quantize(a, hi):
    return np.clip(np.asarray(a) / hi * 255, 0, 255).astype(np.uint8)


def run_one(task):
    part, N, seed, over = task
    rec = json.load(open(RECIPE))
    params = dict(rec['metaparams_init'])
    params.update(over)
    row = {'part': part, 'N': N, 'seed': seed,
           **{f'param_{k}': v for k, v in params.items()}}
    try:
        sim = build_sim(params, N=N, seed=seed, init='halfplane')
        alive_s, W_s, F_s = [], [], []
        cl_alive, cl_W, n_sing, n_pers, n_null, n_rot, tox = [], [], [], [], [], [], 0
        flux0 = None
        frames = {}
        prev = None
        for t in range(1, N_STEPS + 1):
            sim.step()
            if t > STATS_FROM:
                tox += sim.get_tox_deaths_last()
            if t == STATS_FROM:
                flux0 = sim.get_chem_flux(1)
            if t % SAMPLE_EVERY == 0:
                a = sim.get_alive()
                alive_s.append(a.mean())
                W = sim.get_W()
                F = sim.get_F()
                W_s.append(W.mean())
                F_s.append(F.mean())
                if a.mean() == 0:
                    break
            if t > STATS_FROM and (t % SPATIAL_EVERY in (0, SPATIAL_EVERY - PERSIST_LAG)):
                F = sim.get_F()
                W = sim.get_W()
                pos, sgn, n, th = phase_singularities(F, W)
                if t % SPATIAL_EVERY == 0:
                    a = sim.get_alive()
                    cl_alive.append(_correlation_length(a))
                    cl_W.append(_correlation_length(W))
                    n_sing.append(len(pos))
                    if prev is not None:
                        pp, pq, pth = prev
                        n_pers.append(persistent(pos, sgn, pp, pq, n))
                        n_rot.append(rotating(pos, sgn, th, pp, pq, pth, n))
                        n_null.append(persistent(pos, sgn, (pp + n // 2) % n, pq, n))
                else:
                    prev = (pos, sgn, th)
            if part == 'B' and seed == 0 and t in FRAME_STEPS:
                frames[f'F_{t}'] = quantize(sim.get_F(), 1.0)
                frames[f'W_{t}'] = quantize(sim.get_W(), 2.0)
                frames[f'alive_{t}'] = sim.get_alive().astype(np.uint8)
        flux1 = sim.get_chem_flux(1)
        sim.free()
        k = len(alive_s) - (N_STEPS - STATS_FROM) // SAMPLE_EVERY
        tail = slice(max(k, 0), None)
        a_t = np.asarray(alive_s[tail])
        per, peak = period_of(a_t, SAMPLE_EVERY)
        steps = N_STEPS - STATS_FROM
        row.update({
            'extinct': int(alive_s[-1] == 0),
            'alive_mean': float(a_t.mean()) if len(a_t) else 0.0,
            'alive_std': float(a_t.std()) if len(a_t) else 0.0,
            'W_mean': float(np.mean(W_s[tail])) if len(a_t) else float('nan'),
            'W_std': float(np.std(W_s[tail])) if len(a_t) else float('nan'),
            'F_mean': float(np.mean(F_s[tail])) if len(a_t) else float('nan'),
            'period': per, 'acf_peak': peak,
            'corr_len_alive': float(np.mean(cl_alive)) if cl_alive else float('nan'),
            'corr_len_W': float(np.mean(cl_W)) if cl_W else float('nan'),
            'n_sing': float(np.mean(n_sing)) if n_sing else float('nan'),
            'n_sing_persist': float(np.mean(n_pers)) if n_pers else float('nan'),
            # chance level: the earlier frame shifted by half the lattice
            'n_sing_persist_null': float(np.mean(n_null)) if n_null else float('nan'),
            # persistent cores with the phase turning around them: spiral cores
            'n_sing_rotating': float(np.mean(n_rot)) if n_rot else float('nan'),
            # waste production and removal per site per step over the window
            'w_in_rate': (flux1[0] - flux0[0]) / steps / N**2 if flux0 else float('nan'),
            'w_out_rate': (flux1[1] - flux0[1]) / steps / N**2 if flux0 else float('nan'),
            'tox_deaths_per_step': tox / steps,
            'error': '',
        })
        if frames:
            os.makedirs(os.path.join(HERE, 'frames'), exist_ok=True)
            tag = f"B_wd{over['wdiff']}_ws{over['waste_sink']}_s{seed}"
            np.savez_compressed(os.path.join(HERE, 'frames', tag + '.npz'), **frames)
    except Exception as e:
        row['error'] = repr(e)
    return row


def tasks():
    out = []
    for fi in FOOD_INC_A:
        for ws in SINKS_A:
            for s in SEEDS:
                out.append(('A', 256, s, {'food_inc': fi, 'waste_sink': ws}))
    for wd in WDIFF_B:
        for ws in SINKS_B:
            for s in SEEDS:
                out.append(('B', 512, s, {'wdiff': wd, 'waste_sink': ws}))
    # longest first, so the pool's tail is short
    out.sort(key=lambda t: -t[1])
    return out


def main():
    cfg_path = save_scan_config(HERE, N=512, n_steps=N_STEPS,
                                sample_every=SAMPLE_EVERY, seed=0, shadow=False)
    # N and seed vary per row (columns N, seed); record that, not one value
    cfg = json.load(open(cfg_path))
    cfg.update({'N': 'per row: part A 256, part B 512', 'seed': 'per row: 0, 1, 2',
                'recipe': os.path.relpath(RECIPE, HERE)})
    json.dump(cfg, open(cfg_path, 'w'), indent=2)
    ts = tasks()
    print(f"Running {len(ts)} runs across {N_WORKERS} workers...", flush=True)
    t0 = time.perf_counter()
    results = []
    with Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(run_one, ts):
            results.append(r)
            if len(results) % 20 == 0 or len(results) == len(ts):
                print(f"  {len(results)}/{len(ts)}  {time.perf_counter() - t0:5.0f}s", flush=True)
    results.sort(key=lambda r: (r['part'], r.get('param_food_inc', 0),
                                r.get('param_wdiff', 0), r.get('param_waste_sink', 0),
                                r['seed']))
    keys = sorted({k for r in results for k in r})
    with open(os.path.join(HERE, 'results.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(results)
    print(f"Wrote {len(results)} rows; errors: {sum(bool(r['error']) for r in results)}")


if __name__ == '__main__':
    main()
