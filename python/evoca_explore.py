"""
Headless EvoCA harness for parameter exploration.

run_sim(params, n_steps, ...) → dict of summary metrics.

No SDL, no widgets — just sim.step() in a tight loop with periodic
metric sampling. Designed for use from multiprocessing.Pool drivers
(see Scans/<date>_<name>/scan.py for examples).

Metrics are organised on three axes:

  Evolution
    n_distinct_genomes_mean      : mean over time of distinct G-buckets
    n_distinct_genomes_temporal_std
    unique_top_genomes           : count of distinct genomes ever in top-1
    excess_activity_final        : ΣG − ΣN at end of run (Channon)
    excess_activity_slope        : linear slope of ΣG − ΣN over time

  Spatial structure
    F_std_mean                   : mean over time of spatial std of F_food
    correlation_length_mean      : mean over time of v_curr correlation length
    largest_patch_mean           : mean over time of largest alive component
    n_patches_mean               : mean over time of alive patch count

  Temporal variation of spatial structure
    F_std_temporal_std
    correlation_length_temporal_std
    largest_patch_temporal_std

  Survival
    final_pop, min_pop, extinct
"""

import csv
import json
import os
import sys
from collections import defaultdict
from datetime import datetime

import numpy as np
from scipy import ndimage

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from evoca_py import EvoCA   # noqa: E402


# ── Fixed initialisation used by run_sim ─────────────────────────────
# These are the state() keyword arguments equivalent to what run_sim
# does manually. They get snapshotted into each scan's scan_config.json
# so any scan row can be reconstructed as a .evoca recipe later.

SCAN_INIT_DEFAULTS = {
    'lut': 'gol',
    'egenome': 'uniform',
    'egenome_value': 0b000011,   # = 3, plus-shape fiducial
    'v_density': 0.5,            # random binary with this density
    'f_init': 0.5,               # initial f_priv per alive cell
    'F': 'uniform',
    'F_init': 1.0,               # initial F_food everywhere
    'alive': 'halfplane',
    'alive_axis': 0,             # left half alive
}


# ── Spatial metric helpers ────────────────────────────────────────────

def _correlation_length(arr_2d):
    """First-zero-crossing radius of the 2D autocorrelation, averaged
    over angles. Returns float in [0, min(N)/2]."""
    arr = arr_2d.astype(np.float64) - arr_2d.mean()
    if arr.std() == 0:
        return 0.0
    h, w = arr.shape
    fft = np.fft.fft2(arr)
    acf = np.fft.ifft2(fft * np.conj(fft)).real
    norm = acf[0, 0]
    if norm <= 0:
        return 0.0
    acf = acf / norm
    max_r = min(h, w) // 2
    angles = np.linspace(0, 2 * np.pi, 16, endpoint=False)
    cos_a, sin_a = np.cos(angles), np.sin(angles)
    for r in range(1, max_r):
        ii = (np.round(r * cos_a).astype(int)) % h
        jj = (np.round(r * sin_a).astype(int)) % w
        if acf[ii, jj].mean() <= 0:
            return float(r)
    return float(max_r)


def _spatial_metrics(alive_2d, F_2d):
    out = {}
    out['F_std'] = float(F_2d.std())
    out['alive_density'] = float(alive_2d.mean())

    if alive_2d.any():
        labels, ncomp = ndimage.label(alive_2d)
        if ncomp > 0:
            sizes = ndimage.sum_labels(alive_2d, labels,
                                       np.arange(1, ncomp + 1))
            out['largest_patch'] = float(sizes.max())
            out['n_patches'] = int(ncomp)
        else:
            out['largest_patch'] = 0.0
            out['n_patches'] = 0
    else:
        out['largest_patch'] = 0.0
        out['n_patches'] = 0

    out['correlation_length'] = _correlation_length(alive_2d)
    return out


# ── Main entry point ──────────────────────────────────────────────────

_PARAM_KEYS = ('food_inc', 'm_scale', 'gdiff', 'mu_lut', 'mu_egenome',
               'tax', 'restricted_mu')


def run_sim(params, n_steps=5000, sample_every=100, N=256, seed=0,
            shadow=True, init='halfplane'):
    """
    Run one EvoCA configuration headless.

    params: dict with keys from _PARAM_KEYS. Defaults applied if missing.
    n_steps: total ticks.
    sample_every: tick interval at which we snapshot all metrics.
    N: grid side length.
    seed: numpy seed for v_curr initialisation.
    shadow: whether to enable Channon-style neutral shadow + N-activity.
    init: 'halfplane' (left half alive) or 'fraction' (random 50%).

    Returns a dict of summary metrics + parameter echo.
    """
    rng = np.random.default_rng(seed)

    sim = EvoCA()
    sim.init(N, **{k: params[k] for k in _PARAM_KEYS if k in params})
    sim.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    sim.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
    sim.set_f_all(0.5)
    sim.set_F_all(1.0)
    if init == 'halfplane':
        sim.set_alive_halfplane(0)
    elif init == 'fraction':
        sim.set_alive_fraction(0.5)
    else:
        raise ValueError(f"unknown init: {init}")

    if shadow:
        sim.neutral_enable()

    samples = defaultdict(list)
    dominator_history = set()

    for t in range(n_steps + 1):
        if t % sample_every == 0:
            alive_view = np.ctypeslib.as_array(
                sim._lib.evoca_get_alive(),
                shape=(N * N,)).reshape(N, N)
            F_view = np.ctypeslib.as_array(
                sim._lib.evoca_get_F(),
                shape=(N * N,)).reshape(N, N)

            spatial = _spatial_metrics(alive_view, F_view)
            for k, v in spatial.items():
                samples[k].append(v)

            sim._lib.evoca_activity_update()
            G = sim.get_activity(max_n=200000)
            gp = G['pop_count'] > 0
            n_g = int(gp.sum())
            samples['n_distinct_genomes'].append(n_g)
            samples['G_total_activity'].append(
                int(G['activity'][gp].sum()) if n_g else 0)
            if n_g:
                top_idx = int(np.argmax(G['pop_count'][gp]))
                top_hash = int(G['hash'][gp][top_idx])
                dominator_history.add(top_hash)

            if shadow:
                sim.n_activity_update()
                Nt = sim.get_n_activity(max_n=200000)
                npop = Nt['pop_count'] > 0
                samples['N_total_activity'].append(
                    int(Nt['activity'][npop].sum()) if npop.any() else 0)

            samples['t'].append(t)

        if t < n_steps:
            sim.step()

    # ── Compose summary ──────────────────────────────────────────────
    out = {}
    for k in _PARAM_KEYS:
        if k in params:
            out[f'param_{k}'] = params[k]
    out['N'] = N
    out['n_steps'] = n_steps
    out['seed'] = seed
    out['init'] = init
    out['shadow'] = shadow

    def _agg(name, key):
        if samples[key]:
            arr = np.asarray(samples[key], dtype=float)
            out[f'{name}_mean'] = float(arr.mean())
            out[f'{name}_temporal_std'] = float(arr.std())
            out[f'{name}_final'] = float(arr[-1])

    _agg('F_std', 'F_std')
    _agg('correlation_length', 'correlation_length')
    _agg('largest_patch', 'largest_patch')
    _agg('n_patches', 'n_patches')
    _agg('alive_density', 'alive_density')
    _agg('n_distinct_genomes', 'n_distinct_genomes')

    out['unique_top_genomes'] = len(dominator_history)

    if shadow and samples['G_total_activity'] and samples['N_total_activity']:
        G_arr = np.asarray(samples['G_total_activity'], dtype=float)
        N_arr = np.asarray(samples['N_total_activity'], dtype=float)
        excess = G_arr - N_arr
        out['excess_activity_final'] = float(excess[-1])
        if len(excess) >= 2:
            t_arr = np.asarray(samples['t'], dtype=float)
            slope, _ = np.polyfit(t_arr, excess, 1)
            out['excess_activity_slope'] = float(slope)
        else:
            out['excess_activity_slope'] = 0.0

    pop_arr = np.asarray(samples['alive_density'], dtype=float) * (N * N)
    if pop_arr.size:
        out['final_pop'] = int(pop_arr[-1])
        out['min_pop'] = int(pop_arr.min())
        out['extinct'] = bool(out['min_pop'] == 0)
    else:
        out['final_pop'] = 0
        out['min_pop'] = 0
        out['extinct'] = True

    sim.free()
    return out


# ── Scan dir helpers: scan_config.json + .evoca recipe export ────────

def save_scan_config(scan_dir, *, N, n_steps, sample_every, seed, shadow,
                     init_overrides=None):
    """Write scan_config.json into `scan_dir` documenting the fixed
    initialisation used by all configs in this scan. Any keys in
    init_overrides override SCAN_INIT_DEFAULTS for this scan."""
    init = dict(SCAN_INIT_DEFAULTS)
    if init_overrides:
        init.update(init_overrides)
    cfg = {
        'created': datetime.now().isoformat(timespec='seconds'),
        'N': int(N),
        'n_steps': int(n_steps),
        'sample_every': int(sample_every),
        'seed': int(seed),
        'shadow': bool(shadow),
        'initialization': init,
    }
    out_path = os.path.join(scan_dir, 'scan_config.json')
    with open(out_path, 'w') as f:
        json.dump(cfg, f, indent=2)
    return out_path


def _load_scan(scan_dir):
    cfg_path = os.path.join(scan_dir, 'scan_config.json')
    csv_path = os.path.join(scan_dir, 'results.csv')
    if not os.path.exists(cfg_path):
        raise FileNotFoundError(
            f"{cfg_path} not found — the scan dir needs a scan_config.json. "
            f"Re-run scan.py (now writes one), or call save_scan_config() "
            f"manually with the values that were used.")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"{csv_path} not found")
    with open(cfg_path) as f:
        cfg = json.load(f)
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    return cfg, rows


def _row_param_value(row, key):
    """Pull a param_<key> from a CSV row, parsing booleans."""
    raw = row.get(f'param_{key}', '')
    if raw in ('True', 'true'):  return True
    if raw in ('False', 'false'): return False
    try: return float(raw) if '.' in raw or 'e' in raw.lower() else int(raw)
    except (ValueError, TypeError): return raw


def evoca_from_scan(scan_dir, config_idx, descriptor=None,
                    runs_dir=None, probes=None, colormode=0):
    """Read scan_config.json + the row with given config_idx from
    results.csv; write a .evoca recipe to Runs/ (or runs_dir) that
    can be loaded with import_run(). Returns the recipe filepath.

    config_idx matches the 'config_idx' column in results.csv.
    """
    cfg, rows = _load_scan(scan_dir)

    # Find row by config_idx
    row = None
    for r in rows:
        try:
            if int(r.get('config_idx', -1)) == int(config_idx):
                row = r
                break
        except (ValueError, TypeError):
            continue
    if row is None:
        raise ValueError(f"config_idx={config_idx} not found in {scan_dir}/results.csv")

    metaparams = {}
    for k in ('food_inc', 'm_scale', 'gdiff', 'mu_lut', 'mu_egenome',
              'tax', 'restricted_mu'):
        v = _row_param_value(row, k)
        if v != '':
            metaparams[k] = v

    if runs_dir is None:
        runs_dir = os.path.normpath(
            os.path.join(_HERE, '..', 'Runs'))
    os.makedirs(runs_dir, exist_ok=True)

    if descriptor is None:
        scan_name = os.path.basename(os.path.normpath(scan_dir))
        descriptor = f'{scan_name}_cfg{config_idx}'
    safe = descriptor.replace(' ', '_').replace('/', '_')
    filename = f"{datetime.now().strftime('%Y-%m-%d')}_{safe}.evoca"
    filepath = os.path.join(runs_dir, filename)

    recipe = {
        'version': 3,
        'created': datetime.now().isoformat(timespec='seconds'),
        'descriptor': descriptor,
        'source_scan': os.path.basename(os.path.normpath(scan_dir)),
        'source_config_idx': int(config_idx),
        'N': int(cfg['N']),
        'metaparams_init': metaparams,
        'metaparams_final': metaparams,
        'initialization': dict(cfg['initialization']),
        'display': {
            'colormode': int(colormode),
            'probes': probes or {},
        },
    }
    with open(filepath, 'w') as f:
        json.dump(recipe, f, indent=2)
    return filepath


def evoca_from_scan_top(scan_dir, top_k=5, score_keys=None, descriptor_prefix=None,
                        runs_dir=None, probes=None, colormode=0):
    """Convenience: rank rows by a composite of normalised score_keys
    (default: correlation_length_mean + largest_patch_temporal_std +
    n_distinct_genomes_temporal_std + unique_top_genomes), and emit a
    .evoca for each of the top-K. Filters out extinct and >95% saturated
    configurations first.

    Returns list of (config_idx, filepath) tuples in score order."""
    cfg, rows = _load_scan(scan_dir)
    if score_keys is None:
        score_keys = [
            'correlation_length_mean',
            'largest_patch_temporal_std',
            'n_distinct_genomes_temporal_std',
            'unique_top_genomes',
        ]

    def _num(s, d=0.0):
        try: return float(s)
        except (TypeError, ValueError): return d

    candidates = [r for r in rows
                  if r.get('extinct', '').lower() == 'false'
                  and _num(r.get('alive_density_mean', 0)) <= 0.95]
    if not candidates:
        raise ValueError("no non-extinct, non-saturated rows in this scan")

    arrs = {k: np.array([_num(r.get(k, 0)) for r in candidates])
            for k in score_keys}
    def _norm(a):
        if a.max() == a.min(): return np.zeros_like(a)
        return (a - a.min()) / (a.max() - a.min())
    score = sum(_norm(arrs[k]) for k in score_keys) / len(score_keys)
    order = np.argsort(-score)[:top_k]

    out = []
    for rank, i in enumerate(order):
        r = candidates[i]
        idx = int(r['config_idx'])
        desc = (f"{descriptor_prefix or 'top'}_{rank+1}"
                f"_cfg{idx}_score{score[i]:.2f}")
        path = evoca_from_scan(scan_dir, idx, descriptor=desc,
                               runs_dir=runs_dir, probes=probes,
                               colormode=colormode)
        out.append((idx, path))
    return out
