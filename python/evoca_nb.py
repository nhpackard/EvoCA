"""Small helpers shared by the lab notebooks Research/<WS>/<WS>.ipynb.

Construction of the standard headless initialization lives in
evoca_explore.build_sim (used directly by the notebooks); this module
holds the notebook-side glue around recipes, scans and frames:

    rgb, argb_to_rgb         lattice frame -> (N, N, 3) uint8
    make_sim                 bare EvoCA or ChemCA (params['model']): init + set_seed
    load_recipe              import_run + reproducible seeding
    sim_from_scan            scan row -> recipe (temp dir) -> load_recipe
    params_from_row          scan results row -> build_sim params
    export_recipe_once       export a recipe unless that descriptor exists
    slope                    least-squares slope
"""
import inspect
import os
import re
import shutil
import tempfile
import warnings
import json

import numpy as np

import evoca_paths
import sim_model
from evoca_py import EvoCA, import_run


# ── frames ────────────────────────────────────────────────────────────

def argb_to_rgb(a):
    """ARGB int32/uint32 array (..., ) -> (..., 3) uint8 RGB."""
    a = np.asarray(a).astype(np.uint32)
    return np.stack([(a >> 16) & 255, (a >> 8) & 255, a & 255], -1).astype(np.uint8)


def rgb(sim, colormode=0):
    """Headless lattice frame for a colormode, as (N, N, 3) uint8."""
    N = sim._N
    return argb_to_rgb(sim.colorize(np.zeros(N * N, np.int32), colormode).reshape(N, N))


# ── construction ──────────────────────────────────────────────────────

def make_sim(N, seed=0, **params):
    """EvoCA() + init(N, **params) + set_seed(seed); no state() call.
    params['model'] == 'chemca' (sim_model) makes a ChemCA instead."""
    model, params = sim_model.split_model(params)
    if model == 'chemca':
        from chemca_py import ChemCA
        s = ChemCA()
    else:
        s = EvoCA()
    s.init(N, **params)
    s.set_seed(seed)
    return s


def _check_overrides(overrides):
    bad = [k for k in overrides if not hasattr(EvoCA, f'update_{k}')]
    if bad:
        raise TypeError(f"no EvoCA.update_<key>() setter for overrides {bad}")


def load_recipe(name, N=None, seed=0, recipe='init', seed_v=True,
                shadow=False, **overrides):
    """import_run(name, N, recipe), then make it reproducible: set_seed(seed)
    and, if seed_v, a v field from np.random.default_rng(seed) (import_run
    draws v from an unseeded generator). Overrides go through the
    update_<key>() setters; shadow=True calls neutral_enable() last.
    Returns (sim, display_kwargs)."""
    sim, kw = import_run(name, N=N, recipe=recipe)
    if not isinstance(sim, EvoCA):
        # ChemCA: state() already applied the recipe's own seed, so the run
        # is reproducible as loaded; overrides go through update().
        if shadow:
            raise ValueError("shadow is an EvoCA probe")
        if overrides:
            sim.update(**overrides)
        return sim, kw
    _check_overrides(overrides)
    sim.set_seed(seed)
    if seed_v:
        n = sim._N
        sim.set_v(np.random.default_rng(seed).integers(0, 2, (n, n), dtype=np.uint8))
    for k, v in overrides.items():
        getattr(sim, f'update_{k}')(v)
    if shadow:
        sim.neutral_enable()
    return sim, kw


def sim_from_scan(scan, config_idx, N=None, seed=0, runs_dir=None,
                  recipe='init', seed_v=True, shadow=False, **overrides):
    """A scan row as a live sim: evoca_from_scan (into a throwaway temp dir
    unless runs_dir is given, so Runs/ is never littered) -> load_recipe.
    Returns (sim, display_kwargs)."""
    import evoca_explore as ee
    _check_overrides(overrides)
    tmp = None
    if runs_dir is None:
        runs_dir = tmp = tempfile.mkdtemp(prefix='evoca_nb_scan_')
    try:
        path = ee.evoca_from_scan(scan, config_idx, runs_dir=runs_dir)
        return load_recipe(path, N=N, seed=seed, recipe=recipe, seed_v=seed_v,
                           shadow=shadow, **overrides)
    finally:
        if tmp is not None:
            shutil.rmtree(tmp, ignore_errors=True)


def _is_nan(v):
    return isinstance(v, (float, np.floating)) and np.isnan(v)


def params_from_row(row):
    """Scan results row (pandas Series) -> build_sim params: every
    param_<key> column for key in _PARAM_KEYS; restricted_mu as int,
    the rest as float. A missing column -- or a NaN cell, which is how
    pandas fills a column a row never had (error rows, concatenated old and
    new scans) -- is left out (build_sim default), apart from
    evoca_explore.correct_scan_row_params (a pre-A9 row without
    param_food_cap gets the standard food_cap=1)."""
    import evoca_explore as ee
    p = {k: row['param_' + k] for k in ee._PARAM_KEYS
         if 'param_' + k in row.index and not _is_nan(row['param_' + k])}
    p['restricted_mu'] = int(p.get('restricted_mu', 0))
    p = {k: float(v) if k != 'restricted_mu' else v for k, v in p.items()}
    return ee.correct_scan_row_params(p)


# ── recipes ───────────────────────────────────────────────────────────

_DATED = re.compile(r'^\d{4}-\d{2}-\d{2}_(.+\.evoca)$')


def _init_defaults():
    return {k: p.default for k, p in inspect.signature(EvoCA.init).parameters.items()
            if p.default is not inspect.Parameter.empty}


def _recipe_diff(old, new):
    """Fields where importing `old` would differ from `new` (a missing
    metaparam counts as its EvoCA.init default; v_density, a realized
    draw, and 'created' are ignored)."""
    out = []
    dflt = _init_defaults()
    for sec in ('metaparams_init', 'metaparams_final'):
        a, b = old.get(sec, {}), new.get(sec, {})
        for k in sorted(set(a) | set(b)):
            va, vb = a.get(k, dflt.get(k)), b.get(k, dflt.get(k))
            if va != vb:
                out.append(f'{sec}.{k}: {va!r} (existing) vs {vb!r}')
    a, b = old.get('initialization', {}), new.get('initialization', {})
    for k in sorted((set(a) | set(b)) - {'v_density'}):
        if a.get(k) != b.get(k):
            out.append(f'initialization.{k}: {a.get(k)!r} (existing) vs {b.get(k)!r}')
    for k in ('N', 'display'):
        if old.get(k) != new.get(k):
            out.append(f'{k}: {old.get(k)!r} (existing) vs {new.get(k)!r}')
    return out


def export_recipe_once(sim, descriptor, workstream=None, runs_dir=None, **kw):
    """sim.export_recipe(descriptor, ...) unless a recipe with that
    descriptor (any date prefix) is already in the Runs dir; then return
    the newest existing one and write nothing (warns if it would import
    differently from `sim`). kw: probes, colormode. Returns the path."""
    if runs_dir is None:
        runs_dir = evoca_paths.runs_dir(workstream)
    tmp = tempfile.mkdtemp(prefix='evoca_nb_export_')
    try:
        fresh = sim.export_recipe(descriptor, runs_dir=tmp, **kw)
        suffix = _DATED.match(os.path.basename(fresh)).group(1)
        hits = sorted(f for f in (os.listdir(runs_dir) if os.path.isdir(runs_dir) else [])
                      if (m := _DATED.match(f)) and m.group(1) == suffix)
        if not hits:
            os.makedirs(runs_dir, exist_ok=True)
            dest = os.path.join(runs_dir, os.path.basename(fresh))
            shutil.move(fresh, dest)
            return dest
        dest = os.path.join(runs_dir, hits[-1])
        with open(dest) as f:
            old = json.load(f)
        with open(fresh) as f:
            new = json.load(f)
        diff = _recipe_diff(old, new)
        if diff:
            warnings.warn(f"{os.path.basename(dest)} kept, but it differs from the "
                          f"current sim: " + '; '.join(diff), stacklevel=2)
        return dest
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ── misc ──────────────────────────────────────────────────────────────

def slope(t, y):
    """Least-squares slope of y on t (0.0 for fewer than two points)."""
    return float(np.polyfit(t, y, 1)[0]) if len(t) > 1 else 0.0
