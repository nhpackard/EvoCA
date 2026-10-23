"""Galleries: named, indexed lists of runnable EvoCA configurations.

A Gallery is the one interface to a list of recipes or scan rows, so a
lab notebook never keeps a private tuple list or hunts for a recipe name:

    G = Gallery.open('rd_pareto_picks')      # Research/<WS>/Galleries/<name>.json
    G                                        # index table; the integer index is what you pass
    G[3], G['cfg246']                        # Entry by position or by label
    G.show(3, N=512)                         # interactive SDL run (paused=False)
    sim, kw = G.load(3, N=256)               # headless sim, reproducibly seeded
    G.split('layer')[1]                      # sub-gallery: every entry with meta layer == 1
    R = G.run(N=128, n_steps=2000, snaps=(2000,))   # parallel and cached on disk
    R.tiles(t=2000, cm=1); R.strip(); R.summary(t_start=1000)

An Entry has one source: a recipe file (`recipe`, a bare name found in
Research/*/Runs) or a scan row (`scan`, `cfg`; loaded through
evoca_nb.sim_from_scan, so scan-row load-time corrections apply). `overrides`
go through the sim's update_<key>() setters after loading; `meta` holds
anything else (scores, Pareto layer, regime). Figures label each frame
"[i] label", so the index in a title is the index into the gallery.

Runs (Gallery.run) are spread over worker processes (the C library keeps
one lattice per process) and cached in REPO_ROOT/.cache/gallery_runs/,
keyed by the entry's source content, the run arguments, and a fingerprint
of the C library and the Python loaders, so a stale result is never
returned after the model or a load-time correction changes.
"""
import contextlib
import copy
import dataclasses
import hashlib
import io
import json
import os
import pickle
import re
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime
import multiprocessing as mp

import numpy as np
import pandas as pd

import evoca_paths

CACHE_DIR = os.path.join(evoca_paths.REPO_ROOT, '.cache', 'gallery_runs')
_PARAM_COLS = ('food_inc', 'm_scale', 'gdiff', 'tax', 'mu_lut', 'mu_egene',
               'restricted_mu', 'food_cap', 'food_min', 'mouthful_min')


def galleries_dir(workstream):
    return os.path.join(evoca_paths.RESEARCH, evoca_paths._check(workstream), 'Galleries')


def all_galleries(workstream=None):
    """Every saved gallery file (Research/<WS>/Galleries/*.json), one workstream or all."""
    wss = (workstream,) if workstream else evoca_paths.WORKSTREAMS
    out = []
    for ws in wss:
        d = galleries_dir(ws)
        if os.path.isdir(d):
            out += [os.path.join(d, f) for f in sorted(os.listdir(d)) if f.endswith('.json')]
    return out


def scan_label(scan, cfg):
    """Default label for a scan row: '<scan name without date>.cfg<idx>',
    e.g. 'push_edges.cfg246'. Recipes exported from a scan get the same label."""
    name = os.path.basename(os.path.normpath(str(scan)))
    return f"{re.sub(r'^\d{4}-\d{2}-\d{2}_', '', name)}.cfg{int(cfg)}"


def _find_recipe(name):
    if os.path.exists(name):
        return name
    hits = [r for r in evoca_paths.all_runs() if os.path.basename(r) == name]
    if len(hits) != 1:
        raise FileNotFoundError(f"{name!r}: {len(hits)} matches in Research/*/Runs/")
    return hits[0]


# ── entries ───────────────────────────────────────────────────────────

@dataclass
class Entry:
    """One runnable configuration: a recipe file or a scan row, plus overrides."""
    label: str = ''
    recipe: str = None            # bare .evoca name (or a path)
    scan: str = None              # bare scan name or path ...
    cfg: int = None               # ... and its config_idx
    params: str = 'init'          # which recipe metaparams: 'init' or 'final'
    title: str = ''
    overrides: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        if (self.recipe is None) == (self.scan is None):
            raise ValueError("an Entry needs exactly one of recipe= or scan= (+ cfg=)")
        if self.scan is not None:
            if self.cfg is None:
                raise ValueError("a scan Entry needs cfg=")
            self.cfg = int(self.cfg)
        if not self.label:
            self.label = self._default_label()

    def _default_label(self):
        if self.scan is not None:
            return scan_label(self.scan, self.cfg)
        d = self.recipe_json()
        if d.get('source_scan') is not None and d.get('source_config_idx') is not None:
            return scan_label(d['source_scan'], d['source_config_idx'])
        return re.sub(r'^\d{4}-\d{2}-\d{2}_', '', os.path.basename(self.recipe)).replace('.evoca', '')

    @property
    def source(self):
        return self.recipe if self.recipe is not None else f"{os.path.basename(str(self.scan))}#{self.cfg}"

    def recipe_json(self):
        """The recipe file's contents (recipe entries only)."""
        with open(_find_recipe(self.recipe)) as f:
            return json.load(f)

    def scan_row(self):
        """The scan results row behind this entry (pandas Series), or None
        for a hand-found recipe with no source scan."""
        import evoca_explore as ee
        if self.scan is not None:
            scan, cfg = self.scan, self.cfg
        else:
            d = self.recipe_json()
            if d.get('source_scan') is None:
                return None
            scan, cfg = os.path.join(evoca_paths.REPO_ROOT, d['source_scan']), d['source_config_idx']
        df = pd.read_csv(os.path.join(ee._resolve_scan_dir(scan), 'results.csv'))
        return df.set_index('config_idx').loc[int(cfg)]

    def load(self, N=None, seed=0, seed_v=True, shadow=False, quiet=False, **overrides):
        """(sim, display_kwargs): evoca_nb.load_recipe for a recipe entry,
        evoca_nb.sim_from_scan for a scan entry; the entry's overrides,
        then these, go through the update_<key>() setters. import_run's legacy
        food_cap notice is suppressed when food_cap is set here, or quiet=True."""
        from evoca_nb import load_recipe, sim_from_scan
        ov = {**self.overrides, **overrides}
        hush = quiet or 'food_cap' in ov
        with contextlib.redirect_stdout(io.StringIO()) if hush else contextlib.nullcontext():
            if self.recipe is not None:
                return load_recipe(self.recipe, N=N, seed=seed, recipe=self.params,
                                   seed_v=seed_v, shadow=shadow, **ov)
            return sim_from_scan(self.scan, self.cfg, N=N, seed=seed, recipe=self.params,
                                 seed_v=seed_v, shadow=shadow, **ov)

    def effective_params(self):
        """The metaparameters a load actually sets (corrections and overrides applied)."""
        sim, _ = self.load(N=8, quiet=True)
        try:
            return {k: getattr(sim, k) for k in _PARAM_COLS}
        finally:
            sim.free()

    def show(self, N=256, cell_px=2, paused=False, seed=0, **overrides):
        """Open the interactive SDL window (controls.run_with_controls)."""
        from controls import run_with_controls
        sim, kw = self.load(N=N, seed=seed, **overrides)
        return run_with_controls(sim, cell_px=cell_px, paused=paused, **kw)

    def fingerprint(self):
        """What a run of this entry depends on: the source's content (not its
        label or title), the metaparams choice and the overrides."""
        if self.recipe is not None:
            with open(_find_recipe(self.recipe), 'rb') as f:
                src = hashlib.sha256(f.read()).hexdigest()
        else:
            import evoca_explore as ee
            d = ee._resolve_scan_dir(self.scan)
            with open(os.path.join(d, 'scan_config.json'), 'rb') as f:
                cfg_h = hashlib.sha256(f.read()).hexdigest()
            src = [cfg_h, self.scan_row().to_json()]
        return {'src': src, 'params': self.params, 'overrides': self.overrides}

    def to_dict(self):
        return {k: v for k, v in dataclasses.asdict(self).items() if v not in (None, '', {})}


# ── galleries ─────────────────────────────────────────────────────────

class Gallery:
    """An ordered, named list of Entries; see the module docstring."""

    def __init__(self, name, entries, description='', workstream=None, owner=''):
        self.name, self.description = name, description
        self.workstream, self.owner = workstream, owner
        self.entries = [e if isinstance(e, Entry) else Entry(**e) for e in entries]
        labels = [e.label for e in self.entries]
        dup = sorted({l for l in labels if labels.count(l) > 1})
        if dup:
            raise ValueError(f"gallery {name!r}: duplicate labels {dup}")

    # -- construction --------------------------------------------------

    @classmethod
    def from_scan(cls, name, scan, cfgs, metas=None, **kw):
        """Entries for config_idx values of one scan (no recipe files written)."""
        metas = metas or [{}] * len(cfgs)
        return cls(name, [Entry(scan=scan, cfg=int(c), meta=dict(m)) for c, m in zip(cfgs, metas)], **kw)

    @classmethod
    def open(cls, name):
        """A saved gallery by name (searched in Research/*/Galleries) or path."""
        path = name if os.path.exists(name) else None
        if path is None:
            fn = name if name.endswith('.json') else name + '.json'
            hits = [p for p in all_galleries() if os.path.basename(p) == fn]
            if len(hits) != 1:
                raise FileNotFoundError(f"gallery {name!r}: {len(hits)} matches in Research/*/Galleries/")
            path = hits[0]
        with open(path) as f:
            d = json.load(f)
        return cls(d['name'], d['entries'], description=d.get('description', ''),
                   workstream=d.get('workstream'), owner=d.get('owner', ''))

    def to_dict(self):
        return {'name': self.name, 'description': self.description,
                'workstream': self.workstream, 'owner': self.owner,
                'entries': [e.to_dict() for e in self.entries]}

    def save(self, workstream=None):
        """Write Research/<WS>/Galleries/<name>.json, only if the content
        changed (so a Run All does not churn the file). Returns the path."""
        ws = workstream or self.workstream
        if ws is None:
            raise ValueError("save() needs a workstream (argument or Gallery(workstream=...))")
        self.workstream = ws
        path = os.path.join(galleries_dir(ws), self.name + '.json')
        new = self.to_dict()
        if os.path.exists(path):
            with open(path) as f:
                old = json.load(f)
            if {k: v for k, v in old.items() if k != 'saved'} == new:
                return path
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            json.dump({**new, 'saved': datetime.now().isoformat(timespec='seconds')}, f, indent=1)
        return path

    # -- access --------------------------------------------------------

    def __len__(self):
        return len(self.entries)

    def __iter__(self):
        return iter(self.entries)

    def index(self, key):
        """Position of a label. Exact match first, then a unique match on a
        trailing '.'-separated part ('cfg246' finds 'push_edges.cfg246')."""
        if isinstance(key, (int, np.integer)):
            k = int(key)
            if not -len(self) <= k < len(self):
                raise IndexError(f"{self.name}: index {k} out of range ({len(self)} entries)")
            return k % len(self)
        labels = [e.label for e in self.entries]
        if key in labels:
            return labels.index(key)
        hits = [i for i, l in enumerate(labels) if l.endswith('.' + key)]
        if len(hits) == 1:
            return hits[0]
        raise KeyError(f"{self.name}: label {key!r} matches {[labels[i] for i in hits] or 'nothing'}")

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self._sub(range(len(self))[key])
        if isinstance(key, (list, tuple, np.ndarray, range, pd.Index)):
            key = list(key)
            if key and all(isinstance(k, (bool, np.bool_)) for k in key):
                key = [i for i, b in enumerate(key) if b]
            return self._sub([self.index(k) for k in key])
        return self.entries[self.index(key)]

    def _sub(self, idx, name=None):
        return Gallery(name or self.name, [copy.deepcopy(self.entries[i]) for i in idx],
                       description=self.description, workstream=self.workstream, owner=self.owner)

    def where(self, fn=None, **meta):
        """Sub-gallery of entries with fn(entry) true and entry.meta[k] == v for each k=v."""
        keep = [i for i, e in enumerate(self.entries)
                if (fn is None or fn(e)) and all(e.meta.get(k) == v for k, v in meta.items())]
        return self._sub(keep)

    def split(self, key):
        """{value: sub-gallery} over the values of meta[key], sorted, e.g.
        split('layer')[1] is the Pareto-layer-1 sub-gallery."""
        vals = sorted({e.meta.get(key) for e in self.entries if key in e.meta})
        return {v: self._sub([i for i, e in enumerate(self.entries) if e.meta.get(key) == v],
                             name=f'{self.name}[{key}={v}]') for v in vals}

    def load(self, key, **kw):
        return self[key].load(**kw)

    def show(self, key, **kw):
        return self[key].show(**kw)

    def labels(self):
        return [e.label for e in self.entries]

    def table(self, params=True):
        """One row per entry, index = gallery position. params=True adds the
        metaparams a load actually sets (corrections and overrides applied)."""
        rows = []
        for e in self.entries:
            r = {'label': e.label, 'title': e.title, 'source': e.source}
            if e.params != 'init':
                r['params'] = e.params
            if e.overrides:
                r['overrides'] = ', '.join(f'{k}={v}' for k, v in e.overrides.items())
            r.update(e.meta)
            if params:
                r.update(e.effective_params())
            rows.append(r)
        df = pd.DataFrame(rows)
        return df.loc[:, [c for c in df.columns if df[c].astype(str).ne('').any()]]

    def __repr__(self):
        return f"Gallery({self.name!r}, {len(self)} entries)\n" + self.table(params=False).to_string()

    def _repr_html_(self):
        head = f"<b>Gallery <code>{self.name}</code></b>, {len(self)} entries"
        if self.description:
            head += f": {self.description}"
        with pd.option_context('display.max_colwidth', 60):
            return head + self.table()._repr_html_()

    # -- running -------------------------------------------------------

    def run(self, N=128, n_steps=2000, every=25, snaps=(), cms=(0, 1), seed=0, seed_v=True,
            metrics='spatial', overrides=None, workers=None, cache=True, progress=True):
        """Run every entry headless; returns a GalleryRun aligned with the gallery.

        metrics: 'spatial' (evoca_explore._spatial_metrics every `every` steps,
        including alive_density) or 'alive' (alive_density only, cheap).
        snaps: steps at which to keep RGB frames for colormodes `cms`.
        overrides: applied to every entry on top of its own.
        workers: processes (default: cpu count - 1); 1 runs in this process."""
        spec = dict(N=int(N), n_steps=int(n_steps), every=int(every), snaps=sorted(int(t) for t in snaps),
                    cms=[int(c) for c in cms], seed=int(seed), seed_v=bool(seed_v), metrics=metrics,
                    overrides=dict(overrides or {}))
        if metrics not in ('spatial', 'alive'):
            raise ValueError(f"metrics={metrics!r}: 'spatial' or 'alive'")
        jobs = [(e.to_dict(), spec) for e in self.entries]
        keys = [_cache_key(e, spec) for e in self.entries]
        tracks = [_cache_get(k) if cache else None for k in keys]
        todo = [i for i, t in enumerate(tracks) if t is None]
        if todo:
            n_w = workers or max(1, (os.cpu_count() or 2) - 1)
            n_w = min(n_w, len(todo))
            if progress:
                print(f"{self.name}: running {len(todo)} of {len(self)} entries "
                      f"({len(self) - len(todo)} cached), {n_w} worker(s)", flush=True)
            if n_w == 1:
                results = [_run_job(jobs[i]) for i in todo]
            else:
                with ProcessPoolExecutor(n_w, mp_context=mp.get_context('spawn')) as ex:
                    results = list(ex.map(_run_job, [jobs[i] for i in todo]))
            for i, r in zip(todo, results):
                tracks[i] = r
                if cache:
                    _cache_put(keys[i], r)
        return GalleryRun(self, spec, tracks)


# ── running (module level, so spawn workers can import it) ───────────

def _run_job(job):
    from evoca_explore import _spatial_metrics
    from evoca_nb import rgb
    ed, spec = job
    e = Entry(**ed)
    sim, _ = e.load(N=spec['N'], seed=spec['seed'], seed_v=spec['seed_v'], quiet=True, **spec['overrides'])
    N, snaps = sim.N, set(spec['snaps'])
    ts, rows, shots = [], [], {}
    try:
        if 0 in snaps:
            shots[0] = {cm: rgb(sim, cm) for cm in spec['cms']}
        for t in range(1, spec['n_steps'] + 1):
            sim.step()
            if t % spec['every'] == 0:
                al = sim.get_alive().reshape(N, N)
                if spec['metrics'] == 'spatial':
                    rows.append(_spatial_metrics(al.copy(), sim.get_F().reshape(N, N).copy()))
                else:
                    rows.append({'alive_density': float(al.mean())})
                ts.append(t)
            if t in snaps:
                shots[t] = {cm: rgb(sim, cm) for cm in spec['cms']}
    finally:
        sim.free()
    return Track(pd.DataFrame(rows, index=pd.Index(ts, name='t')), shots)


_CODE_FP = None


def _code_fingerprint():
    """Hash of the C library and the Python files a run goes through."""
    global _CODE_FP
    if _CODE_FP is None:
        h = hashlib.sha256()
        root = evoca_paths.REPO_ROOT
        lib = [os.path.join(root, 'C', f) for f in ('libevoca.dylib', 'libevoca.so')]
        py = [os.path.join(root, 'python', f) for f in
              ('evoca_py.py', 'evoca_nb.py', 'evoca_explore.py', 'evoca_gallery.py')]
        for p in lib + py:
            if os.path.exists(p):
                with open(p, 'rb') as f:
                    h.update(f.read())
        _CODE_FP = h.hexdigest()
    return _CODE_FP


def _cache_key(entry, spec):
    """Overrides are keyed merged, as _run_job applies them (run-time ones win), so an
    override set on the entry and the same one passed to run() share a key."""
    fp = entry.fingerprint()
    fp['overrides'] = {**fp['overrides'], **spec['overrides']}
    blob = json.dumps({'entry': fp, 'spec': {k: v for k, v in spec.items() if k != 'overrides'},
                       'code': _code_fingerprint()}, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:32]


def _cache_get(key):
    p = os.path.join(CACHE_DIR, key + '.pkl')
    if not os.path.exists(p):
        return None
    try:
        with open(p, 'rb') as f:
            return pickle.load(f)
    except Exception:
        return None


def _cache_put(key, track):
    os.makedirs(CACHE_DIR, exist_ok=True)
    tmp = os.path.join(CACHE_DIR, f'{key}.{os.getpid()}.tmp')
    with open(tmp, 'wb') as f:
        pickle.dump(track, f)
    os.replace(tmp, os.path.join(CACHE_DIR, key + '.pkl'))


def clear_cache():
    """Delete every cached gallery run."""
    if os.path.isdir(CACHE_DIR):
        for f in os.listdir(CACHE_DIR):
            os.remove(os.path.join(CACHE_DIR, f))


# ── results ───────────────────────────────────────────────────────────

@dataclass
class Track:
    ts: pd.DataFrame              # metrics, index = step
    shots: dict                   # {step: {colormode: (N, N, 3) uint8}}


_CM_NAMES = None


def _cm_name(cm):
    global _CM_NAMES
    if _CM_NAMES is None:
        from colormodes import COLOR_MODES
        _CM_NAMES = COLOR_MODES
    return _CM_NAMES[cm] if 0 <= cm < len(_CM_NAMES) else str(cm)


class GalleryRun:
    """Results of Gallery.run, aligned with the gallery: run[i] / run['cfg246'] is a Track."""

    def __init__(self, gallery, spec, tracks):
        self.gallery, self.spec, self.tracks = gallery, spec, tracks

    def __len__(self):
        return len(self.tracks)

    def __getitem__(self, key):
        return self.tracks[self.gallery.index(key)]

    def tag(self, i):
        return f"[{i}] {self.gallery.entries[i].label}"

    def summary(self, t_start=0):
        """Per-entry means over steps > t_start; index = gallery position."""
        rows = []
        for i, tr in enumerate(self.tracks):
            d = tr.ts[tr.ts.index > t_start]
            r = {'label': self.gallery.entries[i].label,
                 'alive_mean': d.alive_density.mean(), 'alive_std': d.alive_density.std(),
                 'alive_final': tr.ts.alive_density.iloc[-1]}
            if 'correlation_length' in d:
                r.update(corr_len_mean=d.correlation_length.mean(), patch_std=d.largest_patch.std(),
                         n_patches_mean=d.n_patches.mean())
            rows.append(r)
        return pd.DataFrame(rows)

    def tiles(self, t=None, cm=1, idx=None, ncols=8, size=2.0, title=None):
        """One frame per entry (step t, default the last snap); titles '[i] label'
        plus title(i, entry) if given."""
        import matplotlib.pyplot as plt
        t = self.spec['snaps'][-1] if t is None else t
        idx = list(range(len(self))) if idx is None else [self.gallery.index(k) for k in idx]
        nrows = -(-len(idx) // ncols)
        fig, axs = plt.subplots(nrows, ncols, figsize=(size * ncols, size * 1.15 * nrows), squeeze=False)
        for ax in axs.flat:
            ax.axis('off')
        for ax, i in zip(axs.flat, idx):
            ax.imshow(self.tracks[i].shots[t][cm])
            extra = '' if title is None else '\n' + title(i, self.gallery.entries[i])
            ax.set_title(self.tag(i) + extra, fontsize=7)
        fig.suptitle(f"{self.gallery.name}: t={t}, {_cm_name(cm)}, N={self.spec['N']}", fontsize=9)
        plt.tight_layout()
        return fig

    def strip(self, idx=None, snaps=None, cms=None, size=1.75):
        """Rows = entries, columns = snaps x colormodes; row label '[i] label' + title."""
        import matplotlib.pyplot as plt
        snaps = self.spec['snaps'] if snaps is None else snaps
        cms = self.spec['cms'] if cms is None else cms
        idx = list(range(len(self))) if idx is None else [self.gallery.index(k) for k in idx]
        nc = len(snaps) * len(cms)
        fig, axs = plt.subplots(len(idx), nc, figsize=(size * nc, size * 1.05 * len(idx)), squeeze=False)
        for r, i in enumerate(idx):
            for j, t in enumerate(snaps):
                for k, cm in enumerate(cms):
                    ax = axs[r, len(cms) * j + k]
                    ax.imshow(self.tracks[i].shots[t][cm]); ax.set_xticks([]); ax.set_yticks([])
                    if r == 0:
                        ax.set_title(f't={t}  {_cm_name(cm)}', fontsize=7)
            e = self.gallery.entries[i]
            axs[r, 0].set_ylabel(self.tag(i) + ('\n' + e.title if e.title else ''), fontsize=6,
                                 rotation=0, ha='right', va='center', labelpad=4)
        plt.tight_layout()
        return fig
