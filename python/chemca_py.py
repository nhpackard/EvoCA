"""ChemCA: Python wrapper for the metabolism-based agent model (C/chemca.c).

Design: Docs/2026-10-23_chem-model-v2.md §4-10; choices and status:
Docs/plan-chemca.md. The interface parallels EvoCA where it makes sense:

    sim = ChemCA()
    sim.init(128, preset='gs3', mu_rgene=0.01)
    sim.state(alive='halfplane', seed=1)
    sim.step(1000)
    sim.get_alive(); sim.get_E('A'); sim.get_store('F'); sim.get_e()

Like EvoCA, the C state is global to the library: one live ChemCA at a
time per process.

Metaparameters. DEFAULTS is the one statement of every default. Scalar
keys are global; the per-compound fields COMPOUND_FIELDS (free energy G,
diffusion, transport u, toxicity tox, open-mode inc and sink) are keyed
'<field>_<compound>' (e.g. 'tox_W') and take DEFAULTS['compound'] unless a
preset overrides them. Every key list is derived from these (param_keys,
defaults); none is maintained by hand.

Presets (PRESETS) fix the compound universe, the light reaction, the
initial genes and the preset's provisional rate values. All values are
PROVISIONAL (chosen so a population persists for a few thousand steps at
N=128); a separate Gray-Scott regime study will supply better ones.
"""
import ctypes
import hashlib
import json
import os
from datetime import datetime

import numpy as np

from evoca_py import _find_lib
import sim_model

# ── Layout (mirrors C/chemca.h; checked against the library at load) ──

NSLOT = 8
GENE_BYTES = 5
NK = 64
C_NONE = 0xFF
KMAX = 8
G_A, G_B, G_C, G_M, G_K = range(GENE_BYTES)

# Energy-budget entries, in the order of the C enum (CHEMCA_EN_*).
ENERGY_KEYS = ('light', 'yield', 'cost', 'death', 'split', 'open_in', 'open_out')

COMPOUND_FIELDS = ('G', 'diff', 'u', 'tox', 'inc', 'sink')   # chemca_set_compound order

# ── Defaults: the one source ─────────────────────────────────────────

DEFAULTS = dict(
    preset='fw2',
    k_min=1e-3,        # k = k_min * 2^(k_idx/4), k_idx in 0..63 (1e-3 .. 55)
    m_e=0.0,           # maintenance energy per step
    c_g=0.0,           # energy cost per active gene per step
    mu_rgene=0.0,      # point mutations per slot per birth (Poisson(mu * 8))
    mu_ractive=0.0,    # active-bit flips per slot per birth
    p_dup=0.5,         # duplication onto a newly active slot
    open_mode=0,       # 1: per-compound uniform inc and proportional sink
    light=0.0,         # rate of the preset's environment (light) reaction
    repro_gate='energy',  # 'energy': e >= 1; 'mass': body mass >= mass_repro
    mass_repro=1.0,    # body-mass threshold of the mass gate (body: u_c == 0)
    compound=dict(G=0.0, diff=0.0, u=0.0, tox=0.0, inc=0.0, sink=0.0),
)
assert tuple(DEFAULTS['compound']) == COMPOUND_FIELDS

# State-initialization keys and their defaults (recorded in recipes).
STATE_DEFAULTS = dict(
    alive='halfplane', alive_fraction=0.5, alive_axis=0,
    mbar=1.0,          # closed mode: mean mass per site, split by the preset
    mass_split=None,   # None: the preset's split
    store_frac=0.5,    # alive cells start with this fraction of their site's E
    e_init=0.5,
    genes='preset',    # 'preset' | 'random'
    k_idx0=None,       # None: the preset's k_idx; int: all preset genes; list: per gene
    E_noise=0.0,       # multiplicative uniform noise on the initial E
    seed=None,         # set_seed + numpy draws; recorded for reproducibility
)


def gene(A, B, C=None, m=0, k_idx=16):
    """A gene as a dict of compound names (C None = uncatalyzed)."""
    return dict(A=A, B=B, C=C, m=m if C is not None else 0, k_idx=k_idx)


PRESETS = {
    # Food/waste. Default metabolism F -> W uncatalyzed; the autocatalytic
    # variant F -> W catalyzed by W (m=1) sits in slot 1, inactive, and is
    # reachable by mutation (active-bit or (C, m) change). Light W -> F.
    'fw2': dict(
        compounds=('F', 'W'),
        mass_split=(0.9, 0.1),
        light_rx=('W', 'F'),
        genes=[gene('F', 'W', None, 0, 29), gene('F', 'W', 'W', 1, 29)],
        active=0b01,
        params=dict(light=0.015, m_e=0.008, c_g=0.002,
                    G_F=1.0, G_W=0.0,
                    diff_F=0.1, diff_W=0.1, u_F=0.2, u_W=0.2, tox_W=0.6),
    ),
    # Three-species closed Gray-Scott loop (doc status section): inside
    # agents F + 2A -> 3A (cubic autocatalysis; F -> A catalyzed by A, m=2)
    # and A -> P; light P -> F. Rates: set S3 ("spiral") of the closed-loop
    # regime study (Research/Chem/Experiments/2026-10-24_gs_closed_regimes):
    # k1 = 6.0, k2 = 0.1, r = 0.005, diffusion (0.6, 0.3, 0.3), at the nearest
    # grid points k_idx 50 (k1 = 5.79) and 27 (k2 = 0.108). A = 0 is
    # absorbing, so alive cells start with A in their stores (mass_split).
    # Toxicity on P (tox_P 0.6, rarely binding here; set tox_A for the
    # option on A). m_e 0.005 and u 0.5 are the agent-side provisional values.
    'gs3': dict(
        compounds=('F', 'A', 'P'),
        mass_split=(0.8, 0.1, 0.1),
        light_rx=('P', 'F'),
        genes=[gene('F', 'A', 'A', 2, 50), gene('A', 'P', None, 0, 27)],
        active=0b11,
        params=dict(light=0.005, m_e=0.001, c_g=0.002,
                    G_F=1.0, G_A=0.5, G_P=0.0,
                    diff_F=0.6, diff_A=0.3, diff_P=0.3,
                    u_F=0.5, u_A=0.5, u_P=0.5, tox_P=0.6),
    ),
}


def param_keys(preset):
    """All metaparameter keys of a preset (derived from DEFAULTS)."""
    scal = [k for k in DEFAULTS if k != 'compound']
    comp = [f'{f}_{c}' for c in PRESETS[preset]['compounds'] for f in COMPOUND_FIELDS]
    return scal + comp


def defaults(preset):
    """Full default metaparameters of a preset: DEFAULTS, then the
    per-compound template, then the preset's overrides."""
    d = {k: v for k, v in DEFAULTS.items() if k != 'compound'}
    d['preset'] = preset
    for c in PRESETS[preset]['compounds']:
        for f, v in DEFAULTS['compound'].items():
            d[f'{f}_{c}'] = v
    d.update(PRESETS[preset]['params'])
    assert set(d) == set(param_keys(preset)), "preset overrides an unknown key"
    return d


# ── ctypes signatures ────────────────────────────────────────────────

_u8p = ctypes.POINTER(ctypes.c_uint8)
_f32p = ctypes.POINTER(ctypes.c_float)
_f64p = ctypes.POINTER(ctypes.c_double)
_i32p = ctypes.POINTER(ctypes.c_int32)
_ip = ctypes.POINTER(ctypes.c_int)
_F, _I, _L, _D = ctypes.c_float, ctypes.c_int, ctypes.c_long, ctypes.c_double

# Scalar metaparameters pushed by a set_/get_ pair named after the key.
_SCALAR_C = {'k_min': _F, 'm_e': _F, 'c_g': _F, 'mu_rgene': _F,
             'mu_ractive': _F, 'p_dup': _F, 'open_mode': _I, 'mass_repro': _F}
REPRO_GATES = ('energy', 'mass')   # index = the C value of chemca_set_repro_gate

_SIGS = {
    'init': ([_I, _I], _I), 'free': ([], None), 'set_seed': ([ctypes.c_uint32], None),
    'get_N': ([], _I), 'get_K': ([], _I), 'get_step': ([], _L),
    'kmax': ([], _I), 'nslot': ([], _I), 'gene_bytes': ([], _I), 'nk': ([], _I),
    'energy_count': ([], _I), 'k_of': ([_I], _F),
    'set_compound': ([_I, _F, _F, _F, _F, _F, _F], _I), 'get_compound': ([_I, _f32p], _I),
    'env_rx_clear': ([], None), 'env_rx_add': ([_I, _I, _F], _I),
    'env_rx_set_rate': ([_I, _F], _I), 'env_rx_set_field': ([_I, _f32p], _I),
    'env_rx_count': ([], _I), 'env_rx_get': ([_I, _ip, _ip, _f32p, _ip], _I),
    'set_alive': ([_u8p], None), 'set_E': ([_I, _f32p], None),
    'set_store': ([_I, _f32p], None), 'set_e': ([_f32p], None),
    'set_genes_all': ([_u8p, ctypes.c_uint8], _I), 'set_genes': ([_u8p, _u8p], _I),
    'set_genes_random': ([], None),
    'step': ([], None), 'step_n': ([_I], None),
    'get_alive': ([_u8p], None), 'get_births': ([_u8p], None),
    'get_E': ([_I, _f32p], None), 'get_store': ([_I, _f32p], None),
    'get_e': ([_f32p], None), 'get_genes': ([_u8p], None), 'get_active': ([_u8p], None),
    'get_mass': ([_f64p, _f64p], None), 'get_mass_flux': ([_I, _f64p], _I),
    'get_energy_budget': ([_f64p, _f64p], None), 'get_e_total': ([], _D),
    'get_population': ([], _L), 'get_births_last': ([], _L), 'get_deaths_last': ([], _L),
    'get_tox_deaths_last': ([], _L), 'get_energy_deaths_last': ([], _L),
    'get_evictions_last': ([], _L),
    'colorize': ([_i32p, _I], None), 'colormode_count': ([], _I),
    'set_color_scale': ([_F], None), 'get_color_scale': ([], _F),
    'set_repro_gate': ([_I], None), 'get_repro_gate': ([], _I),
    **{f'set_{k}': ([t], None) for k, t in _SCALAR_C.items()},
    **{f'get_{k}': ([], t) for k, t in _SCALAR_C.items()},
}


def _ptr(a, p):
    return a.ctypes.data_as(p)


class ChemCA:
    """ctypes wrapper around C/chemca.c (see module docstring)."""

    COLORMODE_NAMES = ('env RGB', 'alive + stores', 'energy', 'k_idx (first active gene)')

    def __init__(self, lib_path=None):
        self._lib = ctypes.CDLL(lib_path or _find_lib())
        for name, (args, res) in _SIGS.items():
            fn = getattr(self._lib, f'chemca_{name}')
            fn.argtypes, fn.restype = args, res
        L = self._lib
        layout = (L.chemca_nslot(), L.chemca_gene_bytes(), L.chemca_nk(), L.chemca_kmax(),
                  L.chemca_energy_count())
        assert layout == (NSLOT, GENE_BYTES, NK, KMAX, len(ENERGY_KEYS)), layout
        self._N = 0
        self._params = None
        self._init_params = None
        self._state_params = {}
        self._mass0 = None

    # ── Lifecycle ────────────────────────────────────────────────────

    def init(self, N, preset=None, **params):
        """Allocate an N x N world for a preset and push every
        metaparameter. Unknown keys raise (none is dropped silently)."""
        preset = preset or params.pop('preset', None) or DEFAULTS['preset']
        if preset not in PRESETS:
            raise ValueError(f"unknown preset {preset!r}; known: {sorted(PRESETS)}")
        p = defaults(preset)
        unknown = sorted(set(params) - set(p))
        if unknown:
            raise ValueError(f"unknown ChemCA params for preset {preset!r}: {unknown}")
        p.update(params)
        self._check(p)
        self._N = int(N)
        self._params = p
        self._apply_params()

    def reinit(self):
        """Re-allocate the world from the current metaparameters (empty
        lattice; replay state() afterwards)."""
        self._apply_params()

    @staticmethod
    def _check(p):
        if p['repro_gate'] not in REPRO_GATES:
            raise ValueError(f"repro_gate must be one of {REPRO_GATES}, got {p['repro_gate']!r}")
        if p['repro_gate'] == 'mass' and not ChemCA._body(p):
            raise ValueError("repro_gate='mass' needs a body compound (some u_<c> = 0)")

    @staticmethod
    def _body(p):
        return [c for c in PRESETS[p['preset']]['compounds'] if float(p[f'u_{c}']) == 0.0]

    @property
    def body_compounds(self):
        """Compounds that count as body mass (transport rate 0)."""
        return self._body(self._params)

    def _apply_params(self):
        L, p = self._lib, self._params
        if L.chemca_init(self._N, len(self.compounds)) != 0:
            raise RuntimeError(f"chemca_init({self._N}, {len(self.compounds)}) failed")
        self._push_params()
        self._init_params = dict(p)
        self._state_params = {}
        self._mass0 = None

    def _push_params(self):
        """Push every metaparameter in self._params to C (no reallocation)."""
        L, p, pre = self._lib, self._params, self.preset_spec
        for k in _SCALAR_C:
            getattr(L, f'chemca_set_{k}')(p[k])
        L.chemca_set_repro_gate(REPRO_GATES.index(p['repro_gate']))
        for ci, c in enumerate(self.compounds):
            if L.chemca_set_compound(ci, *[float(p[f'{f}_{c}']) for f in COMPOUND_FIELDS]):
                raise RuntimeError(f"chemca_set_compound({c}) failed")
        L.chemca_env_rx_clear()
        A, B = (self.cidx(x) for x in pre['light_rx'])
        if L.chemca_env_rx_add(A, B, float(p['light'])) != 0:
            raise RuntimeError("chemca_env_rx_add failed")

    def update(self, **kw):
        """Change metaparameters live (no reallocation; the lattice keeps
        its state). 'preset' cannot change here: use init()."""
        bad = sorted(set(kw) - (set(self._params) - {'preset'}))
        if bad:
            raise ValueError(f"cannot update {bad}")
        self._check({**self._params, **kw})
        self._params.update(kw)
        light_field = self._light_field
        self._push_params()
        if light_field is not None:
            self.set_light_field(light_field)

    _light_field = None

    def set_light_field(self, r):
        """Per-site light rate r(x) (N x N), replacing the scalar 'light';
        None returns to the scalar."""
        if r is None:
            self._light_field = None
            self._lib.chemca_env_rx_set_field(0, None)
            return
        a = np.ascontiguousarray(r, np.float32).reshape(self._N * self._N)
        self._light_field = a.reshape(self._N, self._N).copy()
        self._lib.chemca_env_rx_set_field(0, _ptr(a, _f32p))

    def free(self):
        self._lib.chemca_free()
        self._N = 0

    def set_seed(self, seed):
        """Reseed ChemCA's own C PRNG (independent of EvoCA's)."""
        self._lib.chemca_set_seed(int(seed) & 0xFFFFFFFF)

    # ── Introspection ────────────────────────────────────────────────

    @property
    def N(self):
        return self._N

    @property
    def preset(self):
        return self._params['preset']

    @property
    def preset_spec(self):
        return PRESETS[self.preset]

    @property
    def compounds(self):
        return PRESETS[self._params['preset']]['compounds']

    def cidx(self, c):
        """Compound index from a name or an int."""
        return c if isinstance(c, (int, np.integer)) else self.compounds.index(c)

    def params(self):
        """Current metaparameters, suitable for init(N, **d) (with model)."""
        return {'model': 'chemca', 'N': self._N, **self._params}

    def colormode_count(self):
        return int(self._lib.chemca_colormode_count())

    def k_of(self, k_idx):
        return float(self._lib.chemca_k_of(int(k_idx)))

    # ── State ────────────────────────────────────────────────────────

    @property
    def state_params(self):
        return dict(self._state_params)

    def _preset_gene_bytes(self, k_idx0=None):
        pre = self.preset_spec
        gs = list(pre['genes'])
        if k_idx0 is not None:
            ks = list(k_idx0) if isinstance(k_idx0, (list, tuple)) else [k_idx0] * len(gs)
            if len(ks) != len(gs):
                raise ValueError(f"k_idx0 needs {len(gs)} entries for preset {self.preset}")
            gs = [dict(g, k_idx=int(k)) for g, k in zip(gs, ks)]
        slots = [gs[s] if s < len(gs) else gs[0] for s in range(NSLOT)]
        out = np.zeros((NSLOT, GENE_BYTES), np.uint8)
        for s, g in enumerate(slots):
            out[s] = self.encode_gene(g)
        return out, pre['active']

    def encode_gene(self, g):
        C = C_NONE if g['C'] is None else self.cidx(g['C'])
        return np.array([self.cidx(g['A']), self.cidx(g['B']), C,
                         g['m'] if g['C'] is not None else 0, g['k_idx']], np.uint8)

    def decode_gene(self, b):
        """5 gene bytes -> dict with compound names."""
        nm = self.compounds
        C = None if b[G_C] == C_NONE else nm[b[G_C]]
        return dict(A=nm[b[G_A]], B=nm[b[G_B]], C=C, m=int(b[G_M]), k_idx=int(b[G_K]))

    def state(self, **kw):
        """Initialize the lattice (see STATE_DEFAULTS). Closed-mode mass is
        set here: E_c = mbar * split_c at every site, and each alive cell
        moves store_frac of its site's E into its stores, so the mean mass
        per site is exactly mbar (before E_noise). alive may be
        'halfplane' | 'fraction' | 'all' | an N x N array."""
        bad = sorted(set(kw) - set(STATE_DEFAULTS))
        if bad:
            raise ValueError(f"unknown state params {bad}")
        s = {**STATE_DEFAULTS, **kw}
        N, L, K = self._N, self._lib, len(self.compounds)
        if s['seed'] is not None:
            self.set_seed(s['seed'])
        rng = np.random.default_rng(s['seed'])

        al = s['alive']
        if isinstance(al, str):
            if al == 'halfplane':
                a = np.zeros((N, N), np.uint8)
                if int(s['alive_axis']) == 0:
                    a[:, :N // 2] = 1
                else:
                    a[:N // 2, :] = 1
            elif al == 'fraction':
                a = (rng.random((N, N)) < s['alive_fraction']).astype(np.uint8)
            elif al == 'all':
                a = np.ones((N, N), np.uint8)
            else:
                raise ValueError(f"unknown alive pattern {al!r}")
        else:
            a = (np.asarray(al).reshape(N, N) != 0).astype(np.uint8)
            s['alive'] = 'array'          # recipes record the kind, not the array
        a = np.ascontiguousarray(a.ravel())
        L.chemca_set_alive(_ptr(a, _u8p))

        split = s['mass_split'] or self.preset_spec['mass_split']
        if len(split) != K:
            raise ValueError(f"mass_split needs {K} entries")
        for ci in range(K):
            E = np.full(N * N, s['mbar'] * split[ci], np.float32)
            if s['E_noise']:
                E *= (1.0 + s['E_noise'] * (rng.random(N * N) - 0.5)).astype(np.float32)
            st = np.where(a == 1, np.float32(s['store_frac']) * E, 0).astype(np.float32)
            E = (E - st).astype(np.float32)
            L.chemca_set_E(ci, _ptr(np.ascontiguousarray(E), _f32p))
            L.chemca_set_store(ci, _ptr(np.ascontiguousarray(st), _f32p))
        e = np.where(a == 1, s['e_init'], 0).astype(np.float32)
        L.chemca_set_e(_ptr(e, _f32p))

        if s['genes'] == 'random':
            L.chemca_set_genes_random()
        elif s['genes'] == 'preset':
            g, act = self._preset_gene_bytes(s['k_idx0'])
            if L.chemca_set_genes_all(_ptr(np.ascontiguousarray(g.ravel()), _u8p), act):
                raise RuntimeError("invalid preset genes")
        else:
            raise ValueError(f"unknown genes {s['genes']!r}")
        self._state_params = s
        self._mass0 = self.mass()['total']

    # Low-level setters (constructed fixtures, tests).
    def set_alive(self, a):
        a = np.ascontiguousarray(np.asarray(a).reshape(-1) != 0, np.uint8)
        self._lib.chemca_set_alive(_ptr(a, _u8p))

    def set_E(self, c, E):
        E = np.ascontiguousarray(np.broadcast_to(np.asarray(E, np.float32),
                                                 (self._N, self._N)).reshape(-1))
        self._lib.chemca_set_E(self.cidx(c), _ptr(E, _f32p))

    def set_store(self, c, s):
        s = np.ascontiguousarray(np.broadcast_to(np.asarray(s, np.float32),
                                                 (self._N, self._N)).reshape(-1))
        self._lib.chemca_set_store(self.cidx(c), _ptr(s, _f32p))

    def set_e(self, e):
        e = np.ascontiguousarray(np.broadcast_to(np.asarray(e, np.float32),
                                                 (self._N, self._N)).reshape(-1))
        self._lib.chemca_set_e(_ptr(e, _f32p))

    def set_genes_all(self, genes, active):
        """genes: list of gene dicts (<= 8; unused slots copy slot 0)."""
        g = [genes[s] if s < len(genes) else genes[0] for s in range(NSLOT)]
        b = np.ascontiguousarray(np.concatenate([self.encode_gene(x) for x in g]))
        if self._lib.chemca_set_genes_all(_ptr(b, _u8p), int(active)):
            raise ValueError("invalid gene or empty active mask")

    def set_genes(self, genes, active):
        """Full arrays: genes (N, N, 8, 5) uint8, active (N, N) uint8."""
        g = np.ascontiguousarray(genes, np.uint8).reshape(-1)
        a = np.ascontiguousarray(active, np.uint8).reshape(-1)
        if self._lib.chemca_set_genes(_ptr(g, _u8p), _ptr(a, _u8p)):
            raise ValueError("invalid gene or empty active mask on an alive cell")

    def set_genes_random(self):
        self._lib.chemca_set_genes_random()

    # ── Dynamics ─────────────────────────────────────────────────────

    def step(self, n=1):
        self._lib.chemca_step_n(int(n))

    # ── Getters ──────────────────────────────────────────────────────

    def _getf(self, fn, *args):
        out = np.empty(self._N * self._N, np.float32)
        fn(*args, _ptr(out, _f32p))
        return out.reshape(self._N, self._N)

    def _getu8(self, fn, n_per=1):
        out = np.empty(self._N * self._N * n_per, np.uint8)
        fn(_ptr(out, _u8p))
        return out

    def get_alive(self):
        return self._getu8(self._lib.chemca_get_alive).reshape(self._N, self._N)

    def get_births(self):
        return self._getu8(self._lib.chemca_get_births).reshape(self._N, self._N)

    def get_E(self, c):
        return self._getf(self._lib.chemca_get_E, self.cidx(c))

    def get_store(self, c):
        return self._getf(self._lib.chemca_get_store, self.cidx(c))

    def get_e(self):
        return self._getf(self._lib.chemca_get_e)

    def get_rgenes(self):
        """(genes (N, N, 8, 5) uint8, active (N, N) uint8)."""
        g = self._getu8(self._lib.chemca_get_genes, NSLOT * GENE_BYTES)
        a = self._getu8(self._lib.chemca_get_active)
        return g.reshape(self._N, self._N, NSLOT, GENE_BYTES), a.reshape(self._N, self._N)

    def population(self):
        return int(self._lib.chemca_get_population())

    def get_step(self):
        return int(self._lib.chemca_get_step())

    def mass(self):
        """{'env': [K], 'store': [K], 'total': float} as doubles."""
        K = len(self.compounds)
        me, ms = np.zeros(K), np.zeros(K)
        self._lib.chemca_get_mass(_ptr(me, _f64p), _ptr(ms, _f64p))
        return {'env': me, 'store': ms, 'total': float(me.sum() + ms.sum())}

    def mass_flux(self):
        """Open-mode cumulative (in, out) per compound, shape (K, 2)."""
        K = len(self.compounds)
        out = np.zeros((K, 2))
        for c in range(K):
            self._lib.chemca_get_mass_flux(c, _ptr(out[c], _f64p))
        return out

    def free_energy(self):
        """Phi = sum_c G_c * (mass of c in env and stores)."""
        m = self.mass()
        G = np.array([self._params[f'G_{c}'] for c in self.compounds])
        return float(G @ (m['env'] + m['store']))

    def energy_budget(self, cum=True):
        """Dict over ENERGY_KEYS, cumulative since init or the last step."""
        last, total = np.zeros(len(ENERGY_KEYS)), np.zeros(len(ENERGY_KEYS))
        self._lib.chemca_get_energy_budget(_ptr(last, _f64p), _ptr(total, _f64p))
        return dict(zip(ENERGY_KEYS, total if cum else last))

    def e_total(self):
        return float(self._lib.chemca_get_e_total())

    def events_last(self):
        L = self._lib
        return dict(births=int(L.chemca_get_births_last()), deaths=int(L.chemca_get_deaths_last()),
                    tox_deaths=int(L.chemca_get_tox_deaths_last()),
                    energy_deaths=int(L.chemca_get_energy_deaths_last()),
                    evictions=int(L.chemca_get_evictions_last()))

    def state_hash(self):
        """sha256 hex of alive, E, stores, e, genes, active."""
        h = hashlib.sha256()
        g, a = self.get_rgenes()
        arrs = [self.get_alive(), self.get_e(), g, a]
        arrs += [self.get_E(c) for c in self.compounds]
        arrs += [self.get_store(c) for c in self.compounds]
        for x in arrs:
            h.update(np.ascontiguousarray(x).tobytes())
        return h.hexdigest()

    # ── Gene summaries ───────────────────────────────────────────────

    def active_genes(self):
        """Active genes of alive cells: (cell index, slot, 5 bytes) arrays."""
        g, a = self.get_rgenes()
        al = self.get_alive().reshape(-1)
        g = g.reshape(-1, NSLOT, GENE_BYTES)
        a = a.reshape(-1)
        cells, slots = np.nonzero(((a[:, None] >> np.arange(NSLOT)) & 1).astype(bool)
                                  & (al[:, None] == 1))
        return cells, slots, g[cells, slots]

    def gene_type_name(self, b):
        d = self.decode_gene(b)
        cat = '' if d['C'] is None else f" [{d['C']}^{d['m']}]"
        return f"{d['A']}->{d['B']}{cat}"

    def k_idx_by_type(self):
        """{gene type 'A->B [C^m]': (count, mean k_idx)} over active genes."""
        _, _, gb = self.active_genes()
        out = {}
        if len(gb) == 0:
            return out
        keys = gb[:, :G_K].astype(np.int64) @ np.array([1 << 24, 1 << 16, 1 << 8, 1])
        for key in np.unique(keys):
            sel = keys == key
            out[self.gene_type_name(gb[sel][0])] = (int(sel.sum()),
                                                    float(gb[sel, G_K].mean()))
        return out

    def n_distinct_genesets(self):
        """Distinct active gene sets (as sorted multisets of 5-byte genes)
        among alive cells."""
        g, a = self.get_rgenes()
        al = self.get_alive().reshape(-1) == 1
        g = g.reshape(-1, NSLOT, GENE_BYTES)[al]
        a = a.reshape(-1)[al]
        seen = set()
        for gi, ai in zip(g, a):
            seen.add(tuple(sorted(gi[s].tobytes() for s in range(NSLOT) if ai >> s & 1)))
        return len(seen)

    # ── Display ──────────────────────────────────────────────────────

    def colorize(self, pixels=None, colormode=0, scale=None):
        """Fill an int32 ARGB buffer (N*N) like EvoCA.colorize. scale sets
        the amount mapped to full intensity (default: mbar of the state)."""
        if pixels is None:
            pixels = np.zeros(self._N * self._N, np.int32)
        if scale is None:
            scale = self._state_params.get('mbar', 1.0) or 1.0
        self._lib.chemca_set_color_scale(float(scale))
        self._lib.chemca_colorize(pixels.ctypes.data_as(_i32p), int(colormode))
        return pixels

    def rgb(self, colormode=0, scale=None):
        """(N, N, 3) uint8 frame, headless."""
        a = self.colorize(None, colormode, scale).astype(np.uint32).reshape(self._N, self._N)
        return np.stack([(a >> 16) & 255, (a >> 8) & 255, a & 255], -1).astype(np.uint8)

    def save_frames(self, path, modes=None, scale=None, dpi=100):
        """Render colormodes side by side into a PNG (matplotlib Agg). Returns path."""
        import matplotlib
        matplotlib.use('Agg', force=False)
        import matplotlib.pyplot as plt
        modes = range(self.colormode_count()) if modes is None else modes
        modes = list(modes)
        fig, axs = plt.subplots(1, len(modes), figsize=(3 * len(modes), 3.2))
        for ax, m in zip(np.atleast_1d(axs), modes):
            ax.imshow(self.rgb(m, scale), interpolation='nearest')
            ax.set_title(self.COLORMODE_NAMES[m], fontsize=9)
            ax.axis('off')
        fig.suptitle(f"chemca {self.preset}  N={self._N}  step {self.get_step()}  "
                     f"alive {self.population() / self._N ** 2:.2f}", fontsize=10)
        fig.tight_layout()
        fig.savefig(path, dpi=dpi)
        plt.close(fig)
        return path

    # ── Recipes ──────────────────────────────────────────────────────

    def recipe(self, descriptor='', colormode=0):
        """The recipe dict (model 'chemca')."""
        return {
            'version': 1,
            'model': 'chemca',
            'created': datetime.now().isoformat(timespec='seconds'),
            'descriptor': descriptor,
            'N': self._N,
            'metaparams_init': dict(self._init_params),
            'metaparams_final': dict(self._params),
            'initialization': dict(self._state_params),
            'light_field': self._light_field is not None,
            'display': {'colormode': colormode},
        }

    def export_recipe(self, descriptor, colormode=0, workstream=None, runs_dir=None):
        """Write <runs_dir>/<date>_<descriptor>.evoca (as EvoCA.export_recipe).
        A per-site light field is not serialized (the recipe says so)."""
        import evoca_paths
        if runs_dir is None:
            runs_dir = evoca_paths.runs_dir(workstream)
        os.makedirs(runs_dir, exist_ok=True)
        safe = descriptor.replace(' ', '_')
        path = os.path.join(runs_dir, f"{datetime.now().strftime('%Y-%m-%d')}_{safe}.evoca")
        with open(path, 'w') as f:
            json.dump(self.recipe(descriptor, colormode), f, indent=2)
        return path


def load_recipe(data, recipe='init', N=None, lib_path=None):
    """A chemca recipe (dict or path) -> (sim, display_kwargs), state applied."""
    if not isinstance(data, dict):
        with open(data) as f:
            data = json.load(f)
    if sim_model.model_of(data) != 'chemca':
        raise ValueError("not a chemca recipe")
    mp = dict(data['metaparams_init'] if recipe == 'init' else data['metaparams_final'])
    sim = ChemCA(lib_path=lib_path)
    sim.init(int(N) if N is not None else int(data['N']), **mp)
    init = dict(data.get('initialization', {}))
    if init.get('alive') == 'array':
        raise ValueError("recipe was made from an alive array, which is not serialized")
    sim.state(**init)
    disp = data.get('display', {})
    kw = {'colormode': disp['colormode']} if disp.get('colormode') else {}
    return sim, kw


# ── Headless construction and runs ───────────────────────────────────

def split_params(params):
    """(init kwargs, state kwargs) from a flat params dict that may carry
    'model', metaparameters and STATE_DEFAULTS keys (e.g. mbar). The
    caller has chosen ChemCA, so a missing 'model' is fine here (the
    'evoca' default applies at the dispatch points); a different one raises."""
    model, p = sim_model.split_model({sim_model.MODEL_KEY: 'chemca', **params})
    if model != 'chemca':
        raise ValueError(f"not a chemca params dict (model={model!r})")
    p.pop('N', None)
    state = {k: p.pop(k) for k in list(p) if k in STATE_DEFAULTS}
    return p, state


def build_chemca(params, N=128, seed=0, init='halfplane', lib_path=None):
    """Construct and initialize one ChemCA headless (the counterpart of
    evoca_explore.build_sim, which dispatches here on model='chemca').
    init: 'halfplane' | 'fraction' (alive pattern); seed reseeds the C PRNG
    and the numpy draws."""
    ip, st = split_params(params)
    sim = ChemCA(lib_path=lib_path)
    sim.init(N, **ip)
    st.setdefault('alive', init)
    st.setdefault('seed', seed)
    sim.state(**st)
    return sim


def run_chemca(params, n_steps=3000, N=128, seed=0, sample_every=100,
               init='halfplane', kymo_steps=1500, comp_window=100):
    """Run one ChemCA configuration headless; return a summary dict.

    Checked every step: the mass identity (closed: total mass constant;
    open: change = in - out) and the energy identity
    d(sum e) = yield - cost - death - split. Kymographs of the middle row
    (alive and E of compound 0) over the last kymo_steps steps give
    evoca_waves.site_period. Over the last comp_window steps the stores are
    averaged per cell (cells alive throughout), and composition_clusters() of the
    final stores, with those averages appended, is the isologous-
    diversification readout (meaningful with genes frozen:
    mu_rgene = mu_ractive = 0)."""
    import evoca_waves
    sim = build_chemca(params, N=N, seed=seed, init=init)
    try:
        M0 = sim.mass()['total']
        e_prev = sim.e_total()
        mass_err = energy_err = 0.0
        alive_t, samples = [], []
        kym_a, kym_E = [], []
        K = len(sim.compounds)
        st_sum = np.zeros((K, N, N))
        alive_win = np.ones((N, N), bool)
        tot = dict(births=0, deaths=0, tox_deaths=0, energy_deaths=0, evictions=0)
        for t in range(1, n_steps + 1):
            sim.step()
            m = sim.mass()['total']
            fl = sim.mass_flux()
            mass_err = max(mass_err, abs(m - M0 - (fl[:, 0].sum() - fl[:, 1].sum())) / M0)
            b = sim.energy_budget(cum=False)
            e_now = sim.e_total()
            de = e_now - e_prev
            energy_err = max(energy_err, abs(de - (b['yield'] - b['cost'] - b['death'] - b['split']))
                             / max(1.0, abs(e_now)))
            e_prev = e_now
            for k, v in sim.events_last().items():
                tot[k] += v
            pop = sim.population()
            alive_t.append(pop / N ** 2)
            if t > n_steps - kymo_steps:
                kym_a.append(sim.get_alive()[N // 2].astype(float))
                kym_E.append(sim.get_E(0)[N // 2].astype(float))
            if t > n_steps - comp_window:
                alive_win &= sim.get_alive() == 1
                for c in range(K):
                    st_sum[c] += sim.get_store(c)
            if t % sample_every == 0:
                samples.append(dict(t=t, alive=pop / N ** 2, k_idx=sim.k_idx_by_type()))
            if pop == 0:
                break
        alive_t = np.array(alive_t)
        tail = alive_t[len(alive_t) // 2:]
        pa, pa_h = evoca_waves.site_period(np.array(kym_a)) if len(kym_a) > 10 else (float('nan'), 0.0)
        pe, pe_h = evoca_waves.site_period(np.array(kym_E)) if len(kym_E) > 10 else (float('nan'), 0.0)
        comp = composition_clusters(sim, time_avg=st_sum if len(alive_t) == n_steps else None,
                                    time_avg_mask=alive_win, seed=seed)
        out = dict(
            model='chemca', preset=sim.preset, N=N, seed=seed, n_steps=n_steps,
            steps_run=len(alive_t),
            alive_mean=float(tail.mean()), alive_std=float(tail.std()),
            alive_final=float(alive_t[-1]), extinct=bool(alive_t[-1] == 0),
            mass_err_max=float(mass_err), energy_err_max=float(energy_err),
            energy_budget=sim.energy_budget(cum=True),
            k_idx_by_type=sim.k_idx_by_type(),
            n_genesets_final=sim.n_distinct_genesets(),
            site_period_alive=pa, site_period_alive_peak=pa_h,
            site_period_E=pe, site_period_E_peak=pe_h,
            n_comp_clusters=comp['n_clusters'], comp_silhouette=comp['silhouette'],
            comp_fractions=comp['fractions'], comp_n_cells=comp['n_cells'],
            samples=samples, alive_series=alive_t,
            params=sim.params(), state=sim.state_params,
            **{f'{k}_total': v for k, v in tot.items()},
        )
        return out
    finally:
        sim.free()


# ── Isologous diversification readout ────────────────────────────────

def _kmeans(X, k, rng, n_init=4, iters=50):
    """Plain k-means++ / Lloyd; returns (labels, centers, inertia)."""
    best = None
    for _ in range(n_init):
        C = [X[rng.integers(len(X))]]
        for _j in range(1, k):
            d2 = np.min(((X[:, None, :] - np.array(C)[None]) ** 2).sum(-1), axis=1)
            tot = d2.sum()
            C.append(X[rng.choice(len(X), p=d2 / tot)] if tot > 0 else X[rng.integers(len(X))])
        C = np.array(C)
        for _it in range(iters):
            lab = np.argmin(((X[:, None, :] - C[None]) ** 2).sum(-1), axis=1)
            newC = np.array([X[lab == j].mean(0) if np.any(lab == j) else C[j] for j in range(k)])
            if np.allclose(newC, C):
                break
            C = newC
        lab = np.argmin(((X[:, None, :] - C[None]) ** 2).sum(-1), axis=1)
        inertia = float(((X - C[lab]) ** 2).sum())
        if best is None or inertia < best[2]:
            best = (lab, C, inertia)
    return best


def _silhouette(X, lab):
    """Mean silhouette (numpy, O(n^2) memory; callers subsample)."""
    D = np.sqrt(((X[:, None, :] - X[None]) ** 2).sum(-1))
    ks = np.unique(lab)
    if len(ks) < 2:
        return 0.0
    M = np.stack([D[:, lab == j].sum(1) for j in ks], 1)      # n x k sums
    cnt = np.array([(lab == j).sum() for j in ks])
    own = np.searchsorted(ks, lab)
    n_own = cnt[own] - 1
    a = np.where(n_own > 0, M[np.arange(len(X)), own] / np.maximum(n_own, 1), 0.0)
    other = M / cnt[None]
    other[np.arange(len(X)), own] = np.inf
    b = other.min(1)
    s = np.where(n_own > 0, (b - a) / np.maximum(np.maximum(a, b), 1e-300), 0.0)
    return float(s.mean())


def composition_clusters(sim, time_avg=None, time_avg_mask=None, k_max=6,
                         max_n=1500, seed=0, sil_min=0.5, spread_min=1e-3):
    """Do (genetically identical) cells fall into distinct compositional
    types? Isologous diversification (Furusawa & Kaneko 1998/2000).

    Features: each alive cell's store vector normalized by its total; if
    time_avg (K x N x N store sums over a window) is given, the cell's
    window-averaged normalized composition is appended (cells outside
    time_avg_mask are dropped). k-means for k = 2..k_max on a subsample of
    at most max_n cells; k is chosen by mean silhouette. One cluster is
    returned when the composition spread (largest per-feature std) is below
    spread_min, or when the best silhouette is below sil_min.

    Returns dict(n_clusters, silhouette (best, k >= 2), fractions (sorted,
    over all used cells), centers, silhouettes {k: s}, n_cells)."""
    alive = sim.get_alive() == 1
    if time_avg is not None and time_avg_mask is not None:
        alive = alive & time_avg_mask
    S = np.stack([sim.get_store(c) for c in sim.compounds])[:, alive].T.astype(float)
    tot = S.sum(1)
    keep = tot > 0
    X = S[keep] / tot[keep, None]
    if time_avg is not None:
        T = np.asarray(time_avg, float)[:, alive].T[keep]
        tt = T.sum(1)
        X = np.hstack([X, T / np.where(tt > 0, tt, 1.0)[:, None]])
    out = dict(n_clusters=1, silhouette=0.0, fractions=[1.0], centers=None,
               silhouettes={}, n_cells=int(len(X)))
    if len(X) < 4:
        out['fractions'] = [1.0] if len(X) else []
        return out
    out['centers'] = X.mean(0)
    if X.std(0).max() < spread_min:
        return out
    rng = np.random.default_rng(seed)
    Xs = X[rng.choice(len(X), max_n, replace=False)] if len(X) > max_n else X
    best = None
    for k in range(2, min(k_max, len(Xs) - 1) + 1):
        lab, C, _ = _kmeans(Xs, k, rng)
        sil = _silhouette(Xs, lab)
        out['silhouettes'][k] = sil
        if best is None or sil > best[0]:
            best = (sil, k, C)
    if best is None:
        return out
    out['silhouette'] = best[0]
    if best[0] < sil_min:
        return out
    C = best[2]
    lab = np.argmin(((X[:, None, :] - C[None]) ** 2).sum(-1), axis=1)
    fr = np.bincount(lab, minlength=len(C)) / len(X)
    order = np.argsort(-fr)
    out.update(n_clusters=int(best[1]), fractions=[float(f) for f in fr[order]],
               centers=C[order])
    return out
