"""
evoca_py.py — Python ctypes wrapper for the EvoCA C library.

LUT indexing
------------
The LUT conditions on (v_x, n1, n2, n3):
  v_x ∈ {0,1}  — current cell state
  n1  ∈ {0..4} — active cells at distance 1   (4 orthogonal Moore nbrs)
  n2  ∈ {0..4} — active cells at distance √2  (4 diagonal Moore nbrs)
  n3  ∈ {0..4} — active cells at distance 2   (4 axis-aligned)

Flat bit index:  v_x*125 + n1*25 + n2*5 + n3
Total bits:  2*5*5*5 = 250  →  32 bytes bit-packed per cell.

GoL is exactly encodable: it conditions on n1+n2 and ignores n3.
"""

import ctypes
import os

import numpy as np

# ── LUT constants ──────────────────────────────────────────────────────
LUT_BITS  = 250   # 2*5*5*5
LUT_BYTES = 32    # ceil(250/8)

# Ring strides for flat index
_S = (125, 25, 5, 1)   # (v_x, n1, n2, n3)


def lut_bit_index(v_x, n1, n2, n3):
    return v_x*125 + n1*25 + n2*5 + n3


# D4-symmetric orbit map for the 5×5 fiducial pattern.
# 6 independent orbits (bits 0-5 of egenome).
ORBIT_MAP = np.array([
    [4, 5, 2, 5, 4],
    [5, 3, 1, 3, 5],
    [2, 1, 0, 1, 2],
    [5, 3, 1, 3, 5],
    [4, 5, 2, 5, 4],
], dtype=np.uint8)


def egenome_to_pattern(eg):
    """Expand a 6-bit egenome value to a 5×5 binary numpy array."""
    return ((eg >> ORBIT_MAP) & 1).astype(np.uint8)


def _find_lib():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.join(here, "..")
    for name in ("libevoca.dylib", "libevoca.so"):
        p = os.path.realpath(os.path.join(root, "C", name))
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        "libevoca not found. Run `make` (or the gcc command in CLAUDE.md) first.")


class EvoCA:
    """Thin ctypes wrapper around the EvoCA shared library."""

    def __init__(self, lib_path=None):
        self._lib        = ctypes.CDLL(lib_path or _find_lib())
        self._N          = 0
        self.food_inc    = 0.0
        self.m_scale     = 0.0
        self.gdiff       = 0
        self.mu_lut      = 0.0
        self.mu_egenome   = 0.0
        self.tax         = 0.0
        self.restricted_mu = 0
        self._n_ent      = 2
        self.egenome      = 0
        self._state_params = {}
        self._setup_signatures()

    def _setup_signatures(self):
        L = self._lib
        L.evoca_init.argtypes           = [ctypes.c_int, ctypes.c_float,
                                            ctypes.c_float]
        L.evoca_init.restype            = None
        L.evoca_free.argtypes           = []
        L.evoca_free.restype            = None
        L.evoca_set_food_inc.argtypes   = [ctypes.c_float]
        L.evoca_set_food_inc.restype    = None
        L.evoca_set_m_scale.argtypes    = [ctypes.c_float]
        L.evoca_set_m_scale.restype     = None
        L.evoca_set_gdiff.argtypes      = [ctypes.c_int]
        L.evoca_set_gdiff.restype       = None
        L.evoca_get_gdiff.argtypes      = []
        L.evoca_get_gdiff.restype       = ctypes.c_int
        L.evoca_set_mu_lut.argtypes     = [ctypes.c_float]
        L.evoca_set_mu_lut.restype      = None
        L.evoca_set_mu_egenome.argtypes  = [ctypes.c_float]
        L.evoca_set_mu_egenome.restype   = None
        L.evoca_get_mu_lut.argtypes     = []
        L.evoca_get_mu_lut.restype      = ctypes.c_float
        L.evoca_get_mu_egenome.argtypes  = []
        L.evoca_get_mu_egenome.restype   = ctypes.c_float
        L.evoca_set_tax.argtypes        = [ctypes.c_float]
        L.evoca_set_tax.restype         = None
        L.evoca_get_tax.argtypes        = []
        L.evoca_get_tax.restype         = ctypes.c_float
        L.evoca_set_restricted_mu.argtypes = [ctypes.c_int]
        L.evoca_set_restricted_mu.restype  = None
        L.evoca_get_restricted_mu.argtypes = []
        L.evoca_get_restricted_mu.restype  = ctypes.c_int
        L.evoca_set_diag.argtypes          = [ctypes.c_int]
        L.evoca_set_diag.restype           = None
        L.evoca_get_diag.argtypes          = []
        L.evoca_get_diag.restype           = ctypes.c_int
        L.evoca_set_v_all.argtypes      = [ctypes.POINTER(ctypes.c_uint8),
                                            ctypes.c_int]
        L.evoca_set_v_all.restype       = None
        L.evoca_set_lut_all.argtypes    = [ctypes.POINTER(ctypes.c_uint8)]
        L.evoca_set_lut_all.restype     = None
        L.evoca_set_lut.argtypes        = [ctypes.c_int,
                                            ctypes.POINTER(ctypes.c_uint8)]
        L.evoca_set_lut.restype         = None
        L.evoca_set_egenome_all.argtypes = [ctypes.c_uint8]
        L.evoca_set_egenome_all.restype  = None
        L.evoca_set_f_all.argtypes      = [ctypes.c_float]
        L.evoca_set_f_all.restype       = None
        L.evoca_set_F_all.argtypes      = [ctypes.c_float]
        L.evoca_set_F_all.restype       = None
        L.evoca_set_env_mask.argtypes   = [ctypes.POINTER(ctypes.c_uint8)]
        L.evoca_set_env_mask.restype    = None
        L.evoca_get_env_mask.argtypes   = []
        L.evoca_get_env_mask.restype    = ctypes.POINTER(ctypes.c_uint8)
        # Alive data plane
        L.evoca_set_alive.argtypes      = [ctypes.POINTER(ctypes.c_uint8)]
        L.evoca_set_alive.restype       = None
        L.evoca_get_alive.argtypes      = []
        L.evoca_get_alive.restype       = ctypes.POINTER(ctypes.c_uint8)
        L.evoca_set_alive_all.argtypes  = []
        L.evoca_set_alive_all.restype   = None
        L.evoca_set_alive_fraction.argtypes = [ctypes.c_float]
        L.evoca_set_alive_fraction.restype  = None
        L.evoca_set_alive_patch.argtypes    = [ctypes.c_int]
        L.evoca_set_alive_patch.restype     = None
        L.evoca_set_alive_halfplane.argtypes = [ctypes.c_int]
        L.evoca_set_alive_halfplane.restype  = None
        L.evoca_step.argtypes           = []
        L.evoca_step.restype            = None
        L.evoca_colorize.argtypes       = [ctypes.POINTER(ctypes.c_int32),
                                            ctypes.c_int]
        L.evoca_colorize.restype        = None
        L.evoca_get_v.argtypes          = []
        L.evoca_get_v.restype           = ctypes.POINTER(ctypes.c_uint8)
        L.evoca_get_F.argtypes          = []
        L.evoca_get_F.restype           = ctypes.POINTER(ctypes.c_float)
        L.evoca_get_f.argtypes          = []
        L.evoca_get_f.restype           = ctypes.POINTER(ctypes.c_float)
        L.evoca_get_egenome.argtypes     = []
        L.evoca_get_egenome.restype      = ctypes.POINTER(ctypes.c_uint8)
        L.evoca_get_lut.argtypes        = []
        L.evoca_get_lut.restype         = ctypes.POINTER(ctypes.c_uint8)
        L.evoca_get_births.argtypes     = []
        L.evoca_get_births.restype      = ctypes.POINTER(ctypes.c_uint8)
        L.evoca_get_N.argtypes          = []
        L.evoca_get_N.restype           = ctypes.c_int
        L.evoca_get_cell_px.argtypes    = []
        L.evoca_get_cell_px.restype     = ctypes.c_int
        L.evoca_set_act_ymax.argtypes   = [ctypes.c_int]
        L.evoca_set_act_ymax.restype    = None
        L.evoca_get_act_ymax.argtypes   = []
        L.evoca_get_act_ymax.restype    = ctypes.c_int
        L.evoca_get_repro_age_hist.argtypes = []
        L.evoca_get_repro_age_hist.restype  = ctypes.POINTER(ctypes.c_uint32)
        L.evoca_get_repro_age_max.argtypes  = []
        L.evoca_get_repro_age_max.restype   = ctypes.c_int
        L.evoca_get_step.argtypes           = []
        L.evoca_get_step.restype            = ctypes.c_uint32
        L.evoca_set_repro_age_t0.argtypes   = [ctypes.c_uint32]
        L.evoca_set_repro_age_t0.restype    = None
        L.evoca_get_repro_age_t0.argtypes   = []
        L.evoca_get_repro_age_t0.restype    = ctypes.c_uint32
        L.evoca_reset_repro_age_hist.argtypes = []
        L.evoca_reset_repro_age_hist.restype  = None
        L.evoca_activity_update.argtypes  = []
        L.evoca_activity_update.restype   = None
        L.evoca_activity_render_col.argtypes = [
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int]
        L.evoca_activity_render_col.restype  = None
        L.evoca_activity_get.argtypes   = [
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_uint64),
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_int32),
            ctypes.c_int]
        L.evoca_activity_get.restype    = ctypes.c_int
        # Egenome activity
        L.evoca_eg_activity_update.argtypes  = []
        L.evoca_eg_activity_update.restype   = None
        L.evoca_eg_activity_render_col.argtypes = [
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int]
        L.evoca_eg_activity_render_col.restype  = None
        L.evoca_eg_activity_get.argtypes = [
            ctypes.POINTER(ctypes.c_uint64),
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_int32)]
        L.evoca_eg_activity_get.restype  = ctypes.c_int
        L.evoca_set_eg_act_ymax.argtypes = [ctypes.c_int]
        L.evoca_set_eg_act_ymax.restype  = None
        L.evoca_get_eg_act_ymax.argtypes = []
        L.evoca_get_eg_act_ymax.restype  = ctypes.c_int
        # LUT complexity
        L.evoca_lut_complexity_counts.argtypes = [
            ctypes.POINTER(ctypes.c_uint32)]
        L.evoca_lut_complexity_counts.restype  = None
        L.evoca_lut_complexity_render_col.argtypes = [
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int]
        L.evoca_lut_complexity_render_col.restype  = None
        # Egenome population banded chart
        L.evoca_eg_pop_render_col.argtypes = [
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int]
        L.evoca_eg_pop_render_col.restype  = None
        # Local-pattern activity & entropy
        L.evoca_set_n_ent.argtypes = [ctypes.c_int]
        L.evoca_set_n_ent.restype  = None
        L.evoca_get_n_ent.argtypes = []
        L.evoca_get_n_ent.restype  = ctypes.c_int
        L.evoca_pat_update.argtypes = []
        L.evoca_pat_update.restype  = None
        L.evoca_get_entropy.argtypes = []
        L.evoca_get_entropy.restype  = ctypes.c_float
        L.evoca_pat_activity_render_col.argtypes = [
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int]
        L.evoca_pat_activity_render_col.restype  = None
        L.evoca_set_pat_act_ymax.argtypes = [ctypes.c_int]
        L.evoca_set_pat_act_ymax.restype  = None
        L.evoca_get_pat_act_ymax.argtypes = []
        L.evoca_get_pat_act_ymax.restype  = ctypes.c_int

    # ── Lifecycle ──────────────────────────────────────────────────────

    def init(self, N, food_inc=0.0, m_scale=1.0, gdiff=0,
             mu_lut=0.0, mu_egenome=0.0, tax=0.0, restricted_mu=0, n_ent=2):
        stop = getattr(self, '_stop_display', None)
        if stop is not None:
            stop()
            self._stop_display = None
        self._N         = N
        self.food_inc   = float(food_inc)
        self.m_scale    = float(m_scale)
        self.gdiff      = int(gdiff)
        self.mu_lut     = float(mu_lut)
        self.mu_egenome  = float(mu_egenome)
        self.tax        = float(tax)
        self.restricted_mu = int(restricted_mu)
        self._n_ent     = int(n_ent)
        self._lib.evoca_init(N, self.food_inc, self.m_scale)
        self._lib.evoca_set_gdiff(self.gdiff)
        self._lib.evoca_set_mu_lut(self.mu_lut)
        self._lib.evoca_set_mu_egenome(self.mu_egenome)
        self._lib.evoca_set_tax(self.tax)
        self._lib.evoca_set_restricted_mu(self.restricted_mu)
        self._lib.evoca_set_n_ent(self._n_ent)
        self._state_params = {}
        self._init_metaparams = {
            'food_inc': self.food_inc,
            'm_scale': self.m_scale,
            'gdiff': self.gdiff,
            'mu_lut': self.mu_lut,
            'mu_egenome': self.mu_egenome,
            'tax': self.tax,
            'restricted_mu': self.restricted_mu,
        }

    def free(self):
        stop = getattr(self, '_stop_display', None)
        if stop is not None:
            stop()
            self._stop_display = None
        self._lib.evoca_free()
        self._N = 0

    def __del__(self):
        try:
            if self._N:
                self.free()
        except Exception:
            pass

    # ── Metaparam setters ──────────────────────────────────────────────

    def update_food_inc(self, f):
        self.food_inc = float(f)
        self._lib.evoca_set_food_inc(self.food_inc)

    def update_m_scale(self, m):
        self.m_scale = float(m)
        self._lib.evoca_set_m_scale(self.m_scale)

    def update_gdiff(self, d):
        self.gdiff = int(d)
        self._lib.evoca_set_gdiff(self.gdiff)

    def update_mu_lut(self, m):
        self.mu_lut = float(m)
        self._lib.evoca_set_mu_lut(self.mu_lut)

    def update_mu_egenome(self, m):
        self.mu_egenome = float(m)
        self._lib.evoca_set_mu_egenome(self.mu_egenome)

    def update_tax(self, t):
        self.tax = float(t)
        self._lib.evoca_set_tax(self.tax)

    def update_restricted_mu(self, r):
        self.restricted_mu = int(r)
        self._lib.evoca_set_restricted_mu(self.restricted_mu)

    def update_act_ymax(self, y):
        self._lib.evoca_set_act_ymax(int(y))

    def update_eg_act_ymax(self, y):
        self._lib.evoca_set_eg_act_ymax(int(y))

    def update_pat_act_ymax(self, y):
        self._lib.evoca_set_pat_act_ymax(int(y))

    # ── Params export ─────────────────────────────────────────────────

    _DEFAULTS = dict(food_inc=0.0, m_scale=1.0,
                      gdiff=0, mu_lut=0.0, mu_egenome=0.0, tax=0.0,
                      restricted_mu=0)

    def params(self):
        """Return current metaparameters as a dict suitable for init(**d)."""
        return dict(N=self._N, food_inc=self.food_inc, m_scale=self.m_scale,
                    gdiff=self.gdiff, mu_lut=self.mu_lut,
                    mu_egenome=self.mu_egenome, tax=self.tax,
                    restricted_mu=self.restricted_mu)

    def params_str(self):
        """Return a copy-pasteable sim.init(...) call with defaults annotated."""
        p = self.params()
        lines = ["sim.init("]
        for k, v in p.items():
            val = repr(v)
            if k in self._DEFAULTS and v != self._DEFAULTS[k]:
                lines.append(f"    {k}={val},   # default: {self._DEFAULTS[k]!r}")
            else:
                lines.append(f"    {k}={val},")
        lines.append(")")
        return "\n".join(lines)

    # ── State initialization ─────────────────────────────────────────

    @property
    def state_params(self):
        """Return a copy of the current state-initialization parameters."""
        return dict(self._state_params)

    def state(self, lut='gol', lut_n_init=3, egenome='uniform', egenome_value=0,
              v_density=0.5, f_init=0.0, F='uniform', F_init=0.0, F_range=None,
              alive='all', alive_fraction=0.5, alive_radius=64, alive_axis=0):
        """Initialize lattice state from parameters.

        See available_state_init() for parameter descriptions.
        """
        self._state_params = {}

        # LUT
        if lut == 'random':
            self.set_lut_random(n_init=lut_n_init)
        else:
            self.set_lut_all(make_gol_lut())

        # Egenome
        if egenome == 'random':
            self.set_egenome_random()
        else:
            self.set_egenome_all(egenome_value)

        # v(x) — random binary with given density
        self.set_v(np.random.default_rng().integers(
            0, 2, (self._N, self._N), dtype=np.uint8), density=v_density)

        # Private food
        self.set_f_all(f_init)

        # Environmental food
        if F == 'random':
            lo, hi = F_range or [0.0, 1.0]
            self.set_F_random(lo, hi)
        else:
            self.set_F_all(F_init)

        # Alive mask (must come after LUT/egenome — dead cells' data is zeroed)
        if alive == 'fraction':
            self.set_alive_fraction(alive_fraction)
        elif alive == 'patch':
            self.set_alive_patch(alive_radius)
        elif alive == 'halfplane':
            self.set_alive_halfplane(alive_axis)
        else:
            self.set_alive_all()

    # ── Grid setters ──────────────────────────────────────────────────

    def set_v(self, v_array, density=None):
        arr = np.ascontiguousarray(v_array.ravel(), dtype=np.uint8)
        self._lib.evoca_set_v_all(
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)), len(arr))
        self._state_params['v_density'] = float(density) if density is not None else float(arr.mean())

    def set_lut_all(self, lut_bytes):
        arr = np.ascontiguousarray(lut_bytes, dtype=np.uint8)
        assert len(arr) == LUT_BYTES, f"need {LUT_BYTES} bytes, got {len(arr)}"
        self._lib.evoca_set_lut_all(
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)))
        self._state_params['lut'] = 'gol'
        self._state_params.pop('lut_n_init', None)

    def set_lut(self, idx, lut_bytes):
        arr = np.ascontiguousarray(lut_bytes, dtype=np.uint8)
        assert len(arr) == LUT_BYTES
        self._lib.evoca_set_lut(
            int(idx), arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)))

    def set_egenome_all(self, eg):
        self.egenome = int(eg) & 0x3F
        self._lib.evoca_set_egenome_all(self.egenome)
        self._state_params['egenome'] = 'uniform'
        self._state_params['egenome_value'] = self.egenome

    def set_egenome_random(self):
        """Set each cell's egenome to a random value in [0, 63].
        No wild-type: all 64 egenomes get distinct hash-based colors."""
        self._lib.evoca_set_egenome_all(0xFF)   # wt=0xFF → no match → all colored
        ptr = self._lib.evoca_get_egenome()
        arr = np.ctypeslib.as_array(ptr, shape=(self._N * self._N,))
        arr[:] = np.random.randint(0, 64, self._N * self._N, dtype=np.uint8)
        self._state_params['egenome'] = 'random'
        self._state_params.pop('egenome_value', None)

    def set_lut_random(self, n_init=3):
        """Set each cell's LUT to an independent random rule.

        n_init: number of rings the rule conditions on (1, 2, or 3).
          1: depends on (v_x, n1) — 10 independent bits per LUT
          2: depends on (v_x, n1, n2) — 50 independent bits
          3: depends on (v_x, n1, n2, n3) — all 250 bits independent
        """
        N2 = self._N * self._N

        if n_init >= 3:
            bits = np.random.randint(0, 2, (N2, 2, 5, 5, 5), dtype=np.uint8)
        elif n_init == 2:
            core = np.random.randint(0, 2, (N2, 2, 5, 5, 1), dtype=np.uint8)
            bits = np.tile(core, (1, 1, 1, 1, 5))
        else:  # n_init == 1
            core = np.random.randint(0, 2, (N2, 2, 5, 1, 1), dtype=np.uint8)
            bits = np.tile(core, (1, 1, 1, 5, 5))

        flat = bits.reshape(N2, LUT_BITS)

        # Bit-pack to (N2, 32) bytes
        padded = np.zeros((N2, 256), dtype=np.uint8)
        padded[:, :LUT_BITS] = flat
        grouped = padded.reshape(N2, 32, 8)
        packed = np.zeros((N2, 32), dtype=np.uint8)
        for b in range(8):
            packed |= grouped[:, :, b] << b

        packed = np.ascontiguousarray(packed)
        for i in range(N2):
            self.set_lut(i, packed[i])
        self._state_params['lut'] = 'random'
        self._state_params['lut_n_init'] = n_init

    def set_f_all(self, f):
        self._lib.evoca_set_f_all(float(f))
        self._state_params['f_init'] = float(f)

    def set_F_all(self, F):
        self._lib.evoca_set_F_all(float(F))
        self._state_params['F'] = 'uniform'
        self._state_params['F_init'] = float(F)
        self._state_params.pop('F_range', None)

    def set_F_random(self, lo=0.0, hi=1.0):
        """Set env food F(x) to uniform random values in [lo, hi]."""
        ptr = self._lib.evoca_get_F()
        arr = np.ctypeslib.as_array(ptr, shape=(self._N * self._N,))
        arr[:] = np.random.uniform(lo, hi, self._N * self._N).astype(np.float32)
        self._state_params['F'] = 'random'
        self._state_params.pop('F_init', None)
        self._state_params['F_range'] = [float(lo), float(hi)]

    def set_env_mask(self, mask):
        """Set the environment food-regeneration mask.
        mask: (N,N) or flat uint8 array. 1=regenerate, 0=no regen."""
        arr = np.ascontiguousarray(np.asarray(mask, dtype=np.uint8).ravel())
        assert len(arr) == self._N * self._N
        self._lib.evoca_set_env_mask(
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)))

    def get_env_mask(self):
        """Return (N, N) uint8 array of environment mask (copy)."""
        ptr = self._lib.evoca_get_env_mask()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    # ── Alive data plane ──────────────────────────────────────────────

    def set_alive(self, arr):
        """Set alive array. arr: (N,N) or flat uint8. Dead cells' data is zeroed."""
        a = np.ascontiguousarray(np.asarray(arr, dtype=np.uint8).ravel())
        assert len(a) == self._N * self._N
        self._lib.evoca_set_alive(
            a.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)))

    def get_alive(self):
        """Return (N, N) uint8 copy of alive array."""
        ptr = self._lib.evoca_get_alive()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def set_alive_all(self):
        """Set all cells alive."""
        self._lib.evoca_set_alive_all()
        self._state_params['alive'] = 'all'
        for k in ('alive_fraction', 'alive_radius', 'alive_axis'):
            self._state_params.pop(k, None)

    def set_alive_fraction(self, frac):
        """Set a random fraction of cells alive; dead cells' data is zeroed."""
        self._lib.evoca_set_alive_fraction(float(frac))
        self._state_params['alive'] = 'fraction'
        self._state_params['alive_fraction'] = float(frac)
        for k in ('alive_radius', 'alive_axis'):
            self._state_params.pop(k, None)

    def set_alive_patch(self, radius):
        """Set a square patch of radius cells at center alive; rest dead."""
        self._lib.evoca_set_alive_patch(int(radius))
        self._state_params['alive'] = 'patch'
        self._state_params['alive_radius'] = int(radius)
        for k in ('alive_fraction', 'alive_axis'):
            self._state_params.pop(k, None)

    def set_alive_halfplane(self, axis=0):
        """Set half the grid alive. axis=0: left half, axis=1: top half."""
        self._lib.evoca_set_alive_halfplane(int(axis))
        self._state_params['alive'] = 'halfplane'
        self._state_params['alive_axis'] = int(axis)
        for k in ('alive_fraction', 'alive_radius'):
            self._state_params.pop(k, None)

    # ── Step and colorize ─────────────────────────────────────────────

    def step(self):
        self._lib.evoca_step()

    def colorize(self, pixels, colormode=0):
        """Fill a (N*N,) int32 numpy array with ARGB values in-place."""
        arr = np.ascontiguousarray(pixels, dtype=np.int32)
        self._lib.evoca_colorize(
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            int(colormode))
        return arr

    # ── Getters ───────────────────────────────────────────────────────

    def get_v(self):
        ptr = self._lib.evoca_get_v()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def get_F(self):
        ptr = self._lib.evoca_get_F()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def get_f(self):
        ptr = self._lib.evoca_get_f()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def get_egenome(self):
        """Return (N, N) uint8 array of fiducial genomes (6-bit values)."""
        ptr = self._lib.evoca_get_egenome()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def get_births(self):
        """Return (N, N) uint8 array of birth events from last step."""
        ptr = self._lib.evoca_get_births()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def get_lut(self, idx):
        """Return LUT_BYTES-length uint8 array for cell at flat index idx."""
        ptr = self._lib.evoca_get_lut()
        arr = np.ctypeslib.as_array(ptr, shape=(self._N * self._N * LUT_BYTES,))
        off = int(idx) * LUT_BYTES
        return arr[off:off + LUT_BYTES].copy()

    def get_lut_complexity(self):
        """Return LUT complexity counts: {n1: count, n2: count, n3: count}."""
        counts = np.zeros(3, dtype=np.uint32)
        self._lib.evoca_lut_complexity_counts(
            counts.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        return {'n1': int(counts[0]), 'n2': int(counts[1]), 'n3': int(counts[2])}

    def get_repro_age_hist(self):
        """Return reproduction age histogram as a numpy array (copy)."""
        n = self._lib.evoca_get_repro_age_max()
        ptr = self._lib.evoca_get_repro_age_hist()
        return np.ctypeslib.as_array(ptr, shape=(n,)).copy()

    def set_repro_age_t0(self, t):
        """Set the step at which reproduction age accumulation begins."""
        self._lib.evoca_set_repro_age_t0(int(t))

    def reset_repro_age_hist(self):
        """Clear the reproduction age histogram."""
        self._lib.evoca_reset_repro_age_hist()

    def get_step(self):
        """Return the current global step counter."""
        return int(self._lib.evoca_get_step())

    def get_eg_activity(self):
        """Return egenome activity table as dict of arrays (64 entries)."""
        acts = np.zeros(64, dtype=np.uint64)
        pops = np.zeros(64, dtype=np.uint32)
        cols = np.zeros(64, dtype=np.int32)
        self._lib.evoca_eg_activity_get(
            acts.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            pops.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            cols.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)))
        return {'activity': acts, 'pop_count': pops, 'color': cols}

    def get_activity(self, max_n=4096):
        """Return activity table as dict of arrays.

        Returns dict with keys: 'hash', 'activity', 'pop_count', 'color',
        each a numpy array of length n (actual number of entries).
        """
        keys = np.zeros(max_n, dtype=np.uint32)
        acts = np.zeros(max_n, dtype=np.uint64)
        pops = np.zeros(max_n, dtype=np.uint32)
        cols = np.zeros(max_n, dtype=np.int32)
        n = self._lib.evoca_activity_get(
            keys.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            acts.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            pops.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            cols.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            max_n)
        return {'hash': keys[:n], 'activity': acts[:n],
                'pop_count': pops[:n], 'color': cols[:n]}

    def plot_eg(self, n=8):
        """Plot the top-n egenomes by population fraction.

        Each subplot shows the 5x5 fiducial pattern on a 0-1 grid.
        Title: integer value (binary) and population fraction.
        """
        import matplotlib.pyplot as plt

        eg = self.get_eg_activity()
        pops = eg['pop_count']
        total = pops.sum()
        if total == 0:
            print("No alive cells.")
            return

        fracs = pops.astype(np.float64) / total
        top_idx = np.argsort(fracs)[::-1][:n]
        # skip zero-population entries
        top_idx = top_idx[pops[top_idx] > 0]
        n_show = len(top_idx)
        if n_show == 0:
            print("No alive cells.")
            return

        cols = min(n_show, 8)
        rows = (n_show + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(cols * 1.4, rows * 1.6))
        if n_show == 1:
            axes = np.array([axes])
        axes = np.atleast_2d(axes)

        for i, eg_val in enumerate(top_idx):
            ax = axes[i // cols, i % cols]
            pat = egenome_to_pattern(eg_val)
            ax.imshow(pat, cmap='Greys', vmin=0, vmax=1, interpolation='nearest')
            for edge in range(6):
                ax.axhline(edge - 0.5, color='gray', linewidth=0.3)
                ax.axvline(edge - 0.5, color='gray', linewidth=0.3)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{eg_val} (0b{eg_val:06b})\n{fracs[eg_val]:.3f}", fontsize=7)

        # hide unused axes
        for i in range(n_show, rows * cols):
            axes[i // cols, i % cols].set_visible(False)

        fig.tight_layout()
        plt.show()

    # ── Recipe export ────────────────────────────────────────────────

    def export_recipe(self, descriptor, probes=None, colormode=0):
        """Export current initialization recipe to Runs/<date>_<descriptor>.evoca.

        Returns the filepath written.
        """
        import json, os
        from datetime import datetime

        runs_dir = os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), 'Runs')
        os.makedirs(runs_dir, exist_ok=True)

        safe_desc = descriptor.replace(' ', '_')
        filename = f"{datetime.now().strftime('%Y-%m-%d')}_{safe_desc}.evoca"
        filepath = os.path.join(runs_dir, filename)

        recipe = {
            'version': 3,
            'created': datetime.now().isoformat(timespec='seconds'),
            'descriptor': descriptor,
            'N': self._N,
            'metaparams_init': dict(self._init_metaparams),
            'metaparams_final': {
                'food_inc': self.food_inc,
                'm_scale': self.m_scale,
                'gdiff': self.gdiff,
                'mu_lut': self.mu_lut,
                'mu_egenome': self.mu_egenome,
                'tax': self.tax,
                'restricted_mu': self.restricted_mu,
            },
            'initialization': dict(self._state_params),
            'display': {
                'colormode': colormode,
                'probes': probes or {},
            },
        }

        with open(filepath, 'w') as f:
            json.dump(recipe, f, indent=2)
        return filepath

    @property
    def N(self):
        return self._N

    @property
    def n_ent(self):
        return self._n_ent

    @property
    def cell_px(self):
        return self._lib.evoca_get_cell_px()


# ── State initialization catalogue ────────────────────────────────────

_AVAILABLE_STATE_INIT = {
    'lut':            "'gol' (GoL rule for all) | 'random' (independent random LUT per cell)",
    'lut_n_init':     "int 1-3: ring depth for random LUT (1=10 bits, 2=50, 3=all 250). Default 3",
    'egenome':        "'uniform' (same value for all) | 'random' (random in 0..63)",
    'egenome_value':  "int 0-63: egenome value when egenome='uniform'. Default 0",
    'v_density':      "float 0-1: fraction of cells with v=1. Default 0.5",
    'f_init':         "float: initial private food per cell. Default 0.0",
    'F':              "'uniform' (constant) | 'random' (uniform in F_range)",
    'F_init':         "float: env food level when F='uniform'. Default 0.0",
    'F_range':        "[lo, hi]: range when F='random'. Default [0.0, 1.0]",
    'alive':          "'all' | 'fraction' | 'patch' | 'halfplane'",
    'alive_fraction': "float 0-1: fraction alive (alive='fraction'). Default 0.5",
    'alive_radius':   "int: half-side of central square (alive='patch'). Default 64",
    'alive_axis':     "int 0|1: 0=left half, 1=top half (alive='halfplane'). Default 0",
}


def available_state_init():
    """Return a dict mapping state-init parameter names to descriptions."""
    return dict(_AVAILABLE_STATE_INIT)


# ── Recipe import ─────────────────────────────────────────────────────

def import_run(filepath=None, recipe='init', lib_path=None):
    """Load a .evoca recipe file and return (sim, display_kwargs).

    If called with no arguments, lists available recipes in Runs/.

    recipe: 'init' uses the initial metaparams, 'final' uses the
            metaparams at export time (after any slider adjustments).

    Usage:
        import_run()                    # list available recipes
        sim, kw = import_run('Runs/2026-03-15_my_run.evoca')
        run_with_controls(sim, **kw)
    """
    import json
    from pathlib import Path

    if filepath is None:
        runs_dir = Path(__file__).resolve().parent.parent / 'Runs'
        if not runs_dir.is_dir():
            print("No Runs/ directory found.")
            return []
        recipes = sorted(runs_dir.glob('*.evoca'))
        if not recipes:
            print("No .evoca files in Runs/.")
            return []
        for r in recipes:
            print(r.name)
        return [str(r) for r in recipes]

    with open(filepath) as f:
        data = json.load(f)

    # version 1 had a single 'metaparams' key; version 2 has init/final
    if 'metaparams_init' in data:
        mp = data['metaparams_init'] if recipe == 'init' else data['metaparams_final']
    else:
        mp = data['metaparams']  # v1 fallback

    sim = EvoCA(lib_path=lib_path)
    sim.init(data['N'], **mp)

    # Build state params — map v2 key names to v3 if needed
    init_raw = data.get('initialization', {})
    _v2_keys = {'lut_method': 'lut', 'egenome_method': 'egenome',
                'F_method': 'F', 'alive_method': 'alive'}
    state_kw = {_v2_keys.get(k, k): v for k, v in init_raw.items()}
    state_kw.pop('v_method', None)   # removed in v3
    sim.state(**state_kw)

    # Build display kwargs
    disp = data.get('display', {})
    display_kwargs = {}
    if disp.get('colormode', 0) != 0:
        display_kwargs['colormode'] = disp['colormode']
    if disp.get('probes'):
        display_kwargs['probes'] = disp['probes']

    return sim, display_kwargs


# ── LUT construction helpers ───────────────────────────────────────────

def pack_lut(bits):
    """Pack a LUT_BITS-length uint8 array of 0/1 values into LUT_BYTES bytes."""
    packed = np.zeros(LUT_BYTES, dtype=np.uint8)
    for i, b in enumerate(bits[:LUT_BITS]):
        if b:
            packed[i >> 3] |= np.uint8(1 << (i & 7))
    return packed


def unpack_lut(packed):
    """Unpack LUT_BYTES bytes into a LUT_BITS-length uint8 0/1 array."""
    bits = np.zeros(LUT_BITS, dtype=np.uint8)
    for i in range(LUT_BITS):
        bits[i] = (packed[i >> 3] >> (i & 7)) & 1
    return bits


def make_gol_lut():
    """
    Build a bit-packed LUT implementing exact Conway's Game of Life.

    GoL rule (B3/S23) on the Moore neighbourhood:
      dead  cell (v_x=0): birth  iff n1+n2 == 3
      alive cell (v_x=1): survive iff n1+n2 in {2, 3}
      n3 is ignored — GoL does not use it.

    Returns a LUT_BYTES-length uint8 array (bit-packed).
    """
    bits = np.zeros(LUT_BITS, dtype=np.uint8)
    for v_x in range(2):
        for n1 in range(5):
            for n2 in range(5):
                moore = n1 + n2
                if v_x == 0:
                    new_state = 1 if moore == 3 else 0
                else:
                    new_state = 1 if moore in (2, 3) else 0
                for n3 in range(5):
                    bits[lut_bit_index(v_x, n1, n2, n3)] = new_state
    return pack_lut(bits)
