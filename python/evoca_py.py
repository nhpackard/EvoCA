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

# Egene representation: 6 D4 orbits, ternary 0/1/wildcard per orbit.
# EGENE_VALUE_COUNT = 64  (legacy 6-bit value-only space; still used by
#                         set_egenome_all and the value-only egenome
#                         scratch view).
# EGENE_KEY_COUNT  = 729  (3^6, the space of distinct ternary egenes).
EGENE_VALUE_COUNT = 64
EGENE_KEY_COUNT   = 729

# Dyn-activity: 250 LUT inputs × 2 outputs = 500 transition buckets.
DYN_ACT_COUNT = LUT_BITS * 2   # 500


def egene_decode_ternary(key):
    """Inverse of the C-side egene_ternary_key. Returns (value, mask)
    as 6-bit ints. Wildcard digits zero both bits in the value, so
    two egenes that differ only at wildcard positions decode
    identically — useful when interpreting a 729-key colour swatch."""
    v = m = 0
    for i in range(6):
        digit = key % 3
        key //= 3
        if digit:
            m |= 1 << i
            if digit == 2:
                v |= 1 << i
    return v, m


def egene_ternary_key(value, mask):
    """Inverse of egene_decode_ternary."""
    key = 0
    p = 1
    for i in range(6):
        if not ((mask >> i) & 1):
            d = 0
        else:
            d = 2 if ((value >> i) & 1) else 1
        key += d * p
        p *= 3
    return key

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


# Default state planted at every stamped cell. v=0/f=0 mirrors the
# convention used by the alive setters (dead cells carry zeroed v/f);
# a freshly stamped organism therefore starts from a neutral substrate
# state and earns its own food, exactly as a reproduction child would.
_STAMP_DEFAULT_V = 0
_STAMP_DEFAULT_F = 0.0


class PatchGenomes:
    """A genome-level snapshot of a square lattice region.

    Holds, in patch-local coordinates (0..side-1), every datum needed to
    reconstruct each organism's heritable genome: the 32-byte LUT, all
    NEGENOME_MAX egene value/mask slot bytes, the per-cell active mask,
    and which cells were alive. Pixel-level state (v, f) is intentionally
    NOT carried — transfer is of recipes, not of running configurations
    (Docs/EvoCA_research_directions.md §11; egene_discussion.md
    "evolutionary test").
    """

    __slots__ = ("side", "src_r0", "src_c0", "src_N",
                 "lut", "egenes", "egenes_mask", "active", "alive")

    def __init__(self, side, src_r0, src_c0, src_N,
                 lut, egenes, egenes_mask, active, alive):
        self.side        = int(side)
        self.src_r0      = int(src_r0)   # provenance (self-transfer check)
        self.src_c0      = int(src_c0)
        self.src_N       = int(src_N)
        self.lut         = lut           # (side, side, LUT_BYTES) uint8
        self.egenes      = egenes        # (side, side, 8)         uint8
        self.egenes_mask = egenes_mask   # (side, side, 8)         uint8
        self.active      = active        # (side, side)            uint8
        self.alive       = alive         # (side, side)            uint8

    @property
    def n_alive(self):
        return int(self.alive.sum())

    @property
    def alive_fraction(self):
        return float(self.alive.mean())

    def __repr__(self):
        return (f"PatchGenomes(side={self.side}, "
                f"src=({self.src_r0},{self.src_c0}) of N={self.src_N}, "
                f"n_alive={self.n_alive})")


class EvoCA:
    """Thin ctypes wrapper around the EvoCA shared library."""

    def __init__(self, lib_path=None):
        self._lib        = ctypes.CDLL(lib_path or _find_lib())
        self._N          = 0
        self.food_inc    = 0.0
        self.m_scale     = 0.0
        self.gdiff       = 0
        self.mu_lut      = 0.0
        self.mu_egene   = 0.0
        self.mu_egenome  = 0.0      # per-active-bit flip rate
        self.p_dup_egene = 1.0
        self.tax         = 0.0
        self.tax_per_egene = 0.0
        self.tax_lut       = 0.0
        self.tax_ring      = 0.0
        self.restricted_mu = 0
        self.food_per      = 1      # env food-regen period (>=1); 1 = no-op
        self.n_chem        = 1      # chem compound-universe size (m1: only 1)
        self.food_cap      = 1      # 1 = mouthful capped at available F (A9 fix); 0 = legacy
        self._n_ent      = 2
        self.egenome      = 0
        # Diagnostic logger sample interval (Python-side; not a C
        # parameter). 1 = log every tick. Used by run_with_controls
        # to throttle ProbeLogs CSV writes.
        self.N_log_interval = 1
        # Origin recipe basename if loaded via import_run; tags log files.
        self._origin_recipe = None
        self._state_params = {}
        self._setup_signatures()

    def _setup_signatures(self):
        L = self._lib
        L.evoca_init.argtypes           = [ctypes.c_int, ctypes.c_float,
                                            ctypes.c_float]
        L.evoca_init.restype            = None
        L.evoca_free.argtypes           = []
        L.evoca_free.restype            = None
        L.evoca_set_seed.argtypes       = [ctypes.c_uint32]
        L.evoca_set_seed.restype        = None
        L.evoca_set_food_inc.argtypes   = [ctypes.c_float]
        L.evoca_set_food_inc.restype    = None
        L.evoca_set_m_scale.argtypes    = [ctypes.c_float]
        L.evoca_set_m_scale.restype     = None
        # Chem machinery (milestone 1: only N_chem=1 supported).
        L.evoca_set_n_chem.argtypes     = [ctypes.c_int]
        L.evoca_set_n_chem.restype      = ctypes.c_int
        L.evoca_get_n_chem.argtypes     = []
        L.evoca_get_n_chem.restype      = ctypes.c_int
        L.evoca_set_chem_energy.argtypes = [ctypes.c_int, ctypes.c_float]
        L.evoca_set_chem_energy.restype  = None
        L.evoca_get_chem_energy.argtypes = [ctypes.c_int]
        L.evoca_get_chem_energy.restype  = ctypes.c_float
        # Chem milestone 2: environment fields, exchange, metabolism.
        L.evoca_set_chem_env.argtypes      = [ctypes.c_int, ctypes.c_float, ctypes.c_float, ctypes.c_float]
        L.evoca_set_chem_env.restype       = None
        L.evoca_set_chem_exchange.argtypes = [ctypes.c_int, ctypes.c_float, ctypes.c_float]
        L.evoca_set_chem_exchange.restype  = None
        L.evoca_set_metab_product.argtypes = [ctypes.c_int]
        L.evoca_set_metab_product.restype  = None
        L.evoca_get_metab_product.argtypes = []
        L.evoca_get_metab_product.restype  = ctypes.c_int
        L.evoca_get_chem_params.argtypes   = [ctypes.c_int, ctypes.POINTER(ctypes.c_float)]
        L.evoca_get_chem_params.restype    = None
        L.evoca_get_chem_env.argtypes      = [ctypes.c_int]
        L.evoca_get_chem_env.restype       = ctypes.POINTER(ctypes.c_float)
        L.evoca_get_chem_priv.argtypes     = [ctypes.c_int, ctypes.POINTER(ctypes.c_float)]
        L.evoca_get_chem_priv.restype      = None
        L.evoca_get_chem_flux.argtypes     = [ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
        L.evoca_get_chem_flux.restype      = None
        L.evoca_get_tox_deaths_last.argtypes = []
        L.evoca_get_tox_deaths_last.restype  = ctypes.c_long
        L.evoca_set_chem_priv.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_float)]
        L.evoca_set_chem_priv.restype  = None
        L.evoca_set_gdiff.argtypes      = [ctypes.c_float]
        L.evoca_set_gdiff.restype       = None
        L.evoca_get_gdiff.argtypes      = []
        L.evoca_get_gdiff.restype       = ctypes.c_float
        L.evoca_set_mu_lut.argtypes     = [ctypes.c_float]
        L.evoca_set_mu_lut.restype      = None
        L.evoca_set_mu_egene.argtypes  = [ctypes.c_float]
        L.evoca_set_mu_egene.restype   = None
        L.evoca_set_mu_egenome.argtypes = [ctypes.c_float]
        L.evoca_set_mu_egenome.restype  = None
        L.evoca_set_p_dup_egene.argtypes = [ctypes.c_float]
        L.evoca_set_p_dup_egene.restype  = None
        L.evoca_get_mu_lut.argtypes     = []
        L.evoca_get_mu_lut.restype      = ctypes.c_float
        L.evoca_get_mu_egene.argtypes  = []
        L.evoca_get_mu_egene.restype   = ctypes.c_float
        L.evoca_get_mu_egenome.argtypes = []
        L.evoca_get_mu_egenome.restype  = ctypes.c_float
        L.evoca_get_p_dup_egene.argtypes = []
        L.evoca_get_p_dup_egene.restype  = ctypes.c_float
        L.evoca_set_tax.argtypes        = [ctypes.c_float]
        L.evoca_set_tax.restype         = None
        L.evoca_set_tax_per_egene.argtypes = [ctypes.c_float]
        L.evoca_set_tax_per_egene.restype  = None
        L.evoca_set_tax_lut.argtypes    = [ctypes.c_float]
        L.evoca_set_tax_lut.restype     = None
        L.evoca_set_tax_ring.argtypes   = [ctypes.c_float]
        L.evoca_set_tax_ring.restype    = None
        L.evoca_get_tax.argtypes        = []
        L.evoca_get_tax.restype         = ctypes.c_float
        L.evoca_get_tax_per_egene.argtypes = []
        L.evoca_get_tax_per_egene.restype  = ctypes.c_float
        L.evoca_get_tax_lut.argtypes    = []
        L.evoca_get_tax_lut.restype     = ctypes.c_float
        L.evoca_get_tax_ring.argtypes   = []
        L.evoca_get_tax_ring.restype    = ctypes.c_float
        L.evoca_set_restricted_mu.argtypes = [ctypes.c_int]
        L.evoca_set_restricted_mu.restype  = None
        L.evoca_get_restricted_mu.argtypes = []
        L.evoca_get_restricted_mu.restype  = ctypes.c_int
        L.evoca_set_food_per.argtypes      = [ctypes.c_int]
        L.evoca_set_food_per.restype       = None
        L.evoca_get_food_per.argtypes      = []
        L.evoca_get_food_per.restype       = ctypes.c_int
        L.evoca_set_food_cap.argtypes      = [ctypes.c_int]
        L.evoca_set_food_cap.restype       = None
        L.evoca_get_food_cap.argtypes      = []
        L.evoca_get_food_cap.restype       = ctypes.c_int
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
        L.evoca_set_egenome_random_all.argtypes = [ctypes.c_uint8]
        L.evoca_set_egenome_random_all.restype  = None
        L.evoca_set_egenome_pair_all.argtypes = [ctypes.c_uint8,
                                                  ctypes.c_uint8]
        L.evoca_set_egenome_pair_all.restype  = None
        L.evoca_get_egenes.argtypes      = []
        L.evoca_get_egenes.restype       = ctypes.POINTER(ctypes.c_uint8)
        L.evoca_get_active.argtypes      = []
        L.evoca_get_active.restype       = ctypes.POINTER(ctypes.c_uint8)
        L.evoca_set_f_all.argtypes      = [ctypes.c_float]
        L.evoca_set_f_all.restype       = None
        L.evoca_chem_sync_from_f.argtypes = []
        L.evoca_chem_sync_from_f.restype  = None
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
        L.evoca_colormode_count.argtypes = []
        L.evoca_colormode_count.restype  = ctypes.c_int
        L.evoca_get_v.argtypes          = []
        L.evoca_get_v.restype           = ctypes.POINTER(ctypes.c_uint8)
        L.evoca_get_F.argtypes          = []
        L.evoca_get_F.restype           = ctypes.POINTER(ctypes.c_float)
        L.evoca_get_f.argtypes          = []
        L.evoca_get_f.restype           = ctypes.POINTER(ctypes.c_float)
        L.evoca_get_mu_lut_gene.argtypes = []
        L.evoca_get_mu_lut_gene.restype  = ctypes.POINTER(ctypes.c_float)
        L.evoca_set_mu_lut_gene_all.argtypes = [ctypes.c_float]
        L.evoca_set_mu_lut_gene_all.restype  = None
        L.evoca_set_mu_lut_gene.argtypes = [ctypes.POINTER(ctypes.c_float)]
        L.evoca_set_mu_lut_gene.restype  = None
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
        # Effective-rule activity (eff_activity)
        L.evoca_eff_activity_update.argtypes = []
        L.evoca_eff_activity_update.restype  = None
        L.evoca_eff_activity_render_col.argtypes = [
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int]
        L.evoca_eff_activity_render_col.restype  = None
        L.evoca_eff_activity_get.argtypes = [
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_uint64),
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_int32),
            ctypes.c_int]
        L.evoca_eff_activity_get.restype  = ctypes.c_int
        L.evoca_eff_count_distinct.argtypes = []
        L.evoca_eff_count_distinct.restype  = ctypes.c_int
        L.evoca_get_eff_hashes.argtypes = [ctypes.POINTER(ctypes.c_uint32)]
        L.evoca_get_eff_hashes.restype  = None
        L.evoca_get_eff_used_mask.argtypes = [ctypes.POINTER(ctypes.c_uint8)]
        L.evoca_get_eff_used_mask.restype  = None
        L.evoca_get_eff_union_mask.argtypes = [ctypes.POINTER(ctypes.c_uint8)]
        L.evoca_get_eff_union_mask.restype  = None
        L.evoca_set_eff_act_ymax.argtypes = [ctypes.c_int]
        L.evoca_set_eff_act_ymax.restype  = None
        L.evoca_get_eff_act_ymax.argtypes = []
        L.evoca_get_eff_act_ymax.restype  = ctypes.c_int
        # Activity quantile probe
        L.evoca_q_activity_deciles.argtypes = [ctypes.POINTER(ctypes.c_float)]
        L.evoca_q_activity_deciles.restype  = None
        # Per-step demography counters
        L.evoca_get_births_last.argtypes = []
        L.evoca_get_births_last.restype  = ctypes.c_int
        L.evoca_get_deaths_last.argtypes = []
        L.evoca_get_deaths_last.restype  = ctypes.c_int
        # Realised per-component mutation flux (fixed-space shadows)
        L.evoca_get_eg_bitflips_last.argtypes  = []
        L.evoca_get_eg_bitflips_last.restype   = ctypes.c_long
        L.evoca_get_lut_bitflips_last.argtypes = []
        L.evoca_get_lut_bitflips_last.restype  = ctypes.c_long
        L.evoca_get_repro_last.argtypes        = []
        L.evoca_get_repro_last.restype         = ctypes.c_long
        # Fixed-space neutral baselines (eg_activity / dyn_activity)
        L.evoca_eg_excess_pc.argtypes      = []
        L.evoca_eg_excess_pc.restype       = ctypes.c_double
        L.evoca_dyn_excess_pc.argtypes     = []
        L.evoca_dyn_excess_pc.restype      = ctypes.c_double
        L.evoca_eg_shadow_total.argtypes   = []
        L.evoca_eg_shadow_total.restype    = ctypes.c_double
        L.evoca_dyn_shadow_total.argtypes  = []
        L.evoca_dyn_shadow_total.restype   = ctypes.c_double
        L.evoca_eg_zscores.argtypes        = [ctypes.POINTER(ctypes.c_double)]
        L.evoca_eg_zscores.restype         = None
        L.evoca_dyn_zscores.argtypes       = [ctypes.POINTER(ctypes.c_double)]
        L.evoca_dyn_zscores.restype        = None
        L.evoca_eg_nbuckets.argtypes       = []
        L.evoca_eg_nbuckets.restype        = ctypes.c_int
        L.evoca_dyn_nbuckets.argtypes      = []
        L.evoca_dyn_nbuckets.restype       = ctypes.c_int
        # Neutral shadow population (Channon-style activity calibration)
        L.evoca_neutral_enable.argtypes        = []
        L.evoca_neutral_enable.restype         = None
        L.evoca_neutral_disable.argtypes       = []
        L.evoca_neutral_disable.restype        = None
        L.evoca_neutral_is_enabled.argtypes    = []
        L.evoca_neutral_is_enabled.restype     = ctypes.c_int
        L.evoca_neutral_get_population.argtypes = []
        L.evoca_neutral_get_population.restype  = ctypes.c_int
        # N-activity probe (shadow population, same hash as G-activity)
        L.evoca_n_activity_update.argtypes     = []
        L.evoca_n_activity_update.restype      = None
        L.evoca_n_activity_render_col.argtypes = [
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int]
        L.evoca_n_activity_render_col.restype  = None
        L.evoca_n_activity_get.argtypes        = [
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_uint64),
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_int32),
            ctypes.c_int]
        L.evoca_n_activity_get.restype         = ctypes.c_int
        L.evoca_set_n_act_ymax.argtypes        = [ctypes.c_int]
        L.evoca_set_n_act_ymax.restype         = None
        L.evoca_get_n_act_ymax.argtypes        = []
        L.evoca_get_n_act_ymax.restype         = ctypes.c_int
        L.evoca_set_n_overlay_alpha.argtypes   = [ctypes.c_int]
        L.evoca_set_n_overlay_alpha.restype    = None
        L.evoca_get_n_overlay_alpha.argtypes   = []
        L.evoca_get_n_overlay_alpha.restype    = ctypes.c_int
        L.evoca_get_n_p90.argtypes             = []
        L.evoca_get_n_p90.restype              = ctypes.c_double
        L.evoca_nq_activity_deciles.argtypes   = [ctypes.POINTER(ctypes.c_float)]
        L.evoca_nq_activity_deciles.restype    = None
        # Activity flux probe
        L.evoca_activity_crossings.argtypes = [
            ctypes.c_uint64, ctypes.c_uint64, ctypes.c_int,
            ctypes.POINTER(ctypes.c_float), ctypes.c_int]
        L.evoca_activity_crossings.restype  = ctypes.c_int
        L.evoca_count_distinct_genomes.argtypes = []
        L.evoca_count_distinct_genomes.restype  = ctypes.c_int
        L.evoca_get_population.argtypes = []
        L.evoca_get_population.restype  = ctypes.c_int
        L.evoca_get_ages.argtypes       = [ctypes.POINTER(ctypes.c_int32)]
        L.evoca_get_ages.restype        = None
        # Opt-in lineage record (research-directions §1)
        L.evoca_enable_lineage.argtypes = [ctypes.c_int]
        L.evoca_enable_lineage.restype  = None
        L.evoca_lineage_enabled.argtypes = []
        L.evoca_lineage_enabled.restype  = ctypes.c_int
        L.evoca_get_lineage_parent_hash.argtypes = [
            ctypes.POINTER(ctypes.c_uint32)]
        L.evoca_get_lineage_parent_hash.restype  = None
        L.evoca_get_lut_hashes.argtypes = [
            ctypes.POINTER(ctypes.c_uint32)]
        L.evoca_get_lut_hashes.restype  = None
        L.evoca_get_lineage_birth_id.argtypes = [
            ctypes.POINTER(ctypes.c_uint64)]
        L.evoca_get_lineage_birth_id.restype  = None
        L.evoca_get_lineage_parent_id.argtypes = [
            ctypes.POINTER(ctypes.c_uint64)]
        L.evoca_get_lineage_parent_id.restype  = None
        # Per-lineage net metabolic return reducer (checkpoint §D2)
        L.evoca_lineage_return_update.argtypes = []
        L.evoca_lineage_return_update.restype  = None
        L.evoca_lineage_return_count.argtypes  = []
        L.evoca_lineage_return_count.restype   = ctypes.c_int
        L.evoca_get_lineage_return.argtypes = [
            ctypes.POINTER(ctypes.c_uint32),   # hash
            ctypes.POINTER(ctypes.c_double),   # intake
            ctypes.POINTER(ctypes.c_double),   # tax
            ctypes.POINTER(ctypes.c_uint64),   # births
            ctypes.POINTER(ctypes.c_uint32),   # pop
            ctypes.POINTER(ctypes.c_uint32),   # first_step
            ctypes.POINTER(ctypes.c_uint32)]   # last_step
        L.evoca_get_lineage_return.restype  = ctypes.c_int
        L.evoca_get_lineage_intake_total.argtypes = []
        L.evoca_get_lineage_intake_total.restype  = ctypes.c_double
        L.evoca_get_lineage_tax_total.argtypes = []
        L.evoca_get_lineage_tax_total.restype  = ctypes.c_double
        L.evoca_get_lineage_intake.argtypes = [
            ctypes.POINTER(ctypes.c_double)]
        L.evoca_get_lineage_intake.restype  = None
        L.evoca_get_lineage_tax.argtypes = [
            ctypes.POINTER(ctypes.c_double)]
        L.evoca_get_lineage_tax.restype  = None
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
        # Egene food intake probe
        L.evoca_eg_food_render_col.argtypes = [
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int]
        L.evoca_eg_food_render_col.restype  = None
        L.evoca_eg_food_get.argtypes = [
            ctypes.POINTER(ctypes.c_uint64),
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_int32),
        ]
        L.evoca_eg_food_get.restype  = ctypes.c_int
        L.evoca_set_eg_food_ymax.argtypes = [ctypes.c_int]
        L.evoca_set_eg_food_ymax.restype  = None
        L.evoca_get_eg_food_ymax.argtypes = []
        L.evoca_get_eg_food_ymax.restype  = ctypes.c_int
        # Egenome scalar stats: out[0..4] =
        #   [mean Negene, std Negene, distinct egene values,
        #    mean max-match, frac at Negene_max]
        L.evoca_egenome_stats.argtypes = [ctypes.POINTER(ctypes.c_float)]
        L.evoca_egenome_stats.restype  = None
        # Egene cognitive stats: out[0..2] =
        #   [mean cognitive specificity (non-* cell-positions per
        #    active egene), mean per-cell cognitive load,
        #    mean intake per alive eater this step]
        L.evoca_egene_stats.argtypes  = [ctypes.POINTER(ctypes.c_float)]
        L.evoca_egene_stats.restype   = None
        # Dissipation / energy-flux accounting (checkpoint §L): out[0..15]
        L.evoca_get_dissipation.argtypes = [ctypes.POINTER(ctypes.c_double)]
        L.evoca_get_dissipation.restype  = None
        # Dyn-activity probe (LUT-entry-level): 500 buckets
        L.evoca_dyn_activity_render_col.argtypes = [
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int]
        L.evoca_dyn_activity_render_col.restype  = None
        L.evoca_dyn_activity_get.argtypes = [
            ctypes.POINTER(ctypes.c_uint64),
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_int32),
        ]
        L.evoca_dyn_activity_get.restype  = ctypes.c_int
        L.evoca_set_dyn_act_ymax.argtypes = [ctypes.c_int]
        L.evoca_set_dyn_act_ymax.restype  = None
        L.evoca_get_dyn_act_ymax.argtypes = []
        L.evoca_get_dyn_act_ymax.restype  = ctypes.c_int
        L.evoca_get_dyn_distinct.argtypes = []
        L.evoca_get_dyn_distinct.restype  = ctypes.c_int
        # Egene mask accessor: per-slot wildcard-mask byte (bit i = 1
        # means orbit i is non-wildcard).
        L.evoca_get_egenes_mask.argtypes = []
        L.evoca_get_egenes_mask.restype  = ctypes.POINTER(ctypes.c_uint8)
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

    def init(self, N, food_inc=0.0, m_scale=1.0, gdiff=0.0,
             mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0,
             p_dup_egene=1.0,
             tax=0.0, tax_per_egene=0.0, tax_lut=0.0, tax_ring=0.0,
             restricted_mu=0, food_per=1, n_ent=2,
             N_log_interval=1, lineage=False, n_chem=1, food_cap=1,
             wdiff=0.0, waste_sink=0.0, waste_disp=0.0):
        # Opt-in lineage record (research-directions §1). The C flag must
        # be set before evoca_init so the arrays get allocated. Default
        # OFF -> zero added work in the Phase-1/Phase-4 hot path.
        self._lib.evoca_enable_lineage(1 if lineage else 0)
        stop = getattr(self, '_stop_display', None)
        if stop is not None:
            stop()
            self._stop_display = None
        self._N         = N
        self.food_inc   = float(food_inc)
        self.m_scale    = float(m_scale)
        self.gdiff      = float(gdiff)
        self.mu_lut     = float(mu_lut)
        self.mu_egene  = float(mu_egene)
        self.mu_egenome = float(mu_egenome)
        self.p_dup_egene = float(p_dup_egene)
        self.tax        = float(tax)
        self.tax_per_egene = float(tax_per_egene)
        self.tax_lut       = float(tax_lut)
        self.tax_ring      = float(tax_ring)
        self.restricted_mu = int(restricted_mu)
        self.food_per   = max(1, int(food_per))
        self.food_cap   = 1 if int(food_cap) else 0
        # Food/waste chemistry (chem milestone 2; Docs/plan-chem.md §7).
        # n_chem=2 is the food/waste preset; the waste metaparams need it.
        self.wdiff      = float(wdiff)
        self.waste_sink = float(waste_sink)
        self.waste_disp = float(waste_disp)
        if int(n_chem) not in (1, 2):
            raise ValueError(f"n_chem={n_chem}: only 1 (single food) and 2 (food/waste) "
                             f"are configured; larger universes are chem milestone 3")
        if int(n_chem) == 1 and (self.wdiff or self.waste_sink or self.waste_disp):
            raise ValueError("wdiff / waste_sink / waste_disp need n_chem=2 (food/waste)")
        self._n_ent     = int(n_ent)
        self.N_log_interval = max(1, int(N_log_interval))
        self._lib.evoca_init(N, self.food_inc, self.m_scale)
        # Chem compound-universe size: 1 = single food (the standard model),
        # 2 = food/waste, configured through the general chem API.
        self.n_chem = int(n_chem)
        self._lib.evoca_set_n_chem(self.n_chem)
        self.n_chem = int(self._lib.evoca_get_n_chem())  # reflect any clamp
        if self.n_chem == 2:
            self._configure_food_waste()
        self._lib.evoca_set_gdiff(self.gdiff)
        self._lib.evoca_set_mu_lut(self.mu_lut)
        self._lib.evoca_set_mu_egene(self.mu_egene)
        self._lib.evoca_set_mu_egenome(self.mu_egenome)
        self._lib.evoca_set_p_dup_egene(self.p_dup_egene)
        self._lib.evoca_set_tax(self.tax)
        self._lib.evoca_set_tax_per_egene(self.tax_per_egene)
        self._lib.evoca_set_tax_lut(self.tax_lut)
        self._lib.evoca_set_tax_ring(self.tax_ring)
        self._lib.evoca_set_restricted_mu(self.restricted_mu)
        self._lib.evoca_set_food_per(self.food_per)
        # Mouthful cap (A9 fix). 1 = capped at available food (default,
        # corrected); 0 = legacy uncapped rule, for reproducing pre-fix runs.
        self._lib.evoca_set_food_cap(self.food_cap)
        # Seed the per-cell mutation-rate gene (mutation-of-mutation Wave-1
        # storage scaffold) to this run's mu_lut. The C library also seeds
        # it at evoca_init, but from the *previous* global gmu_lut (set_mu_lut
        # runs after evoca_init), so re-seed here for a clean default. The
        # gene is inherited on reproduction but read nowhere that affects
        # dynamics, so this is a pure no-op.
        self._lib.evoca_set_mu_lut_gene_all(ctypes.c_float(self.mu_lut))
        self._lib.evoca_set_n_ent(self._n_ent)
        self._state_params = {}
        self._init_metaparams = {
            'food_inc': self.food_inc,
            'm_scale': self.m_scale,
            'gdiff': self.gdiff,
            'mu_lut': self.mu_lut,
            'mu_egene': self.mu_egene,
            'mu_egenome': self.mu_egenome,
            'p_dup_egene': self.p_dup_egene,
            'tax': self.tax,
            'tax_per_egene': self.tax_per_egene,
            'tax_lut': self.tax_lut,
            'tax_ring': self.tax_ring,
            'restricted_mu': self.restricted_mu,
            'food_per': self.food_per,
            'n_chem': self.n_chem,
            'food_cap': self.food_cap,
            'wdiff': self.wdiff,
            'waste_sink': self.waste_sink,
            'waste_disp': self.waste_disp,
        }

    # ── Food/waste chemistry (chem milestone 2) ───────────────────────

    WASTE = 1   # compound id of waste in the food/waste preset

    def _configure_food_waste(self):
        """The food/waste special case of the general chemistry: compound 1
        (waste) carries no energy, kills at a private amount of 1, is what the
        tax converts food into, is excreted at rate waste_disp * w * (1 - W),
        and in the environment has no source, a uniform sink waste_sink and
        diffusion wdiff (the counterparts of food_inc and gdiff)."""
        W = self.WASTE
        self._lib.evoca_set_chem_energy(W, ctypes.c_float(0.0))
        self._lib.evoca_set_chem_exchange(W, ctypes.c_float(self.waste_disp), ctypes.c_float(1.0))
        self._lib.evoca_set_chem_env(W, ctypes.c_float(0.0), ctypes.c_float(self.waste_sink),
                                     ctypes.c_float(self.wdiff))
        self._lib.evoca_set_metab_product(W)

    def _require_waste(self, what):
        if self.n_chem != 2:
            raise ValueError(f"{what} needs the food/waste chemistry (init with n_chem=2)")

    def update_wdiff(self, d):
        """Waste diffusion (same rule as gdiff)."""
        self._require_waste('wdiff')
        self.wdiff = max(0.0, float(d)); self._configure_food_waste()

    def update_waste_sink(self, r):
        """Uniform removal of environmental waste per step (floored at 0)."""
        self._require_waste('waste_sink')
        self.waste_sink = max(0.0, float(r)); self._configure_food_waste()

    def update_waste_disp(self, r):
        """Waste excretion rate: d = waste_disp * w * max(0, 1 - W) per step."""
        self._require_waste('waste_disp')
        self.waste_disp = max(0.0, float(r)); self._configure_food_waste()

    def get_W(self):
        """Environmental waste field W (N x N copy); zeros when n_chem == 1."""
        if self.n_chem != 2:
            return np.zeros((self._N, self._N), dtype=np.float32)
        ptr = self._lib.evoca_get_chem_env(self.WASTE)
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)).copy().reshape(self._N, self._N)

    def get_w(self):
        """Private waste of every cell (N x N); zeros when n_chem == 1."""
        out = np.zeros(self._N * self._N, dtype=np.float32)
        if self.n_chem == 2:
            self._lib.evoca_get_chem_priv(self.WASTE, out.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
        return out.reshape(self._N, self._N)

    def set_w(self, w):
        """Seed the private waste of every alive cell (scalar or N x N)."""
        self._require_waste('set_w')
        a = np.ascontiguousarray(np.broadcast_to(np.float32(w) if np.isscalar(w) else np.asarray(w, np.float32),
                                                 (self._N, self._N)), dtype=np.float32).ravel()
        self._lib.evoca_set_chem_priv(self.WASTE, a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))

    def get_chem_params(self, c):
        """(inc, sink, diff, disp, tox) of compound c >= 1 as set in C."""
        out = (ctypes.c_float * 5)()
        self._lib.evoca_get_chem_params(int(c), out)
        return tuple(float(x) for x in out)

    def get_metab_product(self):
        return int(self._lib.evoca_get_metab_product())

    def get_chem_flux(self, c):
        """Cumulative (in, out) mass of compound c >= 1: in = environmental
        source + metabolic production, out = environmental sink."""
        out = (ctypes.c_double * 2)()
        self._lib.evoca_get_chem_flux(int(c), out)
        return float(out[0]), float(out[1])

    def get_tox_deaths_last(self):
        """Cells killed by a toxic compound (waste) in the last step."""
        return int(self._lib.evoca_get_tox_deaths_last())

    def free(self):
        stop = getattr(self, '_stop_display', None)
        if stop is not None:
            stop()
            self._stop_display = None
        self._lib.evoca_free()
        self._N = 0

    def set_seed(self, seed):
        """Deterministically reset the C PRNG state.

        Pool workers are reused, so the C RNG would otherwise carry
        state across configs. Reseeding identically per run makes
        results reproducible. seed==0 is remapped C-side."""
        self._lib.evoca_set_seed(ctypes.c_uint32(int(seed) & 0xFFFFFFFF))

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
        self.gdiff = float(d)
        self._lib.evoca_set_gdiff(self.gdiff)

    def update_mu_lut(self, m):
        self.mu_lut = float(m)
        self._lib.evoca_set_mu_lut(self.mu_lut)

    def update_mu_egene(self, m):
        self.mu_egene = float(m)
        self._lib.evoca_set_mu_egene(self.mu_egene)

    def update_mu_egenome(self, m):
        self.mu_egenome = float(m)
        self._lib.evoca_set_mu_egenome(self.mu_egenome)

    def update_p_dup_egene(self, p):
        self.p_dup_egene = float(p)
        self._lib.evoca_set_p_dup_egene(self.p_dup_egene)

    def update_tax(self, t):
        self.tax = float(t)
        self._lib.evoca_set_tax(self.tax)

    def update_tax_per_egene(self, t):
        self.tax_per_egene = float(t)
        self._lib.evoca_set_tax_per_egene(self.tax_per_egene)

    def update_tax_lut(self, t):
        self.tax_lut = float(t)
        self._lib.evoca_set_tax_lut(self.tax_lut)

    def update_tax_ring(self, t):
        self.tax_ring = float(t)
        self._lib.evoca_set_tax_ring(self.tax_ring)

    def update_N_log_interval(self, n):
        """Set the ProbeLogs CSV sampling interval (Python-side; the C
        library is unaware). 1 logs every tick; 5 logs every 5th, etc."""
        self.N_log_interval = max(1, int(n))

    def update_restricted_mu(self, r):
        self.restricted_mu = int(r)
        self._lib.evoca_set_restricted_mu(self.restricted_mu)

    def update_food_per(self, p):
        """Set the environmental food-regeneration period (>=1).

        p=1 regenerates every step (bit-for-bit the pre-food_per
        behaviour); p=K fires only every K-th step, adding
        food_inc*K when it fires so the mean input rate is unchanged
        but temporal lumpiness rises."""
        self.food_per = max(1, int(p))
        self._lib.evoca_set_food_per(self.food_per)

    def update_food_cap(self, c):
        """Toggle the Phase-3 mouthful cap (A9 fix).

        c=1 (default) caps each mouthful at the food present, max(F,0),
        so F_food never goes negative; c=0 restores the legacy uncapped
        rule (overdraws F when m_scale*score > 25) to reproduce pre-fix
        runs. At m_scale <= 1 the two are identical."""
        self.food_cap = 1 if int(c) else 0
        self._lib.evoca_set_food_cap(self.food_cap)

    def update_act_ymax(self, y):
        self._lib.evoca_set_act_ymax(int(y))

    def update_eg_act_ymax(self, y):
        self._lib.evoca_set_eg_act_ymax(int(y))

    def update_eg_food_ymax(self, y):
        self._lib.evoca_set_eg_food_ymax(int(y))

    def update_dyn_act_ymax(self, y):
        self._lib.evoca_set_dyn_act_ymax(int(y))

    def update_pat_act_ymax(self, y):
        self._lib.evoca_set_pat_act_ymax(int(y))

    # ── Params export ─────────────────────────────────────────────────

    _DEFAULTS = dict(food_inc=0.0, m_scale=1.0,
                      gdiff=0, mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0,
                      p_dup_egene=1.0, tax=0.0, tax_per_egene=0.0,
                      tax_lut=0.0, tax_ring=0.0,
                      restricted_mu=0, food_per=1, n_chem=1, food_cap=1,
                      wdiff=0.0, waste_sink=0.0, waste_disp=0.0)

    def params(self):
        """Return current metaparameters as a dict suitable for init(**d)."""
        return dict(N=self._N, food_inc=self.food_inc, m_scale=self.m_scale,
                    gdiff=self.gdiff, mu_lut=self.mu_lut,
                    mu_egene=self.mu_egene, mu_egenome=self.mu_egenome,
                    p_dup_egene=self.p_dup_egene, tax=self.tax,
                    tax_per_egene=self.tax_per_egene,
                    tax_lut=self.tax_lut, tax_ring=self.tax_ring,
                    restricted_mu=self.restricted_mu, food_per=self.food_per,
                    n_chem=self.n_chem, food_cap=self.food_cap,
                    wdiff=self.wdiff, waste_sink=self.waste_sink, waste_disp=self.waste_disp)

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

    # ── Env / agent metaparam grouping ────────────────────────────────
    #
    # Two *views* over the flat metaparameter namespace (checkpoint §A /
    # §I' A1). They are DERIVED from params() — the single source of
    # truth (the per-attribute storage) — so there is exactly one owner
    # of the metaparameters; do NOT fork storage into separate grouped
    # dicts. `params()`, `params_str()`, and the .evoca recipe remain
    # the flat, canonical API and are untouched.
    #
    # The niche metaphor (arctic ↔ jungle as points in env-metaparam
    # space) assigns most knobs cleanly. The one place the boundary is
    # non-obvious is the TAX FAMILY, which splits by *causal owner* — who
    # imposes the cost — not by whether the word "tax" appears:
    #   - baseline `tax` is an ENVIRONMENTAL cost-of-living: the niche
    #     imposes it (a harsh/arctic niche = high tax) -> env group.
    #   - `tax_per_egene` / `tax_lut` / `tax_ring` are AGENT-caused
    #     metabolic costs: the organism pays for the complexity of the
    #     genome it carries -> agent group. (These are the ones that
    #     would become per-cell evolvable genes; baseline tax stays an
    #     environmental control.)
    _ENV_KEYS = ('N', 'food_inc', 'gdiff', 'tax')
    _AGENT_KEYS = ('m_scale', 'mu_lut', 'mu_egene', 'mu_egenome',
                   'p_dup_egene', 'tax_per_egene', 'tax_lut', 'tax_ring',
                   'restricted_mu')

    def env_params(self):
        """Return the environment-metaparam view as a dict.

        Derived from params(); {N, food_inc, gdiff, tax}. Baseline `tax`
        is the niche's cost-of-living (see class comment). These are the
        knobs one varies to move between niches (arctic ↔ jungle)."""
        p = self.params()
        return {k: p[k] for k in self._ENV_KEYS}

    def agent_params(self):
        """Return the agent-metaparam view as a dict.

        Derived from params(); {m_scale, mu_lut, mu_egene, mu_egenome,
        p_dup_egene, tax_per_egene, tax_lut, tax_ring, restricted_mu}.
        These are the candidates for being carried per-cell and mutated
        as part of the evolutionary process."""
        p = self.params()
        return {k: p[k] for k in self._AGENT_KEYS}

    def grouped_params(self):
        """Return {'env': env_params(), 'agent': agent_params()}.

        The two views together partition params() exactly (verified by
        the union check below) — no metaparam is dropped or duplicated."""
        env, agent = self.env_params(), self.agent_params()
        # Self-consistency guard: the grouping must partition params()
        # exactly. Catches a future metaparam that gets added to
        # params() but not to _ENV_KEYS/_AGENT_KEYS.
        flat = set(self.params())
        grouped = set(env) | set(agent)
        assert grouped == flat and len(env) + len(agent) == len(flat), (
            "env/agent grouping does not partition params(): "
            f"missing={flat - grouped}, extra={grouped - flat}")
        return {'env': env, 'agent': agent}

    # ── State initialization ─────────────────────────────────────────

    @property
    def state_params(self):
        """Return a copy of the current state-initialization parameters."""
        return dict(self._state_params)

    def state(self, lut='gol', lut_n_init=3, lut_density=0.5,
              lut_uniform=False, lut_seed=None, lut_to_dead=False,
              egenome='uniform', egenome_value=0,
              v_density=0.5, f_init=0.0, F='uniform', F_init=0.0, F_range=None,
              alive='all', alive_fraction=0.5, alive_radius=64, alive_axis=0,
              v_seed=None):
        """Initialize lattice state from parameters.

        See available_state_init() for parameter descriptions. v_seed seeds
        the numpy draw of the initial v(x) field; None (default) draws it
        unseeded, so the run is not reproducible even after set_seed (the C
        PRNG does not touch v). Pass an int for a reproducible run; it is
        recorded in the recipe.
        """
        self._state_params = {}

        # LUT
        if lut == 'random':
            self.set_lut_random(n_init=lut_n_init,
                                density=lut_density,
                                uniform=lut_uniform,
                                seed=lut_seed)
        else:
            self.set_lut_all(make_gol_lut())

        # Egenome
        if egenome == 'random':
            self.set_egenome_random()
        else:
            self.set_egenome_all(egenome_value)

        # v(x) — random binary with given density
        rng = np.random.default_rng(v_seed)
        self.set_v((rng.random((self._N, self._N)) < v_density).astype(np.uint8),
                   density=v_density)
        if v_seed is not None:
            self._state_params['v_seed'] = int(v_seed)

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

        # Optional: re-apply LUT (and egenome) to all cells AFTER the
        # alive setter zeroed dead cells' rules. Gives a uniform CA
        # substrate regardless of which cells start "alive" — dead
        # cells still don't eat/reproduce but they DO participate in
        # the CA via their own LUT (e.g. v_curr can transition there).
        if lut_to_dead:
            if lut == 'random':
                self.set_lut_random(n_init=lut_n_init,
                                    density=lut_density,
                                    uniform=lut_uniform,
                                    seed=lut_seed)
            else:
                self.set_lut_all(make_gol_lut())
            if egenome == 'random':
                self.set_egenome_random()
            else:
                self.set_egenome_all(egenome_value)
        self._state_params['lut_to_dead'] = bool(lut_to_dead)

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
        """Random egene init: each cell gets 8 independent random 6-bit
        egenes and exactly one random active slot (Negene = 1).
        wt=0xFF so all 64 egene values get distinct hash-based colors."""
        self._lib.evoca_set_egenome_random_all(0xFF)
        self._state_params['egenome'] = 'random'
        self._state_params.pop('egenome_value', None)

    def set_lut_random(self, n_init=3, density=0.5, uniform=False, seed=None):
        """Set each cell's LUT to a random rule.

        n_init: number of rings the rule conditions on (1, 2, or 3).
          1: depends on (v_x, n1) — 10 independent bits per LUT
          2: depends on (v_x, n1, n2) — 50 independent bits
          3: depends on (v_x, n1, n2, n3) — all 250 bits independent
        density: probability that each independent bit is 1 (default 0.5).
          Use lower values to bias toward "stay dead" rules (sparse
          activity), higher to bias toward "stay alive".
        uniform: if True, generate one rule and assign it to every
          cell (a true uniform CA substrate). Use for "fixed-LUT
          population" experiments where you want one global dynamics
          and only the egenome to evolve.
        seed: numpy seed for reproducibility (None = use global RNG)."""
        N2 = self._N * self._N
        n_dup = 1 if uniform else N2
        rng = np.random.default_rng(seed) if seed is not None else np.random

        # Each independent bit is 1 with probability `density`.
        if n_init >= 3:
            bits = (rng.random((n_dup, 2, 5, 5, 5)) < density).astype(np.uint8)
        elif n_init == 2:
            core = (rng.random((n_dup, 2, 5, 5, 1)) < density).astype(np.uint8)
            bits = np.tile(core, (1, 1, 1, 1, 5))
        else:  # n_init == 1
            core = (rng.random((n_dup, 2, 5, 1, 1)) < density).astype(np.uint8)
            bits = np.tile(core, (1, 1, 1, 5, 5))

        if uniform:
            bits = np.tile(bits, (N2, 1, 1, 1, 1))

        flat = bits.reshape(N2, LUT_BITS)

        # Bit-pack to (N2, 32) bytes
        padded = np.zeros((N2, 256), dtype=np.uint8)
        padded[:, :LUT_BITS] = flat
        grouped = padded.reshape(N2, 32, 8)
        packed = np.zeros((N2, 32), dtype=np.uint8)
        for b in range(8):
            packed |= grouped[:, :, b] << b

        packed = np.ascontiguousarray(packed)
        if uniform:
            # All cells share row 0 of `packed`; one set_lut_all is faster
            # than N2 set_lut calls.
            self._lib.evoca_set_lut_all(
                packed[0].ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)))
        else:
            for i in range(N2):
                self.set_lut(i, packed[i])
        self._state_params['lut']        = 'random'
        self._state_params['lut_n_init'] = n_init
        self._state_params['lut_density']= float(density)
        self._state_params['lut_uniform']= bool(uniform)
        self._state_params['lut_seed']   = None if seed is None else int(seed)

    def set_f_all(self, f):
        self._lib.evoca_set_f_all(float(f))
        self._state_params['f_init'] = float(f)

    # ── Chem machinery (milestone 1: single-food parity) ───────────────
    def set_n_chem(self, n):
        """Set the compound-universe size (C: 1..8; resets every compound
        c >= 1 to inert and empties their stores). Low-level: init(n_chem=2)
        also applies the food/waste preset. Returns the effective value."""
        self._lib.evoca_set_n_chem(int(n))
        self.n_chem = int(self._lib.evoca_get_n_chem())
        return self.n_chem

    def get_food_cap(self):
        """Effective Phase-3 mouthful-cap flag from the C library."""
        return int(self._lib.evoca_get_food_cap())

    def get_n_chem(self):
        return int(self._lib.evoca_get_n_chem())

    def set_chem_energy(self, compound, energy):
        """Reduction weight for a compound in the vector->f mapping
        (default 1.0 for every compound; compound 0 is the energy 'food')."""
        self._lib.evoca_set_chem_energy(int(compound), float(energy))

    def get_chem_energy(self, compound):
        return float(self._lib.evoca_get_chem_energy(int(compound)))

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

    def get_mu_lut_gene(self):
        """Return (N, N) float32 array of the per-cell LUT mutation-rate gene.

        Mutation-of-mutation Wave-1 storage scaffold. The gene is inherited
        on reproduction but READ NOWHERE that affects dynamics (the mutation
        draw still uses the global mu_lut), and is not itself mutated — so a
        run with genes present is bit-for-bit identical to one without. Seed
        a distribution via set_mu_lut_gene()/set_mu_lut_gene_all() for the
        zero-mutation distribution-seeding assay."""
        ptr = self._lib.evoca_get_mu_lut_gene()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def set_mu_lut_gene_all(self, g):
        """Seed every cell's LUT mutation-rate gene to the scalar `g`."""
        self._lib.evoca_set_mu_lut_gene_all(ctypes.c_float(g))

    def set_mu_lut_gene(self, arr):
        """Seed the per-cell LUT mutation-rate gene from an (N, N) array.

        Accepts any shape broadcastable to (N, N); it is cast to float32,
        made C-contiguous, and copied into the C array. Use this to install
        a spatial or random distribution of rate-gene values for the
        distribution-seeding selection assay."""
        buf = np.ascontiguousarray(
            np.broadcast_to(arr, (self._N, self._N)), dtype=np.float32)
        self._lib.evoca_set_mu_lut_gene(
            buf.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))

    def get_egenome(self):
        """Return (N, N) uint8 array of fiducial genomes (6-bit values).
        Each entry is the lowest-active-slot's value byte for that cell —
        a backward-compat single-byte view of the multi-egene genome.
        For the full per-slot picture, use get_egenes() and
        get_egenes_mask() instead."""
        ptr = self._lib.evoca_get_egenome()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def get_egenes(self):
        """Return (N, N, NEGENOME_MAX) uint8 array of egene value bytes,
        one per slot. Pair with get_egenes_mask() and get_active() to
        reconstruct the full ternary egenome of any cell."""
        ptr = self._lib.evoca_get_egenes()
        return np.ctypeslib.as_array(
            ptr, shape=(self._N * self._N * 8,)).copy().reshape(
                self._N, self._N, 8)

    def get_active(self):
        """Return (N, N) uint8 array — bit i = 1 means slot i is active.
        Negene = popcount of this byte. NB: this is the per-cell active
        mask (for egene slots), distinct from get_alive() which marks
        which cells in the lattice host an organism."""
        ptr = self._lib.evoca_get_active()
        return np.ctypeslib.as_array(
            ptr, shape=(self._N * self._N,)).copy().reshape(
                self._N, self._N)

    def cell_inspect(self, row, col):
        """Convenience: return a dict describing one cell's egenome.
        Includes alive flag, Negene, the list of active (value, mask)
        pairs as ints in [0, 64), and a list of orbit-strings for each
        active slot using '.', '0', '1' for wildcard / off / on."""
        i = int(row) * self._N + int(col)
        alive = int(np.ctypeslib.as_array(
            self._lib.evoca_get_alive(),
            shape=(self._N * self._N,))[i])
        active = int(np.ctypeslib.as_array(
            self._lib.evoca_get_active(),
            shape=(self._N * self._N,))[i])
        vals = np.ctypeslib.as_array(
            self._lib.evoca_get_egenes(),
            shape=(self._N * self._N * 8,))[i*8:(i+1)*8]
        masks = np.ctypeslib.as_array(
            self._lib.evoca_get_egenes_mask(),
            shape=(self._N * self._N * 8,))[i*8:(i+1)*8]
        slots = []
        for s in range(8):
            if (active >> s) & 1:
                v, m = int(vals[s]) & 0x3F, int(masks[s]) & 0x3F
                # Build per-orbit string using '.', '0', '1'.
                orbits = ''.join(
                    '.' if not ((m >> i) & 1) else
                    ('1' if (v >> i) & 1 else '0')
                    for i in range(6))
                slots.append({'slot': s, 'value': v, 'mask': m,
                              'orbits_b0..b5': orbits})
        return {'alive': bool(alive), 'active_byte': active,
                'Negene': int(bin(active).count('1')),
                'slots': slots}

    def set_egenome_pair_all(self, value, mask):
        """Set every cell to a single active egene with the given
        6-bit (value, mask) ternary pair. Negene = 1 everywhere; only
        slot 0 active. Use this to seed runs at any starting cognition
        level — e.g. (0, 0x01) for a minimal "centre = dead" one-orbit
        specification from which evolution can grow specificity."""
        v = int(value) & 0x3F
        m = int(mask)  & 0x3F
        self._lib.evoca_set_egenome_pair_all(v, m)
        self._state_params['egenome'] = 'pair'
        self._state_params['egenome_value'] = v
        self._state_params['egenome_mask']  = m

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
        """Return egene-activity table as dict of arrays (EGENE_KEY_COUNT
        = 729 ternary keys)."""
        acts = np.zeros(EGENE_KEY_COUNT, dtype=np.uint64)
        pops = np.zeros(EGENE_KEY_COUNT, dtype=np.uint32)
        cols = np.zeros(EGENE_KEY_COUNT, dtype=np.int32)
        self._lib.evoca_eg_activity_get(
            acts.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            pops.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            cols.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)))
        return {'activity': acts, 'pop_count': pops, 'color': cols}

    def get_eg_food(self):
        """Return egene food-intake table as dict of arrays
        (EGENE_KEY_COUNT = 729 ternary keys). `food` values are uint64
        scaled by 1e6 (divide by 1e6 to get food units)."""
        food = np.zeros(EGENE_KEY_COUNT, dtype=np.uint64)
        pops = np.zeros(EGENE_KEY_COUNT, dtype=np.uint32)
        cols = np.zeros(EGENE_KEY_COUNT, dtype=np.int32)
        self._lib.evoca_eg_food_get(
            food.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            pops.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            cols.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)))
        return {'food': food, 'pop_count': pops, 'color': cols}

    def get_dyn_activity(self):
        """Return per-(LUT input, output) activity table as dict of
        arrays of length DYN_ACT_COUNT = 500. Each bucket counts how
        often the alive population took the corresponding transition.
        Bucket index i = input * 2 + output where input is
        LUT_IDX(v_x, n1, n2, n3) and output ∈ {0, 1}.
        - 'activity'  : cumulative count since sim init (uint64)
        - 'pop_count' : count this step (uint32) — sum = alive cells
                        going into Phase 1
        - 'color'     : ARGB per bucket (int32)"""
        acts = np.zeros(DYN_ACT_COUNT, dtype=np.uint64)
        pops = np.zeros(DYN_ACT_COUNT, dtype=np.uint32)
        cols = np.zeros(DYN_ACT_COUNT, dtype=np.int32)
        self._lib.evoca_dyn_activity_get(
            acts.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            pops.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            cols.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)))
        return {'activity': acts, 'pop_count': pops, 'color': cols}

    def egenome_stats(self):
        """Return scalar Negene-distribution stats over alive cells:
        {'mean_negene', 'std_negene', 'distinct_egene_values',
         'mean_max_match', 'frac_at_max'}."""
        out = np.zeros(5, dtype=np.float32)
        self._lib.evoca_egenome_stats(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
        return {
            'mean_negene':           float(out[0]),
            'std_negene':            float(out[1]),
            'distinct_egene_values': int(out[2]),
            'mean_max_match':        float(out[3]),
            'frac_at_max':           float(out[4]),
        }

    def egene_stats(self):
        """Return scalar cognitive stats over alive cells:
        {'mean_specificity', 'mean_load', 'mean_intake'}.
        - mean_specificity: non-wildcard cell-positions per active egene (0..25).
        - mean_load: per-cell sum of cell-positions across active egenes (0..200).
        - mean_intake: per-eater mouthful this step (0..1)."""
        out = np.zeros(3, dtype=np.float32)
        self._lib.evoca_egene_stats(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
        return {
            'mean_specificity': float(out[0]),
            'mean_load':        float(out[1]),
            'mean_intake':      float(out[2]),
        }

    def get_dissipation(self):
        """Energy-flux / dissipation readout (checkpoint §L, action L1).

        EvoCA is a dissipative structure: food (free energy) is injected
        each step (source = food_inc over regen cells) and destroyed
        through three sinks — (a) tax removed from private food, (b) the
        remaining food of a cell killed by eviction, (c) food_inc lost to
        the F=1 regeneration cap. Eating (mouthful) is an internal
        F_food→f_priv transfer, not a sink.

        Returns a dict (all values in food units):
          Per-step (last sim.step()):
            'source_step', 'diss_step', 'tax_step', 'death_step',
            'overflow_step', 'mouthful_step' (agent intake this step)
          Cumulative (since init):
            'source_cum', 'diss_cum', 'tax_cum', 'death_cum',
            'overflow_cum'
          Instantaneous:
            'alive'         alive count (Φ denominator)
            'phi'           per-capita energy throughput
                            = diss_step / alive (Chaisson energy-rate-
                            density; 0 when alive==0)
            'storage_env'   current Σ F_food
            'storage_priv'  current Σ f_priv
            'storage'       storage_env + storage_priv
            'step'          current step counter

        Whole-system energy identity across one step:
            source_step == d(storage) + diss_step
        (checked in metric-flux-test.ipynb)."""
        out = np.zeros(16, dtype=np.float64)
        self._lib.evoca_get_dissipation(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        storage_env  = float(out[13])
        storage_priv = float(out[14])
        return {
            'source_step':   float(out[0]),
            'diss_step':     float(out[1]),
            'tax_step':      float(out[2]),
            'death_step':    float(out[3]),
            'overflow_step': float(out[4]),
            'source_cum':    float(out[5]),
            'diss_cum':      float(out[6]),
            'tax_cum':       float(out[7]),
            'death_cum':     float(out[8]),
            'overflow_cum':  float(out[9]),
            'alive':         int(out[10]),
            'mouthful_step': float(out[11]),
            'phi':           float(out[12]),
            'storage_env':   storage_env,
            'storage_priv':  storage_priv,
            'storage':       storage_env + storage_priv,
            'step':          int(out[15]),
        }

    def get_egenes_mask(self):
        """Return (N, N, NEGENOME_MAX) uint8 array of egene wildcard masks.
        Bit i of mask byte = 1 means orbit i is non-wildcard."""
        ptr = self._lib.evoca_get_egenes_mask()
        return np.ctypeslib.as_array(
            ptr, shape=(self._N * self._N * 8,)).copy().reshape(
                self._N, self._N, 8)

    def population(self):
        """Total alive cells."""
        return int(self._lib.evoca_get_population())

    def n_distinct_genomes(self):
        """Count of distinct LUT hashes among alive cells (after
        evoca_activity_update has been called this tick)."""
        return int(self._lib.evoca_count_distinct_genomes())

    def get_ages(self):
        """Return (N, N) int32 array of per-cell ages (g_step -
        last_event_step). Dead cells get -1."""
        out = np.zeros(self._N * self._N, dtype=np.int32)
        self._lib.evoca_get_ages(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)))
        return out.reshape(self._N, self._N)

    # ── Opt-in lineage record (research-directions §1) ───────────────
    def enable_lineage(self, on=True):
        """Toggle the opt-in per-cell lineage record on an initialized
        lattice (cells alive now become founders). To record from step 0,
        pass init(..., lineage=True) instead. Calling this before init()
        raises: init() sets the flag from its own `lineage` kwarg, so a
        pre-call would be silently undone. Default OFF -> zero added cost
        in the CA hot path."""
        if not getattr(self, '_N', 0):
            raise RuntimeError(
                "enable_lineage() before init() would be overridden by "
                "init(lineage=False); use init(..., lineage=True)")
        self._lib.evoca_enable_lineage(1 if on else 0)

    def lineage_enabled(self):
        """True if the lineage record is currently enabled."""
        return bool(self._lib.evoca_lineage_enabled())

    def get_lineage_parent_hash(self):
        """Return (N, N) uint32 copy of each cell's recorded parent
        full-genome FNV-1a hash (same hash as the activity tracker).
        0 for founders / cells never born under tracking. All zeros if
        lineage was not enabled."""
        out = np.zeros(self._N * self._N, dtype=np.uint32)
        self._lib.evoca_get_lineage_parent_hash(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        return out.reshape(self._N, self._N)

    def get_lut_hashes(self):
        """Return (N, N) uint32 of each cell's LUT-only FNV-1a hash
        (hash of just the 32 LUT bytes), NOT the combined LUT‖egene
        genome hash used by the activity tracker / lineage_parent_hash.
        Lets the #8 lineage metric decompose parent→child change into
        ΔLUT independently of Δegene. Dead cells hash their zeroed LUT
        to a fixed value; gate on get_alive()."""
        out = np.zeros(self._N * self._N, dtype=np.uint32)
        self._lib.evoca_get_lut_hashes(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        return out.reshape(self._N, self._N)

    def get_lineage_birth_id(self):
        """Return (N, N) uint64 copy of each cell's monotonic birth id
        (0 == cell never born under tracking, e.g. a founder)."""
        out = np.zeros(self._N * self._N, dtype=np.uint64)
        self._lib.evoca_get_lineage_birth_id(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)))
        return out.reshape(self._N, self._N)

    def get_lineage_parent_id(self):
        """Return (N, N) uint64 copy of each cell's parent birth id
        (0 == parent was a founder / born before tracking). Walk
        child.parent_id == parent's birth_id to follow lineage chains."""
        out = np.zeros(self._N * self._N, dtype=np.uint64)
        self._lib.evoca_get_lineage_parent_id(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)))
        return out.reshape(self._N, self._N)

    # ── Per-lineage net metabolic return reducer (checkpoint §D2) ─────
    def lineage_return_update(self):
        """Refresh each lineage's live pop_count and extinction step.
        Call once per step (or periodically) to keep last_step current.
        No-op when lineage is disabled."""
        self._lib.evoca_lineage_return_update()

    def get_lineage_return(self, update=True):
        """Return the per-lineage net-return reducer as a structured
        numpy array, one row per distinct genome-hash lineage seen so far
        (including extinct ones — the table is never compacted).

        Fields: hash (uint32 full-genome FNV-1a), intake, tax, net
        (= intake − tax, all float64 cumulative over the run), births
        (uint64), pop (uint32 current live count), first_step, last_step
        (uint32; last_step is the extinction-step proxy — the most recent
        step the lineage had a live member as of the last update).

        Net return of a lineage = cumulative mouthful eaten by its members
        − cumulative tax paid by its members. Sampling this each step (and
        keying by hash) builds the per-lineage net-return trajectory the
        §D punctuation detector and §H adaptive-vs-drift program need.

        With update=True (default) pop_count/last_step are refreshed first.
        Returns an empty array when lineage is disabled."""
        if update:
            self._lib.evoca_lineage_return_update()
        n = int(self._lib.evoca_lineage_return_count())
        dt = np.dtype([('hash', np.uint32), ('intake', np.float64),
                       ('tax', np.float64), ('net', np.float64),
                       ('births', np.uint64), ('pop', np.uint32),
                       ('first_step', np.uint32), ('last_step', np.uint32)])
        if n == 0:
            return np.empty(0, dtype=dt)
        h  = np.zeros(n, dtype=np.uint32)
        ic = np.zeros(n, dtype=np.float64)
        tx = np.zeros(n, dtype=np.float64)
        bi = np.zeros(n, dtype=np.uint64)
        pp = np.zeros(n, dtype=np.uint32)
        fs = np.zeros(n, dtype=np.uint32)
        ls = np.zeros(n, dtype=np.uint32)
        m = int(self._lib.evoca_get_lineage_return(
            h.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            ic.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            tx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            bi.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            pp.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            fs.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            ls.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))))
        out = np.empty(m, dtype=dt)
        out['hash'] = h[:m]; out['intake'] = ic[:m]; out['tax'] = tx[:m]
        out['net'] = ic[:m] - tx[:m]
        out['births'] = bi[:m]; out['pop'] = pp[:m]
        out['first_step'] = fs[:m]; out['last_step'] = ls[:m]
        return out

    def lineage_return_totals(self):
        """Return (intake_total, tax_total): run-cumulative intake and tax
        summed over all lineages. Independent counters that must equal the
        column sums of get_lineage_return() (a completeness cross-check)."""
        return (float(self._lib.evoca_get_lineage_intake_total()),
                float(self._lib.evoca_get_lineage_tax_total()))

    def get_lineage_intake(self):
        """Return (N, N) float64 of each cell's cumulative intake over its
        *current* life (reset to 0 when the cell is (re)born). All zeros
        when lineage is disabled."""
        out = np.zeros(self._N * self._N, dtype=np.float64)
        self._lib.evoca_get_lineage_intake(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return out.reshape(self._N, self._N)

    def get_lineage_tax(self):
        """Return (N, N) float64 of each cell's cumulative tax paid over
        its *current* life (reset to 0 when the cell is (re)born). All
        zeros when lineage is disabled."""
        out = np.zeros(self._N * self._N, dtype=np.float64)
        self._lib.evoca_get_lineage_tax(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return out.reshape(self._N, self._N)

    def activity_crossings(self, decile_lo=2, decile_hi=3, window=8,
                           max_n=65536):
        """Return activity-wave slope samples across a decile band.

        Reconstructs the last `window` ticks of activity for every live
        G-activity bucket and collects one slope sample per tick whose
        activity increment overlaps [decile_lo..decile_hi) of the current
        live-bucket distribution. Each sample is the bucket's population
        at that tick.
        """
        # Deciles (p10..p90) are normalized by D = distinct live genomes.
        deciles = np.zeros(9, dtype=np.float32)
        self._lib.evoca_q_activity_deciles(
            deciles.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
        D = max(self.n_distinct_genomes(), 1)
        # decile_lo=2 -> p20 (index 1); decile_hi=3 -> p30 (index 2). Clamp.
        lo_i = max(0, min(8, int(decile_lo) - 1))
        hi_i = max(lo_i, min(9, int(decile_hi) - 1))
        a_lo = int(deciles[lo_i] * D) if lo_i < 9 else 0
        a_hi = int(deciles[hi_i] * D) if hi_i < 9 else 0
        if a_hi <= a_lo:
            return np.zeros(0, dtype=np.float32)
        buf = np.zeros(max_n, dtype=np.float32)
        n = self._lib.evoca_activity_crossings(
            ctypes.c_uint64(a_lo), ctypes.c_uint64(a_hi),
            int(window),
            buf.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            int(max_n))
        return buf[:n].copy()

    def activity_flux(self, decile_lo=2, decile_hi=3, window=8):
        """Mean activity flow per tick through the decile band."""
        slopes = self.activity_crossings(decile_lo=decile_lo,
                                         decile_hi=decile_hi, window=window)
        if len(slopes) == 0:
            return 0.0
        return float(slopes.sum()) / float(window)

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

    # ── Effective-rule activity (eff_activity) ───────────────────────

    def eff_activity_update(self):
        """Accumulate the effective-rule activity table for this tick.

        Keys each alive cell on its LUT MASKED to the entries the
        population has exercised (the union of the per-cell used-masks
        accumulated during Phase 1), folded with the same egene stream as
        the full-genome activity key. Two genomes differing only in
        entries never exercised by any cell collapse to one effective-
        rule species, so silent LUT churn does not mint a new species
        (unlike full-genome activity). Scans alive cells. The alive
        eff-species count is <= the alive full-genome species count each
        step (the map is a coarsening of genome space)."""
        self._lib.evoca_eff_activity_update()

    def get_eff_activity(self, max_n=4096):
        """Return the effective-rule activity table as a dict of arrays
        with keys 'hash', 'activity', 'pop_count', 'color'. Parallel to
        get_activity() but over the effective-rule species key. Call
        eff_activity_update() first to populate pop_count/activity."""
        keys = np.zeros(max_n, dtype=np.uint32)
        acts = np.zeros(max_n, dtype=np.uint64)
        pops = np.zeros(max_n, dtype=np.uint32)
        cols = np.zeros(max_n, dtype=np.int32)
        n = self._lib.evoca_eff_activity_get(
            keys.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            acts.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            pops.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            cols.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            max_n)
        return {'hash': keys[:n], 'activity': acts[:n],
                'pop_count': pops[:n], 'color': cols[:n]}

    def eff_count_distinct(self):
        """Count of distinct effective-rule species currently alive
        (after eff_activity_update this tick). Always <= n_distinct
        full-genome species."""
        return int(self._lib.evoca_eff_count_distinct())

    def get_eff_hashes(self):
        """Return (N, N) uint32 of each cell's effective-rule hash: the
        LUT masked by the population-union used-mask, folded with the
        cell's egenes. Flipping a LUT bit never exercised by any cell
        leaves this unchanged (while get_lut_hashes / the full-genome key
        would change). Dead cells hash their zeroed LUT; gate on
        get_alive()."""
        out = np.zeros(self._N * self._N, dtype=np.uint32)
        self._lib.evoca_get_eff_hashes(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        return out.reshape(self._N, self._N)

    def get_eff_used_mask(self):
        """Return (N, N, LUT_BYTES) uint8 of each cell's accumulated
        per-cell used-mask: bit-packed like the LUT, a set bit = that LUT
        entry was exercised by the cell during Phase 1 over its lifetime.
        This is the diagnostic substrate; the eff_activity species key
        uses the population UNION of these (see get_eff_union_mask)."""
        out = np.zeros(self._N * self._N * LUT_BYTES, dtype=np.uint8)
        self._lib.evoca_get_eff_used_mask(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)))
        return out.reshape(self._N, self._N, LUT_BYTES)

    def get_eff_union_mask(self):
        """Return (LUT_BYTES,) uint8 population-union used-mask: the OR of
        every cell's exercised entries since init, bit-packed like the
        LUT. A 0 bit = a LUT entry never exercised by any cell; flipping
        it is behaviourally silent and mints no eff-species. This is the
        mask the eff_activity species key applies."""
        out = np.zeros(LUT_BYTES, dtype=np.uint8)
        self._lib.evoca_get_eff_union_mask(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)))
        return out

    def set_eff_act_ymax(self, y):
        self._lib.evoca_set_eff_act_ymax(int(y))

    def get_eff_act_ymax(self):
        return int(self._lib.evoca_get_eff_act_ymax())

    # ── Per-step demography counters ─────────────────────────────────

    def get_births_last(self):
        return int(self._lib.evoca_get_births_last())

    def get_deaths_last(self):
        return int(self._lib.evoca_get_deaths_last())

    # ── Fixed-space neutral baselines (eg_activity / dyn_activity) ────

    def get_realised_flux(self):
        """Realised per-component mutation flux from the most recent
        step (measured bit-flip counts applied to children, not the
        nominal mu_* metaparameter). Returns a dict:
          'eg_bitflips'  total egene value+mask + presence flips
          'lut_bitflips' total LUT bit flips
          'repro'        reproduction events this step
        Used internally to drive the fixed-space drift baselines."""
        return {
            'eg_bitflips':  int(self._lib.evoca_get_eg_bitflips_last()),
            'lut_bitflips': int(self._lib.evoca_get_lut_bitflips_last()),
            'repro':        int(self._lib.evoca_get_repro_last()),
        }

    def eg_excess_pc(self):
        """Per-component-normalised cumulative excess for eg_activity
        (729 ternary egene keys): observed ΣG/D_G − drift-shadow
        ΣN/D_N. ≈ 0 under neutral egene drift, clearly > 0 under
        selection on the eating genome. Consistent with
        evoca_explore.py's excess_pc_*. Call eg_activity_update()
        first so eg_pop reflects the current step."""
        return float(self._lib.evoca_eg_excess_pc())

    def dyn_excess_pc(self):
        """Per-component-normalised cumulative excess for dyn_activity
        (500 LUT (input,output) buckets): observed ΣG/D_G − drift-
        shadow ΣN/D_N. ≈ 0 under neutral LUT-entry drift, > 0 when
        effective rule transitions are under selection. dyn_pop is
        kept current by Phase 1 every step (no update call needed)."""
        return float(self._lib.evoca_dyn_excess_pc())

    def eg_shadow_total(self):
        """Σ cumulative activity of the eg_activity drift shadow."""
        return float(self._lib.evoca_eg_shadow_total())

    def dyn_shadow_total(self):
        """Σ cumulative activity of the dyn_activity drift shadow."""
        return float(self._lib.evoca_dyn_shadow_total())

    def get_eg_zscores(self):
        """Per-bucket z-scores for the eg_activity fixed space (729
        ternary-egene keys) — the well-posed per-component ("per-gene")
        neutral shadow.  For each key k:

            z[k] = (observed cumulative activity[k] − shadow drift
                    mean[k]) / closed-form WF drift SD[k]

        The drift mean is the shadow's cumulative activity and the drift
        SD is sqrt of the accumulated Wright-Fisher drift variance
        T²·qn(1−qn)/Np.  z ≈ 0 for every key under a null/frozen (zero
        realised-flux) arm; z ≫ 0 on keys the real run concentrates
        activity on faster than fitness-blind drift (selection), z < 0
        where drift stays warmer than selection.  Keys with no
        accumulated drift variance return 0.

        Returns a float64 array of length EGENE_KEY_COUNT (729)."""
        out = np.zeros(EGENE_KEY_COUNT, dtype=np.float64)
        self._lib.evoca_eg_zscores(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return out

    def get_dyn_zscores(self):
        """Per-bucket z-scores for the dyn_activity fixed space (500
        LUT (input,output) transition buckets) — the well-posed per-
        component neutral shadow.  Same definition as get_eg_zscores:

            z[k] = (observed cumulative activity[k] − shadow drift
                    mean[k]) / closed-form WF drift SD[k]

        z ≈ 0 under a null/frozen arm, z ≫ 0 on transitions the real run
        concentrates on faster than drift.  NOTE: dyn z-scores are
        interpretable only where mu_lut>0 — under a frozen rich LUT
        substrate (mu_lut=0) the shadow freezes onto the realised
        composition, giving z ≈ 0 with no selection signal (this does
        avoid the dyn_excess frozen-rich-substrate artefact, but the
        arm carries no dyn signal either).

        Returns a float64 array of length DYN_ACT_COUNT (500)."""
        out = np.zeros(DYN_ACT_COUNT, dtype=np.float64)
        self._lib.evoca_dyn_zscores(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return out

    # ── N-activity probe (Channon-style neutral shadow) ──────────────

    def neutral_enable(self):
        """Enable the neutral shadow population, sized 1:1 with the
        current alive count and seeded by copying current alive cells'
        LUTs. Each subsequent step() mirrors real births/deaths into the
        shadow with random selection (uniform parent / uniform kill,
        same gmu_lut). N-activity counters are reset."""
        self._lib.evoca_neutral_enable()

    def neutral_disable(self):
        """Disable the neutral shadow and free its state."""
        self._lib.evoca_neutral_disable()

    def neutral_is_enabled(self):
        return bool(self._lib.evoca_neutral_is_enabled())

    def neutral_population(self):
        """Live shadow population (matches alive count while enabled)."""
        return int(self._lib.evoca_neutral_get_population())

    def n_activity_update(self):
        self._lib.evoca_n_activity_update()

    def get_n_activity(self, max_n=4096):
        """Return N-activity table (shadow). Same key/value layout as
        get_activity — both use the same FNV-1a LUT-content hash."""
        keys = np.zeros(max_n, dtype=np.uint32)
        acts = np.zeros(max_n, dtype=np.uint64)
        pops = np.zeros(max_n, dtype=np.uint32)
        cols = np.zeros(max_n, dtype=np.int32)
        n = self._lib.evoca_n_activity_get(
            keys.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            acts.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            pops.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            cols.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            max_n)
        return {'hash': keys[:n], 'activity': acts[:n],
                'pop_count': pops[:n], 'color': cols[:n]}

    def nq_activity_deciles(self):
        """Return 9-element float array with N-activity deciles p10..p90."""
        out = np.zeros(9, dtype=np.float32)
        self._lib.evoca_nq_activity_deciles(
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
        return out

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

    def export_recipe(self, descriptor, probes=None, colormode=0,
                      workstream=None, runs_dir=None):
        """Export current initialization recipe to
        <runs_dir>/<date>_<descriptor>.evoca. runs_dir defaults to
        Research/<workstream>/Runs; workstream defaults to the
        EVOCA_WORKSTREAM env var, else Research/_unsorted/Runs.

        Returns the filepath written.
        """
        import json, os
        from datetime import datetime
        import evoca_paths

        if runs_dir is None:
            runs_dir = evoca_paths.runs_dir(workstream)
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
                'mu_egene': self.mu_egene,
                'mu_egenome': self.mu_egenome,
                'p_dup_egene': self.p_dup_egene,
                'tax': self.tax,
                'tax_per_egene': self.tax_per_egene,
                'tax_lut': self.tax_lut,
                'tax_ring': self.tax_ring,
                'restricted_mu': self.restricted_mu,
                'food_per': self.food_per,
                'n_chem': self.n_chem,
                'food_cap': self.food_cap,
                'wdiff': self.wdiff,
                'waste_sink': self.waste_sink,
                'waste_disp': self.waste_disp,
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

    # ── Patch transfer (genome-level) ─────────────────────────────────
    #
    # extract_patch / stamp_patch implement the §11 transplant API of
    # Docs/EvoCA_research_directions.md. They move *genomes* (LUT + all
    # egene slots + active masks), never pixel state, so a transferred
    # patch is a recipe re-seeded into a host field, not a frozen
    # configuration. Built purely on existing per-cell accessors: the C
    # getters return live array pointers, and evoca_set_lut() is the
    # per-cell hook that refreshes the genome-hash cache (it folds the
    # current egenes/active into the hash), so stamping writes egene
    # planes first and finalises each cell with set_lut().

    def _patch_bounds(self, r0, c0, side):
        N = self._N
        side = int(side)
        if side <= 0:
            raise ValueError(f"side must be positive, got {side}")
        if side > N:
            raise ValueError(f"side {side} exceeds lattice N={N}")
        r0, c0 = int(r0), int(c0)
        if not (0 <= r0 and r0 + side <= N and 0 <= c0 and c0 + side <= N):
            raise ValueError(
                f"patch r0={r0} c0={c0} side={side} out of bounds for "
                f"N={N} (no periodic wrap; keep the region in-bounds)")
        return r0, c0, side

    def extract_patch(self, r0, c0, side):
        """Snapshot the genomes of the square region whose top-left
        corner is (r0, c0) and whose edge length is `side`.

        Returns a PatchGenomes with all positions normalised to
        patch-local coords (0..side-1). Pure read of the live C arrays;
        the sim is left untouched."""
        r0, c0, side = self._patch_bounds(r0, c0, side)
        N = self._N
        rs, cs = slice(r0, r0 + side), slice(c0, c0 + side)

        lut = np.ctypeslib.as_array(
            self._lib.evoca_get_lut(),
            shape=(N * N * LUT_BYTES,)).reshape(N, N, LUT_BYTES)
        egv = np.ctypeslib.as_array(
            self._lib.evoca_get_egenes(),
            shape=(N * N * 8,)).reshape(N, N, 8)
        egm = np.ctypeslib.as_array(
            self._lib.evoca_get_egenes_mask(),
            shape=(N * N * 8,)).reshape(N, N, 8)
        act = np.ctypeslib.as_array(
            self._lib.evoca_get_active(),
            shape=(N * N,)).reshape(N, N)
        alv = np.ctypeslib.as_array(
            self._lib.evoca_get_alive(),
            shape=(N * N,)).reshape(N, N)

        return PatchGenomes(
            side=side, src_r0=r0, src_c0=c0, src_N=N,
            lut=lut[rs, cs, :].copy(),
            egenes=egv[rs, cs, :].copy(),
            egenes_mask=egm[rs, cs, :].copy(),
            active=act[rs, cs].copy(),
            alive=alv[rs, cs].copy())

    def stamp_patch(self, patch, r0, c0, mode='overwrite'):
        """Write a PatchGenomes into the square region of this sim whose
        top-left corner is (r0, c0).

        mode='overwrite' : every patch cell replaces the host cell.
        mode='invade'    : only patch-alive cells that land on a *dead*
                           host cell are written (propagule into empty
                           space); occupied host cells are left intact.

        Stamped organisms are set alive with default v/f
        (_STAMP_DEFAULT_V / _STAMP_DEFAULT_F); patch-dead cells in
        'overwrite' mode kill the host cell. Returns the number of cells
        actually written."""
        if not isinstance(patch, PatchGenomes):
            raise TypeError("patch must be a PatchGenomes")
        if mode not in ('overwrite', 'invade'):
            raise ValueError(f"mode must be 'overwrite' or 'invade', "
                             f"got {mode!r}")
        side = patch.side
        r0, c0, side = self._patch_bounds(r0, c0, side)
        N = self._N

        lut = np.ctypeslib.as_array(
            self._lib.evoca_get_lut(),
            shape=(N * N * LUT_BYTES,)).reshape(N, N, LUT_BYTES)
        egv = np.ctypeslib.as_array(
            self._lib.evoca_get_egenes(),
            shape=(N * N * 8,)).reshape(N, N, 8)
        egm = np.ctypeslib.as_array(
            self._lib.evoca_get_egenes_mask(),
            shape=(N * N * 8,)).reshape(N, N, 8)
        act = np.ctypeslib.as_array(
            self._lib.evoca_get_active(),
            shape=(N * N,)).reshape(N, N)
        alv = np.ctypeslib.as_array(
            self._lib.evoca_get_alive(),
            shape=(N * N,)).reshape(N, N)
        vv = np.ctypeslib.as_array(
            self._lib.evoca_get_v(),
            shape=(N * N,)).reshape(N, N)
        ff = np.ctypeslib.as_array(
            self._lib.evoca_get_f(),
            shape=(N * N,)).reshape(N, N)

        n_written = 0
        for pr in range(side):
            tr = r0 + pr
            for pc in range(side):
                tc = c0 + pc
                src_alive = bool(patch.alive[pr, pc])
                if mode == 'invade':
                    # Only a live propagule cell may settle, and only on
                    # empty (dead) ground. This is the null-respecting
                    # invasion: it never displaces a resident.
                    if not src_alive or alv[tr, tc]:
                        continue
                    _alive = 1
                else:  # overwrite
                    _alive = 1 if src_alive else 0

                if _alive:
                    egv[tr, tc, :] = patch.egenes[pr, pc, :]
                    egm[tr, tc, :] = patch.egenes_mask[pr, pc, :]
                    act[tr, tc]    = patch.active[pr, pc]
                    vv[tr, tc]     = _STAMP_DEFAULT_V
                    ff[tr, tc]     = _STAMP_DEFAULT_F
                    alv[tr, tc]    = 1
                    # set_lut() last: it rehashes the genome cache using
                    # the egenes/active just written above.
                    self._lib.evoca_set_lut(
                        int(tr) * N + int(tc),
                        patch.lut[pr, pc, :].ctypes.data_as(
                            ctypes.POINTER(ctypes.c_uint8)))
                else:
                    # Patch said this cell is empty: clear the host cell
                    # (mirrors evoca_set_alive's dead-cell zeroing) so
                    # the transplanted region is faithful, then rehash.
                    alv[tr, tc]    = 0
                    vv[tr, tc]     = 0
                    ff[tr, tc]     = 0.0
                    act[tr, tc]    = 0
                    lut[tr, tc, :] = 0
                    self._lib.evoca_set_lut(
                        int(tr) * N + int(tc),
                        lut[tr, tc, :].ctypes.data_as(
                            ctypes.POINTER(ctypes.c_uint8)))
                n_written += 1
        # f was written straight into f_priv above; rebuild the chem
        # store from it, or the next tax/eat would restore the host
        # cell's old private food into the stamped organism.
        self._lib.evoca_chem_sync_from_f()
        return n_written

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
    'lut':            "'gol' (GoL rule for all) | 'random' (random LUT)",
    'lut_n_init':     "int 1-3: ring depth for random LUT (1=10 bits, 2=50, 3=all 250). Default 3",
    'lut_density':    "float 0-1: probability each independent LUT bit is 1 (lut='random'). Default 0.5",
    'lut_uniform':    "bool: lut='random' shares one rule across all cells if True, else independent per-cell. Default False",
    'lut_seed':       "int|None: numpy seed for the random LUT (None = global RNG). Default None",
    'lut_to_dead':    "bool: if True, re-apply the chosen LUT (and egenome) to all cells AFTER the alive setter has zeroed dead cells. Gives a uniform CA substrate regardless of alive mask — dead cells run the rule too (their v_curr can transition). Default False",
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

def import_run(filepath=None, recipe='init', lib_path=None, N=None):
    """Load a .evoca recipe file and return (sim, display_kwargs).

    If called with no arguments, lists available recipes in
    Research/*/Runs/. A bare filename (no directory) is looked up there.

    recipe: 'init' uses the initial metaparams, 'final' uses the
            metaparams at export time (after any slider adjustments).
    N:      override the grid size from the recipe. Useful for promoting
            a small-grid scan winner to a larger interactive run, e.g.
            sim, kw = import_run(path, N=512). All other init parameters
            (lut, egenome, density, alive layout, food, metaparams) are
            unchanged.

    Usage:
        import_run()                    # list available recipes
        sim, kw = import_run('2026-03-15_my_run.evoca')   # bare name: searched in Research/*/Runs
        sim, kw = import_run(path, N=512)
        run_with_controls(sim, **kw)
    """
    import json
    from pathlib import Path

    import evoca_paths
    if filepath is None:
        recipes = evoca_paths.all_runs()
        if not recipes:
            print("No .evoca files in Research/*/Runs/.")
            return []
        root = evoca_paths.RESEARCH
        for r in recipes:
            print(os.path.relpath(r, root))
        return recipes
    if not os.path.exists(filepath) and os.sep not in str(filepath):
        hits = [r for r in evoca_paths.all_runs()
                if os.path.basename(r) == filepath]
        if len(hits) != 1:
            raise FileNotFoundError(
                f"{filepath!r}: {len(hits)} matches in Research/*/Runs/")
        filepath = hits[0]

    with open(filepath) as f:
        data = json.load(f)

    # version 1 had a single 'metaparams' key; version 2 has init/final
    if 'metaparams_init' in data:
        mp = data['metaparams_init'] if recipe == 'init' else data['metaparams_final']
    else:
        mp = data['metaparams']  # v1 fallback

    mp = correct_legacy_mu_egenome(dict(mp))
    _food_cap_legacy_notice(mp, filepath)

    sim = EvoCA(lib_path=lib_path)
    grid_N = int(N) if N is not None else int(data['N'])
    sim.init(grid_N, **mp)
    # Tag origin so the probe logger can label log files.
    sim._origin_recipe = os.path.basename(filepath)

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


def correct_legacy_mu_egenome(mp):
    """Load-time correction for the 2026-04-28 rename (4bd3fd9): the per-bit
    egene mutation rate was called `mu_egenome` before it, and that name now
    means the Negene slot-growth rate. A metaparam dict with `mu_egenome` but
    no `mu_egene` predates the rename (rows of the 2026-04-27 scans, recipes
    exported from them; export_recipe always writes both keys), so its value
    is mu_egene. Corrects in place; returns mp. Shared by import_run and
    evoca_explore.correct_scan_row_params."""
    if 'mu_egenome' in mp and 'mu_egene' not in mp:
        mp['mu_egene'] = mp.pop('mu_egenome')
    return mp


def _food_cap_legacy_notice(mp, filepath):
    """Recipes written before the A9 fix (2026-10-15) carry no `food_cap`
    key. They are imported with the standard (capped) rule; when m_scale > 1
    the dynamics can differ from the run that produced the recipe, so say
    how to get the legacy behavior back. Returns True if it printed."""
    if 'food_cap' in mp:
        return False
    try:
        m = float(mp.get('m_scale', 1.0))
    except (TypeError, ValueError):
        return False
    if m > 1.0:
        print(f"[food_cap] {os.path.basename(str(filepath))}: no food_cap key "
              f"and m_scale={m:g} > 1; running the standard capped rule. The "
              f"legacy uncapped rule (pre-A9 fix) is food_cap=0.")
        return True
    return False


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


def make_random_lut(n_init=3, density=0.5, seed=None):
    """Build one bit-packed LUT with the given ring-restriction and
    bit-density. Returns a LUT_BYTES-length uint8 array, suitable for
    sim.set_lut_all(lut).

    n_init: 1 (n1 only), 2 (n1+n2), 3 (full n1+n2+n3).
    density: probability that each independent bit is 1.
    seed: numpy seed for reproducibility (None = global RNG).

    Use to inspect or reproduce a particular substrate independently
    of the simulation:

        from python.evoca_py import make_random_lut, unpack_lut
        lut = make_random_lut(n_init=2, density=0.4, seed=7)
        sim.set_lut_all(lut)
        # later: print/inspect  unpack_lut(lut)
    """
    rng = np.random.default_rng(seed) if seed is not None else np.random
    if n_init >= 3:
        bits = (rng.random((2, 5, 5, 5)) < density).astype(np.uint8)
    elif n_init == 2:
        core = (rng.random((2, 5, 5, 1)) < density).astype(np.uint8)
        bits = np.tile(core, (1, 1, 1, 5))
    else:
        core = (rng.random((2, 5, 1, 1)) < density).astype(np.uint8)
        bits = np.tile(core, (1, 1, 5, 5))
    return pack_lut(bits.flatten())


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
