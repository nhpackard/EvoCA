"""Synthetic drift-shadow runs for calibration and power (Docs/shadows.md
§7 (a), (b)).

A synthetic run drives the C drift core (evoca_fs_drift_synth_*) with an
explicit event stream and no lattice: a constant population of Np
individuals, b deaths and b births per step, per-birth realized mutation
counts drawn here from Poisson at the chosen per-bit rates, and fixed
query weights n_i = Np·w_i. R+1 replicates are run; replicate 0 plays the
real run (neutral, or with imposed selection / without-replacement
parents), replicates 1..R are the shadows. The result has the
EvoCA.fs_drift_get() layout, so python/drift_stats.py scores it exactly
as it scores an EvoCA run.

Founders are monoclonal (the run_sim founder): LUT all zero, every egene
slot (value 3, mask 0x3F), slot 0 active, marker 0.
"""
import ctypes
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from evoca_py import EvoCA   # noqa: E402

_LIB = None


def _lib():
    global _LIB
    if _LIB is None:
        _LIB = EvoCA()
    return _LIB._lib


def weights(kind, L=250):
    """Query-mass shares over inputs: 'uniform', 'geo0.7', 'geo0.9'."""
    if kind == 'uniform':
        w = np.ones(L)
    elif kind.startswith('geo'):
        w = float(kind[3:]) ** np.arange(L)
    else:
        raise ValueError(kind)
    return w / w.sum()


def founder(Np):
    lut = np.zeros((Np, 32), np.uint8)
    eg = np.full((Np, 8), 3, np.uint8)
    egm = np.full((Np, 8), 0x3F, np.uint8)
    act = np.ones(Np, np.uint8)
    mk = np.zeros(Np, np.uint16)
    return lut, eg, egm, act, mk


def event_stream(rng, steps, b, mu_lut=0.0, mu_act=0.0, mu_eg=0.0, mu_mk=0.0):
    """Per-birth realized counts (nf, na, ne, nm) for steps*b births."""
    n = steps * b
    recs = np.zeros((n, 4), np.int32)
    recs[:, 0] = rng.poisson(mu_lut * 250, n)
    recs[:, 1] = rng.poisson(mu_act * 8, n)
    recs[:, 2] = rng.poisson(mu_eg * 96, n)
    recs[:, 3] = rng.poisson(mu_mk * 12, n)
    return recs


def synth_run(seed, R=8, Np=500, b=25, steps=1500, mu_lut=0.0, mu_act=0.0,
              mu_eg=0.0, mu_mk=0.0, w='geo0.7', sel_mode=0, s=0.0, loci=(),
              p_dup=1.0):
    """One synthetic run. Returns an fs_drift_get()-layout dict with
    replicate 0 as 'real' and replicates 1..R as shadows.

    sel_mode: 0 neutral; 1 dyn purifying weight (1-s)^(# allele-1 at
    `loci`); 2 eg advantage (1+s) for an active egene with orbit-0 digit 2;
    3 neutral with without-replacement parents."""
    L = _lib()
    rng = np.random.default_rng(seed)
    lut, eg, egm, act, mk = founder(Np)
    U8 = ctypes.POINTER(ctypes.c_uint8)
    L.evoca_set_p_dup_egene(ctypes.c_float(p_dup))
    rc = L.evoca_fs_drift_synth_init(
        R + 1, Np + b + 1, Np, lut.ctypes.data_as(U8), eg.ctypes.data_as(U8),
        egm.ctypes.data_as(U8), act.ctypes.data_as(U8),
        mk.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16)),
        int(rng.integers(1, 2**32 - 1)))
    if rc != R + 1:
        raise RuntimeError(f"synth_init returned {rc}")
    loci = np.asarray(loci, dtype=np.int32)
    CI = ctypes.POINTER(ctypes.c_int)
    L.evoca_fs_drift_synth_set_selection(int(sel_mode), float(s),
                                         loci.ctypes.data_as(CI), int(loci.size))
    D = np.full(steps, b, np.int32)
    B = np.full(steps, b, np.int32)
    recs = event_stream(rng, steps, b, mu_lut, mu_act, mu_eg, mu_mk)
    nq = np.ascontiguousarray(np.tile(Np * weights(w), (steps, 1)))
    rc = L.evoca_fs_drift_synth_run(steps, D.ctypes.data_as(CI),
                                    B.ctypes.data_as(CI),
                                    recs.ctypes.data_as(CI),
                                    nq.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                    1)
    if rc != 0:
        raise RuntimeError(f"synth_run returned {rc}")
    P = ctypes.POINTER(ctypes.c_double)
    Rt = R + 1
    eg_sh = np.zeros((Rt, 729)); dyn_sh = np.zeros((Rt, 500))
    mk_sh = np.zeros((Rt, 729)); neg = np.zeros(Rt + 1)
    L.evoca_fs_drift_get(None, eg_sh.ctypes.data_as(P), None, None,
                         dyn_sh.ctypes.data_as(P), None, mk_sh.ctypes.data_as(P),
                         neg.ctypes.data_as(P))
    steps_done = int(L.evoca_fs_drift_steps())
    L.evoca_fs_drift_disable()
    return {'eg_real': eg_sh[0], 'eg_shadow': eg_sh[1:],
            'dyn_real': dyn_sh[0], 'dyn_queried': None, 'dyn_shadow': dyn_sh[1:],
            'marker_real': mk_sh[0], 'marker_shadow': mk_sh[1:],
            'negene_real': neg[1], 'negene_shadow': neg[2:],
            'steps': steps_done, 'R': R}
