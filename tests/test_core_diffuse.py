"""Shared-core diffusion (C/core.c core_diffuse, CHEM-P1).

The core extraction replaced evoca.c's three box-blur copies (food full
pass, food partial pass, compound-field pass) by one function. The older
goldens (chem_parity, metrics regression) all run at diffusion < 1, so they
guard only the partial pass; a mutation of the full-pass branch passed the
whole suite. These tests close that gap:

  * golden trajectory hashes at gdiff >= 1 and wdiff >= 1, computed with the
    pre-extraction library (main @ a9854b4, C/evoca.c alone);
  * the C function, called directly, against a float32 reference: a full
    pass is exactly sum/9, a partial pass matches to rounding, mass is
    conserved, and the buffers swap.
"""
import ctypes
import hashlib

import numpy as np
import pytest

from evoca_py import EvoCA, _find_lib

# (params, sha256[:16]) from the pre-extraction library; see module docstring.
GOLDEN = [
    (dict(food_inc=0.02, m_scale=1.0, gdiff=2.37, mu_lut=0.02, mu_egene=0.01,
          mu_egenome=0.01, tax=0.02, restricted_mu=1), 'be03556fea696559'),
    (dict(food_inc=0.016, m_scale=3.5, gdiff=1.0, mu_lut=0.03, mu_egene=0.003,
          mu_egenome=0.003, tax=0.035, n_chem=2, wdiff=2.0, waste_sink=0.015,
          waste_disp=0.5, waste_exchange=1), 'd26edd5a7b6b3592'),
]


def trajectory_hash(params, lib_path=None, N=32, steps=150, seed=11):
    """sha256[:16] of alive, F, f, egenes (and W, w at n_chem=2) every 25 steps."""
    s = EvoCA(lib_path=lib_path)
    s.init(N, **params)
    s.set_seed(seed)
    s.state(lut='gol', egenome='uniform', egenome_value=3, v_seed=seed,
            f_init=0.5, F_init=1.0, alive='halfplane')
    h = hashlib.sha256()
    for t in range(steps):
        s.step()
        if t % 25 == 24:
            arrs = [s.get_alive(), s.get_F(), s.get_f(), s.get_egenes()]
            if params.get('n_chem', 1) == 2:
                arrs += [s.get_W(), s.get_w()]
            for a in arrs:
                h.update(np.ascontiguousarray(a).tobytes())
    s.free()
    return h.hexdigest()[:16]


@pytest.mark.parametrize('params,golden', GOLDEN)
def test_full_pass_diffusion_golden(params, golden):
    assert trajectory_hash(params) == golden


def _core():
    lib = ctypes.CDLL(_find_lib())
    fp = ctypes.POINTER(ctypes.c_float)
    lib.core_diffuse.argtypes = [ctypes.POINTER(fp), ctypes.POINTER(fp),
                                 ctypes.c_int, ctypes.c_float]
    lib.core_diffuse.restype = None
    return lib, fp


def _reference(F, d):
    """float32 reference with the C summation order (di, dj row-major)."""
    F = F.astype(np.float32)
    n_full = int(d)
    passes = [np.float32(1.0)] * n_full
    frac = np.float32(d) - np.float32(n_full)
    if frac > 0:
        passes.append(frac)
    for s in passes:
        acc = np.zeros_like(F)
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                acc = (acc + np.roll(np.roll(F, -di, 0), -dj, 1)).astype(np.float32)
        t = np.float32(1.0) - s
        w9 = s * np.float32(1.0 / 9.0)
        F = (t * F + w9 * acc).astype(np.float32)
    return F


def _run_core(F, d):
    lib, fp = _core()
    N = F.shape[0]
    a = np.ascontiguousarray(F, np.float32).copy()
    b = np.zeros_like(a)
    pa, pb = a.ctypes.data_as(fp), b.ctypes.data_as(fp)
    lib.core_diffuse(ctypes.byref(pa), ctypes.byref(pb), N, d)
    out_addr = ctypes.cast(pa, ctypes.c_void_p).value
    out = a if out_addr == a.ctypes.data else b
    return out.copy(), out_addr, a.ctypes.data, b.ctypes.data


def test_core_diffuse_full_pass_is_exact_mean():
    F = np.random.default_rng(3).random((17, 17)).astype(np.float32)
    out, *_ = _run_core(F, 2.0)
    assert np.array_equal(out, _reference(F, 2.0))


def test_core_diffuse_partial_and_mass():
    F = np.random.default_rng(4).random((16, 16)).astype(np.float32)
    for d in (0.3, 1.7):
        out, *_ = _run_core(F, d)
        assert np.allclose(out, _reference(F, d), rtol=0, atol=1e-6)
        assert abs(out.sum(dtype=np.float64) - F.sum(dtype=np.float64)) < 1e-4


def test_core_diffuse_swaps_and_noop():
    F = np.random.default_rng(5).random((8, 8)).astype(np.float32)
    out, addr, a, b = _run_core(F, 1.0)        # one pass: result in the scratch
    assert addr == b
    out, addr, a, b = _run_core(F, 2.0)        # two passes: back in the field
    assert addr == a
    out, addr, a, b = _run_core(F, 0.0)        # no-op
    assert addr == a and np.array_equal(out, F)
