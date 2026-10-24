"""Model dispatch (sim_model): recipes, build_sim, make_sim, load_recipe.

A params dict or recipe without 'model' is EvoCA (load-time default).
Each door is entered as a caller would, and asserted on what reached C.
"""
import glob
import json
import os

import numpy as np
import pytest

import evoca_build
import evoca_explore as ee
import evoca_nb
import evoca_py
from chemca_py import ChemCA, build_chemca
from evoca_py import EvoCA

ROOT = evoca_build.ROOT


def test_recipe_round_trip_identical_hash(tmp_path):
    p = dict(model='chemca', preset='gs3', mu_rgene=0.05, mu_ractive=0.05, light=0.025)
    sim = build_chemca(p, N=32, seed=9)
    path = sim.export_recipe('rt', runs_dir=str(tmp_path))
    sim.step(200)
    h1 = sim.state_hash()
    sim.free()
    data = json.load(open(path))
    assert data['model'] == 'chemca' and data['initialization']['seed'] == 9
    sim2, _ = evoca_py.import_run(path)
    assert isinstance(sim2, ChemCA)
    sim2.step(200)
    assert sim2.state_hash() == h1
    sim2.free()
    sim3, _ = evoca_nb.load_recipe(path)               # notebook door
    sim3.step(200)
    assert sim3.state_hash() == h1
    sim3.free()


def test_recipe_final_params(tmp_path):
    sim = build_chemca(dict(model='chemca', preset='fw2'), N=16, seed=1)
    sim.update(light=0.04)
    path = sim.export_recipe('fin', runs_dir=str(tmp_path))
    sim.free()
    s2, _ = evoca_py.import_run(path, recipe='final')
    rx = s2._lib
    import ctypes
    A, B, r, f = ctypes.c_int(), ctypes.c_int(), ctypes.c_float(), ctypes.c_int()
    rx.chemca_env_rx_get(0, ctypes.byref(A), ctypes.byref(B), ctypes.byref(r), ctypes.byref(f))
    assert r.value == pytest.approx(0.04)
    s2.free()


def test_evoca_recipe_without_model_is_evoca():
    path = os.path.join(ROOT, 'Research/Chem/Runs/2026-10-23_chem_fw_spiralpoint.evoca')
    assert 'model' not in json.load(open(path))
    sim, _ = evoca_py.import_run(path, N=32)
    assert isinstance(sim, EvoCA) and sim.n_chem == 2
    sim.free()


def test_build_sim_dispatch():
    s = ee.build_sim(dict(model='chemca', preset='gs3', mu_rgene=0.01, mbar=1.5), N=16, seed=1)
    assert isinstance(s, ChemCA)
    assert s._lib.chemca_get_mu_rgene() == pytest.approx(0.01)       # reached C
    assert s.mass()['total'] / 16 ** 2 == pytest.approx(1.5, rel=1e-6)  # state key reached
    s.free()
    e = ee.build_sim(dict(food_inc=0.01, tax=0.02), N=16, seed=1)
    assert isinstance(e, EvoCA) and e.tax == pytest.approx(0.02)
    e.free()
    e = ee.build_sim(dict(model='evoca', tax=0.03), N=16, seed=1)     # explicit evoca
    assert isinstance(e, EvoCA) and e.tax == pytest.approx(0.03)
    e.free()
    with pytest.raises(ValueError):
        ee.build_sim(dict(model='bugs'), N=16)


def test_run_sim_refuses_chemca():
    with pytest.raises(ValueError, match='run_chemca'):
        ee.run_sim(dict(model='chemca'), n_steps=1, N=8)


def test_make_sim_dispatch():
    import ctypes
    s = evoca_nb.make_sim(16, seed=3, model='chemca', preset='fw2', light=0.03)
    assert isinstance(s, ChemCA)
    A, B, r, f = ctypes.c_int(), ctypes.c_int(), ctypes.c_float(), ctypes.c_int()
    s._lib.chemca_env_rx_get(0, ctypes.byref(A), ctypes.byref(B), ctypes.byref(r), ctypes.byref(f))
    assert r.value == pytest.approx(0.03) and (A.value, B.value) == (1, 0)   # W -> F
    s.free()
    e = evoca_nb.make_sim(16, seed=3, tax=0.02)
    assert isinstance(e, EvoCA)
    e.free()


def test_unknown_chemca_param_raises():
    with pytest.raises(ValueError, match='unknown'):
        ChemCA().init(8, preset='fw2', tox_P=0.3)          # P is not an fw2 compound


def test_build_statement_is_consistent():
    """Makefile, CLAUDE.md and python/evoca_build.py name the same sources."""
    srcs = ' '.join(evoca_build.SOURCES)
    mk = open(os.path.join(ROOT, 'Makefile')).read()
    assert f'SRCS   = {srcs}' in mk
    cmd = f"gcc -O2 -Wall -fPIC -dynamiclib -o C/libevoca.dylib {srcs}"
    assert cmd in mk and cmd in open(os.path.join(ROOT, 'CLAUDE.md')).read()
