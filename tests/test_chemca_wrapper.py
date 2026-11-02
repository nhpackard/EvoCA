"""ChemCA Python-wrapper defects from the 2026-11-02 ultrareview of chemca:
stale light field, light field lost on reload, N committed before C accepts
it, unchecked array sizes, hidden gcc errors, params() round trip."""
import os
import shutil
import subprocess

import numpy as np
import pytest

import evoca_build
from chemca_py import ChemCA, build_chemca, load_recipe


def _fresh(N, seed=3, **params):
    s = ChemCA()
    s.init(N, **params)
    s.state(seed=seed)
    return s


# ── light field ──────────────────────────────────────────────────────

def test_reinit_same_N_drops_old_light_field():
    """A field set before init() must not come back on the next update()."""
    a = ChemCA()
    a.init(8)
    a.set_light_field(np.zeros((8, 8)))
    a.init(8)
    a.state(seed=3)
    a.update(light=0.05)
    a.step(20)
    assert a.recipe()['light_field'] is False
    b = _fresh(8)
    b.update(light=0.05)
    b.step(20)
    assert a.state_hash() == b.state_hash()


def test_reinit_other_N_then_update_does_not_raise():
    s = ChemCA()
    s.init(8)
    s.set_light_field(np.zeros((8, 8)))
    s.free()
    s.init(10)
    s.state(seed=1)
    s.update(light=0.05)                     # raised: 8x8 field reshaped to 10x10
    assert s.recipe()['light_field'] is False


def test_recipe_with_light_field_refuses_to_load(tmp_path):
    s = _fresh(8)
    s.set_light_field(np.full((8, 8), 0.01))
    p = s.export_recipe('lf', runs_dir=str(tmp_path))
    with pytest.raises(ValueError, match='light field'):
        load_recipe(p)


# ── failed init keeps the old world consistent ──────────────────────

def test_failed_init_keeps_old_N():
    s = _fresh(8)
    h = s.state_hash()
    with pytest.raises(RuntimeError):
        s.init(2)                            # chemca_init rejects N < 3
    assert s.N == 8
    assert s.get_E(0).shape == (8, 8)        # was a buffer sized for N=2
    assert s.state_hash() == h


# ── array sizes ──────────────────────────────────────────────────────

def test_set_alive_wrong_size_raises():
    s = _fresh(8)
    with pytest.raises(ValueError):
        s.set_alive(np.ones((4, 4)))


def test_set_genes_wrong_size_raises():
    s = _fresh(8)
    g, a = s.get_rgenes()
    with pytest.raises(ValueError):
        s.set_genes(g[:4], a)
    with pytest.raises(ValueError):
        s.set_genes(g, a[:4])


# ── params() round trip ──────────────────────────────────────────────

def test_params_round_trip_through_build_chemca():
    s = _fresh(8, light=0.07)
    t = build_chemca(s.params(), N=8, seed=3)
    assert t.params() == s.params()


# ── build errors are visible ─────────────────────────────────────────

def test_build_error_shows_gcc_stderr(tmp_path):
    shutil.copytree(os.path.join(evoca_build.ROOT, 'C'), tmp_path / 'C',
                    ignore=shutil.ignore_patterns('*.dylib', '*.so'))
    with open(tmp_path / 'C' / 'core.c', 'a') as f:
        f.write('\nint broken_on_purpose( {\n')
    with pytest.raises(subprocess.CalledProcessError) as ei:
        evoca_build.build(root=str(tmp_path), force=True)
    assert 'broken_on_purpose' in str(ei.value)
