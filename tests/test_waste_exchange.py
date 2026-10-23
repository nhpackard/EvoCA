"""Waste exchange modes (CHEM-D4). waste_exchange=0 is the original one-way
excretion d = disp * w * max(0, 1 - W) and stays the default, bit for bit;
waste_exchange=1 is two-way passive membrane transport d = disp * (w - W)
(d > 0 excretes, d < 0 takes up), capped by the source. Constructed
fixtures; every capability entered through the public API."""
import hashlib
import json
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'python'))

from evoca_py import EvoCA, import_run      # noqa: E402

DISP = 0.5


def _still(make_sim, mode, N=8, **kw):
    """All alive, no regrowth, diffusion, eating, mutation or tax."""
    p = dict(food_inc=0.0, m_scale=0.0, gdiff=0.0, n_chem=2, waste_disp=DISP,
             waste_exchange=mode)
    p.update(kw)
    s = make_sim(N, **p)
    s.state(lut='gol', egenome='uniform', egenome_value=3)
    s.set_F_all(0.0); s.set_alive_all(); s.set_f_all(0.5)
    return s


def _env(s):
    return np.ctypeslib.as_array(s._lib.evoca_get_chem_env(EvoCA.WASTE), shape=(s.N * s.N,))


def _fixture(make_sim, mode):
    """Private w = 0.3 everywhere; site W a ramp 0..2 (both directions, uptake to toxicity and
    the one-way room cutoff W >= 1 are exercised)."""
    s = _still(make_sim, mode)
    s.set_w(0.3)
    E = _env(s); E[:] = np.linspace(0.0, 2.0, 64, dtype=np.float32)
    return s, E.copy()


def test_default_is_one_way(make_sim):
    s = make_sim(8, n_chem=2)
    assert s.waste_exchange == 0
    assert s._lib.evoca_get_chem_exchange_mode(EvoCA.WASTE) == 0


def test_one_way_step(make_sim):
    s, E0 = _fixture(make_sim, 0)
    s.step()
    w0 = np.float32(0.3)
    d = np.where(E0 < 1, DISP * w0 * (1 - E0), 0).astype(np.float32)
    assert np.allclose(s.get_w().ravel(), w0 - d, atol=1e-6)
    assert np.allclose(_env(s), E0 + d, atol=1e-6)


def test_two_way_step(make_sim):
    s, E0 = _fixture(make_sim, 1)
    s.step()
    w0 = np.float32(0.3)
    d = DISP * (w0 - E0)                         # > 0 out, < 0 in; caps do not bind here
    w1 = w0 - d
    alive = s.get_alive().ravel() == 1
    assert np.allclose(s.get_w().ravel()[alive], w1[alive], atol=1e-6)
    assert np.all(w1[~alive] >= 1.0 - 1e-6), "only toxic cells may have died"
    assert np.allclose(_env(s)[alive], (E0 + d)[alive], atol=1e-6)
    assert s.get_tox_deaths_last() == int((w1 >= 1.0).sum()) > 0   # uptake kills at W high


def test_two_way_uptake_into_a_clean_cell(make_sim):
    """A cell holding no waste at all (no store slot) still takes up waste."""
    s = _still(make_sim, 1)
    assert np.all(s.get_w() == 0)
    E = _env(s); E[:] = 0.4
    s.step()
    assert np.allclose(s.get_w(), DISP * 0.4, atol=1e-6)
    assert np.allclose(_env(s), 0.4 - DISP * 0.4, atol=1e-6)


def test_two_way_dead_cells_do_not_exchange(make_sim):
    s = _still(make_sim, 1)
    a = np.ones((8, 8), np.uint8); a[0, :] = 0
    s.set_alive(a)
    E = _env(s); E[:] = 0.4
    s.step()
    assert np.allclose(_env(s).reshape(8, 8)[0], 0.4)
    assert np.all(s.get_w()[0] == 0)


def _fp(s):
    h = hashlib.sha1()
    for x in (s.get_F(), s.get_f(), s.get_W(), s.get_w(), s.get_alive(), s.get_v()):
        h.update(np.ascontiguousarray(x).tobytes())
    return h.hexdigest()


def _evolving(make_sim, **kw):
    N = 32
    rng = np.random.default_rng(3)
    p = dict(food_inc=0.016, m_scale=3.5, gdiff=0.08, tax=0.035, mu_lut=0.03, mu_egene=0.003,
             n_chem=2, wdiff=0.01, waste_sink=0.015, waste_disp=0.5)
    p.update(kw)
    s = make_sim(N, **p); s.set_seed(3)
    s.state(lut='gol', egenome='uniform', egenome_value=3)
    s.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.5); s.set_F_all(1.0); s.set_alive_halfplane(0)
    for _ in range(200):
        s.step()
    return s


def test_explicit_one_way_is_bit_identical_to_default(make_sim):
    assert _fp(_evolving(make_sim)) == _fp(_evolving(make_sim, waste_exchange=0))


def test_two_way_changes_the_run(make_sim):
    assert _fp(_evolving(make_sim)) != _fp(_evolving(make_sim, waste_exchange=1))


def test_update_and_validation(make_sim):
    s = make_sim(8, n_chem=2)
    s.update_waste_exchange(1)
    assert s.waste_exchange == 1 and s._lib.evoca_get_chem_exchange_mode(EvoCA.WASTE) == 1
    with pytest.raises(ValueError):
        s.update_waste_exchange(2)
    with pytest.raises(ValueError, match="waste_exchange"):
        EvoCA().init(8, waste_exchange=1)                 # needs n_chem=2
    with pytest.raises(ValueError):
        EvoCA().init(8, n_chem=2, waste_exchange=0.5)


def test_restart_keeps_two_way(make_sim):
    s = make_sim(8, n_chem=2, waste_exchange=1)
    s.reinit()
    assert s._lib.evoca_get_chem_exchange_mode(EvoCA.WASTE) == 1


def test_recipe_round_trip(make_sim, tmp_path):
    s = make_sim(16, n_chem=2, waste_exchange=1)
    s.state(lut='gol')
    assert s.params()['waste_exchange'] == 1
    path = s.export_recipe('exch', runs_dir=str(tmp_path))
    rec = json.load(open(path))
    assert rec['metaparams_init']['waste_exchange'] == 1
    r, _ = import_run(path)
    assert r.waste_exchange == 1 and r._lib.evoca_get_chem_exchange_mode(EvoCA.WASTE) == 1
    r.free()
    for k in ('metaparams_init', 'metaparams_final'):     # a recipe from before the flag
        rec[k].pop('waste_exchange')
    old = tmp_path / 'old.evoca'; old.write_text(json.dumps(rec))
    r, _ = import_run(str(old))
    assert r.waste_exchange == 0 and r._lib.evoca_get_chem_exchange_mode(EvoCA.WASTE) == 0
    r.free()
