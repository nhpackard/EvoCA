"""EvoCA.state(): the initial v(x) field (2026-10-23).

Found by the Ring-complexity notebook pass: state() drew v from an unseeded
numpy generator, so make_sim + set_seed + state() runs were not reproducible
(the C PRNG does not touch v), and v_density was only recorded, never used
(the draw was always 50/50). Guards: v_seed makes the run reproducible and
round-trips through a recipe; v_density sets the density.
"""
import hashlib
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))


def _run(make_sim, v_seed, steps=50):
    s = make_sim(64, food_inc=0.015, m_scale=2.0, gdiff=0.08, tax=0.02, mu_lut=0.01)
    s.set_seed(3)
    s.state(lut='gol', egenome='uniform', egenome_value=3, alive='halfplane', F_init=1.0, f_init=0.5, v_seed=v_seed)
    for _ in range(steps):
        s.step()
    h = hashlib.sha1()
    for a in (s.get_F(), s.get_v(), s.get_alive()):
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


def test_v_seed_makes_runs_reproducible(make_sim):
    assert _run(make_sim, 7) == _run(make_sim, 7)
    assert _run(make_sim, 7) != _run(make_sim, 8)


@pytest.mark.parametrize('density', [0.1, 0.5, 0.8])
def test_v_density_sets_the_density(make_sim, density):
    s = make_sim(128)
    s.state(v_density=density, v_seed=1)
    assert s.get_v().mean() == pytest.approx(density, abs=0.02)


def test_v_seed_round_trips_through_a_recipe(tmp_path):
    from evoca_py import EvoCA, import_run
    s = EvoCA(); s.init(32, m_scale=1.0)
    s.state(v_seed=11, v_density=0.3)
    v0 = s.get_v().copy()
    path = s.export_recipe('_pytest_vseed', runs_dir=str(tmp_path)); s.free()
    init = json.load(open(path))['initialization']
    assert init['v_seed'] == 11
    sim, _ = import_run(path)
    try:
        assert np.array_equal(sim.get_v(), v0)
    finally:
        sim.free()
