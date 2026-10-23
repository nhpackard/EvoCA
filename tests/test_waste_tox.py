"""waste_tox, the toxicity threshold of the food/waste chemistry, is a
metaparameter: it reaches C through init() and update_waste_tox(), is
recorded by params() and recipes, and a recipe without it loads the
previous fixed value 1.0."""
import json
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'python'))

from evoca_py import EvoCA, import_run      # noqa: E402

P = dict(food_inc=0.0, m_scale=0.0, gdiff=0.0, n_chem=2, waste_disp=0.0)


def _tox_c(s):
    return s.get_chem_params(EvoCA.WASTE)[4]


def _poisoned(make_sim, **kw):
    """All alive, private waste 0.35, no excretion, no tax: who dies?"""
    s = make_sim(8, **P, **kw)
    s.state(lut='gol', egenome='uniform', egenome_value=3)
    s.set_F_all(0.0); s.set_alive_all(); s.set_f_all(0.5)
    s.set_w(0.35)
    s.step()
    return s


def test_default_is_one_and_reaches_c(make_sim):
    s = make_sim(8, **P)
    assert s.waste_tox == 1.0 and _tox_c(s) == pytest.approx(1.0)
    assert _poisoned(make_sim).get_alive().all()          # 0.35 < 1: survive


def test_init_value_reaches_c_and_kills(make_sim):
    s = _poisoned(make_sim, waste_tox=0.3)
    assert _tox_c(s) == pytest.approx(0.3)
    assert not s.get_alive().any()                        # 0.35 >= 0.3: all die
    assert s.get_tox_deaths_last() == 64


def test_update_reaches_c(make_sim):
    s = make_sim(8, **P)
    s.update_waste_tox(0.25)
    assert s.waste_tox == 0.25 and _tox_c(s) == pytest.approx(0.25)
    with pytest.raises(ValueError):
        s.update_waste_tox(0.0)


def test_refused_without_food_waste():
    s = EvoCA()
    with pytest.raises(ValueError, match="waste_tox"):
        s.init(8, waste_tox=0.5)
    with pytest.raises(ValueError):
        s.init(8, n_chem=2, waste_tox=-1.0)


def test_params_and_recipe_round_trip(make_sim, tmp_path):
    s = make_sim(16, **P, waste_tox=0.4)
    s.state(lut='gol')
    assert s.params()['waste_tox'] == 0.4
    s.update_waste_tox(0.2)
    path = s.export_recipe('tox', runs_dir=str(tmp_path))
    rec = json.load(open(path))
    assert rec['metaparams_init']['waste_tox'] == 0.4
    assert rec['metaparams_final']['waste_tox'] == 0.2
    for which, want in (('init', 0.4), ('final', 0.2)):
        r, _ = import_run(path, recipe=which)
        assert r.waste_tox == want and _tox_c(r) == pytest.approx(want)
        r.free()
    # a recipe written before waste_tox existed loads the old fixed value 1.0
    for k in ('metaparams_init', 'metaparams_final'):
        rec[k].pop('waste_tox')
    old = tmp_path / 'old.evoca'
    old.write_text(json.dumps(rec))
    r, _ = import_run(str(old))
    assert r.waste_tox == 1.0 and _tox_c(r) == pytest.approx(1.0)
    r.free()
