"""Tests for the load-time correction of the 2026-04-28 rename mu_egenome -> mu_egene.

Before 4bd3fd9 the per-bit egene mutation rate was called `mu_egenome`;
the name now means the Negene slot-growth rate. The 2026-04-27 scans
(initial, refined, push_edges, evo_focus) record `param_mu_egenome` with
the old meaning, and recipes exported from them carry `mu_egenome` and no
`mu_egene`. Without the correction those load with mu_egene=0 and a
nonzero Negene growth rate (found 2026-10-21: push_edges cfg246's
population waves vanish).

Guards, each through a public door:
  1. evoca_from_scan on a pre-rename row writes mu_egene, not mu_egenome.
  2. import_run of a pre-rename recipe (2026-04-29_top_evo_*) sets
     sim.mu_egene to the recorded value and leaves sim.mu_egenome at 0.
  3. evoca_nb.params_from_row on a pre-rename row gives mu_egene.
  4. A post-rename recipe (both keys) is left as written.
"""
import json
import os
import sys

import pandas as pd
import pytest

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))

RD = os.path.join(_ROOT, 'Research', 'Reaction-diffusion')
PUSH = os.path.join(RD, 'Scans', '2026-04-27_push_edges')
OLD_RECIPE = os.path.join(RD, 'Runs', '2026-04-29_top_evo_15_cfg71_score0.66.evoca')


def test_evoca_from_scan_pre_rename_row_writes_mu_egene(tmp_path):
    import evoca_explore as ee
    p = ee.evoca_from_scan(PUSH, 246, runs_dir=str(tmp_path))
    mp = json.load(open(p))['metaparams_init']
    assert mp.get('mu_egene') == pytest.approx(0.01)
    assert 'mu_egenome' not in mp


def test_import_run_pre_rename_recipe_sets_mu_egene():
    from evoca_py import import_run
    raw = json.load(open(OLD_RECIPE))['metaparams_init']
    assert 'mu_egenome' in raw and 'mu_egene' not in raw   # the fixture is pre-rename
    sim, _ = import_run(OLD_RECIPE, N=16)
    try:
        assert sim.mu_egene == pytest.approx(raw['mu_egenome'])
        assert sim.mu_egenome == 0.0
    finally:
        sim.free()


def test_params_from_row_pre_rename_row_gives_mu_egene():
    from evoca_nb import params_from_row
    row = pd.read_csv(os.path.join(PUSH, 'results.csv')).set_index('config_idx').loc[246]
    p = params_from_row(row)
    assert p.get('mu_egene') == pytest.approx(0.01)
    assert 'mu_egenome' not in p


def test_post_rename_recipe_untouched(tmp_path):
    from evoca_py import EvoCA, import_run
    s = EvoCA(); s.init(16, mu_egene=0.002, mu_egenome=0.005); s.state()
    p = s.export_recipe('both_keys', runs_dir=str(tmp_path)); s.free()
    sim, _ = import_run(p)
    try:
        assert sim.mu_egene == pytest.approx(0.002)
        assert sim.mu_egenome == pytest.approx(0.005)
    finally:
        sim.free()
