"""Research/<Workstream>/{Runs,Scans} layout, exercised through the public
entry points (import_run, evoca_from_scan, export_recipe) — 2026-10-15."""
import json
import os

import pytest

import evoca_paths
from evoca_explore import evoca_from_scan, _resolve_scan_dir
from evoca_py import EvoCA, import_run


def test_no_top_level_runs_or_scans():
    for d in ('Runs', 'Scans'):
        p = os.path.join(evoca_paths.REPO_ROOT, d)
        assert not os.path.isdir(p) or not any(
            f.endswith(('.evoca', '.csv', '.py')) for _, _, fs in os.walk(p) for f in fs), p


def test_import_run_lists_workstream_recipes(capsys):
    recipes = import_run()
    assert len(recipes) >= 87
    assert all(os.sep + 'Runs' + os.sep in r and r.startswith(evoca_paths.RESEARCH) for r in recipes)


def test_import_run_bare_name_loads():
    sim, kw = import_run('2026-04-27_top_1_cfg236_score0.90.evoca', N=32)
    sim.step()
    sim.free()


def test_every_recipe_source_scan_resolves():
    for r in evoca_paths.all_runs():
        src = json.load(open(r)).get('source_scan')
        if src:
            assert os.path.isfile(os.path.join(evoca_paths.REPO_ROOT, src, 'results.csv')), (r, src)


def test_scan_name_resolves_to_its_workstream():
    d = _resolve_scan_dir('2026-04-27_initial')
    assert d.endswith(os.path.join('Reaction-diffusion', 'Scans', '2026-04-27_initial'))
    assert evoca_paths.runs_dir_for_scan(d) == evoca_paths.runs_dir('Reaction-diffusion')


def test_evoca_from_scan_writes_relative_provenance(tmp_path):
    p = evoca_from_scan('2026-04-27_initial', config_idx=49, descriptor='t', runs_dir=str(tmp_path))
    rec = json.load(open(p))
    assert rec['source_scan'] == 'Research/Reaction-diffusion/Scans/2026-04-27_initial'


def test_export_recipe_workstream_and_env(tmp_path, monkeypatch):
    s = EvoCA(); s.init(16)
    p = s.export_recipe('x', runs_dir=str(tmp_path))
    assert os.path.dirname(p) == str(tmp_path)
    with pytest.raises(ValueError):
        s.export_recipe('x', workstream='Nope')
    monkeypatch.setenv('EVOCA_WORKSTREAM', 'Ring-complexity')
    assert evoca_paths.runs_dir() == os.path.join(evoca_paths.RESEARCH, 'Ring-complexity', 'Runs')
    monkeypatch.delenv('EVOCA_WORKSTREAM')
    assert evoca_paths.runs_dir().endswith(os.path.join('_unsorted', 'Runs'))
    s.free()
