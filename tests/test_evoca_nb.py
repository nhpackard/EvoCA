"""evoca_nb: the lab-notebook helpers reproduce the private copies they
replaced, bit for bit, entering through the public functions.

Each OLD_* below is the notebook code as it stood on main 5437015 (copied
verbatim, modulo names). Only one EvoCA lattice exists per process (the C
state is global), so each comparison builds, fingerprints and frees one
sim before building the other.
"""
import hashlib
import json
import os
import tempfile
import warnings

import numpy as np
import pandas as pd
import pytest

import evoca_nb as nb
import evoca_explore as ee
from evoca_py import EvoCA, import_run

RECIPE = '2026-03-16_2_branch_activity.evoca'          # Activity-metrics/Runs
RD_RECIPE = '2026-04-25_waves_in_blank_half.evoca'      # Reaction-diffusion/Runs
R1_SCAN = '2026-05-18_R1_mscale_bracket'                # Resource-driving/Scans
N, STEPS = 32, 50


def fingerprint(sim, steps=STEPS):
    """Step, then hash every per-cell field plus the step counter (and the
    lineage ids / shadow table when present); frees the sim."""
    for _ in range(steps):
        sim.step()
    h = hashlib.sha256()
    for a in (sim.get_v(), sim.get_alive(), sim.get_f(), sim.get_F(),
              sim.get_egenome(), sim.get_egenes(), sim.get_active(),
              sim.get_lut_hashes(), np.array([sim.get_step()])):
        h.update(np.ascontiguousarray(a).tobytes())
    if sim.lineage_enabled():
        h.update(sim.get_lineage_birth_id().tobytes())
    S = sim.get_n_activity(max_n=100000)
    h.update(np.sort(S['hash']).tobytes())
    sim.free()
    return h.hexdigest()


def pin():
    """Pin the generators import_run consumes *before* the helpers seed
    (the global numpy RNG for an unseeded random LUT, and the C PRNG), so
    old and new start from the same pre-seed state."""
    np.random.seed(123)
    EvoCA().set_seed(123)


# ── old notebook code (verbatim) ──────────────────────────────────────

def OLD_activity_load_recipe(name, N=96, seed=0):
    sim, kw = import_run(name, N=N)
    sim.set_seed(seed)
    sim.set_v(np.random.default_rng(seed).integers(0, 2, (N, N), dtype=np.uint8))
    return sim, kw


def OLD_activity_argb_to_rgb(a):
    a = a.astype(np.uint32)
    return np.stack([(a >> 16) & 255, (a >> 8) & 255, a & 255], -1).astype(np.uint8)


def OLD_rd_load(name, N, recipe='init', seed=0):
    sim, kw = import_run(name, N=N, recipe=recipe); sim.set_seed(seed)
    return sim


def OLD_rd_rgb(sim, N, cm):
    buf = np.zeros(N*N, np.int32); sim.colorize(buf, cm)
    a = buf.view(np.uint32).reshape(N, N)
    return np.stack([(a>>16)&255, (a>>8)&255, a&255], -1).astype(np.uint8)


def OLD_resource_rgb(sim, colormode):
    N = sim._N
    a = sim.colorize(np.zeros(N * N, np.int32), colormode).reshape(N, N).astype(np.int64)
    return np.stack([(a >> 16) & 255, (a >> 8) & 255, a & 255], axis=-1).astype(np.uint8)


def OLD_r1_sim(config_idx, N=128, seed=0, food_inc=None, shadow=False, runs=None):
    path = ee.evoca_from_scan(R1_SCAN, config_idx, runs_dir=runs)
    sim, kw = import_run(path, N=N)
    sim.set_seed(seed)
    sim.set_v(np.random.default_rng(seed).integers(0, 2, (N, N), dtype=np.uint8))
    if food_inc is not None:
        sim.update_food_inc(food_inc)
    if shadow:
        sim.neutral_enable()
    return sim, path


def OLD_params_from_row(r):
    p = {k: r['param_' + k] for k in ee._PARAM_KEYS if 'param_' + k in r.index}
    p['restricted_mu'] = int(p.get('restricted_mu', 0))
    return {k: float(v) if k != 'restricted_mu' else v for k, v in p.items()}


REG3C = dict(m_scale=3.0, food_inc=0.015, tax=0.035, gdiff=0.06)


def OLD_new_3c(N, seed, mu_lut, lineage=True):
    s = EvoCA()
    s.init(N, mu_lut=mu_lut, mu_egene=0.0, mu_egenome=0.0, tax_per_egene=0.0, lineage=lineage, **REG3C)
    s.set_seed(seed)
    s.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    s.set_v(np.random.default_rng(seed).integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.5); s.set_F_all(1.0); s.set_alive_halfplane(0)
    s.neutral_enable()
    return s


def OLD_edge_table_setup(mu_lut, N=64, seed=3):
    s = EvoCA()
    s.init(N, mu_lut=mu_lut, mu_egene=0.01, mu_egenome=0.01, lineage=True, **REG3C)
    s.set_seed(seed); s.state(lut='gol', egenome='uniform', egenome_value=3)
    s.set_v(np.random.default_rng(seed).integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.5); s.set_F_all(1.0); s.set_alive_halfplane(0)
    return s


def OLD_run_lineage_setup(N, seed, mu_lut, mu_eg, mu_egn, egenome):
    s = EvoCA()
    s.init(N, mu_lut=mu_lut, mu_egene=mu_eg, mu_egenome=mu_egn, tax_per_egene=0.0, lineage=True, **REG3C)
    s.set_seed(seed); s.state(lut='gol', egenome=egenome, egenome_value=3)
    s.set_v(np.random.default_rng(seed).integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.5); s.set_F_all(1.0); s.set_alive_halfplane(0)
    return s


def OLD_ring_make_sim(N, seed, **kw):
    s = EvoCA()
    s.init(N, **kw)
    s.set_seed(seed)
    return s


# ── the fingerprint discriminates ─────────────────────────────────────

def test_fingerprint_discriminates_seeds():
    pin(); a = fingerprint(nb.load_recipe(RECIPE, N=N, seed=0)[0])
    pin(); b = fingerprint(nb.load_recipe(RECIPE, N=N, seed=1)[0])
    assert a != b


# ── equivalence: recipes and scans ────────────────────────────────────

@pytest.mark.parametrize('seed', [0, 3])
def test_load_recipe_matches_activity_helper(seed):
    pin(); old = fingerprint(OLD_activity_load_recipe(RECIPE, N=N, seed=seed)[0])
    pin(); new = fingerprint(nb.load_recipe(RECIPE, N=N, seed=seed)[0])
    assert old == new


@pytest.mark.parametrize('recipe', ['init', 'final'])
def test_load_recipe_unseeded_v_matches_rd_load(recipe):
    # RD's load() left v from import_run's unseeded draw: compare the C side
    # with the same v forced afterwards, and the metaparams.
    v = np.random.default_rng(9).integers(0, 2, (N, N), dtype=np.uint8)
    pin(); s = OLD_rd_load(RD_RECIPE, N, recipe=recipe, seed=2); s.set_v(v)
    mp_old = dict(s._init_metaparams, food_inc=s.food_inc, gdiff=s.gdiff)
    old = fingerprint(s)
    pin(); s, _ = nb.load_recipe(RD_RECIPE, N=N, recipe=recipe, seed=2, seed_v=False); s.set_v(v)
    mp_new = dict(s._init_metaparams, food_inc=s.food_inc, gdiff=s.gdiff)
    new = fingerprint(s)
    assert old == new and mp_old == mp_new


@pytest.mark.parametrize('kw', [dict(), dict(food_inc=0.016), dict(food_inc=0.05, shadow=True)])
def test_sim_from_scan_matches_r1_sim(kw, tmp_path):
    pin(); old = fingerprint(OLD_r1_sim(208, N=N, seed=0, runs=str(tmp_path), **kw)[0])
    pin(); new = fingerprint(nb.sim_from_scan(R1_SCAN, 208, N=N, seed=0, **kw)[0])
    assert old == new


def test_sim_from_scan_leaves_no_files(tmp_path):
    before = set(os.listdir(tempfile.gettempdir()))
    sim, _ = nb.sim_from_scan(R1_SCAN, 45, N=16)
    sim.free()
    assert not [d for d in set(os.listdir(tempfile.gettempdir())) - before
                if d.startswith('evoca_nb_')]


def test_override_without_setter_raises():
    with pytest.raises(TypeError, match='bogus'):
        nb.load_recipe(RECIPE, N=16, bogus=1)
    with pytest.raises(TypeError, match='bogus'):
        nb.sim_from_scan(R1_SCAN, 45, N=16, bogus=1)


def test_params_from_row_matches():
    d = pd.read_csv(os.path.join(ee.evoca_paths.scans_dir('Lineage-co-evolution'),
                                 '2026-05-17_pure_evo_egene_on_winners', 'results.csv')).set_index('config_idx')
    # The old helper plus the pre-A9 load-time correction (this scan has no
    # param_food_cap column; its rows load the standard capped rule).
    assert 'param_food_cap' not in d.columns
    for cfg in (236, 237, 252, 253):
        assert nb.params_from_row(d.loc[cfg]) == dict(OLD_params_from_row(d.loc[cfg]), food_cap=1)
    p = nb.params_from_row(d.loc[252])
    old = fingerprint(ee.build_sim(dict(OLD_params_from_row(d.loc[252]), food_cap=1), N=N, seed=0))
    new = fingerprint(ee.build_sim(p, N=N, seed=0))
    assert old == new


def test_params_from_row_agrees_with_evoca_from_scan(tmp_path):
    """Both scan-row loaders apply the same load-time corrections (one
    source: ee.correct_scan_row_params), so a row means the same run
    whichever path loads it."""
    import json
    scan = '2026-05-18_R1_mscale_bracket'          # pre-A9, has scan_config.json
    d = pd.read_csv(os.path.join(ee.evoca_paths.scans_dir('Resource-driving'),
                                 scan, 'results.csv')).set_index('config_idx')
    assert 'param_food_cap' not in d.columns
    rec = json.load(open(ee.evoca_from_scan(scan, 208, runs_dir=str(tmp_path))))
    mp = rec['metaparams_init']
    p = nb.params_from_row(d.loc[208])
    assert p['food_cap'] == mp['food_cap'] == 1             # standard rule (2026-10-23)
    for k, v in p.items():
        assert float(mp[k]) == float(v), k


# ── equivalence: build_sim replaces the Lineage hand setups ───────────

@pytest.mark.parametrize('lineage', [True, False])
def test_build_sim_matches_new_3c(lineage):
    old = fingerprint(OLD_new_3c(N, 1, 0.03, lineage=lineage))
    s = ee.build_sim(dict(REG3C, mu_lut=0.03, mu_egene=0.0, mu_egenome=0.0, tax_per_egene=0.0),
                     N=N, seed=1, lineage=lineage)
    s.neutral_enable()
    assert old == fingerprint(s)


def test_build_sim_matches_edge_table_setup():
    old = fingerprint(OLD_edge_table_setup(0.01, N=N, seed=3))
    new = fingerprint(ee.build_sim(dict(REG3C, mu_lut=0.01, mu_egene=0.01, mu_egenome=0.01),
                                   N=N, seed=3, lineage=True))
    assert old == new


@pytest.mark.parametrize('egenome', ['uniform', 'random'])
def test_build_sim_matches_run_lineage_setup(egenome):
    old = fingerprint(OLD_run_lineage_setup(N, 2, 0.0003, 0.003, 0.003, egenome))
    new = fingerprint(ee.build_sim(dict(REG3C, mu_lut=0.0003, mu_egene=0.003, mu_egenome=0.003,
                                        tax_per_egene=0.0), N=N, seed=2, lineage=True, egenome=egenome))
    assert old == new


def test_build_sim_defaults_unchanged():
    # the scan harness path: no lineage, uniform egenome
    p = dict(food_inc=0.013, m_scale=1.2, gdiff=0.06, tax=0.035, mu_lut=0.01)
    s = EvoCA(); s.init(N, **p); s.set_seed(4)
    s.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    s.set_v(np.random.default_rng(4).integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.5); s.set_F_all(1.0); s.set_alive_halfplane(0)
    old = fingerprint(s)
    s = ee.build_sim(p, N=N, seed=4)
    assert not s.lineage_enabled()
    assert old == fingerprint(s)


def test_make_sim_matches_ring_helper():
    kw = dict(food_inc=0.01, m_scale=2.5, mu_lut=0.002, tax_ring=0.002)
    st = dict(lut='gol', egenome='uniform', egenome_value=3, f_init=0.5, F_init=1.0, alive='halfplane')
    v = np.random.default_rng(5).integers(0, 2, (N, N), dtype=np.uint8)
    s = OLD_ring_make_sim(N, 1, **kw); s.state(**st); s.set_v(v)
    old = fingerprint(s)
    s = nb.make_sim(N, 1, **kw); s.state(**st); s.set_v(v)
    assert old == fingerprint(s)


# ── frames ────────────────────────────────────────────────────────────

def test_rgb_matches_all_old_variants():
    sim, _ = nb.load_recipe(RECIPE, N=N, seed=0)
    for _ in range(20):
        sim.step()
    for cm in range(6):
        # colorize() is stateful in colormode 4 (age), so every variant
        # converts the same frame: replay one buffer through colorize().
        frame = sim.colorize(np.zeros(N * N, np.int32), cm).copy()
        real = sim.colorize
        sim.colorize = lambda px, mode=0: (np.copyto(px, frame), px)[1]
        new = nb.rgb(sim, cm)
        assert new.shape == (N, N, 3) and new.dtype == np.uint8
        assert np.array_equal(new, OLD_rd_rgb(sim, N, cm))
        assert np.array_equal(new, OLD_resource_rgb(sim, cm))
        assert np.array_equal(nb.argb_to_rgb(frame.reshape(N, N)), new)
        assert np.array_equal(OLD_activity_argb_to_rgb(frame.reshape(N, N)), new)
        sim.colorize = real
    sim.free()


def test_slope():
    t = np.arange(10.0)
    assert nb.slope(t, 3 * t + 1) == pytest.approx(3.0)
    assert nb.slope([1.0], [2.0]) == 0.0


# ── export_recipe_once ────────────────────────────────────────────────

def _sim_for_export():
    s = nb.make_sim(24, 1, food_inc=0.01, m_scale=2.5, mu_lut=0.002, tax_ring=0.002)
    s.state(lut='gol', egenome='uniform', egenome_value=3, f_init=0.5, F_init=1.0, alive='halfplane')
    return s


def _naive_export(sim, descriptor, runs_dir, **kw):
    return sim.export_recipe(descriptor, runs_dir=runs_dir, **kw)


def _rerun_twice(export, runs_dir):
    """Simulate two Run Alls on different days: the first file is
    back-dated so a dated re-export would get a new name."""
    s = _sim_for_export()
    p1 = export(s, 'nb_once_test', runs_dir=runs_dir, probes={'ts': True})
    old = os.path.join(runs_dir, '2026-01-01_nb_once_test.evoca')
    os.rename(p1, old)
    p2 = export(s, 'nb_once_test', runs_dir=runs_dir, probes={'ts': True})
    s.free()
    return old, p2, sorted(os.listdir(runs_dir))


def test_export_recipe_once_is_idempotent(tmp_path):
    old, p2, files = _rerun_twice(nb.export_recipe_once, str(tmp_path))
    assert files == ['2026-01-01_nb_once_test.evoca'] and p2 == old


def test_naive_export_writes_a_second_copy(tmp_path):
    # the bug export_recipe_once fixes: the naive export litters Runs/
    old, p2, files = _rerun_twice(_naive_export, str(tmp_path))
    assert len(files) == 2 and p2 != old


def test_export_recipe_once_writes_when_absent_and_matches_export(tmp_path):
    a, b = tmp_path / 'a', tmp_path / 'b'
    s = _sim_for_export()
    pa = nb.export_recipe_once(s, 'x y', runs_dir=str(a), colormode=1)
    pb = s.export_recipe('x y', runs_dir=str(b), colormode=1)
    s.free()
    assert os.path.basename(pa) == os.path.basename(pb)
    ja, jb = json.load(open(pa)), json.load(open(pb))
    ja.pop('created'); jb.pop('created')
    assert ja == jb


def test_export_recipe_once_ignores_other_descriptors(tmp_path):
    (tmp_path / '2026-01-01_other_nb_once_test.evoca').write_text('{}')
    s = _sim_for_export()
    p = nb.export_recipe_once(s, 'nb_once_test', runs_dir=str(tmp_path))
    s.free()
    assert os.path.basename(p).endswith('_nb_once_test.evoca') and len(os.listdir(tmp_path)) == 2


def test_export_recipe_once_warns_when_existing_differs(tmp_path):
    s = _sim_for_export()
    p1 = nb.export_recipe_once(s, 'nb_once_test', runs_dir=str(tmp_path))
    d = json.load(open(p1))
    for sec in ('metaparams_init', 'metaparams_final'):
        del d[sec]['food_cap']                 # missing key == init default (1): no warning
    json.dump(d, open(p1, 'w'))
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        assert nb.export_recipe_once(s, 'nb_once_test', runs_dir=str(tmp_path)) == p1
    s.update_food_inc(0.02)
    with pytest.warns(UserWarning, match='food_inc'):
        assert nb.export_recipe_once(s, 'nb_once_test', runs_dir=str(tmp_path)) == p1
    s.free()
    assert len(os.listdir(tmp_path)) == 1
