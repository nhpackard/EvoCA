"""Tests for `food_cap` — the Phase-3 mouthful cap (A9 fix, 2026-10-15).

Defect A9 (Docs/Oct2026/06-seeding-findings.md): the mouthful
`(m_scale/25)*score*F` was capped only by private headroom, never by the
food present, so F_food went negative whenever m_scale*score > 25, and a
negative F then made later mouthfuls negative (live cells with negative
f_priv). The fix caps the mouthful at max(F, 0) and max(headroom, 0).
`food_cap=1` (default) is the corrected rule; `food_cap=0` is the legacy
rule, kept so pre-fix runs can be reproduced bit-for-bit.

Guards:
  1. food_cap=1 at m_scale 3 and 4: F_food >= 0 and alive f_priv >= 0 at
     every step of a viable run (population alive at the end). This is
     the regression test for A9; it fails on the pre-fix library.
  2. The energy identity source == d(storage) + diss still closes, capped
     and legacy.
  3. food_cap=0 reproduces main (a0798d1) bit-for-bit: golden SHA1
     fingerprints recorded on main.
  4. At m_scale <= 1 the cap never binds: food_cap=1 == food_cap=0.
  5. Default is capped; params()/recipe round-trip; legacy-recipe notice.
  6. Reachability through run_sim (the scan door, via _PARAM_KEYS).
"""
import hashlib
import json
import os
import sys
import warnings

import numpy as np
import pytest

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))


def _sim(make_sim, *, m_scale, N=48, seed=1, **kw):
    """Viable evolving fixture: GoL LUT, plus-shape egenome, halfplane
    start, mutation on. At m_scale 3-4 the legacy rule drives F negative
    within a few hundred steps (A9)."""
    rng = np.random.default_rng(seed)
    p = dict(food_inc=0.015, m_scale=m_scale, gdiff=0.08, tax=0.01,
             mu_lut=0.01, mu_egenome=0.03)
    p.update(kw)
    s = make_sim(N, **p)
    s.set_seed(seed)
    s.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    s.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.5)
    s.set_F_all(1.0)
    s.set_alive_halfplane(0)
    return s


# ── 1. A9 regression: no negative food, capped ────────────────────────

@pytest.mark.parametrize("m_scale", [3.0, 4.0])
def test_capped_food_never_negative(make_sim, m_scale):
    s = _sim(make_sim, m_scale=m_scale)
    assert s.get_food_cap() == 1
    for t in range(300):
        s.step()
        F = s.get_F()
        alive = s.get_alive().astype(bool)
        f = s.get_f()
        assert F.min() >= 0.0, f"F_food < 0 at step {t+1}: {F.min()}"
        if alive.any():
            assert f[alive].min() >= 0.0, \
                f"alive f_priv < 0 at step {t+1}: {f[alive].min()}"
    assert s.get_alive().mean() > 0.1, "fixture went extinct; test is vacuous"


@pytest.mark.parametrize("m_scale", [3.0, 4.0])
def test_legacy_rule_still_overdraws(make_sim, m_scale):
    """food_cap=0 must keep the defect (otherwise it is not the legacy
    rule). Also shows the fixture discriminates the two rules."""
    s = _sim(make_sim, m_scale=m_scale, food_cap=0)
    minF = min(float((s.step(), s.get_F().min())[1]) for _ in range(300))
    assert minF < -0.01


# ── 2. energy identity closes ─────────────────────────────────────────

@pytest.mark.parametrize("food_cap", [1, 0])
def test_energy_identity_closes(make_sim, food_cap):
    s = _sim(make_sim, m_scale=4.0, food_cap=food_cap)
    d0 = s.get_dissipation()
    prev = d0['storage_env'] + d0['storage_priv']
    for _ in range(150):
        s.step()
        d = s.get_dissipation()
        cur = d['storage_env'] + d['storage_priv']
        assert d['source_step'] == pytest.approx(cur - prev + d['diss_step'],
                                                 abs=1e-3)
        prev = cur


def test_mouthful_accounting_matches_food_removed(make_sim):
    """With all sinks on the env side known, the realised mouthful the
    accountant reports equals what left F_food: env storage change ==
    source - overflow - mouthful (gdiff=0 so diffusion is exact)."""
    s = _sim(make_sim, m_scale=4.0, gdiff=0.0)
    prevF = float(s.get_F().astype(np.float64).sum())
    for _ in range(100):
        s.step()
        d = s.get_dissipation()
        curF = float(s.get_F().astype(np.float64).sum())
        assert curF - prevF == pytest.approx(
            d['source_step'] - d['overflow_step'] - d['mouthful_step'],
            abs=1e-3)
        prevF = curF


# ── 3. food_cap=0 is bit-for-bit main ─────────────────────────────────

# SHA1 of (F, f, v, alive, lut_hashes) after 200 steps of the N=32
# fixture below, recorded on main a0798d1 (pre-fix library), macOS arm64.
_MAIN_GOLDEN = {
    3.0: '0598f4e59b3a145fb7c150d3c3cd2f5cc1be37ab',
    4.0: '81246ef5e4b205e9fcc873533c495d41b8af23d4',
}


def _golden_fp(make_sim, m_scale, **kw):
    N = 32
    rng = np.random.default_rng(5)
    s = make_sim(N, food_inc=0.015, m_scale=m_scale, gdiff=0.08, tax=0.01,
                 mu_lut=0.01, mu_egenome=0.03, **kw)
    s.set_seed(5)
    s.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    s.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.5); s.set_F_all(1.0); s.set_alive_halfplane(0)
    for _ in range(200):
        s.step()
    h = hashlib.sha1()
    for a in (s.get_F(), s.get_f(), s.get_v(), s.get_alive(),
              s.get_lut_hashes()):
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


@pytest.mark.parametrize("m_scale", [3.0, 4.0])
def test_food_cap0_matches_main_golden(make_sim, m_scale):
    assert _golden_fp(make_sim, m_scale, food_cap=0) == _MAIN_GOLDEN[m_scale]


@pytest.mark.parametrize("m_scale", [3.0, 4.0])
def test_default_differs_from_main_above_one(make_sim, m_scale):
    """The fix is a real model change for m_scale > 1."""
    assert _golden_fp(make_sim, m_scale) != _MAIN_GOLDEN[m_scale]


# ── 4. cap never binds at m_scale <= 1 ────────────────────────────────

@pytest.mark.parametrize("m_scale", [0.6, 1.0])
def test_cap_is_noop_at_m_scale_le_one(make_sim, m_scale):
    a = _sim(make_sim, m_scale=m_scale, food_cap=1, N=32)
    b = _sim(make_sim, m_scale=m_scale, food_cap=0, N=32)
    for _ in range(300):
        a.step(); b.step()
    for get in ('get_F', 'get_f', 'get_v', 'get_alive', 'get_lut_hashes'):
        assert np.array_equal(getattr(a, get)(), getattr(b, get)()), get
    assert a.get_alive().mean() > 0.1


# ── 5. defaults, params, recipes ──────────────────────────────────────

def test_default_is_capped(make_sim):
    s = make_sim(16)
    assert s.food_cap == 1 and s.get_food_cap() == 1
    assert s.params()['food_cap'] == 1
    s.update_food_cap(0)
    assert s.get_food_cap() == 0 and s.params()['food_cap'] == 0
    t = make_sim(16, food_cap=0)
    assert t.get_food_cap() == 0
    # re-init without the kwarg restores the default (no carry-over
    # of the C global between runs)
    u = make_sim(16)
    assert u.get_food_cap() == 1


@pytest.mark.parametrize("recipe_kind", ["init", "final"])
def test_food_cap_recipe_roundtrip(recipe_kind):
    from evoca_py import EvoCA, import_run
    s = EvoCA()
    s.init(16, m_scale=3.0, food_cap=0)
    s.state(lut='gol', egenome='uniform')
    path = s.export_recipe(f"_pytest_food_cap_{recipe_kind}")
    s.free()
    try:
        with open(path) as fh:
            data = json.load(fh)
        assert data[f'metaparams_{recipe_kind}']['food_cap'] == 0
        sim, _kw = import_run(path, recipe=recipe_kind)
        try:
            assert sim.food_cap == 0 and sim.get_food_cap() == 0
        finally:
            sim.free()
    finally:
        if os.path.exists(path):
            os.remove(path)


@pytest.mark.parametrize("m_scale,expect_notice", [(2.5, True), (1.0, False)])
def test_legacy_recipe_defaults_capped_with_notice(tmp_path, capsys,
                                                   m_scale, expect_notice):
    """A pre-fix recipe has no food_cap key: import capped, and print a
    one-line notice iff m_scale > 1 (where the rules differ)."""
    from evoca_py import import_run
    mp = dict(food_inc=0.015, m_scale=m_scale, tax=0.01)
    recipe = {'version': 3, 'descriptor': 'legacy', 'N': 16,
              'metaparams_init': mp, 'metaparams_final': mp,
              'initialization': {'lut': 'gol', 'egenome': 'uniform'},
              'display': {}}
    path = tmp_path / 'legacy.evoca'
    path.write_text(json.dumps(recipe))
    sim, _ = import_run(str(path))
    try:
        assert sim.get_food_cap() == 1
    finally:
        sim.free()
    out = capsys.readouterr().out
    if expect_notice:
        assert 'food_cap=0' in out and len(out.strip().splitlines()) == 1
    else:
        assert 'food_cap' not in out


# ── 6. reachability through run_sim ───────────────────────────────────

def test_food_cap_reaches_sim_via_run_sim():
    """run_sim filters params through evoca_explore._PARAM_KEYS. A
    dropped food_cap key would warn (an error here) and run capped in both
    calls, so the two runs would be identical."""
    from evoca_explore import run_sim
    P = dict(food_inc=0.015, m_scale=4.0, gdiff=0.08, tax=0.01,
             mu_lut=0.01, mu_egenome=0.03)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        a = run_sim(dict(P, food_cap=1), n_steps=200, sample_every=50,
                    N=32, seed=1, shadow=False)
        b = run_sim(dict(P, food_cap=0), n_steps=200, sample_every=50,
                    N=32, seed=1, shadow=False)
    assert a['param_food_cap'] == 1 and b['param_food_cap'] == 0
    assert a['alive_density_final'] > 0.1
    # The legacy rule moves different amounts of food (overdraw, then
    # negative mouthfuls where F < 0), so intake must differ.
    assert b['mean_intake_mean'] != a['mean_intake_mean']


def test_run_sim_always_records_effective_food_cap():
    from evoca_explore import run_sim
    r = run_sim(dict(food_inc=0.015, m_scale=2.0, tax=0.01), n_steps=10,
                sample_every=10, N=16, seed=1, shadow=False)
    assert r['param_food_cap'] == 1


def test_evoca_from_scan_pre_fix_rows_load_standard_rule(tmp_path, capsys):
    """A scan row with no param_food_cap ran uncapped, but the capped rule
    is the standard model (2026-10-23): its recipe carries food_cap=1
    explicitly, with a notice when the cap can bind (m_scale > 1). Rows
    that recorded the column keep their value, including the legacy 0."""
    import csv, json
    from evoca_explore import evoca_from_scan
    scan = tmp_path / '2026-01-01_old'
    scan.mkdir()
    import shutil, evoca_paths
    real = next(d for d in evoca_paths.all_scan_dirs()
                if os.path.exists(os.path.join(d, "scan_config.json")))
    shutil.copy(os.path.join(real, 'scan_config.json'), scan / 'scan_config.json')
    base = {'config_idx': 0, 'param_food_inc': '0.01', 'param_m_scale': '3.0',
            'param_tax': '0.03'}
    rows = [base, dict(base, config_idx=1, param_food_cap='1'),
            dict(base, config_idx=2, param_food_cap='0')]
    with open(scan / 'results.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[1]))
        w.writeheader()
        w.writerow(dict(rows[0], param_food_cap=''))
        w.writerows(rows[1:])
    out = tmp_path / 'runs'
    capsys.readouterr()
    r0 = json.load(open(evoca_from_scan(str(scan), 0, runs_dir=str(out))))
    assert 'measured under the legacy rule' in capsys.readouterr().out
    r1 = json.load(open(evoca_from_scan(str(scan), 1, runs_dir=str(out))))
    r2 = json.load(open(evoca_from_scan(str(scan), 2, runs_dir=str(out))))
    assert capsys.readouterr().out == ''               # recorded rows: no notice
    assert r0['metaparams_init']['food_cap'] == 1
    assert r1['metaparams_init']['food_cap'] == 1
    assert r2['metaparams_init']['food_cap'] == 0
