"""Tests for the refractory eating gates `food_min` and `mouthful_min` (2026-10-22).

N's rule: M'(x) = 0 when F(x) < food_min, otherwise the usual capped mouthful
M(x) = min((m_scale/25)*max(0,score)*F, F, 1-f). The second gate, mouthful_min,
refuses any mouthful below the threshold, so the food level at which a site can
be eaten again depends on the eater's score. Both default to 0 (off).

Guards:
  1. Gates off (default) is bit-for-bit main b374612, capped and legacy rule:
     golden SHA1 fingerprints recorded on main.
  2. food_min: in one step where only eating can change F, no site with
     F < food_min loses food, and sites at or above it are eaten as before.
  3. mouthful_min: every realised mouthful is 0 or >= mouthful_min; the
     ungated step has mouthfuls in (0, mouthful_min) (fixture discriminates).
  4. F >= 0 and the energy identity still close with the gates on.
  5. Defaults, re-init reset, params(), recipe round trip, run_sim and the
     gallery door (Entry.load overrides; table columns).
  6. Under the legacy rule (food_cap=0) the gates are still well defined.
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


# ── 1. gates off == main ──────────────────────────────────────────────

# SHA1 of (F, f, v, alive, lut_hashes) after 200 steps of the N=32 fixture
# below, recorded on main b374612 (before the gates), macOS arm64. The
# food_cap=0 values equal test_food_cap._MAIN_GOLDEN (pre-A9 main).
_MAIN_GOLDEN = {
    (1, 3.0): 'b98b1e1822a343a3a5c5f8c260c8b0c762e77cdb',
    (1, 4.0): '7128314e205c0f9b9a2a74fd18f0c92f633f5453',
    (0, 3.0): '0598f4e59b3a145fb7c150d3c3cd2f5cc1be37ab',
    (0, 4.0): '81246ef5e4b205e9fcc873533c495d41b8af23d4',
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


@pytest.mark.parametrize("food_cap,m_scale", sorted(_MAIN_GOLDEN))
def test_gates_off_match_main_golden(make_sim, food_cap, m_scale):
    assert _golden_fp(make_sim, m_scale, food_cap=food_cap) == \
        _MAIN_GOLDEN[(food_cap, m_scale)]


@pytest.mark.parametrize("gate", ['food_min', 'mouthful_min'])
def test_gate_on_changes_the_run(make_sim, gate):
    """The fingerprint fixture is sensitive to each gate (else test 1 is vacuous)."""
    assert _golden_fp(make_sim, 3.0, **{gate: 0.3}) != _MAIN_GOLDEN[(1, 3.0)]


# ── 2-3. one constructed step: only eating changes F ──────────────────

def _eat_once(make_sim, **kw):
    """All cells alive, f=0.3 (headroom 0.7), F uniform random in [0,1];
    no regrowth, no diffusion, no tax, no mutation: after one step,
    F_before - F_after is exactly each site's realised mouthful."""
    N = 24
    s = make_sim(N, food_inc=0.0, m_scale=3.5, gdiff=0.0, tax=0.0, **kw)
    s.set_seed(2)
    s.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    s.set_v(np.random.default_rng(2).integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.3)
    s.set_alive_all()
    np.random.seed(3)
    s.set_F_random(0.0, 1.0)
    F0 = s.get_F().astype(np.float64).copy()
    s.step()
    return F0, F0 - s.get_F().astype(np.float64)


def test_food_min_blocks_sites_below_threshold(make_sim):
    th = 0.4
    F0, eaten = _eat_once(make_sim, food_min=th)
    _, eaten0 = _eat_once(make_sim)
    below = F0 < th
    assert (eaten0[below] > 1e-6).sum() > 20, "ungated step must eat below threshold"
    assert np.all(eaten[below] == 0.0)
    # at or above the threshold the gate does nothing
    assert np.array_equal(eaten[~below], eaten0[~below])


def test_mouthful_min_refuses_small_mouthfuls(make_sim):
    th = 0.25
    _, eaten = _eat_once(make_sim, mouthful_min=th)
    _, eaten0 = _eat_once(make_sim)
    small = (eaten0 > 1e-6) & (eaten0 < th - 1e-6)
    assert small.sum() > 20, "ungated step must take mouthfuls below threshold"
    took = eaten > 0
    assert np.all(eaten[took] >= th - 1e-6)
    assert np.all(eaten[small] == 0.0)
    assert np.array_equal(eaten[eaten0 >= th], eaten0[eaten0 >= th])


# ── 4. F >= 0 and energy identity with gates on ───────────────────────

def _evolving(make_sim, N=48, seed=1, **kw):
    rng = np.random.default_rng(seed)
    p = dict(food_inc=0.015, m_scale=3.5, gdiff=0.08, tax=0.01,
             mu_lut=0.01, mu_egenome=0.03)
    p.update(kw)
    s = make_sim(N, **p)
    s.set_seed(seed)
    s.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    s.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.5); s.set_F_all(1.0); s.set_alive_halfplane(0)
    return s


@pytest.mark.parametrize("gate", [dict(food_min=0.3), dict(mouthful_min=0.1)])
def test_gated_run_keeps_food_nonnegative_and_balances(make_sim, gate):
    s = _evolving(make_sim, **gate)
    d0 = s.get_dissipation()
    prev = d0['storage_env'] + d0['storage_priv']
    for _ in range(200):
        s.step()
        d = s.get_dissipation()
        cur = d['storage_env'] + d['storage_priv']
        assert d['source_step'] == pytest.approx(cur - prev + d['diss_step'], abs=1e-3)
        prev = cur
        assert s.get_F().min() >= 0.0
    assert s.get_alive().mean() > 0.05, "fixture went extinct; test is vacuous"


# ── 5. defaults, recipes, run_sim, gallery ───────────────────────────

def test_defaults_and_reinit_reset(make_sim):
    s = make_sim(16)
    assert (s.food_min, s.mouthful_min) == (0.0, 0.0)
    assert (s.get_food_min(), s.get_mouthful_min()) == (0.0, 0.0)
    s.update_food_min(0.3); s.update_mouthful_min(0.05)
    assert s.get_food_min() == pytest.approx(0.3)
    assert s.get_mouthful_min() == pytest.approx(0.05)
    assert s.params()['food_min'] == pytest.approx(0.3)
    s.update_food_min(-1)                                 # clamped at 0
    assert s.get_food_min() == 0.0
    t = make_sim(16, food_min=0.2)
    assert t.get_food_min() == pytest.approx(0.2)
    u = make_sim(16)                                      # no carry-over of the C global
    assert (u.get_food_min(), u.get_mouthful_min()) == (0.0, 0.0)


@pytest.mark.parametrize("recipe_kind", ["init", "final"])
def test_recipe_roundtrip(recipe_kind, tmp_path):
    from evoca_py import EvoCA, import_run
    s = EvoCA()
    s.init(16, m_scale=3.0, food_min=0.3, mouthful_min=0.05)
    s.state(lut='gol', egenome='uniform')
    path = s.export_recipe('_pytest_food_thresh', runs_dir=str(tmp_path))
    s.free()
    data = json.load(open(path))
    assert data[f'metaparams_{recipe_kind}']['food_min'] == pytest.approx(0.3)
    sim, _ = import_run(path, recipe=recipe_kind)
    try:
        assert sim.get_food_min() == pytest.approx(0.3)
        assert sim.get_mouthful_min() == pytest.approx(0.05)
    finally:
        sim.free()


def test_gates_reach_sim_via_run_sim():
    from evoca_explore import run_sim
    P = dict(food_inc=0.015, m_scale=3.5, gdiff=0.08, tax=0.01,
             mu_lut=0.01, mu_egenome=0.03)
    with warnings.catch_warnings():
        warnings.simplefilter('error')                    # a dropped key warns
        a = run_sim(P, n_steps=200, sample_every=50, N=32, seed=1, shadow=False)
        b = run_sim(dict(P, food_min=0.3), n_steps=200, sample_every=50, N=32, seed=1, shadow=False)
        c = run_sim(dict(P, mouthful_min=0.1), n_steps=200, sample_every=50, N=32, seed=1, shadow=False)
    assert (a['param_food_min'], a['param_mouthful_min']) == (0.0, 0.0)
    assert b['param_food_min'] == pytest.approx(0.3)
    assert c['param_mouthful_min'] == pytest.approx(0.1)
    assert b['mean_intake_mean'] != a['mean_intake_mean']
    assert c['mean_intake_mean'] != a['mean_intake_mean']


def test_gates_reach_sim_via_gallery():
    from evoca_gallery import Entry, Gallery
    rec = '2026-10-15_R1cfg208_intermediate_fi0p016.evoca'
    e = Entry(recipe=rec, params='final', overrides={'food_cap': 1, 'food_min': 0.3})
    sim, _ = e.load(N=16)
    try:
        assert sim.get_food_min() == pytest.approx(0.3) and sim.get_food_cap() == 1
    finally:
        sim.free()
    t = Gallery('t', [e]).table()
    assert t.food_min.iloc[0] == pytest.approx(0.3) and t.mouthful_min.iloc[0] == 0.0


# ── 6. legacy rule with a gate ────────────────────────────────────────

def test_food_min_under_legacy_rule_blocks_debt_sites(make_sim):
    """With food_cap=0 the gate still applies: no site with F < food_min
    (including negative F) is eaten, so no negative mouthfuls occur there."""
    s = _evolving(make_sim, food_cap=0, m_scale=4.0, food_min=0.2, gdiff=0.0, food_inc=0.0)
    for _ in range(100):
        F0 = s.get_F().astype(np.float64).copy()
        s.step()
        eaten = F0 - s.get_F().astype(np.float64)
        assert np.all(eaten[F0 < 0.2] == 0.0)
