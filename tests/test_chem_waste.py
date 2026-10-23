"""Chem milestone 2: food/waste as a special case of the general chemistry.

Waste (compound 1) carries no energy, is what the tax converts food into, is
excreted into its own environmental field W at d = waste_disp * w * (1 - W),
kills a cell when its private amount reaches 1, and in the environment has no
source, a uniform sink waste_sink and diffusion wdiff (Docs/plan-chem.md §7).

Guards:
  1. The general N_chem > 1 code paths with every compound inert reproduce
     the single-food model bit for bit (same fingerprint as N_chem = 1).
  2. Constructed one-step fixtures check metabolism, excretion, the sink,
     toxic death (compounds returned to the site), and the 50/50 split.
  3. Conservation on evolving runs: the waste mass identity
     d(sum W + sum w) = in - out every step (covers eviction, tax death,
     toxic death, reproduction, diffusion), and the energy identity.
  4. Public doors: init, run_sim (_PARAM_KEYS), recipes, gallery overrides.
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

# A setting where the population persists with heavy toxic-death turnover (alive ~0.6 at
# N=48); at waste_sink >= tax the sink removes waste as fast as it is made.
WASTE = dict(n_chem=2, wdiff=0.3, waste_sink=0.01, waste_disp=0.2)


def _fingerprint(s):
    h = hashlib.sha1()
    for a in (s.get_F(), s.get_f(), s.get_v(), s.get_alive(), s.get_lut_hashes()):
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


def _evolving(make_sim, N=48, seed=1, n_chem_after_init=None, **kw):
    rng = np.random.default_rng(seed)
    p = dict(food_inc=0.015, m_scale=3.0, gdiff=0.08, tax=0.02, mu_lut=0.01, mu_egene=0.03)
    p.update(kw)
    s = make_sim(N, **p)
    if n_chem_after_init is not None:
        s._lib.evoca_set_n_chem(n_chem_after_init)  # C level: every compound c >= 1 inert
    s.set_seed(seed)
    s.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    s.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.5); s.set_F_all(1.0); s.set_alive_halfplane(0)
    return s


# ── 1. inert general path == single food ──────────────────────────────

def test_inert_n_chem_2_is_bit_for_bit_single_food(make_sim):
    """N_chem=2 with no metabolism, exchange or environment dynamics runs the
    general multi-compound paths (split, excretion loop, eviction return) but
    must give exactly the single-food run."""
    fp = {}
    for n in (None, 2):
        s = _evolving(make_sim, n_chem_after_init=n)
        assert s.get_n_chem() == (n or 1) and s.get_metab_product() == -1
        for _ in range(300):
            s.step()
        fp[n] = _fingerprint(s)
    assert fp[None] == fp[2]


def test_food_waste_changes_the_run(make_sim):
    fp = []
    for kw in ({}, WASTE):                      # one C lattice at a time: run sequentially
        s = _evolving(make_sim, **kw)
        for _ in range(200):
            s.step()
        fp.append(_fingerprint(s))
        if kw:
            assert s.get_W().sum() > 0 and s.get_w().max() > 0
    assert fp[0] != fp[1]


# ── 2. constructed one-step fixtures ──────────────────────────────────

def _still(make_sim, N=16, **kw):
    """All cells alive, no regrowth, diffusion, eating or mutation; f = 0.5."""
    p = dict(food_inc=0.0, m_scale=0.0, gdiff=0.0, n_chem=2)
    p.update(kw)
    s = make_sim(N, **p)
    s.state(lut='gol', egenome='uniform', egenome_value=3)
    s.set_F_all(0.0); s.set_alive_all(); s.set_f_all(0.5)
    return s


def test_metabolism_and_excretion_one_step(make_sim):
    t, disp = 0.03, 0.25
    s = _still(make_sim, tax=t, waste_disp=disp)
    s.step()
    w, W, f = s.get_w(), s.get_W(), s.get_f()
    assert np.allclose(f, 0.5 - t)                         # food -> waste ...
    assert np.allclose(W, t * disp, atol=1e-7)              # ... excreted into W = 0
    assert np.allclose(w, t * (1 - disp), atol=1e-7)
    assert np.allclose(w + W, t, atol=1e-7)                 # nothing lost
    s.step()                                                # second step: (1 - W) slows excretion
    w1 = t * (1 - disp) + t
    assert np.allclose(s.get_W(), t * disp + disp * w1 * (1 - t * disp), atol=1e-6)


def test_excretion_stops_when_environment_is_saturated(make_sim):
    s = _still(make_sim, tax=0.0, waste_disp=0.5)
    s.set_w(0.4)
    W = np.ones((16, 16), np.float32)
    ptr = s._lib.evoca_get_chem_env(1)
    np.ctypeslib.as_array(ptr, shape=(256,))[:] = W.ravel()
    s.step()
    assert np.allclose(s.get_w(), 0.4) and np.allclose(s.get_W(), 1.0)


def test_waste_sink_removes_uniformly_floored_at_zero(make_sim):
    s = _still(make_sim, tax=0.0, waste_sink=0.1)
    env = np.ctypeslib.as_array(s._lib.evoca_get_chem_env(1), shape=(256,))
    env[:] = np.linspace(0, 0.3, 256, dtype=np.float32)
    before = env.copy()
    s.step()
    assert np.allclose(s.get_W().ravel(), np.maximum(before - 0.1, 0), atol=1e-7)
    assert s.get_chem_flux(1)[1] == pytest.approx(np.minimum(before, 0.1).sum(), rel=1e-5)


def test_toxic_death_returns_food_and_waste_to_the_site(make_sim):
    s = _still(make_sim, tax=0.02)
    w0 = np.full((16, 16), 0.5, np.float32); w0[3, 4] = 0.99    # one cell crosses 1 this step
    s.set_w(w0)
    s.step()
    alive = s.get_alive().reshape(16, 16)
    assert alive[3, 4] == 0 and alive.sum() == 255
    assert s.get_tox_deaths_last() == 1
    assert s.get_F()[3, 4] == pytest.approx(0.48)               # f after tax, returned (not lost)
    assert s.get_W()[3, 4] == pytest.approx(1.01)               # its waste, returned
    assert s.get_w()[3, 4] == 0 and s.get_f()[3, 4] == 0


def test_reproduction_splits_waste_50_50(make_sim):
    s = _still(make_sim, N=8, tax=0.0)
    a = np.zeros((8, 8), np.uint8); a[4, 4] = 1
    s.set_alive(a); s.set_f_all(1.0); s.set_w(0.4)
    s.step()
    alive = s.get_alive().reshape(8, 8)
    assert alive.sum() == 2
    assert np.allclose(s.get_w()[alive == 1], 0.2) and np.allclose(s.get_f()[alive == 1], 0.5)


# ── 3. conservation on evolving runs ──────────────────────────────────

@pytest.mark.parametrize('extra', [dict(), dict(wdiff=0.0), dict(waste_disp=0.0),
                                   dict(waste_exchange=1), dict(waste_exchange=1, wdiff=0.0)])
def test_waste_mass_identity_every_step(make_sim, extra):
    """sum W + sum w changes by exactly (in - out) each step: nothing is
    created or lost by eviction, tax death, toxic death, the split or diffusion."""
    s = _evolving(make_sim, **{**WASTE, **extra})
    prev = s.get_W().astype(np.float64).sum() + s.get_w().astype(np.float64).sum()
    pin, pout = s.get_chem_flux(1)
    tox = 0
    for _ in range(300):
        s.step()
        tox += s.get_tox_deaths_last()
        cur = s.get_W().astype(np.float64).sum() + s.get_w().astype(np.float64).sum()
        cin, cout = s.get_chem_flux(1)
        assert cur - prev == pytest.approx((cin - pin) - (cout - pout), abs=2e-3)
        prev, pin, pout = cur, cin, cout
    if extra.get('waste_disp') == 0.0:
        assert tox > 0, "no excretion must poison cells (toxic deaths exercised)"
    else:
        assert s.get_alive().mean() > 0.02, "fixture went extinct; identity is vacuous"


@pytest.mark.parametrize('extra', [dict(), dict(waste_exchange=1)])
def test_energy_identity_holds_with_waste(make_sim, extra):
    """Waste carries no energy, and a toxic death returns food to F (a transfer,
    not a sink), so the single-food dissipation identity is unchanged."""
    s = _evolving(make_sim, **{**WASTE, **extra})
    d0 = s.get_dissipation(); prev = d0['storage_env'] + d0['storage_priv']
    tox = 0
    for _ in range(300):
        s.step(); tox += s.get_tox_deaths_last()
        d = s.get_dissipation(); cur = d['storage_env'] + d['storage_priv']
        assert d['source_step'] == pytest.approx(cur - prev + d['diss_step'], abs=2e-3)
        prev = cur
    assert tox > 0


# ── 4. public doors ───────────────────────────────────────────────────

def test_init_applies_the_food_waste_preset_and_refuses_misuse(make_sim):
    s = make_sim(16, **WASTE)
    assert s.get_n_chem() == 2 and s.get_metab_product() == 1
    assert s.get_chem_energy(1) == 0.0
    inc, sink, diff, disp, tox = s.get_chem_params(1)
    assert (inc, tox) == (0.0, 1.0)
    assert (sink, diff, disp) == pytest.approx((0.01, 0.3, 0.2))
    s.update_waste_disp(0.5)
    assert s.get_chem_params(1)[3] == pytest.approx(0.5)
    with pytest.raises(ValueError, match='need n_chem=2'):
        make_sim(16, wdiff=0.1)
    with pytest.raises(ValueError, match='n_chem=3'):
        make_sim(16, n_chem=3)
    t = make_sim(16)
    with pytest.raises(ValueError, match='n_chem=2'):
        t.update_wdiff(0.1)
    assert (t.get_n_chem(), t.get_metab_product()) == (1, -1)


def test_waste_reaches_sim_via_run_sim():
    from evoca_explore import run_sim
    P = dict(food_inc=0.015, m_scale=3.0, gdiff=0.08, tax=0.02, mu_lut=0.01)
    with warnings.catch_warnings():
        warnings.simplefilter('error')                      # a dropped key warns
        a = run_sim(P, n_steps=300, sample_every=50, N=32, seed=1, shadow=False)
        b = run_sim(dict(P, **WASTE), n_steps=300, sample_every=50, N=32, seed=1, shadow=False)
    assert (b['param_n_chem'], b['param_wdiff'], b['param_waste_sink'], b['param_waste_disp']) == \
        (2, 0.3, 0.01, 0.2)
    assert a['alive_density_mean'] != b['alive_density_mean']


@pytest.mark.parametrize('recipe_kind', ['init', 'final'])
def test_recipe_roundtrip(recipe_kind, tmp_path):
    from evoca_py import EvoCA, import_run
    s = EvoCA(); s.init(16, m_scale=3.0, **WASTE)
    s.state(lut='gol', egenome='uniform')
    path = s.export_recipe('_pytest_waste', runs_dir=str(tmp_path)); s.free()
    assert json.load(open(path))[f'metaparams_{recipe_kind}']['waste_disp'] == pytest.approx(0.2)
    sim, _ = import_run(path, recipe=recipe_kind)
    try:
        assert sim.get_n_chem() == 2 and sim.get_chem_params(1)[1:4] == pytest.approx((0.01, 0.3, 0.2))
    finally:
        sim.free()


def test_gallery_entry_override(tmp_path, monkeypatch):
    from evoca_py import EvoCA
    from evoca_gallery import Entry, Gallery
    s = EvoCA(); s.init(16, m_scale=3.0, **WASTE); s.state(lut='gol', egenome='uniform')
    path = s.export_recipe('_pytest_waste_gal', runs_dir=str(tmp_path)); s.free()
    e = Entry(recipe=path, overrides={'waste_sink': 0.03})
    sim, _ = e.load(N=16)
    try:
        assert sim.get_chem_params(1)[1] == pytest.approx(0.03)
    finally:
        sim.free()
    t = Gallery('t', [e]).table()
    assert t.n_chem.iloc[0] == 2 and t.waste_sink.iloc[0] == pytest.approx(0.03)


def _colormode_fixture(make_sim):
    """Row 0 dead, the rest alive; env F = 0.5, env W a ramp 0..1.5 (to
    test the clip), private f = 0.5, private w a ramp 0..0.6."""
    s = _still(make_sim, N=8, tax=0.0)
    a = np.ones((8, 8), np.uint8); a[0, :] = 0
    s.set_alive(a)
    s.set_F_all(0.5)
    env = np.ctypeslib.as_array(s._lib.evoca_get_chem_env(1), shape=(64,))
    env[:] = np.linspace(0, 1.5, 64, dtype=np.float32)
    s.set_w(np.linspace(0, 0.6, 64, dtype=np.float32).reshape(8, 8))
    w = s.get_w().ravel()
    f = s.get_f().ravel()
    return s, a.ravel() == 1, env.copy(), w, f


def _channels(s, cm):
    px = s.colorize(np.zeros(64, np.int32), cm).view(np.uint32)
    return (px >> 16) & 255, (px >> 8) & 255, px & 255


def _u8(x):
    return (np.clip(x, 0, 1) * 255).astype(np.uint32)


def test_food_waste_colormode_channels(make_sim):
    """Colormode 6, built like env-food (1): dead = env F green + env W blue;
    alive = red 180 + private f green + private w blue."""
    s, alive, env, w, f = _colormode_fixture(make_sim)
    assert w[alive].max() > 0 and np.all(w[~alive] == 0)
    r, g, b = _channels(s, 6)
    assert np.array_equal(r, np.where(alive, 180, 0))
    assert np.array_equal(g, np.where(alive, _u8(f), _u8(np.full(64, 0.5))))
    assert np.array_equal(b, np.where(alive, _u8(w), _u8(env)))


def test_env_waste_colormode_channels(make_sim):
    """Colormode 7, the waste analog of env-food: dead = env W blue;
    alive = red 180 + private w blue; no green."""
    s, alive, env, w, f = _colormode_fixture(make_sim)
    r, g, b = _channels(s, 7)
    assert np.array_equal(r, np.where(alive, 180, 0))
    assert np.all(g == 0)
    assert np.array_equal(b, np.where(alive, _u8(w), _u8(env)))
