"""Tests for `food_per` — the temporal food-regeneration PERIOD
(checkpoint 3 July 2026 §C / action C1).

Phase 2 environmental regeneration fires only every `food_per`-th step;
when it fires it adds `food_inc * food_per`, so the long-run mean input
rate is held fixed while temporal lumpiness (and the lifetime of
self-generated depletion structure) rises with `food_per`.

Guards:
  1. food_per=1 reproduces the pre-change rule (`F += food_inc` every
     step) exactly, elementwise, in float32.
  2. The default init path (no food_per kwarg) is bit-for-bit identical
     to explicit food_per=1 under full dynamics (agents, eating, repro,
     mutation).
  3. food_per=K produces the intended sawtooth: F is flat between fires
     and jumps by food_inc*K on firing steps (g_step % K == 0).
  4. food_per=K is mean-preserving vs food_per=1 over a full period.
  5. food_per carries through params() and .evoca recipe round-trips.
  6. food_per < 1 clamps to 1 (Python + C guard).
"""
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))


# ── helpers ───────────────────────────────────────────────────────────

def _dead_regen_sim(make_sim, *, food_inc, food_per, N=16, F0=0.0):
    """A sim whose F field is driven ONLY by Phase-2 regeneration:
    all cells dead (no eating/tax/repro), no diffusion. This isolates
    the regen schedule so the arithmetic can be checked exactly."""
    s = make_sim(N, food_inc=food_inc, gdiff=0.0, food_per=food_per)
    s.state(lut='gol', egenome='uniform', egenome_value=0)
    s.set_alive_fraction(0.0)          # every cell dead
    s.set_f_all(0.0)
    s.set_F_all(F0)
    return s


def _full_sim(make_sim, *, food_per, N=32, seed=123):
    """A full evolving sim (eating, repro, mutation) parameterised only
    by food_per; used to check the food_per=1 no-op invariant."""
    rng = np.random.default_rng(seed)
    kw = dict(food_inc=0.02, m_scale=1.2, gdiff=0.0, tax=0.01,
              mu_lut=0.01, mu_egenome=0.01)
    if food_per is not None:
        kw['food_per'] = food_per
    s = make_sim(N, **kw)
    s.set_seed(seed)
    s.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    s.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.5)
    s.set_F_all(1.0)
    s.set_alive_halfplane(0)
    return s


def _fingerprint(s):
    return (s.get_F().tobytes(), s.get_v().tobytes(),
            s.get_alive().tobytes(), s.get_f().tobytes())


# ── 1. food_per=1 == pre-change rule (F += food_inc every step) ────────

def test_food_per1_matches_add_every_step(make_sim):
    food_inc, steps = 0.01, 50
    s = _dead_regen_sim(make_sim, food_inc=food_inc, food_per=1)

    # Reference: exactly `F += food_inc` every step, float32, clamped —
    # the literal pre-food_per Phase-2 body.
    ref = np.zeros((16, 16), dtype=np.float32)
    for _ in range(steps):
        ref = np.minimum(ref + np.float32(food_inc), np.float32(1.0))
        s.step()
        assert np.array_equal(s.get_F(), ref), "food_per=1 diverged from pre-change rule"


# ── 2. default path is bit-for-bit identical to explicit food_per=1 ────

def test_default_equals_explicit_food_per1(make_sim):
    a = _full_sim(make_sim, food_per=None)   # rely on C/Python default (1)
    b = _full_sim(make_sim, food_per=1)      # explicit
    for _ in range(120):
        a.step(); b.step()
    assert _fingerprint(a) == _fingerprint(b), \
        "default food_per path is not bit-for-bit == explicit food_per=1"


# ── 3. food_per=K sawtooth: flat between fires, jump on firing steps ───

def test_food_per_k_sawtooth(make_sim):
    food_inc, K = 0.01, 5
    s = _dead_regen_sim(make_sim, food_inc=food_inc, food_per=K)
    jump = np.float32(food_inc) * np.float32(K)   # 0.05, stays < 1

    # g_step increments at the top of each step(); regen fires when
    # g_step % K == 0, i.e. on the K-th, 2K-th, ... step() call.
    means = []
    for step_idx in range(1, 3 * K + 1):
        s.step()
        means.append(float(s.get_F().mean()))

    # Between fires the field is flat (no source, no sink).
    for i in range(1, K):
        assert means[i - 1] == pytest.approx(0.0), \
            f"F changed on non-firing step {i}"
    # Fires on steps K, 2K, 3K -> cumulative jump*n.
    assert means[K - 1] == pytest.approx(float(jump))
    assert means[2 * K - 1] == pytest.approx(float(2 * jump))
    assert means[3 * K - 1] == pytest.approx(float(3 * jump))


# ── 4. mean-preserving: after a full period, per-1 == per-K ────────────

def test_food_per_mean_preserving(make_sim):
    food_inc, K, periods = 0.005, 4, 3
    a = _dead_regen_sim(make_sim, food_inc=food_inc, food_per=1)
    b = _dead_regen_sim(make_sim, food_inc=food_inc, food_per=K)
    for _ in range(K * periods):
        a.step(); b.step()
    # No clamp in this regime (food_inc*K*periods = 0.06 < 1).
    assert a.get_F().mean() == pytest.approx(b.get_F().mean()), \
        "food_per=K not mean-preserving vs food_per=1"


# ── 5. round-trip through params() and .evoca recipe ───────────────────

def test_params_includes_food_per(make_sim):
    s = make_sim(16, food_per=7)
    assert s.params().get('food_per') == 7
    assert s._lib.evoca_get_food_per() == 7


@pytest.mark.parametrize("recipe_kind", ["init", "final"])
def test_food_per_recipe_roundtrip(recipe_kind):
    from evoca_py import EvoCA, import_run
    s = EvoCA()
    s.init(16, food_per=9)
    s.state(lut='gol', egenome='uniform')
    path = s.export_recipe(f"_pytest_food_per_{recipe_kind}")
    s.free()
    try:
        sim, _kw = import_run(path, recipe=recipe_kind)
        try:
            assert sim.food_per == 9, f"food_per lost in metaparams_{recipe_kind}"
            assert sim._lib.evoca_get_food_per() == 9
        finally:
            sim.free()
    finally:
        if os.path.exists(path):
            os.remove(path)


# ── 6. food_per < 1 clamps to 1 (bit-for-bit safety) ───────────────────

def test_food_per_clamps_to_one(make_sim):
    s = make_sim(16, food_inc=0.01, food_per=0)
    assert s.food_per == 1
    assert s._lib.evoca_get_food_per() == 1
    s.update_food_per(-3)
    assert s.food_per == 1
    assert s._lib.evoca_get_food_per() == 1


# ── 7. reachability through the public run_sim door ────────────────────

def test_food_per_reaches_sim_via_run_sim():
    """run_sim must forward food_per to the sim (it filters params through
    evoca_explore._PARAM_KEYS). A dropped key would warn and run at
    food_per=1, giving the same cumulative source as the default run.
    300 steps at K=8 fire 37 times (296 steps' worth of input), so the
    cumulative source must differ from the per-step run's 300 steps."""
    import warnings
    from evoca_explore import run_sim
    P = dict(food_inc=0.01, m_scale=1.5, gdiff=0.08, mu_lut=0.01,
             mu_egene=0.003, mu_egenome=0.03, tax=0.03)
    a = run_sim(P, n_steps=300, sample_every=50, N=32, seed=1)
    with warnings.catch_warnings():
        warnings.simplefilter('error')        # a dropped key would raise
        b = run_sim(dict(P, food_per=8), n_steps=300, sample_every=50,
                    N=32, seed=1)
    n_cells = 32 * 32
    assert a['diss_source_cum_final'] == pytest.approx(0.01 * n_cells * 300, rel=1e-5)
    assert b['diss_source_cum_final'] == pytest.approx(0.01 * 8 * n_cells * 37, rel=1e-5)


# ── 8. energy identity still closes with lumpy regeneration ───────────

def test_food_per_energy_identity(make_sim):
    """source_step == d(storage) + diss_step every step, food_per=K."""
    s = _full_sim(make_sim, food_per=6)
    d0 = s.get_dissipation()
    prev = d0['storage_env'] + d0['storage_priv']
    for _ in range(60):
        s.step()
        d = s.get_dissipation()
        cur = d['storage_env'] + d['storage_priv']
        assert d['source_step'] == pytest.approx(cur - prev + d['diss_step'], abs=1e-3)
        prev = cur
