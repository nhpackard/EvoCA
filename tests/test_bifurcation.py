"""Bifurcation-harness reducer + Feigenbaum-estimator correctness.

These tests exercise the *analysis* of the harness (the reducer and the
Feigenbaum helper) on synthetic series with known recurrent sets, so
they are fast and do not run the CA at all.  A separate light smoke test
runs one tiny real sweep to confirm the driver wiring is intact.
"""
import numpy as np
import pytest

import bifurcation as bif


# ── recurrent_extrema: fixed point ────────────────────────────────────

def test_reducer_fixed_point():
    """A series that decays to a constant has a single recurrent value."""
    t = np.arange(2000)
    series = 5.0 + 10.0 * np.exp(-t / 50.0)   # transient -> fixed pt 5.0
    pts = bif.recurrent_extrema(series, transient_frac=0.5)
    assert pts.size == 1
    assert pts[0] == pytest.approx(5.0, abs=1e-3)


def test_reducer_fixed_point_exact_constant():
    series = np.full(500, 42.0)
    pts = bif.recurrent_extrema(series, transient_frac=0.4)
    assert pts.size == 1
    assert pts[0] == pytest.approx(42.0)


# ── recurrent_extrema: period-2 ───────────────────────────────────────

def test_reducer_period_two():
    """An alternating signal has exactly two recurrent values."""
    base = np.tile([1.0, 9.0], 1000)          # period-2: extrema 1 and 9
    pts = np.sort(bif.recurrent_extrema(base, transient_frac=0.5))
    assert pts.size == 2
    assert pts[0] == pytest.approx(1.0, abs=1e-6)
    assert pts[1] == pytest.approx(9.0, abs=1e-6)


def test_reducer_period_two_with_transient_and_noise():
    rng = np.random.default_rng(0)
    n = 4000
    t = np.arange(n)
    # decaying transient then a clean period-2 oscillation
    osc = np.where(t % 2 == 0, 2.0, 6.0).astype(float)
    transient = 8.0 * np.exp(-t / 80.0)
    series = osc + transient + rng.normal(0, 1e-4, n)
    pts = np.sort(bif.recurrent_extrema(series, transient_frac=0.6))
    assert pts.size == 2
    assert pts[0] == pytest.approx(2.0, abs=5e-3)
    assert pts[1] == pytest.approx(6.0, abs=5e-3)


def test_reducer_period_four():
    """Period-4 signal -> four distinct recurrent values."""
    base = np.tile([1.0, 4.0, 2.0, 7.0], 800)
    pts = np.sort(bif.recurrent_extrema(base, transient_frac=0.5))
    # turning points of a 4-cycle are its local maxima/minima; the
    # repeating pattern yields the distinct extreme values.
    assert pts.size == 4
    assert np.allclose(pts, [1.0, 2.0, 4.0, 7.0], atol=1e-6)


# ── feigenbaum_delta: synthetic geometric cascade ─────────────────────

def test_feigenbaum_on_synthetic_cascade():
    """Thresholds whose gaps shrink geometrically by exactly delta must
    yield delta back."""
    delta = 4.669201609
    p0, gap0 = 1.0, 1.0
    thresholds = [p0]
    gap = gap0
    for _ in range(8):
        p0 = p0 + gap
        thresholds.append(p0)
        gap /= delta
    res = bif.feigenbaum_delta(thresholds)
    assert res['delta_estimate'] == pytest.approx(delta, rel=1e-6)
    assert res['delta_mean'] == pytest.approx(delta, rel=1e-6)
    assert res['n_thresholds'] == len(thresholds)


def test_feigenbaum_too_few_thresholds():
    res = bif.feigenbaum_delta([1.0, 2.0])
    assert res['delta_estimate'] is None
    assert res['delta_mean'] is None
    assert res['n_thresholds'] == 2


def test_doubling_thresholds_and_delta_roundtrip():
    """A fabricated sweep whose recurrent-set size doubles at
    geometrically converging param values recovers delta end-to-end."""
    delta = 4.669201609
    # doubling param locations whose successive gaps shrink by exactly
    # delta (a clean geometric cascade).
    locs = [0.0]
    gap = 1.0
    for _ in range(5):
        locs.append(locs[-1] + gap)
        gap /= delta
    # Build a sweep whose recurrent-set size is 1 before the first
    # location and doubles 1->2->4->... exactly at each location.  One
    # result row per cascade point plus a leading size-1 row, so
    # doubling_thresholds recovers the locations exactly.
    results = [{'value': -0.5, 'n_recurrent': 1, 'points': np.zeros(1)}]
    size = 1
    for loc in locs:
        size *= 2
        results.append({'value': float(loc), 'n_recurrent': size,
                        'points': np.zeros(size)})
    th = bif.doubling_thresholds(results)
    # thresholds recovered should be exactly the cascade locations
    assert th.size == len(locs)
    assert np.allclose(th, locs)
    res = bif.feigenbaum_delta(th)
    assert res['delta_estimate'] == pytest.approx(delta, rel=1e-6)


# ── to_scatter shape ──────────────────────────────────────────────────

def test_to_scatter_flattens():
    results = [
        {'value': 0.1, 'points': np.array([1.0, 2.0])},
        {'value': 0.2, 'points': np.array([3.0])},
    ]
    xs, ys = bif.to_scatter(results)
    assert xs.tolist() == [0.1, 0.1, 0.2]
    assert ys.tolist() == [1.0, 2.0, 3.0]


# ── light real-driver smoke test (fast: tiny N, few ticks) ────────────

def test_sweep_driver_smoke():
    """One tiny real sweep: confirms the driver/sim wiring is intact and
    reproducible at fixed seed.  Kept very small (seconds)."""
    vals = [0.05, 0.15]
    r1 = bif.sweep('food_inc', vals, observable='pop',
                   n_steps=40, N=32, seed=0, transient_frac=0.5)
    assert len(r1) == 2
    for r in r1:
        assert r['n_recurrent'] >= 1
        assert np.isfinite(r['mean'])
        assert r['points'].size == r['n_recurrent']
    # fixed seed -> identical recurrent sets on a repeat run
    r2 = bif.sweep('food_inc', vals, observable='pop',
                   n_steps=40, N=32, seed=0, transient_frac=0.5)
    for a, b in zip(r1, r2):
        assert np.array_equal(a['points'], b['points'])
