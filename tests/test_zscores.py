"""Per-bucket z-scores for the fixed-space neutral shadows.

Docs/checkpoint 3 July 2026.md §G / §I′ action G2: the fixed-space
neutral baselines (eg = 729 ternary keys, dyn = 500 LUT input/output
buckets) previously exposed only an AGGREGATE per-component-normalised
excess scalar.  This adds PER-BUCKET z-scoring — the well-posed
"per-gene neutral shadow" the user asked for.  For each bucket k:

    z[k] = (observed cumulative activity[k] − shadow drift mean[k])
           / sqrt(closed-form Wright-Fisher drift variance[k])

where the drift mean is the shadow's cumulative activity and the drift
variance is the T²·qn(1−qn)/Np WF form the realised-flux shadow already
implies, accumulated over steps.

Headline properties asserted here:

  * API surface: get_eg_zscores() / get_dyn_zscores() return finite
    float64 arrays of length 729 / 500.
  * Null/frozen arm (zero realised mutation flux): the shadow tracks the
    realised composition exactly, so EVERY bucket's z is ~0.
  * Selection arm: some buckets carry clearly large |z| — the real run
    concentrates cumulative activity on far fewer buckets than
    fitness-blind drift, and leaves drift-warm buckets it never visits.

Small and fast: N=40, a few hundred steps.
"""
import math
import os
import sys

import numpy as np

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))

EG_NB = 729
DYN_NB = 500


def _run(make_sim, *, mu_lut, mu_egene, mu_egenome, tax, tax_per_egene,
         food_inc, gdiff, n_steps, N=40, seed=0):
    """Drive a sim directly; return (alive, eg_zscores, dyn_zscores)."""
    sim = make_sim(N, food_inc=food_inc, m_scale=1.2, gdiff=gdiff,
                   mu_lut=mu_lut, mu_egene=mu_egene,
                   mu_egenome=mu_egenome, tax=tax,
                   tax_per_egene=tax_per_egene)
    sim.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    rng = np.random.default_rng(seed)
    sim.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
    sim.set_f_all(0.5)
    sim.set_F_all(1.0)
    sim.set_alive_fraction(0.6)
    for _ in range(n_steps):
        sim.step()
    alive = int(np.ctypeslib.as_array(
        sim._lib.evoca_get_alive(), shape=(N * N,)).sum())
    return alive, sim.get_eg_zscores(), sim.get_dyn_zscores()


def test_zscore_api_surface(make_sim):
    """get_eg_zscores/get_dyn_zscores return finite float64 arrays of
    the fixed-space bucket counts (729 / 500)."""
    alive, egz, dyz = _run(
        make_sim,
        mu_lut=0.01, mu_egene=0.004, mu_egenome=0.004,
        tax=0.035, tax_per_egene=0.0,
        food_inc=0.013, gdiff=0.06,
        n_steps=200, seed=0)
    assert egz.shape == (EG_NB,) and egz.dtype == np.float64
    assert dyz.shape == (DYN_NB,) and dyz.dtype == np.float64
    assert np.all(np.isfinite(egz)) and np.all(np.isfinite(dyz))


def test_null_frozen_zscores_are_zero(make_sim):
    """Genuinely neutral/frozen arm: ZERO mutation (no realised flux) +
    abundant uniform food so the population persists.  With no flux the
    drift shadow tracks the realised composition EXACTLY, so the
    per-bucket numerator (observed − drift mean) is 0 and every z is 0.
    This is the per-bucket analogue of the excess≈0 neutral property."""
    for seed in (3, 5, 7):
        alive, egz, dyz = _run(
            make_sim,
            mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0,
            tax=0.002, tax_per_egene=0.0,
            food_inc=0.5, gdiff=0.0,        # food saturates → persists
            n_steps=300, seed=seed)
        assert alive > 0, f"neutral run went extinct (seed {seed})"
        assert np.all(np.isfinite(egz)) and np.all(np.isfinite(dyz))
        # Zero realised flux ⇒ shadow == observed ⇒ every z ≈ 0.
        assert np.abs(egz).max() < 1e-6, (
            f"frozen eg z not ~0 (seed {seed}): {np.abs(egz).max()}")
        assert np.abs(dyz).max() < 1e-6, (
            f"frozen dyn z not ~0 (seed {seed}): {np.abs(dyz).max()}")


def test_selection_zscores_clearly_large(make_sim):
    """Strong-selection regime: scarce diffused food + a per-egene tax
    put eating/rule genomes under real selection, so the real run
    concentrates cumulative activity on far fewer buckets than
    fitness-blind drift.  Some buckets must therefore carry clearly
    large |z| — far above the frozen (all-zero) arm.  This is the
    discriminating property the per-bucket shadow exists to provide."""
    for seed in (3, 5, 7):
        a_s, egz, dyz = _run(
            make_sim,
            mu_lut=0.005, mu_egene=0.004, mu_egenome=0.004,
            tax=0.02, tax_per_egene=0.006,
            food_inc=0.06, gdiff=0.5,
            n_steps=400, seed=seed)
        assert a_s > 0, f"selective run went extinct (seed {seed})"
        assert np.all(np.isfinite(egz)) and np.all(np.isfinite(dyz))
        # Under selection at least one bucket in each space is a strong
        # outlier vs its drift shadow (|z| well above any plausible
        # neutral fluctuation).  Thresholds are conservative — observed
        # maxima are ~15 (eg) and ~600 (dyn).
        assert np.abs(egz).max() > 3.0, (
            f"selected eg z not clearly large (seed {seed}): "
            f"{np.abs(egz).max()}")
        assert np.abs(dyz).max() > 3.0, (
            f"selected dyn z not clearly large (seed {seed}): "
            f"{np.abs(dyz).max()}")
        # And a positive concentration signal specifically (selected
        # buckets accrue activity faster than drift).
        assert egz.max() > 2.0, f"no positive eg z outlier (seed {seed})"
        assert dyz.max() > 2.0, f"no positive dyn z outlier (seed {seed})"
