"""Fixed-space neutral baselines for eg_activity / dyn_activity.

Docs/EvoCA_research_directions.md §2: eg_activity (729 ternary egene
keys) and dyn_activity (500 LUT (input,output) buckets) live in FIXED
finite bucket spaces, so their neutral models are closed-form
Wright-Fisher drift baselines driven by the *realised* per-component
mutation flux.  This closes the gap the LUT-hash Channon shadow leaves
for egene-driven runs: get_activity hashes the FULL genome while the
Channon shadow models LUT bytes only, so excess activity was undefined
for egene-driven runs.

The headline properties asserted here:

  * Under genuinely neutral drift (zero realised mutation flux) the
    per-component-normalised excess is ~0 within noise -- the shadow
    tracks the realised composition exactly, so observed - baseline
    vanishes.
  * Under strong selection the excess is clearly > 0 and far above
    the neutral case -- the real run concentrates cumulative activity
    on fewer buckets than fitness-blind drift does.

Small and fast: N=40, a few hundred steps.
"""
import math
import os
import sys

import numpy as np

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))


# ── API surface + finiteness ────────────────────────────────────────

def test_run_sim_emits_fixed_space_excess():
    """run_sim must emit eg_/dyn_ excess_pc final + slope, finite,
    mirroring the excess_pc_* outputs (additive, never dropped)."""
    from evoca_explore import run_sim
    out = run_sim(
        dict(food_inc=0.013, m_scale=1.2, gdiff=0.06, tax=0.035,
             mu_lut=0.01, mu_egene=0.004, mu_egenome=0.004),
        n_steps=300, sample_every=50, N=40, seed=0, shadow=True)

    for k in ('eg_excess_pc_final', 'eg_excess_pc_slope',
              'dyn_excess_pc_final', 'dyn_excess_pc_slope'):
        assert k in out, f"run_sim no longer emits {k}"
        assert math.isfinite(out[k]), f"{k} not finite: {out[k]}"


def test_excess_pc_independent_of_channon_shadow_flag():
    """The fixed-space drift shadows run inside sim.step(), so the
    eg_/dyn_ excess must be emitted even with shadow=False (LUT-hash
    Channon shadow disabled). This is exactly the gap §2 calls out:
    egene-driven excess must be defined regardless of the LUT shadow."""
    from evoca_explore import run_sim
    out = run_sim(
        dict(food_inc=0.013, m_scale=1.2, gdiff=0.06, tax=0.035,
             mu_lut=0.01, mu_egene=0.004, mu_egenome=0.004),
        n_steps=200, sample_every=50, N=40, seed=1, shadow=False)
    for k in ('eg_excess_pc_final', 'dyn_excess_pc_final'):
        assert k in out and math.isfinite(out[k])


# ── Direct C-level drift vs selection contrast ──────────────────────

def _run(make_sim, *, mu_lut, mu_egene, mu_egenome, tax, tax_per_egene,
         food_inc, gdiff, n_steps, N=40, seed=0):
    """Drive a sim directly; return (alive, eg_excess, dyn_excess)."""
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
    return alive, sim.eg_excess_pc(), sim.dyn_excess_pc()


def test_neutral_drift_excess_is_zero(make_sim):
    """Genuinely neutral regime: ZERO mutation (no realised flux) and
    abundant uniform food so the population persists. With no flux the
    drift shadow tracks the realised bucket composition exactly, so
    observed - baseline must be ~0 for BOTH fixed spaces. This is the
    defining property of a valid neutral model (excess ~= 0 within
    noise under neutral drift)."""
    for seed in (3, 5, 7):
        alive, eg_x, dyn_x = _run(
            make_sim,
            mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0,
            tax=0.002, tax_per_egene=0.0,
            food_inc=0.5, gdiff=0.0,        # food saturates → persists
            n_steps=300, seed=seed)
        assert alive > 0, f"neutral run went extinct (seed {seed})"
        assert math.isfinite(eg_x) and math.isfinite(dyn_x)
        # Zero realised flux ⇒ shadow == observed ⇒ excess ≈ 0. Allow a
        # tiny tolerance for FP/rounding-carry; the value is 0 exactly
        # in practice.
        assert abs(eg_x) < 1.0, (
            f"neutral eg excess not ~0 (seed {seed}): {eg_x}")
        assert abs(dyn_x) < 1.0, (
            f"neutral dyn excess not ~0 (seed {seed}): {dyn_x}")


def test_selection_excess_clearly_positive_and_above_neutral(make_sim):
    """Strong-selection regime vs the neutral regime: with scarce
    diffused food + a per-egene tax, eating/rule genomes are under
    real selection, so the real run concentrates cumulative activity
    on fewer buckets than fitness-blind drift. The per-component
    excess must be clearly positive and far above the neutral (zero-
    flux) value. This is the discriminating property the baseline
    exists to provide."""
    for seed in (3, 5, 7):
        # Neutral reference (zero flux, persists): excess ≈ 0.
        a_n, eg_n, dyn_n = _run(
            make_sim,
            mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0,
            tax=0.002, tax_per_egene=0.0,
            food_inc=0.5, gdiff=0.0,
            n_steps=400, seed=seed)

        # Selective regime: scarce diffused food + per-egene tax →
        # differential survival of eating/rule genomes.
        a_s, eg_s, dyn_s = _run(
            make_sim,
            mu_lut=0.005, mu_egene=0.004, mu_egenome=0.004,
            tax=0.02, tax_per_egene=0.006,
            food_inc=0.06, gdiff=0.5,
            n_steps=400, seed=seed)

        assert a_s > 0, f"selective run went extinct (seed {seed})"
        assert all(math.isfinite(v)
                   for v in (eg_n, dyn_n, eg_s, dyn_s))

        # Under selection the population concentrates activity onto
        # fewer high-activity buckets faster than the drift baseline,
        # so the per-component excess is clearly positive and far
        # above the near-zero neutral value.
        assert eg_s > 0.0, f"selected eg excess not > 0: {eg_s}"
        assert dyn_s > 0.0, f"selected dyn excess not > 0: {dyn_s}"
        assert eg_s > eg_n + 100.0, (
            f"selection eg excess ({eg_s:.1f}) not clearly above "
            f"neutral ({eg_n:.1f}), seed {seed}")
        assert dyn_s > dyn_n + 100.0, (
            f"selection dyn excess ({dyn_s:.1f}) not clearly above "
            f"neutral ({dyn_n:.1f}), seed {seed}")
