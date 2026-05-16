"""Regression guard for the roadmap-§2 metric fix.

run_sim must emit BOTH the legacy raw excess and the new
component-normalised excess, and the new one must be finite. This
locks in the fix so a future refactor can't silently drop it (the
exact failure mode the original bug exploited).
"""
import math
import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))


def test_run_sim_emits_component_normalised_excess():
    from evoca_explore import run_sim
    out = run_sim(
        dict(food_inc=0.013, m_scale=1.2, gdiff=0.06, tax=0.035,
             mu_lut=0.01, mu_egene=0.003, mu_egenome=0.005),
        n_steps=300, sample_every=100, N=48, seed=0, shadow=True)

    for k in ('excess_activity_slope', 'excess_activity_final',
              'excess_pc_slope', 'excess_pc_final'):
        assert k in out, f"run_sim no longer emits {k}"
        assert math.isfinite(out[k]), f"{k} is not finite: {out[k]}"
