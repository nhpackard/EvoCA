"""Paired total-variation distance for the fixed-space drift shadows
(N, 2026-11-03: compare each key with its own shadow instead of
ΣG/D_G − ΣN/D_N).

    TV = ½ Σ_k |G_k − N_k| / Σ_k G_k  ∈ [0, 1]

G_k, N_k are the real and drift-shadow cumulative activities of key k
(eg: 729 ternary egene keys; dyn: 500 LUT transition buckets). The
shadow receives the same total mass as the real run every step, so the
plain paired mean Σ(G_k − N_k) is identically 0 and the absolute value
is what carries the signal.

Asserted here:
  * get_eg_shadow / get_dyn_shadow return the per-key arrays, with equal
    totals (the conservation the definition relies on);
  * eg_tv / dyn_tv (computed in C) equal TV recomputed from those arrays;
  * a frozen arm (zero flux) gives TV = 0;
  * a selective arm gives TV clearly above 0;
  * run_sim emits eg_tv_* and dyn_tv_* (public-boundary reachability).
"""
import math
import os
import sys

import numpy as np

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))

EG_NB = 729
DYN_NB = 500

FROZEN = dict(mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0, tax=0.002,
              tax_per_egene=0.0, food_inc=0.5, gdiff=0.0)
SELECTIVE = dict(mu_lut=0.005, mu_egene=0.004, mu_egenome=0.004, tax=0.02,
                 tax_per_egene=0.006, food_inc=0.06, gdiff=0.5)


def _run(make_sim, p, n_steps, N=40, seed=0):
    sim = make_sim(N, m_scale=1.2, **p)
    sim.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    rng = np.random.default_rng(seed)
    sim.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
    sim.set_f_all(0.5)
    sim.set_F_all(1.0)
    sim.set_alive_fraction(0.6)
    for _ in range(n_steps):
        sim.step()
    assert int(sim.get_alive().sum()) > 0, 'run went extinct'
    return sim


def _tv(pair):
    g, n = pair['obs'], pair['shadow']
    return 0.5 * np.abs(g - n).sum() / g.sum() if g.sum() > 0 else 0.0


def test_shadow_arrays_and_conservation(make_sim):
    sim = _run(make_sim, SELECTIVE, 300, seed=3)
    for pair, K in ((sim.get_eg_shadow(), EG_NB), (sim.get_dyn_shadow(), DYN_NB)):
        assert pair['obs'].shape == (K,) and pair['shadow'].shape == (K,)
        assert pair['obs'].sum() > 0
        # Same mass per step on both sides; only the sub-unit carry and
        # the renormalization rounding separate the totals.
        assert abs(pair['obs'].sum() - pair['shadow'].sum()) <= 1e-6 * pair['obs'].sum() + K
    assert math.isclose(pair['shadow'].sum(), sim.dyn_shadow_total(), rel_tol=0, abs_tol=K)


def test_c_tv_matches_arrays(make_sim):
    sim = _run(make_sim, SELECTIVE, 300, seed=5)
    assert math.isclose(sim.eg_tv(), _tv(sim.get_eg_shadow()), rel_tol=1e-12, abs_tol=1e-15)
    assert math.isclose(sim.dyn_tv(), _tv(sim.get_dyn_shadow()), rel_tol=1e-12, abs_tol=1e-15)


def test_frozen_tv_is_zero(make_sim):
    for seed in (3, 5, 7):
        sim = _run(make_sim, FROZEN, 300, seed=seed)
        assert sim.eg_tv() < 1e-6, f'frozen eg TV {sim.eg_tv()} (seed {seed})'
        assert sim.dyn_tv() < 1e-6, f'frozen dyn TV {sim.dyn_tv()} (seed {seed})'


def test_selective_tv_positive(make_sim):
    for seed in (3, 5, 7):
        sim = _run(make_sim, SELECTIVE, 400, seed=seed)
        assert 0.01 < sim.eg_tv() <= 1.0, f'eg TV {sim.eg_tv()} (seed {seed})'
        assert 0.01 < sim.dyn_tv() <= 1.0, f'dyn TV {sim.dyn_tv()} (seed {seed})'


def test_run_sim_emits_tv():
    from evoca_explore import run_sim
    out = run_sim(dict(food_inc=0.013, m_scale=1.2, gdiff=0.06, tax=0.035,
                       mu_lut=0.01, mu_egene=0.004, mu_egenome=0.004),
                  n_steps=300, sample_every=50, N=40, seed=0, shadow=False)
    for k in ('eg_tv_final', 'eg_tv_mean', 'dyn_tv_final', 'dyn_tv_mean'):
        assert k in out, f'run_sim does not emit {k}'
        assert 0.0 <= out[k] <= 1.0, f'{k} = {out[k]}'
    assert out['dyn_tv_final'] > 0.0
