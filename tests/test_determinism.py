"""Determinism guard for the S1 RNG-reseed fix.

The C-side xorshift32 PRNG is a per-process global that was never
reset. multiprocessing.Pool workers are reused, so the C RNG state
carried across configs and per-config results were not reproducible.
run_sim now calls sim.set_seed(seed) so the C RNG is reseeded
identically every call.

These tests pin that invariant:
  - same seed (same process) => identical key metrics
  - different seeds          => different results (seed actually bites)
"""
import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))

_KEYS = (
    'n_distinct_genomes_mean',
    'cog_specificity_mean',
    'final_pop',
    'excess_pc_slope',
)

_CFG = dict(food_inc=0.013, m_scale=1.2, gdiff=0.06, tax=0.035,
            mu_lut=0.01, mu_egene=0.003, mu_egenome=0.005)


def _run(seed):
    from evoca_explore import run_sim
    out = run_sim(dict(_CFG), n_steps=300, sample_every=100,
                  N=48, seed=seed, shadow=True)
    return {k: out[k] for k in _KEYS}


def test_same_seed_is_reproducible():
    a = _run(7)
    b = _run(7)
    assert a == b, f"same-seed runs diverged: {a} != {b}"


def test_different_seeds_differ():
    a = _run(1)
    b = _run(2)
    assert a != b, f"different seeds gave identical metrics: {a}"
