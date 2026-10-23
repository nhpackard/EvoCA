"""Milestone-1 regression gate for the `chem` branch.

Contract: at ``N_chem == 1`` the vector-valued-resource ("chem") build must
reproduce the pre-change single-food dynamics **bit for bit**. The golden
snapshot in ``tests/data/chem_parity_golden.npz`` was captured from the
pre-change build (via ``chem_parity.regen_golden()``); this test replays the
identical fixed-seed initial state on the current build and asserts exact
equality of the private food ``f``, environmental food ``F``, ``alive`` mask,
and per-cell LUT hashes at every checkpoint.

If this fails, the N_chem=1 code path has diverged from single-food behaviour
(most likely: an extra PRNG draw in the hot path, or a non-identity
vector->f reduction).
"""
import pytest
import os
import sys

_HERE = os.path.dirname(__file__)
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, '..', 'python'))

import chem_parity as cp  # noqa: E402


def test_n_chem_1_is_bit_for_bit_parity():
    from evoca_py import EvoCA
    golden = cp.load_golden()
    got = cp.run_and_capture(EvoCA())
    ok, lines = cp.compare(golden, got)
    assert ok, "N_chem=1 chem build diverged from single-food golden:\n" + \
        "\n".join(lines)


def test_n_chem_defaults_to_one_and_clamps():
    from evoca_py import EvoCA
    sim = EvoCA()
    sim.init(N=8, food_inc=0.01, m_scale=1.0)
    assert sim.get_n_chem() == 1
    # Milestone 2, C level: 1..CHEM_MAX_N (8) are accepted; larger sizes
    # clamp to 1 (surfaced).
    assert sim._lib.evoca_set_n_chem(2) == 1 and sim.get_n_chem() == 2
    assert sim._lib.evoca_set_n_chem(9) == 0
    assert sim.get_n_chem() == 1
    # Python level: only the configured chemistries (1, 2); others refused.
    assert sim.set_n_chem(2) == 2
    with pytest.raises(ValueError):
        sim.set_n_chem(9)
    assert sim.get_n_chem() == 2


def test_chem_energy_default_is_unit():
    from evoca_py import EvoCA
    sim = EvoCA()
    sim.init(N=8, food_inc=0.01, m_scale=1.0)
    assert sim.get_chem_energy(0) == 1.0


def test_stamp_patch_overwrite_keeps_chem_store_in_sync():
    """Regression: stamp_patch writes f_priv directly through the
    evoca_get_f() pointer. The chem store must be resynced from it, or
    the next tax silently restores the host cell's old private food
    into the stamped organism. Pre-chem behavior: stamped cells start
    at _STAMP_DEFAULT_F (0.0) and are taxed to death on the next step."""
    import numpy as np
    from evoca_py import EvoCA, _STAMP_DEFAULT_F
    N = 32
    src = EvoCA()
    src.init(N, food_inc=0.0, tax=0.01)
    src.state(lut='gol', egenome='uniform', egenome_value=3)
    src.set_f_all(0.5)
    src.set_alive_halfplane(0)
    patch = src.extract_patch(0, 0, 8)
    src.free()

    host = EvoCA()
    host.init(N, food_inc=0.0, m_scale=0.0, tax=0.01)
    host.state(lut='gol', egenome='uniform', egenome_value=3)
    host.set_F_all(0.0)
    host.set_f_all(0.8)                       # fed living host
    host.stamp_patch(patch, 10, 10, mode='overwrite')
    host.step()
    f = host.get_f()[10:18, 10:18]
    alive = host.get_alive().reshape(N, N)[10:18, 10:18]
    expected_f = max(0.0, _STAMP_DEFAULT_F - 0.01)
    host.free()
    assert np.allclose(f, expected_f), "stamped cells kept the host's food"
    assert int(alive.sum()) == (64 if expected_f > 0 else 0)


def test_n_chem_reaches_sim_via_run_sim(capfd):
    """run_sim must forward n_chem (filtered through _PARAM_KEYS): an
    unconfigured value reaches EvoCA.init, which refuses it; a dropped key
    would warn instead and the value would never arrive."""
    import warnings
    from evoca_explore import run_sim
    with warnings.catch_warnings():
        warnings.simplefilter('error')        # a dropped key would raise
        with pytest.raises(ValueError, match='n_chem=4'):
            run_sim(dict(food_inc=0.01, m_scale=1.5, tax=0.03, n_chem=4),
                    n_steps=20, sample_every=10, N=16, seed=1)


def test_n_chem_recipe_roundtrip():
    from evoca_py import EvoCA
    s = EvoCA()
    s.init(8, n_chem=1)
    assert s.params()['n_chem'] == 1
    assert s._init_metaparams['n_chem'] == 1
    s.free()


def test_chem_energy_non_unit_rejected_at_milestone_1(capfd):
    """Regression: a weight != 1 on compound 0 let f exceed 1 and broke
    the dissipation energy identity (intake/tax/headroom move amounts,
    not weighted energy). Compound 0 must keep weight 1 (milestone 2 lets
    other compounds be 0 or 1)."""
    from evoca_py import EvoCA
    sim = EvoCA()
    sim.init(N=8, food_inc=0.01, m_scale=1.0)
    sim.set_chem_energy(0, 2.0)
    assert sim.get_chem_energy(0) == 1.0
    assert 'unsupported for compound 0' in capfd.readouterr().err
    sim.free()
