"""Tax / survival thresholds.

The death rule is: f(x) -= tax each step, die when f reaches 0, and
food_inc=0 means no eating can offset it. So with initial private
food f0 and constant tax t (and no other tax terms), every cell must
be dead by ceil(f0 / t) steps, and must still be alive strictly
before that. tax=0 with no death pressure keeps the whole population
alive indefinitely. These bracket the survival arithmetic.
"""
import math

import numpy as np


def _world(make_sim, N, tax, f0):
    s = make_sim(N, mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0,
                 tax=tax, food_inc=0.0, m_scale=0.0)
    from evoca_py import make_gol_lut
    s.set_lut_all(make_gol_lut())
    s.set_alive_all()
    s.set_v(np.zeros((N, N), dtype=np.uint8))
    s.set_f_all(f0)
    return s


def test_high_tax_kills_all_by_threshold(make_sim):
    N = 16
    tax, f0 = 0.1, 0.5
    s = _world(make_sim, N, tax, f0)
    deadline = math.ceil(f0 / tax)          # = 5

    # Strictly before the deadline, the population is still fully alive.
    for _ in range(deadline - 1):
        s.step()
    assert s.get_alive().sum() == N * N, \
        "cells died earlier than the tax arithmetic allows"

    # By the deadline (one or two more steps), everything is dead.
    for _ in range(2):
        s.step()
    assert s.get_alive().sum() == 0, \
        "cells survived past the tax-depletion deadline"


def test_zero_tax_population_persists(make_sim):
    N = 16
    s = _world(make_sim, N, tax=0.0, f0=0.5)
    for _ in range(200):
        s.step()
    assert s.get_alive().sum() == N * N, \
        "cells died with tax=0 and no death pressure"


def test_dead_cell_state_is_zeroed(make_sim):
    """A dead cell must have v=0 and a zeroed LUT (structural
    invariant relied on by neighbour counting and activity)."""
    N = 16
    s = _world(make_sim, N, tax=0.1, f0=0.05)   # all die within ~1 step
    for _ in range(3):
        s.step()
    assert s.get_alive().sum() == 0
    assert s.get_v().sum() == 0, "dead cells must have v=0"
    lut = s.get_lut(0)
    assert int(np.asarray(lut).sum()) == 0, "dead cell LUT must be zeroed"
