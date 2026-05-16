"""Conway's Game of Life correctness.

With the GoL LUT loaded uniformly and all mutation/tax/food off,
EvoCA must reproduce exact B3/S23 on the v-plane. These are the
strongest "is the CA core correct?" anchors: oscillator period,
still-life invariance, and glider translation are independent
failure modes.
"""
import numpy as np
import pytest


def _gol_world(make_sim, N):
    from evoca_py import make_gol_lut
    s = make_sim(N, mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0,
                 tax=0.0, food_inc=0.0)
    s.set_lut_all(make_gol_lut())
    s.set_alive_all()          # all cells alive so v evolves as pure GoL
    return s


def _set_cells(s, N, coords):
    v = np.zeros((N, N), dtype=np.uint8)
    for (r, c) in coords:
        v[r, c] = 1
    s.set_v(v)


def test_blinker_period_2(make_sim):
    """A 3-cell blinker oscillates with period 2."""
    N = 16
    s = _gol_world(make_sim, N)
    horiz = [(8, 7), (8, 8), (8, 9)]
    vert = [(7, 8), (8, 8), (9, 8)]
    _set_cells(s, N, horiz)

    s.step()
    v1 = s.get_v()
    assert sorted(zip(*np.where(v1 == 1))) == sorted(vert), \
        "blinker did not flip horizontal->vertical after 1 step"

    s.step()
    v2 = s.get_v()
    assert sorted(zip(*np.where(v2 == 1))) == sorted(horiz), \
        "blinker did not return to horizontal after 2 steps"


def test_block_still_life(make_sim):
    """A 2x2 block is invariant under GoL for many steps."""
    N = 16
    s = _gol_world(make_sim, N)
    block = [(7, 7), (7, 8), (8, 7), (8, 8)]
    _set_cells(s, N, block)
    for _ in range(10):
        s.step()
    v = s.get_v()
    assert sorted(zip(*np.where(v == 1))) == sorted(block), \
        "block still-life was not preserved"


def test_glider_translates(make_sim):
    """A glider returns to its shape translated by (1,1) after 4 steps."""
    N = 24
    s = _gol_world(make_sim, N)
    glider = [(5, 6), (6, 7), (7, 5), (7, 6), (7, 7)]
    _set_cells(s, N, glider)
    for _ in range(4):
        s.step()
    v = s.get_v()
    got = sorted(zip(*np.where(v == 1)))
    expected = sorted([(r + 1, c + 1) for (r, c) in glider])
    assert got == expected, f"glider did not translate by (1,1): {got}"


def test_gol_population_stable_without_mutation(make_sim):
    """With mu=0, every cell keeps the wild-type GoL LUT: the distinct
    LUT-genome count stays exactly 1 over a long run (guards against
    silent LUT corruption / accidental mutation paths)."""
    N = 32
    s = _gol_world(make_sim, N)
    rng = np.random.default_rng(0)
    s.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
    for _ in range(50):
        s.step()
    s._lib.evoca_activity_update()        # populate pop_count (as run_sim does)
    g = s.get_activity(max_n=4096)
    alive_buckets = int((g['pop_count'] > 0).sum())
    assert alive_buckets == 1, \
        f"expected 1 LUT genome with mu=0, got {alive_buckets}"
