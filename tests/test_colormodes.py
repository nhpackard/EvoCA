"""Color modes: the shared Python list matches the C side, and the
'neighborhood' mode (5) colors each cell by its 3x3 Moore pattern of v."""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'python'))

from colormodes import COLOR_MODES          # noqa: E402
from evoca_py import EvoCA                  # noqa: E402

BLACK = np.int32(-16777216)                 # 0xFF000000


def _sim(N):
    sim = EvoCA()
    sim.init(N, food_inc=0.0, m_scale=1.0)
    return sim


def _patterns(v):
    """Expected 9-bit Moore pattern per cell, bit (dy+1)*3+(dx+1), periodic."""
    p = np.zeros_like(v, dtype=np.int64)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            bit = (dy + 1) * 3 + (dx + 1)
            p |= np.roll(np.roll(v, -dy, 0), -dx, 1).astype(np.int64) << (8 - bit)
    return p


def test_python_list_matches_c():
    sim = _sim(8)
    assert len(COLOR_MODES) == sim._lib.evoca_colormode_count()
    assert COLOR_MODES[0] == 'genome' and COLOR_MODES[5] == 'neighborhood'


def test_neighborhood_single_cell_bit_order_and_wrap():
    N = 8
    sim = _sim(N)
    v = np.zeros((N, N), np.uint8)
    v[0, 0] = 1                              # corner: exercises periodic wrap
    sim.set_v(v)
    px = sim.colorize(np.zeros(N * N, np.int32), 5).reshape(N, N)
    exp = _patterns(v)
    assert set(np.unique(exp[exp > 0])) == {1 << b for b in range(9)}
    assert np.all(px[exp == 0] == BLACK)     # all-zero pattern is black
    lit = px[exp > 0]
    assert np.all(lit != BLACK) and len(np.unique(lit)) == 9


def test_neighborhood_color_is_a_function_of_pattern_and_injective():
    N = 64
    sim = _sim(N)
    v = np.random.default_rng(0).integers(0, 2, (N, N), dtype=np.uint8)
    sim.set_v(v)
    px = sim.colorize(np.zeros(N * N, np.int32), 5).reshape(N, N)
    pat = _patterns(v)
    color_of = {}
    for p, c in zip(pat.ravel(), px.ravel()):
        assert color_of.setdefault(int(p), int(c)) == int(c)
    assert len(color_of) > 450               # most of the 512 patterns present
    assert len(set(color_of.values())) == len(color_of)
    assert all(c != BLACK for p, c in color_of.items() if p != 0)
