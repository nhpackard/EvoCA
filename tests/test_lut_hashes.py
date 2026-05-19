"""Accessor guard for #8's enabling API: sim.get_lut_hashes().

The activity tracker / lineage_parent_hash expose the *combined*
LUT‖egene FNV-1a genome hash (lut_hash_cache). #8 (direct lineage
co-evolution metric) must decompose a parent→child change into ΔLUT
*independently* of Δegene, so it needs a per-cell hash of the LUT
bytes alone. evoca_get_lut_hashes() provides exactly that.

These tests pin the defining properties:
  - egene-INVARIANT: changing only the egenome must not change any
    lut-hash (this is the property the combined hash fails — the whole
    reason the accessor exists).
  - LUT-SENSITIVE: changing the LUT changes the hash.
  - CORRECT: equals an independent Python re-implementation of the
    C lut_hash_fn (FNV-1a over LUT_BYTES/4 little-endian uint32 words).
"""
import os
import sys

import numpy as np

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))

from evoca_py import LUT_BYTES, make_gol_lut   # noqa: E402


def _fnv1a_lut(lut_bytes):
    """Python mirror of C lut_hash_fn: FNV-1a over LUT_BYTES/4 uint32
    little-endian words, 32-bit wrap."""
    w = np.frombuffer(bytes(lut_bytes[:LUT_BYTES]), dtype='<u4')
    h = 0x811c9dc5
    for x in w:
        h = ((h ^ int(x)) * 0x01000193) & 0xFFFFFFFF
    return h


def test_lut_hash_is_egene_invariant_and_correct(make_sim):
    N = 8
    s = make_sim(N)
    gol = make_gol_lut()
    s.set_lut_all(gol)
    s.set_egenome_all(0b000011)

    h0 = s.get_lut_hashes()
    assert h0.shape == (N, N)
    # uniform LUT => one distinct hash, matching the Python FNV mirror
    expect = _fnv1a_lut(gol)
    assert np.all(h0 == expect), "uniform GoL LUT must hash uniformly to the FNV mirror"

    # CORE PROPERTY: change only the egenome -> lut-hash must not move
    s.set_egenome_all(0b111111)
    h1 = s.get_lut_hashes()
    assert np.array_equal(h0, h1), "lut-hash leaked egene state (not LUT-only)"


def test_lut_hash_tracks_lut_change(make_sim):
    N = 8
    s = make_sim(N)
    s.set_lut_all(make_gol_lut())
    s.set_egenome_all(0)
    h_gol = s.get_lut_hashes()

    rng = np.random.default_rng(0)
    rnd = rng.integers(0, 256, size=LUT_BYTES, dtype=np.uint8)
    s.set_lut_all(rnd)
    h_rnd = s.get_lut_hashes()

    assert np.all(h_rnd != h_gol), "lut-hash did not respond to a LUT change"
    assert np.all(h_rnd == _fnv1a_lut(rnd)), "lut-hash != independent FNV mirror"
