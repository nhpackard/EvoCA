"""Opt-in per-cell lineage record (research-directions §1).

The lineage field records, for every Phase-4 birth, the parent's
full-genome FNV-1a hash plus monotonic birth ids. These tests pin two
properties:

1. Default OFF: accessors return all-zero and nothing is recorded,
   so the feature is genuinely opt-in (zero observable effect).
2. With lineage ON and mutation OFF (one shared genome), every recorded
   child parent-hash equals the single live genome hash reported by the
   activity tracker — i.e. the recorded parent identity is correct.
3. A controlled single reproduction (one alive founder, dead
   neighbourhood) produces exactly one child whose parent-hash is the
   founder's genome hash, whose birth id is 1, and whose parent id is 0
   (the founder was never born under tracking).

Additivity guard: enabling lineage must not perturb the CA — a GoL
oscillator still has the right period with lineage ON.
"""
import numpy as np


def _gol_world(make_sim, N, lineage=False, **kw):
    from evoca_py import make_gol_lut
    s = make_sim(N, mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0,
                 tax=0.0, food_inc=0.0, lineage=lineage, **kw)
    s.set_lut_all(make_gol_lut())
    return s


def test_default_off_is_zero_cost_and_zero_data(make_sim):
    """Without lineage=True the accessors return all zeros and the
    feature reports disabled."""
    s = _gol_world(make_sim, 32, lineage=False)
    assert s.lineage_enabled() is False
    s.set_alive_all()
    s.set_egenome_all(0)
    s.set_f_all(2.0)            # force reproduction everywhere
    for _ in range(3):
        s.step()
    ph = s.get_lineage_parent_hash()
    bid = s.get_lineage_birth_id()
    pid = s.get_lineage_parent_id()
    assert ph.shape == (32, 32)
    assert not ph.any(), "parent-hash nonzero with lineage OFF"
    assert not bid.any(), "birth-id nonzero with lineage OFF"
    assert not pid.any(), "parent-id nonzero with lineage OFF"


def test_parent_hash_matches_single_live_genome(make_sim):
    """Lineage ON, mutation OFF, one uniform genome: every recorded
    child parent-hash equals the one live-genome hash from the activity
    tracker."""
    N = 24
    s = _gol_world(make_sim, N, lineage=True)
    assert s.lineage_enabled() is True
    s.set_alive_all()
    s.set_egenome_all(0)
    s.set_f_all(2.0)           # everyone over reproduction threshold
    s.step()                   # one wave of births, no mutation

    s._lib.evoca_activity_update()
    act = s.get_activity()
    live = act['pop_count'] > 0
    assert live.sum() == 1, \
        "mutation off should leave exactly one live genome bucket"
    genome_hash = np.uint32(act['hash'][live][0])

    ph = s.get_lineage_parent_hash()
    bid = s.get_lineage_birth_id()
    born = bid > 0
    assert born.any(), "no births recorded with f=2.0 everywhere"
    rec = np.unique(ph[born])
    assert rec.tolist() == [genome_hash], \
        f"recorded parent hashes {rec} != live genome {genome_hash}"


def test_known_single_reproduction_link(make_sim):
    """One alive founder in a dead field, f above threshold: exactly
    one child is born; its parent-hash is the founder's genome hash,
    its birth id is 1, and its parent id is 0 (founder never tracked)."""
    N = 16
    s = _gol_world(make_sim, N, lineage=True)

    alive = np.zeros((N, N), dtype=np.uint8)
    alive[8, 8] = 1
    s.set_lut_all(__import__('evoca_py').make_gol_lut())
    s.set_egenome_all(0)
    s.set_alive(alive)         # zeroes dead cells' data
    s.set_f_all(2.0)           # founder (and dead cells) get f=2.0;
                               # only the alive founder can reproduce

    # Founder genome hash, before any birth.
    s._lib.evoca_activity_update()
    act = s.get_activity()
    live = act['pop_count'] > 0
    assert live.sum() == 1
    founder_hash = np.uint32(act['hash'][live][0])

    s.step()

    bid = s.get_lineage_birth_id()
    ph = s.get_lineage_parent_hash()
    pid = s.get_lineage_parent_id()

    born_idx = np.argwhere(bid > 0)
    assert born_idx.shape[0] == 1, \
        f"expected exactly one birth, got {born_idx.shape[0]}"
    r, c = born_idx[0]

    # Child must be a Moore neighbour of the founder (8,8).
    dr = min((r - 8) % N, (8 - r) % N)
    dc = min((c - 8) % N, (8 - c) % N)
    assert dr <= 1 and dc <= 1 and not (dr == 0 and dc == 0), \
        f"child at ({r},{c}) is not a Moore neighbour of founder"

    assert bid[r, c] == 1, f"first birth id should be 1, got {bid[r, c]}"
    assert pid[r, c] == 0, \
        f"founder parent id should be 0, got {pid[r, c]}"
    assert np.uint32(ph[r, c]) == founder_hash, \
        f"child parent-hash {ph[r, c]} != founder hash {founder_hash}"


def test_lineage_on_does_not_perturb_ca(make_sim):
    """Additivity: a GoL blinker still has period 2 with lineage ON
    (the lineage write must not touch CA state)."""
    N = 16
    s = _gol_world(make_sim, N, lineage=True)
    s.set_alive_all()          # all alive so v evolves as pure GoL
    v = np.zeros((N, N), dtype=np.uint8)
    for (rr, cc) in [(8, 7), (8, 8), (8, 9)]:
        v[rr, cc] = 1
    s.set_v(v)
    s.step()
    v1 = s.get_v()
    s.step()
    v2 = s.get_v()
    expect1 = {(7, 8), (8, 8), (9, 8)}
    assert set(map(tuple, np.argwhere(v1 == 1))) == expect1, \
        "blinker did not flip with lineage ON"
    assert set(map(tuple, np.argwhere(v2 == 1))) == {(8, 7), (8, 8), (8, 9)}, \
        "blinker did not return to phase-0 with lineage ON"
