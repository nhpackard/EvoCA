"""Effective-rule activity (eff_activity) — checkpoint §I' action F1.

Full-genome `activity` keys on an FNV-1a hash of the whole 32-byte LUT,
so flipping a LUT bit for an input the cell never encounters mints a
"new species" with no behavioural change (silent churn). `eff_activity`
keys instead on the LUT MASKED to the entries the population actually
exercises (the union of the per-cell used-masks accumulated in Phase 1),
folded with the same egene stream. This makes the species key a pure
coarsening of genome space:

  1. silent churn (flipping an entry no cell ever exercises) mints NO
     new effective-rule species, though it does mint a full-genome one;
  2. the alive eff-species count is <= the alive full-genome count;
  3. the CA dynamics are bit-identical (additive no-op).

Design note: the key uses the POPULATION-UNION mask, not each cell's own
per-cell mask. A per-cell key would split one genome across cells with
divergent neighbourhood histories, making eff-species LARGER than
full-species and defeating the collapse goal (see C/evoca.c
eff_union_mask comment). The per-cell masks are still accumulated and
exposed via get_eff_used_mask as the diagnostic substrate.
"""
import numpy as np
from collections import defaultdict


def _gol_world(make_sim, N, **kw):
    from evoca_py import make_gol_lut
    s = make_sim(N, mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0,
                 tax=0.0, food_inc=0.0, **kw)
    s.set_lut_all(make_gol_lut())
    s.set_alive_all()
    s.set_egenome_all(0)
    return s


def test_gol_dynamics_unchanged(make_sim):
    """Additivity guard: with the eff_activity machinery compiled in, a
    GoL blinker still oscillates with period 2 (the used-mask
    accumulation must not perturb the CA)."""
    N = 16
    s = _gol_world(make_sim, N)
    horiz = [(8, 7), (8, 8), (8, 9)]
    vert = [(7, 8), (8, 8), (9, 8)]
    v = np.zeros((N, N), np.uint8)
    for (r, c) in horiz:
        v[r, c] = 1
    s.set_v(v)
    s.step()
    assert sorted(zip(*np.where(s.get_v() == 1))) == sorted(vert)
    s.step()
    assert sorted(zip(*np.where(s.get_v() == 1))) == sorted(horiz)


def test_existing_activity_unchanged_and_eff_collapses(make_sim):
    """Under GoL with mu=0 the full-genome activity still reports exactly
    one species (existing behaviour), and eff_activity likewise collapses
    the uniform genome to a single effective-rule species."""
    N = 32
    s = _gol_world(make_sim, N)
    rng = np.random.default_rng(0)
    s.set_v(rng.integers(0, 2, (N, N), np.uint8))
    for _ in range(50):
        s.step()
    s._lib.evoca_activity_update()
    assert s.n_distinct_genomes() == 1
    s.eff_activity_update()
    assert s.eff_count_distinct() == 1


def test_silent_unused_bit_churn_mints_no_eff_species(make_sim):
    """The headline property: randomising each cell's globally-unused LUT
    entries mints many NEW full-genome species but ZERO new effective-
    rule species — every alive cell's eff hash is unchanged, so silent
    churn is invisible to eff_activity."""
    from evoca_py import LUT_BITS
    N = 32
    s = _gol_world(make_sim, N)
    rng = np.random.default_rng(3)
    s.set_v(rng.integers(0, 2, (N, N), np.uint8))
    for _ in range(40):            # accumulate the population union mask
        s.step()
    s._lib.evoca_activity_update()
    s.eff_activity_update()

    eff_before = s.get_eff_hashes().copy()
    full_before = s.n_distinct_genomes()
    eff_cd_before = s.eff_count_distinct()

    union = s.get_eff_union_mask()
    unused = [b for b in range(LUT_BITS)
              if not ((union[b >> 3] >> (b & 7)) & 1)]
    assert len(unused) > 0, "test needs some globally-unexercised LUT entries"

    alive = s.get_alive().astype(bool)
    for i in range(N * N):
        r, c = divmod(i, N)
        if not alive[r, c]:
            continue
        lb = s.get_lut(i).copy()
        cr = np.random.default_rng(10_000 + i)
        for b in unused:            # randomise ONLY globally-unused entries
            if cr.integers(0, 2):
                lb[b >> 3] ^= np.uint8(1 << (b & 7))
        s.set_lut(i, lb)           # NB: no step() — keeps those entries unused

    s._lib.evoca_activity_update()
    s.eff_activity_update()

    am = alive
    eff_after = s.get_eff_hashes()
    full_after = s.n_distinct_genomes()
    eff_cd_after = s.eff_count_distinct()

    # 1. silent churn changes NO cell's effective-rule hash
    assert np.array_equal(eff_after[am], eff_before[am])
    # 2. so the eff-species count is unchanged by the churn
    assert eff_cd_after == eff_cd_before
    # 3. the SAME churn minted many new full-genome species
    assert full_after > full_before
    # 4. eff-species count <= full-genome species count (by a wide margin)
    assert eff_cd_after <= full_after


def test_eff_is_a_coarsening_of_genome(make_sim):
    """Rigorous (hash-aliasing-immune) form of eff <= full: in a churning
    population, every set of cells sharing an EXACT genome maps to a
    single effective-rule hash. So eff is a function of the genome and
    the alive eff-species count cannot exceed the alive genome count."""
    N = 40
    s = make_sim(N, mu_lut=0.02, mu_egene=0.0, mu_egenome=0.0,
                 tax=0.0, food_inc=0.0)
    from evoca_py import make_gol_lut
    s.set_lut_all(make_gol_lut())
    s.set_alive_all()
    s.set_egenome_all(0)           # uniform egenome => genome identity == LUT bytes
    s.set_f_all(2.0)
    rng = np.random.default_rng(0)
    s.set_v(rng.integers(0, 2, (N, N), np.uint8))
    for _ in range(25):
        s.step()
        s.set_f_all(2.0)           # keep everyone reproducing (mutation churn)
    s.eff_activity_update()

    alive = s.get_alive().astype(bool).ravel()
    eff = s.get_eff_hashes().ravel()
    luts = [s.get_lut(i).tobytes() for i in range(N * N)]

    groups = defaultdict(set)
    for i in range(N * N):
        if alive[i]:
            groups[luts[i]].add(int(eff[i]))

    # eff hash is constant within each exact-genome group ⇒ eff = f(genome)
    assert all(len(v) == 1 for v in groups.values())
    # therefore distinct eff-species <= distinct genomes
    assert s.eff_count_distinct() <= len(groups)


def test_per_cell_used_mask_is_subset_of_union(make_sim):
    """Sanity on the exposed diagnostic substrate: every per-cell
    used-mask bit is also set in the population-union mask, and the union
    is non-empty after a run."""
    N = 24
    s = _gol_world(make_sim, N)
    rng = np.random.default_rng(7)
    s.set_v(rng.integers(0, 2, (N, N), np.uint8))
    for _ in range(20):
        s.step()
    per_cell = s.get_eff_used_mask()          # (N, N, LUT_BYTES)
    union = s.get_eff_union_mask()            # (LUT_BYTES,)
    assert union.any()
    or_of_cells = np.bitwise_or.reduce(per_cell.reshape(-1, per_cell.shape[-1]), axis=0)
    # union must cover every per-cell bit (union == OR of per-cell masks
    # for alive-ever cells; dead cells only ever add via their alive life)
    assert np.array_equal(np.bitwise_or(or_of_cells, union), union)
