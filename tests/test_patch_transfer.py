"""Patch-transfer API (Docs/EvoCA_research_directions.md §11).

Three properties of the genome-level transplant API:

  1. extract -> stamp is a faithful genome round-trip: every alive
     cell's heritable genome (LUT bytes, active egene value/mask slots,
     active mask) is reproduced exactly, so the species hash is
     preserved. Dead-cell scratch bytes are NOT heritable and are
     deliberately not transferred — equality is asserted on alive cells.
  2. Self-transfer is approximately neutral: a population stamped back
     into a field initialised from its own recipe neither collapses nor
     explodes (the §11 null control behaves like a null).
  3. mode='invade' only writes cells that land on DEAD host cells; it
     never displaces a resident.

Kept fast: small N, few ticks.
"""
import numpy as np

from patch_transfer import run_assay, self_transfer


def _seed(s):
    """A heterogeneous, reproducible donor population."""
    s.state(lut='random', lut_n_init=3, lut_density=0.5, lut_seed=12345,
            egenome='random', alive='fraction', alive_fraction=0.6)


# ── 1. round-trip ─────────────────────────────────────────────────────

def test_extract_stamp_roundtrip_genomes_exact(make_sim):
    N, side, r0, c0 = 32, 12, 8, 8
    src = make_sim(N, food_inc=0.02, m_scale=1.2)
    _seed(src)
    patch = src.extract_patch(r0, c0, side)

    # Stamp into a deliberately DIFFERENT host so any leakage shows up.
    dst = make_sim(N, food_inc=0.02, m_scale=1.2)
    dst.state(lut='gol', egenome='uniform', alive='all')
    n = dst.stamp_patch(patch, r0, c0, mode='overwrite')
    assert n == side * side

    back = dst.extract_patch(r0, c0, side)
    am = patch.alive.astype(bool)

    assert np.array_equal(patch.alive, back.alive), \
        "alive mask not reproduced"
    assert np.array_equal(patch.lut[am], back.lut[am]), \
        "LUT bytes of alive cells changed across round-trip"
    assert np.array_equal(patch.active[am], back.active[am]), \
        "active mask of alive cells changed across round-trip"
    assert np.array_equal(patch.egenes[am], back.egenes[am]), \
        "egene value bytes of alive cells changed across round-trip"
    assert np.array_equal(patch.egenes_mask[am], back.egenes_mask[am]), \
        "egene mask bytes of alive cells changed across round-trip"


def test_roundtrip_preserves_species_hash(make_sim):
    """Identical alive genomes => identical species count. Stamping the
    patch onto itself must not change how many distinct genomes the
    activity machinery sees in that region."""
    N, side, r0, c0 = 24, 9, 6, 6
    src = make_sim(N, food_inc=0.02, m_scale=1.2)
    _seed(src)
    src._lib.evoca_activity_update()
    before = src.n_distinct_genomes()

    patch = src.extract_patch(r0, c0, side)
    src.stamp_patch(patch, r0, c0, mode='overwrite')
    src._lib.evoca_activity_update()
    after = src.n_distinct_genomes()
    assert before == after, \
        f"distinct genome count drifted on self-stamp: {before} -> {after}"


# ── 2. self-transfer ≈ neutral ────────────────────────────────────────

def test_self_transfer_approximately_neutral(make_sim):
    """The §11 null: a population transplanted into its own field must
    not systematically collapse. We require the patch alive-fraction to
    stay within a loose band of t=0 over a short run — i.e. it does not
    crash to extinction (the failure mode the control exists to rule
    out)."""
    N, side, ticks = 32, 10, 25
    recipe = dict(lut='random', lut_n_init=2, lut_density=0.5,
                  lut_seed=7, egenome='uniform', alive='fraction',
                  alive_fraction=0.5)
    meta = dict(food_inc=0.05, m_scale=1.2, tax=0.0)
    res = self_transfer(recipe, N, side, ticks, meta, mode='overwrite')

    f0 = res['alive_frac_t0']
    ff = res['alive_frac_final']
    assert f0 > 0.1, "degenerate seed — control is uninformative"
    # Not a collapse: final fraction is a meaningful share of t=0.
    assert ff > 0.4 * f0, \
        f"self-transfer collapsed: {f0:.3f} -> {ff:.3f} (not neutral)"
    # Not a blow-up beyond the lattice either.
    assert ff <= 1.0 + 1e-9


# ── 3. invade only writes dead cells ──────────────────────────────────

def test_invade_only_writes_dead_cells(make_sim):
    N, side, r0, c0 = 28, 14, 7, 7

    # Donor: a fully-alive patch with a distinctive uniform genome so
    # we can detect any cell it overwrote.
    donor = make_sim(N, food_inc=0.0, m_scale=0.0)
    donor.state(lut='random', lut_n_init=3, lut_density=0.7,
                lut_seed=99, egenome='uniform', alive='all')
    patch = donor.extract_patch(r0, c0, side)
    assert patch.alive.all(), "donor patch should be fully alive"

    # Host: a known checkerboard of alive/dead in the target region.
    host = make_sim(N, food_inc=0.0, m_scale=0.0)
    host.state(lut='gol', egenome='uniform', alive='all')
    a = host.get_alive()
    rr, cc = np.mgrid[0:N, 0:N]
    a[(rr + cc) % 2 == 0] = 0          # half the cells dead
    host.set_alive(a)

    host_alive_before = host.get_alive().astype(bool)
    host_lut_before = host.extract_patch(r0, c0, side).lut.copy()

    n = host.stamp_patch(patch, r0, c0, mode='invade')

    host_alive_after = host.get_alive().astype(bool)
    rs, cs = slice(r0, r0 + side), slice(c0, c0 + side)
    region_before = host_alive_before[rs, cs]
    region_after = host_alive_after[rs, cs]
    back = host.extract_patch(r0, c0, side)

    # Every host cell that was ALIVE before is byte-identical after:
    # invade must never displace a resident.
    occ = region_before
    assert np.array_equal(back.lut[occ], host_lut_before[occ]), \
        "invade overwrote a living resident's genome"

    # Every host cell that was DEAD before is now alive and carries the
    # donor genome (the donor patch was fully alive everywhere).
    empty = ~region_before
    assert region_after[empty].all(), \
        "invade failed to settle into available dead cells"
    assert np.array_equal(back.lut[empty], patch.lut[empty]), \
        "settled cells did not receive the donor genome"

    # Outside the original-dead set nothing new was born.
    assert np.array_equal(region_after[occ], region_before[occ])
    assert n == int(empty.sum())

    # And outside the patch the host is completely untouched.
    outside = host_alive_before.copy()
    outside[rs, cs] = host_alive_after[rs, cs]
    assert np.array_equal(outside, host_alive_after), \
        "stamp_patch touched cells outside the target region"


# ── harness smoke ─────────────────────────────────────────────────────

def test_run_assay_returns_well_formed_dict(make_sim):
    N, side, ticks = 24, 8, 12
    donor = dict(lut='random', lut_n_init=2, lut_seed=1,
                 egenome='uniform', alive='fraction', alive_fraction=0.6)
    host = dict(lut='gol', egenome='uniform', alive='fraction',
                alive_fraction=0.3)
    meta = dict(food_inc=0.05, m_scale=1.2)
    r = run_assay(donor, host, N, side, ticks, meta, mode='overwrite')

    assert len(r['alive_frac']) == ticks + 1
    assert len(r['boundary_flux']) == ticks
    assert r['alive_frac_rel'][0] == 0.0
    assert isinstance(r['success'], bool)
    assert r['perimeter_cells'] > 0
