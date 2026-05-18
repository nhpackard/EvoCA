"""Ring-dependence tax (research-directions §8).

`tax_ring` adds, per alive cell per step, an extra private-food
decrement of `tax_ring * level`, where `level in {1,2,3}` is the
minimum ring set the cell's LUT actually conditions on — the exact
quantity the `lut_complexity` probe computes.

Two properties are pinned here:

1. Differential: under the same `tax_ring > 0`, a population whose
   LUTs are level-1 (condition on the n1 ring only) is taxed strictly
   less per step than a population whose LUTs are level-3 (condition
   on the n3 ring), so its private food decays slower and it survives
   strictly longer.

2. No-op: `tax_ring = 0.0` reproduces the previous dynamics
   bit-for-bit. With the S1 deterministic-PRNG fix (`set_seed`), a
   short stochastic run with `tax_ring=0` must be byte-identical to
   the same run built without the parameter ever being touched.
"""
import numpy as np

from evoca_py import (LUT_BITS, lut_bit_index, make_gol_lut, pack_lut)


def _level1_lut():
    """LUT that conditions on the n1 ring only (and v_x): the output
    varies with n1 but is constant across n2 and n3. lut_complexity
    must classify this as level 1."""
    bits = np.zeros(LUT_BITS, dtype=np.uint8)
    for v_x in range(2):
        for n1 in range(5):
            val = n1 & 1                       # depends on n1 only
            for n2 in range(5):
                for n3 in range(5):
                    bits[lut_bit_index(v_x, n1, n2, n3)] = val
    return pack_lut(bits)


def _level3_lut():
    """LUT that conditions on the n3 ring (output = n3 parity), which
    forces lut_complexity to return 3."""
    bits = np.zeros(LUT_BITS, dtype=np.uint8)
    for v_x in range(2):
        for n1 in range(5):
            for n2 in range(5):
                for n3 in range(5):
                    bits[lut_bit_index(v_x, n1, n2, n3)] = n3 & 1
    return pack_lut(bits)


def _run_world(make_sim, N, lut, tax_ring, f0, n_steps):
    """Build a fresh world (all cells alive, given LUT, no eating / no
    baseline tax / no mutation, so the only private-food change is the
    ring tax), run it n_steps, and return (decrement_after_1_step,
    extinction_step or None).

    The C core has a single global state, so only one EvoCA instance
    may be active at a time — each world is built and consumed fully
    before the next is created (same pattern as test_determinism)."""
    s = make_sim(N, mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0,
                 tax=0.0, food_inc=0.0, m_scale=0.0, tax_ring=tax_ring)
    s.set_lut_all(lut)
    s.set_alive_all()
    s.set_v(np.zeros((N, N), dtype=np.uint8))
    s.set_f_all(f0)

    s.step()
    dec1 = f0 - float(s.get_f()[0, 0])
    assert np.all(s.get_alive() == 1), \
        "no cell should die in one step at f0=1.0"

    extinct_step = None
    for k in range(2, n_steps + 1):
        s.step()
        if s.get_alive().sum() == 0:
            extinct_step = k
            break
    return dec1, extinct_step


def test_level1_taxed_less_than_level3(make_sim):
    """Same tax_ring > 0: a level-1 LUT world loses strictly less
    private food per step than a level-3 LUT world, and so survives
    strictly longer. Each world is run to completion before the next
    is created (single global C state)."""
    N = 16
    tax_ring, f0 = 0.02, 1.0

    dec1, ext1 = _run_world(make_sim, N, _level1_lut(),
                            tax_ring, f0, n_steps=200)
    dec3, ext3 = _run_world(make_sim, N, _level3_lut(),
                            tax_ring, f0, n_steps=200)

    # Per-cell decrement is exactly tax_ring * level.
    assert abs(dec1 - tax_ring * 1) < 1e-6, \
        f"level-1 must lose exactly tax_ring*1; got {dec1}"
    assert abs(dec3 - tax_ring * 3) < 1e-6, \
        f"level-3 must lose exactly tax_ring*3; got {dec3}"
    assert dec1 < dec3, \
        "level-1 must be taxed strictly less per step than level-3"

    assert ext3 is not None, "level-3 world never went extinct"
    assert ext1 is not None, "level-1 world never went extinct"
    assert ext3 < ext1, \
        "level-3 (taxed 3x) must go extinct strictly before level-1"


def test_tax_ring_zero_is_bitwise_noop(make_sim):
    """tax_ring=0.0 reproduces prior dynamics exactly. With the S1
    deterministic PRNG reseed, a short stochastic run with the
    parameter at its default must be byte-identical to a run that
    never sets it.

    The config is chosen so the snapshot carries non-trivial
    deterministic state that a tax-path regression would perturb:
    a random GoL v-field evolves (v-sum changes every step) while the
    baseline tax depletes the population toward a staggered extinction,
    so the Phase 2c tax branch (the code path tax_ring lives in) fires
    on every step. Snapshot is taken mid-trajectory (step 40), with the
    population still alive and v actively evolving, so all four arrays
    are non-uniform fingerprints rather than a frozen lattice."""
    N = 24

    def _run(set_tax_ring):
        kw = dict(mu_lut=0.05, mu_egene=0.0, mu_egenome=0.0,
                  tax=0.01, food_inc=0.05, m_scale=2.0)
        if set_tax_ring:
            kw['tax_ring'] = 0.0
        s = make_sim(N, **kw)
        s.set_seed(20260516)
        s.set_lut_all(make_gol_lut())
        s.set_alive_all()
        rng = np.random.default_rng(0)
        s.set_v((rng.random((N, N)) < 0.4).astype(np.uint8))
        s.set_f_all(0.5)
        for _ in range(40):
            s.step()
        a, v, f = s.get_alive().copy(), s.get_v().copy(), s.get_f().copy()
        # Guard: the fingerprint must actually be non-trivial, else the
        # equality below would pass even for a broken tax path.
        assert a.sum() == N * N, "snapshot must be taken while alive"
        assert 0 < int(v.sum()) < N * N, "v-field must be evolving"
        return (a, v, f, s.get_lut(0).copy())

    a0, v0, f0, l0 = _run(set_tax_ring=False)   # parameter never touched
    a1, v1, f1, l1 = _run(set_tax_ring=True)    # explicitly tax_ring=0.0

    assert np.array_equal(a0, a1), "alive grid changed with tax_ring=0"
    assert np.array_equal(v0, v1), "v grid changed with tax_ring=0"
    assert np.array_equal(f0, f1), "private food changed with tax_ring=0"
    assert np.array_equal(l0, l1), "LUT changed with tax_ring=0"
