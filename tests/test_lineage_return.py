"""Per-lineage net metabolic return reducer (checkpoint §D2/§B/§H).

The reducer sits on top of the opt-in lineage field. It attributes, per
genome-hash lineage, the cumulative mouthful its members eat and the
cumulative tax they pay, so net return = intake − tax is measurable per
lineage over the whole run (extinct lineages retained). These tests pin:

1. Default OFF: the reducer is empty, the totals are zero, the per-cell
   accumulators are all-zero — genuinely opt-in.
2. Zero perturbation: with lineage ON vs OFF (identical seed/setup) the CA
   state and alive mask are bit-identical after a run — the instrument does
   not touch dynamics.
3. Exact tax accounting: in a no-eating / no-birth / no-death world every
   alive cell pays `tax` each step, so the reducer's cumulative tax equals
   K·tax·N_alive exactly and the per-cell accumulators match.
4. Completeness: in an evolving world (eating + tax + mutation) the reducer
   column sums equal the independent run-cumulative intake/tax totals, and
   net == intake − tax.
5. Birth attribution + per-cell reset: a single controlled reproduction
   registers one birth in the reducer and the child's per-cell accumulators
   start at zero.
"""
import numpy as np


def _world(make_sim, N, lineage=False, **kw):
    from evoca_py import make_gol_lut
    s = make_sim(N, lineage=lineage, **kw)
    s.set_lut_all(make_gol_lut())
    s.set_egenome_all(0)           # value 0, mask 0x3F (all orbits set)
    s.set_alive_all()
    return s


def test_default_off_reducer_is_empty_and_zero(make_sim):
    """Lineage OFF: no reducer rows, zero totals, zero per-cell arrays."""
    s = _world(make_sim, 16, lineage=False,
               food_inc=0.1, m_scale=0.5, tax=0.02)
    s.set_v(np.zeros((16, 16), dtype=np.uint8))
    s.set_F_all(1.0)
    s.set_f_all(0.5)
    for _ in range(4):
        s.step()
    red = s.get_lineage_return()
    assert red.size == 0, "reducer should be empty with lineage OFF"
    assert s.lineage_return_totals() == (0.0, 0.0)
    assert not s.get_lineage_intake().any(), "per-cell intake nonzero (OFF)"
    assert not s.get_lineage_tax().any(), "per-cell tax nonzero (OFF)"


def test_lineage_return_does_not_perturb_dynamics(make_sim):
    """Enabling the reducer must not change the CA. Identical seed +
    setup with lineage ON vs OFF -> identical v and alive after a run."""
    N, K = 20, 12

    def run(lineage):
        s = _world(make_sim, N, lineage=lineage,
                   food_inc=0.1, m_scale=0.5, tax=0.02,
                   mu_lut=0.002, mu_egenome=0.02)
        s.set_seed(4242)
        s.set_v(np.zeros((N, N), dtype=np.uint8))
        s.set_F_all(1.0)
        s.set_f_all(0.5)
        for _ in range(K):
            s.step()
        return s.get_v().copy(), s.get_alive().copy()

    v_off, a_off = run(False)
    v_on,  a_on  = run(True)
    assert np.array_equal(v_off, v_on), "lineage ON perturbed CA state"
    assert np.array_equal(a_off, a_on), "lineage ON perturbed alive mask"


def test_reducer_tax_accounting_exact(make_sim):
    """No eating (m_scale=0), f in (0,1) so nobody dies or reproduces:
    every alive cell pays `tax` each step. The reducer's cumulative tax
    must equal K·tax·N_alive exactly, intake must be zero, and the single
    live genome must own all of it."""
    N, K = 8, 5
    tax = 0.01
    s = _world(make_sim, N, lineage=True,
               food_inc=0.0, m_scale=0.0, tax=tax)
    s.set_v(np.zeros((N, N), dtype=np.uint8))
    s.set_f_all(0.9)               # stays in (0,1): no death, no repro
    for _ in range(K):
        s.step()
        s.lineage_return_update()

    assert int(s.get_alive().sum()) == N * N, "a cell died unexpectedly"
    intake_total, tax_total = s.lineage_return_totals()
    assert intake_total == 0.0, "no eating expected"

    # C accumulates (double)(float32 tax) per alive cell per step.
    expect = K * N * N * float(np.float32(tax))
    assert abs(tax_total - expect) < 1e-9 * expect, \
        f"tax_total {tax_total} != {expect}"

    red = s.get_lineage_return(update=False)
    assert red.size == 1, "mutation off should leave one genome lineage"
    row = red[0]
    assert abs(row['tax'] - tax_total) < 1e-9 * expect
    assert row['intake'] == 0.0
    assert row['net'] == row['intake'] - row['tax']
    assert row['net'] < 0.0, "tax-only lineage should have negative return"
    assert row['pop'] == N * N
    assert row['births'] == 0, "no reproduction, so no births"
    # per-cell accumulator: each alive cell paid K·tax over its life.
    pc_tax = s.get_lineage_tax()
    assert np.allclose(pc_tax, K * float(np.float32(tax)))
    assert not s.get_lineage_intake().any()


def test_reducer_accounting_completeness(make_sim):
    """Evolving world (eat + tax + mutation): the reducer's per-lineage
    columns must sum to the independent run-cumulative totals, and
    net == intake − tax for every lineage."""
    N, K = 24, 20
    s = _world(make_sim, N, lineage=True,
               food_inc=0.1, m_scale=0.5, tax=0.02,
               mu_lut=0.002, mu_egenome=0.02)
    s.set_seed(1234)
    s.set_v(np.zeros((N, N), dtype=np.uint8))
    s.set_F_all(1.0)
    s.set_f_all(0.5)
    for _ in range(K):
        s.step()
        s.lineage_return_update()

    red = s.get_lineage_return(update=False)
    intake_total, tax_total = s.lineage_return_totals()
    assert red.size >= 1
    assert intake_total > 0.0, "cells should have eaten"
    assert tax_total > 0.0, "cells should have paid tax"
    assert abs(red['intake'].sum() - intake_total) < 1e-6 * intake_total, \
        "per-lineage intake does not sum to the run total"
    assert abs(red['tax'].sum() - tax_total) < 1e-6 * tax_total, \
        "per-lineage tax does not sum to the run total"
    assert np.allclose(red['net'], red['intake'] - red['tax'])
    # Lifecycle fields are sane: appeared and last-seen within the run.
    assert (red['first_step'] >= 1).all()
    assert (red['first_step'] <= K).all()
    assert (red['last_step'] <= K).all()
    # Some births must have occurred (f started at 0.5 and cells ate).
    assert red['births'].sum() >= 1


def test_birth_registers_and_resets_child(make_sim):
    """One founder in a dead field reproduces once: the reducer records at
    least one birth and the child's per-cell accumulators start at zero
    (it is born in Phase 4, after that step's tax/eating)."""
    from evoca_py import make_gol_lut
    N = 16
    s = make_sim(N, food_inc=0.0, m_scale=0.0, tax=0.0, lineage=True)
    s.set_lut_all(make_gol_lut())
    s.set_egenome_all(0)
    alive = np.zeros((N, N), dtype=np.uint8)
    alive[8, 8] = 1
    s.set_alive(alive)
    s.set_f_all(2.0)               # only the alive founder can reproduce
    s.step()
    s.lineage_return_update()

    bid = s.get_lineage_birth_id()
    born = np.argwhere(bid > 0)
    assert born.shape[0] == 1, "expected exactly one child"
    r, c = born[0]
    intake = s.get_lineage_intake()
    tax = s.get_lineage_tax()
    assert intake[r, c] == 0.0, "freshly born child should have 0 intake"
    assert tax[r, c] == 0.0, "freshly born child should have 0 tax"

    red = s.get_lineage_return(update=False)
    assert red['births'].sum() >= 1, "reducer did not register the birth"
