"""ChemCA constructed one-step fixtures (Docs/plan-chemca.md, tests).

Each test switches on one mechanism in an otherwise inert world
(chemca_helpers.isolated) and compares with a hand calculation.
"""
import math

import numpy as np
import pytest

from chemca_helpers import isolated, cell
from chemca_py import gene

RTOL = 2e-6   # float32 stores


def k_of(sim, k_idx):
    return sim._params['k_min'] * 2 ** (k_idx / 4)


@pytest.mark.parametrize('C,m', [(None, 0), ('A', 1), ('A', 2)])
def test_conversion_matches_exponential_integrator(C, m):
    # F -> A, uncatalyzed / catalyzed by A at order 1 or 2.
    kidx = 20
    sim = isolated(genes=[gene('F', 'A', C, m, kidx)], stores=dict(F=0.8, A=0.3, P=0.1))
    sim.step()
    lam = k_of(sim, kidx) * (0.3 ** m if C else 1.0)
    conv = 0.8 * (1 - math.exp(-lam))
    assert cell(sim.get_store('F')) == pytest.approx(0.8 - conv, rel=RTOL)
    assert cell(sim.get_store('A')) == pytest.approx(0.3 + conv, rel=RTOL)
    # yield G_F - G_A = 0.5 per unit, booked and credited
    assert cell(sim.get_e()) == pytest.approx(0.5 * conv, rel=RTOL)
    assert sim.energy_budget()['yield'] == pytest.approx(0.5 * conv, rel=RTOL)
    sim.free()


def test_two_genes_share_substrate_in_proportion():
    k1, k2 = 22, 14
    sim = isolated(genes=[gene('F', 'A', None, 0, k1), gene('F', 'P', None, 0, k2)],
                   active=0b11, stores=dict(F=0.9))
    sim.step()
    l1, l2 = k_of(sim, k1), k_of(sim, k2)
    conv = 0.9 * (1 - math.exp(-(l1 + l2)))
    assert cell(sim.get_store('F')) == pytest.approx(0.9 - conv, rel=RTOL)
    assert cell(sim.get_store('A')) == pytest.approx(conv * l1 / (l1 + l2), rel=RTOL)
    assert cell(sim.get_store('P')) == pytest.approx(conv * l2 / (l1 + l2), rel=RTOL)
    # slot order does not matter
    sim2 = isolated(genes=[gene('F', 'P', None, 0, k2), gene('F', 'A', None, 0, k1)],
                    active=0b11, stores=dict(F=0.9))
    sim2.step()
    assert cell(sim2.get_store('A')) == pytest.approx(conv * l1 / (l1 + l2), rel=RTOL)


def test_uphill_gene_stops_when_energy_runs_out():
    # P -> F costs G_F - G_P = 1 per unit; e = 0.01 pays for 0.01 units.
    sim = isolated(genes=[gene('P', 'F', None, 0, 31)], stores=dict(P=0.5), e=0.01)
    sim.step()
    assert cell(sim.get_store('F')) == pytest.approx(0.01, rel=1e-5)
    assert cell(sim.get_e()) == pytest.approx(0.0, abs=1e-9)
    f1 = cell(sim.get_store('F'))
    sim.step()                         # no energy left: nothing more runs
    assert cell(sim.get_store('F')) == pytest.approx(f1, abs=1e-9)
    assert cell(sim.get_e()) >= 0.0


@pytest.mark.parametrize('a,E', [(0.2, 1.0), (1.0, 0.2)])
def test_transport_both_directions(a, E):
    u = 0.25
    sim = isolated(stores=dict(F=a), E=dict(F=E), u_F=u, genes=[gene('A', 'P')])
    sim.step()
    d = u * (E - a)
    assert cell(sim.get_store('F')) == pytest.approx(a + d, rel=RTOL)
    assert cell(sim.get_E('F')) == pytest.approx(E - d, rel=RTOL)
    sim.free()


def test_transport_only_for_alive_cells():
    sim = isolated(stores=dict(F=0.2), E=dict(F=1.0), u_F=0.25, genes=[gene('A', 'P')])
    sim.step()
    assert sim.get_E('F')[0, 0] == pytest.approx(1.0)   # dead site untouched


def test_toxic_death_returns_stores():
    st = dict(F=0.3, A=0.2, P=0.35)
    sim = isolated(stores=st, e=0.4, tox_P=0.3, genes=[gene('A', 'P', None, 0, 0)],
                   k_min=1e-9)
    sim.step()
    assert sim.population() == 0
    assert sim.events_last()['tox_deaths'] == 1
    for c, v in st.items():
        assert cell(sim.get_E(c)) == pytest.approx(v, rel=1e-5)
    assert sim.energy_budget()['death'] == pytest.approx(0.4, rel=1e-6)


def test_energy_death_when_costs_cannot_be_paid():
    sim = isolated(stores=dict(F=0.3), e=0.05, m_e=0.1, genes=[gene('A', 'P')])
    sim.step()
    assert sim.population() == 0
    assert sim.events_last()['energy_deaths'] == 1
    assert cell(sim.get_E('F')) == pytest.approx(0.3)
    b = sim.energy_budget()
    assert b['cost'] == pytest.approx(0.05, rel=1e-6) and b['death'] == 0.0


def test_costs_per_active_gene():
    sim = isolated(e=0.5, m_e=0.01, c_g=0.02, active=0b111,
                   genes=[gene('A', 'P')] * 3)
    sim.step()
    assert cell(sim.get_e()) == pytest.approx(0.5 - 0.01 - 3 * 0.02, rel=1e-6)


def _repro_world(neighbor_e, dead=()):
    """Parent at the center (e = 1.5) with all 8 neighbors alive at the
    given energies except `dead`; no conversion (genes use F, stores
    hold only A and P)."""
    cells = [(r, c) for r in (1, 2, 3) for c in (1, 2, 3) if (r, c) not in dead]
    N = 5
    st = {c: np.zeros((N, N), np.float32) for c in 'FAP'}
    e = np.zeros((N, N), np.float32)
    for k, (r, c) in enumerate(cells):
        st['A'][r, c] = 0.1 + 0.01 * k
        st['P'][r, c] = 0.05 + 0.003 * k
    st['A'][2, 2], st['P'][2, 2] = 0.37, 0.21
    for (r, c), v in neighbor_e.items():
        e[r, c] = v
    e[2, 2] = 1.5
    sim = isolated(cells=cells, stores=st, e=e, genes=[gene('F', 'A', 'A', 2, 0)])
    return sim, st, e


def test_eviction_returns_evictee_stores_and_split_is_exact():
    ne = {(r, c): 0.2 + 0.05 * k for k, (r, c) in
          enumerate([(1, 1), (1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2), (3, 3)])}
    ne[(3, 2)] = 0.01                                    # the lowest e
    sim, st, e = _repro_world(ne)
    sim.step()
    ev = sim.events_last()
    assert ev['births'] == 1 and ev['evictions'] == 1
    # evictee's stores went to its site; its e was lost
    assert cell(sim.get_E('A'), (3, 2)) == pytest.approx(float(st['A'][3, 2]))
    assert cell(sim.get_E('P'), (3, 2)) == pytest.approx(float(st['P'][3, 2]))
    assert sim.energy_budget()['death'] == pytest.approx(0.01, rel=1e-6)
    # 50/50 split, exact in float
    for c in 'AP':
        x = np.float32(st[c][2, 2])
        assert sim.get_store(c)[3, 2] == x * np.float32(0.5)
        assert sim.get_store(c)[2, 2] == x - x * np.float32(0.5)
    assert sim.get_e()[2, 2] == np.float32(0.75) and sim.get_e()[3, 2] == np.float32(0.75)
    assert sim.energy_budget()['split'] == 0.0
    g, a = sim.get_rgenes()
    assert np.array_equal(g[3, 2], g[2, 2]) and a[3, 2] == a[2, 2]


def test_dead_neighbor_preferred_over_lowest_e():
    ne = {(r, c): 0.3 for r in (1, 2, 3) for c in (1, 2, 3)}
    ne[(1, 1)] = 0.0                                     # alive with e = 0
    sim, st, e = _repro_world(ne, dead=((3, 3),))
    sim.step()
    ev = sim.events_last()
    assert ev['births'] == 1 and ev['evictions'] == 0
    assert sim.get_alive()[3, 3] == 1 and sim.get_alive()[1, 1] == 1


def test_cycle_yield_is_zero():
    # fw2: F -> W and W -> F at equal rates. At F = W the two fluxes cancel
    # and the yield is exactly paid back by the uphill cost.
    k = 20
    sim = isolated('fw2', genes=[gene('F', 'W', None, 0, k), gene('W', 'F', None, 0, k)],
                   active=0b11, stores=dict(F=0.5, W=0.5), e=0.5)
    sim.step()
    assert cell(sim.get_e()) == pytest.approx(0.5, abs=1e-7)
    assert cell(sim.get_store('F')) == pytest.approx(0.5, abs=1e-7)
    # Away from equilibrium: e is a state function, e - e0 = -sum G (a - a0),
    # so any path that returns the stores returns e (zero yield per cycle).
    sim2 = isolated('fw2', genes=[gene('F', 'W', None, 0, k), gene('W', 'F', None, 0, k - 3)],
                    active=0b11, stores=dict(F=0.8, W=0.2), e=0.5)
    for _ in range(40):
        sim2.step()
        dF = cell(sim2.get_store('F')) - 0.8
        assert cell(sim2.get_e()) - 0.5 == pytest.approx(-1.0 * dF, abs=2e-6)
