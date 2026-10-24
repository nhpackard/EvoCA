"""ChemCA reproduction gate 'mass' (Kaneko style; Docs/plan-chemca.md).

Body compounds are those with transport rate u_c = 0. With
repro_gate='mass' a cell divides when its body mass reaches mass_repro (a
newborn holds about mass_repro / 2, so division is at a doubling of birth
mass) and evicts the neighbor with the lowest body mass.
"""
import numpy as np
import pytest

from chemca_helpers import isolated, cell
from chemca_py import ChemCA, gene, build_chemca
import test_chemca_identities as ident

MASS = dict(repro_gate='mass', mass_repro=0.4)   # gs3 with u_A = 0: body = A


def _grower(gate):
    # One cell, body A = 0.2 (half the threshold); F and P are transported
    # (not body), so food F comes in from the site to make more A.
    return isolated(genes=[gene('F', 'A', 'A', 2, 31)], stores=dict(F=0.3, A=0.2),
                    E=dict(F=0.3), u_F=0.5, u_P=0.5, e=0.0,
                    repro_gate=gate, mass_repro=0.4)


def test_division_fires_exactly_at_doubling():
    sim = _grower('mass')
    assert sim.body_compounds == ['A']
    prev = cell(sim.get_store('A'))
    for t in range(200):
        sim.step()
        if sim.events_last()['births']:
            break
        prev = cell(sim.get_store('A'))
        assert prev < 0.4                        # not yet doubled: no division
    else:
        pytest.fail("never divided")
    a_parent = cell(sim.get_store('A'))
    assert 2 * a_parent >= 0.4                   # body at division >= 2 x 0.2
    assert prev < 0.4                            # ... and it was below one step earlier
    assert sim.population() == 2
    assert cell(sim.get_e()) < 1.0               # energy did not gate it
    sim.free()


def test_energy_gate_does_not_divide_on_mass():
    sim = _grower('energy')
    for _ in range(200):
        sim.step()
        if cell(sim.get_store('A')) > 0.45:
            break
    assert cell(sim.get_store('A')) > 0.45 and sim.population() == 1
    assert cell(sim.get_e()) < 1.0
    sim.free()


def test_eviction_targets_lowest_body_mass():
    N = 5
    cells = [(r, c) for r in (1, 2, 3) for c in (1, 2, 3)]
    A = np.zeros((N, N), np.float32)
    e = np.zeros((N, N), np.float32)
    for k, (r, c) in enumerate(cells):
        A[r, c] = 0.15 + 0.01 * k
        e[r, c] = 0.9 - 0.05 * k             # lowest e at (3, 3)
    A[1, 2] = 0.05                           # lowest body mass
    A[2, 2] = 0.5                            # the parent (>= 0.4)
    sim = isolated(cells=cells, stores=dict(A=A), e=e, genes=[gene('F', 'A', 'A', 2, 0)],
                   u_F=0.5, u_P=0.5, **MASS)
    sim.step()
    ev = sim.events_last()
    assert ev['births'] == 1 and ev['evictions'] == 1
    assert cell(sim.get_E('A'), (1, 2)) == pytest.approx(0.05)   # evictee's body returned
    assert cell(sim.get_store('A'), (1, 2)) == pytest.approx(0.25)
    sim.free()


def test_mass_gate_needs_a_body_compound():
    with pytest.raises(ValueError, match='body compound'):
        ChemCA().init(8, preset='gs3', repro_gate='mass')
    with pytest.raises(ValueError, match='repro_gate'):
        ChemCA().init(8, preset='gs3', repro_gate='size')


def test_mass_mode_identities():
    p = dict(model='chemca', preset='gs3', u_A=0.0, tox_P=0.15,
             mu_rgene=0.05, mu_ractive=0.05, **MASS)
    seen = ident._check_run(p)
    assert all(v > 0 for v in seen.values()), seen


def test_reachable_through_init_build_and_recipe(tmp_path):
    import evoca_py
    sim = build_chemca(dict(model='chemca', preset='gs3', u_A=0.0, **MASS), N=16, seed=2)
    assert sim._lib.chemca_get_repro_gate() == 1
    assert sim._lib.chemca_get_mass_repro() == pytest.approx(0.4)
    path = sim.export_recipe('gate', runs_dir=str(tmp_path))
    sim.free()
    sim2, _ = evoca_py.import_run(path)
    assert sim2._params['repro_gate'] == 'mass'
    assert sim2._lib.chemca_get_repro_gate() == 1
    sim2.update(repro_gate='energy')
    assert sim2._lib.chemca_get_repro_gate() == 0
    sim2.free()
