"""ChemCA mass and energy identities, every step, on evolving runs.

Closed mode: total mass sum_c (sum E_c + sum store_c) is constant to float
rounding (the light reaction and the stores are float32; diffusion uses the
unbiased core_diffuse_precise). Energy: d(sum e) = yield - cost - death -
split, and the compound free energy Phi = sum_c G_c M_c changes by
light - yield (+ open_in - open_out). The runs are chosen so that every
death path (toxic, energy, eviction) and reproduction occur, and the
assertions say so.
"""
import numpy as np
import pytest

from chemca_py import build_chemca

MASS_TOL = 1e-6      # relative, per step
ENERGY_TOL = 1e-9    # absolute, per step (C books the realized float changes of e)

RUNS = {
    'fw2': dict(model='chemca', preset='fw2', mu_rgene=0.05, mu_ractive=0.05,
                m_e=0.018, light=0.03),
    'gs3': dict(model='chemca', preset='gs3', mu_rgene=0.05, mu_ractive=0.05),
}


def _check_run(params, steps=400, N=32, open_mode=False):
    sim = build_chemca(params, N=N, seed=3, init='fraction')
    M0 = sim.mass()['total']
    e0, P0 = sim.e_total(), sim.free_energy()
    seen = dict(births=0, tox_deaths=0, energy_deaths=0, evictions=0)
    for _ in range(steps):
        sim.step()
        fl = sim.mass_flux()
        M = sim.mass()['total']
        assert abs(M - M0 - (fl[:, 0].sum() - fl[:, 1].sum())) / M0 < MASS_TOL
        b = sim.energy_budget(cum=False)
        e1, P1 = sim.e_total(), sim.free_energy()
        assert abs((e1 - e0) - (b['yield'] - b['cost'] - b['death'] - b['split'])) < ENERGY_TOL
        assert abs((P1 - P0) - (b['light'] - b['yield'] + b['open_in'] - b['open_out'])) \
            < MASS_TOL * M0 + ENERGY_TOL
        e0, P0 = e1, P1
        for k in seen:
            seen[k] += sim.events_last()[k]
        if sim.population() == 0:
            break
    sim.free()
    return seen


@pytest.mark.parametrize('preset', sorted(RUNS))
def test_closed_mass_and_energy_identities(preset):
    seen = _check_run(RUNS[preset])
    assert all(v > 0 for v in seen.values()), seen     # every path exercised


def test_open_mode_mass_balance():
    p = dict(RUNS['gs3'], open_mode=1, inc_F=0.002, sink_P=0.01, sink_A=0.001)
    _check_run(p, steps=200)


def test_closed_mass_conserved_over_long_run():
    # Drift check: 3000 steps at N=64, total relative error stays < 1e-5.
    sim = build_chemca(dict(model='chemca', preset='gs3'), N=64, seed=1)
    M0 = sim.mass()['total']
    sim.step(3000)
    assert abs(sim.mass()['total'] - M0) / M0 < 1e-5
    sim.free()
