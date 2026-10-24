"""ChemCA mutation reachability: requested through the public doors
(init(..., mu_rgene, mu_ractive) and run_chemca), asserted on the genes."""
import numpy as np
import pytest

from chemca_py import ChemCA, C_NONE, G_A, G_B, G_C, G_M, G_K, NK, PRESETS, run_chemca


# Small-lattice settings that keep a reproducing population at N=32.
SMALL = {'fw2': dict(m_e=0.02, light=0.03), 'gs3': {}}


def _evolve(preset, mu_rgene, mu_ractive, steps=400):
    sim = ChemCA()
    sim.init(32, preset=preset, mu_rgene=mu_rgene, mu_ractive=mu_ractive, **SMALL[preset])
    sim.state(alive='fraction', seed=5)
    g0, a0 = sim.get_rgenes()
    r0, c0 = np.argwhere(sim.get_alive() == 1)[0]
    mutant_births = 0
    for _ in range(steps):
        sim.step()
        mutant_births += int((sim.get_births() == 2).sum())
    g, a = sim.get_rgenes()
    al = sim.get_alive() == 1
    out = dict(g0=g0[r0, c0].copy(), g=g[al], a=a[al], mutant_births=mutant_births,
               a0=int(a0[r0, c0]), pop=int(al.sum()))
    sim.free()
    return out


def _valid(g):
    C, m = g[..., G_C], g[..., G_M]
    return (np.all(g[..., G_A] != g[..., G_B])
            and np.all((C == C_NONE) == (m == 0))
            and np.all(g[..., G_K] < NK) and np.all(m <= 2))


@pytest.mark.parametrize('preset', ['fw2', 'gs3'])
def test_mutation_reachable_through_init(preset):
    r = _evolve(preset, 0.1, 0.1)
    assert r['pop'] > 0 and r['mutant_births'] > 0
    g = r['g']
    assert _valid(g)                                   # A = B never appears
    assert np.all(r['a'] != 0)
    assert len(np.unique(g[..., G_K])) > len(np.unique(r['g0'][:, G_K]))   # k_idx moved
    cm0 = {(x[G_C], x[G_M]) for x in r['g0']}
    cm = {(x[G_C], x[G_M]) for x in g.reshape(-1, 5)}
    assert cm - cm0                                    # C / m toggled
    ab0 = {(x[G_A], x[G_B]) for x in r['g0']}
    ab = {(x[G_A], x[G_B]) for x in g.reshape(-1, 5)}
    if preset == 'gs3':
        assert ab - ab0                                # A or B redrawn
    assert len(np.unique(r['a'])) > 1                  # active bits flipped


def test_no_mutation_at_zero_rates():
    r = _evolve('gs3', 0.0, 0.0)
    assert r['mutant_births'] == 0
    assert np.all(r['g'] == r['g0'][None]) and np.all(r['a'] == r['a0'])


def test_mutation_reachable_through_run_chemca():
    p = dict(model='chemca', preset='fw2', mu_rgene=0.1, mu_ractive=0.1, **SMALL['fw2'])
    r = run_chemca(p, n_steps=400, N=32, seed=5, sample_every=100, init='fraction')
    r0 = run_chemca(dict(p, mu_rgene=0.0, mu_ractive=0.0), n_steps=400, N=32, seed=5, init='fraction',
                    sample_every=100)
    assert r0['n_genesets_final'] == 1 and list(r0['k_idx_by_type']) == ['F->W']
    assert r['n_genesets_final'] > 1
    k0 = PRESETS['fw2']['genes'][0]['k_idx']
    assert r0['k_idx_by_type']['F->W'][1] == k0
    assert len(r['k_idx_by_type']) > 1 or r['k_idx_by_type']['F->W'][1] != k0
    assert r['params']['mu_rgene'] == 0.1
