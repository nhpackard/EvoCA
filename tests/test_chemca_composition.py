"""composition_clusters: the isologous-diversification readout, on
constructed populations (homogeneous -> 1, planted compositions -> 2/3)."""
import numpy as np
import pytest

from chemca_py import ChemCA, composition_clusters, run_chemca


def _pop(blocks, noise=0.01, N=32, seed=0):
    """All cells alive; columns split into len(blocks) bands, each with the
    composition blocks[j] = (F, A, P) times (1 + noise * normal)."""
    sim = ChemCA()
    sim.init(N, preset='gs3')
    sim.state(alive='all', seed=seed)
    rng = np.random.default_rng(seed)
    band = np.arange(N) * len(blocks) // N
    for ci, c in enumerate('FAP'):
        v = np.array([b[ci] for b in blocks])[band][None, :].repeat(N, 0)
        sim.set_store(c, v * (1 + noise * rng.standard_normal((N, N))))
    return sim


def test_homogeneous_is_one_cluster():
    sim = _pop([(0.4, 0.1, 0.1)])
    r = composition_clusters(sim)
    assert r['n_clusters'] == 1 and r['fractions'] == [1.0]
    sim.free()


def test_exactly_uniform_is_one_cluster():
    sim = _pop([(0.4, 0.1, 0.1)], noise=0.0)
    assert composition_clusters(sim)['n_clusters'] == 1
    sim.free()


@pytest.mark.parametrize('blocks', [
    [(0.4, 0.1, 0.1), (0.1, 0.3, 0.2)],
    [(0.4, 0.1, 0.1), (0.1, 0.3, 0.2), (0.05, 0.05, 0.5)],
])
def test_planted_compositions_are_found(blocks):
    sim = _pop(blocks)
    r = composition_clusters(sim)
    assert r['n_clusters'] == len(blocks)
    assert r['silhouette'] > 0.8
    assert np.allclose(sorted(r['fractions']), sorted([1 / len(blocks)] * len(blocks)), atol=0.04)
    sim.free()


def test_scale_invariant():
    # Same composition, different totals: one type (stores are normalized).
    sim = _pop([(0.4, 0.1, 0.1), (0.8, 0.2, 0.2)])
    assert composition_clusters(sim)['n_clusters'] == 1
    sim.free()


def test_run_chemca_reports_composition():
    r = run_chemca(dict(model='chemca', preset='gs3'), n_steps=300, N=32, seed=1,
                   kymo_steps=100, comp_window=50)
    for k in ('n_comp_clusters', 'comp_silhouette', 'comp_fractions', 'comp_n_cells'):
        assert k in r
    assert r['comp_n_cells'] > 0 and abs(sum(r['comp_fractions']) - 1) < 1e-9
