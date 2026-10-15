"""Tests for the mu-genome (A2) Wave-1 storage scaffold: the per-cell
LUT mutation-rate gene `mu_lut_gene`.

In this scope the gene is seeded to the run's mu_lut, inherited exactly
on reproduction, never mutated, and read nowhere that affects dynamics.
"""
import numpy as np

N = 48


def _sim(make_sim, mu_lut=0.01):
    s = make_sim(N, food_inc=0.013, m_scale=1.2, gdiff=0.06, tax=0.035,
                 mu_lut=mu_lut, mu_egene=0.003, mu_egenome=0.005)
    s.set_seed(3)
    rng = np.random.default_rng(3)
    s.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    s.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.5)
    s.set_F_all(1.0)
    s.set_alive_halfplane(0)
    return s


def test_gene_seeded_to_mu_lut(make_sim):
    s = _sim(make_sim, mu_lut=0.01)
    g = s.get_mu_lut_gene()
    assert g.shape == (N, N)
    assert np.all(g == np.float32(0.01))


def test_seeded_gene_is_dynamically_inert(make_sim):
    """A wild gene field must leave the whole state bit-for-bit unchanged."""
    a = _sim(make_sim)
    b = _sim(make_sim)
    b.set_mu_lut_gene(np.random.default_rng(9).uniform(0, 0.5, (N, N)))
    for _ in range(200):
        a.step(); b.step()
    for get in ('get_v', 'get_alive', 'get_f', 'get_F', 'get_lut_hashes'):
        assert np.array_equal(getattr(a, get)(), getattr(b, get)()), get


def test_gene_inherited_not_mutated(make_sim):
    s = _sim(make_sim)
    band = np.zeros((N, N), np.float32)
    band[:, :N // 3] = 0.02
    band[:, N // 3:2 * N // 3] = 0.09
    band[:, 2 * N // 3:] = 0.15
    s.set_mu_lut_gene(band)
    births = 0
    for _ in range(200):
        s.step()
        births += int((s.get_births() > 0).sum())
    g = s.get_mu_lut_gene()
    alive = s.get_alive().reshape(N, N) == 1
    assert births > 0
    assert int((g != band).sum()) > 0               # inheritance moved genes
    assert set(np.unique(g[alive])) <= {np.float32(0.02), np.float32(0.09),
                                        np.float32(0.15)}
