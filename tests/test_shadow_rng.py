"""The Channon neutral shadow must be a pure observer.

Before the `shadow-rng` branch the shadow's demography and mutation
(neutral_apply_demography / neut_mutate_inplace in C/evoca.c) drew from
the SAME xorshift32 stream as the dynamics, so enabling the shadow
advanced the dynamics RNG and changed the real trajectory: same seed,
shadow on vs off, gave different lattices. The shadow now draws from its
own stream (g_shadow_rng), derived from the seed in evoca_set_seed.

All tests enter through public boundaries: the EvoCA object API and
evoca_explore.run_sim. The fixture is the viable pure-evo config from
tools/parity_probe.py; each test checks the population is alive at the
end, because an extinct fixture would pass trivially.
"""
import hashlib

import numpy as np

# Viable pure-evo config (tools/parity_probe.py), LUT mutation on so the
# shadow actually consumes RNG draws every step.
_P = dict(food_inc=0.01, m_scale=1.5, gdiff=0.08, mu_lut=0.01,
          mu_egene=0.003, mu_egenome=0.03, p_dup_egene=1.0,
          tax=0.03, tax_per_egene=0.0, restricted_mu=False)
_N = 64


def _h(a):
    return hashlib.sha1(np.ascontiguousarray(a).tobytes()).hexdigest()[:16]


def _trajectory(make_sim, *, shadow, seed=7, n_steps=600, every=100):
    sim = make_sim(_N, **{k: _P[k] for k in
                          ('food_inc', 'm_scale', 'gdiff', 'mu_lut',
                           'mu_egenome', 'tax')})
    sim.set_seed(seed)
    sim.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    sim.set_v(np.random.default_rng(seed).integers(0, 2, (_N, _N),
                                                    dtype=np.uint8))
    sim.set_f_all(0.5)
    sim.set_F_all(1.0)
    sim.set_alive_halfplane(0)
    if shadow:
        sim.neutral_enable()
        assert sim.neutral_is_enabled()
    traj = []
    for t in range(n_steps):
        sim.step()
        if t % every == every - 1:
            traj.append((_h(sim.get_v()), _h(sim.get_alive()),
                         _h(sim.get_f()), _h(sim.get_F()),
                         _h(sim.get_egenome()), _h(sim.get_lut_hashes())))
    alive = int(sim.get_alive().sum())
    shadow_pop = sim.neutral_population() if shadow else None
    return traj, alive, shadow_pop


def test_shadow_does_not_perturb_trajectory(make_sim):
    """Same seed, shadow on vs off -> bit-identical lattice trajectory.
    Fails on main (530c334), where shadow draws advanced the dynamics RNG."""
    off, alive_off, _ = _trajectory(make_sim, shadow=False)
    on, alive_on, shadow_pop = _trajectory(make_sim, shadow=True)
    assert alive_off > 0, "fixture went extinct; test would pass trivially"
    # The shadow really ran (and so really drew random numbers): its
    # population mirrors the real one.
    assert shadow_pop == alive_on > 0
    for i, (a, b) in enumerate(zip(off, on)):
        assert a == b, f"shadow perturbed dynamics at sample {i}: {a} != {b}"


def _run_sim(shadow, seed=3):
    from evoca_explore import run_sim
    return run_sim(dict(_P), n_steps=600, sample_every=100, N=_N,
                   seed=seed, shadow=shadow)


def test_run_sim_scalars_independent_of_shadow():
    """Every run_sim scalar that exists with shadow=False is identical
    with shadow=True (only the shadow-only keys are added)."""
    off = _run_sim(False)
    on = _run_sim(True)
    assert off['final_pop'] > 0, "fixture went extinct"
    diffs = {k: (off[k], on.get(k)) for k in off
             if k != 'shadow' and off[k] != on.get(k)}
    assert not diffs, f"shadow changed non-shadow metrics: {diffs}"
    # Shadow-only keys are present (the shadow did run).
    assert 'excess_pc_final' in on and 'excess_pc_final' not in off


def test_shadow_outputs_reproducible_same_seed():
    """Two shadow-on runs in the same process with the same seed agree
    on everything, including shadow-derived metrics. This also pins that
    set_seed reseeds the shadow stream: without it the second run would
    continue the first run's shadow stream and its shadow metrics differ."""
    a = _run_sim(True, seed=5)
    b = _run_sim(True, seed=5)
    assert a['final_pop'] > 0
    assert a == b, {k: (a[k], b[k]) for k in a if a[k] != b[k]}


def test_shadow_tables_reproducible_object_api(make_sim):
    """Object-API check of the full N-activity table (not only summary
    scalars): same seed -> identical shadow hash/activity/pop tables."""
    def shadow_table():
        sim = make_sim(_N, **{k: _P[k] for k in
                              ('food_inc', 'm_scale', 'gdiff', 'mu_lut',
                               'mu_egenome', 'tax')})
        sim.set_seed(11)
        sim.state(lut='gol', egenome='uniform', egenome_value=0b000011)
        sim.set_v(np.random.default_rng(11).integers(0, 2, (_N, _N),
                                                     dtype=np.uint8))
        sim.set_f_all(0.5); sim.set_F_all(1.0); sim.set_alive_halfplane(0)
        sim.neutral_enable()
        for _ in range(300):
            sim.step()
        sim.n_activity_update()
        Nt = sim.get_n_activity(max_n=200000)
        order = np.argsort(Nt['hash'])
        return ({k: Nt[k][order].tolist() for k in ('hash', 'activity',
                                                     'pop_count')},
                int(sim.get_alive().sum()))
    (t1, alive), (t2, _) = shadow_table(), shadow_table()
    assert alive > 0
    assert len(t1['hash']) > 1, "shadow never mutated; test is vacuous"
    assert t1 == t2


def test_frozen_shadow_is_clonal(make_sim):
    """Null sanity: with LUT mutation frozen and a uniform GoL LUT, the
    shadow stays a single wild-type bucket whose population tracks the
    real alive count (demography mirrored, no spurious variation from
    the new stream)."""
    sim = make_sim(_N, food_inc=_P['food_inc'], m_scale=_P['m_scale'],
                   gdiff=_P['gdiff'], mu_lut=0.0, mu_egenome=0.0,
                   tax=_P['tax'])
    sim.set_seed(2)
    sim.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    sim.set_v(np.random.default_rng(2).integers(0, 2, (_N, _N),
                                                dtype=np.uint8))
    sim.set_f_all(0.5); sim.set_F_all(1.0); sim.set_alive_halfplane(0)
    sim.neutral_enable()
    for _ in range(300):
        sim.step()
    alive = int(sim.get_alive().sum())
    assert alive > 0
    sim.n_activity_update()
    Nt = sim.get_n_activity(max_n=4096)
    live = Nt['pop_count'] > 0
    assert int(live.sum()) == 1
    assert int(Nt['pop_count'][live][0]) == alive == sim.neutral_population()
