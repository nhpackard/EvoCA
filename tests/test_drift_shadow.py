"""Drift shadow + neutral marker (Docs/shadows.md §5-7).

Every test enters through a public boundary: evoca_explore.run_sim, the
EvoCA object API, or the synthetic hook the drift statistics are
calibrated with (python/drift_synth.py). Fixtures are checked to be alive
at the end, because an extinct fixture passes trivially.

  observer  : fs_shadow modes and the marker leave the real trajectory and
              every tracking / Channon key bit-identical
  reach     : fs_shadow, fs_reps, fs_resync, marker_mu reach the run
  repro     : same seed -> identical drift keys
  demography: every replicate mirrors the real births/deaths each step,
              keeps n_r == alive count, and replays the real realized
              flip counts; incremental tallies == recount
  parity    : C replicates == pure-Python reference, tally-for-tally
  wiring    : clonal founder, all mu = 0 -> all drift statistics 0
"""
import ctypes
import hashlib
import json
import math
import os

import numpy as np
import pytest

# pure-evo joint winner (Research/Activity-metrics/Runs/
# 2026-10-15_pure_evo_joint_winner.evoca): survives its founder bottleneck
# at N=64 for every seed tried (min alive >= 24 over 1500 steps, 12 seeds),
# where smaller lattices go extinct for some seeds.
_P = dict(food_inc=0.018, m_scale=2.5, gdiff=0.08, mu_lut=0.01, mu_egene=0.003,
          mu_egenome=0.003, p_dup_egene=1.0, tax=0.035, tax_per_egene=0.0008,
          restricted_mu=False)
_N = 64
_DRIFT_PREFIXES = ('eg_drift_', 'dyn_drift_', 'marker_drift_', 'negene_drift',
                   'dyn_expression_bias', 'drift_')


def _h(a):
    return hashlib.sha1(np.ascontiguousarray(a).tobytes()).hexdigest()[:16]


def _build(seed=7, params=_P, N=_N):
    from evoca_explore import build_sim
    return build_sim(params, N=N, seed=seed)


def _genome_hashes(sim):
    N = sim._N
    egv = np.ctypeslib.as_array(sim._lib.evoca_get_egenes(), shape=(N * N * 8,))
    egm = np.ctypeslib.as_array(sim._lib.evoca_get_egenes_mask(), shape=(N * N * 8,))
    return (_h(sim.get_v()), _h(sim.get_alive()), _h(sim.get_f()), _h(sim.get_F()),
            _h(sim.get_lut_hashes()), _h(egv), _h(egm))


def _traj(*, drift_R=0, marker_mu=0.0, channon=False, seed=7, steps=500):
    sim = _build(seed)
    if channon:
        sim.neutral_enable()
    sim.set_marker_mu(marker_mu)
    if drift_R:
        sim.fs_drift_enable(drift_R)
    out = []
    for t in range(steps):
        sim.step()
        if t % 50 == 49:
            out.append(_genome_hashes(sim))
    alive = int(sim.get_alive().sum())
    sim.free()
    return out, alive


def _is_drift_key(k):
    return k.startswith(_DRIFT_PREFIXES) or k in ('fs_shadow', 'fs_reps',
                                                   'fs_resync', 'marker_mu')


# ── observer (P4) ─────────────────────────────────────────────────────

def test_drift_and_marker_do_not_perturb_trajectory():
    base, alive = _traj()
    assert alive > 50, "fixture went extinct; test would pass trivially"
    for kw in (dict(drift_R=8), dict(drift_R=3, marker_mu=0.05),
               dict(marker_mu=0.05), dict(drift_R=8, channon=True)):
        got, _ = _traj(**kw)
        if kw.get('channon'):
            ref, _ = _traj(channon=True)
        else:
            ref = base
        assert got == ref, f"real trajectory changed with {kw}"


def test_run_sim_modes_bit_identical_on_common_keys():
    from evoca_explore import run_sim
    kw = dict(n_steps=600, sample_every=100, N=_N, seed=3, shadow=True)
    tr = run_sim(_P, fs_shadow='tracking', **kw)
    bo = run_sim(_P, fs_shadow='both', **kw)
    dr = run_sim(_P, fs_shadow='drift', **kw)
    bm = run_sim(_P, fs_shadow='both', marker_mu=0.02, **kw)
    assert not tr['extinct']
    for k in tr:
        if not _is_drift_key(k):
            assert tr[k] == bo[k] == bm[k], k
            if not (k.startswith(('eg_excess', 'dyn_excess',
                                  'eg_concentration', 'dyn_concentration'))):
                assert dr[k] == tr[k], k
    for k in bo:
        if (_is_drift_key(k) and not k.startswith('marker')
                and k not in ('fs_shadow', 'fs_reps', 'fs_resync')):
            assert dr[k] == bo[k] or (isinstance(dr[k], float) and
                                      math.isnan(dr[k]) and math.isnan(bo[k])), k


# ── reachability through run_sim ──────────────────────────────────────

def test_fs_shadow_modes_reach_run_sim():
    from evoca_explore import run_sim, _PARAM_KEYS
    kw = dict(n_steps=200, sample_every=50, N=_N, seed=1, shadow=False)
    tr = run_sim(_P, fs_shadow='tracking', **kw)
    dr = run_sim(_P, fs_shadow='drift', fs_reps=5, **kw)
    bo = run_sim(_P, **kw)
    assert bo['fs_shadow'] == 'both' and bo['fs_reps'] == 8 and bo['drift_R'] == 8
    assert 'eg_excess_pc_final' in tr and 'eg_drift_concentration' not in tr
    assert tr['fs_reps'] == 0
    assert 'eg_excess_pc_final' not in dr and 'eg_concentration_final' not in dr
    assert dr['drift_R'] == 5 and dr['fs_reps'] == 5
    assert dr['drift_steps'] == 200 and dr['drift_steps_last'] == 100
    assert 'eg_excess_pc_final' in bo and 'dyn_drift_consistency' in bo
    assert bo['eg_concentration_final'] == bo['eg_excess_pc_final']
    assert bo['dyn_concentration_slope'] == bo['dyn_excess_pc_slope']
    with pytest.raises(ValueError):
        run_sim(_P, fs_shadow='drfit', **kw)
    for k in ('fs_shadow', 'fs_reps', 'fs_resync', 'marker_mu'):
        assert k not in _PARAM_KEYS


def test_fs_resync_and_marker_mu_reach_run_sim():
    from evoca_explore import run_sim
    kw = dict(n_steps=300, sample_every=100, N=_N, seed=2, shadow=False,
              fs_shadow='drift', fs_reps=4)
    r0 = run_sim(_P, **kw)
    assert not r0['extinct']
    rs = run_sim(_P, fs_resync=70, **kw)
    assert rs['fs_resync'] == 70
    assert rs['drift_steps_last'] == 300 - 280        # last window [280, 300)
    assert r0['drift_steps_last'] == 150
    assert rs['eg_drift_tv'] != r0['eg_drift_tv']      # reseeding took effect
    # marker: clonal at mu 0 (every marker statistic exactly 0) ...
    assert r0['marker_mu'] == 0.0
    assert r0['marker_drift_tv'] == 0.0 and r0['marker_drift_concentration'] == 0.0
    # ... and varying when mu > 0
    rm = run_sim(_P, marker_mu=0.02, **kw)
    assert rm['marker_mu'] == 0.02 and rm['marker_drift_tv'] > 0.0


# ── reproducibility (P5) ──────────────────────────────────────────────

def test_drift_keys_reproducible_with_seed():
    from evoca_explore import run_sim
    kw = dict(n_steps=300, sample_every=100, N=_N, shadow=False,
              fs_shadow='drift', fs_reps=4, marker_mu=0.01)
    a = run_sim(_P, seed=5, **kw)
    assert not a['extinct']
    b = run_sim(_P, seed=5, **kw)
    c = run_sim(_P, seed=6, **kw)
    dk = [k for k in a if _is_drift_key(k)]
    assert len(dk) > 30
    for k in dk:
        assert a[k] == b[k] or (math.isnan(a[k]) and math.isnan(b[k])), k
    assert a['eg_drift_tv'] != c['eg_drift_tv']


def test_replicate_streams_distinct_and_seeded():
    sim = _build(seed=9)
    sim.fs_drift_enable(4)
    for _ in range(200):
        sim.step()
    tl = sim.fs_drift_tallies()
    assert tl['n'].min() > 50
    assert len({_h(tl['c'][r]) for r in range(4)}) == 4, "replicates share a stream"
    sim.free()


# ── demography and flip-count replay ─────────────────────────────────

@pytest.mark.parametrize('restricted', [False, True])
def test_demography_and_flip_replay_match_real(restricted):
    P = dict(_P, restricted_mu=restricted)
    sim = _build(seed=4, params=P)
    sim.set_marker_mu(0.02)
    sim.fs_drift_enable(5)
    nm_tot = 0
    for t in range(400):
        sim.step()
        L = sim._lib
        B, D = L.evoca_get_births_last(), L.evoca_get_deaths_last()
        ev = sim.fs_drift_last_events()
        assert np.all(ev[:, 0] == B) and np.all(ev[:, 1] == D), t
        assert np.all(ev[:, 2] == L.evoca_get_lut_bitflips_last()), t
        assert np.all(ev[:, 3] + ev[:, 4] == L.evoca_get_eg_bitflips_last()), t
        assert np.all(ev[:, 5] == ev[0, 5]), t
        nm_tot += int(ev[0, 5])
        assert np.all(sim.fs_drift_tallies()['n'] == int(sim.get_alive().sum())), t
    assert nm_tot > 0, "marker flips were not replayed"
    assert sim.fs_drift_check_tallies() == 0
    assert int(sim.get_alive().sum()) > 50
    sim.free()


def test_real_allele_counts_match_recount():
    sim = _build(seed=8)
    sim.fs_drift_enable(2)
    for _ in range(150):
        sim.step()
    ro = sim.fs_drift_real_obs()
    N = sim._N
    lut = np.ctypeslib.as_array(sim._lib.evoca_get_lut(), shape=(N * N, 32))
    al = sim.get_alive().ravel() == 1
    bits = np.unpackbits(lut[al], axis=1, bitorder='little')[:, :250].sum(0)
    assert ro['Np'] == al.sum() > 0
    assert np.array_equal(bits, ro['c'])
    mk = sim.get_marker()[al]
    from drift_ref import marker_key
    keys = np.bincount([marker_key(int(m)) for m in mk], minlength=729)
    assert np.array_equal(keys, ro['m'])
    sim.free()


# ── C vs pure-Python reference (bit parity) ───────────────────────────

def test_c_replicates_match_python_reference():
    from drift_ref import Replicate, drift_seed_for
    from evoca_py import EvoCA
    lib = EvoCA()._lib
    rng = np.random.default_rng(12)
    n0, R, steps, seed = 40, 3, 1000, 777
    lut = rng.integers(0, 256, (n0, 32), dtype=np.uint8)
    lut[:, 31] &= 0x03                                  # bits >= 250 unused
    eg = rng.integers(0, 64, (n0, 8), dtype=np.uint8)
    egm = rng.integers(0, 64, (n0, 8), dtype=np.uint8)
    act = rng.integers(1, 256, n0, dtype=np.uint8)
    mk = (rng.integers(0, 64, n0) | (rng.integers(0, 64, n0) << 8)).astype(np.uint16)
    lib.evoca_set_p_dup_egene(ctypes.c_float(0.5))
    U8 = ctypes.POINTER(ctypes.c_uint8)
    CI = ctypes.POINTER(ctypes.c_int)
    cap = 200
    assert lib.evoca_fs_drift_synth_init(
        R, cap, n0, lut.ctypes.data_as(U8), eg.ctypes.data_as(U8),
        egm.ctypes.data_as(U8), act.ctypes.data_as(U8),
        mk.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16)), seed) == R
    refs = [Replicate(lut, eg, egm, act, mk, drift_seed_for(seed, r), p_dup=0.5)
            for r in range(R)]
    c = np.zeros((R, 250), np.int32); e = np.zeros((R, 729), np.int32)
    m = np.zeros((R, 729), np.int32); n = np.zeros(R, np.int32)
    T = np.zeros(R, np.int64)
    I32 = ctypes.POINTER(ctypes.c_int32)
    nq = np.ones(250)
    npop = n0
    for t in range(steps):
        D = int(rng.integers(0, 5))
        B = int(np.clip(D + rng.integers(-1, 2), 0, None))
        if npop - min(D, npop) + B < 20:
            B += 2
        if npop - min(D, npop) + B > 120:
            B = max(B - 2, 0)
        npop = npop - min(D, npop) + B
        recs = np.stack([rng.poisson(1.5, B), rng.poisson(0.4, B),
                         rng.poisson(0.8, B), rng.poisson(0.5, B)], 1).astype(np.int32)
        recs = np.ascontiguousarray(recs)
        Da, Ba = np.array([D], np.int32), np.array([B], np.int32)
        assert lib.evoca_fs_drift_synth_run(
            1, Da.ctypes.data_as(CI), Ba.ctypes.data_as(CI), recs.ctypes.data_as(CI),
            nq.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 0) == 0
        for ref in refs:
            ref.step(D, [tuple(int(x) for x in row) for row in recs])
        lib.evoca_fs_drift_tallies(c.ctypes.data_as(I32), e.ctypes.data_as(I32),
                                   m.ctypes.data_as(I32), n.ctypes.data_as(I32),
                                   T.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)))
        for r, ref in enumerate(refs):
            assert n[r] == len(ref.lut), (t, r)
            assert np.array_equal(c[r], ref.c), (t, r, 'c')
            assert np.array_equal(e[r], ref.e), (t, r, 'e')
            assert np.array_equal(m[r], ref.m), (t, r, 'm')
            assert T[r] == ref.T, (t, r, 'T')
    assert n.min() > 5, "population collapsed; parity would be trivial"
    lib.evoca_fs_drift_disable()
    lib.evoca_set_p_dup_egene(ctypes.c_float(1.0))


# ── clonal wiring and statistics behavior ─────────────────────────────

def test_clonal_frozen_run_reads_exactly_zero():
    from evoca_explore import run_sim
    P = dict(_P, mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0)
    out = run_sim(P, n_steps=300, sample_every=100, N=_N, seed=1, shadow=False,
                  fs_shadow='drift', fs_reps=4)
    assert not out['extinct']
    for k in ('eg_drift_concentration', 'eg_drift_z', 'eg_drift_tv', 'eg_drift_nhi',
              'dyn_drift_consistency', 'dyn_drift_z', 'dyn_drift_concentration',
              'dyn_drift_tv', 'dyn_drift_nhi', 'marker_drift_concentration',
              'marker_drift_tv', 'negene_drift', 'eg_drift_tv_z', 'dyn_drift_tv_z',
              'dyn_expression_bias'):
        assert out[k] == 0.0, (k, out[k])
    for k in ('eg_drift_p', 'dyn_drift_p', 'eg_drift_tv_p', 'dyn_drift_tv_p'):
        assert out[k] == 1.0, k


def test_synthetic_strong_selection_is_detected_and_neutral_is_not():
    from drift_synth import synth_run
    from drift_stats import drift_statistics
    R = 9
    sel = drift_statistics(synth_run(1, R=R, Np=300, b=15, steps=1500, mu_lut=0.01,
                                     sel_mode=1, s=0.3, loci=[0, 1]))
    assert sel['dyn_drift_p'] == pytest.approx(1 / (R + 1))
    assert sel['dyn_drift_consistency'] > 3
    eg = drift_statistics(synth_run(2, R=R, Np=300, b=15, steps=1500, mu_eg=0.01,
                                    sel_mode=2, s=0.5))
    assert eg['eg_drift_p'] == pytest.approx(1 / (R + 1))
    assert eg['eg_drift_concentration'] > 0


def _forced_run(params, force, rank=0, seed=11):
    """Run 100 steps, pick the rank-th most-queried input i*, then
    (force=None) leave it, or force output `force(major)` at i* in every
    birth; re-baseline the drift shadow at the switch and run 1500 steps.
    Returns (stats, zeta at i*)."""
    from drift_stats import drift_statistics, consistency_stats
    sim = _build(seed=seed, params=params)
    sim.fs_drift_enable(8)
    for _ in range(100):
        sim.step()
    acc = sim.fs_drift_get()
    q = acc['dyn_queried'].reshape(250, 2).sum(1)
    i_star = int(np.argsort(-q)[rank])
    major = int(np.argmax(acc['dyn_real'].reshape(250, 2)[i_star]))
    if force is not None:
        sim.test_force_bit(i_star, force(major))
    sim.fs_drift_enable(8)                 # re-baseline at the switch
    for _ in range(1500):
        sim.step()
    acc = sim.fs_drift_get()
    sim.test_force_bit(-1, 0)
    alive = int(sim.get_alive().sum())
    sim.free()
    assert alive > 50
    zeta = consistency_stats(acc['dyn_real'], acc['dyn_shadow'])[2]
    return drift_statistics(acc), float(zeta[i_star])


def test_forced_output_fixture_purifying_detected_at_target():
    """Constructed EvoCA fixture (doc §7(b)): holding the majority output at
    the most-queried input in every birth is non-neutral (purifying-like).
    The target input's consistency deviation must jump (zeta > 3) relative
    to the same seed unforced (zeta < 2). DEVIATION from the doc's
    criterion (aggregate dyn_drift_z > 5): in an evolving base the other
    inputs diverge from the founder faster than neutral (negative zeta),
    which outweighs one purified input in the aggregate."""
    st_f, z_f = _forced_run(_P, force=lambda major: major)
    st_u, z_u = _forced_run(_P, force=None)
    assert z_f > 3, z_f
    assert z_u < 2, z_u
    assert st_f['dyn_drift_consistency'] > st_u['dyn_drift_consistency']


def test_forced_output_fixture_sweep_detected_two_sided():
    """Otherwise-neutral fixture: frozen clonal genome, forcing the MINORITY
    output at a queried input (a directional sweep inside the window).
    Cumulative homozygosity at that input drops, so the one-sided
    consistency statistic reads NEGATIVE (design finding); the two-sided
    TV test catches it at the rank-p floor."""
    P = dict(_P, mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0)
    st, z = _forced_run(P, force=lambda major: 1 - major)
    assert z < 0 and st['dyn_drift_consistency'] < 0
    assert st['dyn_drift_tv_p'] == pytest.approx(1 / 9)
    assert st['dyn_drift_tv_z'] > 5, st['dyn_drift_tv_z']
    assert st['eg_drift_tv'] == 0.0          # nothing else departs


# ── recipe / scan measurement block ──────────────────────────────────

def test_measurement_block_in_recipes(tmp_path):
    from evoca_explore import save_scan_config
    from evoca_py import import_run
    p = save_scan_config(str(tmp_path), N=32, n_steps=10, sample_every=5, seed=0,
                         shadow=True, fs_reps=19, marker_mu=0.003)
    cfg = json.load(open(p))
    assert cfg['measurement'] == {'fs_shadow': 'both', 'fs_reps': 19,
                                  'fs_resync': None, 'marker_mu': 0.003}
    sim = _build(seed=1, N=32)
    sim.set_marker_mu(0.004)
    sim.fs_drift_enable(6)
    path = sim.export_recipe('drift_rt', runs_dir=str(tmp_path))
    sim.free()
    rec = json.load(open(path))
    assert rec['version'] == 4
    assert rec['measurement']['fs_shadow'] == 'both'
    assert rec['measurement']['fs_reps'] == 6
    assert rec['measurement']['marker_mu'] == pytest.approx(0.004)
    sim2, _ = import_run(path)
    assert sim2.get_marker_mu() == pytest.approx(0.004)
    assert sim2._measurement['fs_reps'] == 6
    assert not sim2.fs_drift_is_enabled()
    sim2.free()


def test_evoca_from_scan_writes_measurement(tmp_path):
    from evoca_explore import save_scan_config, evoca_from_scan
    sd = tmp_path / 'scan'
    sd.mkdir()
    save_scan_config(str(sd), N=32, n_steps=10, sample_every=5, seed=0,
                     shadow=True, fs_shadow='drift', fs_reps=8)
    (sd / 'results.csv').write_text('config_idx,param_food_inc,param_food_cap\n0,0.01,1\n')
    path = evoca_from_scan(str(sd), config_idx=0, runs_dir=str(tmp_path))
    rec = json.load(open(path))
    assert rec['version'] == 4 and rec['measurement']['fs_shadow'] == 'drift'
    os.remove(sd / 'scan_config.json')
    cfg = {'N': 32, 'initialization': {}}
    (sd / 'scan_config.json').write_text(json.dumps(cfg))
    path = evoca_from_scan(str(sd), config_idx=0, runs_dir=str(tmp_path),
                           descriptor='old')
    assert json.load(open(path))['measurement']['fs_shadow'] == 'tracking'
