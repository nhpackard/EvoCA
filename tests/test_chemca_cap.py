"""ChemCA open-mode source cap cap_<c> (EvoCA's F <= 1 food rule), 2026-11-02.

After the source fires a site holds at most cap; the excess is discarded and
booked as outflow. cap = 0 (the default) is off and bit-identical to the
engine before the cap existed (goldens recorded on that engine)."""
import numpy as np
import pytest

import evoca_explore
from chemca_py import ChemCA, run_chemca, load_recipe

PORT = dict(preset='fw2', open_mode=1, light=0.0, m_e=0.05, u_F=0.75, u_W=0.09, tox_W=1.6,
            inc_F=0.012, sink_F=0.004, sink_W=0.1, diff_F=0.08, diff_W=0.01,
            mu_rgene=0.02, mu_ractive=0.02)


def _hash(N, steps, state, **p):
    s = ChemCA()
    s.init(N, **p)
    s.state(seed=7, **state)
    s.step(steps)
    pop, h = s.population(), s.state_hash()
    s.free()
    return pop, h


# ── off is unchanged ─────────────────────────────────────────────────

def test_cap_off_matches_pre_cap_engine():
    pop, h = _hash(64, 25, dict(mass_split=(1.0, 0.0), store_frac=0.0, k_idx0=63), **PORT)
    assert pop == 1911                                   # alive and reproducing
    assert h == '0b99e3e7b573a9345d9986f805e4c0b6ed1da58d0d1b1719b7a63cf5ac7650d4'
    _, h = _hash(32, 300, {}, preset='gs3', mu_rgene=0.02, mu_ractive=0.02)
    assert h == '1ba146f92c467e1f15dbda608223c8a2b5483f6f5692adf90ab18cf0c8ae595a'


# ── the rule, on a constructed field ─────────────────────────────────

def _empty_open(N=4, **p):
    """Open mode, no cells, no diffusion, no light: only source/sink/cap act."""
    s = ChemCA()
    s.init(N, preset='fw2', open_mode=1, light=0.0, diff_F=0.0, diff_W=0.0, **p)
    s.set_alive(np.zeros((N, N)))
    return s


def test_cap_clips_after_the_source_and_books_the_excess():
    s = _empty_open(inc_F=0.012, cap_F=1.0)
    E = np.full((4, 4), 0.3, np.float32)
    E[0, 0], E[0, 1] = 0.995, 1.5               # just under the cap; above it already
    s.set_E('F', E)
    s.step()
    F = s.get_E('F')
    assert F[0, 0] == pytest.approx(1.0) and F[0, 1] == pytest.approx(1.0)
    assert F[1:, :] == pytest.approx(0.312)
    excess = (0.995 + 0.012 - 1.0) + (1.5 + 0.012 - 1.0)
    assert s.mass_flux()[s.cidx('F'), 1] == pytest.approx(excess, rel=1e-5)
    assert s.energy_budget()['open_out'] == pytest.approx(excess * 1.0, rel=1e-5)   # G_F = 1


def test_cap_acts_without_a_source():
    """As in EvoCA (the clip runs even at food_inc = 0)."""
    s = _empty_open(inc_F=0.0, cap_F=1.0)
    s.set_E('F', 2.0)
    s.step()
    assert s.get_E('F') == pytest.approx(1.0)


def test_cap_then_sink():
    s = _empty_open(inc_F=0.5, cap_F=1.0, sink_F=0.1)
    s.set_E('F', 0.8)
    s.step()
    assert s.get_E('F') == pytest.approx(1.0 * np.exp(-0.1), rel=1e-6)


# ── through the public doors ─────────────────────────────────────────

def test_cap_reaches_C_through_build_sim_and_update():
    s = evoca_explore.build_sim(dict(model='chemca', **PORT, cap_F=1.0), N=16)
    assert s.compound_params('F')['cap'] == 1.0
    assert s.compound_params('W')['cap'] == 0.0
    s.update(cap_F=0.5)
    assert s.compound_params('F')['cap'] == 0.5
    s.free()


def test_cap_round_trips_through_a_recipe(tmp_path):
    s = ChemCA()
    s.init(16, **PORT, cap_F=1.0)
    s.state(seed=3)
    p = s.export_recipe('cap', runs_dir=str(tmp_path))
    s.free()
    t, _ = load_recipe(p)
    assert t.params()['cap_F'] == 1.0 and t.compound_params('F')['cap'] == 1.0
    t.free()


def test_mass_identity_holds_with_cap():
    r = run_chemca(dict(model='chemca', **{**PORT, 'sink_F': 0.0}, cap_F=1.0), n_steps=300, N=32,
                   seed=1, kymo_steps=50)
    assert r['mass_err_max'] < 1e-5
