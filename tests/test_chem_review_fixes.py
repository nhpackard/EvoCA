"""Regression tests for the chem-branch review fixes (Research/Chem/plans/plan-chem.md,
"Review fixes"). Each enters through the public door a caller uses
(stamp_patch, set_f_all, set_alive*, set_w, set_n_chem, set_chem_energy,
init, the controls Restart path, colorize, params_from_row) on a constructed
fixture, and asserts the waste-mass identity

    d(sum W + sum w) == d(in) - d(out)        (get_chem_flux(1))

and, where energy moves, the energy identity

    source_step == d(storage_env + storage_priv) + diss_step.

The still world (_still): every cell alive, f = 0.5, no food regrowth,
eating, diffusion or mutation, no tax unless asked; waste_disp = 0 so a
seeded private waste stays put until a test moves it.
"""
import hashlib
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))

WASTE = dict(n_chem=2, wdiff=0.3, waste_sink=0.01, waste_disp=0.2)


def _still(make_sim, N=8, **kw):
    p = dict(food_inc=0.0, m_scale=0.0, gdiff=0.0, n_chem=2)
    p.update(kw)
    s = make_sim(N, **p)
    s.state(lut='gol', egenome='uniform', egenome_value=3)
    s.set_F_all(0.0); s.set_alive_all(); s.set_f_all(0.5)
    return s


def _mass(s):
    return s.get_W().astype(np.float64).sum() + s.get_w().astype(np.float64).sum()


def _net_in(s):
    cin, cout = s.get_chem_flux(1)
    return cin - cout


class _Balance:
    """Waste-mass balance across an operation: check() asserts that the
    change in sum W + sum w equals the booked in - out since the last mark."""
    def __init__(self, s):
        self.s = s; self.mark()

    def mark(self):
        self.m, self.n = _mass(self.s), _net_in(self.s)

    def check(self, abs_tol=1e-4):
        dm, dn = _mass(self.s) - self.m, _net_in(self.s) - self.n
        assert dm == pytest.approx(dn, abs=abs_tol), f"waste mass moved {dm:+.6f}, booked {dn:+.6f}"
        self.mark()


def _fingerprint(s):
    h = hashlib.sha1()
    for a in (s.get_F(), s.get_f(), s.get_v(), s.get_alive(), s.get_lut_hashes(), s.get_W(), s.get_w()):
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


_STATE = dict(lut='gol', egenome='uniform', egenome_value=3, v_seed=7,
              f_init=0.5, F='uniform', F_init=1.0, alive='halfplane', alive_axis=0)
_RUN = dict(food_inc=0.015, m_scale=3.0, gdiff=0.08, tax=0.02)


# ── stamp_patch (review: cloud 3 / local 1, local 5) ──────────────────

def test_stamp_overwrite_dead_cells_return_host_waste(make_sim):
    """Patch cells marked empty: the host's private waste returns to its site
    (not left in the dead cell to vanish under the next birth)."""
    s = _still(make_sim)
    s.set_f_all(1.0)                    # alive cells will reproduce next step
    s.set_w(0.3)
    patch = s.extract_patch(0, 0, 4)
    patch.alive[:] = 0
    b = _Balance(s)
    s.stamp_patch(patch, 0, 0, mode='overwrite')
    b.check()
    w, W = s.get_w(), s.get_W()
    assert np.all(w[:4, :4] == 0), "a dead cell kept private waste"
    assert np.allclose(W[:4, :4], 0.3)
    s.step()                            # births land on the cleared cells
    assert s.get_alive().reshape(8, 8)[:4, :4].sum() > 0
    b.check()


def test_stamp_overwrite_live_cells_start_without_host_waste(make_sim):
    s = _still(make_sim)
    s.set_w(0.3)
    patch = s.extract_patch(4, 4, 4)    # all alive
    b = _Balance(s)
    s.stamp_patch(patch, 0, 0, mode='overwrite')
    b.check()
    w, W = s.get_w(), s.get_W()
    assert np.all(s.get_alive().reshape(8, 8)[:4, :4] == 1)
    assert np.all(w[:4, :4] == 0), "a stamped organism inherited the host's waste"
    assert np.allclose(W[:4, :4], 0.3) and np.allclose(w[4:, 4:], 0.3)


# ── set_f_all (local 4) ───────────────────────────────────────────────

def test_set_f_all_keeps_private_waste(make_sim):
    s = _still(make_sim)
    s.set_w(0.3)
    b = _Balance(s)
    s.set_f_all(0.7)
    b.check()
    assert np.allclose(s.get_w(), 0.3) and np.allclose(s.get_f(), 0.7)


# ── set_alive* (local 6) ──────────────────────────────────────────────

_MASK = np.zeros((8, 8), np.uint8); _MASK[2:6, 1:7] = 1


@pytest.mark.parametrize('call', [
    lambda s: s.set_alive(_MASK),
    lambda s: s.set_alive_patch(2),
    lambda s: s.set_alive_halfplane(0),
    lambda s: s.set_alive_fraction(0.0),
], ids=['set_alive', 'patch', 'halfplane', 'fraction'])
def test_set_alive_returns_dying_cells_waste(make_sim, call):
    s = _still(make_sim)
    s.set_w(0.3)
    b = _Balance(s)
    call(s)
    b.check()
    dead = s.get_alive().reshape(8, 8) == 0
    assert dead.any()
    assert np.all(s.get_w()[dead] == 0)
    assert np.allclose(s.get_W()[dead], 0.3)


# ── set_w (local 7) ───────────────────────────────────────────────────

def test_set_w_is_booked_in_the_flux(make_sim):
    s = _still(make_sim)
    b = _Balance(s)
    s.set_w(0.3)                        # a rise: booked as input
    b.check()
    w = np.full((8, 8), 0.3, np.float32); w[:4] = 0.1; w[4:, :2] = 0.5
    s.set_w(w)                          # falls and rises: out and in
    b.check()
    assert np.allclose(s.get_w(), w)


# ── energy weight of compounds c >= 1 (cloud 2 / local 2) ─────────────

def _energy_identity_steps(s, n):
    d = s.get_dissipation(); prev = d['storage_env'] + d['storage_priv']
    for _ in range(n):
        s.step()
        d = s.get_dissipation(); cur = d['storage_env'] + d['storage_priv']
        assert d['source_step'] == pytest.approx(cur - prev + d['diss_step'], abs=1e-4)
        prev = cur


def test_set_chem_energy_refuses_an_energy_carrying_waste(make_sim, capfd):
    s = _still(make_sim, N=16, tax=0.05)
    s.set_chem_energy(1, 1.0)
    assert s.get_chem_energy(1) == 0.0
    assert 'unsupported for compound 1' in capfd.readouterr().err
    _energy_identity_steps(s, 5)
    assert s.get_f().mean() == pytest.approx(0.5 - 5 * 0.05, abs=1e-5), "the tax must cost energy"


def test_low_level_metabolism_defaults_to_energy_free_product(make_sim):
    """The C path set_n_chem(2) + set_metab_product(1), with no preset: the
    product must default to energy 0 (it was 1.0, making the tax free)."""
    s = make_sim(16, food_inc=0.0, m_scale=0.0, tax=0.05)
    s.state(lut='gol', egenome='uniform', egenome_value=3)
    s.set_F_all(0.0); s.set_alive_all(); s.set_f_all(0.5)
    s._lib.evoca_set_n_chem(2)
    assert s.get_chem_energy(1) == 0.0
    s._lib.evoca_set_metab_product(1)
    _energy_identity_steps(s, 5)
    assert s.get_f().mean() == pytest.approx(0.25, abs=1e-5)


# ── set_n_chem keeps the Python view faithful (cloud 4 / local 3) ─────

def test_set_n_chem_1_resets_waste_params_and_round_trips(make_sim):
    from evoca_py import EvoCA
    s = make_sim(16, **WASTE)
    s.set_n_chem(1)
    p = s.params()
    assert (p['n_chem'], p['wdiff'], p['waste_sink'], p['waste_disp']) == (1, 0.0, 0.0, 0.0)
    s.free()
    s2 = EvoCA(); s2.init(**p)              # was: ValueError (waste params need n_chem=2)
    assert s2.get_n_chem() == 1
    s2.free()


def test_set_n_chem_2_record_replays_the_same_run(make_sim):
    """A run switched to n_chem=2 after init must replay from params() and
    state_params into the same dynamics (was: inert C chemistry recorded as
    the food/waste preset, so the replay had metabolism and toxicity)."""
    from evoca_py import EvoCA
    s = make_sim(24, **_RUN)
    s.set_n_chem(2)
    assert s.get_metab_product() == 1 and s.get_chem_energy(1) == 0.0
    s.set_seed(3); s.state(**_STATE)
    b = _Balance(s)
    for _ in range(40):
        s.step(); b.check(abs_tol=2e-3)
    assert s.get_W().sum() > 0
    fp, p, sp = _fingerprint(s), s.params(), s.state_params
    s.free()
    r = EvoCA(); r.init(**p); r.set_seed(3); r.state(**sp)
    for _ in range(40):
        r.step()
    assert _fingerprint(r) == fp
    r.free()


def test_set_n_chem_2_keeps_the_configured_waste_params(make_sim):
    s = make_sim(16, **WASTE)
    s.set_n_chem(2)
    assert s.get_chem_params(1) == pytest.approx((0.0, 0.01, 0.3, 0.2, 1.0))
    assert (s.wdiff, s.waste_sink, s.waste_disp) == pytest.approx((0.3, 0.01, 0.2))


def test_set_n_chem_refuses_unconfigured_sizes(make_sim):
    s = make_sim(8, **WASTE)
    with pytest.raises(ValueError):
        s.set_n_chem(3)
    assert s.get_n_chem() == 2 and s.n_chem == 2 and s.get_metab_product() == 1


# ── init validates before any side effect (cloud 5) ───────────────────

@pytest.mark.parametrize('bad', [dict(n_chem=3), dict(wdiff=0.1),
                                 dict(n_chem=2, wdiff=float('nan'))],
                         ids=['n_chem=3', 'waste-without-n_chem', 'nan'])
def test_refused_init_leaves_the_instance_untouched(make_sim, bad):
    s = make_sim(16, **WASTE, **_RUN)
    s.state(**_STATE)
    before, fp = s.params(), _fingerprint(s)
    with pytest.raises(ValueError):
        s.init(32, **_RUN, **bad)
    assert s.params() == before and s.N == 16
    assert s.get_alive().size == 256 and _fingerprint(s) == fp
    b = _Balance(s)
    s.step(); b.check()


# ── controls Restart re-applies the full parameter set (cloud 1) ──────

def test_restart_reapplies_every_metaparameter(make_sim):
    from controls import restart_sim
    from evoca_py import EvoCA
    P = dict(**WASTE, **_RUN, food_cap=0, food_per=3, lineage=True)
    s = make_sim(16, **P)
    s.set_seed(5); s.state(**_STATE)
    for _ in range(20):
        s.step()
    s.set_seed(5); restart_sim(s)
    assert (s.get_n_chem(), s.get_metab_product(), s.get_food_cap(), s._lib.evoca_get_food_per()) == (2, 1, 0, 3)
    assert s.get_chem_params(1) == pytest.approx((0.0, 0.01, 0.3, 0.2, 1.0))
    assert s.get_chem_energy(1) == 0.0
    b = _Balance(s)
    for _ in range(30):
        s.step(); b.check(abs_tol=2e-3)
    fp = _fingerprint(s)
    s.free()
    r = EvoCA(); r.init(16, **P); r.set_seed(5); r.state(**_STATE)
    for _ in range(30):
        r.step()
    assert _fingerprint(r) == fp, "a restarted run must equal a fresh one"
    r.free()


# ── colormode 6 shows the waste field without metabolism (local 8) ────

def test_food_waste_colormode_shows_seeded_waste_without_metabolism(make_sim):
    s = _still(make_sim)
    s._lib.evoca_set_metab_product(-1)
    s.set_w(0.5)
    s.set_alive_halfplane(0)            # the dead half returns its waste to W
    W = s.get_W()
    assert W.max() > 0
    px = s.colorize(np.zeros(64, np.int32), 6).view(np.uint32).reshape(8, 8)
    dead = s.get_alive() == 0
    blue = px & 0xFF                     # colormode 6: env waste is blue on dead cells
    assert np.all(blue[dead & (W > 0)] > 0) and np.all(blue[dead & (W == 0)] == 0)


# ── params_from_row with NaN cells (local 9) ──────────────────────────

def test_params_from_row_drops_nan_cells():
    import pandas as pd
    from evoca_nb import params_from_row
    from evoca_py import EvoCA
    row = pd.Series({'param_food_inc': 0.01, 'param_m_scale': 1.0, 'param_tax': 0.01,
                     'param_n_chem': np.nan, 'param_wdiff': np.nan,
                     'param_waste_sink': np.nan, 'param_waste_disp': np.nan})
    p = params_from_row(row)
    assert not any(k in p for k in ('n_chem', 'wdiff', 'waste_sink', 'waste_disp'))
    s = EvoCA(); s.init(8, **p)
    assert s.get_n_chem() == 1
    s.free()


# ── per-compound accessors (local 12) ─────────────────────────────────

def test_per_compound_accessors(make_sim):
    s = make_sim(8)
    assert np.all(s.get_chem_env(1) == 0) and s.get_chem_env(1).shape == (8, 8)   # absent: zeros, not a NULL read
    assert np.all(s.get_W() == 0) and np.all(s.get_w() == 0)
    s.set_F_all(0.25)
    assert np.allclose(s.get_chem_env(0), 0.25)
    s._lib.evoca_set_n_chem(3)
    s.set_alive_all()
    a = np.zeros(64, np.float32); a[9] = 0.4
    s._lib.evoca_set_chem_priv(2, a.ctypes.data_as(__import__('ctypes').POINTER(__import__('ctypes').c_float)))
    assert s.get_chem_priv(2)[1, 1] == pytest.approx(0.4) and s.get_chem_priv(2).sum() == pytest.approx(0.4)
    assert s.get_chem_env(2).shape == (8, 8) and np.all(s.get_chem_env(7) == 0)


# ── environmental source booked once per step (local 11) ──────────────

def test_env_source_is_booked_as_inc_times_cells(make_sim):
    """Refactor guard (booking moved out of the per-cell loop): behavior-
    preserving, so this pins the value rather than failing without the fix."""
    s = _still(make_sim)
    s._lib.evoca_set_chem_env(1, 0.01, 0.0, 0.0)
    b = _Balance(s)
    for _ in range(5):
        s.step()
    assert s.get_chem_flux(1)[0] == pytest.approx(0.01 * 64 * 5, rel=1e-6)
    b.check()
