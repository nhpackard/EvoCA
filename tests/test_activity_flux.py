"""Activity tables and the activity flux are never truncated (2026-11-03).

Before: the C table readers and the crossings probe stopped at their buffer
(max_n) and returned what fit, with nothing to signal it. The flux summed
the crossings buffer (65536), so at N = 512 it pinned at 65536 / 8 = 8192
in every condition (Research/Chem/Scans/2026-10-24_N7_evolution_waves).
Now the readers return the total count and the Python wrappers re-read
with a buffer that fits; the flux is summed in C with no buffer."""
import numpy as np
import pytest


@pytest.fixture
def diverse(make_sim):
    """A world with many genomes and the activity table updated each tick,
    as the GUI and the flux runs do."""
    s = make_sim(48, mu_lut=0.02, food_inc=0.02, tax=0.01)
    s.set_seed(1)
    s.set_lut_random(2)
    s.set_egenome_random()
    s.set_alive_all()
    s.set_f_all(0.8)
    for _ in range(60):
        s.step()
        s._lib.evoca_activity_update()
    return s


def test_get_activity_returns_the_whole_table(diverse):
    full = diverse.get_activity(max_n=1 << 20)
    assert len(full['hash']) > 50
    small = diverse.get_activity(max_n=5)
    assert len(small['hash']) == len(full['hash'])
    assert np.array_equal(np.sort(small['hash']), np.sort(full['hash']))


def test_shadow_and_eff_tables_return_whole(make_sim):
    s = make_sim(48, mu_lut=0.02, food_inc=0.02, tax=0.01)
    s.set_seed(1)
    s.set_lut_random(2)
    s.set_egenome_random()
    s.set_alive_all()
    s.set_f_all(0.8)
    s.neutral_enable()
    for _ in range(60):
        s.step()
        s.n_activity_update()
        s.eff_activity_update()
    for get in (s.get_n_activity, s.get_eff_activity):
        full = get(max_n=1 << 20)
        assert len(full['hash']) > 20
        assert len(get(max_n=5)['hash']) == len(full['hash'])


def test_crossings_are_not_capped_by_the_buffer(diverse):
    full = diverse.activity_crossings(max_n=1 << 22)
    assert len(full) > 20
    assert len(diverse.activity_crossings(max_n=3)) == len(full)


def test_flux_is_the_uncapped_sum_of_the_crossings(diverse):
    for window in (1, 4, 8):
        c = diverse.activity_crossings(window=window, max_n=1 << 22)
        assert diverse.activity_flux(window=window) == pytest.approx(c.sum() / window, rel=1e-6)
        assert c.sum() > 0


def test_flux_beyond_the_old_65536_cap(make_sim):
    """The N7 failure: more samples than the old buffer held. The old flux
    read at most 65536 / 8 = 8192 here."""
    s = make_sim(256, mu_lut=0.02, food_inc=0.02, tax=0.01)
    s.set_seed(1)
    s.set_lut_random(2)
    s.set_egenome_random()
    s.set_alive_all()
    s.set_f_all(0.8)
    for _ in range(50):
        s.step()
        s._lib.evoca_activity_update()
    c = s.activity_crossings()
    assert len(c) > 2 * 65536
    assert s.activity_flux() == pytest.approx(c.sum() / 8, rel=1e-6)
    assert s.activity_flux() > 2 * 8192
