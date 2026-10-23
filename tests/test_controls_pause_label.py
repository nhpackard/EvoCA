"""The Run/Pause button label must match the sim's pause state.

Regression (N, 2026-10-23): Gallery.show opened running with the button
reading "Run". The front end snaps a slider value to its step
(m_scale 3.5 -> 3.5000000000000004) and echoes it back; the slider
callback took that as a move, paused (label "Run"), and its 200 ms timer
resumed the sim without updating the label.

run_with_controls is driven headless: the SDL subprocess and the notebook
display are stubbed, and the widgets are read from what it displays."""
import io
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
pytest.importorskip("ipywidgets")

import controls                    # noqa: E402
from evoca_py import EvoCA         # noqa: E402


class _FakeProc:
    stdout = io.StringIO("")
    stdin = None
    pid = 0
    returncode = None

    def poll(self):
        return None

    def wait(self, *a, **k):
        return 0

    def terminate(self):
        pass

    def kill(self):
        pass


def _walk(w):
    yield w
    for c in getattr(w, "children", ()):
        yield from _walk(c)


@pytest.fixture
def session(monkeypatch):
    shown = []
    monkeypatch.setattr(controls.subprocess, "Popen", lambda *a, **k: _FakeProc())
    monkeypatch.setattr(controls, "ipy_display", lambda w: shown.append(w))
    sim = EvoCA()
    sim.init(32, food_inc=0.016, m_scale=3.5, tax=0.035)
    sim.state(lut="gol")
    sim.set_alive_halfplane(0)
    controls.run_with_controls(sim, cell_px=1, paused=False)
    widgets = [w for x in shown for w in _walk(x)]
    btn = next(w for w in widgets if type(w).__name__ == "ToggleButton")
    sl_m = next(w for w in widgets
                if type(w).__name__ == "FloatSlider" and w.value == pytest.approx(3.5))
    yield sim, btn, sl_m
    controls._active_stop()


def _settle():
    time.sleep(controls._SLIDER_RESUME_S + 0.2)


def test_starts_running_with_pause_label(session):
    sim, btn, _ = session
    assert (btn.description, btn.value) == ("Pause", False)
    t0 = sim.get_step(); time.sleep(0.2)
    assert sim.get_step() > t0


def test_step_snap_echo_is_not_a_slider_move(session):
    sim, btn, sl_m = session
    sl_m.value = 3.5000000000000004          # the front end's snapped echo
    assert (btn.description, btn.value) == ("Pause", False)
    assert sim.m_scale == 3.5                # the sim keeps the exact value
    _settle()
    assert (btn.description, btn.value) == ("Pause", False)


def test_real_slider_move_resumes_with_pause_label(session):
    sim, btn, sl_m = session
    sl_m.value = 3.0
    assert sim.m_scale == pytest.approx(3.0)
    _settle()
    assert (btn.description, btn.value) == ("Pause", False)
    t0 = sim.get_step(); time.sleep(0.2)
    assert sim.get_step() > t0


# ── Food/waste sliders (n_chem=2 only) ────────────────────────────────

WASTE = ("waste_sink", "wdiff", "waste_disp", "waste_tox")


def _panel(monkeypatch, **init):
    shown = []
    monkeypatch.setattr(controls.subprocess, "Popen", lambda *a, **k: _FakeProc())
    monkeypatch.setattr(controls, "ipy_display", lambda w: shown.append(w))
    sim = EvoCA()
    sim.init(32, food_inc=0.016, m_scale=1.0, tax=0.035, **init)
    sim.state(lut="gol")
    controls.run_with_controls(sim, cell_px=1, paused=True)
    sliders = {w.description.rstrip(":"): w for x in shown for w in _walk(x)
               if type(w).__name__ == "FloatSlider"}
    return sim, sliders


def test_waste_sliders_absent_for_single_food(monkeypatch):
    sim, sliders = _panel(monkeypatch)
    try:
        assert not set(WASTE) & set(sliders)
    finally:
        controls._active_stop()


def test_waste_sliders_reach_the_sim(monkeypatch):
    sim, sliders = _panel(monkeypatch, n_chem=2, wdiff=0.01, waste_sink=0.015,
                          waste_disp=0.5, waste_tox=1.0)
    try:
        assert set(WASTE) <= set(sliders)
        assert sliders["waste_tox"].min > 0
        new = dict(waste_sink=0.02, wdiff=0.05, waste_disp=0.3, waste_tox=0.4)
        for k, v in new.items():
            sliders[k].value = v
        for k, v in new.items():
            assert getattr(sim, k) == pytest.approx(v, rel=1e-6), k
        inc, sink, diff, disp, tox = sim.get_chem_params(EvoCA.WASTE)
        assert (sink, diff, disp, tox) == pytest.approx((0.02, 0.05, 0.3, 0.4), rel=1e-6)
    finally:
        controls._active_stop()
