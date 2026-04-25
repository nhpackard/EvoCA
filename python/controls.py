"""
controls.py — ipywidgets control panel + SDL2 display for EvoCA.

Architecture
------------
SDL2 requires macOS's actual main thread (thread 0), which belongs to the
Jupyter kernel's event loop.  Running SDL2 in any other thread crashes the
kernel.  The fix: run SDL2 in a *subprocess* that has its own main thread.

  Main process (Jupyter kernel)
  ├── Jupyter event loop  →  ipywidgets callbacks fire here
  ├── Sim thread          →  sim.step() + sim.colorize() → shared memory
  ├── Reader thread       →  relays SDL worker stdout to terminal
  └── subprocess (sdl_worker.py)
      └── SDL2 main thread  →  reads shared memory, renders window

Pixel data flows via POSIX shared memory with no copying:
  sim.colorize(pixels_shm, …) writes directly into shared memory;
  the subprocess reads the same buffer.

Usage (Jupyter cell)
--------------------
    from python.controls import run_with_controls
    run_with_controls(sim_sdl)
    # Cell returns immediately; widgets appear below; SDL2 window opens.
    # Click Quit or press Q/Esc in the SDL2 window to stop.
    # Then call sim_sdl.free() in the next cell.
"""

import atexit
import os
import sys
import subprocess
import threading
import time

import numpy as np
import ipywidgets as widgets
import matplotlib.pyplot as plt
from IPython.display import display as ipy_display
from multiprocessing.shared_memory import SharedMemory

import ctypes

COLOR_MODES      = ["state", "env-food", "priv-food", "births", "age"]
_SLIDER_RESUME_S = 0.20   # seconds after last slider touch before auto-resume
_WORKER          = os.path.join(os.path.dirname(__file__), "sdl_worker.py")

# ctrl_shm layout (5 × int32)
_QUIT, _CMODE, _STEP, _FPS10, _PAUSED = 0, 1, 2, 3, 4

# Probe strip-chart constants
PROBE_W = 512    # pixels = time steps visible
PROBE_H = 128    # pixel height of each strip chart

# Grouped time-series probe: 6 traces in 2 stacked strips (3 each).
# Order below must match _TS_COLORS / _TS_LABELS in sdl_worker.py.
_TS_TRACES = ('pop', 'F_mean', 'f_mean',
              'lut_div', 'eg_ent', 'activity_flux')
TS_N_TRACES = len(_TS_TRACES)
_TS_GROUPS  = ((0, 1, 2), (3, 4, 5))
TS_STRIPS   = len(_TS_GROUPS)

# Map probe name → (C getter function name, ctype element type)
_PROBE_GETTER = {
    'env_food':  ('evoca_get_F',      ctypes.c_float),
    'priv_food': ('evoca_get_f',      ctypes.c_float),
    'births':    ('evoca_get_births', ctypes.c_uint8),
}

_AVAILABLE_PROBES = {
    'env_food':       'Mean +/- std of environmental food F(x)',
    'priv_food':      'Mean +/- std of private food f(x)',
    'births':         'Mean +/- std of births array',
    'activity':       'LUT genome activity (scrolling hash-colored strip)',
    'eg_activity':    'Egenome activity (scrolling hash-colored strip)',
    'lut_complexity': 'Stacked area: LUT ring-dependency level',
    'eg_pop':         'Stacked area: egenome population fractions',
    'entropy':        'Local-pattern Shannon entropy',
    'pat_activity':   'Local-pattern activity (scrolling hash-colored strip)',
    'q_activity':     'Activity quantile profile (decile strip chart)',
    'ts':             ('Grouped time-series: pop, F_mean, f_mean '
                       '(top) / lut_div, eg_ent, activity_flux (bottom)'),
}


def available_probes():
    """Return a dict mapping probe names to short descriptions."""
    return dict(_AVAILABLE_PROBES)


# Module-level handle: stop any previous session before starting a new one.
_active_stop = None


def run_with_controls(sim, cell_px=None, colormode=0, paused=True, probes=None,
                      diag=False):
    """
    Display ipywidgets controls and open an SDL2 simulation window.

    Returns immediately (non-blocking).  The simulation runs in a background
    thread; the SDL2 display runs in a subprocess.

    Parameters
    ----------
    sim       : initialised EvoCA instance
    cell_px   : screen pixels per cell (default: sim.cell_px from CELL_PX #define)
    colormode : initial colour mode (0=state, 1=env-food, 2=priv-food)
    paused    : if True, start in paused state
    probes    : dict of probe names to enable, e.g. {'env_food': True, 'priv_food': True}

    Returns
    -------
    threading.Thread — the simulation thread (can be .join()-ed if desired)
    """
    # ── Tear down any previous session ────────────────────────────
    global _active_stop
    if _active_stop is not None:
        _active_stop()
        _active_stop = None

    N  = sim.N
    px = cell_px if cell_px is not None else sim.cell_px

    # ── Diagnostics ─────────────────────────────────────────────────
    sim._lib.evoca_set_diag(int(diag))

    # ── Shared memory ─────────────────────────────────────────────
    pixel_shm = SharedMemory(create=True, size=N * N * 4)
    ctrl_shm  = SharedMemory(create=True, size=5 * 4)   # 5 int32

    # numpy views into shared memory
    pixels = np.ndarray((N * N,), dtype=np.int32, buffer=pixel_shm.buf)
    ctrl   = np.ndarray((5,),     dtype=np.int32, buffer=ctrl_shm.buf)
    ctrl[:] = [0, colormode, 0, 0, int(paused)]

    # ── Activity probe setup ────────────────────────────────────────
    activity_enabled = bool((probes or {}).get('activity'))
    ACT_H = 2 * PROBE_H  # 256 pixels tall
    activity_shm     = None
    activity_cursor  = None   # int32 view
    activity_pixels  = None   # int32[ACT_H, PROBE_W] view
    activity_col     = None   # temp column buffer

    if activity_enabled:
        act_shm_size = 4 + PROBE_W * ACT_H * 4
        activity_shm = SharedMemory(create=True, size=act_shm_size)
        _abuf = np.ndarray((act_shm_size,), dtype=np.uint8,
                           buffer=activity_shm.buf)
        _abuf[:] = 0
        activity_cursor = np.ndarray((1,), dtype=np.int32,
                                     buffer=activity_shm.buf)
        activity_pixels = np.ndarray((ACT_H, PROBE_W), dtype=np.int32,
                                     buffer=activity_shm.buf, offset=4)
        activity_col = np.zeros(ACT_H, dtype=np.int32)

    # ── Egenome activity probe setup ─────────────────────────────────
    eg_activity_enabled = bool((probes or {}).get('eg_activity'))
    eg_activity_shm     = None
    eg_activity_cursor  = None
    eg_activity_pixels  = None
    eg_activity_col     = None

    if eg_activity_enabled:
        ega_meta_off = 4 + PROBE_W * ACT_H * 4
        ega_shm_size = ega_meta_off + 64*8 + 64*4 + 64*4 + 4  # +acts,pops,cols,ymax
        eg_activity_shm = SharedMemory(create=True, size=ega_shm_size)
        _egabuf = np.ndarray((ega_shm_size,), dtype=np.uint8,
                             buffer=eg_activity_shm.buf)
        _egabuf[:] = 0
        eg_activity_cursor = np.ndarray((1,), dtype=np.int32,
                                        buffer=eg_activity_shm.buf)
        eg_activity_pixels = np.ndarray((ACT_H, PROBE_W), dtype=np.int32,
                                        buffer=eg_activity_shm.buf, offset=4)
        eg_activity_col = np.zeros(ACT_H, dtype=np.int32)
        # Metadata for click → egenome mapping (written by sim thread)
        ega_m_acts = np.ndarray((64,), dtype=np.uint64,
                                buffer=eg_activity_shm.buf, offset=ega_meta_off)
        ega_m_pops = np.ndarray((64,), dtype=np.uint32,
                                buffer=eg_activity_shm.buf, offset=ega_meta_off + 64*8)
        ega_m_cols = np.ndarray((64,), dtype=np.int32,
                                buffer=eg_activity_shm.buf, offset=ega_meta_off + 64*8 + 64*4)
        ega_m_ymax = np.ndarray((1,), dtype=np.int32,
                                buffer=eg_activity_shm.buf, offset=ega_meta_off + 64*8 + 64*4*2)

    # ── Entropy probe setup ──────────────────────────────────────────
    entropy_enabled      = bool((probes or {}).get('entropy'))
    entropy_shm          = None
    entropy_cursor       = None   # int32[1] view
    entropy_buf          = None   # float32[PROBE_W] view

    if entropy_enabled:
        ent_shm_size = 4 + PROBE_W * 4   # cursor + float32[PROBE_W]
        entropy_shm = SharedMemory(create=True, size=ent_shm_size)
        _ebuf = np.ndarray((ent_shm_size,), dtype=np.uint8,
                           buffer=entropy_shm.buf)
        _ebuf[:] = 0
        entropy_cursor = np.ndarray((1,), dtype=np.int32,
                                    buffer=entropy_shm.buf)
        entropy_buf = np.ndarray((PROBE_W,), dtype=np.float32,
                                 buffer=entropy_shm.buf, offset=4)

    # ── Pattern activity probe setup ─────────────────────────────────
    pat_activity_enabled = bool((probes or {}).get('pat_activity'))
    PAT_H                = ACT_H   # 256 px tall (same as activity window)
    pat_activity_shm     = None
    pat_activity_cursor  = None
    pat_activity_pixels  = None
    pat_activity_col     = None

    if pat_activity_enabled:
        pa_shm_size = 4 + PROBE_W * PAT_H * 4
        pat_activity_shm = SharedMemory(create=True, size=pa_shm_size)
        _pabuf = np.ndarray((pa_shm_size,), dtype=np.uint8,
                            buffer=pat_activity_shm.buf)
        _pabuf[:] = 0
        pat_activity_cursor = np.ndarray((1,), dtype=np.int32,
                                         buffer=pat_activity_shm.buf)
        pat_activity_pixels = np.ndarray((PAT_H, PROBE_W), dtype=np.int32,
                                         buffer=pat_activity_shm.buf, offset=4)
        pat_activity_col = np.zeros(PAT_H, dtype=np.int32)

    pat_enabled = entropy_enabled or pat_activity_enabled

    # ── Activity quantile (q_activity) probe setup ──────────────────
    q_activity_enabled  = bool((probes or {}).get('q_activity'))
    QA_N_DECILES        = 9
    q_activity_shm      = None
    q_activity_cursor   = None
    q_activity_deciles  = None   # list of 9 float32[PROBE_W] views
    q_activity_col      = None   # temp float32[9]

    if q_activity_enabled:
        qa_shm_size = 4 + QA_N_DECILES * PROBE_W * 4
        q_activity_shm = SharedMemory(create=True, size=qa_shm_size)
        _qabuf = np.ndarray((qa_shm_size,), dtype=np.uint8,
                             buffer=q_activity_shm.buf)
        _qabuf[:] = 0
        q_activity_cursor = np.ndarray((1,), dtype=np.int32,
                                        buffer=q_activity_shm.buf)
        q_activity_deciles = []
        off = 4
        for _ in range(QA_N_DECILES):
            q_activity_deciles.append(
                np.ndarray((PROBE_W,), dtype=np.float32,
                           buffer=q_activity_shm.buf, offset=off))
            off += PROBE_W * 4
        q_activity_col = np.zeros(QA_N_DECILES, dtype=np.float32)

    # ── Set n_ent on C library ───────────────────────────────────────
    if pat_enabled:
        sim._lib.evoca_set_n_ent(sim.n_ent)

    # ── LUT complexity probe setup ───────────────────────────────────
    lut_complexity_enabled = bool((probes or {}).get('lut_complexity'))
    lut_complexity_shm     = None
    lut_complexity_cursor  = None
    lut_complexity_pixels  = None
    lut_complexity_col     = None

    if lut_complexity_enabled:
        lc_shm_size = 4 + PROBE_W * PROBE_H * 4
        lut_complexity_shm = SharedMemory(create=True, size=lc_shm_size)
        _lcbuf = np.ndarray((lc_shm_size,), dtype=np.uint8,
                            buffer=lut_complexity_shm.buf)
        _lcbuf[:] = 0
        lut_complexity_cursor = np.ndarray((1,), dtype=np.int32,
                                           buffer=lut_complexity_shm.buf)
        lut_complexity_pixels = np.ndarray((PROBE_H, PROBE_W), dtype=np.int32,
                                           buffer=lut_complexity_shm.buf, offset=4)
        lut_complexity_col = np.zeros(PROBE_H, dtype=np.int32)

    # ── Egenome population banded chart setup ───────────────────────
    eg_pop_enabled = bool((probes or {}).get('eg_pop'))
    eg_pop_shm     = None
    eg_pop_cursor  = None
    eg_pop_pixels  = None
    eg_pop_col     = None

    if eg_pop_enabled:
        ep_shm_size = 4 + PROBE_W * PROBE_H * 4
        eg_pop_shm = SharedMemory(create=True, size=ep_shm_size)
        _epbuf = np.ndarray((ep_shm_size,), dtype=np.uint8,
                            buffer=eg_pop_shm.buf)
        _epbuf[:] = 0
        eg_pop_cursor = np.ndarray((1,), dtype=np.int32,
                                    buffer=eg_pop_shm.buf)
        eg_pop_pixels = np.ndarray((PROBE_H, PROBE_W), dtype=np.int32,
                                    buffer=eg_pop_shm.buf, offset=4)
        eg_pop_col = np.zeros(PROBE_H, dtype=np.int32)

    # ── Grouped time-series (ts) probe setup ────────────────────────
    ts_enabled = bool((probes or {}).get('ts'))
    ts_shm     = None
    ts_cursor  = None
    ts_bufs    = None   # list of 6 float32[PROBE_W]
    ts_getters = None   # (F_ptr_fn, f_ptr_fn, alive_ptr_fn)

    if ts_enabled:
        ts_shm_size = 4 + TS_N_TRACES * PROBE_W * 4
        ts_shm = SharedMemory(create=True, size=ts_shm_size)
        _tsbuf = np.ndarray((ts_shm_size,), dtype=np.uint8,
                            buffer=ts_shm.buf)
        _tsbuf[:] = 0
        ts_cursor = np.ndarray((1,), dtype=np.int32, buffer=ts_shm.buf)
        ts_bufs = []
        off = 4
        for _ in range(TS_N_TRACES):
            ts_bufs.append(np.ndarray((PROBE_W,), dtype=np.float32,
                                       buffer=ts_shm.buf, offset=off))
            off += PROBE_W * 4
        ts_getters = (sim._lib.evoca_get_F,
                       sim._lib.evoca_get_f,
                       sim._lib.evoca_get_alive)

    # ── Probe setup ─────────────────────────────────────────────────
    probe_names = [k for k, v in (probes or {}).items() if v and k in _PROBE_GETTER]
    n_probes    = len(probe_names)
    probe_shm   = None
    probe_cursor = None       # int32 view into probe_shm
    probe_means  = []         # list of float32[PROBE_W] views
    probe_stds   = []         # list of float32[PROBE_W] views

    # Build list of (C getter func, numpy dtype) for fast access in sim thread
    probe_getters = []   # list of (getter_fn, np.dtype)
    for pname in probe_names:
        fn_name, ctype = _PROBE_GETTER[pname]
        getter = getattr(sim._lib, fn_name)
        getter.argtypes = []
        getter.restype  = ctypes.POINTER(ctype)
        probe_getters.append((getter, np.dtype(ctype)))

    if n_probes > 0:
        probe_shm_size = 4 + n_probes * 2 * PROBE_W * 4
        probe_shm = SharedMemory(create=True, size=probe_shm_size)
        _pbuf = np.ndarray((probe_shm_size,), dtype=np.uint8, buffer=probe_shm.buf)
        _pbuf[:] = 0
        probe_cursor = np.ndarray((1,), dtype=np.int32, buffer=probe_shm.buf)
        off = 4
        for _ in range(n_probes):
            probe_means.append(np.ndarray((PROBE_W,), dtype=np.float32,
                                          buffer=probe_shm.buf, offset=off))
            off += PROBE_W * 4
            probe_stds.append(np.ndarray((PROBE_W,), dtype=np.float32,
                                         buffer=probe_shm.buf, offset=off))
            off += PROBE_W * 4

    # ── SDL2 subprocess ───────────────────────────────────────────
    cmd = [sys.executable, _WORKER,
           pixel_shm.name, ctrl_shm.name, str(N), str(px)]
    if n_probes > 0:
        cmd += [probe_shm.name, ",".join(probe_names)]
    if activity_enabled:
        cmd += ["--activity=" + activity_shm.name]
    if eg_activity_enabled:
        cmd += ["--eg-activity=" + eg_activity_shm.name]
    if lut_complexity_enabled:
        cmd += ["--lut-complexity=" + lut_complexity_shm.name]
    if entropy_enabled:
        cmd += ["--entropy=" + entropy_shm.name]
    if pat_activity_enabled:
        cmd += ["--pat-activity=" + pat_activity_shm.name]
    if eg_pop_enabled:
        cmd += ["--eg-pop=" + eg_pop_shm.name]
    if q_activity_enabled:
        cmd += ["--q-activity=" + q_activity_shm.name]
    if ts_enabled:
        cmd += ["--ts=" + ts_shm.name]
    sdl_proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )

    # Relay subprocess output to the terminal (for diagnostics)
    def _reader():
        for line in sdl_proc.stdout:
            print(line.decode(errors='replace'), end='', flush=True)

    threading.Thread(target=_reader, name="evoca-sdl-reader", daemon=True).start()

    # ── ts record helper ────────────────────────────────────────────
    def _record_ts():
        if not ts_enabled:
            return
        F_ptr = ts_getters[0]()
        F_arr = np.ctypeslib.as_array(F_ptr, shape=(N * N,))
        f_ptr = ts_getters[1]()
        f_arr = np.ctypeslib.as_array(f_ptr, shape=(N * N,))
        a_ptr = ts_getters[2]()
        a_arr = np.ctypeslib.as_array(a_ptr, shape=(N * N,))
        pop = int(a_arr.sum())
        F_mean = float(F_arr.mean())
        if pop > 0:
            mask = a_arr.astype(bool)
            f_mean = float(f_arr[mask].mean())
        else:
            f_mean = 0.0
        # Activity must be updated at least once per tick for diversity
        # and flux to be current.
        if not (activity_enabled or q_activity_enabled):
            sim._lib.evoca_activity_update()
        lut_div = sim.n_distinct_genomes()
        # Egenome entropy: Shannon entropy (bits) over the 64-bucket egenome
        # distribution. Unlike a raw distinct-count, this stays responsive
        # to mutation pressure even after all 64 buckets are occupied.
        # Needs a recent eg update.
        if not (eg_activity_enabled or eg_pop_enabled):
            sim._lib.evoca_eg_activity_update()
        eg_info = sim.get_eg_activity()
        eg_pops = eg_info['pop_count']
        total = int(eg_pops.sum())
        if total > 0:
            p = eg_pops[eg_pops > 0].astype(np.float64) / total
            eg_ent = float(-(p * np.log2(p)).sum())
        else:
            eg_ent = 0.0
        flux = sim.activity_flux()
        ts_cur = int(ts_cursor[0])
        vals = (float(pop), F_mean, f_mean,
                float(lut_div), eg_ent, float(flux))
        for ti, v in enumerate(vals):
            ts_bufs[ti][ts_cur] = v
        ts_cursor[0] = (ts_cur + 1) % PROBE_W

    # ── Shared sim state ──────────────────────────────────────────
    st    = dict(paused=bool(paused), running=True, colormode=colormode, step_cnt=0)
    # _alive guards widget callbacks against use-after-free:
    # set to False just before shared memory is released.
    _alive        = [True]
    _cleanup_done = [False]

    # ── Guaranteed cleanup (sim thread or atexit, whichever fires first) ──
    def _do_cleanup():
        if _cleanup_done[0]:
            return
        _cleanup_done[0] = True
        _alive[0] = False
        if sdl_proc.poll() is None:
            sdl_proc.terminate()
            try:
                sdl_proc.wait(timeout=2)
            except Exception:
                pass
        all_shm = [pixel_shm, ctrl_shm]
        if probe_shm is not None:
            all_shm.append(probe_shm)
        if activity_shm is not None:
            all_shm.append(activity_shm)
        if eg_activity_shm is not None:
            all_shm.append(eg_activity_shm)
        if lut_complexity_shm is not None:
            all_shm.append(lut_complexity_shm)
        if entropy_shm is not None:
            all_shm.append(entropy_shm)
        if pat_activity_shm is not None:
            all_shm.append(pat_activity_shm)
        if eg_pop_shm is not None:
            all_shm.append(eg_pop_shm)
        if q_activity_shm is not None:
            all_shm.append(q_activity_shm)
        if ts_shm is not None:
            all_shm.append(ts_shm)
        for shm in all_shm:
            try:
                shm.unlink()
            except Exception:
                pass
            try:
                shm.close()
            except Exception:
                pass

    atexit.register(_do_cleanup)

    def _stop():
        """Signal sim thread to quit and wait for cleanup."""
        global _active_stop
        st['running'] = False
        if _alive[0]:
            ctrl[_QUIT] = 1
        # Give sim thread time to finish and call _do_cleanup
        time.sleep(0.1)
        _do_cleanup()
        _active_stop = None

    _active_stop = _stop
    sim._stop_display = _stop

    # ── ipywidgets ─────────────────────────────────────────────────
    btn_pause = widgets.ToggleButton(
        value=bool(paused), description="Run" if paused else "Pause",
        button_style="", layout=widgets.Layout(width="90px"))
    btn_restart = widgets.Button(
        description="Restart", layout=widgets.Layout(width="80px"))
    btn_step  = widgets.Button(
        description="Step",  layout=widgets.Layout(width="70px"))
    btn_step200 = widgets.Button(
        description="Step200", layout=widgets.Layout(width="80px"))
    btn_quit  = widgets.Button(
        description="Quit",  button_style="danger",
        layout=widgets.Layout(width="70px"))
    btn_save  = widgets.Button(
        description="Save Plots", layout=widgets.Layout(width="100px"))
    txt_descriptor = widgets.Text(
        placeholder='run descriptor', layout=widgets.Layout(width="200px"))
    btn_export = widgets.Button(
        description="Export", layout=widgets.Layout(width="70px"))

    sl_kw = dict(continuous_update=True,
                 style={"description_width": "90px"},
                 layout=widgets.Layout(width="440px"))

    # Slider limits scale to the initialized parameter value: max = 2 × init,
    # with a tiny floor so a zero-init slider isn't completely degenerate.
    # Step gives ~100 increments, no finer than the floor.
    def _flt_lims(init, floor):
        mx = max(2.0 * float(init), floor)
        st = max(mx / 100.0, floor)
        return mx, st

    def _int_lims(init, floor):
        return max(2 * int(init), floor)

    fi_max, fi_step = _flt_lims(sim.food_inc,   1e-3)
    sl_food_inc = widgets.FloatSlider(
        value=sim.food_inc,   min=0.0, max=fi_max, step=fi_step,
        description="food_inc:",   readout_format=".4f", **sl_kw)
    ms_max, ms_step = _flt_lims(sim.m_scale,    1e-2)
    sl_m_scale = widgets.FloatSlider(
        value=sim.m_scale,    min=0.0, max=ms_max, step=ms_step,
        description="m_scale:",    readout_format=".3f", **sl_kw)
    gd_max = _int_lims(sim.gdiff, 2)
    sl_gdiff = widgets.IntSlider(
        value=sim.gdiff, min=0, max=gd_max, step=1,
        description="gdiff:",
        style={"description_width": "90px"},
        layout=widgets.Layout(width="440px"))
    mul_max, mul_step = _flt_lims(sim.mu_lut,     1e-6)
    sl_mu_lut = widgets.FloatSlider(
        value=sim.mu_lut, min=0.0, max=mul_max, step=mul_step,
        description="mu_lut:",    readout_format=".6f", **sl_kw)
    mue_max, mue_step = _flt_lims(sim.mu_egenome, 1e-4)
    sl_mu_egenome = widgets.FloatSlider(
        value=sim.mu_egenome, min=0.0, max=mue_max, step=mue_step,
        description="mu_egenome:", readout_format=".4f", **sl_kw)
    tx_max, tx_step = _flt_lims(sim.tax,        1e-4)
    sl_tax = widgets.FloatSlider(
        value=sim.tax, min=0.0, max=tx_max, step=tx_step,
        description="tax:", readout_format=".4f", **sl_kw)

    # ── ymax halve / double buttons ──────────────────────────────────
    def _make_ymax_btns(name, initial, update_fn):
        val = [initial]
        lbl = widgets.Label(value=f"{name}: {initial}",
                            layout=widgets.Layout(width="160px"))
        btn_h = widgets.Button(description="<|",
                               layout=widgets.Layout(width="35px"))
        btn_d = widgets.Button(description="|>",
                               layout=widgets.Layout(width="35px"))
        def halve(_):
            if not _alive[0]: return
            val[0] = max(val[0] // 2, 100)
            update_fn(val[0])
            lbl.value = f"{name}: {val[0]}"
        def double(_):
            if not _alive[0]: return
            val[0] *= 2
            update_fn(val[0])
            lbl.value = f"{name}: {val[0]}"
        btn_h.on_click(halve)
        btn_d.on_click(double)
        return widgets.HBox([btn_h, lbl, btn_d])

    _ymax_btns = []
    if activity_enabled:
        _ymax_btns.append(_make_ymax_btns(
            "act_ymax", 2000, sim.update_act_ymax))
    if eg_activity_enabled:
        _ymax_btns.append(_make_ymax_btns(
            "eg_act_ymax", 2000, sim.update_eg_act_ymax))
    if pat_activity_enabled:
        _ymax_btns.append(_make_ymax_btns(
            "pat_act_ymax", 2000, sim.update_pat_act_ymax))

    cb_restricted_mu = widgets.Checkbox(
        value=bool(sim.restricted_mu),
        description="restricted_mu",
        layout=widgets.Layout(width="200px"))
    color_dd   = widgets.Dropdown(
        options=COLOR_MODES, value=COLOR_MODES[colormode],
        description="Color:",
        layout=widgets.Layout(width="200px"))
    status_lbl = widgets.Label(value="Starting…")

    _rows = [
        widgets.HBox([btn_pause, btn_restart, btn_step, btn_step200,
                      btn_quit, btn_save, txt_descriptor, btn_export]),
        sl_food_inc, sl_m_scale, sl_gdiff,
        sl_mu_lut, sl_mu_egenome, sl_tax,
    ]
    if _ymax_btns:
        _rows.append(widgets.HBox(_ymax_btns))
    _rows.append(widgets.HBox([color_dd, cb_restricted_mu, status_lbl]))
    ipy_display(widgets.VBox(_rows))

    # ── Widget callbacks ──────────────────────────────────────────
    _guard = [False]

    def _set_paused(p):
        if _guard[0] or not _alive[0]:
            return
        _guard[0] = True
        st['paused']          = p
        ctrl[_PAUSED]         = int(p)
        btn_pause.value       = p
        btn_pause.description = "Run" if p else "Pause"
        _guard[0] = False

    def on_pause_toggle(change):
        _set_paused(change['new'])

    def _record_probes():
        if n_probes > 0:
            cur = int(probe_cursor[0])
            for pi, (getter, dt) in enumerate(probe_getters):
                ptr = getter()
                arr = np.ctypeslib.as_array(ptr, shape=(N * N,))
                farr = arr.astype(np.float64) if dt != np.float32 else arr
                probe_means[pi][cur] = farr.mean()
                probe_stds[pi][cur]  = farr.std()
            probe_cursor[0] = (cur + 1) % PROBE_W
        if activity_enabled:
            sim._lib.evoca_activity_update()
            act_cur = int(activity_cursor[0])
            act_col_ptr = activity_col.ctypes.data_as(
                ctypes.POINTER(ctypes.c_int32))
            sim._lib.evoca_activity_render_col(act_col_ptr, ACT_H)
            activity_pixels[:, act_cur] = activity_col
            activity_cursor[0] = (act_cur + 1) % PROBE_W
        if q_activity_enabled:
            if not activity_enabled:
                sim._lib.evoca_activity_update()
            qa_cur = int(q_activity_cursor[0])
            qa_col_ptr = q_activity_col.ctypes.data_as(
                ctypes.POINTER(ctypes.c_float))
            sim._lib.evoca_q_activity_deciles(qa_col_ptr)
            for di in range(QA_N_DECILES):
                q_activity_deciles[di][qa_cur] = q_activity_col[di]
            q_activity_cursor[0] = (qa_cur + 1) % PROBE_W
        if eg_activity_enabled or eg_pop_enabled:
            sim._lib.evoca_eg_activity_update()
        if eg_activity_enabled:
            ega_cur = int(eg_activity_cursor[0])
            ega_col_ptr = eg_activity_col.ctypes.data_as(
                ctypes.POINTER(ctypes.c_int32))
            sim._lib.evoca_eg_activity_render_col(ega_col_ptr, ACT_H)
            eg_activity_pixels[:, ega_cur] = eg_activity_col
            eg_activity_cursor[0] = (ega_cur + 1) % PROBE_W
            sim._lib.evoca_eg_activity_get(
                ega_m_acts.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
                ega_m_pops.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                ega_m_cols.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)))
            ega_m_ymax[0] = sim._lib.evoca_get_eg_act_ymax()
        if lut_complexity_enabled:
            lc_cur = int(lut_complexity_cursor[0])
            lc_col_ptr = lut_complexity_col.ctypes.data_as(
                ctypes.POINTER(ctypes.c_int32))
            sim._lib.evoca_lut_complexity_render_col(lc_col_ptr, PROBE_H)
            lut_complexity_pixels[:, lc_cur] = lut_complexity_col
            lut_complexity_cursor[0] = (lc_cur + 1) % PROBE_W
        if eg_pop_enabled:
            ep_cur = int(eg_pop_cursor[0])
            ep_col_ptr = eg_pop_col.ctypes.data_as(
                ctypes.POINTER(ctypes.c_int32))
            sim._lib.evoca_eg_pop_render_col(ep_col_ptr, PROBE_H)
            eg_pop_pixels[:, ep_cur] = eg_pop_col
            eg_pop_cursor[0] = (ep_cur + 1) % PROBE_W
        if pat_enabled:
            sim._lib.evoca_pat_update()
            if entropy_enabled:
                ec = int(entropy_cursor[0])
                entropy_buf[ec] = sim._lib.evoca_get_entropy()
                entropy_cursor[0] = (ec + 1) % PROBE_W
            if pat_activity_enabled:
                pac = int(pat_activity_cursor[0])
                pa_col_ptr = pat_activity_col.ctypes.data_as(
                    ctypes.POINTER(ctypes.c_int32))
                sim._lib.evoca_pat_activity_render_col(pa_col_ptr, PAT_H)
                pat_activity_pixels[:, pac] = pat_activity_col
                pat_activity_cursor[0] = (pac + 1) % PROBE_W
        _record_ts()

    def on_step(_):
        if not _alive[0]:
            return
        if st['paused']:
            sim.step()
            st['step_cnt'] += 1
            ctrl[_STEP] = st['step_cnt']
            sim.colorize(pixels, st['colormode'])
            _record_probes()
            status_lbl.value = f"t={st['step_cnt']}  (paused)"

    def _step_display(n=1, delay=0.0):
        """Programmatic step loop that advances the SDL window, probes,
        and step counter. Requires the sim to be paused (otherwise the
        background sim thread would race on C state)."""
        if not _alive[0]:
            raise RuntimeError("no display session attached")
        if not st['paused']:
            raise RuntimeError(
                "sim.step_display() requires the sim paused "
                "(click Pause, or set btn_pause.value=True)")
        for _ in range(int(n)):
            sim.step()
            st['step_cnt'] += 1
            ctrl[_STEP] = st['step_cnt']
            sim.colorize(pixels, st['colormode'])
            _record_probes()
            status_lbl.value = f"t={st['step_cnt']}  (stepping)"
            if delay > 0:
                time.sleep(delay)
        status_lbl.value = f"t={st['step_cnt']}  (paused)"

    sim.step_display = _step_display

    def on_step200(_):
        if not _alive[0]:
            return
        if not st['paused']:
            _set_paused(True)
        _step_display(200, delay=0.02)

    def on_quit(_):
        st['running'] = False
        if _alive[0]:
            ctrl[_QUIT] = 1
        status_lbl.value = "Stopped — call sim.free() when ready"

    def on_color(change):
        if not _alive[0]:
            return
        cm = COLOR_MODES.index(change['new'])
        st['colormode'] = cm
        ctrl[_CMODE]    = cm

    def on_save(_):
        if not _alive[0] or n_probes == 0:
            return
        cur = int(probe_cursor[0])
        for i, pname in enumerate(probe_names):
            # Reconstruct time-ordered arrays from ring buffer
            m = np.roll(probe_means[i], -cur)
            s = np.roll(probe_stds[i], -cur)
            t = np.arange(PROBE_W)
            fig_s, ax_s = plt.subplots(figsize=(8, 2.5))
            ax_s.fill_between(t, m - s, m + s, alpha=0.3)
            ax_s.plot(t, m, linewidth=0.8)
            ax_s.set_title(pname)
            ax_s.set_xlabel("t (relative)")
            fig_s.tight_layout()
            fname = f"probe_{pname}.png"
            fig_s.savefig(fname, dpi=150)
            plt.close(fig_s)
        status_lbl.value = f"Saved {n_probes} probe plot(s)"

    def on_export(_):
        desc = txt_descriptor.value.strip() or 'unnamed'
        cm = COLOR_MODES.index(color_dd.value) if color_dd.value in COLOR_MODES else 0
        path = sim.export_recipe(desc, probes=probes, colormode=cm)
        status_lbl.value = f"Exported: {path}"

    def on_restricted_mu(change):
        if not _alive[0]:
            return
        sim.update_restricted_mu(int(change['new']))

    def on_restart(_):
        if not _alive[0]:
            return
        _set_paused(True)
        time.sleep(0.02)   # let sim thread settle into paused state

        saved_state = sim.state_params

        # Re-initialize C library (realloc arrays, reset counters)
        sim._lib.evoca_init(N, sim.food_inc, sim.m_scale)
        sim._lib.evoca_set_gdiff(sim.gdiff)
        sim._lib.evoca_set_mu_lut(sim.mu_lut)
        sim._lib.evoca_set_mu_egenome(sim.mu_egenome)
        sim._lib.evoca_set_tax(sim.tax)
        sim._lib.evoca_set_restricted_mu(sim.restricted_mu)
        if pat_enabled:
            sim._lib.evoca_set_n_ent(sim.n_ent)

        # Replay state initialization
        sim.state(**saved_state)

        # Reset step counter
        st['step_cnt'] = 0
        ctrl[_STEP] = 0

        # Clear probe buffers
        if n_probes > 0:
            probe_cursor[0] = 0
            for pi in range(n_probes):
                probe_means[pi][:] = 0
                probe_stds[pi][:] = 0
        if activity_enabled:
            activity_cursor[0] = 0
            activity_pixels[:] = 0
        if eg_activity_enabled:
            eg_activity_cursor[0] = 0
            eg_activity_pixels[:] = 0
        if lut_complexity_enabled:
            lut_complexity_cursor[0] = 0
            lut_complexity_pixels[:] = 0
        if eg_pop_enabled:
            eg_pop_cursor[0] = 0
            eg_pop_pixels[:] = 0
        if entropy_enabled:
            entropy_cursor[0] = 0
            entropy_buf[:] = 0
        if pat_activity_enabled:
            pat_activity_cursor[0] = 0
            pat_activity_pixels[:] = 0
        if q_activity_enabled:
            q_activity_cursor[0] = 0
            for di in range(QA_N_DECILES):
                q_activity_deciles[di][:] = 0
        if ts_enabled:
            ts_cursor[0] = 0
            for ti in range(TS_N_TRACES):
                ts_bufs[ti][:] = 0

        sim.colorize(pixels, st['colormode'])
        status_lbl.value = "Restarted — t=0  (paused)"

    btn_pause.observe(on_pause_toggle, names='value')
    btn_restart.on_click(on_restart)
    btn_step.on_click(on_step)
    btn_step200.on_click(on_step200)
    btn_quit.on_click(on_quit)
    btn_save.on_click(on_save)
    btn_export.on_click(on_export)
    cb_restricted_mu.observe(on_restricted_mu, names='value')
    color_dd.observe(on_color, names='value')

    # Slider drag: pause on first touch, auto-resume after 200 ms idle
    def _make_slider_cb(attr, slider):
        _timer      = [None]
        _was_paused = [False]

        def on_value(change):
            if not _alive[0]:
                return
            if not st['paused']:
                _was_paused[0] = False
                _set_paused(True)
            else:
                _was_paused[0] = True
            getattr(sim, f"update_{attr}")(change['new'])
            if _timer[0] is not None:
                _timer[0].cancel()

            def _resume():
                if not _alive[0]:
                    return
                if not _was_paused[0]:
                    st['paused']  = False
                    ctrl[_PAUSED] = 0
                    # btn_pause visual won't sync from this thread — acceptable

            _timer[0] = threading.Timer(_SLIDER_RESUME_S, _resume)
            _timer[0].start()

        slider.observe(on_value, names='value')

    _make_slider_cb("food_inc",   sl_food_inc)
    _make_slider_cb("m_scale",    sl_m_scale)
    _make_slider_cb("gdiff",      sl_gdiff)
    _make_slider_cb("mu_lut",     sl_mu_lut)
    _make_slider_cb("mu_egenome",  sl_mu_egenome)
    _make_slider_cb("tax",        sl_tax)

    # ── Simulation thread ─────────────────────────────────────────
    def _sim_thread():
        FPS_ALPHA = 0.05
        fps       = 0.0
        t_last    = time.perf_counter()
        _t_loop_end = time.perf_counter()

        while st['running'] and ctrl[_QUIT] == 0:

            # Check if subprocess died unexpectedly (once, after a few steps)
            if st['step_cnt'] == 5:
                rc = sdl_proc.poll()
                if rc is not None:
                    print(f"EvoCA: SDL worker exited early (rc={rc})", flush=True)
                    status_lbl.value = f"SDL worker crashed (rc={rc}) — check terminal"

            if st['paused']:
                # Colorize so Step button results are shown; then rest.
                sim.colorize(pixels, st['colormode'])
                t_last = time.perf_counter()   # reset timer; avoid fps spike on resume
                _t_loop_end = t_last
                time.sleep(0.01)
                continue

            _t0 = time.perf_counter()
            _gap = _t0 - _t_loop_end  # time between end of last iteration and start of this one

            sim.step()
            st['step_cnt'] += 1

            _t1 = time.perf_counter() if diag else 0

            # colorize writes directly into shared pixel memory
            sim.colorize(pixels, st['colormode'])

            _t2 = time.perf_counter() if diag else 0

            # Record probe data
            if n_probes > 0:
                cur = int(probe_cursor[0])
                for pi, (getter, dt_p) in enumerate(probe_getters):
                    ptr = getter()
                    arr = np.ctypeslib.as_array(ptr, shape=(N * N,))
                    farr = arr.astype(np.float64) if dt_p != np.float32 else arr
                    probe_means[pi][cur] = farr.mean()
                    probe_stds[pi][cur]  = farr.std()
                probe_cursor[0] = (cur + 1) % PROBE_W

            # Record activity data
            if activity_enabled:
                sim._lib.evoca_activity_update()
                act_cur = int(activity_cursor[0])
                act_col_ptr = activity_col.ctypes.data_as(
                    ctypes.POINTER(ctypes.c_int32))
                sim._lib.evoca_activity_render_col(act_col_ptr, ACT_H)
                activity_pixels[:, act_cur] = activity_col
                activity_cursor[0] = (act_cur + 1) % PROBE_W

            # Record activity quantile data
            if q_activity_enabled:
                if not activity_enabled:
                    sim._lib.evoca_activity_update()
                qa_cur = int(q_activity_cursor[0])
                qa_col_ptr = q_activity_col.ctypes.data_as(
                    ctypes.POINTER(ctypes.c_float))
                sim._lib.evoca_q_activity_deciles(qa_col_ptr)
                for di in range(QA_N_DECILES):
                    q_activity_deciles[di][qa_cur] = q_activity_col[di]
                q_activity_cursor[0] = (qa_cur + 1) % PROBE_W

            _t3 = time.perf_counter() if diag else 0

            # Record egenome activity / population
            if eg_activity_enabled or eg_pop_enabled:
                sim._lib.evoca_eg_activity_update()
            if eg_activity_enabled:
                ega_cur = int(eg_activity_cursor[0])
                ega_col_ptr = eg_activity_col.ctypes.data_as(
                    ctypes.POINTER(ctypes.c_int32))
                sim._lib.evoca_eg_activity_render_col(ega_col_ptr, ACT_H)
                eg_activity_pixels[:, ega_cur] = eg_activity_col
                eg_activity_cursor[0] = (ega_cur + 1) % PROBE_W
                sim._lib.evoca_eg_activity_get(
                    ega_m_acts.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
                    ega_m_pops.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                    ega_m_cols.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)))
                ega_m_ymax[0] = sim._lib.evoca_get_eg_act_ymax()

            # Record LUT complexity data
            if lut_complexity_enabled:
                lc_cur = int(lut_complexity_cursor[0])
                lc_col_ptr = lut_complexity_col.ctypes.data_as(
                    ctypes.POINTER(ctypes.c_int32))
                sim._lib.evoca_lut_complexity_render_col(lc_col_ptr, PROBE_H)
                lut_complexity_pixels[:, lc_cur] = lut_complexity_col
                lut_complexity_cursor[0] = (lc_cur + 1) % PROBE_W

            # Record egenome population data
            if eg_pop_enabled:
                ep_cur = int(eg_pop_cursor[0])
                ep_col_ptr = eg_pop_col.ctypes.data_as(
                    ctypes.POINTER(ctypes.c_int32))
                sim._lib.evoca_eg_pop_render_col(ep_col_ptr, PROBE_H)
                eg_pop_pixels[:, ep_cur] = eg_pop_col
                eg_pop_cursor[0] = (ep_cur + 1) % PROBE_W

            # Record pattern entropy / pattern activity data
            if pat_enabled:
                sim._lib.evoca_pat_update()
                if entropy_enabled:
                    ec = int(entropy_cursor[0])
                    entropy_buf[ec] = sim._lib.evoca_get_entropy()
                    entropy_cursor[0] = (ec + 1) % PROBE_W
                if pat_activity_enabled:
                    pac = int(pat_activity_cursor[0])
                    pa_col_ptr = pat_activity_col.ctypes.data_as(
                        ctypes.POINTER(ctypes.c_int32))
                    sim._lib.evoca_pat_activity_render_col(pa_col_ptr, PAT_H)
                    pat_activity_pixels[:, pac] = pat_activity_col
                    pat_activity_cursor[0] = (pac + 1) % PROBE_W

            _record_ts()

            _t4 = time.perf_counter()
            if diag:
                _total = _t4 - _t0
                if _gap > 0.05 or _total > 0.1:  # log gaps>50ms or steps>100ms
                    print(f"DIAG t={st['step_cnt']}: "
                          f"gap={_gap*1000:.1f}ms "
                          f"step={(_t1-_t0)*1000:.1f}ms "
                          f"color={(_t2-_t1)*1000:.1f}ms "
                          f"activity={(_t3-_t2)*1000:.1f}ms "
                          f"probes={(_t4-_t3)*1000:.1f}ms "
                          f"total={_total*1000:.1f}ms",
                          flush=True)
            _t_loop_end = _t4

            # FPS (only meaningful when running)
            t_now = time.perf_counter()
            dt    = t_now - t_last
            t_last = t_now
            if dt > 0:
                fps = FPS_ALPHA * (1.0 / dt) + (1.0 - FPS_ALPHA) * fps

            # Update ctrl_shm so subprocess can update its window title
            sc           = st['step_cnt']
            ctrl[_STEP]  = sc
            ctrl[_FPS10] = int(fps * 10)

            # Update status label (infrequent — tolerate ZMQ non-safety)
            if sc % 100 == 0:
                status_lbl.value = f"t={sc}  fps={fps:.1f}"

        # Cleanup — also called by atexit if this thread doesn't finish in time
        st['running'] = False
        if _alive[0]:
            ctrl[_QUIT] = 1
        sdl_proc.wait(timeout=3)
        _do_cleanup()

    t = threading.Thread(target=_sim_thread, name="evoca-sim", daemon=True)
    t.start()
    return t
