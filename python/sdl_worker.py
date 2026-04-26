"""
sdl_worker.py — SDL2 display subprocess for EvoCA.

Launched by controls.run_with_controls() via subprocess.Popen.
All output goes to the terminal where Jupyter was started.

Usage (internal):
    python sdl_worker.py <pixel_shm> <ctrl_shm> <N> <px> [<probe_shm> <probe_names_csv>]

ctrl_shm layout (5 × int32)
    [0] quit      1 = exit
    [1] colormode 0/1/2
    [2] step_cnt
    [3] fps × 10
    [4] paused    0/1

probe_shm layout (if present)
    [0]      int32   cursor (write position, 0..PROBE_W-1)
    [4..]    per probe: float32[PROBE_W] mean, float32[PROBE_W] std
"""
import sys
import ctypes
import traceback

PROBE_W = 512
PROBE_H = 128
TITLE_BAR_H = 28   # macOS title bar estimate

# Magnifier window
MAG_CELLS = 25      # cells shown in magnifier
MAG_PX    = 8       # pixels per cell in magnifier
MAG_W     = MAG_CELLS * MAG_PX   # 200
MAG_H     = MAG_CELLS * MAG_PX   # 200
MAG_HALF  = MAG_CELLS // 2       # 12

# Egenome overlay
EG_CELL_PX = 5          # pixels per cell in overlay pattern
EG_PAT_N   = 5          # 5x5 pattern
EG_OVL_PX  = EG_PAT_N * EG_CELL_PX   # 25x25

def _c(argb):
    """Convert ARGB uint32 to np.int32 (avoids numpy deprecation warning)."""
    import numpy as np
    return np.int32(argb if argb < 0x80000000 else argb - 0x100000000)

# Probe color palette (ARGB) — converted to int32 once at import
BG_COLOR      = _c(0xFF1A1A1A)
CURSOR_COLOR  = _c(0xFF444444)
_PROBE_COLORS = {
    'env_food':  (_c(0xFF00CC00), _c(0xFF003300)),   # green
    'priv_food': (_c(0xFF4488FF), _c(0xFF112244)),   # blue
    'births':    (_c(0xFFFFFF00), _c(0xFF444400)),   # yellow
}
_DEFAULT_COLORS = (_c(0xFFCCCCCC), _c(0xFF333333))   # grey fallback

# ── Grouped time-series (ts) probe ───────────────────────────────
# Must match controls._TS_TRACES order.
TS_N_TRACES = 6
_TS_GROUPS  = ((0, 1, 2), (3, 4, 5))
TS_STRIPS   = len(_TS_GROUPS)
_TS_COLORS = (
    _c(0xFFFF4488),   # pop           — magenta
    _c(0xFF00CC44),   # F_mean        — green
    _c(0xFF4488FF),   # f_mean        — blue
    _c(0xFFFFAA22),   # lut_div       — orange
    _c(0xFFCC66FF),   # eg_ent        — purple
    _c(0xFF22EEDD),   # activity_flux — cyan/teal
)
_TS_LABELS = ('pop', 'F_mean', 'f_mean', 'lut_div', 'eg_ent', 'flux')
TS_SEP_COLOR = _c(0xFF444444)

# Decile colors for q_activity (p10 blue → p50 green → p90 red)
_QA_COLORS = [
    _c(0xFF3344FF),  # p10 - blue
    _c(0xFF2288FF),  # p20
    _c(0xFF00BBDD),  # p30 - cyan
    _c(0xFF00CC88),  # p40 - teal
    _c(0xFF44DD44),  # p50 - green (median)
    _c(0xFFBBBB00),  # p60 - yellow
    _c(0xFFFF8800),  # p70 - orange
    _c(0xFFFF4422),  # p80
    _c(0xFFFF1144),  # p90 - red
]


def _render_q_activity(dst, decile_bufs, cursor, global_max):
    """Render activity quantile strip chart with log-scaled Y axis.

    decile_bufs: list of 9 float32[PROBE_W] arrays (p10..p90).
    global_max:  all-time observed max (float); updated and returned.
    Returns:     updated global_max.
    """
    import numpy as np
    import math

    dst[:PROBE_H, :PROBE_W] = BG_COLOR

    # Roll each decile buffer so newest is on the right
    rolled = [np.roll(b, -cursor) for b in decile_bufs]

    # Find min/max of positive values in current window
    all_pos = []
    for r in rolled:
        pos = r[r > 0]
        if len(pos) > 0:
            all_pos.append(pos)
    if len(all_pos) == 0:
        return global_max
    cat = np.concatenate(all_pos)
    lo = float(cat.min())
    hi = float(cat.max())

    # Track all-time max so collapses remain visible after scrolling out
    global_max = max(global_max, hi)

    if lo <= 0:
        lo = float(cat[cat > 0].min()) if (cat > 0).any() else 1.0
    hi = global_max
    if hi <= lo:
        hi = lo * 10.0
    log_lo = math.log10(lo)
    log_hi = math.log10(hi)
    # Add 10% margin in log-space so top quantile doesn't touch the ceiling
    span = log_hi - log_lo
    if span < 0.01:
        span = 1.0
    log_hi += span * 0.10
    scale = (PROBE_H - 1) / (log_hi - log_lo)

    xs = np.arange(PROBE_W)
    for di in range(9):
        col = _QA_COLORS[di]
        r = rolled[di]
        mask = r > 0
        if not mask.any():
            continue
        lv = np.log10(r[mask])
        ys = np.clip(((log_hi - lv) * scale).astype(int), 0, PROBE_H - 1)
        dst[ys, xs[mask]] = col

    # Cursor line at right edge
    dst[:PROBE_H, PROBE_W - 1] = CURSOR_COLOR
    return global_max


def _render_probe(dst, pitch_i32, mean_buf, std_buf, cursor, color_mean, color_band,
                   y_fixed=None):
    """Render one strip chart into an SDL surface pixel array.

    y_fixed : optional (lo, hi) tuple to fix the Y axis range.
              If None, auto-scale from data.
    """
    import numpy as np

    # Clear to background
    dst[:PROBE_H, :PROBE_W] = BG_COLOR

    # Build time-ordered view: oldest on left, newest on right
    m = np.roll(mean_buf, -cursor)
    s = np.roll(std_buf, -cursor)

    lo_arr = m - s
    hi_arr = m + s
    # Only scale from non-zero data (buffer may be partially filled)
    filled = (m != 0) | (s != 0)
    if not filled.any():
        return

    if y_fixed is not None:
        y_min, y_max = y_fixed
    else:
        y_min = float(lo_arr[filled].min())
        y_max = float(hi_arr[filled].max())
    if y_max - y_min < 1e-8:
        y_min -= 0.01
        y_max += 0.01
    margin = (y_max - y_min) * 0.05
    y_min -= margin
    y_max += margin
    scale = (PROBE_H - 1) / (y_max - y_min)

    for x in range(PROBE_W):
        mv = m[x]
        sv = s[x]
        if mv == 0.0 and sv == 0.0 and not filled[x]:
            continue
        # Y pixel (0 = top = y_max, PROBE_H-1 = bottom = y_min)
        y_mean = int((y_max - mv) * scale)
        y_lo   = int((y_max - (mv + sv)) * scale)
        y_hi   = int((y_max - (mv - sv)) * scale)
        y_lo   = max(0, min(PROBE_H - 1, y_lo))
        y_hi   = max(0, min(PROBE_H - 1, y_hi))
        y_mean = max(0, min(PROBE_H - 1, y_mean))
        # Draw std band
        for y in range(y_lo, y_hi + 1):
            dst[y, x] = color_band
        # Draw mean line (overwrite band)
        dst[y_mean, x] = color_mean

    # Draw cursor position (right edge of newest data) as dim vertical line
    cx = PROBE_W - 1
    for y in range(PROBE_H):
        dst[y, cx] = CURSOR_COLOR


def _draw_ts_legend(dst, y0, trace_idxs):
    """Paint a small color-coded legend block in the top-left of the
    strip starting at row y0. Each trace gets a 3-row swatch plus a
    label drawn via a minimal 3x5 bitmap font."""
    import numpy as np
    # Minimal font — only chars we need for _TS_LABELS.
    _F = {
        'p': ["###", "# #", "###", "#  ", "#  "],
        'o': ["###", "# #", "# #", "# #", "###"],
        'F': ["###", "#  ", "## ", "#  ", "#  "],
        '_': ["   ", "   ", "   ", "   ", "###"],
        'm': ["# #", "###", "###", "# #", "# #"],
        'e': ["###", "#  ", "###", "#  ", "###"],
        'a': [" ##", "# #", "###", "# #", "# #"],
        'n': ["# #", "###", "###", "# #", "# #"],
        'f': [" ##", "#  ", "###", "#  ", "#  "],
        'l': ["#  ", "#  ", "#  ", "#  ", "###"],
        'u': ["# #", "# #", "# #", "# #", "###"],
        't': ["###", " # ", " # ", " # ", " # "],
        'd': ["## ", "# #", "# #", "# #", "## "],
        'i': ["#", "#", "#", "#", "#"],
        'v': ["# #", "# #", "# #", "# #", " # "],
        'x': ["# #", "# #", " # ", "# #", "# #"],
        'g': ["###", "# #", "###", "  #", "###"],
        ' ': ["  ", "  ", "  ", "  ", "  "],
    }
    ch_w = 4   # per-char width with 1px spacing
    sw_w = 8   # swatch width
    line_h = 8
    x = 3
    y = y0 + 2
    for ti in trace_idxs:
        col = _TS_COLORS[ti]
        label = _TS_LABELS[ti]
        dst[y + 1:y + 4, x:x + sw_w] = col
        lx = x + sw_w + 2
        for ch in label:
            rows = _F.get(ch, _F[' '])
            for rr, bits in enumerate(rows):
                for cc, b in enumerate(bits):
                    if b == '#':
                        dst[y + rr, lx + cc] = col
            lx += ch_w
        x += sw_w + 2 + ch_w * len(label) + 5
        if x > 256:
            break


def _render_ts(dst, bufs, cursor):
    """Render the grouped time-series probe into `dst`.

    `dst` has shape (2*PROBE_H, W+).  `bufs` is a list of 6 float32[PROBE_W]
    views.  Each trace auto-scales independently within its strip; observed
    [min, max] maps to the middle 80% of the strip's height (10% padding
    top and bottom). Adjacent samples are connected by vertical fill so
    each trace renders as a continuous line, not isolated dots. Zero
    samples are treated as uninitialized and excluded from autoscale and
    plotting.
    """
    import numpy as np

    H_total = TS_STRIPS * PROBE_H
    dst[:H_total, :PROBE_W] = BG_COLOR

    rolled = [np.roll(b, -cursor) for b in bufs]
    xs = np.arange(PROBE_W)

    for gi, group in enumerate(_TS_GROUPS):
        y0 = gi * PROBE_H
        strip = dst[y0:y0 + PROBE_H]
        y_top = 0.1 * (PROBE_H - 1)
        y_bot = 0.9 * (PROBE_H - 1)
        H_grid = np.arange(PROBE_H)[:, None]        # (H, 1)

        for ti in group:
            r = rolled[ti]
            mask = r != 0.0
            if not mask.any():
                continue
            vals = r[mask]
            lo = float(vals.min())
            hi = float(vals.max())
            if hi <= lo:
                ys_valid = np.full(int(mask.sum()),
                                   int((y_top + y_bot) * 0.5), dtype=np.int32)
            else:
                scale = (y_bot - y_top) / (hi - lo)
                ys_valid = np.clip(
                    (y_bot - (vals - lo) * scale).astype(np.int32),
                    0, PROBE_H - 1)

            ys = np.full(PROBE_W, -1, dtype=np.int32)
            ys[mask] = ys_valid

            prev_ys = np.roll(ys, 1)
            prev_ys[0] = -1
            both = (ys >= 0) & (prev_ys >= 0)
            y_lo = np.minimum(ys, prev_ys)
            y_hi = np.maximum(ys, prev_ys)
            fill = ((H_grid >= y_lo[None, :])
                    & (H_grid <= y_hi[None, :])
                    & both[None, :])

            col = _TS_COLORS[ti % len(_TS_COLORS)]
            strip[fill] = col
            isolated = mask & ~both
            strip[ys[isolated], xs[isolated]] = col

        # Cursor edge + legend
        strip[:, PROBE_W - 1] = CURSOR_COLOR
        _draw_ts_legend(dst, y0, group)
        if gi < TS_STRIPS - 1:
            dst[y0 + PROBE_H - 1, :PROBE_W] = TS_SEP_COLOR


def main():
    if len(sys.argv) < 5:
        print("EvoCA SDL: bad args", flush=True)
        sys.exit(1)

    pixel_shm_name = sys.argv[1]
    ctrl_shm_name  = sys.argv[2]
    N  = int(sys.argv[3])
    px = int(sys.argv[4])
    W, H = N * px, N * px

    # Optional probe args (positional: argv[5]=probe_shm, argv[6]=names_csv)
    probe_shm_name = sys.argv[5] if len(sys.argv) > 6 else None
    probe_names    = sys.argv[6].split(",") if len(sys.argv) > 6 else []
    n_probes       = len(probe_names)

    # Optional activity probe (--activity=<shm_name>)
    activity_shm_name = None
    eg_activity_shm_name = None
    lut_complexity_shm_name = None
    entropy_shm_name = None
    pat_activity_shm_name = None
    eg_pop_shm_name = None
    q_activity_shm_name = None
    n_activity_shm_name  = None
    nq_activity_shm_name = None
    ts_shm_name = None
    for arg in sys.argv:
        if arg.startswith("--activity="):
            activity_shm_name = arg[len("--activity="):]
        elif arg.startswith("--eg-activity="):
            eg_activity_shm_name = arg[len("--eg-activity="):]
        elif arg.startswith("--lut-complexity="):
            lut_complexity_shm_name = arg[len("--lut-complexity="):]
        elif arg.startswith("--entropy="):
            entropy_shm_name = arg[len("--entropy="):]
        elif arg.startswith("--pat-activity="):
            pat_activity_shm_name = arg[len("--pat-activity="):]
        elif arg.startswith("--eg-pop="):
            eg_pop_shm_name = arg[len("--eg-pop="):]
        elif arg.startswith("--q-activity="):
            q_activity_shm_name = arg[len("--q-activity="):]
        elif arg.startswith("--n-activity="):
            n_activity_shm_name = arg[len("--n-activity="):]
        elif arg.startswith("--nq-activity="):
            nq_activity_shm_name = arg[len("--nq-activity="):]
        elif arg.startswith("--ts="):
            ts_shm_name = arg[len("--ts="):]

    ACT_H = 2 * PROBE_H  # 256

    print(f"EvoCA SDL: starting  N={N} px={px}  probes={probe_names}", flush=True)

    import numpy as np
    from multiprocessing.shared_memory import SharedMemory
    import sdl2

    # Open shared memory
    try:
        pixel_shm = SharedMemory(name=pixel_shm_name)
        ctrl_shm  = SharedMemory(name=ctrl_shm_name)
    except Exception as e:
        print(f"EvoCA SDL: SharedMemory open failed: {e}", flush=True)
        sys.exit(1)

    pixels = np.ndarray((N * N,), dtype=np.int32, buffer=pixel_shm.buf)
    ctrl   = np.ndarray((5,),     dtype=np.int32, buffer=ctrl_shm.buf)

    # Open probe shared memory
    probe_shm    = None
    probe_cursor = None
    probe_means  = []
    probe_stds   = []
    if n_probes > 0 and probe_shm_name:
        try:
            probe_shm = SharedMemory(name=probe_shm_name)
        except Exception as e:
            print(f"EvoCA SDL: probe SharedMemory open failed: {e}", flush=True)
            n_probes = 0
            probe_names = []
        if probe_shm is not None:
            probe_cursor = np.ndarray((1,), dtype=np.int32, buffer=probe_shm.buf)
            off = 4
            for _ in range(n_probes):
                probe_means.append(np.ndarray((PROBE_W,), dtype=np.float32,
                                              buffer=probe_shm.buf, offset=off))
                off += PROBE_W * 4
                probe_stds.append(np.ndarray((PROBE_W,), dtype=np.float32,
                                             buffer=probe_shm.buf, offset=off))
                off += PROBE_W * 4

    # Open activity shared memory
    activity_shm     = None
    activity_cursor  = None
    activity_pixels  = None
    if activity_shm_name:
        try:
            activity_shm = SharedMemory(name=activity_shm_name)
            activity_cursor = np.ndarray((1,), dtype=np.int32,
                                         buffer=activity_shm.buf)
            activity_pixels = np.ndarray((ACT_H, PROBE_W), dtype=np.int32,
                                         buffer=activity_shm.buf, offset=4)
            print(f"EvoCA SDL: activity shm opened ({ACT_H}x{PROBE_W})",
                  flush=True)
        except Exception as e:
            print(f"EvoCA SDL: activity SharedMemory open failed: {e}",
                  flush=True)
            activity_shm_name = None

    # Open egenome activity shared memory
    eg_activity_shm     = None
    eg_activity_cursor  = None
    eg_activity_pixels  = None
    ega_m_acts = None
    ega_m_pops = None
    ega_m_cols = None
    ega_m_ymax = None
    if eg_activity_shm_name:
        try:
            eg_activity_shm = SharedMemory(name=eg_activity_shm_name)
            eg_activity_cursor = np.ndarray((1,), dtype=np.int32,
                                            buffer=eg_activity_shm.buf)
            eg_activity_pixels = np.ndarray((ACT_H, PROBE_W), dtype=np.int32,
                                            buffer=eg_activity_shm.buf, offset=4)
            ega_meta_off = 4 + PROBE_W * ACT_H * 4
            ega_m_acts = np.ndarray((64,), dtype=np.uint64,
                                    buffer=eg_activity_shm.buf, offset=ega_meta_off)
            ega_m_pops = np.ndarray((64,), dtype=np.uint32,
                                    buffer=eg_activity_shm.buf, offset=ega_meta_off + 64*8)
            ega_m_cols = np.ndarray((64,), dtype=np.int32,
                                    buffer=eg_activity_shm.buf, offset=ega_meta_off + 64*8 + 64*4)
            ega_m_ymax = np.ndarray((1,), dtype=np.int32,
                                    buffer=eg_activity_shm.buf, offset=ega_meta_off + 64*8 + 64*4*2)
            print(f"EvoCA SDL: eg_activity shm opened ({ACT_H}x{PROBE_W})",
                  flush=True)
        except Exception as e:
            print(f"EvoCA SDL: eg_activity SharedMemory open failed: {e}",
                  flush=True)
            eg_activity_shm_name = None

    # Open LUT complexity shared memory
    lut_complexity_shm     = None
    lut_complexity_cursor  = None
    lut_complexity_pixels  = None
    if lut_complexity_shm_name:
        try:
            lut_complexity_shm = SharedMemory(name=lut_complexity_shm_name)
            lut_complexity_cursor = np.ndarray((1,), dtype=np.int32,
                                               buffer=lut_complexity_shm.buf)
            lut_complexity_pixels = np.ndarray((PROBE_H, PROBE_W), dtype=np.int32,
                                               buffer=lut_complexity_shm.buf, offset=4)
            print(f"EvoCA SDL: lut_complexity shm opened ({PROBE_H}x{PROBE_W})",
                  flush=True)
        except Exception as e:
            print(f"EvoCA SDL: lut_complexity SharedMemory open failed: {e}",
                  flush=True)
            lut_complexity_shm_name = None

    # Open egenome population shared memory
    eg_pop_shm     = None
    eg_pop_cursor  = None
    eg_pop_pixels  = None
    if eg_pop_shm_name:
        try:
            eg_pop_shm = SharedMemory(name=eg_pop_shm_name)
            eg_pop_cursor = np.ndarray((1,), dtype=np.int32,
                                        buffer=eg_pop_shm.buf)
            eg_pop_pixels = np.ndarray((PROBE_H, PROBE_W), dtype=np.int32,
                                        buffer=eg_pop_shm.buf, offset=4)
            print(f"EvoCA SDL: eg_pop shm opened ({PROBE_H}x{PROBE_W})",
                  flush=True)
        except Exception as e:
            print(f"EvoCA SDL: eg_pop SharedMemory open failed: {e}",
                  flush=True)
            eg_pop_shm_name = None

    # Open activity quantile shared memory
    QA_N_DECILES        = 9
    q_activity_shm      = None
    q_activity_cursor   = None
    q_activity_deciles  = None   # list of 9 float32[PROBE_W]
    if q_activity_shm_name:
        try:
            q_activity_shm = SharedMemory(name=q_activity_shm_name)
            q_activity_cursor = np.ndarray((1,), dtype=np.int32,
                                            buffer=q_activity_shm.buf)
            q_activity_deciles = []
            off = 4
            for _ in range(QA_N_DECILES):
                q_activity_deciles.append(
                    np.ndarray((PROBE_W,), dtype=np.float32,
                               buffer=q_activity_shm.buf, offset=off))
                off += PROBE_W * 4
            print(f"EvoCA SDL: q_activity shm opened ({QA_N_DECILES}x{PROBE_W})",
                  flush=True)
        except Exception as e:
            print(f"EvoCA SDL: q_activity SharedMemory open failed: {e}",
                  flush=True)
            q_activity_shm_name = None

    # Open N-activity (Channon shadow) shared memory
    n_activity_shm    = None
    n_activity_cursor = None
    n_activity_pixels = None
    if n_activity_shm_name:
        try:
            n_activity_shm = SharedMemory(name=n_activity_shm_name)
            n_activity_cursor = np.ndarray((1,), dtype=np.int32,
                                            buffer=n_activity_shm.buf)
            n_activity_pixels = np.ndarray((ACT_H, PROBE_W), dtype=np.int32,
                                            buffer=n_activity_shm.buf, offset=4)
            print(f"EvoCA SDL: n_activity shm opened ({ACT_H}x{PROBE_W})",
                  flush=True)
        except Exception as e:
            print(f"EvoCA SDL: n_activity SharedMemory open failed: {e}",
                  flush=True)
            n_activity_shm_name = None

    # Open Nq-activity (deciles) shared memory
    nq_activity_shm     = None
    nq_activity_cursor  = None
    nq_activity_deciles = None
    if nq_activity_shm_name:
        try:
            nq_activity_shm = SharedMemory(name=nq_activity_shm_name)
            nq_activity_cursor = np.ndarray((1,), dtype=np.int32,
                                             buffer=nq_activity_shm.buf)
            nq_activity_deciles = []
            off = 4
            for _ in range(QA_N_DECILES):
                nq_activity_deciles.append(
                    np.ndarray((PROBE_W,), dtype=np.float32,
                               buffer=nq_activity_shm.buf, offset=off))
                off += PROBE_W * 4
            print(f"EvoCA SDL: nq_activity shm opened "
                  f"({QA_N_DECILES}x{PROBE_W})", flush=True)
        except Exception as e:
            print(f"EvoCA SDL: nq_activity SharedMemory open failed: {e}",
                  flush=True)
            nq_activity_shm_name = None

    # Open ts shared memory
    ts_shm     = None
    ts_cursor  = None
    ts_bufs    = None
    if ts_shm_name:
        try:
            ts_shm = SharedMemory(name=ts_shm_name)
            ts_cursor = np.ndarray((1,), dtype=np.int32, buffer=ts_shm.buf)
            ts_bufs = []
            off = 4
            for _ in range(TS_N_TRACES):
                ts_bufs.append(np.ndarray((PROBE_W,), dtype=np.float32,
                                           buffer=ts_shm.buf, offset=off))
                off += PROBE_W * 4
            print(f"EvoCA SDL: ts shm opened ({TS_N_TRACES} traces)", flush=True)
        except Exception as e:
            print(f"EvoCA SDL: ts SharedMemory open failed: {e}", flush=True)
            ts_shm_name = None

    COLOR_MODES = ["state", "env-food", "priv-food", "births", "age"]

    # ── SDL2 init ─────────────────────────────────────────────────
    if sdl2.SDL_Init(sdl2.SDL_INIT_VIDEO) != 0:
        print(f"EvoCA SDL: SDL_Init failed: {sdl2.SDL_GetError()}", flush=True)
        ctrl[0] = 1
        sys.exit(1)

    print("EvoCA SDL: SDL_Init OK", flush=True)

    # ── Screen geometry for window placement ──────────────────────
    dm = sdl2.SDL_DisplayMode()
    sdl2.SDL_GetCurrentDisplayMode(0, ctypes.byref(dm))
    scr_w, scr_h = dm.w, dm.h

    # Main window: upper-right corner
    main_x = scr_w - W
    main_y = 0

    window_p = sdl2.SDL_CreateWindow(
        b"EvoCA",
        main_x, main_y,
        W, H,
        sdl2.SDL_WINDOW_SHOWN,
    )
    if not window_p:
        print(f"EvoCA SDL: SDL_CreateWindow failed: {sdl2.SDL_GetError()}", flush=True)
        ctrl[0] = 1
        sdl2.SDL_Quit()
        sys.exit(1)

    sdl2.SDL_RaiseWindow(window_p)
    main_window_id = sdl2.SDL_GetWindowID(window_p)
    print("EvoCA SDL: window created and raised", flush=True)

    # Magnifier state (center computed dynamically from window position)
    mag_window_p  = None
    mag_surface_p = None
    mag_dst       = None

    surface_p = sdl2.SDL_GetWindowSurface(window_p)
    if not surface_p:
        print(f"EvoCA SDL: SDL_GetWindowSurface failed: {sdl2.SDL_GetError()}", flush=True)
        ctrl[0] = 1
        sdl2.SDL_DestroyWindow(window_p)
        sdl2.SDL_Quit()
        sys.exit(1)

    sdl2.SDL_SetSurfaceBlendMode(surface_p, sdl2.SDL_BLENDMODE_NONE)

    surf       = surface_p.contents
    pitch_i32  = surf.pitch // 4
    pixels_ptr = ctypes.cast(surf.pixels, ctypes.POINTER(ctypes.c_int32))
    dst_flat   = np.ctypeslib.as_array(pixels_ptr, shape=(H * pitch_i32,))
    dst        = dst_flat.reshape(H, pitch_i32)

    # ── Probe windows ────────────────────────────────────────────
    # Fixed Y ranges for known probes (both food fields are in [0, 1])
    _Y_FIXED = {
        'env_food':  (0.0, 1.0),
        'priv_food': (0.0, 1.0),
        'births':    (0.0, 1.0),
    }

    probe_windows  = []   # SDL_Window pointers
    probe_surfaces = []   # SDL_Surface pointers
    probe_dsts     = []   # numpy pixel arrays
    probe_colors   = []   # (color_mean, color_band) tuples
    probe_yfix     = []   # fixed Y ranges (or None for auto-scale)

    # Stack probe windows top-down.  After creating each window, read back
    # its actual Y position so the next window is placed correctly below it
    # (macOS may adjust positions for the menu bar and window decorations).
    next_probe_y = 0   # requested Y for the next probe window
    real_title_h = TITLE_BAR_H  # updated from first window's actual position

    for i, pname in enumerate(probe_names):
        pw_x = main_x - PROBE_W
        pw = sdl2.SDL_CreateWindow(
            pname.encode(),
            pw_x, next_probe_y,
            PROBE_W, PROBE_H,
            sdl2.SDL_WINDOW_SHOWN,
        )
        if not pw:
            print(f"EvoCA SDL: probe window '{pname}' failed", flush=True)
            continue

        # Read back actual content-area position
        actual_y = ctypes.c_int(0)
        sdl2.SDL_GetWindowPosition(pw, None, ctypes.byref(actual_y))
        # Measure the window's title-bar decoration height (top border)
        if i == 0:
            top_border = ctypes.c_int(0)
            ret = sdl2.SDL_GetWindowBordersSize(
                pw, ctypes.byref(top_border), None, None, None)
            if ret == 0 and top_border.value > 0:
                real_title_h = top_border.value
            print(f"EvoCA SDL: title bar = {real_title_h}px", flush=True)
        print(f"EvoCA SDL: probe '{pname}' y={actual_y.value}", flush=True)
        # Next window: below this content + room for next window's title bar
        next_probe_y = actual_y.value + PROBE_H + real_title_h

        ps = sdl2.SDL_GetWindowSurface(pw)
        if not ps:
            sdl2.SDL_DestroyWindow(pw)
            continue
        sdl2.SDL_SetSurfaceBlendMode(ps, sdl2.SDL_BLENDMODE_NONE)
        psurf     = ps.contents
        pp_i32    = psurf.pitch // 4
        pp_ptr    = ctypes.cast(psurf.pixels, ctypes.POINTER(ctypes.c_int32))
        pd_flat   = np.ctypeslib.as_array(pp_ptr, shape=(PROBE_H * pp_i32,))
        pd        = pd_flat.reshape(PROBE_H, pp_i32)

        probe_windows.append(pw)
        probe_surfaces.append(ps)
        probe_dsts.append(pd)
        probe_colors.append(_PROBE_COLORS.get(pname, _DEFAULT_COLORS))
        probe_yfix.append(_Y_FIXED.get(pname))

    print(f"EvoCA SDL: {len(probe_windows)} probe window(s) created", flush=True)

    # ── Activity window ──────────────────────────────────────────
    act_window_p  = None
    act_surface_p = None
    act_dst       = None
    if activity_shm is not None:
        aw_x = main_x - PROBE_W
        aw = sdl2.SDL_CreateWindow(
            b"activity",
            aw_x, next_probe_y,
            PROBE_W, ACT_H,
            sdl2.SDL_WINDOW_SHOWN,
        )
        if aw:
            actual_y = ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(aw, None, ctypes.byref(actual_y))
            next_probe_y = actual_y.value + ACT_H + real_title_h
            aps = sdl2.SDL_GetWindowSurface(aw)
            if aps:
                sdl2.SDL_SetSurfaceBlendMode(aps, sdl2.SDL_BLENDMODE_NONE)
                asurf   = aps.contents
                ap_i32  = asurf.pitch // 4
                ap_ptr  = ctypes.cast(asurf.pixels,
                                      ctypes.POINTER(ctypes.c_int32))
                ad_flat = np.ctypeslib.as_array(ap_ptr,
                                                shape=(ACT_H * ap_i32,))
                act_dst       = ad_flat.reshape(ACT_H, ap_i32)
                act_window_p  = aw
                act_surface_p = aps
                print("EvoCA SDL: activity window created", flush=True)
            else:
                sdl2.SDL_DestroyWindow(aw)
        else:
            print("EvoCA SDL: activity window creation failed", flush=True)

    # ── Egenome activity window ─────────────────────────────────
    eg_act_window_p  = None
    eg_act_surface_p = None
    eg_act_dst       = None
    eg_act_window_id = 0
    eg_overlay_val   = -1       # -1 = no overlay
    if eg_activity_shm is not None:
        caw_x = main_x - PROBE_W
        caw = sdl2.SDL_CreateWindow(
            b"eg_activity",
            caw_x, next_probe_y,
            PROBE_W, ACT_H,
            sdl2.SDL_WINDOW_SHOWN,
        )
        if caw:
            actual_y = ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(caw, None, ctypes.byref(actual_y))
            next_probe_y = actual_y.value + ACT_H + real_title_h
            caps = sdl2.SDL_GetWindowSurface(caw)
            if caps:
                sdl2.SDL_SetSurfaceBlendMode(caps, sdl2.SDL_BLENDMODE_NONE)
                casurf   = caps.contents
                cap_i32  = casurf.pitch // 4
                cap_ptr  = ctypes.cast(casurf.pixels,
                                       ctypes.POINTER(ctypes.c_int32))
                cad_flat = np.ctypeslib.as_array(cap_ptr,
                                                 shape=(ACT_H * cap_i32,))
                eg_act_dst       = cad_flat.reshape(ACT_H, cap_i32)
                eg_act_window_p  = caw
                eg_act_surface_p = caps
                eg_act_window_id = sdl2.SDL_GetWindowID(caw)
                print("EvoCA SDL: eg_activity window created", flush=True)
            else:
                sdl2.SDL_DestroyWindow(caw)
        else:
            print("EvoCA SDL: eg_activity window creation failed", flush=True)

    # ── LUT complexity window ──────────────────────────────────
    lc_window_p  = None
    lc_surface_p = None
    lc_dst       = None
    if lut_complexity_shm is not None:
        lcw_x = main_x - PROBE_W
        lcw = sdl2.SDL_CreateWindow(
            b"lut_complexity",
            lcw_x, next_probe_y,
            PROBE_W, PROBE_H,
            sdl2.SDL_WINDOW_SHOWN,
        )
        if lcw:
            actual_y = ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(lcw, None, ctypes.byref(actual_y))
            next_probe_y = actual_y.value + PROBE_H + real_title_h
            lps = sdl2.SDL_GetWindowSurface(lcw)
            if lps:
                sdl2.SDL_SetSurfaceBlendMode(lps, sdl2.SDL_BLENDMODE_NONE)
                lsurf   = lps.contents
                lp_i32  = lsurf.pitch // 4
                lp_ptr  = ctypes.cast(lsurf.pixels,
                                      ctypes.POINTER(ctypes.c_int32))
                ld_flat = np.ctypeslib.as_array(lp_ptr,
                                                shape=(PROBE_H * lp_i32,))
                lc_dst       = ld_flat.reshape(PROBE_H, lp_i32)
                lc_window_p  = lcw
                lc_surface_p = lps
                print("EvoCA SDL: lut_complexity window created", flush=True)
            else:
                sdl2.SDL_DestroyWindow(lcw)
        else:
            print("EvoCA SDL: lut_complexity window creation failed", flush=True)

    # ── Egenome population window ──────────────────────────────
    ep_window_p  = None
    ep_surface_p = None
    ep_dst       = None
    ep_window_id = 0
    ep_overlay_val = -1         # -1 = no overlay
    if eg_pop_shm is not None:
        epw_x = main_x - PROBE_W
        epw = sdl2.SDL_CreateWindow(
            b"eg_pop",
            epw_x, next_probe_y,
            PROBE_W, PROBE_H,
            sdl2.SDL_WINDOW_SHOWN,
        )
        if epw:
            actual_y = ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(epw, None, ctypes.byref(actual_y))
            next_probe_y = actual_y.value + PROBE_H + real_title_h
            eps = sdl2.SDL_GetWindowSurface(epw)
            if eps:
                sdl2.SDL_SetSurfaceBlendMode(eps, sdl2.SDL_BLENDMODE_NONE)
                esurf   = eps.contents
                ep_i32  = esurf.pitch // 4
                ep_ptr  = ctypes.cast(esurf.pixels,
                                      ctypes.POINTER(ctypes.c_int32))
                ed_flat = np.ctypeslib.as_array(ep_ptr,
                                                shape=(PROBE_H * ep_i32,))
                ep_dst       = ed_flat.reshape(PROBE_H, ep_i32)
                ep_window_p  = epw
                ep_surface_p = eps
                ep_window_id = sdl2.SDL_GetWindowID(epw)
                print("EvoCA SDL: eg_pop window created", flush=True)
            else:
                sdl2.SDL_DestroyWindow(epw)
        else:
            print("EvoCA SDL: eg_pop window creation failed", flush=True)

    # Create q_activity window
    qa_window_p    = None
    qa_surface_p   = None
    qa_dst         = None
    qa_global_max  = 0.0
    if q_activity_shm is not None:
        qaw_x = main_x - PROBE_W
        qaw = sdl2.SDL_CreateWindow(
            b"q_activity",
            qaw_x, next_probe_y,
            PROBE_W, PROBE_H,
            sdl2.SDL_WINDOW_SHOWN,
        )
        if qaw:
            actual_y = ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(qaw, None, ctypes.byref(actual_y))
            next_probe_y = actual_y.value + PROBE_H + real_title_h
            qaps = sdl2.SDL_GetWindowSurface(qaw)
            if qaps:
                sdl2.SDL_SetSurfaceBlendMode(qaps, sdl2.SDL_BLENDMODE_NONE)
                qasurf  = qaps.contents
                qap_i32 = qasurf.pitch // 4
                qap_ptr = ctypes.cast(qasurf.pixels,
                                       ctypes.POINTER(ctypes.c_int32))
                qad_flat = np.ctypeslib.as_array(qap_ptr,
                                                  shape=(PROBE_H * qap_i32,))
                qa_dst       = qad_flat.reshape(PROBE_H, qap_i32)
                qa_window_p  = qaw
                qa_surface_p = qaps
                print("EvoCA SDL: q_activity window created", flush=True)
            else:
                sdl2.SDL_DestroyWindow(qaw)
        else:
            print("EvoCA SDL: q_activity window creation failed", flush=True)

    # Create N-activity window (Channon shadow)
    n_window_p   = None
    n_surface_p  = None
    n_dst        = None
    if n_activity_shm is not None:
        nw_x = main_x - PROBE_W
        nw = sdl2.SDL_CreateWindow(
            b"n_activity (shadow)",
            nw_x, next_probe_y,
            PROBE_W, ACT_H,
            sdl2.SDL_WINDOW_SHOWN,
        )
        if nw:
            actual_y = ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(nw, None, ctypes.byref(actual_y))
            next_probe_y = actual_y.value + ACT_H + real_title_h
            nps = sdl2.SDL_GetWindowSurface(nw)
            if nps:
                sdl2.SDL_SetSurfaceBlendMode(nps, sdl2.SDL_BLENDMODE_NONE)
                nsurf  = nps.contents
                np_i32 = nsurf.pitch // 4
                np_ptr = ctypes.cast(nsurf.pixels,
                                      ctypes.POINTER(ctypes.c_int32))
                nd_flat = np.ctypeslib.as_array(np_ptr,
                                                 shape=(ACT_H * np_i32,))
                n_dst        = nd_flat.reshape(ACT_H, np_i32)
                n_window_p   = nw
                n_surface_p  = nps
                print("EvoCA SDL: n_activity window created", flush=True)
            else:
                sdl2.SDL_DestroyWindow(nw)
        else:
            print("EvoCA SDL: n_activity window creation failed", flush=True)

    # Create Nq-activity window
    nq_window_p   = None
    nq_surface_p  = None
    nq_dst        = None
    nq_global_max = 0.0
    if nq_activity_shm is not None:
        nqw_x = main_x - PROBE_W
        nqw = sdl2.SDL_CreateWindow(
            b"nq_activity (shadow)",
            nqw_x, next_probe_y,
            PROBE_W, PROBE_H,
            sdl2.SDL_WINDOW_SHOWN,
        )
        if nqw:
            actual_y = ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(nqw, None, ctypes.byref(actual_y))
            next_probe_y = actual_y.value + PROBE_H + real_title_h
            nqps = sdl2.SDL_GetWindowSurface(nqw)
            if nqps:
                sdl2.SDL_SetSurfaceBlendMode(nqps, sdl2.SDL_BLENDMODE_NONE)
                nqsurf  = nqps.contents
                nqp_i32 = nqsurf.pitch // 4
                nqp_ptr = ctypes.cast(nqsurf.pixels,
                                       ctypes.POINTER(ctypes.c_int32))
                nqd_flat = np.ctypeslib.as_array(nqp_ptr,
                                                  shape=(PROBE_H * nqp_i32,))
                nq_dst       = nqd_flat.reshape(PROBE_H, nqp_i32)
                nq_window_p  = nqw
                nq_surface_p = nqps
                print("EvoCA SDL: nq_activity window created", flush=True)
            else:
                sdl2.SDL_DestroyWindow(nqw)
        else:
            print("EvoCA SDL: nq_activity window creation failed", flush=True)

    # Open entropy shared memory + create window
    entropy_shm      = None
    entropy_cursor   = None
    entropy_buf      = None
    entropy_std_buf  = None
    ent_window_p     = None
    ent_surface_p    = None
    ent_dst          = None
    if entropy_shm_name:
        try:
            entropy_shm = SharedMemory(name=entropy_shm_name)
            entropy_cursor = np.ndarray((1,), dtype=np.int32,
                                        buffer=entropy_shm.buf)
            entropy_buf = np.ndarray((PROBE_W,), dtype=np.float32,
                                     buffer=entropy_shm.buf, offset=4)
            entropy_std_buf = np.zeros(PROBE_W, dtype=np.float32)
            print(f"EvoCA SDL: entropy shm opened", flush=True)
            ew_x = main_x - PROBE_W
            ew = sdl2.SDL_CreateWindow(
                b"entropy", ew_x, next_probe_y, PROBE_W, PROBE_H,
                sdl2.SDL_WINDOW_SHOWN)
            if ew:
                actual_y = ctypes.c_int(0)
                sdl2.SDL_GetWindowPosition(ew, None, ctypes.byref(actual_y))
                next_probe_y = actual_y.value + PROBE_H + real_title_h
                eps = sdl2.SDL_GetWindowSurface(ew)
                if eps:
                    sdl2.SDL_SetSurfaceBlendMode(eps, sdl2.SDL_BLENDMODE_NONE)
                    esurf  = eps.contents
                    ep_i32 = esurf.pitch // 4
                    ep_ptr = ctypes.cast(esurf.pixels,
                                         ctypes.POINTER(ctypes.c_int32))
                    ed_flat = np.ctypeslib.as_array(ep_ptr,
                                                    shape=(PROBE_H * ep_i32,))
                    ent_dst      = ed_flat.reshape(PROBE_H, ep_i32)
                    ent_window_p  = ew
                    ent_surface_p = eps
                    print("EvoCA SDL: entropy window created", flush=True)
                else:
                    sdl2.SDL_DestroyWindow(ew)
            else:
                print("EvoCA SDL: entropy window creation failed", flush=True)
        except Exception as e:
            print(f"EvoCA SDL: entropy SharedMemory open failed: {e}", flush=True)
            entropy_shm_name = None

    # Open pattern activity shared memory + create window
    pat_activity_shm     = None
    pat_activity_cursor  = None
    pat_activity_pixels  = None
    pa_window_p  = None
    pa_surface_p = None
    pa_dst       = None
    if pat_activity_shm_name:
        try:
            pat_activity_shm = SharedMemory(name=pat_activity_shm_name)
            pat_activity_cursor = np.ndarray((1,), dtype=np.int32,
                                             buffer=pat_activity_shm.buf)
            pat_activity_pixels = np.ndarray((ACT_H, PROBE_W), dtype=np.int32,
                                             buffer=pat_activity_shm.buf, offset=4)
            print(f"EvoCA SDL: pat_activity shm opened ({ACT_H}x{PROBE_W})", flush=True)
            paw_x = main_x - PROBE_W
            paw = sdl2.SDL_CreateWindow(
                b"pat_activity", paw_x, next_probe_y, PROBE_W, ACT_H,
                sdl2.SDL_WINDOW_SHOWN)
            if paw:
                actual_y = ctypes.c_int(0)
                sdl2.SDL_GetWindowPosition(paw, None, ctypes.byref(actual_y))
                next_probe_y = actual_y.value + ACT_H + real_title_h
                paps = sdl2.SDL_GetWindowSurface(paw)
                if paps:
                    sdl2.SDL_SetSurfaceBlendMode(paps, sdl2.SDL_BLENDMODE_NONE)
                    pasurf   = paps.contents
                    pap_i32  = pasurf.pitch // 4
                    pap_ptr  = ctypes.cast(pasurf.pixels,
                                           ctypes.POINTER(ctypes.c_int32))
                    pad_flat = np.ctypeslib.as_array(pap_ptr,
                                                     shape=(ACT_H * pap_i32,))
                    pa_dst       = pad_flat.reshape(ACT_H, pap_i32)
                    pa_window_p  = paw
                    pa_surface_p = paps
                    print("EvoCA SDL: pat_activity window created", flush=True)
                else:
                    sdl2.SDL_DestroyWindow(paw)
            else:
                print("EvoCA SDL: pat_activity window creation failed", flush=True)
        except Exception as e:
            print(f"EvoCA SDL: pat_activity SharedMemory open failed: {e}", flush=True)
            pat_activity_shm_name = None

    # ── ts window (2×PROBE_H tall) ─────────────────────────────────
    TS_H = TS_STRIPS * PROBE_H
    ts_window_p  = None
    ts_surface_p = None
    ts_dst       = None
    if ts_shm is not None:
        tw_x = main_x - PROBE_W
        tw = sdl2.SDL_CreateWindow(
            b"ts", tw_x, next_probe_y, PROBE_W, TS_H,
            sdl2.SDL_WINDOW_SHOWN)
        if tw:
            actual_y = ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(tw, None, ctypes.byref(actual_y))
            next_probe_y = actual_y.value + TS_H + real_title_h
            tps = sdl2.SDL_GetWindowSurface(tw)
            if tps:
                sdl2.SDL_SetSurfaceBlendMode(tps, sdl2.SDL_BLENDMODE_NONE)
                tsurf  = tps.contents
                tp_i32 = tsurf.pitch // 4
                tp_ptr = ctypes.cast(tsurf.pixels,
                                      ctypes.POINTER(ctypes.c_int32))
                td_flat = np.ctypeslib.as_array(tp_ptr,
                                                 shape=(TS_H * tp_i32,))
                ts_dst       = td_flat.reshape(TS_H, tp_i32)
                ts_window_p  = tw
                ts_surface_p = tps
                print("EvoCA SDL: ts window created", flush=True)
            else:
                sdl2.SDL_DestroyWindow(tw)
        else:
            print("EvoCA SDL: ts window creation failed", flush=True)

    print("EvoCA SDL: entering main loop", flush=True)

    event = sdl2.SDL_Event()
    pending_click = None          # deferred probe click: ('ega'|'ep', x, y)

    while ctrl[0] == 0:

        # Process SDL events
        while sdl2.SDL_PollEvent(ctypes.byref(event)):
            if event.type == sdl2.SDL_QUIT:
                ctrl[0] = 1
            elif event.type == sdl2.SDL_KEYDOWN:
                k = event.key.keysym.sym
                if k in (sdl2.SDLK_q, sdl2.SDLK_ESCAPE):
                    ctrl[0] = 1
            elif event.type == sdl2.SDL_MOUSEBUTTONDOWN:
                _wid = event.button.windowID
                _mx  = event.button.x
                _my  = event.button.y
                if _wid == main_window_id:
                    # Position magnifier centered on click point
                    wx_c, wy_c = ctypes.c_int(0), ctypes.c_int(0)
                    sdl2.SDL_GetWindowPosition(window_p,
                        ctypes.byref(wx_c), ctypes.byref(wy_c))
                    mag_x = wx_c.value + _mx - MAG_W // 2
                    mag_y = wy_c.value + _my - MAG_H // 2
                    if mag_window_p is None:
                        mw = sdl2.SDL_CreateWindow(
                            b"mag", mag_x, mag_y, MAG_W, MAG_H,
                            sdl2.SDL_WINDOW_SHOWN)
                        if mw:
                            ms = sdl2.SDL_GetWindowSurface(mw)
                            if ms:
                                sdl2.SDL_SetSurfaceBlendMode(
                                    ms, sdl2.SDL_BLENDMODE_NONE)
                                msurf  = ms.contents
                                mp_i32 = msurf.pitch // 4
                                mp_ptr = ctypes.cast(msurf.pixels,
                                    ctypes.POINTER(ctypes.c_int32))
                                md_flat = np.ctypeslib.as_array(
                                    mp_ptr, shape=(MAG_H * mp_i32,))
                                mag_dst       = md_flat.reshape(MAG_H, mp_i32)
                                mag_window_p  = mw
                                mag_surface_p = ms
                            else:
                                sdl2.SDL_DestroyWindow(mw)
                    else:
                        sdl2.SDL_SetWindowPosition(mag_window_p,
                            mag_x, mag_y)
                elif (eg_act_window_id and _wid == eg_act_window_id):
                    pending_click = ('ega', _mx, _my)
                elif (ep_window_id and _wid == ep_window_id):
                    pending_click = ('ep', _mx, _my)
            elif event.type == sdl2.SDL_WINDOWEVENT:
                if event.window.event == sdl2.SDL_WINDOWEVENT_CLOSE:
                    if (mag_window_p is not None and
                            event.window.windowID == sdl2.SDL_GetWindowID(
                                mag_window_p)):
                        sdl2.SDL_DestroyWindow(mag_window_p)
                        mag_window_p  = None
                        mag_surface_p = None
                        mag_dst       = None
                    elif event.window.windowID == main_window_id:
                        ctrl[0] = 1

        if ctrl[0]:
            break

        # ── Deferred probe click processing ──────────────────────
        if pending_click is not None:
            _kind, _cx, _cy = pending_click
            pending_click = None
            try:
                if _kind == 'ega' and eg_activity_pixels is not None \
                        and ega_m_cols is not None:
                    cur = int(eg_activity_cursor[0])
                    bx = (_cx + cur) % PROBE_W
                    if 0 <= _cy < ACT_H and 0 <= bx < PROBE_W:
                        col2eg = {}
                        for i in range(64):
                            cv = int(ega_m_cols[i])
                            col2eg[cv] = i
                            ru = (cv >> 16) & 0xFF
                            gu = (cv >>  8) & 0xFF
                            bu =  cv        & 0xFF
                            dim = (0xFF000000 | ((ru*15//100) << 16)
                                   | ((gu*15//100) << 8) | (bu*15//100))
                            if dim >= 0x80000000:
                                dim -= 0x100000000
                            col2eg.setdefault(dim, i)
                        best_eg = -1
                        for sy in range(_cy, -1, -1):
                            spx = int(eg_activity_pixels[sy, bx])
                            if spx in col2eg:
                                best_eg = col2eg[spx]
                                break
                        if best_eg >= 0:
                            eg_overlay_val = best_eg
                            total_pop = max(int(ega_m_pops.sum()), 1)
                            frac = int(ega_m_pops[best_eg]) / total_pop
                            alive = ("alive" if ega_m_pops[best_eg] > 0
                                     else "extinct")
                            print(f"eg_activity click: egenome "
                                  f"{best_eg} (0b{best_eg:06b})  "
                                  f"pop={frac:.3f}  {alive}",
                                  flush=True)
                elif _kind == 'ep' and eg_pop_pixels is not None \
                        and ega_m_cols is not None:
                    cur = int(eg_pop_cursor[0])
                    bx = (_cx + cur) % PROBE_W
                    if 0 <= _cy < PROBE_H and 0 <= bx < PROBE_W:
                        col2eg = {}
                        for i in range(64):
                            col2eg[int(ega_m_cols[i])] = i
                        best_eg = -1
                        for sy in range(_cy, -1, -1):
                            spx = int(eg_pop_pixels[sy, bx])
                            if spx in col2eg:
                                best_eg = col2eg[spx]
                                break
                        if best_eg >= 0:
                            ep_overlay_val = best_eg
                            total_pop = max(int(ega_m_pops.sum()), 1)
                            frac = int(ega_m_pops[best_eg]) / total_pop
                            print(f"eg_pop click: egenome {best_eg} "
                                  f"(0b{best_eg:06b})  "
                                  f"pop={frac:.3f}", flush=True)
            except Exception:
                with open("/tmp/evoca_click.log", "a") as _lf:
                    traceback.print_exc(file=_lf)

        # Render main window
        sdl2.SDL_LockSurface(surface_p)
        src = pixels.reshape(N, N)
        if px == 1:
            dst[:N, :N] = src
        else:
            dst[:H, :W] = np.repeat(np.repeat(src, px, axis=0), px, axis=1)
        sdl2.SDL_UnlockSurface(surface_p)
        sdl2.SDL_UpdateWindowSurface(window_p)

        # Render magnifier — center computed from window position
        if mag_window_p is not None:
            wx_m, wy_m = ctypes.c_int(0), ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(mag_window_p,
                ctypes.byref(wx_m), ctypes.byref(wy_m))
            wx_main, wy_main = ctypes.c_int(0), ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(window_p,
                ctypes.byref(wx_main), ctypes.byref(wy_main))
            # Center of magnifier in main-window pixel coords
            cx_px = wx_m.value + MAG_W // 2 - wx_main.value
            cy_px = wy_m.value + MAG_H // 2 - wy_main.value
            cell_col = max(MAG_HALF, min(N - 1 - MAG_HALF, cx_px // px))
            cell_row = max(MAG_HALF, min(N - 1 - MAG_HALF, cy_px // px))
            r0 = cell_row - MAG_HALF
            c0 = cell_col - MAG_HALF
            region = src[r0:r0+MAG_CELLS, c0:c0+MAG_CELLS]
            sdl2.SDL_LockSurface(mag_surface_p)
            mag_dst[:MAG_H, :MAG_W] = np.repeat(
                np.repeat(region, MAG_PX, axis=0), MAG_PX, axis=1)
            sdl2.SDL_UnlockSurface(mag_surface_p)
            sdl2.SDL_UpdateWindowSurface(mag_window_p)
            sdl2.SDL_SetWindowTitle(mag_window_p,
                f"mag ({cell_row},{cell_col})".encode())

        # Render probe windows
        if probe_cursor is not None:
            cur = int(probe_cursor[0])
            for i in range(len(probe_windows)):
                sdl2.SDL_LockSurface(probe_surfaces[i])
                cm, cb = probe_colors[i]
                _render_probe(probe_dsts[i], probe_dsts[i].shape[1],
                              probe_means[i], probe_stds[i], cur, cm, cb,
                              y_fixed=probe_yfix[i])
                sdl2.SDL_UnlockSurface(probe_surfaces[i])
                sdl2.SDL_UpdateWindowSurface(probe_windows[i])

        # Render activity window (scroll so newest column is on the right)
        if act_window_p is not None and activity_pixels is not None:
            sdl2.SDL_LockSurface(act_surface_p)
            cur_act = int(activity_cursor[0])
            act_dst[:ACT_H, :PROBE_W] = np.roll(activity_pixels, -cur_act,
                                                  axis=1)
            sdl2.SDL_UnlockSurface(act_surface_p)
            sdl2.SDL_UpdateWindowSurface(act_window_p)

        # Render egenome activity window
        if eg_act_window_p is not None and eg_activity_pixels is not None:
            sdl2.SDL_LockSurface(eg_act_surface_p)
            cur_ega = int(eg_activity_cursor[0])
            eg_act_dst[:ACT_H, :PROBE_W] = np.roll(eg_activity_pixels,
                                                     -cur_ega, axis=1)
            # Overlay: draw 5×5 egenome pattern in top-right corner
            if eg_overlay_val >= 0:
                _orbit = np.array([[4,5,2,5,4],[5,3,1,3,5],[2,1,0,1,2],
                                   [5,3,1,3,5],[4,5,2,5,4]], dtype=np.uint8)
                pat = ((eg_overlay_val >> _orbit) & 1).astype(np.uint8)
                col = int(ega_m_cols[eg_overlay_val]) if ega_m_cols is not None \
                    else _c(0xFFFFFFFF)
                ox = PROBE_W - EG_OVL_PX - 3
                oy = 3
                # border (1px white)
                eg_act_dst[oy-1:oy+EG_OVL_PX+1, ox-1:ox+EG_OVL_PX+1] = \
                    _c(0xFFFFFFFF)
                for r in range(EG_PAT_N):
                    for c in range(EG_PAT_N):
                        px_col = col if pat[r, c] else _c(0xFF000000)
                        y0 = oy + r * EG_CELL_PX
                        x0 = ox + c * EG_CELL_PX
                        eg_act_dst[y0:y0+EG_CELL_PX, x0:x0+EG_CELL_PX] = px_col
            sdl2.SDL_UnlockSurface(eg_act_surface_p)
            sdl2.SDL_UpdateWindowSurface(eg_act_window_p)

        # Render LUT complexity window
        if lc_window_p is not None and lut_complexity_pixels is not None:
            sdl2.SDL_LockSurface(lc_surface_p)
            cur_lc = int(lut_complexity_cursor[0])
            lc_dst[:PROBE_H, :PROBE_W] = np.roll(lut_complexity_pixels,
                                                   -cur_lc, axis=1)
            sdl2.SDL_UnlockSurface(lc_surface_p)
            sdl2.SDL_UpdateWindowSurface(lc_window_p)

        # Render egenome population window
        if ep_window_p is not None and eg_pop_pixels is not None:
            sdl2.SDL_LockSurface(ep_surface_p)
            cur_ep = int(eg_pop_cursor[0])
            ep_dst[:PROBE_H, :PROBE_W] = np.roll(eg_pop_pixels,
                                                    -cur_ep, axis=1)
            # Overlay: draw 5×5 egenome pattern in top-right corner
            if ep_overlay_val >= 0:
                _orbit = np.array([[4,5,2,5,4],[5,3,1,3,5],[2,1,0,1,2],
                                   [5,3,1,3,5],[4,5,2,5,4]], dtype=np.uint8)
                pat = ((ep_overlay_val >> _orbit) & 1).astype(np.uint8)
                col = int(ega_m_cols[ep_overlay_val]) if ega_m_cols is not None \
                    else _c(0xFFFFFFFF)
                ox = PROBE_W - EG_OVL_PX - 3
                oy = 3
                ep_dst[oy-1:oy+EG_OVL_PX+1, ox-1:ox+EG_OVL_PX+1] = \
                    _c(0xFFFFFFFF)
                for r in range(EG_PAT_N):
                    for c in range(EG_PAT_N):
                        px_col = col if pat[r, c] else _c(0xFF000000)
                        y0 = oy + r * EG_CELL_PX
                        x0 = ox + c * EG_CELL_PX
                        ep_dst[y0:y0+EG_CELL_PX, x0:x0+EG_CELL_PX] = px_col
            sdl2.SDL_UnlockSurface(ep_surface_p)
            sdl2.SDL_UpdateWindowSurface(ep_window_p)

        # Render entropy window
        if ent_window_p is not None and entropy_buf is not None:
            sdl2.SDL_LockSurface(ent_surface_p)
            cur_ent = int(entropy_cursor[0])
            _render_probe(ent_dst, ent_dst.shape[1],
                          entropy_buf, entropy_std_buf, cur_ent,
                          _c(0xFF00CCCC), _c(0xFF004444))
            sdl2.SDL_UnlockSurface(ent_surface_p)
            sdl2.SDL_UpdateWindowSurface(ent_window_p)

        # Render pattern activity window
        if pa_window_p is not None and pat_activity_pixels is not None:
            sdl2.SDL_LockSurface(pa_surface_p)
            cur_pa = int(pat_activity_cursor[0])
            pa_dst[:ACT_H, :PROBE_W] = np.roll(pat_activity_pixels,
                                                -cur_pa, axis=1)
            sdl2.SDL_UnlockSurface(pa_surface_p)
            sdl2.SDL_UpdateWindowSurface(pa_window_p)

        # Render q_activity window
        if qa_window_p is not None and q_activity_deciles is not None:
            sdl2.SDL_LockSurface(qa_surface_p)
            cur_qa = int(q_activity_cursor[0])
            qa_global_max = _render_q_activity(qa_dst, q_activity_deciles,
                                                cur_qa, qa_global_max)
            sdl2.SDL_UnlockSurface(qa_surface_p)
            sdl2.SDL_UpdateWindowSurface(qa_window_p)

        # Render N-activity window (Channon shadow)
        if n_window_p is not None and n_activity_pixels is not None:
            sdl2.SDL_LockSurface(n_surface_p)
            cur_n = int(n_activity_cursor[0])
            n_dst[:ACT_H, :PROBE_W] = np.roll(n_activity_pixels, -cur_n,
                                                axis=1)
            sdl2.SDL_UnlockSurface(n_surface_p)
            sdl2.SDL_UpdateWindowSurface(n_window_p)

        # Render Nq-activity window (shadow deciles)
        if nq_window_p is not None and nq_activity_deciles is not None:
            sdl2.SDL_LockSurface(nq_surface_p)
            cur_nq = int(nq_activity_cursor[0])
            nq_global_max = _render_q_activity(nq_dst, nq_activity_deciles,
                                                cur_nq, nq_global_max)
            sdl2.SDL_UnlockSurface(nq_surface_p)
            sdl2.SDL_UpdateWindowSurface(nq_window_p)

        # Render ts window
        if ts_window_p is not None and ts_bufs is not None:
            sdl2.SDL_LockSurface(ts_surface_p)
            cur_ts = int(ts_cursor[0])
            _render_ts(ts_dst, ts_bufs, cur_ts)
            sdl2.SDL_UnlockSurface(ts_surface_p)
            sdl2.SDL_UpdateWindowSurface(ts_window_p)

        # Window title
        step   = int(ctrl[2])
        mode   = COLOR_MODES[min(int(ctrl[1]), len(COLOR_MODES) - 1)]
        paused = bool(ctrl[4])
        if paused:
            title = f"EvoCA  PAUSED  t={step}  color={mode}"
            sdl2.SDL_Delay(16)
        else:
            fps   = ctrl[3] / 10.0
            title = f"EvoCA  t={step}  fps={fps:.1f}  color={mode}"
            sdl2.SDL_Delay(4)   # cap render ~250 fps; prevents 100% CPU spin
        sdl2.SDL_SetWindowTitle(window_p, title.encode())

    print("EvoCA SDL: exiting cleanly", flush=True)
    if ts_window_p is not None:
        sdl2.SDL_DestroyWindow(ts_window_p)
    if mag_window_p is not None:
        sdl2.SDL_DestroyWindow(mag_window_p)
    if pa_window_p is not None:
        sdl2.SDL_DestroyWindow(pa_window_p)
    if ent_window_p is not None:
        sdl2.SDL_DestroyWindow(ent_window_p)
    if nq_window_p is not None:
        sdl2.SDL_DestroyWindow(nq_window_p)
    if n_window_p is not None:
        sdl2.SDL_DestroyWindow(n_window_p)
    if qa_window_p is not None:
        sdl2.SDL_DestroyWindow(qa_window_p)
    if ep_window_p is not None:
        sdl2.SDL_DestroyWindow(ep_window_p)
    if lc_window_p is not None:
        sdl2.SDL_DestroyWindow(lc_window_p)
    if eg_act_window_p is not None:
        sdl2.SDL_DestroyWindow(eg_act_window_p)
    if act_window_p is not None:
        sdl2.SDL_DestroyWindow(act_window_p)
    for pw in probe_windows:
        sdl2.SDL_DestroyWindow(pw)
    sdl2.SDL_DestroyWindow(window_p)
    sdl2.SDL_Quit()
    pixel_shm.close()
    ctrl_shm.close()
    if probe_shm is not None:
        probe_shm.close()
    if activity_shm is not None:
        activity_shm.close()
    if eg_activity_shm is not None:
        eg_activity_shm.close()
    if lut_complexity_shm is not None:
        lut_complexity_shm.close()
    if entropy_shm is not None:
        entropy_shm.close()
    if pat_activity_shm is not None:
        pat_activity_shm.close()
    if eg_pop_shm is not None:
        eg_pop_shm.close()
    if q_activity_shm is not None:
        q_activity_shm.close()
    if n_activity_shm is not None:
        n_activity_shm.close()
    if nq_activity_shm is not None:
        nq_activity_shm.close()
    if ts_shm is not None:
        ts_shm.close()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
