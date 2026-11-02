"""N's question (2026-11-02): how to widen the spacing between fronts in
2026-11-02_spiral_wave_ish? Sweep food_inc x gdiff x waste_sink on the
recipe (mutation on, as exported), N = 512.

Spacing estimate: front spacing ~ speed x period. Readouts:
  scale    first zero of the row autocorrelation of F (row means removed),
           both axes, last frames: the half-size of food patches
  peak     autocorrelation at its first maximum after that zero (> 0 only
           for a regular wave train) and its lag
  period   site_period of the W kymograph (middle row)

    python3 s6_wavelength.py <i>        # one config -> out/s6/<i>.json + frame png
    python3 s6_wavelength.py report
"""
import itertools
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', '..', 'python'))
import evoca_nb  # noqa: E402
import evoca_waves as ew  # noqa: E402

GRID = dict(food_inc=[0.004, 0.006, 0.008, 0.012],
            gdiff=[0.02, 0.08, 0.3, 1.0],
            waste_sink=[0.01, 0.019])
CONFIGS = [dict(zip(GRID, v)) for v in itertools.product(*GRID.values())]
N, STEPS, KYMO = 512, 6000, 2000


def autocorr_scale(frames):
    acs = []
    for F in frames:
        for X in (F, F.T):
            Y = X - X.mean(axis=1, keepdims=True)
            acs.append(np.fft.irfft(np.abs(np.fft.rfft(Y, axis=1)) ** 2, axis=1).sum(axis=0))
    ac = np.mean(acs, axis=0)
    ac = ac / ac[0]
    n = len(ac) // 2
    z = int(np.argmax(ac[:n] < 0)) if (ac[:n] < 0).any() else n
    j = z + int(np.argmax(ac[z:n])) if z < n else n
    return z, j, float(ac[j]) if j < len(ac) else float('nan')


def one(i):
    cfg = CONFIGS[i]
    sim, _ = evoca_nb.load_recipe('2026-11-02_spiral_wave_ish.evoca', N=N, seed=0, **cfg)
    kym, frames, alive = [], [], []
    for t in range(1, STEPS + 1):
        sim.step()
        if t > STEPS - KYMO:
            kym.append(sim.get_W()[N // 2].copy())
            if t % 400 == 0:
                frames.append(sim.get_F())
        if t % 50 == 0:
            alive.append(sim.population() / N ** 2)
    rgb = evoca_nb.rgb(sim, 1)
    sim.free()
    z, j, pk = autocorr_scale(frames)
    per, acf = ew.site_period(np.array(kym), b=4)
    a = np.array(alive)
    row = dict(i=i, **cfg, alive=float(a[len(a) // 2:].mean()), extinct=bool(a[-1] == 0),
               scale=z, peak_lag=j, peak=pk, period=per, acf=acf)
    os.makedirs(os.path.join(HERE, 'out', 's6'), exist_ok=True)
    json.dump(row, open(os.path.join(HERE, 'out', 's6', f'{i:02d}.json'), 'w'))
    np.save(os.path.join(HERE, 'out', 's6', f'{i:02d}_rgb.npy'), rgb)
    print(row, flush=True)


def report():
    import glob
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(HERE, 'out', 's6', '*.json')))]
    d = pd.DataFrame(rows)
    pd.set_option('display.width', 200)
    print(d.round(3).to_string(index=False))
    for ws in GRID['waste_sink']:
        fig, axs = plt.subplots(len(GRID['food_inc']), len(GRID['gdiff']),
                                figsize=(3 * len(GRID['gdiff']), 3 * len(GRID['food_inc'])))
        for _, r in d[d.waste_sink == ws].iterrows():
            a = axs[GRID['food_inc'].index(r.food_inc), GRID['gdiff'].index(r.gdiff)]
            a.imshow(np.load(os.path.join(HERE, 'out', 's6', f'{int(r.i):02d}_rgb.npy')), interpolation='nearest')
            a.set_title(f"inc {r.food_inc} gdiff {r.gdiff}\nscale {r.scale}, T {r.period:.0f}, alive {r.alive:.2f}",
                        fontsize=8)
            a.axis('off')
        fig.suptitle(f'spiral_wave_ish, N=512, t=6000, waste_sink={ws}')
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, 'figs', f's6_sweep_wsink{ws}.png'), dpi=70)
        plt.close(fig)


if __name__ == '__main__':
    report() if sys.argv[1] == 'report' else one(int(sys.argv[1]))
