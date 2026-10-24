"""CHEM-N6 analysis of the saved frame stacks (scan.py): layered-front
statistics per run (fronts.csv), seed means per configuration (printed),
and two figures:
  spectra.png  radial spectrum of the alive field (time-averaged over the 51
               frames), one panel per configuration, seeds as separate lines
  frames.png   the colormode-6 frame at t=6000, seed 0 of each configuration
  kymo_view.png  space-time plots (middle row, steps 5400-6000, x 0-511) of
               alive (9-cell running mean) and environmental W, for the wave runs
Frames are 4 x 4 block means, so ring_spectrum runs with b=1 and wavelengths
are converted back to cells."""
import csv
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
while not os.path.exists(os.path.join(ROOT, 'C', 'evoca.h')):
    ROOT = os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(ROOT, 'python'))
import evoca_waves as ew                                     # noqa: E402

PAT = re.compile(r'(\w+)_N(\d+)_ws([\d.]+)_wd([\d.]+)_tox([\d.]+)_x(\d)_s(\d)\.npz')
RING = ('wavelength', 'prominence', 'contrast', 'layers', 'isotropy', 'moving')
KW = ('moving', 'wavelength', 'speed', 'peak_share')


def config_of(f):
    g, n, ws, wd, tox, x, s = PAT.match(os.path.basename(f)).groups()
    return (g, int(n), float(ws), float(wd), float(tox), int(x)), int(s)


def label(c):
    g, n, ws, wd, tox, x = c
    return f"{'two-way' if x else 'one-way'} N={n}\nsink {ws}, wdiff {wd}, tox {tox}"


def main():
    rows, spectra = [], {}
    for f in sorted(glob.glob(os.path.join(HERE, 'frames', '*.npz'))):
        c, s = config_of(f)
        d = np.load(f)
        b = int(d['block'])
        row = dict(zip(('group', 'N', 'waste_sink', 'wdiff', 'waste_tox', 'waste_exchange'), c), seed=s,
                   alive_mean=float(d['alive_series'][-10:].mean()), n_frames=len(d['alive']))
        for name in ('alive', 'W'):
            fr = d[name].astype(float)
            if len(fr) < 2:
                continue
            r = ew.ring_spectrum(fr, b=1)
            r['wavelength'] *= b
            row.update({f'ring_{name}_{k}': r[k] for k in RING})
            k = ew.kw_spectrum(d['kymo_' + name].astype(float), b=4)
            row.update({f'kw_{name}_{q}': k[q] for q in KW})
            row[f'period_{name}'], row[f'period_{name}_peak'] = ew.site_period(d['kymo_' + name].astype(float))
            if name == 'alive':
                spectra.setdefault(c, []).append(r['S'])
        rows.append(row)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(HERE, 'fronts.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    import pandas as pd
    df = pd.DataFrame(rows)
    cols = ['period_alive', 'period_alive_peak', 'period_W', 'period_W_peak',
            'ring_alive_wavelength', 'ring_alive_contrast', 'ring_alive_prominence', 'ring_alive_layers',
            'ring_alive_isotropy', 'ring_alive_moving', 'ring_W_contrast', 'ring_W_moving',
            'kw_alive_peak_share', 'kw_alive_speed', 'alive_mean']
    g = df.groupby(['waste_exchange', 'N', 'waste_sink', 'wdiff', 'waste_tox'])[cols]
    pd.set_option('display.width', 250)
    print(g.mean().round(2).to_string())
    print('\nseed spread (max - min):')
    print((g.max() - g.min()).round(2).to_string())
    figures(spectra)


def figures(spectra):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cfgs = sorted(spectra, key=lambda c: (c[5], c[1], c[3], c[2]))
    nc = 4
    nr = -(-len(cfgs) // nc)
    fig, axes = plt.subplots(nr, nc, figsize=(3.2 * nc, 2.6 * nr), sharey=False, squeeze=False)
    for ax, c in zip(axes.flat, cfgs):
        N = c[1]
        for S, shade in zip(spectra[c], ('#2a78d6', '#5598e7', '#9ec5f4')):
            k = np.arange(1, len(S))
            ax.plot(N / k, S[1:] / np.median(S[1:len(S) // 2]), color=shade, lw=1.5)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_title(label(c), fontsize=8)
        ax.set_xlabel('wavelength (cells)', fontsize=7)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.25, lw=0.5)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    for ax in axes.flat[len(cfgs):]:
        ax.axis('off')
    axes[0, 0].set_ylabel('S / median S', fontsize=7)
    fig.suptitle('CHEM-N6: radial spectrum of alive density (seeds 0-2, darkest = seed 0); '
                 'a ring is a hump above the low-k end', fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'spectra.png'), dpi=130)
    plt.close(fig)
    lasts = {}
    for f in sorted(glob.glob(os.path.join(HERE, 'frames', '*_s0.npz'))):
        lasts[config_of(f)[0]] = np.load(f)['last']
    fig, axes = plt.subplots(nr, nc, figsize=(3.2 * nc, 3.4 * nr), squeeze=False)
    for ax, c in zip(axes.flat, cfgs):
        ax.imshow(lasts[c], interpolation='nearest')
        ax.set_title(label(c), fontsize=8)
        ax.axis('off')
    for ax in axes.flat[len(cfgs):]:
        ax.axis('off')
    fig.suptitle('CHEM-N6: colormode 6 at t = 6000, seed 0', fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'frames.png'), dpi=80)
    plt.close(fig)
    kymo_figure(plt)


KYMO_RUNS = ['oneway_N512_ws0.019_wd0.01_tox0.4_x0_s0', 'oneway_N1024_ws0.019_wd0.01_tox0.4_x0_s0',
             'twoway_N1024_ws0.019_wd0.03_tox0.6_x1_s0', 'twoway_N1024_ws0.015_wd0.08_tox1.0_x1_s0']


def kymo_figure(plt):
    fig, ax = plt.subplots(2, len(KYMO_RUNS), figsize=(4 * len(KYMO_RUNS), 8), squeeze=False)
    for j, r in enumerate(KYMO_RUNS):
        d = np.load(os.path.join(HERE, 'frames', r + '.npz'))
        ka = d['kymo_alive'].astype(float)
        ks = np.apply_along_axis(lambda x: np.convolve(np.r_[x[-4:], x, x[:4]], np.ones(9) / 9, 'valid'), 1, ka)
        ax[0, j].imshow(ks[:, :512], aspect='auto', cmap='Greys', interpolation='nearest')
        ax[0, j].set_title(label(config_of(os.path.join('frames', r + '.npz'))[0]) + '\nalive (9-cell mean)', fontsize=8)
        ax[1, j].imshow(d['kymo_W'].astype(float)[:, :512], aspect='auto', cmap='Blues', interpolation='nearest')
        ax[1, j].set_title('environmental W', fontsize=8)
        for a in ax[:, j]:
            a.set_xlabel('x (cells)', fontsize=7)
            a.set_ylabel('step - 5400', fontsize=7)
            a.tick_params(labelsize=7)
    fig.suptitle('CHEM-N6: space-time plots, seed 0. Fronts leave sources (V), travel at 1-2 cells/step (by eye) '
                 'and annihilate on collision (inverted V): excitation waves', fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'kymo_view.png'), dpi=90)
    plt.close(fig)


if __name__ == '__main__':
    main()
