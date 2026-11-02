"""Summary table and figures for the s1 EvoCA characterization runs."""
import os, sys, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', '..', 'python'))
import evoca_waves as ew

os.makedirs(os.path.join(HERE, 'figs'), exist_ok=True)
runs = sorted(glob.glob(os.path.join(HERE, 'out', 's1_*.npz')))
rows = []
for p in runs:
    tag = os.path.basename(p)[3:-4]
    d = np.load(p)
    late = d['ts_t'] > d['ts_t'].max() - 3000
    per, pk = ew.site_period(d['kym_W'], b=4)
    perA, _ = ew.site_period(d['kym_alive'], b=4)
    r = dict(run=tag, alive=d['ts_alive'][late].mean(), tax_pc=np.nanmean(d['ts_tax_pc'][late]),
             intake_pc=np.nanmean(d['ts_intake_pc'][late]), eat_frac=np.nanmean(d['ts_eat_frac'][late]),
             F=d['ts_F'][late].mean(), W=d['ts_W'][late].mean(), F_alive=np.nanmean(d['ts_F_alive'][late]),
             w_alive=np.nanmean(d['ts_w_alive'][late]),
             tox_frac=d['ts_tox'][late].sum() / max(d['ts_deaths'][late].sum() + d['ts_tox'][late].sum(), 1),
             period_W=per, acf=pk, period_alive=perA)
    rows.append(r)
    # frames
    ft = d['frame_t']; sel = [i for i, t in enumerate(ft) if t in (1000, 2000, 3000, 4000, 5000, 6000)]
    fig, axs = plt.subplots(3, len(sel), figsize=(2.6 * len(sel), 8))
    for j, i in enumerate(sel):
        axs[0, j].imshow(d['fr_rgb'][i], interpolation='nearest'); axs[0, j].set_title(f't={ft[i]}', fontsize=9)
        axs[1, j].imshow(d['fr_F'][i], vmin=0, vmax=1, cmap='Greens', interpolation='nearest')
        axs[2, j].imshow(d['fr_W'][i], vmin=0, vmax=max(0.05, d['fr_W'][i].max()), cmap='Purples', interpolation='nearest')
        for a in axs[:, j]: a.axis('off')
    axs[0, 0].text(-0.1, 0.5, 'colormode 1', transform=axs[0, 0].transAxes, rotation=90, va='center', ha='right')
    axs[1, 0].text(-0.1, 0.5, 'F', transform=axs[1, 0].transAxes, rotation=90, va='center', ha='right')
    axs[2, 0].text(-0.1, 0.5, 'W', transform=axs[2, 0].transAxes, rotation=90, va='center', ha='right')
    fig.suptitle(f'EvoCA {tag}'); fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'figs', f's1_frames_{tag}.png'), dpi=80); plt.close(fig)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    axs[0].plot(d['ts_t'], d['ts_alive'], label='alive'); axs[0].plot(d['ts_t'], d['ts_F'], label='mean F')
    axs[0].plot(d['ts_t'], d['ts_W'], label='mean W'); axs[0].legend(); axs[0].set_xlabel('step')
    axs[1].imshow(d['kym_alive'], aspect='auto', cmap='gray', interpolation='nearest'); axs[1].set_title('alive, middle row (last 2000)')
    axs[2].imshow(d['kym_W'], aspect='auto', cmap='Purples', interpolation='nearest'); axs[2].set_title('W, middle row')
    fig.suptitle(f'EvoCA {tag}'); fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'figs', f's1_ts_{tag}.png'), dpi=80); plt.close(fig)
import pandas as pd
pd.set_option('display.width', 200)
print(pd.DataFrame(rows).round(4).to_string(index=False))
