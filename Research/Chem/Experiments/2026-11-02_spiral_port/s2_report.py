"""Rank the s2 ChemCA port scan against the EvoCA target and plot the best."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
# EvoCA target: mean of base_N256 and nomut_N256 (s1)
TARGET = dict(alive=0.22, F=0.53, W=0.19, period=146.0)


def score(r):
    """Distance to the EvoCA target: log ratios for alive and period,
    absolute differences for F and W (both O(0.1-1)); inf without a period."""
    if not np.isfinite(r.period_W):
        return np.inf
    return (abs(np.log(r.alive / TARGET['alive'])) + abs(r.F - TARGET['F'])
            + abs(r.W - TARGET['W']) / TARGET['W'] + abs(np.log(r.period_W / TARGET['period'])))


d = pd.read_csv(os.path.join(HERE, 'out', 's2_results.csv'))
v = d[~d.extinct].copy()
v['score'] = v.apply(score, axis=1)
v = v.sort_values('score')
pd.set_option('display.width', 250)
cols = ['i', 'u_F', 'tox_W', 'u_W', 'sink_W', 'sink_F', 'alive', 'F', 'W', 'w_alive', 'tox_frac',
        'period_W', 'acf_W', 'score']
print(v[cols].head(15).round(3).to_string(index=False))
top = v.head(8)
fig, axs = plt.subplots(4, len(top), figsize=(2.7 * len(top), 11))
for j, (_, r) in enumerate(top.iterrows()):
    z = np.load(os.path.join(HERE, 'out', 's2', f'{int(r.i):03d}.npz'))
    axs[0, j].imshow(z['fr_rgb'][-1], interpolation='nearest')
    axs[0, j].set_title(f"#{int(r.i)} uF={r.u_F} tox={r.tox_W}\nuW={r.u_W} sW={r.sink_W} sF={r.sink_F}", fontsize=8)
    axs[1, j].imshow(z['fr_F'][-1], vmin=0, vmax=1, cmap='Greens', interpolation='nearest')
    axs[2, j].imshow(z['fr_W'][-1], vmin=0, vmax=max(0.05, float(z['fr_W'][-1].max())), cmap='Purples',
                     interpolation='nearest')
    axs[3, j].imshow(z['kym_alive'], aspect='auto', cmap='gray', interpolation='nearest')
    axs[3, j].set_title(f"period {r.period_W:.0f}, alive {r.alive:.2f}", fontsize=8)
    for a in axs[:3, j]:
        a.axis('off')
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'figs', 's2_top8.png'), dpi=70)
v.to_csv(os.path.join(HERE, 'out', 's2_ranked.csv'), index=False)
