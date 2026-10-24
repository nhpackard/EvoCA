import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, sys
from s2_1d_scan import run1d
cases = [(57.7,0.024),(86.8,0.092),(57.7,0.143),(25.6,0.346),(25.6,0.538),(17,0.838),(17,1.3),(11.3,2.0),(38.4,0.222),(130,0.092),(295,0.059),(17,3.16)]
fig, axs = plt.subplots(3, 4, figsize=(20, 12))
for ax, (a, rho) in zip(axs.flat, cases):
    u, v, p, h = run1d(a, rho, L=400, T=1600, rec=True)
    ax.imshow(h, aspect='auto', origin='lower', extent=[0, 400, 0, 1600], cmap='magma')
    ax.set_title(f'v(x,t): a={a}, rho={rho}'); ax.set_xlabel('x'); ax.set_ylabel('step')
plt.tight_layout(); plt.savefig('figs/f2b_1d_gallery.png', dpi=70)
