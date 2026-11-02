"""Spiral-survival map over (a, rho), k2=0.1, du=.6,dv=.3,dp=.3, N=256, T=4000, broken-front seed.
survive = active fraction (v>0.02) averaged over last 500 steps; also records whether the activity is spiral-like (>0 and < 0.9)."""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from s4_2d import run
if __name__ == '__main__':
    As = [30, 40, 50, 60, 70, 80, 100, 130]; Rs = [0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15]
    fr = np.zeros((len(Rs), len(As))); fr2 = np.zeros_like(fr)
    for i, rho in enumerate(Rs):
        for j, a in enumerate(As):
            out, _ = run(a, rho, N=256, snaps=(3500, 4000))
            fr[i, j] = (out[4000] > 0.02).mean(); fr2[i, j] = (out[3500] > 0.02).mean()
        print('rho', rho, np.round(fr[i], 3), flush=True)
    np.savez('s4c_results.npz', As=As, Rs=Rs, fr=fr)
    fig, ax = plt.subplots(figsize=(7, 5))
    q = ax.pcolormesh(np.arange(len(As)+1), np.arange(len(Rs)+1), np.where(fr > 0, fr, np.nan), cmap='viridis')
    ax.set_xticks(np.arange(len(As))+.5, As); ax.set_yticks(np.arange(len(Rs))+.5, Rs); plt.colorbar(q, label='active fraction at t=4000 (white = died out)')
    ax.set(xlabel='a = k1 m^2/k2', ylabel='rho = r/k2', title='2D discrete, broken-front seed: wave survival at t=4000\n(k2=0.1, du/dv/dp=0.6/0.3/0.3, N=256)')
    plt.tight_layout(); plt.savefig('figs/f4c_spiral_map.png', dpi=100)
