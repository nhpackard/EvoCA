"""Spiral search: broken-front seed, N=384, low rho band. Montage of t=1500, 3500, 6000 and active-fraction table."""
import numpy as np, matplotlib, time; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from s4_2d import run
if __name__ == '__main__':
    As = [40, 60, 90, 130]; Rs = [0.015, 0.03, 0.05]
    snaps = (1500, 3500, 6000)
    fig, axs = plt.subplots(len(Rs)*len(As)//2, 6, figsize=(24, 4*len(Rs)*len(As)//2))
    k = 0; table = []
    for rho in Rs:
        for a in As:
            out, _ = run(a, rho, N=384, snaps=snaps)
            row = axs[k//2]; off = 3*(k % 2)
            for j, (t, f) in enumerate(out.items()):
                row[off+j].imshow(f, cmap='magma', vmin=0, vmax=max(0.05, f.max())); row[off+j].axis('off'); row[off+j].set_title(f'a={a} rho={rho} t={t}', fontsize=9)
            fr = (out[6000] > 0.02).mean(); table.append((a, rho, fr, out[6000].max())); print(a, rho, 'active frac', round(fr, 3), 'vmax', round(out[6000].max(), 3), flush=True)
            k += 1
    plt.tight_layout(); plt.savefig('figs/f4b_spiral_scan.png', dpi=45)
