"""Task 2c: 2D discrete-time simulations with broken-front seeds (spiral test)."""
import numpy as np, matplotlib, sys, time; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from gs2d import make_stepper

def seed_broken(N, m=1.0, kick=0.4, refr=0.5, width=6, back=30):
    u = np.full((N, N), m); v = np.zeros((N, N)); p = np.zeros((N, N)); c = N//2
    # axis0 = y, axis1 = x. front along y in lower half
    v[:c, c:c+width] = kick*m; u[:c, c:c+width] -= kick*m
    p[:c, c-back:c] = refr*m; u[:c, c-back:c] -= refr*m
    return u, v, p

def run(a, rho, k2=0.1, du=0.6, dv=0.3, dp=0.3, N=256, T=3000, snaps=(300, 800, 1500, 3000), m=1.0, seed=seed_broken, tag=''):
    k1 = a*k2/m**2; r = rho*k2
    u, v, p = seed(N, m)
    st = make_stepper((N, N)); out = {}; t = 0; tl = []
    for target in snaps:
        st(u, v, p, k1, k2, r, du, dv, dp, nsteps=target-t); t = target; out[t] = v.copy()
    return out, (u, v, p)

if __name__ == '__main__':
    cases = [('A', 86.8, 0.092), ('B', 141.5, 0.094), ('C', 57.7, 0.024), ('D', 25.6, 0.346), ('E', 79.4, 0.128), ('F', 252, 0.069)]
    snaps = (200, 600, 1200, 2400)
    fig, axs = plt.subplots(len(cases), len(snaps), figsize=(4*len(snaps), 4*len(cases)))
    for row, (nm, a, rho) in zip(axs, cases):
        t0 = time.time(); out, _ = run(a, rho, snaps=snaps)
        for ax, (t, f) in zip(row, out.items()):
            ax.imshow(f, cmap='magma', vmin=0, vmax=max(0.05, f.max())); ax.set_title(f'{nm}: a={a}, rho={rho}, t={t}', fontsize=9); ax.axis('off')
        print(nm, 'done', round(time.time()-t0, 1), 's', 'vmax end', out[snaps[-1]].max(), flush=True)
    plt.tight_layout(); plt.savefig('figs/f4_2d_broken_front.png', dpi=60)
