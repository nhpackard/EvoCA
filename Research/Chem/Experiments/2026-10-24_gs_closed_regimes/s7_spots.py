"""Spot / Turing search in the closed system. Large Du/Dv (>=25 needed for linear Turing instability of the active state).
Seeds: (S) single v-spot in trivial medium; (N) noise around the active fixed point (when it exists)."""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from gs2d import make_stepper
from gslib import fixed_points
def go(a, rho, du, dv, dp, seed, N=160, T=6000, k2=0.1):
    rng = np.random.default_rng(1); k1 = a*k2; r = rho*k2
    if seed == 'S':
        u = np.ones((N, N)); v = np.zeros((N, N)); p = np.zeros((N, N)); c = N//2
        v[c-4:c+4, c-4:c+4] = 0.4; u[c-4:c+4, c-4:c+4] = 0.6
    else:
        pts = fixed_points(a, rho); uu, vv, pp = pts[1]
        u = uu + 0.01*rng.standard_normal((N, N)); v = vv*(1+0.05*rng.standard_normal((N, N))); p = np.full((N, N), pp)
        s = u+v+p; u += 1-s
    st = make_stepper((N, N)); st(u, v, p, k1, k2, r, du, dv, dp, nsteps=T)
    return v, u
if __name__ == '__main__':
    cases = [(30, 3.0, 0.9, 0.03, 0.3, 'N'), (400, 0.2, 0.9, 0.03, 0.3, 'N'), (30, 3.0, 0.9, 0.03, 0.3, 'S'), (60, 1.0, 0.9, 0.05, 0.3, 'S'), (30, 3.0, 0.9, 0.3, 0.3, 'S'), (200, 0.1, 0.9, 0.03, 0.3, 'S')]
    fig, axs = plt.subplots(1, len(cases), figsize=(4*len(cases), 4.4))
    for ax, (a, rho, du, dv, dp, sd) in zip(axs, cases):
        v, u = go(a, rho, du, dv, dp, sd)
        ax.imshow(v, cmap='viridis'); ax.axis('off'); ax.set_title(f'a={a} rho={rho}\ndu/dv/dp={du}/{dv}/{dp} seed {sd}\nv std={v.std():.3g} mean={v.mean():.3g}', fontsize=9)
        print(a, rho, du, dv, dp, sd, 'v mean', v.mean(), 'std', v.std(), 'max', v.max(), flush=True)
    plt.tight_layout(); plt.savefig('figs/f7_spots_search.png', dpi=70)
