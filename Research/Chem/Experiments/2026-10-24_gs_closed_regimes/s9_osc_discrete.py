"""Discrete-time well-mixed oscillation map (k2=0.1) and 2D oscillatory medium (noisy uniform start)."""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from multiprocessing import Pool
from gslib import *
from gs2d import make_stepper

def wm(args):
    a, rho, k2 = args
    best = (0.0, np.nan)
    for y0 in [(0.6, 0.4, 0.0), (0.9, 0.05, 0.05)]:
        pts = fixed_points(a, rho)
        ics = [y0] + ([ (pts[1][0], pts[1][1]*1.05, pts[1][2]-0.05*pts[1][1]) ] if len(pts) == 3 else [])
        for ic in ics:
            tr = well_mixed_discrete(a*k2, k2, rho*k2, ic, 12000)[:, 1][-4000:]
            if tr.min() > 1e-3 and tr.max()-tr.min() > 1e-3:
                idx, _ = find_peaks(tr)
                if len(idx) >= 3 and tr.max()-tr.min() > best[0]: best = (tr.max()-tr.min(), np.median(np.diff(idx)))
    return best

if __name__ == '__main__':
    k2 = 0.1
    Ag = np.logspace(1.2, 3.2, 30); Rg = np.logspace(-2.0, -0.3, 30)
    with Pool(12) as P: out = P.map(wm, [(a, r, k2) for r in Rg for a in Ag])
    amp = np.array([o[0] for o in out]).reshape(30, 30); per = np.array([o[1] for o in out]).reshape(30, 30)
    np.savez('s9_wm.npz', Ag=Ag, Rg=Rg, amp=amp, per=per)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    q = ax[0].pcolormesh(Ag, Rg, np.where(amp > 0, amp, np.nan), shading='auto'); plt.colorbar(q, ax=ax[0]); ax[0].plot(4*(1+1/Rg), Rg, 'k--', lw=1)
    ax[0].set(xscale='log', yscale='log', title='discrete-time well-mixed (k2=0.1): limit-cycle amplitude of v (blank = none)\ndashed: saddle-node a=4(1+1/rho)', xlabel='a', ylabel='rho')
    q = ax[1].pcolormesh(Ag, Rg, per, shading='auto'); plt.colorbar(q, ax=ax[1]); ax[1].set(xscale='log', yscale='log', title='period (steps)', xlabel='a')
    plt.tight_layout(); plt.savefig('figs/f9a_discrete_wm_cycles.png', dpi=90); plt.close()
    ii = np.argwhere(amp > 0)
    for i, j in ii[::max(1, len(ii)//12)]: print(f'a={Ag[j]:.1f} rho={Rg[i]:.3f} amp={amp[i,j]:.3f} period={per[i,j]:.0f}')
    # 2D oscillatory medium
    N = 128; rng = np.random.default_rng(3)
    fig, axs = plt.subplots(2, 4, figsize=(18, 9))
    for row, (a, rho, d) in zip(axs, [(141.5, 0.094, (0.6, 0.3, 0.3)), (141.5, 0.094, (0.2, 0.1, 0.1))]):
        pts = fixed_points(a, rho); uu, vv, pp = pts[1]
        u = np.full((N, N), uu); v = vv*(1+0.2*rng.standard_normal((N, N))); p = np.full((N, N), pp); u = 1-v-p
        st = make_stepper((N, N)); k1 = a*k2; r = rho*k2; mean = []
        snaps = {0: None}
        t = 0
        for tgt in (150, 400, 1500, 4000):
            st(u, v, p, k1, k2, r, *d, nsteps=tgt-t); t = tgt
            row[[150, 400, 1500, 4000].index(tgt)].imshow(v, cmap='magma'); row[[150, 400, 1500, 4000].index(tgt)].set_title(f'a={a} rho={rho} d={d} t={tgt}\nv mean {v.mean():.3f} std {v.std():.3f}', fontsize=9)
    plt.tight_layout(); plt.savefig('figs/f9b_2d_oscillatory_medium.png', dpi=60)
