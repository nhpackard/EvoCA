"""Final figures for recommended sets: well-mixed time series, 1D space-time, 2D spiral (u, v, p fields)."""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from gslib import *
from gs2d import make_stepper
from s4_2d import seed_broken
k2 = 0.1; d = (0.6, 0.3, 0.3)
SETS = {'S1 oscillatory': (141.5, 0.094), 'S2 excitable pulse': (86.8, 0.092), 'S3 spiral': (60, 0.05), 'S4 invasion front': (25.6, 0.346)}
if __name__ == '__main__':
    fig, ax = plt.subplots(2, 4, figsize=(20, 8))
    for k, (nm, (a, rho)) in enumerate(SETS.items()):
        wm = well_mixed_discrete(a*k2, k2, rho*k2, (0.6, 0.4, 0.0), 1500)
        ax[0, k].plot(wm[:, 0], label='u'); ax[0, k].plot(wm[:, 1], label='v'); ax[0, k].plot(wm[:, 2], label='p'); ax[0, k].legend()
        ax[0, k].set(title=f'{nm}: well-mixed discrete\nk1={a*k2:.3g}, k2={k2}, r={rho*k2:.3g}', xlabel='step')
        L = 400; u = np.ones(L); v = np.zeros(L); p = np.zeros(L); c = L//2
        v[c-8:c+8] = 0.4; u[c-8:c+8] = 0.6
        H = []
        for t in range(1600):
            u, v, p = step_discrete(u, v, p, a*k2, k2, rho*k2, *d)
            if t % 2 == 0: H.append(v.copy())
        ax[1, k].imshow(np.array(H), aspect='auto', origin='lower', extent=[0, L, 0, 1600], cmap='magma'); ax[1, k].set(title='1D v(x,t), d=(0.6,0.3,0.3)', xlabel='x', ylabel='step')
    plt.tight_layout(); plt.savefig('figs/f11a_recommended_1d.png', dpi=75); plt.close()
    # front speed
    a, rho = SETS['S4 invasion front']; L = 800; u = np.ones(L); v = np.zeros(L); p = np.zeros(L); c = L//2
    v[c-8:c+8] = 0.4; u[c-8:c+8] = 0.6; pos = {}
    for t in range(1500):
        u, v, p = step_discrete(u, v, p, a*k2, k2, rho*k2, *d)
        if t in (700, 1400): pos[t] = np.max(np.nonzero(v > 0.02)[0])
    print('invasion front speed (sites/step):', (pos[1400]-pos[700])/700, 'plateau v', v[c])
    # spiral 2D
    N = 256; a, rho = SETS['S3 spiral']; u, v, p = seed_broken(N); st = make_stepper((N, N)); t = 0
    fig, axs = plt.subplots(1, 5, figsize=(25, 5))
    for ax_, tgt in zip(axs[:2], (1500, 6000)):
        st(u, v, p, a*k2, k2, rho*k2, *d, nsteps=tgt-t); t = tgt
        ax_.imshow(v, cmap='magma'); ax_.set_title(f'S3 spiral: v, t={tgt} (a=60, rho=0.05, k2=0.1)'); ax_.axis('off')
    for ax_, f, nm in zip(axs[2:], (u, p, u+v+p), ('u (food)', 'p (waste)', 'u+v+p (mass; =1)')):
        q = ax_.imshow(f, cmap='viridis'); plt.colorbar(q, ax=ax_, fraction=0.046); ax_.set_title(nm + ' at t=6000'); ax_.axis('off')
    plt.tight_layout(); plt.savefig('figs/f11b_spiral_S3.png', dpi=65)
    print('mass check', (u+v+p).mean(), (u+v+p).min(), (u+v+p).max())
