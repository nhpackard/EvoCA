"""Task 4: closed two-species food/waste analog.  u + w = m.  Reaction: u + w -> 2w (k u w), w -> u (r w, light).
1D local dynamics: dw = w (k (m-w) - r).  Show: no oscillation (1D flow), and what diffusion can do:
 (a) linear stability of the homogeneous state with D=diag(Du,Dw): max Re lambda(q) over a grid
 (b) 1D simulation: front invasion, final state homogeneous."""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
def eig_scan(k, r, m):
    out = []
    ws = m - r/k
    for (u, w) in [(r/k, ws)] if ws > 0 else [(m, 0.0)]:
        J = np.array([[-k*w, -k*u + r], [k*w, k*u - r]])      # d(u,w)
        for Du in np.logspace(-2, 2, 41):
            for Dw in np.logspace(-2, 2, 41):
                best = max(np.linalg.eigvals(J - q2*np.diag([Du, Dw])).real.max() for q2 in np.logspace(-4, 2, 60))
                out.append((Du, Dw, best, np.linalg.eigvals(J)))
    return out
def sim(k, r, m, Du, Dw, L=400, T=15000, dt=0.02):
    u = np.full(L, m); w = np.zeros(L); w[200:205] = 0.3; u[200:205] -= 0.3
    lap = lambda x: np.roll(x, 1)+np.roll(x, -1)-2*x
    H = []
    for t in range(T):
        R = k*u*w - r*w
        # R = net conversion u -> w (autocatalysis minus light regeneration); mass u+w conserved
        u, w = u + dt*(-R + Du*lap(u)), w + dt*(R + Dw*lap(w))
        if t % 150 == 0: H.append(w.copy())
    return np.array(H), u, w
if __name__ == '__main__':
    k, r, m = 1.0, 0.2, 1.0
    res = eig_scan(k, r, m)
    mx = max(x[2] for x in res); ev = res[0][3]
    print('two-species eigenvalues at homogeneous steady state (k=1,r=0.2,m=1):', ev, '; max Re lambda over all D ratios in [1e-2,1e2]^2 and q:', mx)
    print('any complex eigenvalue for any q,D?', any(np.any(abs(np.linalg.eigvals(np.array([[-k*(m-r/k), -r+r], [k*(m-r/k), 0]])-q2*np.diag([Du, Dw])).imag) > 1e-9) for q2 in np.logspace(-3, 2, 20) for Du in (0.01, 1, 100) for Dw in (0.01, 1, 100)))
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for axx, (Du, Dw) in zip(ax[:2], [(1.0, 1.0), (10.0, 0.1)]):
        H, u, w = sim(k, r, m, Du, Dw)
        axx.imshow(H, aspect='auto', origin='lower', extent=[0, 400, 0, 15000*0.02], cmap='viridis'); axx.set(title=f'2-species closed: w(x,t), Du={Du}, Dw={Dw}\n(Fisher-type invasion front -> homogeneous w*=m-r/k=0.8)', xlabel='x', ylabel='t')
        print(f'Du={Du} Dw={Dw}: final w min/max', w.min(), w.max())
    # local phase line
    ww = np.linspace(0, 1, 200); ax[2].plot(ww, ww*(k*(m-ww)-r)); ax[2].axhline(0, color='k', lw=.5); ax[2].set(xlabel='w', ylabel='dw/dt', title='local dynamics: 1D flow, stable fixed point w*=0.8 (no oscillation possible)')
    plt.tight_layout(); plt.savefig('figs/f6_two_species.png', dpi=90)
