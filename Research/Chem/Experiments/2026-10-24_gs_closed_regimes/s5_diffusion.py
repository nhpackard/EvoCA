"""Which diffusion ratios matter? (linear/Turing part; pulse-speed part lives in s10_pulse_props.py)
 (i) Linear (Turing) analysis of the active fixed point of the 3-species system: max Re lambda(q) over q for D=diag(Du,Dv,Dp).
 (ii) 1D discrete simulations: one-way pulse speed / width / survival over (du, dp) at dv fixed, and over dv."""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
from gslib import *

def jac3(u, v, a, rho):
    # d/d(u,v,p) of (-a u v^2 + rho p, a u v^2 - v, v - rho p)
    return np.array([[-a*v*v, -2*a*u*v, rho], [a*v*v, 2*a*u*v - 1, 0], [0, 1, -rho]])

def max_growth(a, rho, Dr_u, Dr_p, branch=1):
    """Dispersion of upper fixed point with Dv=1 (units absorbed in q): max over q^2 of Re lambda."""
    pts = fixed_points(a, rho)
    if len(pts) < 3: return np.nan
    u, v, p = pts[branch]; J = jac3(u, v, a, rho)
    best = -np.inf
    for q2 in np.logspace(-4, 2, 120):
        ev = np.linalg.eigvals(J - q2*np.diag([Dr_u, 1.0, Dr_p]))
        best = max(best, ev.real.max())
    return best

def pulse(args):
    a, rho, k2, du, dv, dp = args
    k1 = a*k2; r = rho*k2; L = 500; c = 100
    u = np.ones(L); v = np.zeros(L); p = np.zeros(L)
    v[c:c+6] = 0.4; u[c:c+6] -= 0.4; p[c-40:c] = 0.5; u[c-40:c] -= 0.5
    pos = []
    for t in range(1500):
        u, v, p = step_discrete(u, v, p, k1, k2, r, du, dv, dp)
        if t % 50 == 49: pos.append(np.argmax(v) if v.max() > 1e-3 else -1)
    alive = v.max() > 1e-3 and pos[-1] > pos[len(pos)//2] and pos[-1] > 0
    sp = np.nan; w = np.nan
    if alive:
        x = np.arange(L); sp = (pos[-1]-pos[len(pos)//2])/(50*(len(pos)-1-len(pos)//2))
        w = (v > 0.5*v.max()).sum()
    return alive, sp, w, v.max()

if __name__ == '__main__':
    # (i) Turing check
    Dr = np.logspace(-1, 2, 25)
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for k, (a, rho) in enumerate([(30, 3.0), (400, 0.2)]):
        pts = fixed_points(a, rho)
        if len(pts) == 1:  # no active fixed point
            ax[k].set_title(f'a={a}, rho={rho}: no active fixed point'); continue
        Z = np.array([[max_growth(a, rho, du_, dp_, branch=2) for du_ in Dr] for dp_ in Dr])  # branch 2 = lower-v branch (saddle) skip
        Z1 = np.array([[max_growth(a, rho, du_, dp_, branch=1) for du_ in Dr] for dp_ in Dr])
        q = ax[k].pcolormesh(Dr, Dr, Z1, shading='auto', cmap='RdBu_r', vmin=-abs(np.nanmax(Z1)), vmax=abs(np.nanmax(Z1))); plt.colorbar(q, ax=ax[k])
        ax[k].set(xscale='log', yscale='log', xlabel='Du/Dv', ylabel='Dp/Dv', title=f'max Re(lambda) over q, high-v fixed point, a={a}, rho={rho}\n(at q=0 the value is the homogeneous rate; red = unstable)')
    plt.tight_layout(); plt.savefig('figs/f5a_turing_dispersion.png', dpi=90); plt.close()
    # does diffusion change the stability of the upper fixed point vs the q=0 value? print summary
    for (a, rho) in [(20, 0.3), (30, 3.0), (200, 0.05), (400, 0.2)]:
        pts = fixed_points(a, rho)
        if len(pts) < 3: continue
        J = jac3(*pts[1][:2], a, rho)
        e0 = np.linalg.eigvals(J).real.max()
        worst = max(max_growth(a, rho, x, y) for x in Dr for y in Dr)
        print(f'a={a} rho={rho}: q=0 max Re={e0:.4g}; max over all (Du/Dv,Dp/Dv,q) = {worst:.4g}  -> Turing-type gain only if this exceeds q=0 value when q=0 is <0')
    # (ii) removed: pulse propagation vs diffusion weights is measured properly in s10_pulse_props.py
