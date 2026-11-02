"""Task 1: well-mixed analysis of closed Gray-Scott in scaled variables (a=k1 m^2/k2, rho=r/k2)."""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from gslib import *

# --- 1. regime map over (a, rho): stability of the nontrivial fixed points
A = np.logspace(0, 3.5, 300); R = np.logspace(-3, 1.5, 300)
regime = np.zeros((len(R), len(A)))  # 0 none; 1 stable node/focus (upper); 2 unstable focus (osc); 3 bistable-ish
nontrivial = np.zeros_like(regime); hopf_a = []
for i, rho in enumerate(R):
    for j, a in enumerate(A):
        pts = fixed_points(a, rho)
        if len(pts) == 1:
            regime[i, j] = 0; continue
        # upper-v branch is pts[1] (s=+1: larger v, smaller u)
        u, v, p = pts[1]
        ev = np.linalg.eigvals(jac2(u, v, a, rho))
        if np.all(ev.real < 0):
            regime[i, j] = 1 if np.any(abs(ev.imag) > 0) else 2   # stable focus / stable node
        else:
            regime[i, j] = 3 if np.any(abs(ev.imag) > 0) else 4     # 3 unstable focus, 4 unstable node/saddle
print('regime counts', {k: int((regime == k).sum()) for k in np.unique(regime)})

# trivial state (1,0,0): Jacobian [[-1,0],[1,-rho]] -> always stable (v decays: -1, p: -rho). So trivial state always stable.
# --- 2. Numerical: for each (a, rho), integrate from a perturbation of the nontrivial branch and from kicks, classify long-term
def long_term(a, rho, y0, T=400):
    s = simulate(a, rho, y0, T, n=2000)
    v = s.y[1]; tail = v[len(v)//2:]
    return tail.max() - tail.min(), tail.mean(), s
Ag = np.logspace(0.3, 3.3, 40); Rg = np.logspace(-2.5, 1, 40)
amp = np.zeros((len(Rg), len(Ag))); vend = np.zeros_like(amp)
for i, rho in enumerate(Rg):
    for j, a in enumerate(Ag):
        # initial: small v seed in food-rich medium (v=0.05) -> tests whether the active state is reached
        y0 = [0.9, 0.1, 0.0]
        T = 60/min(rho, 1)
        s = simulate(a, rho, y0, T, n=3000)
        tail = s.y[1][-1000:]
        amp[i, j] = tail.max() - tail.min(); vend[i, j] = tail.mean()
fig, ax = plt.subplots(1, 3, figsize=(17, 5))
cmap = matplotlib.colors.ListedColormap(['#dddddd', '#4477aa', '#ee6677', '#ccbb44', '#228833'])
im = ax[0].pcolormesh(A, R, regime, cmap=matplotlib.colors.ListedColormap(['#dddddd', '#4477aa', '#88ccee', '#ee6677', '#ccbb44']), vmin=-0.5, vmax=4.5, shading='auto')
ax[0].set(xscale='log', yscale='log', xlabel='a = k1 m^2 / k2', ylabel='rho = r / k2',
          title='Upper fixed point: grey none, dark blue stable focus,\nlight blue stable node, red unstable focus (limit cycle)')
aa = 4*(1 + 1/R); ax[0].plot(aa, R, 'k--', lw=1, label='saddle-node a=4(1+1/rho)'); ax[0].legend(loc='lower left')
ax[0].set_xlim(A[0], A[-1]); ax[0].set_ylim(R[0], R[-1])
p1 = ax[1].pcolormesh(Ag, Rg, amp, shading='auto', cmap='magma'); plt.colorbar(p1, ax=ax[1])
ax[1].set(xscale='log', yscale='log', xlabel='a', ylabel='rho', title='ODE from (u,v,p)=(0.9,0.1,0): late-time v peak-to-peak\n(>0 sustained oscillation)')
p2 = ax[2].pcolormesh(Ag, Rg, vend, shading='auto', cmap='viridis'); plt.colorbar(p2, ax=ax[2])
ax[2].set(xscale='log', yscale='log', xlabel='a', ylabel='rho', title='late-time mean v (0 = extinct)')
plt.tight_layout(); plt.savefig('figs/f1_wellmixed_regimes.png', dpi=110); plt.close()

# --- Hopf locus: trace(J)=0 at upper branch; scan rho, find a
from scipy.optimize import brentq
def tr_det(a, rho):
    u, v, p = fixed_points(a, rho)[1]; J = jac2(u, v, a, rho); return np.trace(J), np.linalg.det(J)
hopf = []
for rho in np.logspace(-3, 1.5, 60):
    aa_ = np.logspace(np.log10(4*(1+1/rho))+1e-9, 6, 400)
    tr = np.array([tr_det(a, rho)[0] for a in aa_])
    idx = np.where(np.sign(tr[:-1]) != np.sign(tr[1:]))[0]
    hopf.append((rho, [brentq(lambda a: tr_det(a, rho)[0], aa_[k], aa_[k+1]) for k in idx]))
np.save('hopf_locus.npy', np.array([(r, *(h if h else [np.nan])) for r, h in hopf if len(h) <= 2 and len(h)>=1], dtype=object), allow_pickle=True)
print('Hopf locus (rho: a values where trace=0):')
for rho, h in hopf[::6]: print(f'  rho={rho:.4g}  a_hopf={["%.4g"%x for x in h]}  saddle-node a={4*(1+1/rho):.4g}')

# --- example trajectories + excitability
fig, ax = plt.subplots(2, 3, figsize=(16, 8))
# a limit cycle example: pick rho=0.05, a=?
cases = [('oscillatory', 100, 0.03), ('excitable-ish', 20, 0.3), ('stable focus', 30, 3)]
for k, (nm, a, rho) in enumerate(cases):
    pts = classify(a, rho)
    print(nm, 'a', a, 'rho', rho, [(tuple(np.round(q, 4)), np.round(e, 4)) for q, e in pts])
    s = simulate(a, rho, [0.9, 0.1, 0.0], 400/ min(1,rho*5), 6000)
    ax[0, k].plot(s.t, s.y[0], label='u'); ax[0, k].plot(s.t, s.y[1], label='v'); ax[0, k].plot(s.t, s.y[2], label='p')
    ax[0, k].set(title=f'{nm}: a={a}, rho={rho} (time in 1/k2)', xlabel='t'); ax[0, k].legend()
    ax[1, k].plot(s.y[1], s.y[2]); ax[1, k].set(xlabel='v', ylabel='p', title='phase plane (v,p)')
plt.tight_layout(); plt.savefig('figs/f1b_wellmixed_examples.png', dpi=100); plt.close()
