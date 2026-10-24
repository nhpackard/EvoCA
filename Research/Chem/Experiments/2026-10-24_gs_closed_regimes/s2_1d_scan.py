"""Task 2a: 1D discrete-time scan over (a, rho) at fixed per-step k2 and diffusion weights.
Initial: trivial state (u=m,v=0,p=0), one-sided kick (v-bump with refractory p on the left) on a ring."""
import numpy as np, sys, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
from gslib import *

def run1d(a, rho, k2=0.1, m=1.0, du=0.6, dv=0.3, dp=0.3, L=600, T=1500, rec=False, kick=0.4, onesided=True):
    k1 = a*k2/m**2; r = rho*k2
    u = np.full(L, m); v = np.zeros(L); p = np.zeros(L)
    c = L//2
    v[c:c+6] = kick*m; u[c:c+6] -= kick*m
    if onesided:   # refractory (waste-rich) zone on left of kick
        p[c-40:c] = 0.5*m; u[c-40:c] -= 0.5*m
    hist = []
    for t in range(T):
        u, v, p = step_discrete(u, v, p, k1, k2, r, du, dv, dp)
        if rec and t % 2 == 0: hist.append(v.copy())
    return (u, v, p, np.array(hist)) if rec else (u, v, p)

def metrics(args):
    a, rho, kw = args
    L = kw.get('L', 600); T = kw.get('T', 1500)
    k2 = kw.get('k2', 0.1)
    u, v, p, h = run1d(a, rho, rec=True, **kw)
    # activity: fraction of sites with v>0.01 over time (every 2 steps)
    act = (h > 0.02).mean(1)
    frac_end = act[-100:].mean()
    # fate (classified on late-time v(x,t) record, last 200 samples = 400 steps)
    late = h[-200:]
    tstd = late.std(0).mean()          # temporal std per site, averaged
    sstd = late[-1].std()              # spatial std of last frame
    if late.max() < 1e-3: fate = 'extinct'
    elif tstd > 0.01: fate = 'dynamic'                 # sustained waves/pulses/oscillation
    elif sstd < 0.01*max(late.mean(), 1e-9) + 1e-4 and late.mean() > 0.01: fate = 'uniform steady'
    else: fate = 'frozen'                              # static non-uniform pattern or still-invading front
    # front speed: extent of active region at t=100..300 (before wrapping)
    ext = (h > 0.02).sum(1)
    sp = np.polyfit(np.arange(50, 150)*2, ext[50:150], 1)[0]/2 if ext[50:150].max() > 0 else 0   # sites/step (both sides)
    return dict(a=a, rho=rho, fate=fate, frac_end=frac_end, speed=sp, vmean=h[-1].mean(), vstd=h[-50:].std(0).mean())

if __name__ == '__main__':
    Ag = np.logspace(0.7, 3.0, 14); Rg = np.logspace(-2.0, 0.5, 14)
    kw = dict(k2=0.1, du=0.6, dv=0.3, dp=0.3, L=600, T=2500)
    jobs = [(a, r, kw) for r in Rg for a in Ag]
    with Pool(12) as P: out = P.map(metrics, jobs)
    sh = (len(Rg), len(Ag))
    code = {'extinct': 0, 'frozen': 1, 'uniform steady': 2, 'dynamic': 3}
    fate = np.array([code[o['fate']] for o in out]).reshape(sh); fe = np.array([o['frac_end'] for o in out]).reshape(sh)
    sp = np.array([o['speed'] for o in out]).reshape(sh)
    np.savez('s2_scan.npz', Ag=Ag, Rg=Rg, fate=fate, fe=fe, sp=sp)
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    cm = matplotlib.colors.ListedColormap(['#dddddd', '#999944', '#4477aa', '#ee6677'])
    ax[0].pcolormesh(Ag, Rg, fate, cmap=cm, vmin=-.5, vmax=3.5, shading='auto'); ax[0].set(title='1D discrete (k2=0.1,du=.6,dv=.3,dp=.3): fate at T=1500\ngrey extinct, olive frozen, blue uniform active, red dynamic', xscale='log', yscale='log', xlabel='a', ylabel='rho')
    q = ax[1].pcolormesh(Ag, Rg, fe, shading='auto'); plt.colorbar(q, ax=ax[1]); ax[1].set(title='active fraction at end', xscale='log', yscale='log', xlabel='a')
    q = ax[2].pcolormesh(Ag, Rg, sp, shading='auto', cmap='magma'); plt.colorbar(q, ax=ax[2]); ax[2].set(title='front expansion speed (sites/step, 2 fronts)', xscale='log', yscale='log', xlabel='a')
    plt.tight_layout(); plt.savefig('figs/f2a_1d_scan.png', dpi=100)
    from collections import Counter; print(Counter(o['fate'] for o in out))
    for o in out:
        if o['fate'] in ('dynamic','frozen'): print(o['fate'], round(o['a'],2), round(o['rho'],3), round(o['frac_end'],2))
