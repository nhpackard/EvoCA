"""Task 2b: pulse-train emitters (target waves) in 1D discrete model: period, speed, wavelength over (a, rho).
Symmetric center kick, probes at c+60, c+120; peaks of v(t)."""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from multiprocessing import Pool
from gslib import *

def run_probe(a, rho, k2=0.1, m=1.0, du=0.6, dv=0.3, dp=0.3, L=600, T=5000, kick=0.4, probes=(60, 120)):
    k1 = a*k2/m**2; r = rho*k2
    u = np.full(L, m); v = np.zeros(L); p = np.zeros(L); c = L//2
    v[c:c+6] = kick*m; u[c:c+6] -= kick*m
    rec = np.empty((T, len(probes)))
    for t in range(T):
        u, v, p = step_discrete(u, v, p, k1, k2, r, du, dv, dp)
        for i, o in enumerate(probes): rec[t, i] = v[c+o]
    return rec, v

def analyze(a, rho, **kw):
    probes = kw.pop('probes', (60, 120))
    rec, vend = run_probe(a, rho, probes=probes, **kw)
    out = dict(a=a, rho=rho, n=0, period=np.nan, speed=np.nan, wl=np.nan, fate='none')
    pk = []
    for i in range(len(probes)):
        x = rec[:, i]
        idx, _ = find_peaks(x, prominence=0.02 * max(x.max(), 1e-9), height=0.01)
        pk.append(idx)
    out['n'] = len(pk[0])
    if len(pk[0]) == 0: out['fate'] = 'extinct'; return out
    if rec[-1500:, 0].max() < 1e-3 and len(pk[0]) == 1: out['fate'] = 'single pulse'
    if len(pk[0]) >= 3:
        # period from late peaks
        out['period'] = np.median(np.diff(pk[0][len(pk[0])//2:])) if len(pk[0]) >= 4 else np.diff(pk[0]).mean()
        out['fate'] = 'pulse train'
    if len(pk[0]) >= 1 and len(pk[1]) >= 1:
        dt = pk[1][0] - pk[0][0]; out['speed'] = (probes[1]-probes[0])/dt
        if not np.isnan(out['period']): out['wl'] = out['speed']*out['period']
    if len(pk[0]) < 3 and rec[-1000:, 0].mean() > 0.01: out['fate'] = 'steady active'
    return out

def job(args):
    a, rho, kw = args
    return analyze(a, rho, **dict(kw))

if __name__ == '__main__':
    Ag = np.logspace(np.log10(25), np.log10(600), 12); Rg = np.logspace(np.log10(0.02), np.log10(0.6), 12)
    kw = dict(k2=0.1, du=0.6, dv=0.3, dp=0.3, L=600, T=5000)
    jobs = [(a, r, kw) for r in Rg for a in Ag]
    with Pool(12) as P: out = P.map(job, jobs)
    sh = (len(Rg), len(Ag)); g = lambda k: np.array([o[k] for o in out], float).reshape(sh)
    per, sp, wl = g('period'), g('speed'), g('wl')
    code = {'extinct': 0, 'single pulse': 1, 'steady active': 2, 'pulse train': 3, 'none': 0}
    fate = np.array([code[o['fate']] for o in out]).reshape(sh)
    np.savez('s3_results.npz', Ag=Ag, Rg=Rg, per=per, sp=sp, wl=wl, fate=fate)
    fig, ax = plt.subplots(1, 4, figsize=(21, 4.8))
    cm = matplotlib.colors.ListedColormap(['#dddddd', '#ccbb44', '#4477aa', '#ee6677'])
    ax[0].pcolormesh(Ag, Rg, fate, cmap=cm, vmin=-.5, vmax=3.5, shading='auto'); ax[0].set(title='1D discrete, center kick: grey extinct, yellow single pulse,\nblue steady active, red pulse-train source', xscale='log', yscale='log', xlabel='a', ylabel='rho')
    for axx, z, t in zip(ax[1:], (per, sp, wl), ('emission period (steps)', 'pulse speed (sites/step)', 'wavelength = speed x period (sites)')):
        q = axx.pcolormesh(Ag, Rg, z, shading='auto'); plt.colorbar(q, ax=axx); axx.set(title=t, xscale='log', yscale='log', xlabel='a')
    plt.tight_layout(); plt.savefig('figs/f3_sources_1d.png', dpi=100)
    for o in out:
        if o['fate'] == 'pulse train': print(f"a={o['a']:.1f} rho={o['rho']:.3f} period={o['period']:.1f} speed={o['speed']:.3f} wl={o['wl']:.1f}")
