"""Task 1b: limit cycles and excitability in the well-mixed ODE (scaled a, rho). Parallel scan."""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
from gslib import *

def classify_point(args):
    a, rho = args
    pts = fixed_points(a, rho)
    res = dict(a=a, rho=rho, peak=0.0, final='none', amp=0.0, period=np.nan, thr=np.nan)
    # excitability: smallest kick v0 (with u=1-v0,p=0) for which excursion peak > 0.2 (large)
    T = 40/min(rho, 1) + 200
    def run(y0, T=T, n=4000):
        return simulate(a, rho, y0, T, n)
    # threshold scan on kicks
    ks = np.array([0.005, 0.01, 0.02, 0.04, 0.08, 0.15, 0.3, 0.6])
    big = None
    for k in ks:
        s = run([1-k, k, 0.0], T=T, n=1500)
        if s.y[1].max() > 0.2:
            big = k; break
    res['thr'] = big if big is not None else np.nan
    # fate of a large kick: tail behaviour
    s = run([0.7, 0.3, 0.0], T=3*T, n=6000)
    v = s.y[1]; tail = v[-2000:]
    res['peak'] = v.max()
    if tail.max() < 1e-4: res['final'] = 'return'
    elif tail.max() - tail.min() < 1e-5*max(1, tail.max()): res['final'] = 'active'
    else:
        res['final'] = 'osc'; res['amp'] = tail.max()-tail.min()
        # period from upward crossings of mean
        m = tail.mean(); idx = np.where((tail[:-1] < m) & (tail[1:] >= m))[0]
        if len(idx) >= 2: res['period'] = np.mean(np.diff(s.t[-2000:][idx]))
    # oscillation also from perturbed active focus
    if len(pts) == 3 and res['final'] != 'osc':
        u, vv, p = pts[1]
        s = run([u*1.0, vv*1.02, p-0.02*vv], T=3*T, n=6000)
        t2 = s.y[1][-2000:]
        if t2.max()-t2.min() > 1e-4 and t2.min() > 1e-4:
            res['final'] = 'osc'; res['amp'] = t2.max()-t2.min()
            m = t2.mean(); idx = np.where((t2[:-1] < m) & (t2[1:] >= m))[0]
            if len(idx) >= 2: res['period'] = np.mean(np.diff(s.t[-2000:][idx]))
    return res

if __name__ == '__main__':
    Ag = np.logspace(0.5, 3.0, 36); Rg = np.logspace(-2.5, 1.0, 36)
    jobs = [(a, r) for r in Rg for a in Ag]
    with Pool(12) as P: out = P.map(classify_point, jobs, chunksize=8)
    sh = (len(Rg), len(Ag))
    peak = np.array([o['peak'] for o in out]).reshape(sh); thr = np.array([o['thr'] for o in out]).reshape(sh)
    code = {'none': 0, 'return': 1, 'active': 2, 'osc': 3}
    fin = np.array([code[o['final']] for o in out]).reshape(sh)
    per = np.array([o['period'] for o in out]).reshape(sh); amp = np.array([o['amp'] for o in out]).reshape(sh)
    np.savez('s1b_results.npz', Ag=Ag, Rg=Rg, peak=peak, thr=thr, fin=fin, per=per, amp=amp)
    fig, ax = plt.subplots(1, 4, figsize=(21, 4.8))
    cm = matplotlib.colors.ListedColormap(['#dddddd', '#ccbb44', '#4477aa', '#ee6677'])
    ax[0].pcolormesh(Ag, Rg, fin, cmap=cm, vmin=-.5, vmax=3.5, shading='auto')
    ax[0].set(title='fate of kick (u,v,p)=(0.7,0.3,0):\nyellow excursion+return (excitable), blue -> active steady, red sustained oscillation', xscale='log', yscale='log', xlabel='a', ylabel='rho')
    ax[0].plot(4*(1+1/Rg), Rg, 'k--', lw=1); ax[0].axvline(4, color='k', ls=':'); ax[0].set_xlim(Ag[0], Ag[-1]); ax[0].set_ylim(Rg[0], Rg[-1])
    q = ax[1].pcolormesh(Ag, Rg, peak, shading='auto', cmap='viridis'); plt.colorbar(q, ax=ax[1]); ax[1].set(title='peak v excursion after kick', xscale='log', yscale='log', xlabel='a')
    q = ax[2].pcolormesh(Ag, Rg, np.log10(thr), shading='auto', cmap='plasma'); plt.colorbar(q, ax=ax[2]); ax[2].set(title='log10 threshold kick v0 for excursion peak>0.2\n(blank = none up to 0.6)', xscale='log', yscale='log', xlabel='a')
    q = ax[3].pcolormesh(Ag, Rg, per, shading='auto', cmap='cividis'); plt.colorbar(q, ax=ax[3]); ax[3].set(title='limit-cycle period (1/k2 units)', xscale='log', yscale='log', xlabel='a')
    plt.tight_layout(); plt.savefig('figs/f1c_excitability_cycles.png', dpi=100); plt.close()
    ii = np.argwhere(fin == 3)
    print('n osc points', len(ii))
    for i, j in ii[::max(1, len(ii)//15)]: print(f'a={Ag[j]:.4g} rho={Rg[i]:.4g} amp={amp[i,j]:.3g} period={per[i,j]:.4g} (sat-node a={4*(1+1/Rg[i]):.4g})')
