"""Clean 1D pulse measurement (two-sided kick on a long ring, first-arrival times at +-150, +-250) and the (du, dv) propagation map."""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
from gslib import *

def pulse_props(args, kickw=20, L=1000, T=1500):
    a, rho, k2, du, dv, dp = args
    k1 = a*k2; r = rho*k2; c = L//2
    u = np.ones(L); v = np.zeros(L); p = np.zeros(L)
    v[c-kickw//2:c+kickw//2] = 0.4; u[c-kickw//2:c+kickw//2] = 0.6
    probes = [c+150, c+250]; arr = [None, None]; prof = None; vpk = 0
    for t in range(T):
        u, v, p = step_discrete(u, v, p, k1, k2, r, du, dv, dp)
        for i, x in enumerate(probes):
            if arr[i] is None and v[x-1:x+2].max() > 0.05: arr[i] = t
        if arr[1] is not None and prof is None:
            prof = v[probes[1]-40:probes[1]+10].copy()
    if arr[1] is None: return dict(prop=False)
    sp = 100/(arr[1]-arr[0])
    j = np.argmax(prof); half = prof[j]/2; lo = j
    while lo > 0 and prof[lo] > half: lo -= 1
    hi = j
    while hi < len(prof)-1 and prof[hi] > half: hi += 1
    # wake (refractory) length: distance behind peak until v < 0.01 * peak behind (profile truncated at 40)
    return dict(prop=True, speed=sp, fwhm=hi-lo, vpeak=float(prof[j]), second=float(v[c+200:].max() > 0.05))

if __name__ == '__main__':
    k2 = 0.1
    sets = {'osc (79.4,.128)': (79.4, .128), 'osc2 (141.5,.094)': (141.5, .094), 'excit (86.8,.092)': (86.8, .092), 'spiral (60,.05)': (60, .05), 'spiral2 (60,.03)': (60, .03), 'front (25.6,.346)': (25.6, .346)}
    for nm, (a, rho) in sets.items():
        print(nm, pulse_props((a, rho, k2, 0.6, 0.3, 0.3)))
    dus = [0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.9]; dvs = [0.02, 0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.9]
    with Pool(12) as P: res = P.map(pulse_props, [(60, 0.05, k2, du, dv, 0.3) for dv in dvs for du in dus])
    S = np.array([r['speed'] if r['prop'] else np.nan for r in res]).reshape(len(dvs), len(dus))
    print('speed map (rows dv, cols du):'); print(np.round(S, 3))
    fig, ax = plt.subplots(figsize=(7, 5.5))
    q = ax.imshow(S, origin='lower', cmap='viridis'); plt.colorbar(q, label='pulse speed (sites/step); blank = no propagation')
    ax.set_xticks(range(len(dus)), dus); ax.set_yticks(range(len(dvs)), dvs); ax.set(xlabel='du (box-blur weight)', ylabel='dv', title='1D pulse propagation vs diffusion weights\n(a=60, rho=0.05, k2=0.1, dp=0.3, kick width 20)')
    plt.tight_layout(); plt.savefig('figs/f10_pulse_du_dv.png', dpi=100)
    dps = [0.0, 0.05, 0.15, 0.3, 0.6, 0.9]
    with Pool(12) as P: res = P.map(pulse_props, [(60, 0.05, k2, 0.6, 0.3, dp) for dp in dps])
    print('dp sweep:', [(dp, round(r['speed'], 3) if r['prop'] else None, r.get('fwhm')) for dp, r in zip(dps, res)])
