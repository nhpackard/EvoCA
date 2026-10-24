"""Characterize candidate parameter sets: well-mixed discrete dynamics, minimal sustainable stimulation period,
2D spiral (period at a probe, radial spectrum wavelength, speed = wl/period, survival vs lattice size and k2)."""
import numpy as np, matplotlib, sys, json; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from multiprocessing import Pool
from gslib import *
from gs2d import make_stepper
from s4_2d import seed_broken

def min_period(k1, k2, r, du, dv, dp, L=500):
    """periodic stimulation at left end: smallest period P (steps) for which pulses 1..N all reach x=300"""
    out = {}
    for P in (40, 60, 80, 100, 120, 150, 200, 300, 400):
        u = np.ones(L); v = np.zeros(L); p = np.zeros(L); cnt = 0; prev = False
        for t in range(3000):
            if t % P == 0 and t < 3000-P*1:
                v[5:9] = np.maximum(v[5:9], 0.4); u[5:9] = np.minimum(u[5:9], 1 - v[5:9])  # inject only where still excitable? injection adds v taken from u
            u, v, p = step_discrete(u, v, p, k1, k2, r, du, dv, dp)
            now = v[300] > 0.05
            if now and not prev: cnt += 1
            prev = now
        out[P] = cnt
    return out

def spiral(args):
    a, rho, k2, du, dv, dp, N, T = args
    k1 = a*k2; r = rho*k2
    u, v, p = seed_broken(N); st = make_stepper((N, N))
    probes = [(N//4, N//4), (N//4, 3*N//4), (3*N//4, N//2)]
    rec = []
    st(u, v, p, k1, k2, r, du, dv, dp, nsteps=T//2)
    for t in range(T//2):
        st(u, v, p, k1, k2, r, du, dv, dp, nsteps=1); rec.append([v[i, j] for i, j in probes])
    rec = np.array(rec); periods = []
    for i in range(3):
        idx, _ = find_peaks(rec[:, i], height=0.05)
        if len(idx) >= 3: periods.append(np.median(np.diff(idx)))
    f = v - v.mean(); P = np.abs(np.fft.fft2(f))**2
    ky = np.fft.fftfreq(N)[:, None]; kx = np.fft.fftfreq(N)[None, :]; kk = np.sqrt(kx**2+ky**2)
    bins = np.arange(0, 0.3, 1.0/N); h = np.histogram(kk, bins, weights=P)[0]/np.maximum(np.histogram(kk, bins)[0], 1)
    cen = (bins[1:]+bins[:-1])/2; jj = np.argmax(h[2:])+2; w_ = h[jj-1:jj+2]; kpk = (cen[jj-1:jj+2]*w_).sum()/w_.sum()
    alive = (v > 0.02).mean()
    return dict(a=a, rho=rho, k2=k2, N=N, alive=alive, period=float(np.median(periods)) if periods else None,
                wl=float(1/kpk) if alive > 0 else None, vmax=float(v.max()))

if __name__ == '__main__':
    sets = {'osc': (79.4, 0.128), 'osc2': (141.5, 0.094), 'excit': (86.8, 0.092), 'spiral': (60, 0.05), 'spiral2': (60, 0.03), 'front': (25.6, 0.346)}
    k2 = 0.1; d = (0.6, 0.3, 0.3)
    print('--- well-mixed discrete: long-run v after kick (last 1000 steps min/max)')
    for nm, (a, rho) in sets.items():
        wm = well_mixed_discrete(a*k2, k2, rho*k2, (0.6, 0.4, 0.0), 8000)[:, 1]
        t = wm[-2000:]; idx, _ = find_peaks(t, height=0.01)
        print(nm, a, rho, 'v range', t.min().round(4), t.max().round(4), 'period', np.median(np.diff(idx)) if len(idx) > 2 else None)
    # (1D pulse speed/width are measured in s10_pulse_props.py)
    print('--- min sustainable stimulation period (pulses reaching x=300 out of injections), P: count')
    for nm in ('excit', 'spiral', 'osc'):
        a, rho = sets[nm]; print(nm, min_period(a*k2, k2, rho*k2, *d), flush=True)
    print('--- 2D spirals')
    jobs = [(60, 0.05, 0.1, *d, N, 8000) for N in (256, 512)] + [(60, 0.03, 0.1, *d, 256, 10000)] \
         + [(60, 0.05, k, *d, 256, 8000) for k in (0.03, 0.2)]
    for dd in ((0.6, 0.3, 0.0), (0.6, 0.3, 0.9), (0.3, 0.3, 0.3), (0.9, 0.3, 0.3), (0.9, 0.15, 0.3), (0.45, 0.3, 0.3)):
        jobs.append((60, 0.05, 0.1, *dd, 256, 6000))
    for j in jobs: print('d=', j[3:6], {k: v for k, v in spiral(j).items() if k != 'a'}, flush=True)
