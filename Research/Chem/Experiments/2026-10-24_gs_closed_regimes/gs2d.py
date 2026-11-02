"""Fast 2D discrete-time step (numba): reaction (exponential integrator per channel) + 3x3 box blur, periodic.
Same arithmetic as gslib.step_discrete (checked in test_gs.py)."""
import numpy as np
from numba import njit, prange

@njit(parallel=True, cache=True)
def react_nb(u, v, p, k1, k2, r):
    N, M = u.shape
    e2 = 1 - np.exp(-k2); e3 = 1 - np.exp(-r)
    for i in prange(N):
        for j in range(M):
            uu, vv, pp = u[i, j], v[i, j], p[i, j]
            c = uu*(1 - np.exp(-k1*vv*vv)); uu -= c; vv += c
            c = vv*e2; vv -= c; pp += c
            c = pp*e3; pp -= c; uu += c
            u[i, j], v[i, j], p[i, j] = uu, vv, pp

@njit(parallel=True, cache=True)
def blur_nb(x, out, d):
    N, M = x.shape
    for i in prange(N):
        i0 = (i-1) % N; i1 = (i+1) % N
        for j in range(M):
            j0 = (j-1) % M; j1 = (j+1) % M
            s = (x[i0, j0]+x[i0, j]+x[i0, j1]+x[i, j0]+x[i, j]+x[i, j1]+x[i1, j0]+x[i1, j]+x[i1, j1])/9.0
            out[i, j] = (1-d)*x[i, j] + d*s

def make_stepper(shape):
    tmp = np.empty(shape)
    def step(u, v, p, k1, k2, r, du, dv, dp, nsteps=1):
        for _ in range(nsteps):
            react_nb(u, v, p, k1, k2, r)
            for x, d in ((u, du), (v, dv), (p, dp)):
                blur_nb(x, tmp, d); x[:] = tmp
    return step
