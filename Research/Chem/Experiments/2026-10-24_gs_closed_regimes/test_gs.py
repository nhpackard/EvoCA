"""Tests: mass conservation, nonnegativity, numba==numpy, scaled-ODE agreement of the discrete reaction at small rates,
effective diffusion coefficient of the box blur (D = d/3), reachability of diffusion params through public step."""
import numpy as np
from gslib import *
from gs2d import make_stepper
rng = np.random.default_rng(0)
N = 48
u = rng.random((N, N)); v = rng.random((N, N))*.3; p = rng.random((N, N))*.2
m0 = (u+v+p).sum()
u2, v2, p2 = u.copy(), v.copy(), p.copy()
st = make_stepper((N, N))
for _ in range(5):
    u, v, p = step_discrete(u, v, p, 0.7, 0.1, 0.05, 0.6, 0.3, 0.2)
st(u2, v2, p2, 0.7, 0.1, 0.05, 0.6, 0.3, 0.2, nsteps=5)
assert abs((u+v+p).sum()-m0) < 1e-9*m0 and min(u.min(), v.min(), p.min()) >= 0
assert np.allclose(u, u2) and np.allclose(v, v2) and np.allclose(p, p2), 'numba != numpy'
# diffusion: Gaussian spread variance after n steps = n*2*D*... 1D: var growth per step = d*2/3 -> D=d/3
x = np.zeros(401); x[200] = 1
for d in (0.3, 0.6):
    y = x.copy()
    for _ in range(100): y = blur(y, d)
    var = (y*(np.arange(401)-200)**2).sum()
    assert abs(var/100/2 - d/3) < 1e-6, (var/100/2, d/3)
# reachability: changing dv alters result
a = rng.random((N, N)); b = a.copy()
r1 = step_discrete(a, b*.2, a*0, .5, .1, .05, .6, .3, .3)[1]; r2 = step_discrete(a, b*.2, a*0, .5, .1, .05, .6, .5, .3)[1]
assert not np.allclose(r1, r2)
# discrete well-mixed vs ODE at small rates (a=100, rho=0.1, k2=0.01): compare v(t) peak time
k2 = 0.01; a_, rho = 100., 0.1
wm = well_mixed_discrete(a_*k2, k2, rho*k2, (0.9, 0.1, 0), 3000)
s = simulate(a_, rho, [0.9, 0.1, 0], 3000*k2, 3000)
assert abs(wm[:, 1].argmax()*k2 - s.t[s.y[1].argmax()]) < 0.3, (wm[:, 1].argmax()*k2, s.t[s.y[1].argmax()])
print('all tests passed')
