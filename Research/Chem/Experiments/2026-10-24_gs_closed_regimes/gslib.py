"""Shared library: closed-loop Gray-Scott (u,v,p) well-mixed + spatial integrators.
Reactions: u+2v->3v (k1 u v^2), v->p (k2 v), p->u (r p).  Mass u+v+p = mbar conserved.
Scaled form: time in 1/k2, concentration in mbar:  a = k1*mbar^2/k2, rho = r/k2.
  du = -a u v^2 + rho p ; dv = a u v^2 - v ; dp = v - rho p ;  u+v+p = 1.
"""
import numpy as np
from scipy.integrate import solve_ivp

def rhs3(t, y, a, rho):
    u, v, p = y
    c = a*u*v*v
    return [-c + rho*p, c - v, v - rho*p]

def fixed_points(a, rho):
    """nontrivial fixed points (scaled, mbar=1): (1+1/rho) v^2 - v + 1/a = 0 ; plus trivial (1,0,0)."""
    A = 1 + 1/rho
    disc = 1 - 4*A/a
    pts = [(1.0, 0.0, 0.0)]
    if disc >= 0:
        for s in (+1, -1):
            v = (1 + s*np.sqrt(disc))/(2*A)
            pts.append((1/(a*v), v, v/rho))
    return pts

def jac2(u, v, a, rho):
    """Jacobian of reduced (v,p) system, u = 1-v-p."""
    # dv = a(1-v-p)v^2 - v ; dp = v - rho p
    return np.array([[a*(2*u*v - v*v) - 1, -a*v*v], [1.0, -rho]])

def classify(a, rho):
    out = []
    for (u, v, p) in fixed_points(a, rho):
        J = jac2(u, v, a, rho)
        ev = np.linalg.eigvals(J)
        out.append(((u, v, p), ev))
    return out

def simulate(a, rho, y0, T, n=4000):
    return solve_ivp(rhs3, (0, T), y0, args=(a, rho), method='LSODA', rtol=1e-9, atol=1e-12,
                     t_eval=np.linspace(0, T, n))

# ---------------- spatial: continuous-time explicit Euler (dx=1) ----------------
def lap(x):
    if x.ndim == 1:
        return np.roll(x, 1) + np.roll(x, -1) - 2*x
    return np.roll(x, 1, 0) + np.roll(x, -1, 0) + np.roll(x, 1, 1) + np.roll(x, -1, 1) - 4*x

def step_euler(u, v, p, a, rho, Du, Dv, Dp, dt):
    c = a*u*v*v
    nu = u + dt*(-c + rho*p + Du*lap(u))
    nv = v + dt*(c - v + Dv*lap(v))
    np_ = p + dt*(v - rho*p + Dp*lap(p))
    return nu, nv, np_

# ---------------- spatial: discrete-time (reaction by exponential integrator, then 3x3 box blur) ----
def blur(x, d):
    """x <- (1-d) x + d * mean over 3x3 (2D) or 3-point (1D) neighborhood, periodic. D_eff = d/3 (2D) , d/3 (1D)."""
    if x.ndim == 1:
        m = (np.roll(x, 1) + x + np.roll(x, -1))/3
    else:
        r = np.roll(x, 1, 0) + x + np.roll(x, -1, 0)
        m = (np.roll(r, 1, 1) + r + np.roll(r, -1, 1))/9
    return (1-d)*x + d*m

def react_discrete(u, v, p, k1, k2, r):
    """Per-channel exponential integrator (sequential, mass-conserving, non-negative):
       1. autocatalysis: converted = u*(1-exp(-k1 v^2))  (u->v, rate k1 v^2 per unit u)
       2. decay: v -> p, converted = v*(1-exp(-k2))
       3. regeneration: p -> u, converted = p*(1-exp(-r))
    Channels use state after the previous channel (Lie splitting); order 1->2->3."""
    c = u*(1 - np.exp(-k1*v*v)); u = u - c; v = v + c
    c = v*(1 - np.exp(-k2)); v = v - c; p = p + c
    c = p*(1 - np.exp(-r)); p = p - c; u = u + c
    return u, v, p

def step_discrete(u, v, p, k1, k2, r, du, dv, dp):
    u, v, p = react_discrete(u, v, p, k1, k2, r)
    return blur(u, du), blur(v, dv), blur(p, dp)

def well_mixed_discrete(k1, k2, r, y0, nsteps):
    u, v, p = [np.array(float(x)) for x in y0]
    out = np.empty((nsteps, 3))
    for i in range(nsteps):
        u, v, p = react_discrete(u, v, p, k1, k2, r)
        out[i] = (u, v, p)
    return out
