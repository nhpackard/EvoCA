"""Sanity check: open Gray-Scott (Pearson 1993) with explicit Euler, dt=1, dx=1, D_u=0.16, D_v=0.08
(equivalent to Pearson's D_u=2e-5, D_v=1e-5, dx=0.01, dt=1 up to D*dt/dx^2 = 0.2/0.1 -- same ratio 2).
Known patterns: (F,k)=(0.035,0.065) spots; (0.025,0.060) / (0.030,0.062) worms/ stripes-labyrinth; (0.022,0.051) spiral waves/ solitons-ish;
(0.037,0.060) mitosis/ spots splitting."""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numba import njit

@njit(cache=True)
def gs_steps(U, V, F, k, Du, Dv, n):
    N = U.shape[0]
    U2 = U.copy(); V2 = V.copy()
    for _ in range(n):
        for i in range(N):
            im = (i-1) % N; ip = (i+1) % N
            for j in range(N):
                jm = (j-1) % N; jp = (j+1) % N
                u = U[i, j]; v = V[i, j]
                lu = U[im, j]+U[ip, j]+U[i, jm]+U[i, jp]-4*u
                lv = V[im, j]+V[ip, j]+V[i, jm]+V[i, jp]-4*v
                uvv = u*v*v
                U2[i, j] = u + Du*lu - uvv + F*(1-u)
                V2[i, j] = v + Dv*lv + uvv - (F+k)*v
        U, U2 = U2, U; V, V2 = V2, V
    return U, V

def init(N, seed=0):
    rng = np.random.default_rng(seed)
    U = np.ones((N, N)); V = np.zeros((N, N)); c = N//2; w = N//10
    U[c-w:c+w, c-w:c+w] = 0.5; V[c-w:c+w, c-w:c+w] = 0.25
    U += 0.02*rng.standard_normal((N, N)); V += 0.02*rng.standard_normal((N, N))
    return U, V

if __name__ == '__main__':
    pats = [('spots (mitosis)', 0.037, 0.060), ('spots', 0.035, 0.065), ('worms/stripes', 0.030, 0.057), ('waves/solitons', 0.014, 0.054)]
    fig, axs = plt.subplots(1, 4, figsize=(20, 5.4))
    for ax, (nm, F, k) in zip(axs, pats):
        U, V = init(200)
        U, V = gs_steps(U, V, F, k, 0.16, 0.08, 20000)
        ax.imshow(V, cmap='viridis'); ax.set_title(f'open GS: {nm}\nF={F}, k={k}, Du=0.16, Dv=0.08, t=20000'); ax.axis('off')
        print(nm, 'V range', V.min(), V.max(), 'mean', V.mean())
    plt.tight_layout(); plt.savefig('figs/f0_open_gs_check.png', dpi=80)
