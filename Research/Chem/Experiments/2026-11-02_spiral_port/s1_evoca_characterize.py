"""Characterize N's EvoCA spiral recipe (2026-11-02_spiral_wave_ish) as the
target for a ChemCA port: the effective rates a port must match, and whether
the spirals need evolution or the waste rule.

    python3 s1_evoca_characterize.py <variant> [N] [steps]

variants: base (the recipe), nomut (all mutation off), exch1 (two-way waste
exchange, the ChemCA transport rule). Writes out/s1_<variant>_N<N>.npz.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', '..', 'python'))
import evoca_nb  # noqa: E402

RECIPE = '2026-11-02_spiral_wave_ish.evoca'
VARIANTS = {
    'base': {},
    'nomut': dict(mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0),
    'exch1': dict(waste_exchange=1),
}


def main(variant, N=256, steps=5000, kymo=2000, every=10):
    sim, _ = evoca_nb.load_recipe(RECIPE, N=N, seed=0, **VARIANTS[variant])
    ts = {k: [] for k in ('t', 'alive', 'tax_pc', 'intake_pc', 'eat_frac', 'F', 'W',
                          'F_alive', 'w_alive', 'tox', 'births', 'deaths')}
    kym = {k: [] for k in ('alive', 'F', 'W')}
    frames = {}
    tox = births = deaths = 0
    for t in range(1, steps + 1):
        sim.step()
        tox += sim.get_tox_deaths_last()
        births += sim.get_births_last()
        deaths += sim.get_deaths_last()
        need_fields = t % every == 0 or t > steps - kymo or t % 500 == 0
        if not need_fields:
            continue
        A = sim.get_alive().astype(bool)
        F, W = sim.get_F(), sim.get_W()
        if t > steps - kymo:
            kym['alive'].append(A[N // 2].astype(np.uint8))
            kym['F'].append(F[N // 2].astype(np.float32))
            kym['W'].append(W[N // 2].astype(np.float32))
        if t % 500 == 0:
            frames[t] = dict(rgb=evoca_nb.rgb(sim, 1), F=F.astype(np.float32),
                             W=W.astype(np.float32), alive=A)
        if t % every == 0:
            d = sim.get_dissipation()
            n = max(int(A.sum()), 1)
            m = d['mouthful_step']
            ts['t'].append(t)
            ts['alive'].append(A.mean())
            ts['tax_pc'].append(d['tax_step'] / n)
            ts['intake_pc'].append(m / n)
            # F at alive sites before eating = after + eaten (eating is the
            # last F change in the step; births after it are a small error)
            ts['eat_frac'].append(m / max(float(F[A].sum()) + m, 1e-12))
            ts['F'].append(F.mean())
            ts['W'].append(W.mean())
            ts['F_alive'].append(F[A].mean() if A.any() else np.nan)
            ts['w_alive'].append(sim.get_w()[A].mean() if A.any() else np.nan)
            ts['tox'].append(tox); ts['births'].append(births); ts['deaths'].append(deaths)
            tox = births = deaths = 0
        if t % 1000 == 0:
            print(variant, N, t, f"alive {A.mean():.3f}", flush=True)
    sim.free()
    os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
    np.savez_compressed(
        os.path.join(HERE, 'out', f's1_{variant}_N{N}.npz'),
        **{f'ts_{k}': np.array(v) for k, v in ts.items()},
        **{f'kym_{k}': np.array(v) for k, v in kym.items()},
        frame_t=np.array(sorted(frames)),
        **{f'fr_{k}': np.array([frames[t][k] for t in sorted(frames)]) for k in ('rgb', 'F', 'W', 'alive')})


if __name__ == '__main__':
    v = sys.argv[1]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 256
    steps = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
    main(v, N, steps)
