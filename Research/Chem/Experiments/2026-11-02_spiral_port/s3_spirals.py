"""Spiral-core persistence, EvoCA recipe vs ChemCA port candidates, N=512.

A spiral core is a phase defect of theta = atan2(W - W0, F - F0) on 8 x 8
blocks (evoca_waves), with (F0, W0) each block's time average. A core
persists for many periods; noise defects do not. Readout: of the cores at
time t, the fraction with a same-sign core within 2 blocks at t + lag,
against the same with the later frame shifted by half the lattice (chance).

    python3 s3_spirals.py <case> <seed>     # one run -> out/s3/<case>_s<seed>.npz
    python3 s3_spirals.py report            # table from all out/s3/*.npz
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', '..', 'python'))
import evoca_waves as ew  # noqa: E402

N, STEPS, WINDOW, EVERY, B = 512, 5000, 2500, 10, 8
# ChemCA candidates: s2 scan indices (out/s2_ranked.csv)
CHEM = {
    'chem109': dict(u_F=0.35, tox_W=1.2, u_W=0.09, sink_W=0.1, sink_F=0.004),
    'chem307': dict(u_F=0.75, tox_W=1.6, u_W=0.09, sink_W=0.1, sink_F=0.004),
    'chem418': dict(u_F=1.0, tox_W=2.4, u_W=0.09, sink_W=0.05, sink_F=0.004),
    'chem13': dict(u_F=0.2, tox_W=0.792, u_W=0.4, sink_W=0.1, sink_F=0.004),
}
EVO = {'evoca': {}, 'evoca_nomut': dict(mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0)}


def make(case, seed):
    """(step(), fields() -> (F, W, alive)) for a case."""
    if case in EVO:
        import evoca_nb
        sim, _ = evoca_nb.load_recipe('2026-11-02_spiral_wave_ish.evoca', N=N, seed=seed, **EVO[case])
        return sim, lambda: (sim.get_F(), sim.get_W(), sim.get_alive())
    import portlib as P
    from chemca_py import ChemCA
    p, st = P.port_params(m_e=0.0501, **CHEM[case])
    sim = ChemCA()
    sim.init(N, **p)
    sim.state(**st, seed=seed)
    return sim, lambda: (sim.get_E('F'), sim.get_E('W'), sim.get_alive())


def run(case, seed):
    sim, fields = make(case, seed)
    Fb, Wb, Ab = [], [], []
    for t in range(1, STEPS + 1):
        sim.step()
        if t > STEPS - WINDOW and t % EVERY == 0:
            F, W, A = fields()
            Fb.append(ew.block_means(F, B)); Wb.append(ew.block_means(W, B))
            Ab.append(ew.block_means(A.astype(float), B))
    sim.free()
    os.makedirs(os.path.join(HERE, 'out', 's3'), exist_ok=True)
    np.savez_compressed(os.path.join(HERE, 'out', 's3', f'{case}_s{seed}.npz'),
                        Fb=np.array(Fb, np.float32), Wb=np.array(Wb, np.float32), Ab=np.array(Ab, np.float32))


def persistence(Fb, Wb, lags_steps=(150, 300, 600), amp_min=0.05):
    """{lag: (matched fraction, chance fraction, mean cores per frame)}."""
    center = (Fb.mean(axis=0), Wb.mean(axis=0))
    cs = [ew.cores(*ew.phase_field(F, W, 1, center), amp_min=amp_min) for F, W in zip(Fb, Wb)]
    n = Fb.shape[1]
    out = {}
    for lag in lags_steps:
        k = lag // EVERY
        m = c = tot = 0
        for (pa, sa), (pb, sb) in zip(cs[:-k], cs[k:]):
            tot += len(pa)
            m += ew.matched(pa, sa, pb, sb, n)
            c += ew.matched(pa, sa, (pb + n // 2) % n, sb, n)
        out[lag] = (m / max(tot, 1), c / max(tot, 1), tot / max(len(cs) - k, 1))
    return out


def report():
    import glob
    import pandas as pd
    rows = []
    for p in sorted(glob.glob(os.path.join(HERE, 'out', 's3', '*.npz'))):
        z = np.load(p)
        case_seed = os.path.basename(p)[:-4]
        per, pk = ew.period_of(z['Ab'].mean(axis=(1, 2)), EVERY)
        r = dict(run=case_seed, alive=float(z['Ab'].mean()), F=float(z['Fb'].mean()), W=float(z['Wb'].mean()),
                 period=per, acf=pk)
        for lag, (m, c, nc) in persistence(z['Fb'], z['Wb']).items():
            r.setdefault('cores', nc)              # per frame (shortest lag: most frames)
            r[f'match{lag}'] = m
            r[f'chance{lag}'] = c
        rows.append(r)
    pd.set_option('display.width', 250)
    print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == '__main__':
    if sys.argv[1] == 'report':
        report()
    else:
        run(sys.argv[1], int(sys.argv[2]))
