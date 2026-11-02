"""Scan the ChemCA port over the knobs that cannot be ported by number:
u_F (food uptake: EvoCA eats by pattern match, capped by headroom), tox_W
and u_W (ChemCA makes waste when food is converted, i.e. at intake, where
EvoCA makes it as tax is paid), and sink_W (EvoCA's sink is a constant
amount, ChemCA's is proportional), and sink_F (EvoCA discards only the
food above its cap of 1; a proportional sink wastes food at every level).

    python3 s2_scan.py [n_workers] [N] [steps]

Writes out/s2_results.csv and out/s2/<i>.npz (final frames, kymographs).
"""
import itertools
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', '..', 'python'))
import portlib as P  # noqa: E402
import evoca_waves as ew  # noqa: E402

GRID = dict(u_F=[0.2, 0.35, 0.5, 0.75, 1.0],
            tox_W=[0.792, 1.2, 1.6, 2.4, 0.0],      # 0 = no toxicity
            u_W=[0.09, 0.2, 0.4],
            sink_W=[0.1, 0.05],
            sink_F=[0.012, 0.004, 0.0])             # 0.012: F saturates at 1; 0: unbounded


def configs():
    keys = list(GRID)
    return [dict(zip(keys, v)) for v in itertools.product(*GRID.values())]


def one(args):
    i, cfg, N, steps = args
    p, st = P.port_params(m_e=0.0501, **cfg)      # nomut EvoCA's effective tax
    d = P.run(p, st, N=N, steps=steps, kymo=1500, every=10,
              frames_at=range(steps - 1500, steps + 1, 250))
    t = d['ts_t']
    late = t > t.max() - 2000 if len(t) else t
    ext = bool(len(t) == 0 or d['ts_alive'][-1] == 0)
    row = dict(i=i, **cfg, steps_run=int(t.max()) if len(t) else 0, extinct=ext)
    if not ext and late.any():
        per, pk = ew.site_period(d['kym_W'], b=4) if len(d['kym_W']) > 100 else (np.nan, 0)
        perA, pkA = ew.site_period(d['kym_alive'], b=4) if len(d['kym_alive']) > 100 else (np.nan, 0)
        tot = d['ts_deaths'][late].sum()
        row.update(alive=float(d['ts_alive'][late].mean()), alive_sd=float(d['ts_alive'][late].std()),
                   F=float(d['ts_F'][late].mean()), W=float(d['ts_W'][late].mean()),
                   w_alive=float(np.nanmean(d['ts_w_alive'][late])),
                   intake_pc=float(np.nanmean(d['ts_intake_pc'][late])),
                   eat_frac=float(np.nanmean(d['ts_eat_frac'][late])),
                   tox_frac=float(d['ts_tox'][late].sum() / max(tot, 1)),
                   period_W=per, acf_W=pk, period_alive=perA, acf_alive=pkA)
    os.makedirs(os.path.join(HERE, 'out', 's2'), exist_ok=True)
    np.savez_compressed(os.path.join(HERE, 'out', 's2', f'{i:03d}.npz'),
                        **{k: d[k] for k in ('ts_t', 'ts_alive', 'ts_F', 'ts_W', 'kym_alive', 'kym_W',
                                             'frame_t', 'fr_rgb', 'fr_F', 'fr_W')})
    print(i, row.get('alive'), row.get('period_W'), flush=True)
    return row


if __name__ == '__main__':
    import pandas as pd
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 256
    steps = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
    jobs = [(i, c, N, steps) for i, c in enumerate(configs())]
    with Pool(nw) as pool:
        rows = pool.map(one, jobs, chunksize=1)
    pd.DataFrame(rows).to_csv(os.path.join(HERE, 'out', 's2_results.csv'), index=False)
    print('done', len(rows))
