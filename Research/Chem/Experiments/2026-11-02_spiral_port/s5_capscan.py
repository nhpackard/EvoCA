"""s2 again with ChemCA's food cap: food is now EvoCA's rule exactly
(inc_F 0.012, cap_F 1, sink_F 0), leaving u_F, tox_W, u_W, sink_W calibrated.

    python3 s5_capscan.py [n_workers] [N] [steps]

Writes out/s5_results.csv and out/s5/<i>.npz.
"""
import itertools
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s2_scan as S2  # noqa: E402

GRID = dict(u_F=[0.2, 0.35, 0.5, 0.75, 1.0],
            tox_W=[0.792, 1.2, 1.6, 2.4, 0.0],      # 0 = no toxicity
            u_W=[0.09, 0.2, 0.4],
            sink_W=[0.1, 0.05])


def one(args):
    i, cfg, N, steps = args
    row = S2.one((i, dict(cfg, cap_F=1.0), N, steps), sub='s5')
    return row


if __name__ == '__main__':
    import pandas as pd
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 256
    steps = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
    keys = list(GRID)
    jobs = [(i, dict(zip(keys, v)), N, steps) for i, v in enumerate(itertools.product(*GRID.values()))]
    with Pool(nw) as pool:
        rows = pool.map(one, jobs, chunksize=1)
    pd.DataFrame(rows).to_csv(os.path.join(HERE, 'out', 's5_results.csv'), index=False)
    print('done', len(rows))
