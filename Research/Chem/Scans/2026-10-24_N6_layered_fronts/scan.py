"""CHEM-N6, part 2: layered fronts. Recipe 2026-10-23_chem_fw_spiralpoint
with wdiff 0.01 unless set, 6000 steps; seeds 0-2.

  oneway   N's interactive point (one-way exchange, waste_tox 0.4,
           waste_sink 0.019), never measured by the CHEM-N1 statistics:
           wdiff {0.01, 0.08} x N {512, 1024}, and N=256 at wdiff 0.01 (as watched)
  twoway   N=1024, the strongest CHEM-N1 points by ring contrast and by
           defect survival: (sink 0.019, wdiff 0.03, tox 0.6),
           (0.015, 0.08, 1.0), (0.019, 0.01, 0.4)

Saved per run (frames/, not committed): 4 x 4 block means of alive and of
environmental W every 20 steps over 5000-6000 (float16), the middle row of
alive and W every step over 5400-6000, and the colormode-6 frame at 6000.
The statistics are computed offline by analyze.py, so they can change
without a re-run."""

import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np


def _repo_root(p):
    while not os.path.exists(os.path.join(p, 'C', 'evoca.h')):
        if os.path.dirname(p) == p:
            raise FileNotFoundError('EvoCA repo root not found')
        p = os.path.dirname(p)
    return p


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = _repo_root(HERE)
sys.path.insert(0, os.path.join(ROOT, 'python'))

import evoca_waves as ew                                     # noqa: E402
from evoca_explore import build_sim, save_scan_config        # noqa: E402

RECIPE = os.path.join(ROOT, 'Research', 'Chem', 'Runs', '2026-10-23_chem_fw_spiralpoint.evoca')
N_WORKERS = 28
N_STEPS = 6000
FRAMES_FROM, FRAME_EVERY = 5000, 20
KYMO_FROM = 5400
BLOCK = 4
SEEDS = [0, 1, 2]


def rgb6(sim, N):
    c = sim.colorize(np.zeros(N * N, np.int32), 6).view(np.uint32).reshape(N, N)
    return np.dstack([(c >> 16) & 255, (c >> 8) & 255, c & 255]).astype(np.uint8)


def tag(group, N, seed, over):
    return (f"{group}_N{N}_ws{over['waste_sink']}_wd{over['wdiff']}"
            f"_tox{over['waste_tox']}_x{over['waste_exchange']}_s{seed}")


def run_one(task):
    group, N, seed, over = task
    params = dict(json.load(open(RECIPE))['metaparams_init'])
    params.update(over)
    t0 = time.perf_counter()
    sim = build_sim(params, N=N, seed=seed, init='halfplane')
    A, W, ka, kw, alive = [], [], [], [], []
    for t in range(1, N_STEPS + 1):
        sim.step()
        if t >= FRAMES_FROM and t % FRAME_EVERY == 0:
            a = sim.get_alive().reshape(N, N).astype(float)
            A.append(ew.block_means(a, BLOCK).astype(np.float16))
            W.append(ew.block_means(sim.get_W().reshape(N, N), BLOCK).astype(np.float16))
        if t > KYMO_FROM:
            ka.append(sim.get_alive().reshape(N, N)[N // 2].astype(np.uint8))
            kw.append(sim.get_W().reshape(N, N)[N // 2].astype(np.float16))
        if t % 100 == 0:
            alive.append(float(sim.get_alive().mean()))
            if alive[-1] == 0:
                break
    last = rgb6(sim, N)
    sim.free()
    os.makedirs(os.path.join(HERE, 'frames'), exist_ok=True)
    np.savez_compressed(os.path.join(HERE, 'frames', tag(group, N, seed, over) + '.npz'),
                        alive=np.array(A), W=np.array(W), kymo_alive=np.array(ka),
                        kymo_W=np.array(kw), last=last, alive_series=np.array(alive),
                        block=BLOCK, frame_every=FRAME_EVERY, params=json.dumps(params))
    return tag(group, N, seed, over), time.perf_counter() - t0, alive[-1] if alive else 0.0


def tasks():
    out = []
    one = dict(waste_exchange=0, waste_sink=0.019, waste_tox=0.4)
    for wd in (0.01, 0.08):
        for N in (512, 1024):
            out += [('oneway', N, s, dict(one, wdiff=wd)) for s in SEEDS]
    out += [('oneway', 256, s, dict(one, wdiff=0.01)) for s in SEEDS]
    for ws, wd, tox in ((0.019, 0.03, 0.6), (0.015, 0.08, 1.0), (0.019, 0.01, 0.4)):
        out += [('twoway', 1024, s, dict(waste_exchange=1, waste_sink=ws, wdiff=wd, waste_tox=tox))
                for s in SEEDS]
    out.sort(key=lambda t: -t[1])
    return out


def main():
    cfg = save_scan_config(HERE, N=1024, n_steps=N_STEPS, sample_every=FRAME_EVERY, seed=0, shadow=False)
    c = json.load(open(cfg))
    c.update({'N': 'per run: 256, 512, 1024', 'seed': '0, 1, 2', 'recipe': os.path.relpath(RECIPE, HERE)})
    json.dump(c, open(cfg, 'w'), indent=2)
    ts = tasks()
    print(f"Running {len(ts)} runs across {N_WORKERS} workers...", flush=True)
    t0 = time.perf_counter()
    with Pool(N_WORKERS) as pool:
        for name, sec, a in pool.imap_unordered(run_one, ts):
            print(f"  {name}: {sec:5.0f}s alive {a:.3f}  [{time.perf_counter() - t0:5.0f}s]", flush=True)
    print('done')


if __name__ == '__main__':
    main()
