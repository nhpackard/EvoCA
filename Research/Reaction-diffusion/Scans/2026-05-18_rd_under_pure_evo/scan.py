"""Priority #4 — do the 2026-04-27 evo+spatial (RD) winners keep their
spatial structure under a long pure-evolution run?

User's central intuition: robust evolution needs life-death turnover,
so pure-evo selection should NOT collapse to static structure; if it
does, suspect the metrics. #1/#2 hinted correlation_length persists;
this is the direct test — take the actual RD-generating winner
metaparams and run them long, then ask whether spatial structure
*decays over the run* (final ≪ mean ⇒ collapsing toward static/death)
or *persists* (final ≈ mean ⇒ RD robust without spatial optimisation).

Bases: the 30 `Runs/2026-04-27_top_*` balanced-composite winners
(excludes top_evo/top_spatial/auto). Each run records spatial
final-vs-mean plus the corrected pure-evo metrics, multi-seed.
"""
import csv
import glob
import json
import os
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'python'))
from evoca_explore import run_sim, save_scan_config   # noqa: E402

SEEDS = (0, 1, 2)
N_WORKERS = 30
N_STEPS = 12000          # long enough to see structure decay if it will
N_GRID = 256
SAMPLE_EVERY = 100


def load_bases():
    out, seen = [], set()
    for p in sorted(glob.glob(os.path.join(ROOT, 'Runs',
                                           '2026-04-27_top_*.evoca'))):
        b = os.path.basename(p)
        if 'top_evo' in b or 'top_spatial' in b or 'auto' in b:
            continue
        d = json.load(open(p))
        mp = d.get('metaparams_init') or d.get('metaparams') or {}
        key = tuple(sorted(mp.items()))
        if key in seen or not mp:
            continue
        seen.add(key)
        out.append((b, mp))
    return out


def run_one(arg):
    cfg_idx, (label, mp, seed) = arg
    try:
        r = run_sim(mp, n_steps=N_STEPS, sample_every=SAMPLE_EVERY,
                    N=N_GRID, seed=seed, shadow=True, init='halfplane')
        r['error'] = ''
    except Exception as e:
        r = {'error': repr(e), **{f'param_{k}': v for k, v in mp.items()}}
    r['config_idx'] = cfg_idx
    r['base'] = label
    r['seed'] = seed
    # decay diagnostics: <1 ⇒ structure shrank over the run
    for m in ('correlation_length', 'largest_patch'):
        mn, fn = r.get(f'{m}_mean'), r.get(f'{m}_final')
        r[f'{m}_final_over_mean'] = (fn / mn) if (mn not in (None, 0)
                                                  and fn is not None) else ''
    return r


def main():
    save_scan_config(HERE, N=N_GRID, n_steps=N_STEPS,
                     sample_every=SAMPLE_EVERY, seed=0, shadow=True)
    bases = load_bases()
    jobs = [(lbl, mp, s) for (lbl, mp) in bases for s in SEEDS]
    print(f"[#4 RD-under-pure-evo] {len(bases)} winner recipes x "
          f"{len(SEEDS)} seeds = {len(jobs)} runs / {N_WORKERS} workers "
          f"/ {N_STEPS} ticks", flush=True)
    t0 = time.perf_counter()
    args = list(enumerate(jobs))
    results = []
    with Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(run_one, args):
            results.append(r)
            d = len(results)
            if d % 20 == 0 or d == len(args):
                dt = time.perf_counter() - t0
                rate = d / dt if dt else 0
                eta = (len(args) - d) / rate if rate else 0
                print(f"  {d}/{len(args)}  {dt:5.0f}s  ~{eta:.0f}s left",
                      flush=True)
    dt = time.perf_counter() - t0
    print(f"Done in {dt:.0f}s.")
    out_path = os.path.join(HERE, 'results.csv')
    keys = sorted({k for r in results for k in r})
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(results)
    print(f"Wrote {len(results)} rows -> {out_path}")


if __name__ == '__main__':
    main()
