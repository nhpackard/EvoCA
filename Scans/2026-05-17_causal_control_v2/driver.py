"""Causal control v2 — corrected metric, real run length.

The 2026-05-16 v1 run exposed a bug: the Channon shadow modelled
LUT-only while activity hashed the full genome, so excess activity
was undefined for egene-driven evolution (EGENE_only excess ≈ −26).
S2d (`eg-dyn-shadows`, merged) added FIXED-SPACE neutral baselines:
`eg_excess_pc_slope` (729 ternary egene keys) and `dyn_excess_pc_slope`
(500 LUT-entry buckets), per-component-normalised. S1 (`evoca_set_seed`)
makes each run reproducible.

Scientific question this can now actually answer, per arm:
  EGENE_only (LUT frozen): is eg_excess_pc_slope > 0 ?  → egene
      evolution genuinely beats neutral drift through egene space.
  LUT_only  (egene frozen): is dyn_excess_pc_slope > 0 ? → the
      alive population's *effective* rule transitions are under
      selection (not silent hash churn).
  JOINT: both positive ⇒ both subsystems under selection together.

Each (arm, seed) runs in its own process (maxtasksperchild=1) so the
C PRNG is fresh; with S1's set_seed these results are reproducible.
"""
import csv
import os
import sys
from multiprocessing import Pool

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(_HERE, '..', '..', 'python')))
from evoca_explore import run_sim   # noqa: E402

N = 256
N_STEPS = 15000
SAMPLE_EVERY = 250
SEEDS = (0, 1, 2, 3)

BASE = dict(food_inc=0.013, m_scale=1.2, gdiff=0.06,
            tax=0.035, tax_per_egene=8e-4)

ARMS = {
    'JOINT':      dict(BASE, mu_lut=0.01, mu_egene=0.003, mu_egenome=0.005),
    'LUT_only':   dict(BASE, mu_lut=0.01, mu_egene=0.0,   mu_egenome=0.0),
    'EGENE_only': dict(BASE, mu_lut=0.0,  mu_egene=0.003, mu_egenome=0.005),
}

REPORT = ['cog_specificity_mean', 'cog_load_mean', 'mean_intake_mean',
          'n_distinct_genomes_mean',
          'excess_pc_slope', 'eg_excess_pc_slope', 'dyn_excess_pc_slope',
          'final_pop', 'min_pop', 'extinct']


def _one(job):
    arm, seed = job
    out = run_sim(ARMS[arm], n_steps=N_STEPS, sample_every=SAMPLE_EVERY,
                  N=N, seed=seed, shadow=True, init='halfplane')
    out['arm'] = arm
    out['seed'] = seed
    return out


def main():
    jobs = [(a, s) for a in ARMS for s in SEEDS]
    with Pool(processes=6, maxtasksperchild=1) as pool:
        rows = pool.map(_one, jobs)
    cols = ['arm', 'seed'] + REPORT
    with open(os.path.join(_HERE, 'results.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\n=== causal control v2: N={N} steps={N_STEPS} seeds={SEEDS} ===")
    hdr = f"{'arm':<11}" + "".join(f"{c:>20}" for c in
                                   ['cog_spec', 'intake', 'n_distinct',
                                    'excess_pc', 'eg_excess_pc',
                                    'dyn_excess_pc'])
    print(hdr)
    keys = ['cog_specificity_mean', 'mean_intake_mean',
            'n_distinct_genomes_mean', 'excess_pc_slope',
            'eg_excess_pc_slope', 'dyn_excess_pc_slope']
    for arm in ARMS:
        ar = [r for r in rows if r['arm'] == arm]
        line = f"{arm:<11}"
        for k in keys:
            vals = [r[k] for r in ar if isinstance(r.get(k), (int, float))]
            m = sum(vals) / len(vals) if vals else float('nan')
            line += f"{m:>20.4f}"
        print(line)
    print(f"\nfull rows -> {os.path.join(_HERE, 'results.csv')}")


if __name__ == '__main__':
    main()
