"""
Reciprocal causal-control experiment (research-directions §7.3).

Three arms at a fixed productive corner (mirrors the egene_discussion
V3 validation corner so results are directly comparable):

  JOINT      : LUT and egene both evolve            (baseline)
  LUT_only   : egene frozen (mu_egene=mu_egenome=0), LUT evolves
  EGENE_only : LUT frozen (mu_lut=0), egene evolves

run_sim already seeds lut='gol' + egenome='uniform' every run, so the
two freeze conditions need NO code change — just zeroed mutation rates
(this is the "zero new code" claim in §7.3, verified).

Each (arm, seed) runs in its OWN process (Pool maxtasksperchild=1) so
the per-process C xorshift RNG is fresh — this sidesteps the
pool-reuse determinism bug documented in scan_analysis.md and means
these particular results ARE per-config reproducible.
"""
import csv
import os
import sys
from multiprocessing import Pool

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(_HERE, '..', '..', 'python')))
from evoca_explore import run_sim   # noqa: E402

N = 128
N_STEPS = 3000
SAMPLE_EVERY = 100
SEEDS = (0, 1, 2)

# Fixed productive corner (egene_discussion V3): cognition is costly
# (tax_per_egene=8e-4) so "spec stays high" is a real selection
# signal, not a free default.
BASE = dict(food_inc=0.013, m_scale=1.2, gdiff=0.06,
            tax=0.035, tax_per_egene=8e-4)

ARMS = {
    'JOINT':      dict(BASE, mu_lut=0.01, mu_egene=0.003, mu_egenome=0.005),
    'LUT_only':   dict(BASE, mu_lut=0.01, mu_egene=0.0,   mu_egenome=0.0),
    'EGENE_only': dict(BASE, mu_lut=0.0,  mu_egene=0.003, mu_egenome=0.005),
}

REPORT = ['cog_specificity_mean', 'cog_load_mean', 'mean_intake_mean',
          'n_distinct_genomes_mean', 'excess_pc_slope',
          'excess_activity_slope', 'final_pop', 'min_pop', 'extinct']


def _one(job):
    arm, seed = job
    out = run_sim(ARMS[arm], n_steps=N_STEPS, sample_every=SAMPLE_EVERY,
                  N=N, seed=seed, shadow=True, init='halfplane')
    out['arm'] = arm
    out['seed'] = seed
    return out


def main():
    jobs = [(arm, s) for arm in ARMS for s in SEEDS]
    with Pool(processes=4, maxtasksperchild=1) as pool:
        rows = pool.map(_one, jobs)

    cols = ['arm', 'seed'] + REPORT
    csv_path = os.path.join(_HERE, 'results.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # Compact mean-over-seeds summary to stdout.
    print(f"\n=== causal control: N={N} steps={N_STEPS} "
          f"seeds={SEEDS} ===")
    hdr = f"{'arm':<11}" + "".join(f"{c:>22}" for c in REPORT[:5])
    print(hdr)
    for arm in ARMS:
        ar = [r for r in rows if r['arm'] == arm]
        line = f"{arm:<11}"
        for c in REPORT[:5]:
            vals = [r.get(c) for r in ar if isinstance(r.get(c), (int, float))]
            m = sum(vals) / len(vals) if vals else float('nan')
            line += f"{m:>22.4f}"
        print(line)
    print(f"\nfull rows -> {csv_path}")


if __name__ == '__main__':
    main()
