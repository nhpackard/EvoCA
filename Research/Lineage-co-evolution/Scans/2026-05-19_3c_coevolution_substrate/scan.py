"""Campaign #3c — egene evolution on a frozen *evolved* substrate.

Spec: Docs/plan-3c-coevolution-substrate.md (user-approved 2026-05-19).

Single continuous sim per run. Phase-1 evolves LUT-only; at the switch
the rule is either frozen (mu_lut->0) or kept evolving, and egene
mutation is enabled for phase-2. Decisive contrast = A vs C from a
byte-identical shared phase-1 substrate (S1 deterministic replay,
asserted via get_lut_hashes()). B / B' = GoL-frozen calibration
anchors to #3/3b.

Arms
  A  evolved-frozen   : phase-1 mu_lut=p1mu (LUT evolves) -> mu_lut=0
  C  evolved-co-evo   : phase-1 mu_lut=p1mu               -> mu_lut stays
  B  GoL-frozen        : GoL, mu_lut=0, phase-1=PHASE1 (settled GoL)
  B' GoL settled-short : GoL, mu_lut=0, phase-1=GOL_PRERUN (still active)
All arms: phase-2 = PHASE2 ticks with egene mutation ON.

Readout: eg_excess_pc / dyn_excess_pc slope over phase-2 (whole), and
EARLY vs LATE phase-2 windows (the user's long-timescale
"is evolution still proceeding robustly at 8000 ticks" diagnostic);
viability; substrate-evolvedness at the switch (get_lut_hashes());
C's phase-2 rule drift. Slopes are fit over phase-2 samples only, so
the S2d fixed-space shadow is effectively scoped to phase-2 without a
C-side reset (no C change — pure Python on main).

Designed for torque (32 cores); launched detached (nohup; the
recovered #4/R1 lesson — SSH is the only valid torque health check).
"""
import csv
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'python'))

from evoca_py import EvoCA, make_gol_lut, LUT_BYTES        # noqa: E402
from evoca_explore import save_scan_config                  # noqa: E402

N_GRID       = 256
PHASE1       = 3000
PHASE2       = 8000
GOL_PRERUN   = 100            # B' short settling
SAMPLE_EVERY = 100
SEEDS        = tuple(range(10))
P1_MU        = (0.01, 0.03, 0.06)     # phase-1 mu_lut: substrate quality
REGIME       = dict(m_scale=3.0, food_inc=0.015, tax=0.035, gdiff=0.06)
P2_MU_EGENE   = 0.01
P2_MU_EGENOME = 0.01
WIN          = 2000           # early/late window length (phase-2 ticks)
N_WORKERS    = 30


def _gol_hash():
    """FNV-1a mirror of C lut_hash_fn over the GoL LUT bytes."""
    gol = make_gol_lut()
    w = np.frombuffer(bytes(bytearray(np.asarray(gol, np.uint8)[:LUT_BYTES])),
                       dtype='<u4')
    h = 0x811c9dc5
    for x in w:
        h = ((h ^ int(x)) * 0x01000193) & 0xFFFFFFFF
    return h


GOL_HASH = _gol_hash()


def _new_sim(seed, mu_lut):
    """Fresh sim at the #3c regime, GoL-seeded, egene frozen (phase-1)."""
    sim = EvoCA()
    sim.init(N_GRID, mu_lut=mu_lut, mu_egene=0.0, mu_egenome=0.0,
             tax_per_egene=0.0, **REGIME)
    sim.set_seed(seed)
    sim.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    rng = np.random.default_rng(seed)
    sim.set_v(rng.integers(0, 2, (N_GRID, N_GRID), dtype=np.uint8))
    sim.set_f_all(0.5)
    sim.set_F_all(1.0)
    sim.set_alive_halfplane(0)
    sim.neutral_enable()
    return sim


def _snap(sim):
    """One phase-2 sample: the corrected metrics + pop + diversity."""
    sim._lib.evoca_eg_activity_update()
    sim._lib.evoca_activity_update()
    G = sim.get_activity(max_n=200000)
    gp = G['pop_count'] > 0
    alive = np.ctypeslib.as_array(sim._lib.evoca_get_alive(),
                                  shape=(N_GRID * N_GRID,))
    npop = int((alive == 1).sum())
    return dict(
        eg_excess_pc=float(sim.eg_excess_pc()),
        dyn_excess_pc=float(sim.dyn_excess_pc()),
        n_distinct_genomes=int(gp.sum()),
        pop=npop,
    )


def _substrate_diag(sim):
    """Evolvedness of the phase-1 substrate at the switch, via the
    just-merged get_lut_hashes() (#8 enabler)."""
    lh = sim.get_lut_hashes().ravel()
    alive = np.ctypeslib.as_array(sim._lib.evoca_get_alive(),
                                  shape=(N_GRID * N_GRID,)) == 1
    if not alive.any():
        return dict(sub_frac_evolved=0.0, sub_n_lut=0, sub_alive=0,
                    sub_digest=0)
    lha = lh[alive]
    return dict(
        sub_frac_evolved=float((lha != GOL_HASH).mean()),
        sub_n_lut=int(np.unique(lha).size),
        sub_alive=int(alive.sum()),
        sub_digest=int(np.bitwise_xor.reduce(
            (lha.astype(np.uint64) * np.arange(1, lha.size + 1,
                                               dtype=np.uint64)))),
    )


def _run_phase(sim, n_ticks, t0, sample=False):
    """Step n_ticks; if sample, collect _snap()s tagged with abs tick."""
    rows = []
    for i in range(n_ticks):
        sim.step()
        t = t0 + i + 1
        if sample and (t % SAMPLE_EVERY == 0):
            s = _snap(sim)
            s['t'] = t
            rows.append(s)
    return rows


def _phase2_summary(prefix, p2rows, switch_t):
    """Slope over whole phase-2 + EARLY vs LATE windows (the
    long-timescale ongoing-evolution diagnostic)."""
    out = {}
    if not p2rows:
        return out
    t = np.array([r['t'] for r in p2rows], dtype=float)
    rel = t - switch_t                       # 0 .. PHASE2
    for key in ('eg_excess_pc', 'dyn_excess_pc', 'n_distinct_genomes'):
        v = np.array([r[key] for r in p2rows], dtype=float)
        out[f'{prefix}{key}_switch'] = float(v[0])
        out[f'{prefix}{key}_final'] = float(v[-1])

        def _slope(mask):
            if mask.sum() >= 2:
                return float(np.polyfit(t[mask], v[mask], 1)[0])
            return 0.0
        out[f'{prefix}{key}_slope'] = _slope(np.ones_like(rel, bool))
        out[f'{prefix}{key}_slope_early'] = _slope(rel <= WIN)
        out[f'{prefix}{key}_slope_late'] = _slope(rel >= PHASE2 - WIN)
        early = v[rel <= WIN]
        late = v[rel >= PHASE2 - WIN]
        out[f'{prefix}{key}_late_minus_early_mean'] = (
            float(late.mean() - early.mean())
            if early.size and late.size else 0.0)
    pop = np.array([r['pop'] for r in p2rows], dtype=float)
    out[f'{prefix}final_pop'] = int(pop[-1])
    out[f'{prefix}min_pop'] = int(pop.min())
    out[f'{prefix}extinct'] = bool(pop.min() == 0)
    return out


def _one_arm(seed, mu_lut_p1, freeze, phase1_ticks, label):
    """Build sim, run phase-1, switch, run phase-2; return (row, diag,
    switch lut-hash array for the determinism assert)."""
    sim = _new_sim(seed, mu_lut_p1)
    _run_phase(sim, phase1_ticks, 0, sample=False)
    diag = _substrate_diag(sim)
    sw_lh = sim.get_lut_hashes().copy()
    if freeze:
        sim.update_mu_lut(0.0)
    sim.update_mu_egene(P2_MU_EGENE)
    sim.update_mu_egenome(P2_MU_EGENOME)
    p2 = _run_phase(sim, PHASE2, phase1_ticks, sample=True)
    sim.free()
    row = dict(arm=label, seed=seed, param_mu_lut_p1=mu_lut_p1,
               phase1_ticks=phase1_ticks, **REGIME)
    row.update(diag)
    row.update(_phase2_summary('', p2, phase1_ticks))
    return row, diag, sw_lh


def run_one(arg):
    cfg_idx, job = arg
    try:
        kind = job['kind']
        seed = job['seed']
        if kind == 'AC':
            p1 = job['p1mu']
            rowA, dA, hA = _one_arm(seed, p1, True,  PHASE1, 'A')
            rowC, dC, hC = _one_arm(seed, p1, False, PHASE1, 'C')
            det = bool(np.array_equal(hA, hC))
            for r in (rowA, rowC):
                r['determinism_ok'] = det
                r['config_idx'] = cfg_idx
                r['error'] = ''
            if not det:
                # shared-fork premise (D-4) violated — flag loudly
                rowA['error'] = rowC['error'] = 'AC_phase1_not_identical'
            return [rowA, rowC]
        elif kind == 'B':
            row, _, _ = _one_arm(seed, 0.0, True, PHASE1, 'B')
        elif kind == 'Bprime':
            row, _, _ = _one_arm(seed, 0.0, True, GOL_PRERUN, "B'")
        else:
            raise ValueError(kind)
        row['determinism_ok'] = ''
        row['config_idx'] = cfg_idx
        row['error'] = ''
        return [row]
    except Exception as e:                                  # noqa: BLE001
        return [dict(config_idx=cfg_idx, error=repr(e),
                     arm=job.get('kind'), seed=job.get('seed'))]


def main():
    save_scan_config(HERE, N=N_GRID, n_steps=PHASE1 + PHASE2,
                     sample_every=SAMPLE_EVERY, seed=0, shadow=True)
    jobs = []
    for s in SEEDS:
        for p1 in P1_MU:
            jobs.append({'kind': 'AC', 'seed': s, 'p1mu': p1})
        jobs.append({'kind': 'B', 'seed': s})
        jobs.append({'kind': 'Bprime', 'seed': s})
    args = list(enumerate(jobs))
    n_rows = len(SEEDS) * (len(P1_MU) * 2 + 2)
    print(f"[#3c coevolution-substrate] {len(jobs)} jobs -> {n_rows} rows "
          f"/ {N_WORKERS} workers / phase1={PHASE1} phase2={PHASE2} "
          f"@ N={N_GRID}", flush=True)
    t0 = time.perf_counter()
    results = []
    with Pool(N_WORKERS) as pool:
        for rows in pool.imap_unordered(run_one, args):
            results.extend(rows)
            d = len(results)
            if d % 20 == 0 or d >= n_rows:
                dt = time.perf_counter() - t0
                print(f"  {d}/{n_rows}  {dt:5.0f}s", flush=True)
    dt = time.perf_counter() - t0
    print(f"Done in {dt:.0f}s.", flush=True)
    out_path = os.path.join(HERE, 'results.csv')
    keys = sorted({k for r in results for k in r})
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(results)
    print(f"Wrote {len(results)} rows -> {out_path}", flush=True)


if __name__ == '__main__':
    main()
