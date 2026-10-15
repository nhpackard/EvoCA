"""Campaign #8 — direct lineage co-evolution metric.

Spec: Docs/plan-8-lineage-coevolution.md.

ΔLUT and Δegene are drawn by independent Poisson passes at the
mutation step, so they are independent across reproduction edges *by
construction*. Any (ΔLUT,Δegene) association among SURVIVING edges is
selection jointly retaining co-adapted LUT+egene — direct, shadow-free
co-evolution evidence (immune to the #3c viability + dyn-shadow
confounds).

Per run: lineage ON, warm-up, then snapshot (birth_id, parent_id,
lut-hash, egene-key, alive) every S ticks over a post-transient
window; rebuild parent->child edges; 2x2 (ΔLUT × Δegene) over
surviving edges; phi / odds / excess vs independence + a within-run
label-permutation null. Arms: co-evo (mu_lut>0) and frozen-LUT
(mu_lut=0, ΔLUT≡0 — the clean "3b cannot show it" control).

Pure Python on main (analysis-only, no C change). Launched detached
on torque (nohup; SSH-only health check — the #4/R1/#3c lesson).
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

from evoca_py import EvoCA                                  # noqa: E402
from evoca_explore import save_scan_config                  # noqa: E402

N_GRID   = 256
T_WARM   = 4000
T_MEAS   = 6000
S        = 50            # snapshot interval
T_SURV   = 200           # min observed lifespan for a "surviving" edge
N_PERM   = 200
SEEDS    = tuple(range(8))
REGIME   = dict(m_scale=3.0, food_inc=0.015, tax=0.035, gdiff=0.06,
                mu_egene=0.01, mu_egenome=0.01, tax_per_egene=0.0)
ARMS     = {'coevo': 0.03, 'frozenLUT': 0.0}     # phase mu_lut
N_WORKERS = 16
_SLOT = np.arange(8, dtype=np.uint8)


def _egene_keys(sim):
    """Per-cell FNV-1a over active_byte ‖ active-slot (value,mask)
    bytes (inactive slots zeroed → key invariant to junk)."""
    act = sim.get_active().ravel().astype(np.uint8)              # (NN,)
    val = sim.get_egenes().reshape(-1, 8).astype(np.uint8)       # (NN,8)
    msk = sim.get_egenes_mask().reshape(-1, 8).astype(np.uint8)
    on = ((act[:, None] >> _SLOT) & 1).astype(bool)              # (NN,8)
    val = np.where(on, val, 0)
    msk = np.where(on, msk, 0)
    h = np.full(act.shape, 2166136261, dtype=np.uint64)
    for col in (act, *val.T, *msk.T):
        h = ((h ^ col.astype(np.uint64)) * np.uint64(16777619)) \
            & np.uint64(0xFFFFFFFF)
    return h


def _phi(dl, de):
    """phi coefficient of the 2x2 (ΔLUT × Δegene). 0 if degenerate."""
    n11 = int(np.sum(dl & de)); n10 = int(np.sum(dl & ~de))
    n01 = int(np.sum(~dl & de)); n00 = int(np.sum(~dl & ~de))
    a = (n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00)
    if a == 0:
        return 0.0, (n11, n10, n01, n00)
    return (n11 * n00 - n10 * n01) / np.sqrt(a), (n11, n10, n01, n00)


def _run_one(arg):
    cfg_idx, (arm, seed) = arg
    try:
        mu_lut = ARMS[arm]
        sim = EvoCA()
        # lineage MUST be passed through init() — init() resets the C
        # flag from this kwarg (a separate enable_lineage() pre-call is
        # clobbered).
        sim.init(N_GRID, mu_lut=mu_lut, lineage=True, **REGIME)
        sim.set_seed(seed)
        sim.state(lut='gol', egenome='uniform', egenome_value=0b000011)
        rng = np.random.default_rng(seed)
        sim.set_v(rng.integers(0, 2, (N_GRID, N_GRID), dtype=np.uint8))
        sim.set_f_all(0.5)
        sim.set_F_all(1.0)
        sim.set_alive_halfplane(0)

        for _ in range(T_WARM):
            sim.step()

        # bid -> [lut_hash, egene_key, parent_id, first_t, last_t]
        rec = {}
        t = T_WARM
        for _ in range(T_MEAS // S):
            for _ in range(S):
                sim.step()
            t += S
            alive = sim.get_alive().ravel() == 1
            if not alive.any():
                break
            bid = sim.get_lineage_birth_id().ravel()
            pid = sim.get_lineage_parent_id().ravel()
            lh = sim.get_lut_hashes().ravel()
            ek = _egene_keys(sim)
            ai = np.flatnonzero(alive)
            for i in ai:
                b = int(bid[i])
                if b == 0:
                    continue
                r = rec.get(b)
                if r is None:
                    rec[b] = [int(lh[i]), int(ek[i]), int(pid[i]), t, t]
                else:
                    r[4] = t
        sim.free()

        # ── build parent->child edges ──
        dl, de, surv = [], [], []
        for b, (lu, eg, pa, f_t, l_t) in rec.items():
            if pa == 0:
                continue
            pr = rec.get(pa)
            if pr is None:
                continue
            dl.append(lu != pr[0])
            de.append(eg != pr[1])
            surv.append((l_t - f_t) >= T_SURV)
        dl = np.array(dl, bool)
        de = np.array(de, bool)
        surv = np.array(surv, bool)
        n_edges = dl.size
        out = dict(config_idx=cfg_idx, arm=arm, seed=seed,
                   mu_lut=mu_lut, n_bids=len(rec), n_edges=n_edges,
                   n_surv_edges=int(surv.sum()), error='')
        if n_edges == 0 or surv.sum() < 10:
            out['error'] = 'too_few_edges'
            return out

        # all-edge marginals (≈ independent mutation step) for reference
        out['p_dLUT_all'] = float(dl.mean())
        out['p_dEG_all'] = float(de.mean())
        out['p_dLUT_analytic'] = float(1.0 - np.exp(-mu_lut * 250.0))

        sdl, sde = dl[surv], de[surv]
        phi, (n11, n10, n01, n00) = _phi(sdl, sde)
        pL, pE = sdl.mean(), sde.mean()
        p11 = (sdl & sde).mean()
        out.update(surv_p_dLUT=float(pL), surv_p_dEG=float(pE),
                   surv_p11=float(p11),
                   surv_excess_joint=float(p11 - pL * pE),
                   surv_phi=float(phi),
                   surv_n11=n11, surv_n10=n10, surv_n01=n01,
                   surv_n00=n00)
        # odds ratio (Haldane-Anscombe 0.5 correction)
        out['surv_odds_ratio'] = float(
            ((n11 + .5) * (n00 + .5)) / ((n10 + .5) * (n01 + .5)))

        # label-permutation null on phi (cancels survivor-capture bias)
        if pL in (0.0, 1.0) or pE in (0.0, 1.0):
            out['phi_perm_z'] = 0.0
            out['phi_perm_p'] = 1.0
            out['frozen_null'] = (mu_lut == 0.0)
        else:
            pr = np.random.default_rng(seed)
            null = np.empty(N_PERM)
            de_s = sde.copy()
            for k in range(N_PERM):
                pr.shuffle(de_s)
                null[k] = _phi(sdl, de_s)[0]
            mu, sd = null.mean(), null.std()
            out['phi_perm_z'] = float((phi - mu) / sd) if sd > 0 else 0.0
            out['phi_perm_p'] = float((np.sum(null >= phi) + 1)
                                      / (N_PERM + 1))
            out['frozen_null'] = False
        return out
    except Exception as e:                                  # noqa: BLE001
        return dict(config_idx=cfg_idx, arm=arg[1][0],
                    seed=arg[1][1], error=repr(e))


def main():
    save_scan_config(HERE, N=N_GRID, n_steps=T_WARM + T_MEAS,
                     sample_every=S, seed=0, shadow=False)
    jobs = [(a, s) for a in ARMS for s in SEEDS]
    args = list(enumerate(jobs))
    print(f"[#8 lineage-coevolution] {len(jobs)} runs "
          f"(2 arms x {len(SEEDS)} seeds) / {N_WORKERS} workers / "
          f"warm={T_WARM} meas={T_MEAS} S={S} @ N={N_GRID}", flush=True)
    t0 = time.perf_counter()
    results = []
    with Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(_run_one, args):
            results.append(r)
            d = len(results)
            print(f"  {d}/{len(args)}  {time.perf_counter()-t0:5.0f}s  "
                  f"{r.get('arm')}/{r.get('seed')} "
                  f"phi={r.get('surv_phi')} z={r.get('phi_perm_z')} "
                  f"err={r.get('error')!r}", flush=True)
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
