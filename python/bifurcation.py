"""Bifurcation-diagram harness for EvoCA global dynamics.

Research-directions R3, program (c): bifurcations in the temporal
dynamics of global observables as a metaparameter is varied
(fixed point -> limit cycle -> period-doubling cascade -> chaos),
with the Feigenbaum constant delta ~ 4.6692 as the canonical scaling.

This is ANALYSIS only -- no C changes.  The driver sweeps one
metaparameter over a fine grid at fixed N and fixed seed, runs each
configuration long enough to discard a transient, then reduces the
post-transient scalar time series to its *recurrent set* (the local
extrema -- the standard bifurcation-diagram reduction).  The result is
a list of (param_value, recurrent_points) suitable for a scatter plot.

The thin per-tick sim loop here (`scalar_series`) deliberately mirrors
`evoca_explore.run_sim`'s init sequence (same `state`, `set_v`,
`set_f_all`, `set_F_all`, alive init, `set_seed`) so the dynamics match
the rest of the codebase -- but it records the chosen observable *every
tick* rather than at `sample_every`, which a bifurcation reduction
needs.  We do not modify `run_sim`.

Observables (global scalars, one per tick):
  - ``pop``      : number of alive cells
  - ``activity`` : total cumulative LUT activity (ts `activity_flux`-like)
  - ``lut_div``  : number of distinct living LUT genomes (diversity)

CLI / quick example (small N, coarse sweep -- seconds, not a campaign):

    python3 python/bifurcation.py --demo

obstacles/decisions:
  * EvoCA's per-tick activity getter (`get_activity`) is O(table); to
    keep sweeps tractable the activity/lut_div observables sample the
    table each tick but with a modest `max_n`.  `pop` is the cheapest
    and the recommended default observable for quick diagrams.
  * The Feigenbaum estimator is best-effort: it only returns a
    meaningful number if a *clean* period-doubling cascade was supplied;
    it is documented as such and never fabricated for noisy input.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from evoca_py import EvoCA  # noqa: E402

# Reuse the canonical parameter-key set so callers can sweep any
# metaparameter run_sim understands.
from evoca_explore import _PARAM_KEYS  # noqa: E402

# Default sweep candidates, per research-directions R3 (c).
DEFAULT_SWEEP_PARAMS = ('food_inc', 'mu_lut', 'tax')


# ----------------------------------------------------------------------
# Thin per-tick sim loop
# ----------------------------------------------------------------------

def _observe(sim, observable, N, max_n):
    """Return one global scalar for the current sim state."""
    if observable == 'pop':
        return float((sim.get_alive() == 1).sum())
    if observable == 'activity':
        sim._lib.evoca_activity_update()
        a = sim.get_activity(max_n=max_n)
        live = a['pop_count'] > 0
        return float(a['activity'][live].sum()) if live.any() else 0.0
    if observable == 'lut_div':
        sim._lib.evoca_activity_update()
        a = sim.get_activity(max_n=max_n)
        return float((a['pop_count'] > 0).sum())
    raise ValueError(f"unknown observable: {observable!r}")


def scalar_series(params, *, observable='pop', n_steps=2000, N=64,
                  seed=0, init='halfplane', max_n=20000):
    """Run one EvoCA config headless; return the per-tick scalar series.

    Mirrors `evoca_explore.run_sim`'s init sequence exactly (so the
    dynamics are identical to the rest of the codebase) but records
    `observable` *every tick* instead of at `sample_every`.  `run_sim`
    is left untouched.

    Returns a float numpy array of length ``n_steps + 1`` (state is
    recorded before each step and once after the last).
    """
    rng = np.random.default_rng(seed)

    sim = EvoCA()
    try:
        sim.init(N, **{k: params[k] for k in _PARAM_KEYS if k in params})
        sim.set_seed(seed)
        sim.state(lut='gol', egenome='uniform', egenome_value=0b000011)
        sim.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
        sim.set_f_all(0.5)
        sim.set_F_all(1.0)
        if init == 'halfplane':
            sim.set_alive_halfplane(0)
        elif init == 'fraction':
            sim.set_alive_fraction(0.5)
        else:
            raise ValueError(f"unknown init: {init}")

        out = np.empty(n_steps + 1, dtype=float)
        for t in range(n_steps + 1):
            out[t] = _observe(sim, observable, N, max_n)
            if t < n_steps:
                sim.step()
        return out
    finally:
        try:
            sim.free()
        except Exception:
            pass


# ----------------------------------------------------------------------
# Recurrent-set reduction (the bifurcation-diagram reducer)
# ----------------------------------------------------------------------

def recurrent_extrema(series, *, transient_frac=0.5, transient=None,
                      atol=1e-9, rtol=1e-6, max_points=512):
    """Reduce a scalar time series to its recurrent set: the local
    extrema *after* a discarded transient.

    This is the standard bifurcation-diagram reduction.  For a fixed
    point the post-transient signal is (near-)constant and collapses to
    a single value; for a period-2 cycle it yields the two alternating
    extrema; for a period-doubling cascade the number of distinct
    extrema doubles; for chaos it fills an interval.

    Parameters
    ----------
    series : 1-D array of the observable, one sample per tick.
    transient : int, number of leading samples to discard.  If None,
        ``int(transient_frac * len(series))`` is used.
    atol, rtol : tolerances for collapsing near-equal extrema (so a
        genuine fixed point returns exactly one point, not numerical
        dust).
    max_points : safety cap on returned points (chaotic series are
        decimated, not truncated, so the spread is preserved).

    Returns a 1-D float array of recurrent observable values (the y
    coordinates of one column of the bifurcation diagram).  A flat /
    constant tail returns its single steady value.
    """
    s = np.asarray(series, dtype=float).ravel()
    n = s.size
    if n == 0:
        return np.empty(0, dtype=float)

    if transient is None:
        transient = int(transient_frac * n)
    transient = max(0, min(transient, n - 1))
    tail = s[transient:]
    if tail.size == 0:
        tail = s[-1:]

    # Degenerate / constant tail -> single recurrent value (fixed point).
    spread = float(tail.max() - tail.min())
    scale = max(abs(float(tail.mean())), 1.0)
    if tail.size < 3 or spread <= atol + rtol * scale:
        return np.array([float(tail.mean())])

    # Interior local maxima and minima (turning points).  Plateaus are
    # handled by treating a non-strict change after a strict one as an
    # extremum boundary; for bifurcation purposes the standard reduction
    # uses strict turning points, which is sufficient for clean cycles
    # and still fills an interval under chaos.
    d = np.diff(tail)
    # sign of slope, with zeros carried forward so flat steps don't
    # spuriously register as turning points.
    sign = np.sign(d)
    nz = sign != 0
    if not nz.any():
        return np.array([float(tail.mean())])
    # forward-fill zero slopes with the last non-zero slope sign
    idx = np.where(nz, np.arange(sign.size), 0)
    np.maximum.accumulate(idx, out=idx)
    filled = sign[idx]
    turn = np.where(np.diff(filled) != 0)[0] + 1  # indices into tail
    extrema = tail[turn]

    if extrema.size == 0:
        return np.array([float(tail.mean())])

    # Collapse near-equal extrema (cluster) so a period-k cycle returns
    # exactly k points rather than k +/- numerical noise.
    vals = np.sort(extrema)
    tol = atol + rtol * max(abs(vals).max(), 1.0)
    # also use a fraction of the overall spread so genuine but tiny
    # limit-cycle amplitudes still separate, while jitter collapses.
    tol = max(tol, 1e-3 * spread)
    # Single-linkage from each cluster's anchor (first member): a value
    # joins the current cluster only if it is within tol of the anchor,
    # so a genuine period-k spacing wider than tol always splits.
    members = [[vals[0]]]
    for v in vals[1:]:
        if v - members[-1][0] <= tol:
            members[-1].append(v)
        else:
            members.append([v])
    reduced = np.array([np.mean(m) for m in members], dtype=float)

    if reduced.size > max_points:
        sel = np.linspace(0, reduced.size - 1, max_points).astype(int)
        reduced = reduced[np.unique(sel)]
    return reduced


# ----------------------------------------------------------------------
# Sweep driver
# ----------------------------------------------------------------------

def sweep(param, values, *, base_params=None, observable='pop',
          n_steps=2000, N=64, seed=0, init='halfplane',
          transient_frac=0.5, max_n=20000, progress=False):
    """Sweep one metaparameter over `values`; return the bifurcation
    data.

    Returns a list of dicts, one per swept value::

        {'value': float, 'points': np.ndarray, 'n_recurrent': int,
         'mean': float, 'series_tail_std': float}

    `points` is the recurrent set (y values for that column of the
    bifurcation diagram).  Everything is computed at fixed N and fixed
    seed so the only varying input is `param`.
    """
    if param not in _PARAM_KEYS:
        raise ValueError(
            f"{param!r} not a sweepable metaparameter; one of {_PARAM_KEYS}")
    base = dict(base_params or {})
    results = []
    for i, v in enumerate(values):
        p = dict(base)
        p[param] = v
        s = scalar_series(p, observable=observable, n_steps=n_steps,
                          N=N, seed=seed, init=init, max_n=max_n)
        pts = recurrent_extrema(s, transient_frac=transient_frac)
        transient = int(transient_frac * s.size)
        tail = s[transient:]
        results.append({
            'value': float(v),
            'points': pts,
            'n_recurrent': int(pts.size),
            'mean': float(tail.mean()) if tail.size else 0.0,
            'series_tail_std': float(tail.std()) if tail.size else 0.0,
        })
        if progress:
            print(f"  [{i + 1}/{len(values)}] {param}={v:.6g} "
                  f"-> {pts.size} recurrent pts "
                  f"(tail std {results[-1]['series_tail_std']:.4g})",
                  flush=True)
    return results


def to_scatter(results):
    """Flatten sweep results into (xs, ys) arrays for a scatter plot."""
    xs, ys = [], []
    for r in results:
        for y in r['points']:
            xs.append(r['value'])
            ys.append(y)
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


# ----------------------------------------------------------------------
# Feigenbaum-ratio estimator (best-effort)
# ----------------------------------------------------------------------

def feigenbaum_delta(thresholds):
    """Best-effort Feigenbaum-delta estimate from a period-doubling
    cascade.

    Given the parameter values ``p_1, p_2, p_3, ...`` at which successive
    period-doublings occur (1->2, 2->4, 4->8, ...), the Feigenbaum
    constant is the limit of the ratio of successive *gaps*::

        delta_k = (p_k - p_{k-1}) / (p_{k+1} - p_k)   ->   4.66920...

    This estimator returns the per-step ratios and their last value /
    mean.  It is **best-effort and only meaningful if a clean cascade
    was supplied** -- it cannot detect whether the input really is a
    period-doubling route to chaos; that is the caller's responsibility
    (e.g. by inspecting `n_recurrent` doubling: 1,2,4,8,... across the
    sweep).  With <3 thresholds no ratio is defined.

    Returns a dict::

        {'ratios': np.ndarray, 'delta_estimate': float|None,
         'delta_mean': float|None, 'n_thresholds': int}
    """
    th = np.asarray(thresholds, dtype=float).ravel()
    th = th[np.isfinite(th)]
    n = th.size
    if n < 3:
        return {'ratios': np.empty(0), 'delta_estimate': None,
                'delta_mean': None, 'n_thresholds': int(n)}
    gaps = np.diff(th)               # p_{k} - p_{k-1}
    # ratio_k = gap_k / gap_{k+1}
    with np.errstate(divide='ignore', invalid='ignore'):
        ratios = gaps[:-1] / gaps[1:]
    ratios = ratios[np.isfinite(ratios)]
    if ratios.size == 0:
        return {'ratios': ratios, 'delta_estimate': None,
                'delta_mean': None, 'n_thresholds': int(n)}
    return {
        'ratios': ratios,
        'delta_estimate': float(ratios[-1]),   # best as cascade deepens
        'delta_mean': float(ratios.mean()),
        'n_thresholds': int(n),
    }


def doubling_thresholds(results):
    """Heuristic: from sweep results, find the param values where the
    recurrent-set size first reaches 1, 2, 4, 8, ... (a candidate
    period-doubling cascade).  Best-effort -- returns whatever clean
    doublings it can see; feed the result to `feigenbaum_delta`.
    """
    th = []
    target = 1
    for r in results:
        # advance through powers of two as the recurrent set grows
        while r['n_recurrent'] >= target * 2 and target < 1 << 20:
            th.append(r['value'])
            target *= 2
    return np.asarray(th, dtype=float)


# ----------------------------------------------------------------------
# CLI / quick example
# ----------------------------------------------------------------------

def _demo():
    """Tiny, fast diagram (small N, coarse sweep) so the output format
    can be eyeballed.  NOT a campaign -- seconds, not hours."""
    param = 'food_inc'
    values = np.linspace(0.02, 0.20, 8)
    print(f"[demo] sweeping {param} over {len(values)} values, "
          f"N=48, n_steps=400, seed=0 (fast, illustrative only)")
    results = sweep(param, values, observable='pop', n_steps=400,
                    N=48, seed=0, transient_frac=0.5, progress=True)
    xs, ys = to_scatter(results)
    print(f"\n[demo] bifurcation scatter: {xs.size} (param, observable) "
          f"points across {len(results)} columns")
    print(f"[demo] {param:>10s} | n_recurrent | tail_mean | tail_std")
    for r in results:
        print(f"        {r['value']:10.4f} | {r['n_recurrent']:11d} | "
              f"{r['mean']:9.2f} | {r['series_tail_std']:.4g}")
    th = doubling_thresholds(results)
    fb = feigenbaum_delta(th)
    print(f"\n[demo] candidate doubling thresholds: {th.tolist()}")
    print(f"[demo] Feigenbaum (best-effort, needs a CLEAN cascade): "
          f"delta_estimate={fb['delta_estimate']} "
          f"delta_mean={fb['delta_mean']}")
    print("[demo] done (this coarse run is for format inspection, "
          "not for drawing scientific conclusions).")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="EvoCA bifurcation-diagram harness (research R3 c).")
    ap.add_argument('--demo', action='store_true',
                    help="run a tiny fast illustrative sweep and exit")
    ap.add_argument('--param', default='food_inc',
                    choices=_PARAM_KEYS,
                    help="metaparameter to sweep")
    ap.add_argument('--observable', default='pop',
                    choices=('pop', 'activity', 'lut_div'))
    ap.add_argument('--lo', type=float, default=0.02)
    ap.add_argument('--hi', type=float, default=0.20)
    ap.add_argument('--steps', type=int, default=8,
                    help="number of swept parameter values")
    ap.add_argument('--n-steps', type=int, default=600,
                    help="ticks per configuration")
    ap.add_argument('--N', type=int, default=48)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--transient-frac', type=float, default=0.5)
    args = ap.parse_args(argv)

    if args.demo:
        _demo()
        return 0

    values = np.linspace(args.lo, args.hi, args.steps)
    print(f"sweeping {args.param} in [{args.lo}, {args.hi}] "
          f"({args.steps} pts), observable={args.observable}, "
          f"N={args.N}, n_steps={args.n_steps}, seed={args.seed}")
    results = sweep(args.param, values, observable=args.observable,
                    n_steps=args.n_steps, N=args.N, seed=args.seed,
                    transient_frac=args.transient_frac, progress=True)
    xs, ys = to_scatter(results)
    print(f"\nbifurcation scatter: {xs.size} points across "
          f"{len(results)} columns")
    th = doubling_thresholds(results)
    fb = feigenbaum_delta(th)
    print(f"candidate doubling thresholds: {th.tolist()}")
    print(f"Feigenbaum (best-effort): {fb}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
