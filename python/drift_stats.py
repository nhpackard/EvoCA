"""Drift-shadow statistics (Docs/shadows.md §6).

Input: the cumulative activities returned by EvoCA.fs_drift_get() (or the
synthetic driver python/drift_synth.py, which returns the same layout):
one real array and R replicate arrays per channel. Every statistic is
computed for the real run AND for each replicate scored against the
other R units (the real run included), so the R+1 units are exchangeable
under H0 and the rank p-value

    p = (1 + #{r : stat_r >= stat_real}) / (R + 1)

is exact if the shadow is the right neutral model. `*_z` is the real
statistic standardized by the replicates' statistics:
(stat_real − mean_r stat_r) / sd_r stat_r.

Channels and statistics (prefix `<ch>_drift_`):

  eg, marker (729 ternary keys)
    concentration  log Simpson(A_real) − mean_r log Simpson(A_r)
                   (intensive; > 0: real activity more concentrated than
                   neutral diffusion from the founder)
  dyn (500 (input, output) components; real side = whole-population
       allele frequency x real query counts, §5.4 design C)
    consistency    Z_dyn = Σ_i w_i ζ_i / sqrt(Σ w_i²), ζ_i the
                   standardized cumulative per-input homozygosity
                   F_i = (A_i0² + A_i1²)/(A_i0 + A_i1)², w_i = input i's
                   share of query mass (primary dyn statistic)
    concentration  as for eg, over the 500 components (secondary)
  all channels
    tv             ½ Σ_k |A_real,k − mean_r A_r,k| / Σ A_real
    tv_z, tv_p     two-sided departure test on tv (each unit's TV from
                   the mean of the other R). Needed because the one-sided
                   statistics above score a directional sweep that
                   completes inside the window as NEGATIVE: the sweep
                   splits cumulative activity between the old and new
                   component, lowering cumulative homozygosity/Simpson.
    nhi            #{k: z_k > 3} minus the replicates' mean count
    frac_z3        #{k: |z_k| > 3} / #present components (calibration)
  plus
    negene_drift          time-mean real Negene − replicate mean
    dyn_expression_bias   ½ Σ_k |Q_k − A_real,k| / M, Q = real queried
                          transitions (design B), A_real = design C;
                          the dynamical self-selection the C real side
                          removes (descriptor, not a test)
    dyn_expression_bias_signed  Σ_i (Q_i1 − A_real,i1) / M

Per-component z_k = (A_real,k − μ_k)/σ_k with σ_k = max(sd_k, 1,
0.01·μ_k) (doc §6.2 floor). For the per-input homozygosity the analogous
floor is max(sd, 1e-6, 0.01·(1 − mean F)).
"""
import numpy as np

KEYS_CONC = ('concentration', 'z', 'p', 'tv', 'tv_z', 'tv_p', 'nhi', 'frac_z3')


def _safe_div(num, den):
    """num/den; 0 where both are 0, ±inf where only den is 0."""
    num, den = np.broadcast_arrays(np.asarray(num, dtype=float),
                                   np.asarray(den, dtype=float))
    out = np.zeros(num.shape)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    pos = ~ok & (num > 0)
    neg = ~ok & (num < 0)
    out[pos] = np.inf
    out[neg] = -np.inf
    return out


def _loo(U):
    """Leave-one-out reference mean and sd (ddof=1) for each unit.
    U: (R+1, ...) → (mean, sd) of the same shape."""
    n = U.shape[0]
    S1 = U.sum(axis=0)
    S2 = (U * U).sum(axis=0)
    m = (S1[None] - U) / (n - 1)
    var = ((S2[None] - U * U) - (n - 1) * m * m) / (n - 2)
    return m, np.sqrt(np.maximum(var, 0.0))


def _rank_stats(stat):
    """stat[0] = real, stat[1:] = replicates (each scored vs the rest)."""
    real, reps = float(stat[0]), np.asarray(stat[1:], dtype=float)
    R = reps.size
    p = (1.0 + float(np.sum(reps >= real))) / (R + 1.0)
    sd = reps.std(ddof=1) if R > 1 else 0.0
    z = float(_safe_div(real - reps.mean(), sd))
    return z, p


def _log_simpson(U):
    tot = U.sum(axis=-1)
    with np.errstate(divide='ignore', invalid='ignore'):
        s = (U * U).sum(axis=-1) / (tot * tot)
    return np.log(np.where(tot > 0, s, 1.0))


def _component_tails(U):
    """Per-component z of each unit vs the others; returns
    (count z>3 per unit, frac |z|>3 of the real, z of the real)."""
    m, sd = _loo(U)
    sig = np.maximum(np.maximum(sd, 1.0), 0.01 * m)
    z = (U - m) / sig
    present = (U[0] > 0) | (m[0] > 0)
    hi = (z > 3).sum(axis=1)
    npres = max(int(present.sum()), 1)
    frac = float((np.abs(z[0]) > 3)[present].sum()) / npres
    return hi, frac, z[0]


def concentration_stats(real, shadow):
    """eg-type statistics for one channel. real (K,), shadow (R, K)."""
    U = np.vstack([np.asarray(real, float)[None], np.asarray(shadow, float)])
    R = U.shape[0] - 1
    out = {}
    ls = _log_simpson(U)
    if R >= 2:
        stat = ls - (ls.sum() - ls) / R
        out['concentration'] = float(stat[0])
        out['z'], out['p'] = _rank_stats(stat)
        hi, frac, _ = _component_tails(U)
        out['nhi'] = float(hi[0] - hi[1:].mean())
        out['frac_z3'] = frac
    else:
        out['concentration'] = float(ls[0] - ls[1:].mean()) if R else float('nan')
        out['z'] = out['p'] = out['nhi'] = out['frac_z3'] = float('nan')
    mu = U[1:].mean(axis=0) if R else np.zeros_like(U[0])
    tot = U[0].sum()
    out['tv'] = float(0.5 * np.abs(U[0] - mu).sum() / tot) if tot > 0 else 0.0
    if R >= 2:
        # Two-sided departure test: each unit's TV distance from the mean
        # of the other R units, ranked. Catches directional sweeps, which
        # the one-sided concentration / consistency statistics score
        # negative (a sweep inside the window splits cumulative activity).
        ref = (U.sum(axis=0)[None] - U) / R
        tots = U.sum(axis=1)
        tvu = np.where(tots > 0, 0.5 * np.abs(U - ref).sum(axis=1)
                       / np.where(tots > 0, tots, 1.0), 0.0)
        out['tv_z'], out['tv_p'] = _rank_stats(tvu)
    else:
        out['tv_z'] = out['tv_p'] = float('nan')
    return out


def consistency_stats(real, shadow):
    """dyn per-input consistency Z_dyn. real (500,), shadow (R, 500)
    laid out as key = 2*i + o. Returns (stats, F[R+1, 250], zeta_real[250],
    w[250]): F the cumulative per-input homozygosity of every unit, zeta
    the real run's standardized per-input deviation, w the query shares."""
    U = np.vstack([np.asarray(real, float)[None], np.asarray(shadow, float)])
    R = U.shape[0] - 1
    A0, A1 = U[:, 0::2], U[:, 1::2]
    mass = A0 + A1
    with np.errstate(divide='ignore', invalid='ignore'):
        F = np.where(mass > 0, (A0 * A0 + A1 * A1) / (mass * mass), 1.0)
    w = mass[0] / mass[0].sum() if mass[0].sum() > 0 else np.zeros(mass.shape[1])
    out = {}
    if R < 2 or not np.any(w > 0):
        out['consistency'] = out['z'] = out['p'] = float('nan')
        return out, F, np.zeros(F.shape[1]), w
    m, sd = _loo(F)
    sig = np.maximum(np.maximum(sd, 1e-6), 0.01 * (1.0 - m))
    zeta = (F - m) / sig
    Z = (zeta * w[None]).sum(axis=1) / np.sqrt((w * w).sum())
    out['consistency'] = float(Z[0])
    out['z'], out['p'] = _rank_stats(Z)
    return out, F, zeta[0], w


def drift_statistics(acc):
    """All drift-shadow run-level statistics from one accumulator dict
    (fs_drift_get layout). Returns a flat dict of `<ch>_drift_*` keys."""
    out = {}
    R = int(acc.get('R', np.asarray(acc['eg_shadow']).shape[0]))
    steps = int(acc.get('steps', 0))
    for ch in ('eg', 'marker'):
        st = concentration_stats(acc[f'{ch}_real'], acc[f'{ch}_shadow'])
        for k in KEYS_CONC:
            out[f'{ch}_drift_{k}'] = st[k]
    cons = consistency_stats(acc['dyn_real'], acc['dyn_shadow'])[0]
    out['dyn_drift_consistency'] = cons['consistency']
    out['dyn_drift_z'] = cons['z']
    out['dyn_drift_p'] = cons['p']
    st = concentration_stats(acc['dyn_real'], acc['dyn_shadow'])
    out['dyn_drift_concentration'] = st['concentration']
    out['dyn_drift_concentration_z'] = st['z']
    out['dyn_drift_concentration_p'] = st['p']
    for k in ('tv', 'tv_z', 'tv_p', 'nhi', 'frac_z3'):
        out[f'dyn_drift_{k}'] = st[k]
    if steps > 0 and R > 0:
        out['negene_drift'] = float(acc['negene_real'] / steps
                                    - np.mean(acc['negene_shadow']) / steps)
    else:
        out['negene_drift'] = float('nan')
    Q = acc.get('dyn_queried')
    G = np.asarray(acc['dyn_real'], float)
    M = G.sum()
    if Q is not None and M > 0:
        Q = np.asarray(Q, float)
        out['dyn_expression_bias'] = float(0.5 * np.abs(Q - G).sum() / M)
        out['dyn_expression_bias_signed'] = float((Q[1::2] - G[1::2]).sum() / M)
    else:
        out['dyn_expression_bias'] = out['dyn_expression_bias_signed'] = float('nan')
    # Q9 check: the same consistency statistic with RUN-MEAN query weights
    # n̄_i (fixed-weight locus test). A material difference from the
    # contemporaneous version flags a Cov_t(n_i, P_i) term (§5.4).
    if acc.get('psum_real') is not None and steps > 0:
        nbar = (G[0::2] + G[1::2]) / steps
        def _rm(ps):
            ps = np.asarray(ps, float)
            o = np.zeros(ps.shape[:-1] + (500,))
            o[..., 1::2] = nbar * ps
            o[..., 0::2] = nbar * (steps - ps)
            return o
        rm = consistency_stats(_rm(acc['psum_real']), _rm(acc['psum_shadow']))[0]
        out['dyn_drift_consistency_runmean'] = rm['consistency']
        out['dyn_drift_z_runmean'] = rm['z']
    out['drift_steps'] = steps
    out['drift_R'] = R
    return out


_ACC_ARRAYS = ('eg_real', 'eg_shadow', 'dyn_real', 'dyn_queried', 'dyn_shadow',
               'marker_real', 'marker_shadow', 'negene_shadow', 'psum_real',
               'psum_shadow')


def window(acc_end, acc_start):
    """Accumulators over (start, end]: elementwise difference."""
    w = {}
    for k in _ACC_ARRAYS:
        if acc_end.get(k) is not None:
            w[k] = np.asarray(acc_end[k]) - np.asarray(acc_start[k])
    w['negene_real'] = acc_end['negene_real'] - acc_start['negene_real']
    w['steps'] = acc_end['steps'] - acc_start['steps']
    w['R'] = acc_end['R']
    return w
