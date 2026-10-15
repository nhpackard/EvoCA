"""
dyn_activity as the PRIMARY phenotypic-evolution readout.

Checkpoint 3-July-2026 §F / §I' action D1: promote the *used-transition*
signal (``dyn_activity`` + fixed-space ``dyn_excess``) to the primary
phenotypic-evolution readout, and DEMOTE full-genome ``activity`` to a
churn diagnostic.

Why (§F, verbatim intent): the full-genome ``activity`` probe keys on an
FNV-1a hash of the *whole* 250-bit LUT, so flipping a bit for an input the
cell never encounters mints a "new species" without changing behaviour —
**silent churn**. ``dyn_activity`` counts only the 500 ``(LUT-input,
output)`` transitions that *actually fired* for alive cells, so it moves
only when the effective dynamics moves.

The KEY CONTRAST this module computes and names:

    full-genome activity churns  +  dyn_activity stays concentrated
        ==> SILENT CHURN   (bit-shuffling on unused entries; no
                            phenotypic change)

    full-genome activity churns  +  dyn_activity ALSO moves
        ==> FUNCTIONAL CHANGE  (the population is genuinely shifting
                               which transitions it uses / selection
                               on the effective rule)

The empirical lesson (measured, see dyn-promote-demo.ipynb): under any
LUT mutation the full-genome ``activity`` probe reports a continuous
parade of new species in BOTH the silent and the functional regime --
by that measure alone the two are indistinguishable. What separates
them is whether ``dyn_activity`` MOVES: the *effective number of
exercised transitions* (``dyn_eff``, a participation ratio over the 500
buckets) climbs when the change is functional and stays low & flat when
the churn is silent. That recruitment slope -- not the raw
``dyn_excess`` magnitude -- is the robust discriminator (``dyn_excess``
is reported but is artifact-prone in sparse / frozen-substrate regimes;
it corroborates, it does not gate).

This is pure re-interpretation of EXISTING probes (no C change): it reads
``sim.get_dyn_activity()``, ``sim.dyn_excess_pc()`` and
``sim.get_activity()`` that already ship.

Public API
----------
``dyn_report(params, ...)`` -> dict with time-series + summary + verdict.
``print_report(rep)``       -> human-readable one-screen digest.

Metric glossary (all in the returned dict)
------------------------------------------
CHURN track (full-genome ``activity`` -- demoted to diagnostic):
  ``genome_distinct``       D_G(t): distinct alive LUT-genome hashes.
  ``genome_cum_distinct``   cumulative distinct hashes ever seen.
  ``genome_turnover``       cum_distinct[-1] / mean(pop): how many times
                            over the standing population's worth of
                            distinct genomes the run has minted (>1 =>
                            the genome pool has cycled through more than
                            one population's worth of species: churn).
  ``genome_new_rate``       new distinct hashes per 1000 steps (slope of
                            cum_distinct): the PRIMARY churn magnitude --
                            > 0 means the activity probe keeps minting
                            species.

FUNCTION track (``dyn_activity`` -- promoted to primary):
  ``dyn_eff``               effective # used transitions carrying the
                            activity = participation ratio
                            (Sum a)^2 / Sum a^2 over the 500 buckets.
                            Concentration = LOW & FLAT dyn_eff.
  ``dyn_support``           # of the 500 transitions ever exercised.
  ``dyn_excess_pc``         fixed-space per-component neutral-drift
                            excess (observed - drift); >0 => selection
                            on used transitions, ~0 => drift.

The verdict combines a churn score (is full-genome activity turning
over?) with a function score (is dyn_activity moving / selected?).
"""

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from evoca_explore import build_sim   # noqa: E402


# ── concentration helpers ────────────────────────────────────────────

def participation_ratio(a):
    """Effective number of nonzero components of a nonnegative vector.

    PR = (Sum a)^2 / Sum a^2, in [1, support]. For dyn_activity's 500
    cumulative bucket counts this is the effective number of transitions
    that carry the population's activity -- a concentration measure that,
    unlike raw support, does NOT saturate once every bucket has fired
    once (it stays weighted by how often each transition is used)."""
    a = np.asarray(a, dtype=float)
    s2 = float((a * a).sum())
    if s2 <= 0.0:
        return 0.0
    return float(a.sum() ** 2 / s2)


# ── main entry point ─────────────────────────────────────────────────

def dyn_report(params, n_steps=2000, sample_every=50, N=96, seed=0,
               init='halfplane',
               churn_rate_hi=50.0, func_eff_level=60.0):
    """Run one config and report the dyn_activity-vs-activity contrast.

    params : dict of EvoCA metaparams (see evoca_explore._PARAM_KEYS).
    n_steps, sample_every, N, seed, init : run controls.

    Verdict thresholds (tune per study):
      churn_rate_hi   : genome_new_rate (new hashes/1000 steps) above
                        this => full-genome activity is minting species
                        (churn present).
      func_eff_level  : dyn_eff (effective transitions, of 500) above
                        this => the population exercises a SPREAD of the
                        transition space (functional). Below => the used
                        transitions are CONCENTRATED (silent churn).

    Why LEVEL, not slope, gates the verdict: dyn_eff's growth *rate*
    depends on where in the transient you measure (a functional run's
    slope is high early and ~0 once it plateaus), so a slope gate
    mislabels short runs. The converged *level* (how many of the 500
    transitions the population actually exercises) is horizon-robust:
    the silent regime settles near ~25/500, the functional near
    ~200/500. The recruitment slope IS reported (dyn_eff_slope) as a
    corroborating movement signal, and is what the demo figure shows
    visually (flat vs climbing) -- it just does not gate.

    ``dyn_excess_pc`` is likewise reported but NOT a gate: its magnitude
    is artifact-prone in sparse / frozen-substrate regimes (project
    memory: dyn-shadow frozen-rich-substrate artifact).

    NOTE: this level gate targets *recruitment*-type functional change
    (the population bringing new transitions into use). Selection that
    reshapes *which* of a small fixed transition set dominates -- without
    growing the count -- would keep dyn_eff low; that subtler case needs
    per-bucket composition tracking (checkpoint §G2) and is out of scope
    here.

    Returns a dict: params echo, 't' plus the time-series named in the
    module docstring, scalar summaries, and 'verdict' / 'verdict_reason'.
    """
    sim = build_sim(params, N=N, seed=seed, init=init)

    t_hist = []
    g_distinct, g_cum_distinct = [], []
    dyn_eff_hist, dyn_support_hist, dyn_excess_hist = [], [], []
    pop_hist = []
    seen_hashes = set()

    for t in range(n_steps + 1):
        if t % sample_every == 0:
            # ── FUNCTION track: dyn_activity (promoted) ──────────────
            d = sim.get_dyn_activity()
            a = d['activity'].astype(float)          # cumulative, len 500
            dyn_eff_hist.append(participation_ratio(a))
            dyn_support_hist.append(int((a > 0).sum()))
            dyn_excess_hist.append(sim.dyn_excess_pc())

            # ── CHURN track: full-genome activity (demoted) ──────────
            sim._lib.evoca_activity_update()
            G = sim.get_activity(max_n=200000)
            live = G['pop_count'] > 0
            n_live = int(live.sum())
            g_distinct.append(n_live)
            if n_live:
                seen_hashes.update(int(h) for h in G['hash'][live])
            g_cum_distinct.append(len(seen_hashes))

            alive = np.ctypeslib.as_array(
                sim._lib.evoca_get_alive(), shape=(N * N,))
            pop_hist.append(int((alive == 1).sum()))
            t_hist.append(t)

        if t < n_steps:
            sim.step()

    sim.free()

    t_arr = np.asarray(t_hist, dtype=float)
    pop_arr = np.asarray(pop_hist, dtype=float)
    mean_pop = float(pop_arr[pop_arr > 0].mean()) if (pop_arr > 0).any() \
        else 0.0
    final_pop = float(pop_arr[-1]) if pop_arr.size else 0.0
    min_pop = float(pop_arr.min()) if pop_arr.size else 0.0
    cum = np.asarray(g_cum_distinct, dtype=float)
    eff = np.asarray(dyn_eff_hist, dtype=float)

    def _slope_per_1000(y):
        if len(t_arr) >= 2 and np.ptp(t_arr) > 0:
            m = np.polyfit(t_arr, y, 1)[0]
            return float(m * 1000.0)
        return 0.0

    genome_turnover = float(cum[-1] / mean_pop) if mean_pop > 0 else 0.0
    genome_new_rate = _slope_per_1000(cum)
    dyn_eff_slope = _slope_per_1000(eff)
    dyn_excess_final = float(dyn_excess_hist[-1]) if dyn_excess_hist \
        else 0.0

    # ── verdict: cross the churn track with the function track ───────
    # Churn = full-genome activity keeps minting species (new_rate).
    # Recruitment = dyn_activity moves (effective transitions climbing).
    dyn_eff_final = float(eff[-1]) if eff.size else 0.0
    churning = genome_new_rate >= churn_rate_hi
    spread = dyn_eff_final >= func_eff_level      # transitions exercised

    # Extinction gate uses FINAL/MIN population, not mean: a run that
    # dies mid-way freezes its metrics (dyn_eff flat, cum-distinct
    # plateaued) and would otherwise be mis-labelled 'silent-churn'.
    if final_pop <= 0 or min_pop <= 0:
        verdict = 'extinct'
        reason = (f'population reached 0 (min_pop={min_pop:.0f}, '
                  f'final_pop={final_pop:.0f}); metrics frozen post-death, '
                  f'no valid readout')
    elif churning and not spread:
        verdict = 'silent-churn'
        reason = (f'activity minting {genome_new_rate:.0f} species/kstep '
                  f'(churning) but dyn_activity concentrated '
                  f'({dyn_eff_final:.0f}/500 effective < {func_eff_level:.0f}'
                  f', slope {dyn_eff_slope:+.1f}/kstep): genome hashes '
                  f'churn while used transitions stay put => SILENT churn')
    elif churning and spread:
        verdict = 'functional'
        reason = (f'activity minting {genome_new_rate:.0f} species/kstep '
                  f'AND dyn_activity spread ({dyn_eff_final:.0f}/500 '
                  f'effective >= {func_eff_level:.0f}, recruited at '
                  f'{dyn_eff_slope:+.1f}/kstep): both move => real '
                  f'functional change')
    elif spread:
        verdict = 'functional-low-churn'
        reason = (f'dyn_activity spread ({dyn_eff_final:.0f}/500 effective) '
                  f'with low genome churn ({genome_new_rate:.0f} '
                  f'species/kstep): transitions exercised without '
                  f'genome-hash churn')
    else:
        verdict = 'quiescent'
        reason = (f'neither track moving (churn {genome_new_rate:.0f}'
                  f'/kstep, dyn_eff {dyn_eff_final:.0f}/500): '
                  f'effectively static')

    return {
        'params': dict(params),
        'N': N, 'n_steps': n_steps, 'sample_every': sample_every,
        'seed': seed, 'init': init,
        # time series
        't': t_arr,
        'genome_distinct': np.asarray(g_distinct, dtype=float),
        'genome_cum_distinct': cum,
        'dyn_eff': eff,
        'dyn_support': np.asarray(dyn_support_hist, dtype=float),
        'dyn_excess_pc': np.asarray(dyn_excess_hist, dtype=float),
        'pop': pop_arr,
        # scalar summary
        'mean_pop': mean_pop,
        'final_pop': final_pop,
        'min_pop': min_pop,
        'genome_turnover': genome_turnover,
        'genome_new_rate': genome_new_rate,
        'dyn_eff_final': dyn_eff_final,
        'dyn_eff_slope': dyn_eff_slope,
        'dyn_support_final': float(dyn_support_hist[-1])
            if dyn_support_hist else 0.0,
        'dyn_excess_pc_final': dyn_excess_final,
        # thresholds used
        'churn_rate_hi': churn_rate_hi,
        'func_eff_level': func_eff_level,
        # verdict
        'verdict': verdict,
        'verdict_reason': reason,
    }


def print_report(rep):
    """One-screen human-readable digest of a dyn_report() result."""
    p = rep['params']
    print(f"config: {p}")
    print(f"  run: N={rep['N']} steps={rep['n_steps']} seed={rep['seed']} "
          f"pop mean/min/final="
          f"{rep['mean_pop']:.0f}/{rep['min_pop']:.0f}/{rep['final_pop']:.0f}")
    print("  CHURN  (full-genome activity -- demoted to diagnostic):")
    print(f"    genome_new_rate     = {rep['genome_new_rate']:.0f} "
          f"new species / 1000 steps  (>0 => activity churning)")
    print(f"    genome_turnover     = {rep['genome_turnover']:.2f} "
          f"population-worths of distinct genomes minted")
    print("  FUNCTION  (dyn_activity -- PRIMARY readout):")
    print(f"    dyn_eff (final)     = {rep['dyn_eff_final']:.1f} / 500 "
          f"effective transitions  (slope {rep['dyn_eff_slope']:+.2f}/kstep)")
    print(f"    dyn_support (final) = {rep['dyn_support_final']:.0f} / 500 "
          f"transitions ever used")
    print(f"    dyn_excess_pc final = {rep['dyn_excess_pc_final']:.2f}  "
          f"(corroborating only; artifact-prone when sparse)")
    print(f"  VERDICT: {rep['verdict'].upper()}")
    print(f"    {rep['verdict_reason']}")
    return rep


if __name__ == '__main__':
    # tiny smoke run: the two demo regimes at equal LUT mutation rate.
    print("=== FUNCTIONAL (rich) ===")
    print_report(dyn_report(
        dict(food_inc=0.2, m_scale=1.0, tax=0.02, mu_lut=0.03,
             mu_egenome=0.0),
        n_steps=800, sample_every=50, N=80, seed=2))
    print("\n=== SILENT-CHURN (converged) ===")
    print_report(dyn_report(
        dict(food_inc=0.05, m_scale=1.0, tax=0.05, mu_lut=0.03,
             mu_egenome=0.0),
        n_steps=800, sample_every=50, N=80, seed=2))
