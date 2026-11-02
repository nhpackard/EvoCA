"""ChemCA port of N's EvoCA food/waste spiral recipe: the parameter mapping
and a runner that records the same readouts as s1_evoca_characterize.py.

Mapping (see notes.md for the reasoning and what is approximate):
  EvoCA food_inc, F <= 1 cap  -> open-mode inc_F, cap_F = 1, sink_F = 0 (exact;
                                 before the ChemCA cap: a proportional sink_F)
  gdiff, wdiff                -> diff_F, diff_W (same box-blur rule)
  waste_sink (constant)       -> proportional sink_W
  waste_disp (one-way)        -> u_W (two-way transport, = EvoCA waste_exchange=1)
  waste_tox                   -> tox_W
  tax (effective, per cell)   -> m_e (c_g = 0)
  eating (pattern-matched)    -> u_F (transport of F into the cell)
  f (private food = energy)   -> e, with fast F -> W conversion (k_idx 63)
  f >= 1 reproduction         -> e >= 1 (same split and eviction rule)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', '..', 'python'))
from chemca_py import ChemCA  # noqa: E402

# EvoCA recipe 2026-11-02_spiral_wave_ish (metaparams_final), and its
# effective tax: tax + tax_per_egene * 25 (all egene positions specified).
EVOCA = dict(food_inc=0.012, gdiff=0.08, wdiff=0.01, waste_sink=0.019, waste_disp=0.09,
             waste_tox=0.792, tax_eff=0.0301 + 0.0008 * 25)


def port_params(u_F=0.5, sink_F=None, sink_W=None, W_ref=0.1, tox_W=None, u_W=None,
                m_e=None, k_idx=63, cap_F=0.0):
    """(init params, state params) for ChemCA fw2 in open mode. Unset
    arguments take the direct numeric port of the EvoCA recipe:
    sink_F = food_inc (F saturates at inc/sink = 1) unless cap_F is set
    (cap_F = 1, sink_F = 0 is EvoCA's food rule exactly),
    sink_W = waste_sink / W_ref (the constant sink at a typical W)."""
    E = EVOCA
    p = dict(preset='fw2', open_mode=1, light=0.0, mu_rgene=0.0, mu_ractive=0.0,
             m_e=E['tax_eff'] if m_e is None else m_e, c_g=0.0,
             G_F=1.0, G_W=0.0,
             diff_F=E['gdiff'], diff_W=E['wdiff'],
             u_F=u_F, u_W=E['waste_disp'] if u_W is None else u_W,
             tox_F=0.0, tox_W=E['waste_tox'] if tox_W is None else tox_W,
             inc_F=E['food_inc'], cap_F=cap_F,
             sink_F=(0.0 if cap_F else E['food_inc']) if sink_F is None else sink_F,
             inc_W=0.0, sink_W=E['waste_sink'] / W_ref if sink_W is None else sink_W)
    st = dict(alive='halfplane', mbar=1.0, mass_split=(1.0, 0.0), store_frac=0.0,
              e_init=0.5, k_idx0=k_idx)
    return p, st


def run(p, st, N=256, steps=6000, seed=0, kymo=2000, every=10, frames_at=None):
    """Same readouts as s1: dict of time series, kymographs, frames."""
    frames_at = set(frames_at or range(500, steps + 1, 500))
    sim = ChemCA()
    sim.init(N, **p)
    sim.state(**st, seed=seed)
    ts = {k: [] for k in ('t', 'alive', 'tax_pc', 'intake_pc', 'eat_frac', 'F', 'W',
                          'F_alive', 'w_alive', 'tox', 'births', 'deaths')}
    kym = {k: [] for k in ('alive', 'F', 'W')}
    frames = {}
    acc = dict(tox=0, births=0, deaths=0)
    for t in range(1, steps + 1):
        sim.step()
        ev = sim.events_last()
        acc['tox'] += ev['tox_deaths']; acc['births'] += ev['births']; acc['deaths'] += ev['deaths']
        if not (t % every == 0 or t > steps - kymo or t in frames_at):
            continue
        A = sim.get_alive().astype(bool)
        F, W = sim.get_E('F'), sim.get_E('W')
        if t > steps - kymo:
            kym['alive'].append(A[N // 2].astype(np.uint8))
            kym['F'].append(F[N // 2].copy()); kym['W'].append(W[N // 2].copy())
        if t in frames_at:
            frames[t] = dict(rgb=sim.rgb(1), F=F.copy(), W=W.copy(), alive=A)
        if t % every == 0:
            b = sim.energy_budget(cum=False)
            n = max(int(A.sum()), 1)
            y = b['yield']                       # F converted = F taken up (fast k)
            ts['t'].append(t); ts['alive'].append(A.mean())
            ts['tax_pc'].append(b['cost'] / n); ts['intake_pc'].append(y / n)
            ts['eat_frac'].append(y / max(float(F[A].sum()) + y, 1e-12))
            ts['F'].append(F.mean()); ts['W'].append(W.mean())
            ts['F_alive'].append(F[A].mean() if A.any() else np.nan)
            ts['w_alive'].append(sim.get_store('W')[A].mean() if A.any() else np.nan)
            for k in acc:
                ts[k].append(acc[k]); acc[k] = 0
        if not A.any() and t % every == 0:
            break
    sim.free()
    return dict(**{f'ts_{k}': np.array(v) for k, v in ts.items()},
                **{f'kym_{k}': np.array(v) for k, v in kym.items()},
                frame_t=np.array(sorted(frames)),
                **{f'fr_{k}': np.array([frames[t][k] for t in sorted(frames)])
                   for k in ('rgb', 'F', 'W', 'alive')})
