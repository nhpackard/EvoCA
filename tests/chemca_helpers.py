"""Constructed ChemCA fixtures shared by tests/test_chemca_*.py.

isolated(): a small world with every flux switched off (no light, no
diffusion, no transport, no costs, no toxicity) so one mechanism at a time
can be switched back on and checked against a hand calculation.
"""
import numpy as np

from chemca_py import ChemCA, gene

QUIET = dict(light=0.0, m_e=0.0, c_g=0.0, mu_rgene=0.0, mu_ractive=0.0)


def quiet_compounds(preset):
    from chemca_py import PRESETS
    out = {}
    for c in PRESETS[preset]['compounds']:
        out.update({f'diff_{c}': 0.0, f'u_{c}': 0.0, f'tox_{c}': 0.0})
    return out


def isolated(preset='gs3', N=5, cells=((2, 2),), genes=None, active=0b1,
             stores=None, e=0.0, E=None, **params):
    """A ChemCA with only `cells` alive, given genes/stores/e, all fluxes off
    except those re-enabled through params. stores: {compound: value or
    N x N array}; E: {compound: value or array} (default 0)."""
    sim = ChemCA()
    sim.init(N, preset=preset, **{**QUIET, **quiet_compounds(preset), **params})
    sim.set_seed(1)
    a = np.zeros((N, N), np.uint8)
    for r, c in cells:
        a[r, c] = 1
    sim.set_alive(a)
    for c in sim.compounds:
        sim.set_E(c, (E or {}).get(c, 0.0))
        sim.set_store(c, (stores or {}).get(c, 0.0))
    sim.set_e(e)
    sim.set_genes_all(genes or [gene(sim.compounds[0], sim.compounds[1])], active)
    return sim


def cell(arr, rc=(2, 2)):
    return float(arr[rc])
