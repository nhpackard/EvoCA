"""Shared harness for the N_chem=1 chem-branch parity gate (milestone 1).

The `chem` branch introduces vector-valued resource ("chem") machinery.
Milestone 1's contract is that at ``N_chem == 1`` the model reproduces the
pre-change single-food dynamics **bit for bit**. This module builds a fixed
deterministic initial state, steps it, and captures a full-state snapshot
(private food ``f``, environmental food ``F``, ``alive`` mask, and per-cell
LUT hashes) at several checkpoints.

Both ``tests/test_chem_parity.py`` (the pytest gate) and ``chem-test.ipynb``
(the demonstration notebook) call into here, so the golden reference and the
live run are produced by identical code. The golden ``.npz`` was captured from
the pre-change build via ``regen_golden()`` (see ``__main__`` block) and is
committed under ``tests/data/``.

Determinism notes:
  * The LUT init is numpy-seeded (``set_lut_random(seed=...)``), independent of
    the C PRNG.
  * ``set_egenome_random`` / ``set_alive_fraction`` consume the C PRNG, so the
    C PRNG is seeded (``SEED_INIT``) before them for a reproducible init.
  * The C PRNG is re-seeded (``SEED_STEP``) immediately before the step loop so
    the stepping stream is fixed regardless of how many draws init consumed.
  * The N_chem=1 code path must add **zero** PRNG draws relative to the
    pre-change hot path, or the stream diverges and parity fails — this is the
    strongest guard the snapshot provides.
"""
import os
import sys

import numpy as np

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if os.path.join(_ROOT, 'python') not in sys.path:
    sys.path.insert(0, os.path.join(_ROOT, 'python'))

GOLDEN_PATH = os.path.join(os.path.dirname(__file__), 'data',
                           'chem_parity_golden.npz')

# Fixed initial-condition / stepping parameters. Chosen so eating, tax death,
# reproduction, and LUT+egene mutation all fire (a config that exercises every
# phase the chem storage touches).
N          = 48
SEED_INIT  = 20260720   # C PRNG seed for the stochastic init helpers
SEED_LUT   = 1234       # numpy seed for set_lut_random
SEED_STEP  = 4242       # C PRNG seed for the step loop
N_STEPS    = 200
CHECKPOINTS = (50, 100, 150, 200)

# food_cap=0 pins the legacy (pre-A9-fix) uncapped mouthful rule: the golden
# was captured before the 2026-10-15 food_cap fix and m_scale=3 > 1, so the
# default capped rule would (correctly) diverge. This gate tests the chem
# storage path, not the eating rule, so it stays on the rule it was recorded
# with. See tests/test_food_cap.py for the cap itself.
CFG = dict(
    food_inc=0.02, m_scale=3.0, gdiff=0.05, tax=0.02,
    mu_lut=0.01, mu_egene=0.005, mu_egenome=0.01, p_dup_egene=1.0,
    food_cap=0,
)


def build_initial_state(sim):
    """Construct the fixed deterministic initial state on `sim`."""
    sim.init(N=N, **CFG)
    sim.set_seed(SEED_INIT)
    sim.set_lut_random(n_init=3, density=0.5, seed=SEED_LUT)
    sim.set_egenome_random()
    sim.set_alive_fraction(0.6)
    sim.set_f_all(0.3)
    sim.set_F_all(0.5)
    sim.set_seed(SEED_STEP)


def run_and_capture(sim):
    """Step the sim and return a dict of stacked snapshots at CHECKPOINTS.

    Keys: 'f', 'F', 'alive', 'lut_hash' — each shape (len(CHECKPOINTS), N, N).
    """
    build_initial_state(sim)
    fs, Fs, alives, hashes = [], [], [], []
    nxt = list(CHECKPOINTS)
    for step in range(1, N_STEPS + 1):
        sim.step()
        if nxt and step == nxt[0]:
            nxt.pop(0)
            fs.append(sim.get_f())
            Fs.append(sim.get_F())
            alives.append(sim.get_alive())
            hashes.append(sim.get_lut_hashes())
    return {
        'f':        np.stack(fs),
        'F':        np.stack(Fs),
        'alive':    np.stack(alives),
        'lut_hash': np.stack(hashes),
    }


def compare(golden, got):
    """Return (ok, report_lines). Requires bit-for-bit equality on every field.

    Floats (`f`, `F`) are compared with ``array_equal`` (exact bit equality,
    NaN-free by construction), not ``allclose`` — milestone 1 is a *bit-for-bit*
    gate, not an approximate one.
    """
    ok = True
    lines = []
    for key in ('f', 'F', 'alive', 'lut_hash'):
        g, x = golden[key], got[key]
        if g.shape != x.shape:
            ok = False
            lines.append(f"  {key:9s} SHAPE MISMATCH {g.shape} vs {x.shape}")
            continue
        eq = np.array_equal(g, x)
        if not eq:
            ok = False
            n_diff = int(np.count_nonzero(g != x))
            first = np.argwhere(g != x)
            where = tuple(int(v) for v in first[0]) if len(first) else None
            lines.append(f"  {key:9s} DIFFERS  ({n_diff} cells, first @ {where})")
        else:
            lines.append(f"  {key:9s} OK       (identical, {g.size} values)")
    return ok, lines


def load_golden():
    data = np.load(GOLDEN_PATH)
    return {k: data[k] for k in ('f', 'F', 'alive', 'lut_hash')}


def regen_golden():
    """Capture the golden snapshot from the CURRENTLY-BUILT libevoca.

    Run this ONLY against the pre-change (single-food) build to (re)mint the
    committed reference. Running it against the chem build would make the gate
    vacuously pass.
    """
    from evoca_py import EvoCA
    sim = EvoCA()
    snap = run_and_capture(sim)
    os.makedirs(os.path.dirname(GOLDEN_PATH), exist_ok=True)
    np.savez_compressed(GOLDEN_PATH, **snap)
    return snap


if __name__ == '__main__':
    snap = regen_golden()
    print(f"Wrote golden snapshot -> {GOLDEN_PATH}")
    for k, v in snap.items():
        print(f"  {k:9s} shape={v.shape} dtype={v.dtype}")
