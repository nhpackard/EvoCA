"""
patch_transfer.py — §11 patch-transfer success assay.

Implements the protocol of Docs/EvoCA_research_directions.md §11 (and
the "evolutionary test" of Docs/egene_discussion.md): transplant a
genome-level patch of one evolved population into a host field and
measure whether it grows or shrinks.

The science is in the controls, not the transplant. This module
provides:

  * run_assay()        — one directional transplant + its scoring
  * self_transfer()    — the mandatory null: a population stamped back
                         into a field initialised from its own recipe
  * reciprocal_2x2()   — the recommended unit of measurement: the 2×2
                         sign structure of {A-in-B, B-in-A}, EACH
                         reported relative to its own self-transfer
                         control (§11 control 1 + 2)
  * size_series()      — the patch-size sweep (§11 control 3): probes
                         a propagule-pressure / critical-patch-size
                         threshold

Scoring (§11): the leading indicator is the boundary flux across the
patch perimeter (births landing on the perimeter minus deaths of
patch-resident cells); the alive-fraction trajectory relative to t=0 is
the lag indicator. "Success" is sustained positive boundary flux, not
merely "the patch did not shrink in the first K ticks".

Everything is headless and returns plain dicts/lists of floats so the
caller can serialise or aggregate without importing this module's
types.
"""

import numpy as np

from evoca_py import EvoCA


# ── low-level: one directional transplant ─────────────────────────────

def _alive_map(sim):
    return sim.get_alive().astype(bool)


def _make_sim(N, metaparams):
    """Fresh EvoCA at the given lattice size and metaparameters."""
    sim = EvoCA()
    sim.init(N, **metaparams)
    return sim


def _seed_from_recipe(sim, recipe):
    """Apply a state-recipe dict to sim via sim.state(**recipe)."""
    sim.state(**recipe)


def _patch_box(N, side):
    """Centred square of edge `side`; returns (r0, c0)."""
    side = int(side)
    r0 = (N - side) // 2
    return r0, r0


def _perimeter_mask(N, r0, c0, side):
    """Boolean (N,N) mask of the one-cell-thick ring just OUTSIDE the
    patch — the cells a patch organism reproduces into when it expands.
    Clipped to the lattice (no periodic wrap, matching stamp_patch)."""
    m = np.zeros((N, N), dtype=bool)
    r1, c1 = r0 + side, c0 + side
    if r0 - 1 >= 0:
        m[r0 - 1, max(c0 - 1, 0):min(c1 + 1, N)] = True
    if r1 < N:
        m[r1, max(c0 - 1, 0):min(c1 + 1, N)] = True
    if c0 - 1 >= 0:
        m[max(r0 - 1, 0):min(r1 + 1, N), c0 - 1] = True
    if c1 < N:
        m[max(r0 - 1, 0):min(r1 + 1, N), c1] = True
    return m


def run_assay(donor, host, N, side, ticks, metaparams,
              mode='overwrite', t0_extract_ticks=0):
    """Transplant a centred `side`×`side` patch of a `donor` population
    into a `host` population and score the result.

    donor, host : state-recipe dicts (passed to sim.state()).
    metaparams  : kwargs for sim.init() (food_inc, m_scale, tax, mu_*,
                  …) — shared by donor field, host field and the
                  resident dynamics.
    ticks       : number of CA steps to run after the transplant.
    mode        : 'overwrite' | 'invade' (forwarded to stamp_patch).
    t0_extract_ticks : let the donor field settle this many ticks
                  before the patch is extracted (0 = extract from the
                  freshly seeded recipe).

    Returns a dict with the alive-fraction trajectory (patch-local,
    relative to t=0) and the per-tick boundary flux. See module
    docstring for the scoring rationale.
    """
    # 1. Build & (optionally) settle the donor; extract the patch.
    dsim = _make_sim(N, metaparams)
    _seed_from_recipe(dsim, donor)
    for _ in range(int(t0_extract_ticks)):
        dsim.step()
    r0, c0 = _patch_box(N, side)
    patch = dsim.extract_patch(r0, c0, side)

    # 2. Build the host field, stamp the patch in.
    hsim = _make_sim(N, metaparams)
    _seed_from_recipe(hsim, host)
    n_written = hsim.stamp_patch(patch, r0, c0, mode=mode)

    rs, cs = slice(r0, r0 + side), slice(c0, c0 + side)
    perim = _perimeter_mask(N, r0, c0, side)
    perim_n = int(perim.sum())

    # t=0 reference, taken AFTER stamping (the realised propagule).
    a0 = _alive_map(hsim)
    f0 = float(a0[rs, cs].mean())
    n0 = int(a0[rs, cs].sum())

    alive_frac = [f0]                 # patch-local alive fraction
    alive_frac_rel = [0.0]            # relative to t=0
    boundary_flux = []                # births_on_perimeter − patch deaths
    cum_flux = []
    perim_births = []
    patch_deaths = []

    prev_patch_alive = a0[rs, cs].copy()
    running = 0
    for _ in range(int(ticks)):
        hsim.step()
        a = _alive_map(hsim)
        births = hsim.get_births().astype(bool)

        # Leading indicator: expansion across the perimeter ring minus
        # contraction (resident patch cells that died this tick).
        b_out = int((births & perim).sum())
        cur_patch_alive = a[rs, cs]
        died = int((prev_patch_alive & ~cur_patch_alive).sum())
        flux = b_out - died
        running += flux

        perim_births.append(b_out)
        patch_deaths.append(died)
        boundary_flux.append(flux)
        cum_flux.append(running)

        fr = float(cur_patch_alive.mean())
        alive_frac.append(fr)
        alive_frac_rel.append(fr - f0)
        prev_patch_alive = cur_patch_alive.copy()

    n_ticks = int(ticks)
    # "Success" = sustained positive boundary flux: the mean flux over
    # the back half of the run (transient excluded) is > 0.
    half = max(1, n_ticks // 2)
    sustained = float(np.mean(boundary_flux[-half:])) if boundary_flux else 0.0

    return {
        'N': N, 'side': int(side), 'ticks': n_ticks, 'mode': mode,
        't0_extract_ticks': int(t0_extract_ticks),
        'n_written': int(n_written),
        'perimeter_cells': perim_n,
        'patch_alive_t0': n0,
        'alive_frac_t0': f0,
        'alive_frac': alive_frac,                  # len ticks+1
        'alive_frac_rel': alive_frac_rel,          # len ticks+1
        'boundary_flux': boundary_flux,            # len ticks (leading)
        'cum_boundary_flux': cum_flux,             # len ticks
        'perimeter_births': perim_births,
        'patch_deaths': patch_deaths,
        'alive_frac_final': alive_frac[-1],
        'alive_frac_rel_final': alive_frac_rel[-1],
        'sustained_flux': sustained,               # back-half mean flux
        'success': bool(sustained > 0.0),
    }


# ── controls ──────────────────────────────────────────────────────────

def self_transfer(recipe, N, side, ticks, metaparams,
                   mode='overwrite', t0_extract_ticks=0):
    """§11 control 1 — the mandatory null.

    A patch of population P transplanted into a host field initialised
    from P's OWN recipe. Any growth/shrinkage here is a
    boundary/initialisation artefact, not competitive superiority.
    Every invasion result must be read relative to this. Same signature
    as run_assay with donor == host == recipe.
    """
    return run_assay(donor=recipe, host=recipe, N=N, side=side,
                     ticks=ticks, metaparams=metaparams, mode=mode,
                     t0_extract_ticks=t0_extract_ticks)


def reciprocal_2x2(recipe_A, recipe_B, N, side, ticks, metaparams,
                   mode='overwrite', t0_extract_ticks=0):
    """§11 control 2 — the recommended unit of measurement.

    Adaptive-dynamics (Geritz et al. 1998) makes the question precise:
    not "does A grow in B" but the 2×2 sign structure of
    {A-in-B, B-in-A}, each reported *relative to* its own self-transfer
    control. Returns a dict holding all four assays plus a coarse
    classification:

      both invade        → 'coexistence'
      exactly one invades → 'exclusion_A' / 'exclusion_B' (genuine
                            evolutionary success of the invader)
      neither invades     → 'mutual_resistance'

    "Invades" means the cross-transfer's sustained boundary flux
    exceeds its self-transfer control's (the artefact-corrected sign).
    """
    A_in_B = run_assay(recipe_A, recipe_B, N, side, ticks, metaparams,
                        mode, t0_extract_ticks)
    B_in_A = run_assay(recipe_B, recipe_A, N, side, ticks, metaparams,
                        mode, t0_extract_ticks)
    A_self = self_transfer(recipe_A, N, side, ticks, metaparams,
                            mode, t0_extract_ticks)
    B_self = self_transfer(recipe_B, N, side, ticks, metaparams,
                            mode, t0_extract_ticks)

    # Artefact-corrected invasion sign: cross flux minus the invader's
    # own self-transfer flux (the null for its boundary/init artefact).
    A_adv = A_in_B['sustained_flux'] - A_self['sustained_flux']
    B_adv = B_in_A['sustained_flux'] - B_self['sustained_flux']
    A_inv = A_adv > 0.0
    B_inv = B_adv > 0.0
    if A_inv and B_inv:
        cls = 'coexistence'
    elif A_inv and not B_inv:
        cls = 'exclusion_A'      # A is the evolutionary winner
    elif B_inv and not A_inv:
        cls = 'exclusion_B'      # B is the evolutionary winner
    else:
        cls = 'mutual_resistance'

    return {
        'A_in_B': A_in_B, 'B_in_A': B_in_A,
        'A_self': A_self, 'B_self': B_self,
        'A_advantage': A_adv, 'B_advantage': B_adv,
        'A_invades': bool(A_inv), 'B_invades': bool(B_inv),
        'classification': cls,
        'N': N, 'side': int(side), 'ticks': int(ticks), 'mode': mode,
    }


def size_series(donor, host, N, ticks, metaparams,
                sides=None, mode='overwrite', t0_extract_ticks=0,
                with_self_control=True):
    """§11 control 3 — the patch-size sweep.

    Runs the directional assay at patch side ∈ `sides` (default
    {N/6, N/4, N/3}, per §11). A critical patch size below which even a
    superior population fails to establish is itself the key measurable
    (a propagule-pressure threshold; connects to §5 correlation
    length). Each size's cross-transfer is paired with its self-transfer
    control so the threshold is read on the artefact-corrected
    advantage, not the raw fraction.

    Returns a list of per-size dicts, ascending in side.
    """
    if sides is None:
        sides = sorted({max(2, N // 6), max(2, N // 4), max(2, N // 3)})
    out = []
    for side in sides:
        cross = run_assay(donor, host, N, side, ticks, metaparams,
                          mode, t0_extract_ticks)
        rec = {'side': int(side), 'cross': cross}
        if with_self_control:
            ctrl = self_transfer(donor, N, side, ticks, metaparams,
                                  mode, t0_extract_ticks)
            rec['self'] = ctrl
            rec['advantage'] = (cross['sustained_flux']
                                - ctrl['sustained_flux'])
            rec['establishes'] = bool(rec['advantage'] > 0.0)
        out.append(rec)
    return out
