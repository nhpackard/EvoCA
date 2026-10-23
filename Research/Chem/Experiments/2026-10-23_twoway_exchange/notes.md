# Experiment: two-way (passive membrane) waste exchange (2026-10-23)

**Status.** An experiment outside the model. `main` still runs the one-way rule. Whether to adopt the two-way rule is CHEM-D4 (N).

## Why

The food/waste model makes no spirals (CHEM-N5, `Chem.ipynb` "No spirals yet"). Its refractory zone is soft. In the one-way rule a cell excretes `d = disp * a * max(0, 1 - E)`: waste only leaves the cell, and stops leaving once the site holds `E >= 1`. A cell on waste-loaded ground therefore lives about `(tox - a) / tax` steps, and a fast-reproducing lineage can persist there. Waste ground never empties of life, so no clean food lakes form for a front to sweep across.

## The rule tested

Passive diffusion across the membrane, in both directions:

- `d = disp * (a - E)`.
  - `d > 0`: the cell excretes.
  - `d < 0`: the cell takes waste up from its site.
- `d` is capped by what the source holds.
- The steady private load is `E + tax/disp`, so a cell dies wherever the site holds about `E > tox - tax/disp` (0.93 at the defaults). That ground is lethal within a few steps, which makes the refractory zone hard. No extra rule is added.
- Food is unchanged (active uptake through the egenes).
- Mass is conserved: the exchange only moves waste between site and cell. The per-step waste identity error is 8.6e-5 for both rules (float rounding; `run2.py: mass_check`).

## Reproduce

From the repo root:

```
patch -o Research/Chem/Experiments/2026-10-23_twoway_exchange/evoca_patched.c C/evoca.c Research/Chem/Experiments/2026-10-23_twoway_exchange/twoway.patch
gcc -O2 -Wall -fPIC -dynamiclib -I C -o Research/Chem/Experiments/2026-10-23_twoway_exchange/libevoca_2way.dylib Research/Chem/Experiments/2026-10-23_twoway_exchange/evoca_patched.c
cd Research/Chem/Experiments/2026-10-23_twoway_exchange && python3 run2.py && python3 kymo.py && python3 phase.py
```

The patch applies to `C/evoca.c` at commit `0e2ee2a`. Both build products are gitignored.

All runs use recipe `2026-10-23_chem_fw_spiralpoint` (standard rule, mutation on, halfplane start) and seed 0. They are example runs, not campaign results.

## Results

**1. Sweep** (`run2.py`, `grid_sink_wdiff.png`). N=256, t=3000; statistics over steps 1500-3000.

| `waste_sink` | `wdiff` 0.01: alive ± std, toxic share, W | `wdiff` 0.08 |
|---|---|---|
| 0.010 | 0.24 ± 0.012, 97%, 1.34 | 0.22 ± 0.013, 97%, 1.30 |
| 0.013 | 0.30 ± 0.016, 91%, 1.17 | 0.29 ± 0.018, 84%, 1.12 |
| 0.015 | 0.32 ± 0.018, 49%, 0.78 | **0.35 ± 0.037**, 28%, 0.60 |
| 0.017 | 0.33 ± 0.019, 31%, 0.54 | 0.38 ± 0.004, 0%, 0.27 |
| 0.019 | **0.35 ± 0.028**, 33%, 0.45 | 0.38 ± 0.004, 0%, 0.21 |

- **The structured band is wider.** At `wdiff` 0.01 the oscillating three-state mosaic spans `waste_sink` 0.013-0.019. Under the one-way rule it is one or two grid steps (N5: 0.015-0.016 at `wdiff` 0.01).
- **W stays bounded** (at most 1.3, against up to 19 under the one-way rule). Cells take waste back up instead of leaving it to pile up.
- Patches are larger and cleaner than under the one-way rule.

**2. Space-time plots** (`kymo.py`, `kymo.png`). N=512. The middle row is plotted every step over t = 3000-3600.

- **One-way** (sink 0.015, `wdiff` 0.01): short, uncorrelated vertical streaks. Sites cycle locally and nothing travels.
- **Two-way** (sink 0.017 / `wdiff` 0.01, 0.015 / 0.08, 0.019 / 0.01):
  - a coherent oscillation along the whole row, period about 100 steps;
  - **tilted bands**: phase traveling along the row;
  - **chevrons**: sources emitting waves to both sides.

  These are traveling waves. Chevrons arise from pacemakers (target patterns) and from spiral cores alike, so they do not identify spirals by themselves.

**3. Phase maps and core persistence** (`phase.py`, `phase.png`). N=512, 8x8 blocks, t = 3000-4500. A core counts as persistent if a same-sign defect is within 1 block 60 steps later; chance is the same test against the earlier frame shifted by half the lattice.

| Rule | Defects | Persistent | Chance | Ratio |
|---|---|---|---|---|
| one-way 0.015 / 0.01 | 668 | 421 | 392 | 1.08 |
| two-way 0.017 / 0.01 | 354 | 195 | 124 | 1.57 |
| two-way 0.015 / 0.08 | 316 | 189 | 100 | 1.89 |

Two-way gives larger coherent phase domains and defects that persist above chance. But the N5 scan (16x16 blocks) shows the ratio reaching 1.69 in one-way runs where nothing propagates, so a ratio above 1 is not a spiral signature, and no clean pinwheel is visible. The colonies' fine speckle (alive about 0.33, scattered) still breaks up the phase field at block scale.

## Reading

- Two-way exchange turns the one-way model's locally cycling mosaic into a **coherent oscillatory medium with traveling phase waves and wave sources**.
- That is the precondition for spirals, and the step the one-way rule lacked.
- **Spirals themselves are not established.** The next checks:
  - several seeds;
  - lower `waste_tox` under the two-way rule, for a harder zone;
  - N=1024, for wavelength room;
  - a films-based look at defect rotation.
