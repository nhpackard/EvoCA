# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EvoCA is an evolutionary cellular automata simulation. It extends the GeneLife model (see `Docs/genelife-1.pdf` and `Docs/genelife-2.pdf`) with continuous food resources and evolved consumption patterns. Every cell is alive and carries a genome that governs its local CA rule and eating behavior.

## Code Layout

- `C/` — C core: CA lattice, genome data structures, and CA update functions (optimized for speed)
- `python/` — Python layer: `evoca_py.py` (ctypes wrapper), `display.py` (SDL2 window + widgets)
- `Docs/` — Reference papers on the GeneLife model
- `EvoCA.md` — Full model specification (read this for implementation details)
- `evoca_test.ipynb` — Test notebook: GoL verification, benchmarks, SDL2 display launch

## Build & Run

```bash
make                        # builds C/libevoca.dylib (macOS) or C/libevoca.so (Linux)
gcc -O2 -Wall -fPIC -dynamiclib -o C/libevoca.dylib C/evoca.c   # manual macOS build
jupyter notebook evoca_test.ipynb
```

The Python wrapper (`python/evoca_py.py`) loads the shared library via `ctypes`.
The SDL2 display (`python/display.py`) follows the genelife pattern: C fills an int32 ARGB pixel buffer via `evoca_colorize()`, Python reshapes it via numpy and copies to the SDL surface via `pixels2d()`, and software-rendered colored blocks serve as metaparam widgets.

**SDL2 controls** (interactive window):
`SPACE` pause/unpause · `C` cycle color mode · `S` single step (paused) · `Q`/`Esc` quit
Mouse click widget row: left `-` / right `+` to adjust food_inc / m_scale / food_repro.

## Model Architecture

**CA lattice**: Binary 2D grid (periodic boundaries). Every cell carries:
- `v(x)` — binary cell state (0 or 1)
- Rule LUT (bit-packed, 1407 bytes/cell) — maps `(v_x, n1, n2, n3, n4, n5)` to new state
- `c(x)` — fiducial configuration pattern (6-bit D4-symmetric genome for eating)
- `f(x)` — private food store (float)

**LUT indexing — per-ring counts**:
The LUT is indexed by `(v_x, n1, n2, n3, n4, n5)` where `nk` = count of active cells in distance-ring k:

| Ring | Distance | Cells | Max count |
|---|---|---|---|
| n1 | 1 | (±1,0),(0,±1) | 4 |
| n2 | √2 | (±1,±1) | 4 |
| n3 | 2 | (±2,0),(0,±2) | 4 |
| n4 | √5 | (±2,±1),(±1,±2) | 8 |
| n5 | 2√2 | (±2,±2) | 4 |

Flat bit index: `v_x*5625 + n1*1125 + n2*225 + n3*45 + n4*5 + n5`
Total: 2×5×5×5×9×5 = **11 250 bits = 1 407 bytes** (bit-packed) per cell.

**Why per-ring, not weighted sum**: An earlier design considered a Euclidean-norm weighted sum S = Σ v·‖δ‖ over the 5×5 neighbourhood. Because dist(n3) = 2 × dist(n1) and dist(n5) = 2 × dist(n2), these rings' contributions are conflated — a rule that depends on n1 alone (e.g. GoL, which uses n1+n2) cannot be distinguished from one that trades n1 for n3. Separate per-ring counts remove this ambiguity and make GoL exactly encodable.

**GoL initialization** (`make_gol_lut()`): sets new_state = 1 iff Moore count `n1+n2` == 3 (dead cell) or ∈ {2,3} (alive cell), regardless of n3/n4/n5. With mutation=0 (all cells keep the same LUT), the simulation runs exact Conway's Game of Life.

**Food dynamics** each time step:
1. `F(x) += food_inc` (uniform regeneration)
2. Each cell eats: mouthful `M(x) = (m/25) · matches(C(x), c(x))` transferred from `F(x)` to `f(x)`
3. Reproduction: when `f(x) >= food_repro`, copy genome to the Moore-neighbor with lowest `f(x')`; split food 50/50

**Mutation** (applied to child's genome during reproduction):
- `mu_lut`: per-bit flip probability for the 11250-bit LUT. n_flips drawn from Poisson(mu_lut * 11250).
- `mu_cgenom`: per-bit flip probability for the 6-bit cgenom. n_flips drawn from Poisson(mu_cgenom * 6).

**Global metaparameters**: `food_inc`, `m_scale`, `food_repro`, `gdiff`, `mu_lut`, `mu_cgenom`

**Fiducial pattern `c(x)`**: D4-symmetric 5×5 binary pattern. The 25 cells form 6 orbits under D4 (reflections about horizontal/vertical midlines and diagonals), requiring 6 independent bits. The orbit map:

```
4  5  2  5  4
5  3  1  3  5
2  1  0  1  2    (orbit index at each grid position)
5  3  1  3  5
4  5  2  5  4
```

**Performance** (CA step only, no food/repro, GoL LUT):
N=256 → ~150 fps · N=512 → ~30 fps
