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
- Rule LUT (bit-packed, 32 bytes/cell) — maps `(v_x, n1, n2, n3)` to new state
- `c(x)` — fiducial configuration pattern (6-bit D4-symmetric genome for eating)
- `f(x)` — private food store (float)

**LUT indexing — per-ring counts**:
The LUT is indexed by `(v_x, n1, n2, n3)` where `nk` = count of active cells in distance-ring k:

| Ring | Distance | Cells | Max count |
|---|---|---|---|
| n1 | 1 | (±1,0),(0,±1) | 4 |
| n2 | √2 | (±1,±1) | 4 |
| n3 | 2 | (±2,0),(0,±2) | 4 |

Flat bit index: `v_x*125 + n1*25 + n2*5 + n3`
Total: 2×5×5×5 = **250 bits = 32 bytes** (bit-packed) per cell.

GoL is exactly encodable: it conditions on n1+n2 and ignores n3.
The fiducial pattern for eating still uses the full 5×5 neighbourhood.

**GoL initialization** (`make_gol_lut()`): sets new_state = 1 iff Moore count `n1+n2` == 3 (dead cell) or ∈ {2,3} (alive cell), regardless of n3. With mutation=0 (all cells keep the same LUT), the simulation runs exact Conway's Game of Life.

**Food dynamics** each time step:
1. `F(x) += food_inc` (uniform regeneration)
2. Tax: `f(x) -= tax`; if `f(x)` reaches 0, cell's LUT is zeroed (death)
3. Each cell eats: `M(x) = (m/25) · matches · F(x)` (proportional to available food)
4. Reproduction: when `f(x) >= food_repro`, copy genome to the Moore-neighbor with lowest `f(x')`; split food 50/50

**Mutation** (applied to child's genome during reproduction):
- `mu_lut`: per-bit flip probability for the 250-bit LUT. n_flips drawn from Poisson(mu_lut * 250).
- `mu_cgenom`: per-bit flip probability for the 6-bit cgenom. n_flips drawn from Poisson(mu_cgenom * 6).
- `restricted_mu` (toggle, default off): when enabled, LUT mutations are restricted
  to bit positions that were actually queried during the current CA step. A 250-bit
  mask `lut_active` is built during Phase 1; only the `n_active` set bits are
  eligible for mutation (Poisson rate becomes `mu_lut * n_active`). This ensures
  every LUT mutation is immediately phenotypic — no silent mutations that merely
  change the hash without affecting dynamics.

**Random initialization helpers**:
- `set_lut_random(n_init)`: each cell gets an independent random LUT.
  `n_init` controls ring conditioning: 1 = (v_x, n1) only (10 bits),
  2 = (v_x, n1, n2) (50 bits, GoL-level), 3 = all 250 bits.
- `set_cgenom_random()`: each cell gets a random cgenom in [0, 63].

**Global metaparameters**: `food_inc`, `m_scale`, `food_repro`, `gdiff`, `mu_lut`, `mu_cgenom`, `tax`, `restricted_mu`

**Fiducial pattern `c(x)`**: D4-symmetric 5×5 binary pattern. The 25 cells form 6 orbits under D4 (reflections about horizontal/vertical midlines and diagonals), requiring 6 independent bits. The orbit map:

```
4  5  2  5  4
5  3  1  3  5
2  1  0  1  2    (orbit index at each grid position)
5  3  1  3  5
4  5  2  5  4
```

**Colormode 3 (births)**: yellow = normal birth (genome copied exactly),
magenta = mutant birth (LUT or cgenom mutated), dim grey = alive (no birth),
black = dead.

**Activity tracking**: Cumulative presence of each distinct LUT genome.
Each genome is identified by its FNV-1a hash (cached in `lut_hash_cache[N*N]`).
An open-addressing hash table (`act_keys`/`act_vals`) maps hash → `{activity, pop_count, color}`.
`evoca_activity_update()` clears pop_counts, scans alive cells, increments activity.
`evoca_activity_render_col(col, height)` renders one column of the scrolling strip chart.
Alive genomes in full color; extinct genomes dimmed (RGB × 0.15).

**Activity Y-axis — saturation formula** (from genelife):
`y = (H-1) - (H-1) * act / (act + ymax)`, where `act` is the cumulative activity
count and `ymax` (default 2000, tunable via `act_ymax` slider) controls the vertical
scale.  This is a hyperbolic saturation curve: `act=0` maps to `y=H-1` (bottom),
`act=ymax` maps to `y=(H-1)/2` (mid-chart), and `act→∞` approaches `y=0` (top).
New genomes start at the bottom and rise as they accumulate presence; the curve
compresses high-activity genomes toward the top without clipping, giving a natural
logarithmic-like spread.  Lowering `ymax` makes waves rise faster; raising it
spreads out low-activity genomes.

**Cgenom activity tracking**: Mirrors LUT activity for the 6-bit cgenom (fiducial
eating pattern). Since there are only 2^6 = 64 possible cgenoms, uses fixed-size
arrays (`cg_act[64]`, `cg_pop[64]`, `cg_color[64]`) instead of a hash table.
Wild-type cgenom is colored white; mutants get FNV-1a hash colors.
`evoca_cg_activity_update()` / `evoca_cg_activity_render_col()` parallel the LUT
activity functions. Separate `cg_act_ymax` slider. Enable with
`probes={'cg_activity': True}`.

**Reproduction age histogram**: Tracks the distribution of time between successive
reproduction events (or birth-to-first-reproduction). A per-cell timestamp
`last_event_step[N*N]` records the step of each cell's most recent birth or
reproduction. At each reproduction, `age = step - last_event_step[parent]` is
binned into `repro_age_hist[1024]`. A configurable `repro_age_t0` (default 0)
skips transient: only events where both the current step and the parent's last
event are >= t0 are counted. `reset_repro_age_hist()` clears the histogram.

**Performance** (CA step only, no food/repro, GoL LUT):
N=256 → ~150 fps · N=512 → ~30 fps
