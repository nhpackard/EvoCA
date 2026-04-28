# EvoCA Model Reference

Evolutionary Cellular Automata: a spatially inhomogeneous binary 2D CA where every cell is regarded as an organism, carrying a genome governing its local rule.  There is a resource field (food) modeled as a spatial field over a lattice the same size as the CA lattice.  Each organism at each lattice site can interact with the food field by eating a mouthful of food, transferring a fraction of the food from the food field to the organism's private stash of food.  The eating algorithm is controlled by a local pattern match using an additional piece of genetic data carried by each organism.  Every time step, the organism is taxed; its private food stash is diminished by a fixed amount.  The survival of an organism is a competition between the tax and the ability of the organism to evolve a rule genome and an eating genome to get more food than is lost by the tax.  If an organism's food stash reaches 1.0, it reproduces with genetic mutation; the offspring replaces the neighbor with the least food.  If an organism's food stash reaches 0, it dies: its  LUT is set to all zeroes. It no longer eats, its CA state is zero.  Dead cells can only come alive through a reproductive event.

## Table of Contents

- [Grid State](#grid-state)
- [Global Metaparameters](#global-metaparameters)
- [Restricted Mutation](#restricted-mutation)
- [LUT Indexing (Per-Ring Counts)](#lut-indexing-per-ring-counts)
- [Fiducial Pattern c(x) and egenome](#fiducial-pattern-cx-and-egenome)
- [Time Step Phases](#time-step-phases)
- [Visualization (Color Modes)](#visualization-color-modes)
- [Python API](#python-api)
- [Probes, Activity Tracking, and Diagnostics](#probes-activity-tracking-and-diagnostics)
- [Controls and Display](#controls-and-display)
- [LUT Helper Functions](#lut-helper-functions)
- [Build](#build)
- [Architecture Notes](#architecture-notes)

---

## Grid State

The simulation is an N x N lattice with periodic (toroidal) boundaries.
Every cell at position **x** carries:

| Field       | Type    | Description                                     |
|-------------|---------|-------------------------------------------------|
| `alive(x)`  | uint8   | 1 = alive organism, 0 = dead (empty) slot        |
| `v(x)`      | uint8   | Binary CA state (0 or 1) — not life/death        |
| `lut(x)`    | 32 B    | Bit-packed rule LUT (250 bits)                   |
| `egenome(x)` | uint8   | 6-bit fiducial configuration genome (D4-symmetric) |
| `f(x)`      | float   | Private food store, clamped to [0, 1]            |
| `F(x)`      | float   | Environmental food at this location, clamped to [0, 1] |

**Alive vs. CA state**: `alive(x)` tracks whether an organism occupies the
cell. `v(x)` is the binary CA state used by the LUT rule — it is *not*
life/death. Dead cells (`alive=0`) have zeroed LUT, `v=0`, `f=0`, and do
not eat or reproduce. Reproduction is the only way a dead cell becomes alive.

---

## Global Metaparameters

All metaparameters can be set at init or adjusted at runtime via sliders.

| Name            | Type  | Default | Range (slider) | Description                                          |
|-----------------|-------|---------|----------------|------------------------------------------------------|
| `food_inc`      | float | 0.0     | [0, 0.5]       | Environmental food added per cell per step            |
| `m_scale`       | float | 1.0     | [0, 10]        | Mouthful scale factor for eating                      |
| `gdiff`         | int   | 0       | [0, 10]        | Food diffusion passes (3x3 box blur) per step         |
| `mu_lut`        | float | 0.0     | [0, 0.001]     | Per-bit LUT mutation probability on reproduction      |
| `mu_egene`      | float | 0.0     | [0, 0.05]      | Per-egene-bit flip probability on reproduction (across all 8×6 = 48 bits per cell, including inactive slots — those drift as pseudogenes) |
| `mu_egenome`    | float | 0.0     | [0, 0.05]      | Per-active-bit flip probability on reproduction (across the 8 presence bits). Flips that would take Negene to 0 are rejected |
| `p_dup_on_activate` | float | 1.0 | [0, 1]          | Probability that a 0→1 active-bit flip overwrites the new slot's egene byte with a copy from a random currently-active slot (gene duplication). Copy happens before the egene-bit flip pass, so the new copy still receives fresh independent point mutations |
| `tax`           | float | 0.0     | [0, 0.1]       | Constant private-food decrement per step; death if depleted |
| `tax_per_egene` | float | 0.0     | [0, 0.01]      | Additional decrement per active egene per step. Bounds Negene against the unbounded "more is better" pressure of max-match eating |
| `tax_lut`       | float | 0.0     | [0, 0.001]     | Additional decrement per LUT '1' bit per step. Penalises rule complexity |
| `restricted_mu` | int   | 0       | checkbox       | If 1, restrict LUT mutations to dynamically active bits |

---

## Restricted Mutation

When `restricted_mu=1`, LUT mutations are restricted to bit positions that
were actually queried during the current CA step.

During Phase 1, a 250-bit mask `lut_active` records which LUT indices were
looked up.  Before Phase 4, the set bits are extracted into an `active_bits`
array (typically 20-50 of 250 bits are active in typical dynamics).

During reproduction, instead of `Poisson(mu_lut * 250)` flips at random
positions, the mutation draws `Poisson(mu_lut * n_active)` flips from only
the active positions.  This ensures every mutation is immediately phenotypic
-- no silent mutations that change the genome hash without affecting dynamics.

When `restricted_mu=0` (default), all 250 bits are eligible (original behavior).

---

## LUT Indexing (Per-Ring Counts)

The LUT maps the local neighborhood configuration to a new cell state.
It is indexed by `(v_x, n1, n2, n3)` where each `nk` is the count of
active (v=1) cells in distance-ring k:

| Ring | Distance | Offsets                            | Cells | Max count |
|------|----------|------------------------------------|-------|-----------|
| n1   | 1        | (+-1,0), (0,+-1)                   | 4     | 4         |
| n2   | sqrt(2)  | (+-1,+-1)                          | 4     | 4         |
| n3   | 2        | (+-2,0), (0,+-2)                   | 4     | 4         |

**Flat bit index**: `v_x*125 + n1*25 + n2*5 + n3`

**Total**: 2 x 5 x 5 x 5 = **250 bits = 32 bytes** per cell.

**Ring map** (indexed by `[di+2][dj+2]`, -1 = centre or outside LUT):

```
-1  -1   2  -1  -1
-1   1   0   1  -1
 2   0  -1   0   2
-1   1   0   1  -1
-1  -1   2  -1  -1
```

GoL (B3/S23) is exactly encodable: it conditions on n1+n2 and ignores n3.

**GoL initialization**: The function `make_gol_lut()` in `python/evoca_py.py`
constructs a LUT implementing exact Conway's Game of Life (B3/S23).  For a
dead cell (v_x=0), the new state is 1 iff Moore count n1+n2 == 3.  For an
alive cell (v_x=1), the new state is 1 iff n1+n2 in {2, 3}.  The n3 value
is a don't-care: all five n3 entries for each (v_x, n1, n2) combination are
set identically.  With mutation=0 (all cells share this LUT), the simulation
runs exact GoL.

Note: the fiducial pattern for eating still uses the full 5x5
neighbourhood (all 6 D4 orbits).  Only the CA rule LUT is restricted
to 3 rings.

---

## Fiducial Pattern c(x) and egenome

Each cell carries a small **list of egenes** — up to `NEGENOME_MAX = 8`
slots — and an **active mask** marking which slots are live.

- `egenes[NEGENOME_MAX][6 bits]` per cell — each slot is a 6-bit
  D4-symmetric fiducial pattern (same encoding as the legacy single
  egenome).
- `active[8]` per cell — 1 bit per slot. `Negene = popcount(active)`.
  An alive cell always has Negene ≥ 1.

Inactive slots' egene bytes are still mutated at the per-egene-bit
rate, accumulating drift like pseudogenes; if an inactive slot is
later activated by a presence-bit flip, its drifted bits become live
genome (and may be overwritten by a duplicated copy of an existing
active egene first — see `p_dup_on_activate`).

The list of slots is **unordered** for purposes of identity: the
species hash sorts active egene bytes lex before folding them into the
genome hash. Two cells with the same LUT and the same set of active
egenes are the same species, regardless of slot positions.

When the section below talks about "the egenome" it refers to a single
slot. The next subsection covers how the cell uses *all* its slots.

### D4 Orbits

The 25 positions of the 5x5 neighborhood fall into 6 orbits under the
D4 symmetry group (reflections about horizontal, vertical, and both
diagonal axes).  Each orbit is controlled by one bit of `egenome`:

```
 4   5   2   5   4        bit 0: centre               (1 cell)
 5   3   1   3   5        bit 1: axis neighbors        (4 cells)
 2   1   0   1   2        bit 2: axis distance-2       (4 cells)
 5   3   1   3   5        bit 3: diagonal neighbors    (4 cells)
 4   5   2   5   4        bit 4: corners               (4 cells)
                           bit 5: knight-move positions (8 cells)
```

The number at each position is the **orbit index** = the egenome bit
that controls it.  The pattern value at position (i,j) is:

    pattern[i][j] = (egenome >> orbit_map[i][j]) & 1

### egenome Examples

| egenome     | Binary   | Pattern description                      |
|------------|----------|------------------------------------------|
| `0b000000` | 000000   | All zeros (matches dead cells everywhere) |
| `0b111111` | 111111   | All ones (matches alive cells everywhere) |
| `0b000010` | 000010   | Only axis neighbors = 1 (plus sign)       |
| `0b001010` | 001010   | Axis + diagonal neighbors (3x3 filled)    |
| `0b000001` | 000001   | Only centre = 1                           |
| `0b010100` | 010100   | Corners + axis-2 (ring pattern)           |

To visualize any egenome value:

```python
from python.evoca_py import egenome_to_pattern
pat = egenome_to_pattern(0b001010)
print(pat)
# [[0 0 0 0 0]
#  [0 1 1 1 0]
#  [0 1 0 1 0]
#  [0 1 1 1 0]
#  [0 0 0 0 0]]
```

### How egenome affects eating

The **fiducial match count** for a single egene compares the actual
cell states in the 5x5 neighborhood against the fiducial pattern:

    matches(egene) = sum over all 25 positions of (v_actual == c_fiducial)

With multiple egenes per cell, eating uses the **best (maximum)** match
across the cell's active egenes:

    matches(cell) = max over active s of matches(egene[s])

The cell's **mouthful** is:

    M(x) = (m_scale / 25) * matches(cell) * F(x)

(capped to `1 - f_priv` so private food can't exceed 1).

The max-match rule is a "breadth bonus": more active egenes can only
increase mouthful, never decrease it. Without the per-egene tax
(`tax_per_egene`) this bonus would drive `Negene → NEGENOME_MAX` for
free; the tax provides the metabolic counterweight.

With egenome=0 (all-zero fiducial), matches counts **dead** neighbors.
With egenome=0b111111 (all-one), matches counts **alive** neighbors.

---

## Time Step Phases

Each call to `evoca_step()` executes five phases in order:

### Phase 1: CA State Update

Double-buffered.  For each cell, count active neighbors per ring,
compute the LUT bit index, look up the new state from the cell's
private LUT.  For alive cells only, the queried bit index is OR'd into
the `lut_active` mask (for restricted mutation).  Dead cells still
participate in neighbor counting (their zeroed LUT produces v=0).
Then swap `v_curr` and `v_next`.

If `restricted_mu` is enabled, the active bit indices are extracted
into `active_bits[]` for use in Phase 4.

### Phase 2: Environmental Food Regeneration

    F(x) += food_inc,   clamped to [0, 1]

### Phase 2b: Food Diffusion

If `gdiff > 0`, perform `gdiff` passes of 3x3 box blur on F:

    F_new(x) = (1/9) * sum of F over 3x3 Moore neighborhood including centre

Each pass uses a scratch buffer and pointer swap (double-buffered).
Periodic boundary conditions.

### Phase 2c: Tax and Death

For each alive cell, the per-step decrement is the sum of three terms:

    t = tax  +  tax_per_egene * Negene  +  tax_lut * popcount(LUT bytes)
    f(x) -= t,   clamped to 0

- `tax` is the unconditional baseline.
- `tax_per_egene` makes wider egenomes metabolically more expensive,
  which combined with max-match eating gives Negene a finite optimum.
- `tax_lut` taxes "1" bits in the LUT, penalising rule complexity so
  innovations that don't pay for themselves get pruned.

Defaulting `tax_per_egene = tax_lut = 0` reproduces the previous
single-rate tax exactly.

If `f(x)` reaches 0, the cell dies: `alive(x)` is set to 0, its LUT is
zeroed, the active mask is cleared (so `Negene = 0`), and `f(x)` is
cleared. Dead cells are skipped by subsequent phases.

### Phase 3: Eating

For each alive cell:
1. Compute fiducial matches (0-25) between 5x5 neighborhood and egenome
2. `mouthful = (m_scale / 25) * matches * F(x)` (proportional to available food)
3. Clamp to headroom: `mouthful = min(mouthful, 1.0 - f(x))`
4. Transfer: `F(x) -= mouthful`,  `f(x) += mouthful`

Private food f(x) is hard-capped at 1.0.

### Phase 4: Reproduction and Mutation

For each alive cell where `f(x) >= 1.0`:
1. Find the Moore neighbor (8 cells) with the lowest f(x').
   Ties broken by uniform random (xorshift32 PRNG).
2. Set `alive(child) = 1`. Copy parent genome to child: LUT, all 8
   egene slots, active mask. (v_curr is dynamical state, NOT copied.)
3. **Mutate child's LUT**: if `restricted_mu`, draw n_flips ~
   Poisson(mu_lut * n_active) and flip only active bits; otherwise
   draw n_flips ~ Poisson(mu_lut * 250) and flip random bits.
4. **Mutate child's egenome** (three sub-steps run in this order):
   - **Active-bit pass**: draw `na ~ Poisson(mu_egenome * 8)` flips
     across the 8 presence bits. Reject any flip that would make
     Negene = 0.
   - **Dup-on-activate**: for each presence bit that just went 0→1,
     with probability `p_dup_on_activate`, overwrite that slot's
     egene byte with a copy from a uniformly random
     currently-active slot.
   - **Egene-bit pass**: draw `ne ~ Poisson(mu_egene * 8 * 6)` flips,
     each picking a uniformly random slot and a uniformly random
     bit in [0, 6). Inactive slots get mutated too (pseudogenes).
   The duplicated copy from sub-step 2 sees fresh independent
   point mutations from sub-step 3, so duplication is followed by
   divergence.
5. Update child's caches: `lut_color` (LUT-only hash → ARGB) and
   `lut_hash_cache` (full genome hash = FNV-1a over LUT bytes ‖
   sorted active egene bytes).
6. Split food: `f(parent) = f(child) = f(parent) / 2`

The `births` array is set for each child cell (used by colormode 3).

---

## Visualization (Color Modes)

Five modes, selectable via dropdown or the `colormode` parameter. Each
mode uses ARGB pixel values; alpha is always `0xFF`. All `[0, 1]` floats
are clamped before the `× 255` cast.

| Mode | Name         | Encodes |
|------|--------------|---------|
| 0    | `state`      | alive flag × CA dynamical state, with genome-hash colour |
| 1    | `env-food`   | Environmental food on dead cells, private food on alive cells |
| 2    | `priv-food`  | Private food on alive cells (blue), CA state (red) |
| 3    | `births`     | Birth events highlighted, modulated by CA state |
| 4    | `age`        | Cell age via cool→hot log gradient |

### Mode 0: `state` — genome colour by CA state

| | dead | alive, v=0 | alive, v=1 |
|---|---|---|---|
| pixel | `0xFF000000` (black) | `0xFF333333` (dark grey) | `lut_color[i]` |

`lut_color[i]` is cached on each cell from an FNV-1a hash of its 32-byte
LUT. `set_lut_all()` stores the hash of the LUT it loads as the
**wild-type hash**; any cell whose LUT still hashes to that value
displays as white (`0xFFFFFFFF`). Mutant LUTs hash to other values and
get pseudo-random colours from the lower 24 bits of the hash. White
on screen therefore means "unmutated wild-type" — useful for tracking
how fast the original LUT decays under mutation.

### Mode 1: `env-food` — food field with hunger overlay on alive cells

The green channel encodes **environmental food** for dead cells and
**private food (cell's reserve)** for alive cells, since `F_food` is
eaten down to ~0 wherever organisms are sustained and so carries no
information there.

| | f_priv = 0 | f_priv = 0.5 | f_priv = 1 |
|---|---|---|---|
| **alive** (red=180, green=`f_priv`×255) | (180, 0, 0) red — starving | (180, 128, 0) orange | (180, 255, 0) yellow — about to reproduce |

| | F_food = 0 | F_food = 0.5 | F_food = 1 |
|---|---|---|---|
| **dead** (red=0, green=`F_food`×255) | (0, 0, 0) black | (0, 128, 0) mid-green | (0, 255, 0) bright green |

Visual story: dead cells brighten from black to green as environmental
food accumulates (at rate `food_inc`); alive cells flash red→orange→
yellow as their private reserve fills up between reproduction events.

### Mode 2: `priv-food` — private food + CA state on alive cells

Dead cells are black. Alive cells encode private food on the **blue**
channel; the **red** channel marks CA state (red=180 if v=1, else 0).

| | f_priv = 0 | f_priv = 0.5 | f_priv = 1 |
|---|---|---|---|
| **alive, v=0** (red=0, blue=`f_priv`×255) | (0, 0, 0) black | (0, 0, 128) mid-blue | (0, 0, 255) bright blue |
| **alive, v=1** (red=180, blue=`f_priv`×255) | (180, 0, 0) red | (180, 0, 128) magenta | (180, 0, 255) bright magenta |

Dead cells are pure black `0xFF000000` regardless of food.

### Mode 3: `births` — birth events highlighted

Highlights cells that have just reproduced this step. Mutant births
(LUT or egenome bits flipped) get magenta; clean copy births get
yellow; cells alive but with no birth this step are dim grey. The CA
state v_curr modulates brightness (v=1 → bright, v=0 → dim).

| | dead | alive, no birth | normal birth (`births=1`) | mutant birth (`births=2`) |
|---|---|---|---|---|
| **v=0** | `0xFF000000` black | `0xFF222222` dark grey | `0xFF808000` dim yellow | `0xFF800080` dim magenta |
| **v=1** | `0xFF000000` black | `0xFF444444` grey | `0xFFFFFF00` bright yellow | `0xFFFF00FF` bright magenta |

Useful for watching reproduction "fronts" advance and identifying
where mutation is most active.

### Mode 4: `age` — cool → hot log gradient

For each alive cell, age is `g_step − last_event_step[i]`, the number
of ticks since this cell was last replaced (born or reproduced into).
Dead cells are black.

A scale `g_age_scale` tracks the maximum observed alive-cell age,
ratcheting up instantly to any new max and decaying slowly (×0.995 per
tick, floor 10) so the gradient adapts to each run's age range without
wild jumps. The per-cell intensity is

$$v = \frac{\log(1 + \text{age})}{\log(1 + g_{\text{age scale}})} \in [0, 1].$$

Pixel channels are then

$$R = 80 + 175\,v,\quad G = 60 + 560\,v\,(1-v),\quad B = 40 + 180\,(1-v).$$

| age | $v$ | (R, G, B) | colour |
|--:|--:|---|---|
| 0 (just born / replaced) | 0.00 | (80, 60, 220) | cool blue-violet |
| ≈ scale × 0.1 | 0.50 | (167, 200, 130) | yellow-green |
| ≈ scale | 1.00 | (255, 60, 40) | hot red |
| dead | — | (0, 0, 0) | black |

Useful for spotting long-lived structures vs constant-turnover regions
in the same window.

---

## Probes, Activity Tracking, and Diagnostics

See **[Docs/probes.md](probes.md)** for full documentation of all probe
windows (activity charts, strip charts, click-to-identify overlays,
Y-axis scale controls, and the reproduction age histogram).

Call `available_probes()` from `python/controls.py` at runtime to list
all probe names and descriptions:

```python
from python.controls import available_probes
available_probes()
```

---

## Python API

### EvoCA class  (`python/evoca_py.py`)

#### Lifecycle

```python
sim = EvoCA(lib_path=None)
    # Load the shared library (auto-finds C/libevoca.dylib or .so)

sim.init(N, food_inc=0.0, m_scale=1.0, gdiff=0,
         mu_lut=0.0, mu_egene=0.0, mu_egenome=0.0,
         p_dup_on_activate=1.0,
         tax=0.0, tax_per_egene=0.0, tax_lut=0.0,
         restricted_mu=0)
    # Allocate N x N lattice, set metaparameters.
    # All grids initialized to zero.

sim.free()
    # Deallocate C-side memory.
```

#### Metaparameter Setters

```python
sim.update_food_inc(f)       # float
sim.update_m_scale(m)        # float
sim.update_gdiff(d)          # int
sim.update_mu_lut(m)         # float, per-bit LUT mutation rate
sim.update_mu_egene(m)       # float, per-egene-bit flip rate (across all 8×6 bits)
sim.update_mu_egenome(m)     # float, per-active-bit flip rate (across 8 presence bits)
sim.update_p_dup_on_activate(p)  # float in [0,1], dup-on-activate probability
sim.update_tax(t)            # float, constant priv-food decrement per step
sim.update_tax_per_egene(t)  # float, additional decrement per active egene
sim.update_tax_lut(t)        # float, additional decrement per LUT '1' bit
sim.update_restricted_mu(r)  # int (0 or 1)
sim.update_act_ymax(y)       # int, Y-scale for LUT activity chart
sim.update_eg_act_ymax(y)    # int, Y-scale for egenome activity chart
```

Each setter updates both the Python attribute (`sim.food_inc`, etc.)
and the C global.

#### Grid Setters

```python
sim.set_v(v_array)
    # Set cell states from (N, N) or flat uint8 array.

sim.set_lut_all(lut_bytes)
    # Set ALL cells' LUT from a single LUT_BYTES-length uint8 array.

sim.set_lut(idx, lut_bytes)
    # Set one cell's LUT (idx = flat cell index).

sim.set_egenome_all(eg)
    # Set all cells' fiducial genome (6-bit value, masked to 0x3F).

sim.set_egenome_random()
    # Set each cell's egenome to a random value in [0, 63].
    # No wild-type: all 64 egenomes get distinct hash-based colors.

sim.set_lut_random(n_init=3)
    # Set each cell's LUT to an independent random rule.
    # n_init: number of rings the rule conditions on (1, 2, or 3).
    #   1: depends on (v_x, n1) only — 10 independent bits per LUT
    #   2: depends on (v_x, n1, n2) — 50 independent bits
    #   3: depends on (v_x, n1, n2, n3) — all 250 bits independent

sim.set_f_all(f)
    # Set all cells' private food to float f.

sim.set_F_all(F)
    # Set all cells' environmental food to float F.

sim.set_F_random(lo=0.0, hi=1.0)
    # Set env food to uniform random floats in [lo, hi].

sim.set_env_mask(mask)
    # Set environment food-regeneration mask. mask: (N,N) or flat uint8.
    # 1 = regenerate, 0 = no regen. Default: all 1s.

sim.get_env_mask() -> np.ndarray
    # (N, N) uint8 copy of the environment mask.
```

#### Alive Data Plane

```python
sim.set_alive(arr)
    # Set alive array from (N,N) or flat uint8. Dead cells' v, f, LUT are zeroed.
    # Call LUT/egenome setters BEFORE this (alive setter zeroes dead cells' data).

sim.get_alive() -> np.ndarray
    # (N, N) uint8 copy of alive array.

sim.set_alive_all()
    # Set all cells alive (default after init).

sim.set_alive_fraction(frac)
    # Random fraction of cells alive; dead cells' data zeroed.

sim.set_alive_patch(radius)
    # Square patch of side 2*radius centered on grid; outside is dead.

sim.set_alive_halfplane(axis=0)
    # Half the grid alive. axis=0: left half, axis=1: top half.
```

**Initialization sequencing**: LUT and egenome setters must be called
*before* alive setters, because alive setters zero dead cells' LUT/v/f:

```python
# GoL in a central patch
sim.set_lut_all(make_gol_lut())
sim.set_egenome_all(0b000011)
sim.set_alive_patch(64)            # outside patch: alive=0, LUT/v/f zeroed

# Random rules on left half only
sim.set_lut_random(n_init=2)
sim.set_egenome_random()
sim.set_alive_halfplane(0)         # right half: dead
```

#### One-shot lattice initialization with `sim.state()`

`sim.state(**kwargs)` is a convenience that drives the LUT, egenome, v,
f_priv, F_food, and alive setters in the correct order from a single
keyword dict. Useful for parameter-sweep notebooks and as the
canonical state-init format used by `.evoca` recipe files
(`evoca_py.import_run` and the harness in `python/evoca_explore.py`).

```python
sim.state(
    lut='gol',          # 'gol' or 'random'  (anything ≠ 'random' = GoL)
    lut_n_init=3,       # if lut='random': 1=10 bits, 2=50, 3=250
    egenome='uniform',  # 'uniform' or 'random'
    egenome_value=0b000011,
    v_density=0.5,      # fraction of cells with v_curr=1
    f_init=0.5,         # initial private food
    F='uniform',        # 'uniform' or 'random'
    F_init=1.0,         # initial env food (when F='uniform')
    F_range=None,       # [lo, hi] (when F='random')
    alive='halfplane',  # 'all' | 'fraction' | 'patch' | 'halfplane'
    alive_fraction=0.5, # alive='fraction'
    alive_radius=64,    # alive='patch'
    alive_axis=0,       # alive='halfplane', 0=left, 1=top
)
```

All keys are optional and default to the values shown. Calling
`evoca_py.available_state_init()` returns the same descriptions as a
dict at runtime.

A separate state call **resets the lattice**: it overwrites the LUT
on every cell with the chosen rule, so any prior mutations are wiped.
If you want to start from GoL but **keep** the simulation actually
running as GoL, set `mu_lut=0` and `mu_egene=0` first (or switch
them off after `sim.state(...)` via `sim.update_mu_lut(0.0)`); otherwise
the per-birth Poisson mutation will drift the LUT away from the GoL
bytes within tens of generations.

#### Step and Colorize

```python
sim.step()
    # Execute one time step (all 5 phases).

sim.colorize(pixels, colormode=0)
    # Fill a (N*N,) int32 numpy array with ARGB pixel values.
    # colormode: 0=state, 1=env-food, 2=priv-food, 3=births
```

#### Getters

```python
sim.get_v()       -> np.ndarray   # (N, N) uint8, CA states (copy)
sim.get_alive()   -> np.ndarray   # (N, N) uint8, alive flags (copy)
sim.get_F()       -> np.ndarray   # (N, N) float32, env food (copy)
sim.get_f()       -> np.ndarray   # (N, N) float32, private food (copy)
sim.get_egenome()  -> np.ndarray   # (N, N) uint8, fiducial genomes (copy)
sim.get_births()  -> np.ndarray   # (N, N) uint8, birth events last step (copy)
sim.get_lut(idx)  -> np.ndarray   # (LUT_BYTES,) uint8, one cell's LUT (copy)
sim.get_step()    -> int           # current global step counter
sim.N             -> int           # lattice size
sim.cell_px       -> int           # CELL_PX compile-time constant

# Activity and diagnostics
sim.get_activity(max_n=4096) -> dict   # {'hash', 'activity', 'pop_count', 'color'}
sim.get_eg_activity()        -> dict   # {'activity', 'pop_count', 'color'} (64 entries)
sim.get_lut_complexity()     -> dict   # {'n1': count, 'n2': count, 'n3': count}
sim.get_repro_age_hist()     -> np.ndarray  # (1024,) uint32 histogram
sim.set_repro_age_t0(t)                # set step threshold for histogram accumulation
sim.reset_repro_age_hist()             # clear histogram

# Params export
sim.params()      -> dict         # current metaparameters as dict
sim.params_str()  -> str          # copy-pasteable sim.init(...) call
```

#### Recipe Export/Import

**Export** (from a running session):
```python
sim.export_recipe('my_run', probes={'activity': True}, colormode=0)
# -> Runs/2026-03-16_my_run.evoca
    # Records both initial metaparams (from init()) and final metaparams
    # (current values at export time, after any slider adjustments).
```

**Import and run**:
```python
from python.evoca_py import import_run

from python.controls import run_with_controls

import_run()                    # list available recipes in Runs/
sim, kw = import_run('Runs/2026-03-16_my_run.evoca')

sim, kw = import_run('Runs/2026-03-15_my_run.evoca', recipe='init')
    # recipe='init' (default): use the metaparams from sim.init()
    # recipe='final': use the metaparams at export time (after slider tweaks)
    # Returns (sim, display_kwargs).
    # Usage:  run_with_controls(sim, **kw)

run_with_controls(sim, **kw)

```

`import_run()` with no arguments prints available `.evoca` files and returns
a list of paths. With a filepath, it returns a `(sim, display_kwargs)` tuple.
The `recipe` argument selects which metaparams to use:
- `'init'` (default): the metaparams from the original `sim.init()` call
- `'final'`: the metaparams at export time (after any slider adjustments)

The recipe records initialization *methods* (not raw grid data), so each
import produces a new random realization with the same parameters. The
`.evoca` file stores both `metaparams_init` and `metaparams_final` so
either the original or the tweaked configuration can be reproduced.

#### Stored Attributes

After `init()` or the corresponding setter, these Python attributes
reflect current values:

    sim.food_inc, sim.m_scale, sim.gdiff,
    sim.mu_lut, sim.mu_egene, sim.tax, sim.restricted_mu, sim.egenome

---

## Controls and Display

### run_with_controls  (`python/controls.py`)

```python
run_with_controls(sim, cell_px=None, colormode=0, paused=False, probes=None)
    -> threading.Thread
```

Opens an SDL2 window and displays ipywidgets controls below the
notebook cell.  Returns immediately (non-blocking).

**Parameters**:

| Parameter   | Type  | Default        | Description                          |
|-------------|-------|----------------|--------------------------------------|
| `sim`       | EvoCA | (required)     | Initialized EvoCA instance           |
| `cell_px`   | int   | `sim.cell_px`  | Screen pixels per simulation cell    |
| `colormode` | int   | 0              | Initial color mode (0/1/2/3)         |
| `paused`    | bool  | False          | Start in paused state                |
| `probes`    | dict  | None           | Probe names to enable (see [probes.md](probes.md)) |

**Widgets**:
- **Run/Pause** toggle button
- **Restart** button (re-initializes state with current `state_params`)
- **Step** button (advances one step when paused)
- **Quit** button (closes SDL window and stops simulation)
- **Save Plots** button (saves probe strip chart images)
- **Export** text field + button: enter a descriptor and click Export to save a `.evoca` recipe
- **food_inc** slider: [0, 0.5], step 0.001
- **m_scale** slider: [0, 10], step 0.1
- **gdiff** slider: [0, 10], step 1
- **mu_lut** slider: [0, 0.001], step 0.00001
- **mu_egene** slider: [0, 0.05], step 0.001
- **tax** slider: [0, 0.1], step 0.001
- **act_ymax** / **eg_act_ymax** / **pat_act_ymax**: halve/double buttons (`<| name |>`)
- **restricted_mu** checkbox: toggle restricted mutation
- **Color** dropdown: state / env-food / priv-food / births
- **Save Plots** button: saves probe strip charts to PNG (when paused)
- **Export Params** button: prints copy-pasteable `sim.init(...)` call
- **Fiducial pattern**: 5x5 matplotlib grid displayed above widgets

**SDL2 window title** shows: time step, FPS (when running), color mode,
PAUSED indicator.

**Keyboard** (SDL2 window): Q or Esc to quit.

**Magnifier**: clicking on the main lattice window opens a 200×200 pixel
magnifier showing a 25×25 cell region at 8× zoom.  The magnifier acts as
a magnifying glass — dragging it over the lattice updates the displayed
region in real time.  Clicking elsewhere on the lattice jumps the
magnifier to that point.  The title bar shows the current center cell
coordinates `mag (row,col)`.  Close the magnifier via its title-bar close
button.  Near grid edges, the center is clamped so the 25×25 region
always fits without wrapping.

**Slider behavior**: dragging any slider auto-pauses the simulation;
200 ms after the last touch, it auto-resumes (unless already paused).

### Architecture

SDL2 requires the main thread on macOS.  Since Jupyter's kernel runs
cells in a worker thread, SDL2 is launched in a **subprocess**
(`python/sdl_worker.py`).  Pixel data and control signals flow via
POSIX shared memory (zero-copy):

```
Main process (Jupyter kernel)
+-- Jupyter event loop    -> ipywidgets callbacks
+-- Simulation thread     -> sim.step() + sim.colorize() -> pixel_shm
+-- Reader thread         -> relays SDL subprocess stdout
+-- subprocess (sdl_worker.py)
    +-- SDL2 main thread  -> reads pixel_shm, renders window
```

**ctrl_shm layout** (5 x int32):

| Index | Name      | Values                        |
|-------|-----------|-------------------------------|
| 0     | quit      | 0 = running, 1 = exit         |
| 1     | colormode | 0 / 1 / 2                     |
| 2     | step_cnt  | current time step             |
| 3     | fps * 10  | FPS scaled for int storage     |
| 4     | paused    | 0 = running, 1 = paused        |

---

## LUT Helper Functions

### make_gol_lut  (`python/evoca_py.py`)

```python
make_gol_lut() -> np.ndarray   # LUT_BYTES-length uint8
```

Builds Conway's Game of Life (B3/S23) as a bit-packed LUT.
- Dead cell (v_x=0): birth iff Moore count n1+n2 == 3
- Alive cell (v_x=1): survive iff Moore count n1+n2 in {2, 3}
- n3 is a don't-care (all values set identically).

With mutation=0 (all cells share this LUT), the simulation runs exact GoL.

### pack_lut / unpack_lut

```python
pack_lut(bits) -> np.ndarray
    # Pack LUT_BITS-length uint8 0/1 array into LUT_BYTES bytes.

unpack_lut(packed) -> np.ndarray
    # Unpack LUT_BYTES bytes into LUT_BITS-length uint8 0/1 array.
```

### egenome_to_pattern

```python
egenome_to_pattern(eg) -> np.ndarray   # (5, 5) uint8
    # Expand a 6-bit egenome value to a 5x5 binary array.
```

### lut_bit_index

```python
lut_bit_index(v_x, n1, n2, n3) -> int
    # Compute the flat bit index into the LUT.
```

### Constants

```python
LUT_BITS  = 250      # bits per LUT
LUT_BYTES = 32       # bytes per LUT (ceil(250/8))
ORBIT_MAP            # (5, 5) uint8 array of orbit indices
```

---

## Build

```bash
# macOS
gcc -O2 -Wall -fPIC -dynamiclib -o C/libevoca.dylib C/evoca.c

# Linux
gcc -O2 -Wall -fPIC -shared -o C/libevoca.so C/evoca.c

# or just:
make
```

The compile-time constant `CELL_PX` (default 2) in `C/evoca.h` sets
display scaling.  Change it and recompile.

---

## Architecture Notes

### Performance (CA step only, no food/repro, GoL LUT)

- N=256: ~150 fps
- N=512: ~30 fps

### Memory

Per cell: 32 (LUT) + 1 (egenome) + 1 (v_curr) + 1 (v_next) + 4 (f_priv)
+ 4 (F_food) + 4 (F_temp) + 1 (births) + 4 (lut_color) + 4 (lut_hash_cache)
+ 4 (last_event_step) = **60 bytes**.

- N=256: ~3.4 MB
- N=512: ~13.6 MB

### SDL2 on macOS

SDL2 video init (Cocoa/NSApplication) requires the actual main thread
(thread 0).  Jupyter's ipykernel runs sync cells in an executor thread,
so threading.Thread for SDL2 crashes the kernel.  The fix: run SDL2 in a
subprocess, which has its own main thread.
