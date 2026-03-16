# EvoCA To-Do Implementation Plan

## Context

The model overview in `Docs/model.md` has 6 To-Do items that evolve the EvoCA model: renaming, documentation, a new environment mask, simplifying reproduction, adding an alive/dead cell distinction, and export/import. These are ordered to minimize rework: small cleanups first (1-4), then the major architectural change (5), then export/import last (6) since it serializes the final API.

## Implementation Order

1. Rename `cgenom` -> `egenome` (mechanical)
2. Move GoL LUT docs into LUT Indexing section (docs only)
3. Add environment mask `em(x)` (new C array + Python wrapper)
4. Fix reproduction threshold to 1.0, remove `food_repro` (removal)
5. Add `alive(x)` data plane (major architectural change)
6. Export/import recipes to JSON (new feature, depends on final API)

**Build and test after each item.** One commit per item.

---

## Item 1: Rename cgenom -> egenome

**Files**: `C/evoca.c`, `C/evoca.h`, `python/evoca_py.py`, `python/controls.py`, `python/sdl_worker.py`, `Docs/model.md`, `CLAUDE.md`, `evoca_test.ipynb`, `evoca.ipynb`

Global find-and-replace of all identifiers, parameters, attributes, and comments:

| Pattern | Replacement | Scope |
|---------|-------------|-------|
| `cgenom` | `egenome` | C arrays, Python attrs, function names, parameters |
| `mu_cgenom` | `mu_egenome` | C global, Python attr, slider, docs |
| `cg_act`, `cg_pop`, `cg_color` | `eg_act`, `eg_pop`, `eg_color` | C activity arrays |
| `cg_activity` | `eg_activity` | probe name, shm name, controls, sdl_worker |
| `cg_act_ymax` | `eg_act_ymax` | C global, slider, Python method |
| `evoca_*cgenom*` / `evoca_cg_*` | `evoca_*egenome*` / `evoca_eg_*` | C functions, header decls |
| `cgenom_to_pattern` | `egenome_to_pattern` | Python helper |
| `wt_cgenom_val` | `wt_egenome_val` | C wild-type tracker |
| `CGENOM_COUNT` | `EGENOME_COUNT` | C #define |

**Verify**: `grep -ri cgenom` returns zero hits. Build succeeds. GoL test passes.

---

## Item 2: GoL LUT Documentation

**File**: `Docs/model.md`

After the existing line "GoL (B3/S23) is exactly encodable..." in the LUT Indexing section (~line 119), add a paragraph explaining how `make_gol_lut()` constructs the B3/S23 rule (dead: birth iff n1+n2==3; alive: survive iff n1+n2 in {2,3}; n3 don't-care). Keep the existing detailed entry under LUT Helper Functions.

---

## Item 3: Environment Mask em(x)

**Files**: `C/evoca.c`, `C/evoca.h`, `python/evoca_py.py`, `Docs/model.md`

**C side**:
- New global: `static uint8_t *env_mask = NULL;`
- `evoca_init()`: allocate, `memset` to 1 (default: all sites regenerate)
- `evoca_free()`: free it
- Phase 2 food regen: `if (env_mask[i]) { F_food[i] += ...; }`
- New functions: `evoca_set_env_mask(const uint8_t *mask)`, `evoca_get_env_mask(void)`

**Python side**:
- ctypes bindings for both functions
- `sim.set_env_mask(mask)`: takes (N,N) or flat uint8
- `sim.get_env_mask()`: returns (N,N) uint8 copy

**Static mask** — set once, no slider. `em(x)=1` everywhere is default (current behavior).

**Verify**: default behavior unchanged. Set half to 0, confirm food only grows on mask=1 side (colormode 1).

---

## Item 4: Fix food_repro to 1.0

**Files**: `C/evoca.c`, `C/evoca.h`, `python/evoca_py.py`, `python/controls.py`, `Docs/model.md`, notebooks

**Remove**:
- C: `gfood_repro` global, `evoca_set_food_repro()`, `food_repro` param from `evoca_init()`
- Python: `food_repro` param from `init()`, `update_food_repro()`, `self.food_repro`, `_DEFAULTS` entry, `params()`/`params_str()` entries
- Controls: `sl_food_repro` slider and callback
- Notebooks: `food_repro=` from any `sim.init()` calls

**Hardcode**: Phase 4 threshold: `if (f_priv[idx] < 1.0f) continue;`

**C init signature change**: `evoca_init(int N, float food_inc, float m_scale)` (3 params, was 4)

**Verify**: build succeeds, reproduction triggers only when f(x) >= 1.0.

---

## Item 5: Add alive(x) Data Plane

**Files**: `C/evoca.c`, `C/evoca.h`, `python/evoca_py.py`, `Docs/model.md`, `CLAUDE.md`

### New C data

- `static uint8_t *alive = NULL;` — N*N array, 1=alive organism, 0=dead slot
- Allocated in `evoca_init()`, initialized to all 1s (backward compatible)
- Freed in `evoca_free()`
- No double-buffering needed (alive changes don't conflict within a step)

### Phase modifications in evoca_step()

| Phase | Change |
|-------|--------|
| **1 (CA update)** | Only record `lut_active` bits for alive cells. Dead cells still produce v_next=0 from zeroed LUT (no code change needed for v_next). |
| **2 (food regen)** | No change — environmental, happens everywhere (subject to env_mask) |
| **2c (tax/death)** | Skip dead cells (`if (!alive[i]) continue`). On death: set `alive[i]=0`, zero LUT, zero f_priv. |
| **3 (eating)** | Skip dead cells (`if (!alive[idx]) continue`) |
| **4 (reproduction)** | Skip dead parents (`if (!alive[idx]) continue`). Set `alive[child]=1` when reproducing into a cell. |

### Colorization changes

| Mode | Current | New |
|------|---------|-----|
| 0 (state) | v=1: genome color, v=0: black | alive+v=1: genome color, alive+v=0: very dark grey (0xFF111111), dead: black |
| 1 (env-food) | green=F(x), red tint if v=1 | green=F(x) always, red tint if alive (not v) |
| 2 (priv-food) | blue=f(x), red tint if v=1 | dead: black. alive: blue=f(x), red tint |
| 3 (births) | v=0: black | dead: black, alive+v=0: dim variants of birth/no-birth colors |

Mode 0 note: alive+v=0 gets 0xFF111111 (nearly black) so pure GoL runs look almost identical to current.

### Activity tracking

Change `if (!v_curr[i]) continue;` to `if (!alive[i]) continue;` in:
- `evoca_activity_update()`
- `evoca_eg_activity_update()` (renamed from cg_)
- `evoca_lut_complexity_counts()`

Note: `evoca_pat_update()` was kept scanning all sites since it measures lattice-wide neighbourhood pattern entropy, not organism-specific properties.

### New C functions

```c
void     evoca_set_alive(const uint8_t *arr);
uint8_t *evoca_get_alive(void);
void     evoca_set_alive_all(void);          // all cells alive
void     evoca_set_alive_fraction(float f);  // random fraction alive
void     evoca_set_alive_patch(int radius);  // square patch at center
void     evoca_set_alive_halfplane(int axis); // 0=left, 1=top
```

For fraction/patch/halfplane: dead cells get zeroed LUT, v=0, f=0.

### New Python methods

`set_alive(arr)`, `get_alive()`, `set_alive_all()`, `set_alive_fraction(f)`, `set_alive_patch(radius)`, `set_alive_halfplane(axis=0)`

### Documentation

- Grid State table: add `alive(x)` (uint8), change `v(x)` description to "Binary CA state (not life/death)"
- Update all phase descriptions
- Eliminate all "v(x)=0 (dead)" / "v(x)=1 (alive)" language — v(x) is just CA state
- Add new Python API methods
- Memory: +1 byte/cell

### Documentation: initialization sequencing examples

The docs must include examples showing that LUT/egenome setters are called **before** alive setters (which zero dead cells' data):

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

**Verify**: GoL mode looks correct (v=0 alive cells nearly invisible). Food+tax run: dead cells black, reproduction fills dead slots. `set_alive_fraction(0.5)` / `set_alive_patch(64)` / `set_alive_halfplane(0)` all work visually.

---

## Item 6: Export/Import Recipes

**Files**: `python/evoca_py.py`, `python/controls.py`, `Docs/model.md`

### Recipe format (JSON, `.evoca` extension)

```json
{
  "version": 1,
  "created": "2026-03-15T12:34:56",
  "descriptor": "user description",
  "N": 256,
  "metaparams": { "food_inc": ..., "m_scale": ..., "gdiff": ...,
                   "mu_lut": ..., "mu_egenome": ..., "tax": ...,
                   "restricted_mu": 0 },
  "initialization": {
    "lut_method": "gol" | "random",  "lut_n_init": 1|2|3|null,
    "egenome_method": "uniform" | "random",  "egenome_value": int|null,
    "v_method": "random",  "v_density": 0.5,
    "f_init": 0.0,
    "F_method": "uniform" | "random",  "F_init": float,  "F_range": [lo,hi],
    "alive_method": "all" | "fraction" | "patch" | "halfplane",
    "alive_fraction": float|null, "alive_radius": int|null, "alive_axis": int|null,
    "env_mask": null  (or base64-encoded if non-trivial)
  },
  "display": { "colormode": 0, "probes": {...} }
}
```

### Recipe tracking

Add `self._recipe = {}` to EvoCA.__init__. Each setter (`set_lut_all`, `set_lut_random`, `set_egenome_all`, etc.) records its method/args in `_recipe['initialization']`.

### Export

`sim.export_recipe(descriptor, probes=None, colormode=0)` -> writes `Runs/YYYY-MM-DD_descriptor.evoca`

Creates `Runs/` directory if needed.

### Import

Module-level `import_run(filepath)` -> `(sim, display_kwargs)`. Reads JSON, constructs EvoCA, calls init + setters, returns sim and kwargs for `run_with_controls`.

### Controls widget change

Replace "Export Params" button with HBox: `Text(placeholder='descriptor')` + `Button('Export')`. On click, calls `sim.export_recipe(text.value, probes=probes, colormode=current)`.

**Verify**: export a recipe, check JSON in `Runs/`. Import it, run `run_with_controls(sim, **display_kwargs)`, confirm it works. Test GoL and random-LUT recipes.

---

## Critical Files

| File | Items |
|------|-------|
| `C/evoca.c` | 1, 3, 4, 5 |
| `C/evoca.h` | 1, 3, 4, 5 |
| `python/evoca_py.py` | 1, 3, 4, 5, 6 |
| `python/controls.py` | 1, 4, 6 |
| `python/sdl_worker.py` | 1 |
| `Docs/model.md` | 1, 2, 3, 4, 5, 6 |
| `CLAUDE.md` | 1, 5 |
| `evoca_test.ipynb`, `evoca.ipynb` | 1, 4 |

---

## Implementation Summary (completed 2026-03-15)

All 6 items implemented in commit `a449cd2` on branch `softmouth`.

### Item 1 — Rename cgenom → egenome
Global find-replace across all C, Python, docs, and notebooks. All identifiers, function names, parameters, comments, section headers, and probe names updated. `grep -ri cgenom` returns zero hits (except archival `Docs/plan-reduce-lut-3rings.md`).

### Item 2 — GoL LUT docs
Added initialization paragraph in the LUT Indexing section of `Docs/model.md` explaining how `make_gol_lut()` constructs B3/S23.

### Item 3 — Environment mask em(x)
New `env_mask` uint8 array in C (`evoca.c`), allocated in `evoca_init()`, freed in `evoca_free()`. Phase 2 food regen checks `if (!env_mask[i]) continue;`. Python `set_env_mask()`/`get_env_mask()` wrappers added. Default all-1s preserves current behavior.

### Item 4 — Fix food_repro to 1.0
Removed the `food_repro` metaparameter entirely:
- C: removed `gfood_repro` global, `evoca_set_food_repro()`, `food_repro` param from `evoca_init()`. Phase 4 threshold hardcoded to `1.0f`.
- Python: removed from `init()`, `_DEFAULTS`, `params()`, `update_food_repro()`.
- Controls: removed `sl_food_repro` slider and callback.
- Notebooks: removed `food_repro=` from all `sim.init()` calls.
- Docs/CLAUDE.md: all references updated.
- `evoca_init` signature: `(int N, float food_inc, float m_scale)` (3 params, was 4).
- `python/display.py`: removed `repro=` from window title.

### Item 5 — alive(x) data plane
New `uint8 alive[N*N]` array separating organism existence from CA state `v(x)`:
- Dead cells (`alive=0`) have zeroed LUT, `v=0`, `f=0`.
- Dead cells don't eat, are skipped by tax, and cannot reproduce.
- Reproduction is the only way a dead cell becomes alive (`alive[child]=1`).
- Tax death sets `alive[i]=0`.
- Phase 1: only alive cells contribute to `lut_active` mask.
- Activity tracking (`evoca_activity_update`, `evoca_eg_activity_update`, `evoca_lut_complexity_counts`): changed to check `alive[i]` instead of `v_curr[i]`.
- `evoca_pat_update` left unchanged (scans all sites for lattice-wide pattern entropy).
- Colorization updated: mode 0 shows alive+v=0 as dark grey (#111111), dead as black. Mode 1 uses alive (not v) for red tint. Mode 2 shows dead as black. Mode 3 shows dead as black, alive no-birth v=0 as dark grey (#222222).
- New initialization helpers: `set_alive_all()`, `set_alive_fraction(f)`, `set_alive_patch(radius)`, `set_alive_halfplane(axis)`. All zero dead cells' LUT/v/f.
- Python wrappers: `set_alive(arr)`, `get_alive()`, plus the four helpers above.
- Documentation updated: Grid State table, all phase descriptions, initialization sequencing examples.

### Item 6 — Recipe export/import
- `self._recipe = {}` in `EvoCA.__init__`. Each setter (`set_lut_all`, `set_lut_random`, `set_egenome_all`, `set_egenome_random`, `set_v`, `set_f_all`, `set_F_all`, `set_F_random`, `set_alive_*`) records its method/args.
- `sim.export_recipe(descriptor, probes, colormode)` writes JSON to `Runs/YYYY-MM-DD_descriptor.evoca`.
- Module-level `import_run(filepath)` reads JSON, constructs EvoCA, calls init + setters, returns `(sim, display_kwargs)`.
- Controls widget: replaced "Export Params" button with text field + "Export" button. On click, calls `sim.export_recipe()` and shows filepath in status label.
- Tested round-trip: GoL recipe and random-LUT + alive_patch recipe both export and import correctly.
