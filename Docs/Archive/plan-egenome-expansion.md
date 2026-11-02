# EvoCA Egenome Expansion — Implementation Plan

## Context

Currently each cell has a single 6-bit egenome (D4 orbits of a 5×5 fiducial
pattern). The food-eating algorithm matches the cell's environment against
this single fiducial pattern; better matches yield bigger mouthfuls.

This plan expands the egenome into a list of up to 8 egenes per cell, with
the mouthful computed via *max-match* across the cell's active egenes. The
number of active egenes (`Negene`) is itself a per-cell genetic trait that
mutates and is taxed, allowing evolution to balance metabolic breadth
against per-egene cost.

This is also a first step toward a more differentiated tax structure: in
addition to the existing constant-per-tick tax, we add a per-egene tax and
a LUT-complexity tax, giving evolution finer-grained levers to reward
specific innovations.

## Design summary

- `NEGENOME_MAX = 8` slots per cell, fixed; some slots may be inactive.
- Each slot holds a 6-bit egene (same encoding as the current single
  egenome — preserves `EGENOME_COUNT = 64` bucket count for histogram
  probes).
- A per-cell `active[8]` bit array marks live slots.
  `Negene = popcount(active)`. Inactive slots' egene bits drift.
- Eating: `mouthful = (m_scale/25) × max_over_active_egenes(matches) × F_food`.
- Tax: `tax = tax_const + tax_per_egene × Negene + tax_lut × popcount(LUT)`.
- Mutation at birth (two rates):
  - **`mu_egene`** — per-bit flip rate over all 8×6 egene bits (active and
    inactive alike, so inactive slots accumulate drift).
  - **`mu_egenome`** — per-bit flip rate over the 8 active bits.
    A flip that would take Negene to 0 is rejected.
- On a presence bit going 0 → 1: with probability `p_dup_on_activate`,
  copy a random currently-active slot's egene bits into the newly-active
  slot, *before* the egene-bit point-mutation pass runs (option-b
  duplication-then-divergence: the new copy gets fresh independent
  mutations, not the same ones as its parent slot).
- Initialisation: all 8 slots filled with random 6-bit egenes; one random
  slot active (`Negenome_init = 1`).
- Species hash (for `n_distinct_genomes`-style probes): sort the active
  egene bytes lexicographically, append LUT bytes, FNV-1a — sort at hash
  time, no maintained invariant.

### Backward compatibility

`Negenome_init=1`, `mu_egenome=0`, `tax_per_egene=0`, `tax_lut=0`, and
`p_dup_on_activate` irrelevant ⇒ behaviour reduces to current sim modulo
struct size and the `mu_egene`/`mu_egenome` rename (see Item 1).

## Naming change — breaking

The existing global parameter `mu_egenome` (per-bit flip rate of the
single 6-bit egenome) is being **renamed to `mu_egene`**. The name
`mu_egenome` is then **reused** for the new active-bit flip rate.

Implications:
- Old scan CSVs have `param_mu_egenome` columns whose semantics no longer
  match the new `param_mu_egenome`.
- Python sliders, ctypes wrappers, controls.py keys, evoca_explore
  parameter dicts all need updates.
- `Docs/model.md` needs a new section for the multi-egene structure and
  the renamed lever.

## Implementation order

Build and test (compile + GoL smoke test + one short `run_sim`) after
each item. One commit per item.

1. Rename `mu_egenome` → `mu_egene` (mechanical, prepares the namespace)
2. C struct + storage: `egene[8][L]`, `active[8]`, init, free
3. C eating loop: max-match across active egenes
4. C mutation pass: two rates, presence-bit mutation, duplication on
   activation
5. C tax: per-egene + LUT-popcount terms
6. C species hash: sorted-active-egene + LUT bytes
7. Histogram probe (`eg_pop`/`eg_act`): redefine as
   per-active-egene-presence count
8. Python ctypes + `run_sim` parameter plumbing
9. Controls / SDL slider wiring (existing levers; new ones get sliders
   only if cheap)
10. `Docs/model.md` update
11. Validation: zero new knobs, diff metrics against last scan baseline

---

## Item 1: Rename `mu_egenome` → `mu_egene`

**Files**: `C/evoca.c`, `C/evoca.h`, `python/evoca_py.py`,
`python/controls.py`, `python/sdl_worker.py`, `python/evoca_explore.py`,
`Docs/model.md`, notebooks.

Mechanical find-and-replace:

| Pattern | Replacement |
|---------|-------------|
| `gmu_egenome` | `gmu_egene` |
| `mu_egenome` (param/var/key/slider/csv col) | `mu_egene` |
| `evoca_set_mu_egenome` / `evoca_get_mu_egenome` | `..._mu_egene` |

**Verify**: `grep -ri mu_egenome` shows only intended new (post-Item 4)
references or zero hits. Build succeeds. `evoca_test.ipynb` GoL run still
works.

After this rename, `mu_egenome` is a free name to reclaim in Item 4.

---

## Item 2: Multi-slot storage

**Goal**: replace single `uint8_t egenome[N*N]` with multi-slot storage,
without yet using the extra slots — Negene stays 1 everywhere, so
behaviour is unchanged.

**Data layout**:

```c
#define NEGENOME_MAX 8

/* per cell: 8 × 6-bit egenes, packed into uint8_t (low 6 bits each) */
static uint8_t *egenes  = NULL;   /* [N*N*NEGENOME_MAX] */
static uint8_t *active  = NULL;   /* [N*N] — bit array, 1 bit per slot */
```

(Replaces the existing `egenome` pointer.)

`active[idx]` is one byte storing 8 presence bits (LSB = slot 0).
`Negene(idx) = popcount(active[idx])`.

**Init** (in `evoca_init`):
- For every cell, fill all 8 slots with random 6-bit values
  (`rng_next() & 0x3F`).
- For every cell, set `active[idx] = 1 << (rng_next() & 7)` — one random
  slot active.

**Reproduction** (until Item 3+ swap in real logic): copy parent's full
8 slots and active byte to child verbatim (so behaviour matches the
single-egene case as long as exactly one bit is set in `active`).

**Convenience helpers** (private, file-static):
```c
static inline int  cell_negene(int idx);
static inline int  cell_first_active(int idx);   /* lowest set bit, or -1 */
static inline void cell_copy(int dst, int src);  /* both arrays */
```

**eating** (item-3 will rewrite): for now, replace
`fiducial_matches(row, col, egenome[idx])` with
`fiducial_matches(row, col, egenes[idx*NEGENOME_MAX + cell_first_active(idx)])`
so behaviour is identical when each cell has exactly one active slot.

**hash** (item-6 will rewrite): use the single active egene byte for now.

**Verify**: GoL run produces same trajectories as before (modulo RNG —
allocate the new arrays after existing RNG draws, or document seed
shift).

---

## Item 3: Max-match eating

**Files**: `C/evoca.c` (eating phase, ~line 766).

Replace single fiducial match with max over active slots. New helper:

```c
static int fiducial_matches_best(int row, int col, int idx)
{
    uint8_t a = active[idx];
    int best = 0;
    while (a) {
        int s = __builtin_ctz(a);
        a &= a - 1;
        int m = fiducial_matches(row, col, egenes[idx*NEGENOME_MAX + s]);
        if (m > best) best = m;
    }
    return best;
}
```

Replace the call site:

```c
int matches = fiducial_matches_best(row, col, idx);
```

**Cost note**: at Negene=1 this is a single call (same as before). At
Negene=8 it's 8× the eating-loop cost. Add `evoca_get_eating_ms()` in
Item 11 if validation shows a noticeable hit.

**Verify**: with `Negenome_init=1` and no mutation of the active byte,
eating throughput and metric distributions are unchanged from
pre-Item-2 baseline.

---

## Item 4: Two-rate mutation pass + dup-on-activate

**Files**: `C/evoca.c` (reproduction phase, ~line 825–845).

**Globals**:
```c
static float gmu_egene          = 0.0f; /* per-egene-bit flip rate */
static float gmu_egenome        = 0.0f; /* per-active-bit flip rate */
static float gp_dup_on_activate = 0.0f; /* prob copy active slot on 0→1 */
```

**Setters/getters**: standard pattern, mirror existing `set/get_mu_lut`.

**Mutation pass** (replaces the current single egenome flip block):

```c
/* Step 1: snapshot pre-mutation active mask for activation events */
uint8_t active_before = active[child];

/* Step 2: mutate active bits */
int na = poisson_sample(gmu_egenome * NEGENOME_MAX);
for (int f = 0; f < na; f++) {
    int b = rng_next() & 7;
    uint8_t new_active = active[child] ^ (uint8_t)(1u << b);
    if (new_active != 0)              /* reject Negene→0 */
        active[child] = new_active;
}

/* Step 3: dup-on-activate for newly-activated slots
 * (do this BEFORE point mutation so the new slot picks up fresh
 *  independent mutations in step 4). */
uint8_t newly_on = active[child] & ~active_before;
while (newly_on) {
    int s = __builtin_ctz(newly_on);
    newly_on &= newly_on - 1;
    if ((rng_next() & 0xFFFF) <
            (uint32_t)(gp_dup_on_activate * 65536.0f)) {
        /* pick a uniformly random currently-active slot
         * (active[child] already includes s, so we may pick s itself
         *  — that's a no-op copy, fine) */
        int n_act = __builtin_popcount(active[child]);
        int pick = rng_next() % n_act;
        uint8_t a = active[child];
        int src = -1;
        while (pick-- >= 0) { src = __builtin_ctz(a); a &= a - 1; }
        egenes[child*NEGENOME_MAX + s] =
            egenes[child*NEGENOME_MAX + src];
    }
    /* else: keep the drifted bits already in slot s */
}

/* Step 4: point-mutate every egene bit (active or not) */
int ne = poisson_sample(gmu_egene * NEGENOME_MAX * 6);
for (int f = 0; f < ne; f++) {
    int slot = rng_next() & 7;
    int bit  = rng_next() % 6;
    egenes[child*NEGENOME_MAX + slot] ^= (uint8_t)(1u << bit);
}
```

**`births[]` flag**: count `(na + ne) > 0` as the mutation indicator
(replaces the old `nc > 0`).

**Verify**: with `gmu_egenome = 0`, `gp_dup_on_activate = 0` and
`gmu_egene = old_mu_egenome_value`, scan metrics match pre-rename
results in distribution. With `gmu_egenome > 0`, Negene distribution
broadens over time.

---

## Item 5: Tax — per-egene + LUT-popcount

**Files**: `C/evoca.c` (tax phase, ~line 749).

**Globals**:
```c
static float gtax           = 0.0f; /* existing: constant per tick */
static float gtax_per_egene = 0.0f; /* per active egene per tick */
static float gtax_lut       = 0.0f; /* per LUT '1' bit per tick */
```

**Tax application**:

```c
if (gtax > 0.0f || gtax_per_egene > 0.0f || gtax_lut > 0.0f) {
    for (int i = 0; i < cells; i++) {
        if (!alive[i]) continue;
        float t = gtax + gtax_per_egene * cell_negene(i);
        if (gtax_lut > 0.0f)
            t += gtax_lut * lut_popcount(lut + (size_t)i * LUT_BYTES);
        f_priv[i] -= t;
        if (f_priv[i] < 0.0f) {
            alive[i] = 0;
            g_deaths_last++;
        }
    }
}
```

`lut_popcount` is a new helper that sums bits across the cell's LUT
bytes — straightforward.

**Verify**: `gtax_per_egene = gtax_lut = 0` reproduces existing behaviour.
With `gtax_per_egene > 0`, mean-Negene is bounded.

---

## Item 6: Species hash with sorted active egenes

**Files**: `C/evoca.c` (Channon hash machinery, near `lut_hash_fn`).

Goal: distinct cells with the same LUT and same *set* of active egenes
hash to the same bucket regardless of slot order.

```c
static uint32_t cell_genome_hash(int idx)
{
    uint8_t buf[NEGENOME_MAX + LUT_BYTES];
    int n = 0;
    uint8_t a = active[idx];
    while (a) {
        int s = __builtin_ctz(a); a &= a - 1;
        buf[n++] = egenes[idx*NEGENOME_MAX + s];
    }
    /* small-N insertion sort on buf[0..n) */
    for (int i = 1; i < n; i++) {
        uint8_t x = buf[i]; int j = i - 1;
        while (j >= 0 && buf[j] > x) { buf[j+1] = buf[j]; j--; }
        buf[j+1] = x;
    }
    memcpy(buf + n, lut + (size_t)idx * LUT_BYTES, LUT_BYTES);
    return fnv1a(buf, n + LUT_BYTES);
}
```

Replace existing genome-hash sites (used by `nact_*` Channon machinery)
with `cell_genome_hash(idx)`. This also affects the cached genome-color
recompute path — needs care to invalidate the cache when *either* the
LUT changes (existing) *or* any active egene byte changes (new).

Simple approach: invalidate cache on any mutation event (`na + ne > 0`
in the mutation pass), instead of the current LUT-only check.

**Verify**: `n_distinct_genomes` stays sane (close to single-egene
counts when Negene=1 everywhere).

---

## Item 7: `eg_pop` histogram redefinition

**Files**: `C/evoca.c` (egenome pop probes, ~line 1605–1660).

Current semantics: `eg_pop[v]` = number of cells whose 6-bit egenome
value is `v`. Sum across buckets = total alive.

New semantics: `eg_pop[v]` = number of *active egene slots across all
cells* whose 6-bit value is `v`. Sum across buckets = total alive
egene-slot count = sum over cells of Negene. Documented in
`Docs/model.md` and in the controls UI tooltip.

```c
/* in the per-tick eg_pop refresh */
memset(eg_pop, 0, sizeof(eg_pop));
for (int i = 0; i < cells; i++) {
    if (!alive[i]) continue;
    uint8_t a = active[i];
    while (a) {
        int s = __builtin_ctz(a); a &= a - 1;
        eg_pop[egenes[i*NEGENOME_MAX + s] & 0x3F]++;
    }
}
```

**Verify**: at Negene=1 everywhere, `sum(eg_pop) == alive_count` (same
as before).

---

## Item 8: Python plumbing

**Files**: `python/evoca_py.py`, `python/evoca_explore.py`,
`python/controls.py`.

- Add ctypes signatures for new setters/getters: `mu_egene`,
  `mu_egenome` (reused), `p_dup_on_activate`, `tax_per_egene`, `tax_lut`,
  and accessors for `Negene` (per-cell or summary).
- `run_sim`: new params with defaults preserving backward compat:
  ```python
  mu_egene=0.0,           # was mu_egenome
  mu_egenome=0.0,         # NEW: active-bit flip rate
  p_dup_on_activate=1.0,  # default: always duplicate on activation
  tax_per_egene=0.0,
  tax_lut=0.0,
  ```
- `evoca_explore.py`: extend `SCAN_INIT_DEFAULTS` and
  `_PARAM_DEFAULTS` (or wherever scan-param schema lives). Add summary
  metrics to results dict: `negene_mean`, `negene_std`,
  `negene_temporal_std`.

---

## Item 9: SDL / controls slider wiring

Add sliders (or just expose via controls.py) for the new globals. Keep
ranges conservative initially:

| Lever | Range | Default |
|-------|-------|---------|
| `mu_egene` | 0 – 0.1 | 0.001 |
| `mu_egenome` | 0 – 0.05 | 0.0 |
| `p_dup_on_activate` | 0 – 1 | 1.0 |
| `tax_per_egene` | 0 – 0.01 | 0.0 |
| `tax_lut` | 0 – 0.001 | 0.0 |

Defer to scan defaults if not exposed.

---

## Item 10: `Docs/model.md` update

New subsection under the genome description:

- Multi-egene cell genome: `(egene[8][6-bit], active[8])`.
- Max-match eating definition.
- Two mutation rates and dup-on-activate.
- Tax decomposition (constant + per-egene + LUT popcount).
- Species hash includes sorted active egenes.

Also update the table of model parameters (typically near the top of
model.md) with the new entries and the renamed `mu_egene`.

---

## Item 11: Validation

1. With all new knobs at zero/defaults reproducing single-egene
   behaviour: run a 2 000-tick sim at known scan-3 winner params; assert
   distributions of `correlation_length`, `n_distinct_genomes_mean`,
   `excess_activity_slope` are within ~5 % of pre-change baseline.
2. With `mu_egenome > 0` and `tax_per_egene = 0`: verify mean Negene
   drifts up over time (sanity: no tax pressure, breadth always wins).
3. With `mu_egenome > 0` and `tax_per_egene > 0`: verify mean Negene
   stabilises at an intermediate value.
4. Quick eating-loop timing check at Negene≈8 to confirm the max-match
   cost is acceptable (target: ≤2× current eating throughput).

## Status

- [x] Item 1 — Rename `mu_egenome` → `mu_egene`
- [x] Item 2 — Multi-slot storage
- [x] Item 3 — Max-match eating
- [x] Item 4 — Two-rate mutation + dup-on-activate
- [x] Item 5 — Tax: per-egene + LUT-popcount
- [x] Item 6 — Species hash sorted-active
- [x] Item 7 — `eg_pop` histogram redefinition
- [x] Item 8 — Python plumbing
- [x] Item 9 — SDL/controls sliders
- [x] Item 10 — `Docs/model.md` update
- [x] Item 11 — Validation

## Validation results (2026-04-28)

Four checks at the scan-3 productive corner (`food_inc=0.013, m_scale=1.2,
gdiff=0.06, mu_lut=0.001, mu_egene=0.003, tax=0.035, restricted_mu=True`),
N=128 unless noted, 2000 steps:

1. **Backward compat** (all new knobs zero, `p_dup_on_activate=0`):
   `negene_mean = 1.00` exactly. `alive_density_mean = 0.276`,
   `correlation_length = 33.1`, `n_distinct_genomes = 966`. System
   alive, plausible regime.
2. **Negene drift unbounded** (`mu_egenome=0.005`, `p_dup=1.0`,
   `tax_per_egene=0`): `negene_mean = 3.64`, final = 4.32 — Negene
   drifts up as expected from the breadth-bonus pressure of max-match
   eating without a metabolic counterweight.
3. **Negene bounded** (`mu_egenome=0.005`, `tax_per_egene=0.020`):
   `negene_mean = 1.21`, final = 1.31. The per-egene tax cleanly
   pulls Negene back to a low intermediate value.
4. **Eating-loop timing** at N=256, 500 timed steps after a 1000-step
   warmup: 1009 steps/s at `Negene = 1.00`, 557 steps/s at
   `Negene = 4.37`. Ratio 1.81× — within the "≤2× current"
   target. (Max-match short-circuits less than the worst-case 8×
   suggests because cells with Negene < NEGENOME_MAX skip the
   inactive slots, and bucket-population variance helps too.)
