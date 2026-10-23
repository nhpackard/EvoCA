# plan-chem.md — the vector-valued resource ("chem") branch

Semi-permanent, collaborator-facing feature branch `chem`. Implements the
checkpoint §E / §I′ actions **E1** (branch + plan), **E3** (design: sparse
storage + vector→`f` mapping pinned up front), and **L3** (branch + single-field
parity gate). Source anchor: Farmer–Kauffman–Packard 1986 (autocatalytic
networks); see `Docs/references/references.md`. `~/Projects/swarmchem` is **not**
a source — it has no reaction network (checkpoint §E LB4).

This document pins the two design decisions that had to be settled **before**
any code, and records exactly what milestone 1 does and does not do.

---

## 0. Scope of this commit (milestone 1 only)

Milestone 1 = **the general `N_chem` machinery skeleton + the N_chem=1 parity
gate.** It lays down the storage and the vector→scalar mapping, and proves that
at `N_chem == 1` the model reproduces the current single-food results
**bit-for-bit**. It does **not** implement reactions, food/waste, or the
egene→compound-profile redesign — those are milestones 2 and 3.

Milestone roadmap (from checkpoint §K′ / §E′ note, N's call 2026-07-19):

| Milestone | Content | Status |
|---|---|---|
| **1 — parity** | `N_chem` machinery + sparse store + vector→`f`; N_chem=1 == single-food bit-for-bit | **this commit** |
| 2 — food/waste | `N_chem=2`; minimal reaction-execution engine (food→waste); waste `W_env` field, disposal ∝ (1−W_env); §B structured-gradient experiment | **implemented 2026-10-23 (branch `chem`; §7)** |
| 3 — larger spaces | autocatalytic reaction-network design (FKP-style); egene→compound-profile eating | deferred |

Each milestone gates the next (deliberately serial — checkpoint §K″).

---

## 1. PINNED DECISION A — sparse agent-local resource storage

**Decision: store each agent's resource as a SPARSE per-cell active-compound
slot list, mirroring the egene slots — NOT a dense per-compound plane.**

This is the tractability decision (checkpoint §E, §J′, LB "single most
consequential structural decision"). A dense scheme costs `O(N_chem · N²)`
memory *touched every step*; at `N_chem=1000, N=512` that is a quarter-billion
cells streamed per step and misses the fps target. But under the "seed sparse,
let reactions populate a fixed universe" dynamics, **most compounds are absent
at most sites almost all the time** — occupancy is sparse — so the store is
`O(active compounds per cell)`, independent of the universe size `N_chem`.

### Layout (C/evoca.c, C/evoca.h)

```
#define CHEM_MAX 8            /* per-cell active-compound slot capacity   */
chem_ns [N*N]                 /* uint8  live slot count for the cell      */
chem_id [N*N*CHEM_MAX]        /* uint8  compound id in each slot          */
chem_amt[N*N*CHEM_MAX]        /* float  amount     in each slot           */
```

`CHEM_MAX` (per-cell slot cap) is deliberately decoupled from `N_chem`
(universe size): memory is `O(N² · CHEM_MAX)`, never `O(N² · N_chem)`. At
N=512 the store is ≈ 10.5 MB regardless of whether the universe has 2 or 10 000
compounds. This mirrors the egene design exactly (`NEGENOME_MAX=8` slots, an
8-bit active mask), which is the same sparse-slot pattern the model already
uses — so the branch adds no new structural idiom.

Optional **dense per-compound planes for visualization** (checkpoint §E, N's
reply (ii)) are a *viz projection* rendered on demand, not the hot-path
representation, and are deferred to when there is >1 compound to show.

### At N_chem=1

Every live cell has exactly one slot: `chem_id=0`, `chem_amt=` the agent's
food. `chem_ns=1` for alive cells, `0` for dead. The slot machinery is fully
exercised (allocated, deposited into, split, cleared, freed) even though the
universe has one compound.

---

## 2. PINNED DECISION B — the vector→scalar `f` mapping

**Decision: the scalar private food `f` that gates reproduction is DEFINED as
the linear reduction**

```
f(cell) = Σ_slots  chem_energy[compound_id] · amount
```

with a per-compound weight vector `chem_energy[]` (default `1.0` for every
compound). This is the "how does multi-compound intake become the private food"
mapping the checkpoint (§E) flagged as the actual scientific core and the one
thing to pin before coding.

### Why linear, and why an energy weight

- **Parity is exact and provable.** At `N_chem=1`, one slot (id 0) with
  `chem_energy[0]=1.0`, so `f = 1.0 · amount`. The reduction is the identity;
  `f` tracks the pre-change scalar to the last bit (`1.0f * x == x` exactly in
  IEEE-754).
- **Milestone 2 is a one-line generalization.** food+waste sets
  `chem_energy = [1, 0]` (food is energy, waste is inert) — waste accumulates
  in a slot but contributes nothing to the reproduction gate, which is exactly
  the intended semantics.
- The weight vector is the natural home for "compound X is high-free-energy"
  (checkpoint §L, the ATP/activation-compound idea) when milestone 3 arrives.

### `f_priv` becomes a cached reduction (not a second source of truth)

`f_priv[N*N]` is retained as the **hot-read cache** of the reduction: every
place that currently reads private food (eating headroom `1−f`, the repro gate
`f≥1`, the 50/50 split, tax death `f≤0`, every probe, and the Python `get_f()`
getter) keeps reading `f_priv` with no per-read reduction cost. The store is the
**storage of record for compound amounts**; `f_priv` is recomputed
(`f_priv[idx] = chem_reduce_to_f(idx)`) immediately after any store mutation.
At `N_chem=1` the cache equals the amount, so it is redundant — but at
`N_chem≥2` the reduction is non-trivial and the cache earns its keep, so this is
forward-looking structure, not dead duplication. The invariant
`f_priv[idx] == chem_reduce_to_f(idx)` holds after every phase.

### Eating at milestone 1

Unchanged from single-food: `mouthful = (m_scale/25) · matches · F_food`,
clamped to headroom, subtracted from `F_food`, deposited into compound 0. The
egene→compound-profile redesign (does an egene target a compound profile or a
spatial pattern?) is deferred to milestone 3; milestones 1–2 reuse the current
spatial-pattern eating on the food component. Resource communication uses the
existing (enlarged Moore) eating neighbourhood — no new messaging layer
(checkpoint §E reply (iii)).

---

## 3. The `N_chem` guard (surfacing the limitation, not hiding it)

`evoca_set_n_chem(n)` implements only `n == 1`. Any other value is **clamped to
1 and reported** on stderr (returns 0), rather than silently pretending to
support a universe whose reactions are not yet built. This follows the
project rule: when the spec asks for more than is implemented, surface it. The
guard is relaxed milestone-by-milestone as the reaction engine lands.

---

## 4. Parity gate (milestone 1 deliverable)

- `tests/chem_parity.py` — shared harness: fixed deterministic initial state
  (numpy-seeded LUT, C-seeded egene/alive init, C PRNG re-seeded before the
  step loop), 200 steps on N=48 with a config that fires eating, tax death,
  reproduction, and LUT+egene mutation. Captures `f`, `F`, `alive`, and per-cell
  LUT hashes at steps 50/100/150/200.
- `tests/data/chem_parity_golden.npz` — golden snapshot captured from the
  **pre-change** build via `chem_parity.regen_golden()`.
- `tests/test_chem_parity.py` — asserts the current (chem) build reproduces the
  golden **bit-for-bit** (`np.array_equal`, not `allclose`), plus that N_chem
  defaults to 1, clamps unsupported values, and `chem_energy` defaults to unit.
- `chem-test.ipynb` — the same parity check as a runnable PASS/FAIL demo on a
  fixed seed.

The single strongest guard: the N_chem=1 hot path adds **zero PRNG draws**
relative to the pre-change path (the deposit/tax/split helpers are pure
arithmetic). Any extra draw would desynchronize the reproduction RNG stream and
fail the ~2300-cell / 200-step bit-exact snapshot immediately.

**Result: parity PASSES** — `f`, `F`, `alive`, `lut_hash` all bit-identical at
every checkpoint.

---

## 5. What is intentionally NOT here (deferred, by design)

- No reactions / reaction-execution engine (milestone 2).
- No food/waste second compound or `W_env` field (milestone 2).
- No autocatalytic-network specification (milestone 3).
- No egene→compound-profile eating (milestone 3).
- No dense per-compound visualization planes (add with milestone 2+).
- `chem_energy` non-unit weights are wired but exercised only from milestone 2.

## 6. Files changed on this branch

- `C/evoca.h` — `CHEM_MAX`/`CHEM_MAX_N`; `evoca_set/get_n_chem`,
  `evoca_set/get_chem_energy` prototypes.
- `C/evoca.c` — chem globals + store; `chem_reduce_to_f` / `chem_cell_set` /
  `chem_cell_add` / `chem_cell_clear` helpers; alloc/free; setters/getters;
  routed Phase 2c (tax), Phase 3 (eating), Phase 4 (split), `set_f_all`, and the
  `set_alive_*` dead-cell zeroing through the store.
- `python/evoca_py.py` — ctypes bindings; `n_chem` init param; `set/get_n_chem`,
  `set/get_chem_energy` methods.
- `tests/chem_parity.py`, `tests/test_chem_parity.py`,
  `tests/data/chem_parity_golden.npz` — the parity gate.
- `chem-test.ipynb` — PASS/FAIL parity demo.
- `Docs/plan-chem.md` — this document.

---

## 7. Milestone 2 — food/waste as a special case of the general chemistry (2026-10-23)

N's spec (checkpoint, LB6) with the `W_env` refinement and N's decisions of
2026-10-23: environmental waste mirrors environmental food, as a **sink**
instead of a source, with its own diffusion; on death a cell's food and waste
return to its site.

### General mechanism (C, every compound c)

| | compound 0 (energy, "food") | compound c >= 1 |
|---|---|---|
| environment field | `F_food` | `chem_env[c]` (allocated by `evoca_set_n_chem`) |
| source per step | `food_inc` on `env_mask`, `F <= 1` cap | `chem_env_inc[c]`, uniform |
| sink per step | none | `chem_env_sink[c]`, uniform, removal floored at 0 |
| diffusion | `gdiff` | `chem_env_diff[c]` (same box-blur rule) |
| intake | eating (egene match, capped mouthful) | none (milestone 3) |
| excretion | none | `d = chem_disp[c] * a * max(0, 1 - E)` into the cell's own site |
| toxicity | none | private amount `>= chem_tox[c]` (> 0) kills |
| energy weight in `f` | 1 (fixed) | 0 or 1 |

**Metabolism, the one reaction:** the tax converts compound 0 into
`g_metab_product` (the energy actually consumed: `t`, or the whole store if
smaller). `-1` = into nothing, which is the single-food model.

**Death:**
- *Tax death:* the cell's compounds `c >= 1` return to its site; its energy is `<= 0`.
- *Toxic death:* everything returns, energy into `F_food`. This is a transfer, not a dissipation sink.
- *Eviction:* the `c >= 1` compounds return; the energy is discarded as before (sink (b)).

**Reproduction** splits every compound 50/50.

**Step order:**

1. CA.
2. Food regeneration and diffusion.
3. Environment fields `c >= 1`: source, then sink, then diffusion.
4. Tax and metabolism (tax death).
5. Excretion and toxicity (toxic death).
6. Eating.
7. Reproduction.

`N_chem` may be 1..`CHEM_MAX_N` (8) in C. At `N_chem == 1` none of the new code runs: the single-food paths are textually unchanged.

### Food/waste preset (Python, `EvoCA.init(n_chem=2, wdiff, waste_sink, waste_disp)`)

| waste (compound 1) | value |
|---|---|
| energy weight | 0 |
| toxic at | 1 (N's spec: "penalty for waste reaching max (1.0): agent dies") |
| metabolic product of the tax | yes (N's spec: tax changes food to waste) |
| environment: source, sink, diffusion | 0, `waste_sink`, `wdiff` |
| excretion rate | `waste_disp` |

The three waste metaparameters go through recipes, `run_sim` (`_PARAM_KEYS`) and
gallery overrides like any other; they require `n_chem=2`. Getters: `get_W()`
(environment), `get_w()` (private), `get_chem_flux(1)`, `get_tox_deaths_last()`;
`set_w()` seeds private waste. Colormode 6 "food-waste" shows F (green), W (red)
and alive (blue).

### Tests (`tests/test_chem_waste.py`, 17)

- **Inert path.** `N_chem=2` with every compound inert is bit-for-bit the single-food run, which exercises the general split, excretion and eviction paths.
- **Constructed one-step fixtures.** Metabolism and excretion, excretion stopping at saturation, the sink, toxic death returning food and waste, and the 50/50 split.
- **Conservation every step on evolving runs.** The waste mass identity `d(ΣW + Σw) = in − out`, and the unchanged energy identity.
- **Public doors.** Init (and misuse), `run_sim`, recipe round trip, gallery override and table, and the colormode channels.

Each guard was shown to fail with the corresponding piece removed:

- no eviction return → 3 identity tests fail;
- toxic death losing food → 2 tests fail;
- waste keys missing from `_PARAM_KEYS` → the `run_sim` test fails.

The milestone-1 parity golden and all earlier goldens still pass.

### Not in milestone 2

- intake of compounds other than food;
- wgenes (egene-like disposal genes);
- source/sink masks for waste;
- a reaction network (milestone 3);
- energy weights other than 0/1.
