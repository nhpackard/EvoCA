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
| energy weight in `f` | 1 (fixed) | 0 (fixed until milestone 3) |

**Metabolism, the one reaction:** the tax converts compound 0 into
`g_metab_product` (the energy actually consumed: `t`, or the whole store if
smaller). `-1` = into nothing, which is the single-food model.

**Death:**
- *Tax death:* the cell's compounds `c >= 1` return to its site; its energy is `<= 0`.
- *Toxic death:* everything returns, energy into `F_food`. This is a transfer, not a dissipation sink.
- *Eviction:* the `c >= 1` compounds return; the energy is discarded as before (sink (b)).
- *Outside the step* (`set_alive*`, `stamp_patch` clearing or replacing an occupant): the `c >= 1` compounds return to the site, as on any death.

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
- energy-carrying compounds `c >= 1` (weights other than 0).

## 8. Review fixes (2026-10-23)

**Summary.** Two reviews of milestone 2 ran on 2026-10-23: a local single-agent
review of `main...chem`, and the cloud `/code-review ultra` of
`food-cap-standard...chem`. Fixed:

- Every setter or stamp that kills a cell, replaces its occupant or rewrites its store now keeps the waste-mass balance `d(ΣW + Σw) = in − out`.
- Compounds `c >= 1` carry no energy.
- The Python view (`params()`, recipes, Restart) stays in agreement with what C runs.

There are 13 regression tests (22 cases) in `tests/test_chem_review_fixes.py`,
all through public entry points with constructed fixtures. Each test failed
with its fix reverted. Section 5 of `tests/notebooks/chem-m2-test.ipynb` holds
the matching PASS/FAIL cells. Dynamics of an ordinary food/waste or
single-food run are unchanged: the goldens pass and the Chem notebook outputs
are unchanged. The shared `diffuse()` refactor was left for a separate change.

| # | finding (L = local, C = cloud) | fix | test |
|---|---|---|---|
| 1 | L1 / C3: `stamp_patch` overwrite left a cleared cell's private waste in place, to vanish under the next birth | `evoca_chem_release(mask)` returns the host's `c >= 1` compounds to the site before the stamp; `evoca_chem_sync_from_f` also returns them for any cell marked dead | `test_stamp_overwrite_dead_cells_return_host_waste` |
| 2 | L5: a stamped organism inherited the host's waste | same release | `test_stamp_overwrite_live_cells_start_without_host_waste` |
| 3 | L4: `set_f_all` erased private waste at `N_chem > 1` | rewrites only the energy compound (`chem_cell_set_energy`); bit-exact at `N_chem = 1` | `test_set_f_all_keeps_private_waste` |
| 4 | L6: `set_alive*` deleted waste instead of returning it | the four setters settle dying cells through `chem_cell_vacate` | `test_set_alive_returns_dying_cells_waste` (4 cases) |
| 5 | L7: `set_w` (`evoca_set_chem_priv`) was not booked | a rise is booked as input, a fall as output | `test_set_w_is_booked_in_the_flux` |
| 6 | L2 / C2: weight 1 was accepted for `c >= 1` (and was the C default), making the tax free while booking it as dissipated | compound 0 is fixed at 1 and `c >= 1` at 0, the default in `evoca_init` and `chem_env_reset`; other requests are clamped with a notice | `test_set_chem_energy_refuses_an_energy_carrying_waste`, `test_low_level_metabolism_defaults_to_energy_free_product` |
| 7 | L3 / C4: `set_n_chem` left `wdiff`/`waste_sink`/`waste_disp` and the preset stale | Python `set_n_chem` configures as `init(n_chem=n)` would: 1 resets the waste params, 2 applies the preset, other sizes are refused. The inert `N_chem=2` parity probe calls the C entry point `_lib.evoca_set_n_chem` | 4 `test_set_n_chem_*` tests, including a record → replay fingerprint check |
| 8 | C5: `init` raised for a bad chemistry after it had changed `N` and the metaparams | `_check_chem` runs first (and refuses non-finite waste params) | `test_refused_init_leaves_the_instance_untouched` (3 cases) |
| 9 | C1: Restart (`controls.on_restart`) re-ran `evoca_init` (which resets to `N_chem = 1`) plus a fixed list of setters, dropping the chemistry | `EvoCA.init` = validate + store, then `_apply_metaparams()`, the single place C is configured from; `EvoCA.reinit()` and `controls.restart_sim()` reuse it | `test_restart_reapplies_every_metaparameter` (restart equals a fresh run) |
| 10 | L8: colormode 6 showed no waste when metabolism was off | shows the metabolic product, or compound 1 | `test_food_waste_colormode_shows_seeded_waste_without_metabolism` |
| 11 | L9: `params_from_row` passed NaN cells, which defeated the truthiness check | NaN cells dropped at load time (defaults apply) | `test_params_from_row_drops_nan_cells` |
| 12 | L12: `get_W`/`get_w` hard-coded `n_chem == 2` and read a NULL field when absent | `get_chem_env(c)` / `get_chem_priv(c)` (zeros when absent); `get_W`/`get_w` are the compound-1 conveniences; the gallery's `W_mean` keys on C's `N_chem` | `test_per_compound_accessors` |
| 13 | L11: the env source was booked once per cell in the hot loop | booked once as `inc × cells` | `test_env_source_is_booked_as_inc_times_cells` (pins the value; behavior-preserving, so it cannot fail with the change reverted) |
| 14 | C7: stale `evoca_set_n_chem` header comment | rewritten (also the `CHEM_MAX_N` and energy-weight comments) | — |

**Skipped:** L10 / C6, the third copy of the 3×3 periodic blur. The shared
`diffuse()` is a separate refactor.

**Energy across a stamp** is not claimed. `stamp_patch` writes `f` directly (an
experimenter injection). Only the waste mass is balanced across it.

**Structural notes:**

- `_apply_metaparams` turns the Restart re-initialization from a hand-maintained setter list into a derived one: a new metaparameter now reaches Restart by being added to `_apply_metaparams`.
- `evoca_chem_release` and `chem_cell_vacate` add outside-the-step call sites of `chem_return_to_env` (it was called only from tax death, toxic death and eviction).


## 9. Display follow-ups (2026-10-23, N)

- **Colormodes.** Colormode 6 (food-waste) is now built like env-food (1): a dead cell shows the environmental fields, food F in green and waste W in blue; an alive cell shows red 180 plus its private stores, food f in green and waste w in blue. The earlier encoding was red = W, green = F, blue = alive. New colormode 7 (env-waste) is the waste analog of env-food: waste only, in blue. Tests: `test_food_waste_colormode_channels` and `test_env_waste_colormode_channels`, both with constructed fields; together with the colormode count check they fail on the previous build.
- **Run/Pause label.** `Gallery.show` opened running with the button reading "Run". The front end snaps a slider value to its step (`m_scale` 3.5 becomes 3.5000000000000004) and sends it back. The slider callback took that as a move and paused, and its 200 ms timer then resumed the sim without updating the label. Now a change below float rounding is ignored, and the resume goes through `_set_paused(False)`. `tests/test_controls_pause_label.py` drives `run_with_controls` headless; 2 of its 3 tests fail without the fix.
- `Chem.ipynb` was re-executed; its text outputs are unchanged, and the film is in the new colors.
