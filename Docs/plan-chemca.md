# chemca: plan and implementation summary

For N, 2026-10-24. Branch `chemca` (not merged). Spec: `Docs/2026-10-23_chem-model-v2.md`
(status section, §4-10) and the coordinator's brief, with three additions N approved during
the work: a mass reproduction gate, an isologous-diversification readout, and the gs3 rates
from the closed-loop Gray-Scott study (`Research/Chem/Experiments/2026-10-24_gs_closed_regimes`,
on main).

**Summary.**
- One library now holds both models: `gcc -O2 -Wall -fPIC -dynamiclib -o C/libevoca.dylib C/evoca.c C/core.c C/chemca.c` (no warnings).
- The shared core extraction is behavior-neutral. Trajectory hashes are identical to main, and the 250 existing tests pass unchanged.
- The old goldens never covered diffusion ≥ 1, so a full-pass mutation passed the whole suite. New goldens close that gap.
- chemca is a separate engine with its own PRNG state. It keeps exact mass bookkeeping (closed mode) and an energy budget that closes to double precision.
- ChemCA (Python) comes with presets fw2 and gs3. Recipes carry `model`, and a recipe without it loads as EvoCA.
- `build_sim`, `make_sim`, `import_run` and `load_recipe` dispatch on `model`. `run_chemca` gives the summary metrics.
- The full suite has 298 tests and all pass. 31 mutations (21 C, 10 Python) and 3 core mutations were each caught.
- `tests/notebooks/chemca-test.ipynb` runs headless: 14/14 PASS.

**Dynamics at the provisional presets (N = 128-256).**
- Both presets persist with genes frozen and with mutation on (mu 0.02).
- **fw2:** patchy, toxicity-limited colonies at about 0.37 alive, rearranging slowly. No oscillation or waves.
- **gs3 at the exact S3 rates:** persists at about 0.2 alive, but irregular and patchy. No clean waves, no measurable site period.
- **gs3 with both rates x10** (same a = k1/k2, light unchanged): clean traveling excitation wave trains with a site period of about 450 steps and a global period of about 490, at N = 128 and 256.
- Spirals are not established at either setting.

---

## 1. Plan followed

| Step | Content | Status |
|---|---|---|
| 0 | `C/core.c`: PRNG step, Poisson, one box-blur diffusion (CHEM-P1); evoca.c calls them; build stated once | done, commit c41e4e0 |
| 1 | `C/chemca.c/.h`: state, genes, step, bookkeeping, API, colorize | done, 1dc4869 |
| 2 | `python/chemca_py.py`, presets, recipes, dispatch, `run_chemca`, headless frames | done, 35df159 |
| 3 | tests `tests/test_chemca_*.py`, `tests/test_core_diffuse.py`; mutation checks | done |
| 4 | `tests/notebooks/chemca-test.ipynb`, this document | done |
| + | repro_gate 'mass'; `composition_clusters`; k grid widened to 6 bits; gs3 = S3 | done (N's additions) |

## 2. What was built

### Step 0: shared core (`C/core.c`, `C/core.h`)

- `core_xorshift32(state*)` (inline) and `core_poisson(next, lambda)`. The stream is passed in as a step function, so each model keeps its own PRNG state. EvoCA's dynamics and shadow streams are unchanged.
- `core_diffuse(field**, tmp**, N, d)`: the one copy of the 3x3 periodic box blur. It replaces the food full pass, the food partial pass and the compound-field pass in evoca.c.
  - A full pass is the partial rule at s = 1. There t = 0 and w9 = 1/9, so the result equals sum/9 bit for bit.
- `core_diffuse_precise`: the same stencil, accumulated in double with the unbiased update F + s(mean - F). See §5, deviation 1.
- `core_wrap(i, N)`. Moore-neighbor selection was not extracted: the two models rank neighbors differently (EvoCA: f with a 1e-7 tolerance; chemca: e or body mass, dead = -inf, exact ties).
- `python/evoca_build.py` is the one Python statement of the sources and the gcc command. `tests/conftest.py`, the five top-level notebooks' build cells and three test notebooks call it. The Makefile and CLAUDE.md carry the same command, and `test_build_statement_is_consistent` asserts that they agree.
- **Makefile and zsh:** on N's machine `make` is an autoloaded zsh function ("function definition file not found"). That is a shell setup issue, not the Makefile. `command make` or `/usr/bin/make` works. The note is in the Makefile and in CLAUDE.md.

### Step 1: `C/chemca.c`

- **State:**
  - dense stores `[N²·K]` (K ≤ 8); see the comment there on departing from pinned decision A;
  - `e`, genes `[N²·8·5]`, the active mask, the env fields `E[K]`;
  - per compound: G, diff, u, tox, inc, sink;
  - up to 8 environment reactions (scalar rate or a per-site field);
  - its own PRNG state.
- **Gene:** 5 bytes `{A, B, C, m, k_idx}`. The rate is `k = k_min·2^(k_idx/4)` with **k_idx 0..63** (§5, deviation 2).
- **Step order** (as in the brief):
  1. environment reactions, then the open-mode source and sink;
  2. diffusion (precise);
  3. two-way transport;
  4. reactions;
  5. costs;
  6. toxicity;
  7. reproduction and mutation.
- **Reactions:** the rates come from a snapshot of the stores. For each substrate, λ = Σ k_j c_j^m_j and the amount converted is a(1 - e^-λ), split in proportion λ_j/λ (order-independent, never overdraws). Downhill genes apply first; uphill genes then run in slot order while e lasts.
- **Deaths** (toxic, energy, eviction) return all stores to the site. e is lost and booked.
- **Energy booking:** the budget books the *realized* float change of e for yield and cost, so `d(Σe) = yield - cost - death - split` holds to double summation. Float rounding of e shows up only in the compound free-energy identity `dΦ = light - yield (+ open_in - open_out)`, at float precision.
- **Mutation:**
  1. active-bit pass: Poisson(mu_ractive·8) flips, rejecting an empty mask;
  2. duplication onto newly active slots with probability p_dup;
  3. point mutations: Poisson(mu_rgene·8), with field weights 0.7 / 0.15 / 0.15 for k_idx±1 (clamped), a new (C, m) pair, and A or B redrawn (rejecting A = B).
  - A birth is flagged "mutant" only if the genome actually changed.
- **Colormodes:**
  - 0: env RGB of compounds 0-2;
  - 1: alive red, with stores of compounds 0 and 1 as green and blue;
  - 2: energy heat map;
  - 3: k_idx of the first active gene.

### Step 2: Python (`python/chemca_py.py`, `python/sim_model.py`)

- `DEFAULTS` is the one source of defaults. Per-compound keys `<field>_<compound>` are derived from `DEFAULTS['compound']`, and `param_keys(preset)` and `defaults(preset)` are derived too; there are no parallel key lists. Unknown keys raise.
- `ChemCA`:
  - `init(N, preset=..., **params)`, `reinit`, `update(**kw)` (live), `state(...)`, `step(n)`, `set_seed`;
  - getters for alive, E, store, e and rgenes; `mass`, `mass_flux`, `free_energy`, `energy_budget`, `events_last`, `state_hash`;
  - `k_idx_by_type`, `n_distinct_genesets`, `colorize`, `rgb`, `save_frames` (PNG, headless);
  - `recipe`, `export_recipe`, `set_light_field` (r(x)).
- `state()` options:
  - alive pattern: halfplane, fraction, all, or an array;
  - `mbar` with the preset's mass split; each alive cell moves `store_frac` of its site's E into its stores, so the mean mass per site is exactly `mbar`;
  - `e_init`, preset or random genes, `k_idx0` (an int or a per-gene list), `E_noise`;
  - `seed`, which is recorded in the recipe, so a recipe reproduces its run.
- Recipes: `version 1`, `"model": "chemca"`, init/final metaparams and the initialization. A per-site light field is not serialized; the recipe says whether one was set.
- `sim_model`: the `model` key, its default `'evoca'` (applied at load time, never written back) and `split_model`.
- Dispatch:
  - `evoca_explore.build_sim` (chemca → `build_chemca`);
  - `evoca_nb.make_sim`;
  - `evoca_py.import_run` (chemca → `chemca_py.load_recipe`);
  - `evoca_nb.load_recipe` (no EvoCA reseeding or v redraw for ChemCA; the recipe's own seed applies).
  - `run_sim` refuses chemca with a pointer to `run_chemca`. EvoCA paths are unchanged.
- `run_chemca(params, n_steps, N, seed, sample_every, init, kymo_steps=1500, comp_window=100)` returns:
  - alive mean, std and final, and extinct;
  - the per-step mass-identity max error (open mode: against in - out) and the per-step energy-identity max error;
  - the cumulative energy budget;
  - k_idx by gene type, and the distinct active gene-set count;
  - `evoca_waves.site_period` of the middle row (alive, and E of compound 0);
  - births and deaths by kind;
  - the composition-cluster readout;
  - samples, and the parameter echo.
  - It is re-exported from `evoca_explore`.
- `composition_clusters(sim, ...)`: numpy only. Features are the stores normalized by their total, plus window-averaged compositions if given. It runs k-means++ for k = 2..6 and picks k by mean silhouette. It returns 1 cluster when the spread is below 1e-3 or the best silhouette is below 0.5. Outputs are the number of clusters, the silhouette, the fractions and the centers.
- `repro_gate='energy' | 'mass'` and `mass_repro`. See §4.
- Not done (per the brief): no change to `controls.py`. Use `save_frames` or `rgb` to look at a run.

## 3. Provisional preset values

All are **provisional**. The gs3 rates come from the regime study; everything agent-side (m_e, c_g, u, tox, store_frac) was chosen so that populations persist at N = 128.

| | fw2 | gs3 |
|---|---|---|
| compounds, G | F 1, W 0 | F 1, A 0.5, P 0 |
| genes (slot: gene, k_idx → k) | 0: F→W, none, 29 → 0.152 (active); 1: F→W cat W m=1, 29 (inactive, reachable by mutation) | 0: F→A cat A m=2, 50 → 5.79 (S3 k1 = 6.0); 1: A→P, 27 → 0.108 (S3 k2 = 0.1); both active |
| light | W→F, r = 0.015 | P→F, r = 0.005 (S3) |
| diffusion | F 0.1, W 0.1 | F 0.6, A 0.3, P 0.3 (S3) |
| transport u | 0.2, 0.2 | 0.5, 0.5, 0.5 |
| toxicity | W 0.6 | P 0.6 (rarely binding at these values; tox_A is the option on A) |
| costs | m_e 0.008, c_g 0.002 | m_e 0.001, c_g 0.002 |
| k_min | 1e-3 | 1e-3 |
| state | mbar 1, split (0.9, 0.1), store_frac 0.5, e_init 0.5 | mbar 1, split (0.8, 0.1, 0.1) (A seeded: A = 0 is absorbing), store_frac 0.5, e_init 0.5 |

- **Cost split.** The total cost of each default genome (one active gene in fw2, two in gs3) equals the cost used while tuning, so frozen-gene dynamics are unchanged by the split. The gene cost c_g is needed once mutation is on (§6).

## 4. The mass reproduction gate (N's addition 1)

- **Body compounds** are those with transport rate u_c = 0.
- With `repro_gate='mass'`, a cell divides when its body mass reaches `mass_repro`. Eviction then targets the neighbor with the lowest body mass; dead neighbors are still preferred.
- Energy is still produced and spent (costs, uphill genes) but does not gate division.
- **Threshold chosen over a per-cell birth mass.** A fixed threshold is cleaner: no extra per-cell state, and the same shape as the energy gate. With the exact 50/50 split, a newborn holds about mass_repro / 2, so division happens at a doubling of birth mass. The constructed test shows it firing on the first step with body ≥ 2 x birth mass (0.394 → 0.401 against 0.4).
- **Body compounds per preset.** Neither preset has a body compound by default (all u > 0), so `repro_gate='mass'` with the defaults raises.
  - **gs3:** set `u_A = 0`; the body is A. This is the natural choice, since A is the autocatalyst that is inherited at division. Persists at mass_repro 0.4 (N = 128: alive about 0.45; e accumulates, to a mean of 3.4, because it is no longer spent on division).
  - **fw2:** the only candidate is W (F is food and must be taken up). With `u_W = 0` waste stays inside cells, so `tox_W` has to be raised or disabled. Not tuned; flagged as an open item.

## 5. Deviations from the brief and resolutions

1. **Diffusion: one stencil, two precisions.** The float box blur has a systematic mass bias: t + 9·w9 ≠ 1 in float, so partial passes lose about 2e-8 of mass per pass and full passes gain mass. That is invisible to EvoCA but breaks a closed-mass identity, at about 2e-5 relative over 1000 steps in chemca.
   - EvoCA must stay bit-exact, so `core_diffuse` keeps the float rule.
   - chemca uses `core_diffuse_precise`: the same stencil in one function, with a `precise` flag in double.
   - Residual closed-mass drift is about 1e-10 per step relative. It comes from float32 storage in the light reaction (3e-7 over 3000 steps at N = 128), and is tested at a 1e-6 per-step and 1e-5 per-3000-step tolerance.
2. **k grid: 6 bits** (k_idx 0..63, coordinator's instruction). The brief's 0..31 at k_min 1e-3 tops out at 0.2 and cannot reach the S3 rate k1 = 6. Range now 1e-3 to 55, still quarter-octave steps; k_min remains a metaparameter.
3. **`run_chemca` kymograph window is 1500 steps, not 600** (parameter `kymo_steps`). gs3 periods are about 450-500 steps, and `site_period` searches only up to T/2, so 600 steps cannot see them. The composition window is separate (`comp_window` = 100): in fw2 almost no cell lives 600 steps.
4. **Open-mode sink is proportional** (E·(1 - e^-sink), Gray-Scott-like washout) rather than EvoCA's uniform floored removal. The brief left the law open.
5. **Environment reactions with a shared substrate** apply in the order added, not split by rate (each has its own exponential integrator). Both presets have one.
6. **Duplication onto a newly active slot** copies a random *other* active slot (EvoCA's egene code can pick the new slot itself, a no-op).
7. **A or B redraw** is uniform over the compounds, rejecting A = B. In fw2 (K = 2) that redraw can only reproduce the current value, so F→W never becomes W→F by one mutation.
8. **Uphill genes** run after all downhill genes in the cell, in slot order. Only when e is binding does the slot order matter.
9. **`set_alive`** zeroes the stores, e and genes of cells set dead (an init-time configuration, not a death: nothing is returned or booked).
10. **Costs:** a cell that cannot pay its cost dies with e = 0, and what it had counts as spent. The death booking is 0. This keeps the e identity without negative "lost" energy.

## 6. Results

### Test suite

- 298 pass: 250 existing + 6 `test_core_diffuse` + 42 chemca.
- **Core extraction:**
  - trajectory hashes over 5 EvoCA configurations (gdiff 0.08, 1.0, 2.37, 3.0; one-way and two-way waste; wdiff 1.3 and 2.0; shadow on) are identical before and after;
  - `test_core_diffuse` holds goldens computed with the pre-extraction library (main @ a9854b4).

### Mutation checks (each guarded piece broken, guarding tests run, restored)

| Mutation | Caught by |
|---|---|
| core: full-pass weight 0.999 | 4 (both new goldens, direct full-pass test); the old suite missed it (250 pass) |
| core: compound diffusion clipped at 1 | 1 (wdiff golden) |
| core: Poisson L = e^-1.01λ | 6 (existing goldens) |
| M1 eviction does not return evictee stores | eviction fixture, all identity tests (5) |
| M2 toxic death does not return stores | toxic fixture, identities (4) |
| M3 energy death: remaining e not booked | energy-death fixture, identities (4) |
| M4 e lost at death not booked | toxic and eviction fixtures, identities (5) |
| M5 split 0.49 / 0.51 | split fixture |
| M6 shared substrate not split | sharing fixture, identities (4) |
| M7 uphill gene not clipped | uphill fixture |
| M8 m = 2 treated as m = 1 | conversion[A-2] |
| M9 linear instead of exponential integrator | 4 conversion fixtures |
| M10 export-only transport | transport[uptake] |
| M11 transport for dead cells | dead-site fixture |
| M12 A = B allowed | mutation reachability (2) |
| M13 no point mutations | mutation reachability (2) |
| M14 (C, m) branch off | mutation reachability (2) |
| M15 dead neighbors not preferred | dead-neighbor fixture |
| M16 mass gate ranks eviction by e | body-mass eviction fixture |
| M17 mass threshold 0.9x | doubling fixture |
| M18 light energy not booked | identities (3) |
| M19 cost booked nominal, not realized | identities (3) |
| M20 open-mode sink not booked | open-mode balance |
| M21 chemca back to float diffusion | identities, long-run drift (3) |
| P1 missing model → chemca | EvoCA-recipe test, build_sim, make_sim (3) |
| P2-P4 build_sim, make_sim, import_run without dispatch | each caught (P4 by 3 recipe tests) |
| P5 recipe drops seed | round-trip hash |
| P6 mu_rgene not pushed to C | build_sim (reads C), mutation reachability (3) |
| P7 repro_gate not pushed | 4 gate tests |
| P8 no silhouette floor | homogeneous → 1, scale invariance |
| P9 no normalization | scale invariance |
| P10 light not pushed | final-params recipe, make_sim (reads C) |

All 31 chemca mutations were caught (run twice, before and after the final preset change). The restored suite passes.

### Dynamics (test notebook, N = 128-256)

- **fw2 (genes frozen):** alive 0.37 ± 0.03 over 3 seeds, 4000 steps. Turnover is mostly toxic deaths (about 85k per 4000 steps against about 8k energy deaths). Colonies are patchy and rearrange slowly. `site_period` is nan.
- **gs3 at the S3 rates:** persists 6000 steps at N = 128 and 256, alive about 0.2 (std 0.02). The pattern is irregular and patchy.
  - Global alive period: 956 at N = 128 (ACF 0.36, weak); none at N = 256.
  - Site period: nan.
  - This is **not** the study's spiral regime: inside agents, the chemistry is coupled through transport (u 0.5), and the colony turnover changes the effective time scales.
- **gs3, both rates x10** (k_min 0.01, same k_idx; k1 = 58, k2 = 1.08, a = 54 unchanged): regular traveling excitation wave trains.
  - N = 128: site period 447, global period 499.
  - N = 256: site period 447, global period 488.
  - The kymographs show curved fronts and V-shaped source and collision structures.
  - This was found by accident: an intermediate draft of the preset still had k_min = 0.01.
  - Recommendation: the regime study's map should be redone with agent coupling, and the time-scale ratio between reactions and transport looks like the control parameter.
- **Spirals:** not established in either setting. Single frames show curved bands; a tip-tracking analysis was not done.

### Evolution (mu_rgene = mu_ractive = 0.02, N = 128)

- **Without a gene cost, fw2 collapses** by step 600-750, even at mu 0.005. The cause is the active-bit channel: activation and duplication give cells extra F→W copies, rates add, and those cells reproduce faster. That raises W until toxicity kills the population: runaway selection on gene copy number.
  - mu_rgene alone persists.
  - The doc's per-gene cost c_g (§8b) stops it: c_g 0.002 persists, while 0.005 and above kill through energy.
  - The presets now carry c_g 0.002, with m_e lowered so the default genome's cost is unchanged.
- **With costs, both persist** (alive about 0.22, with 350-500 distinct gene sets after 4000 steps).
  - fw2 evolves F→W catalyzed by **F** (m = 1, 2), i.e. a rate ∝ F² or F³, more than the W-autocatalytic variant.
  - gs3 evolves the shortcut **F→P catalyzed by A** (the full yield of 1 per unit in one gene), and k(A→P) drifts down (27 → 25.5).
  - First observations from 2 seeds, not results.

### Isologous diversification, genes frozen (N's addition 2)

- `run_chemca` reports **2 compositional clusters** in both presets:
  - fw2: silhouette 0.63-0.65, fractions about 0.78 / 0.22;
  - gs3: silhouette 0.51-0.52 (just above the 0.5 floor), fractions about 0.6 / 0.4.
- **But in fw2 the split is positional.** The minority cluster (more F) has 3.7 alive neighbors on average against 6.3, i.e. colony edge against interior.
- An earlier gs3 draft gave the same picture (4.5 against 7.0 neighbors).
- So this is environment-driven composition, not Kaneko-style cell-intrinsic differentiation. A test for the latter needs the types to persist in a cell and be inherited across divisions (track composition along lineages, or transplant cells), which is not implemented.

## 7. Known gaps and open issues

- **Preset values are provisional.** gs3 at S3 gives no clean waves inside agents. The agent-coupled regime map (k scale, u, m_e) is the next study.
- **fw2 with the mass gate** needs W as body (u_W = 0) and a raised or disabled tox_W. Not tuned.
- **Isologous diversification:** the readout fires on positional composition. A lineage-based persistence test is needed before calling anything differentiation.
- **No live display:** no `controls.py` integration (out of scope). `save_frames` and `rgb` only.
- **Not serialized in recipes:** the per-site light field r(x) and `alive` arrays (the recipe flags the first; loading an array-made recipe raises).
- **One instance per process:** the ChemCA state is global in C, as for EvoCA.
- **Spiral detection:** not done (no core tracking); the visual evidence is curved bands only.
- **Closed-mass drift:** about 1e-10 relative per step (float32 storage in the light reaction). Fine for 1e5 steps; a compensated-sum light reaction would remove it if needed.

## 8. Structural-risk points (flagged per the standing rule)

- **A new model in the shared library, and code moved across files** (core extraction: PRNG step, Poisson, diffusion moved out of evoca.c). Guarded by bit-exact hashes and the new full-pass goldens.
- **A new pluggable point:** the `model` key with its load-time default. There are four dispatch doors (`build_sim`, `make_sim`, `import_run`, `load_recipe`), all calling `sim_model`, each with a reachability test that asserts on what reached C.
- **A second call site for the diffusion** (chemca uses the precise variant of the same stencil).
- **Configured → derived:** body compounds are derived from u_c = 0 rather than configured.
- **Before merging, N may want `/code-review ultra`** (user-triggered and billed; I cannot launch it). The dispatch and the bookkeeping are where an independent reader would add the most.

## 9. Structural and scaling review

- **Artifacts:** each has one owner.
  - The build: `evoca_build.py`, asserted against the Makefile and CLAUDE.md.
  - Defaults: `DEFAULTS` + `PRESETS`.
  - The model key: `sim_model`.
  - The gene layout: `chemca.h`, checked at load by `ChemCA.__init__`.
- **Scale:**
  - `run_chemca` calls `mass()` and `e_total()` each step: O(N²K) in C. Measured at N = 128 (gs3): `mass()` 0.04 ms against 0.3 ms per step, so the checks are cheap. Its per-step `get_alive`/`get_E` copies inside the kymograph window cost more, but only in that window.
  - `n_distinct_genesets` is a Python loop over alive cells (4 ms at N = 128); it runs at the end only.
  - `composition_clusters` subsamples to 1500 cells (an O(n²) silhouette).
- **Opportunity (not done):** a C-side `chemca_step_checked` returning the per-step mass and energy residuals would make the identities free in scans.
- **Opportunity (not done):** `controls.py` → sim methods plus a slider registry, which the doc (§9c) already plans. It is the next pluggable point, and the one most worth `/code-review ultra`.
