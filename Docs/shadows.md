# Neutral shadows in EvoCA

A reference for every neutral shadow in the code, and the design of the stateful
**drift shadow** that replaces the memoryless fixed-space construction.

Written for N, 2026-10-15. Code references are to `main` at `a0798d1` unless marked
`shadow-rng` (branch `shadow-rng`, commit `12f4775`, based on `530c334`). Line numbers
drift. The function names are the stable anchors. Claims are tagged **[code]** (read in
the source), **[run]** (checked by a run this pass, scripts in section 12), **[recorded]**
(in a committed document or CSV), or **[reasoned]** (argument only, not yet tested).

Working names (N may rename; section 9): the existing fixed-space construction is the
**tracking shadow** (`fs_shadow='tracking'`). The new one is the **drift shadow**
(`fs_shadow='drift'`). The LUT-hash shadow is the **Channon shadow**.

---

## 0. Summary

**What a shadow is for.** Bedau–Packard activity statistics count how much each heritable
component is present over time. A *neutral shadow* is the same population, with the same
births, deaths and mutation supply, reproducing *fitness-blind*. Excess activity is real
minus shadow. It should be about zero, in distribution, whenever there is no selection,
whatever the mutation flux.

**What exists.**

| Shadow | Space | State | Valid when | Main defect |
|---|---|---|---|---|
| Channon (`neutral_*`) | open: FNV hash of the LUT | Np individual LUTs | LUT-only runs at low `mu_lut` | Real side hashes LUT plus egenes, shadow LUT only. Genotypes are near-singletons at productive `mu_lut`. On `main` it shares the dynamics RNG |
| tracking (`fs_*`), eg and dyn | 729 ternary egene keys, 500 (input, output) transitions | none carried between steps | never as a selection test | The excess reduces to a concentration index. The frozen null passes by construction. The z-scores count buckets above uniform share |

**The drift shadow, in one paragraph.** Keep R replicate shadow *populations* of
individuals. Each individual carries a full genome: a 250-bit LUT and an 8-slot egenome.
Each step every replicate receives exactly the real run's death count and birth count:
uniform-random victims, then uniform-random parents among the survivors. Each shadow
child replays one real child's *realized* flip counts (`nf`, `na`, `ne`), at positions
drawn from a shadow RNG stream. Nothing selects. The shadow population is then projected
into each fixed space:

- **eg:** the ternary key of every active egene, exactly as `fs_eg_observe` projects the
  real run.
- **dyn:** per input `i`, the shadow's allele frequency `p_i` of output 1, weighted by
  the real run's observed query count `n_i(t)`. The real side is scored the same way,
  from the *population* allele frequency `P_i(t)` over all alive cells, not from the
  queriers.

The shadow's mutation kernel is therefore the real one by construction, so nothing has to
be derived and then kept in sync. Linkage between the 250 loci and the egenome is
preserved. The replicates give a Monte Carlo null for any statistic.

**What changes in the null tests.** A frozen run no longer reads zero by construction. If
the frozen population has standing variation, selection on it is real adaptive activity,
and a positive reading is *correct*. The right null is "no selection", not "no mutation".
It is checked three ways:

1. a synthetic neutral population (exact calibration);
2. an EvoCA **frozen-clonal + neutral-marker** run (real spatial demography, nothing to
   select);
3. the same marker inside evolving runs, which measures the hitchhiking (draft) floor.

**Two findings that shape the design.**

- **[run]** At EvoCA's operating point (N=256, `mu_lut`=0.01) a locus sees
  θ = 2·Np·u ≈ 300 mutations per generation. Panmictic neutral drift is then tiny, and the
  replicate spread is narrow. So model misspecification (space, draft), not sampling
  noise, sets the real floor. Effect sizes against the marker null matter more than
  replicate z-scores.
- **[run]** In fixed spaces with θ ≫ 1 every component is present all the time
  (D = 79.5 of 80 in the stylized test). The Bedau–Packard mean activity per present
  component, Ā = ΣA/D, therefore cannot discriminate. The drift shadow must be read
  through per-component or distributional statistics (section 6).

**Recommended name:** keep **drift shadow** (section 9).

**Open questions for N** are in section 8. The main ones:

- the dyn real side: population allele frequencies (recommended) or queried transitions;
- a separate RNG sub-stream rather than literally sharing the Channon stream
  (recommended);
- R = 8 replicates by default;
- emit tracking and drift together (`fs_shadow='both'` default in `run_sim`);
- add a neutral marker trait for calibration.

---

## 1. Purpose: what a neutral shadow is, and what "excess" should mean

### 1.1 The construction in the literature

- **Bedau and Packard** (1992, *Artificial Life II*) defined evolutionary activity: a
  component's activity is its cumulative usage or presence. **Bedau, Snyder and Packard**
  (1998, *Artificial Life VI*, "A classification of long-term evolutionary dynamics")
  introduced the neutral shadow as the calibration that separates adaptive persistence
  from persistence that happens anyway. They used it to classify dynamics by diversity
  `D`, mean cumulative activity `Ācum`, and new activity `Anew`.
- **Channon** (2001, ECAL, "Passing the ALife test", `Docs/references/channon_ad_ecal2001.pdf`)
  states the operational recipe. The shadow "mirrors the real run in every detail except
  that selection is random. Whenever a real organism is killed, a randomly chosen shadow
  organism is also killed. Whenever a real organism is born … a new shadow organism is
  born as the product of two randomly chosen shadow organisms, using the same
  reproduction procedure with the same rates of crossover and mutation." Two further
  details from the same paper:
  - `Anew(t) = (1/D(t)) · Σ_{i: a_i(t) ∈ [a0, a1]} a_i(t)`: activity in a *band* chosen
    from the real and shadow activity distributions, normalized by diversity. It is not a
    count above a percentile.
  - "isolated occurrences" (a component present in this snapshot but not the previous
    one) are screened out on both sides.
- **Channon** (2006, *GP&EM* 7:253–281) applies the same machinery to a larger system.
  I am not confident of its normalization details and do not cite them here.

### 1.2 The estimand

Let `X_t` be the real run's state and `Y_t` a component-frequency readout (a bucket
count vector). Let `S` be the shadow process: given the realized demography and mutation
supply of `X`, it produces `Y^S_t` with parent choice independent of genotype.

> **Excess** = a functional of the cumulative activity `A(t) = Σ_τ Y_τ`, contrasted
> between real and shadow, such that under H0 (genotype does not affect reproduction or
> survival) the real `A` is a draw from the same distribution as the shadow's `A`.

Three consequences decide most of the design.

1. **What is conditioned on must be stated.** The shadow conditions on the real run's
   realized population size, births, deaths and mutation counts, and, for dyn, on the
   real query mass `n_i(t)`. Any selection that acts *through* those quantities is
   absorbed into the null and is not detected. Examples: selection that raises the
   population size, or makes some inputs more frequent. The test is therefore
   conservative about those routes. Section 5.9 lists them.
2. **Zero flux is not the null; zero selection is.** A no-mutation run with standing
   variation and selection *should* read positive, because selection on standing
   variation is adaptive activity. The tracking shadow reads exactly zero there (section
   4). That is a defect, not a validation.
3. **Excess must be population-comparable.** Cumulative activity scales with the mass
   per step `T` (alive cells for dyn, active egenes for eg). A ranking objective needs an
   intensive form (section 6). Guide 07 §2.1 showed that the current excess is extensive
   in `T`.

### 1.3 Properties a valid excess must have (the acceptance targets)

| # | Property | Test (section 7) |
|---|---|---|
| P1 | Mean about 0 and calibrated tails under no selection, at any flux, including skewed component accessibility | Synthetic neutral; frozen-clonal marker |
| P2 | Positive under selection, with known power | Synthetic imposed selection |
| P3 | Intensive (comparable across `T`, `N`, run length) | Rescaling test |
| P4 | Does not perturb the real run | Bit-identical trajectory, drift shadow on vs off |
| P5 | Reproducible per seed | Two runs, same seed, identical shadow outputs |

---

## 2. Anatomy of one step: where demography, flux and buckets come from

All shadows read counters filled inside `evoca_step` (`C/evoca.c` ~1622). Order matters.

| Phase | What happens | Counters touched | Where |
|---|---|---|---|
| top of step | counters reset | `g_births_last`, `g_deaths_last`, `g_eg_bitflips_last`, `g_lut_bitflips_last`, `g_repro_last` = 0 | ~1626–1632 |
| 1 (CA update) | each cell computes `bit = compute_lut_bit(row,col)` (the input index `i = v_x·125 + n1·25 + n2·5 + n3`) and `new_state = LUT_x[i]`. For **alive** cells: `dyn_act[2i+new_state]++`, `dyn_pop[2i+new_state]++` | `dyn_pop` (per step), `dyn_act` (cumulative since init) | ~1641–1665 |
| 2c (tax) | `f -= tax`. If `f ≤ 0`: dies, LUT zeroed, `active=0` | `g_deaths_last++` (tax death) | ~1709–1785 |
| 4 (reproduction) | each alive cell with `f ≥ 1` copies its genome into the min-`f` Moore neighbor. An alive occupant is evicted. Child LUT gets `nf ~ Poisson(mu_lut·250)` flips (or `·n_active` under `restricted_mu`). Egenome gets `na ~ Poisson(mu_egenome·8)` presence toggles (Negene→0 rejected), dup-on-activate with `p_dup_egene`, then `ne ~ Poisson(mu_egene·96)` value/mask flips over all 8 slots × 12 bits | `g_deaths_last++` (eviction), `g_births_last++`, `g_lut_bitflips_last += nf`, `g_eg_bitflips_last += na+ne`, `g_repro_last++` | ~1838–2003; flux at ~1991–2000 |
| end | `neutral_apply_demography(g_births_last, g_deaths_last)` (Channon). Then `fs_shadow_step(alive_pop)` (tracking) | | ~2043–2055 |

### 2.1 The two fixed bucket spaces, precisely

**dyn bucket** (`DYN_ACT_COUNT = 500`, ~474) **[code]**

- `key = 2·i + o`, with `i ∈ [0,250)` the LUT input index the cell queried in Phase 1
  and `o = LUT_x[i]` the output it produced.
- Every alive cell adds exactly one unit per step, so `T_dyn = Σ_k dyn_pop[k]` is the
  number of cells alive *at Phase 1*.
- **The `i` part is ecological.** It is set by the cell's own `v` and its 12 ring
  neighbors, that is, by the v-field. **Only `o` is heritable**: it is bit `i` of the
  genome.
- Mutation cannot move a cell from input `i` to input `j`. It can only flip which output
  a carrier of input `i` produces.
- So "the dyn composition" mixes a heritable allele frequency with an ecological input
  distribution `n_i(t) = dyn_pop[2i] + dyn_pop[2i+1]`.

**eg bucket** (`EGENE_KEY_COUNT = 729`, `egene_ternary_key` ~421) **[code]**

- `key = Σ_j d_j·3^j` over the 6 D4 orbits, with digit `d_j = 0` if the mask bit is 0
  (wildcard), 1 if mask=1 and value=0, 2 if mask=1 and value=1.
- Each **active** egene of each alive cell adds one unit (`fs_eg_observe` ~3537, run at
  the *end* of the step). So `T_eg = Σ_cells Negene ≈ 1.1·Np` at the causal-control base
  **[run]**.
- The key is heritable. But the genome under it is richer than the key, in two ways:
  - A **value bit under a wildcard is hidden state.** Flipping it does not change the
    key. A later mask flip at that orbit reveals it, as digit 1 or 2 depending on the
    hidden bit.
  - **Inactive slots** (pseudogenes) mutate too. They enter the key multiset only when a
    presence toggle activates them, possibly overwritten by a duplicate of an active slot
    (`p_dup_egene`).
- So the 729-key composition is a *projection* of a Markov process on a larger state.
  The key alone is not Markov.

### 2.2 Flux counters

- `g_lut_bitflips_last` sums `nf` over this step's births.
- `g_eg_bitflips_last` sums `na + ne`. That includes `ne` flips that land on inactive
  slots or on wildcard value bits, both of which leave every active key unchanged, and
  `na` toggles, which change *which* keys a cell carries, not any key's identity.
- **Measured at the causal-control JOINT base** **[run]** (`food_inc=0.013, m_scale=1.2,
  gdiff=0.06, tax=0.035, tax_per_egene=8e-4, mu_lut=0.01, mu_egene=0.003,
  mu_egenome=0.005`, steps 1000–3000):

  | N | alive Np | births/step | deaths/step | LUT flips/birth | eg flips/birth | Negene | ms/step |
  |---|---|---|---|---|---|---|---|
  | 128 | 3,924 | 45.4 | 44.7 | 2.49 | 0.327 | 1.05 | 0.34 |
  | 256 | 15,388 | 181 | 178 | 2.50 | 0.326 | 1.11 | 1.54 |

  Births per step are about 1.2% of Np, so one generation (Np births) is about 85 steps.
  A 15,000-step run is about 180 generations.

---

## 3. Catalog of shadows and their consumers

### 3.1 The Channon shadow (LUT-hash, open space)

| Item | Detail |
|---|---|
| Space | Open: `lut_hash_fn(LUT)`, an FNV-1a hash of the 32 LUT bytes. The real side (`act_*`, `get_activity`) hashes LUT ‖ sorted active egenes (`lut_hash_cache`) |
| State | `neut_luts[neut_n·32]`, `neut_hashes[]`, plus the hash table `nact_*` (activity, pop_count) (`C/evoca.c` ~2619–2700) |
| Seeding | `evoca_neutral_enable` (~2824): copy every alive cell's LUT 1:1. `run_sim` calls it at t=0 when `shadow=True` (`python/evoca_explore.py` ~196) |
| Update | `neutral_apply_demography(B, D)` (~2773), every step after Phase 4. (1) `D` uniform-random victims, swap-remove. (2) `B` births: parent `p ~ U{0..neut_n−1}`, copy LUT, `neut_mutate_inplace` (~2738), which is the same operator as the real run: `nf ~ Poisson(gmu_lut·250)` uniform positions, or `Poisson(gmu_lut·n_active)` over the real run's `active_bits` under `restricted_mu` |
| Borrowed from real run | `B`, `D` (counts only); nominal `gmu_lut`; under `restricted_mu`, the real run's current `active_bits` set |
| RNG | **`main`:** the dynamics stream `rng_next()`. Enabling it changes the real trajectory. **`shadow-rng`:** own stream `shadow_rng_next()`, seeded `shadow_seed_from(seed)` (a splitmix64 finalizer of seed XOR "SHADWRNG"), with `poisson_sample_from(next, λ)`. Shadow on and off then give bit-identical lattices (`tests/test_shadow_rng.py` on that branch) |
| Activity update | `evoca_n_activity_update` (~2897): `activity += pop_count` for live buckets. In `run_sim` it is called **once per sample**, like the real `evoca_activity_update`, so both sides count in samples, not steps |
| Accessors | `sim.neutral_enable/disable/is_enabled/neutral_population`, `n_activity_update`, `get_n_activity` (`python/evoca_py.py` ~1817–1850); `evoca_get_n_p90` (p90 of present shadow activities, cached by the strip renderer ~2960) |
| `run_sim` keys | `excess_activity_final/_slope` (ΣG − ΣN over present components); `excess_pc_final/_slope` (ΣG/D_G − ΣN/D_N); `a_new_final/_slope`, `a_new_thresh_final` (section 3.3) |

**Known pathologies**

1. **Scope mismatch.** The real hash covers LUT plus egenes, the shadow LUT only. Excess
   is undefined whenever egenes vary. Evidence: memory note
   `shadow-models-lut-only-mismatch`; v1 EGENE_only `excess_pc_slope` ≈ −26; guide 07 §3.
2. **Near-singleton genotype space at productive `mu_lut`.** P(no LUT flip in a birth) =
   e^(−250·`mu_lut`) ≈ 8% at 0.01. Both sides mint a new hash at almost every birth, so
   genotype activity measures persistence, not adaptation. Guide 07 §2.3: `A_new`
   threshold of 2–5 units; eff-species keys stay at 86–100% of full species.

   Research-directions §2 calls the same-operator shadow a "no-heredity model" and
   proposes "copy without mutation, inject at the realized novelty rate". I read it
   differently. The same-operator shadow *is* the textbook construction: Channon 2001
   uses "the same rates of crossover and mutation". What fails at high μ is the
   **component space**, which has no persistent units, not the shadow. A novelty-injection
   shadow would have heredity that the real run lacks, and so would not be neutral for
   this population either.
3. **Dynamics-RNG coupling on `main`.** Fixed on `shadow-rng`, not yet merged. Until
   then, `shadow=True` and `shadow=False` runs with the same seed are different
   trajectories.
4. **Minor demographic deviations** **[code]**:
   - Children born earlier in the same step can be chosen as parents (`neut_n` grows
     inside the birth loop). In EvoCA a newborn has f ≈ 0.5 and cannot reproduce in its
     birth step.
   - Victims are drawn from the pre-birth population. Real evictions can hit cells born
     earlier in the same step.
   - When the shadow empties, it reseeds with a *random* LUT.

   All three are small at Np ≫ B, but they are deviations from "mirror every detail".
5. **Mutation is nominal, not realized.** It is equal in distribution today. It would
   become wrong if per-cell `mu_lut_gene` (merged A2 scaffold, ~222, currently read
   nowhere) is ever made to drive mutation: the shadow would still use the global
   `gmu_lut`.
6. **No isolated-occurrence screen** (Channon 2001 §2). Singletons enter both sides, and
   that inflates the low end of the shadow activity distribution that `A_new`'s p90 reads.

**Valid when.** LUT-only evolution (`mu_egene = mu_egenome = 0`, egenome uniform) at low
`mu_lut`, where genotypes persist (roughly 250·`mu_lut` ≲ 0.5), with the `shadow-rng` fix
merged, read through `excess_pc` or a properly banded `A_new`. Elsewhere use it as a
churn diagnostic only.

### 3.2 The tracking shadow (fixed-space, memoryless), eg and dyn

| Item | Detail |
|---|---|
| Spaces | eg (729 keys), dyn (500 transitions), section 2.1 |
| State | `eg_q[729]`, `dyn_q[500]` (overwritten every step, so not state in any dynamical sense); `*_nact` (cumulative shadow activity, `uint64` plus `*_nact_frac` carry); `*_nact_var` (accumulated variance); `fs_eg_obs_pop/act` (private real-side eg mirror) (~3471–3512) |
| Seeding | none needed; "primed" is bookkeeping only (~3614) |
| Update | `fs_shadow_step` (~3649), every step: `rate = flips/(births·bits)` with `bits = 12` (eg) or `1` (dyn), then `fs_space_step` (~3564) per space, below |
| Borrowed from real run | **this step's observed composition** `obs_f`, mass `T`, flip totals, birth count, end-of-step `alive_pop` as "Np" |
| RNG | none (deterministic) |
| Accessors | `sim.eg_excess_pc()`, `sim.dyn_excess_pc()`, `sim.eg_shadow_total()`, `sim.dyn_shadow_total()`, `sim.get_eg_zscores()`, `sim.get_dyn_zscores()`, `sim.get_realised_flux()` (`python/evoca_py.py` ~1732–1813); C `evoca_eg/dyn_excess_pc` (~3703–3713), `evoca_eg/dyn_zscores` (~3776–3786), `evoca_eg/dyn_nbuckets` |
| `run_sim` keys | `eg_excess_pc_final/_slope`, `dyn_excess_pc_final/_slope` (sampled every `sample_every`, slope by `np.polyfit`; `evoca_explore.py` ~284–295, ~422–447). The z-scores are **not** emitted by `run_sim`; they are used in notebooks (`tests/notebooks/metric-zscore-test.ipynb`, Activity-metrics §4b) |

**Update rule, from the code** (`fs_space_step`, per bucket `k`, per step) **[code]**:

```
T      = Σ_k obs_pop[k]                     (mass this step)
occ    = #{k : obs_pop[k] > 0}
m      = 1 − exp(−rate)                      (rate = flips / (births · bits))
q_k    = (1 − m)·obs_pop[k]/T + m/K          (uniform smear over all K)
q_k   += (m/Np)·([obs_pop[k] > 0]/occ − q_k) (mean-field "drift" nudge)
q_k   ← q_k / Σ q                            (renormalize)
nact_k += q_k·T                              (with fractional carry)
nvar_k += m·T²·q_k(1 − q_k)/Np
```

`fs_excess_pc` (~3687):
`ΣG/D_G − ΣN/D_N`, where

- `ΣG` and `D_G` are taken over buckets with `obs_pop > 0` (observed cumulative activity);
- `ΣN` and `D_N` are taken over buckets with `q > 1e−12` and `nact > 0`.

`fs_zscores` (~3761): `z_k = (obs_act_k − nact_k)/sqrt(nvar_k)`, or 0 if `nvar_k ≤ 1e−12`.

**Known pathologies** (derivations in section 4)

1. Memoryless. The excess reduces to a concentration index `A_occ/D_G − M/K`, extensive in
   `T` (guide 07 §2.1, verified there to within 0.6 counts at the final sample).
2. The frozen null passes by construction (§4.2). Evidence: `tests/test_fixed_space_shadows.py`
   `test_neutral_drift_excess_is_zero` and `tests/test_zscores.py`
   `test_null_frozen_zscores_are_zero` both use zero mutation, so they test the code path,
   not the statistic. The paired "selection" tests contrast a zero-flux arm with a
   flux-plus-selection arm, so they confound flux with selection.
3. The z-tail count is about the number of buckets above uniform share. In a stylized
   *neutral* population, 113–173 of 500 buckets reach z>3 (guide 07 §2.2).
4. The baselines are wrong even as smears:
   - dyn: uniform 1/500 ignores input frequencies, and a mutation cannot even move mass
     between inputs.
   - eg: the neutral stationary key distribution is not uniform (§4.3).
5. The rates are not per-component rates applied per birth. `rate` for dyn is flips per
   birth over 1 bit, about 250·`mu_lut` = 2.5, so `m` ≈ 0.92. `m` is then applied to the
   **whole composition every step**, as if every individual were a newborn. For eg, `bits
   = 12` divides a count that includes presence toggles and silent flips. This affects
   only `nvar`, since the mean excess does not depend on `m` once `m > 0`.
6. Timing mismatch:
   - dyn's real side is counted in Phase 1 (cells alive at the start of the step);
   - eg's real side is counted at the end of the step;
   - `Np` is the end-of-step alive count.
7. **The frozen-rich-substrate artifact** (memory note `dyn-shadow-frozen-rich-substrate`;
   #3c arm A `dyn` +335 vs ≈+1). This is a hypothesis, unverified. Before a mid-run
   freeze, every bucket accumulates smear mass, so `D_N = K`. After the freeze `m = 0`
   and `q = obs_f`, so `D_N` drops to the occupied set, and the comparison becomes `A_occ`
   against `N_occ`. `N_occ` is the shadow mass that landed on now-occupied buckets, which
   is less than observed by the pre-freeze smear deficit. So the excess jumps to that
   frozen-in deficit per occupied bucket.

**Valid when.** As a *descriptor of concentration* of the bucket distribution (rename it
that). Never as a selection test.

### 3.3 Consumers of shadows

| Consumer | Reads | Definition in code | Status |
|---|---|---|---|
| `excess_activity_*` | Channon | ΣG − ΣN over present components, per sample (`run_sim` ~363–381) | Legacy; inverts at high `mu_lut`; scope bug |
| `excess_pc_*` | Channon | ΣG/D_G − ΣN/D_N (~383–401) | Scope bug; about 0 everywhere in campaign #1 |
| `a_new`, `a_new_thresh` | Channon | number of present real genotypes with cumulative activity > p90 of present shadow activities (~249–282) | Not the BSP/Channon statistic. Channon's `Anew` sums activities inside a band `[a0,a1]` and divides by `D`. Here the threshold collapses to 2–5 units in evolving runs, and it is blind (0) when the LUT is frozen (guide 07 §2.3) |
| `eg_excess_pc_*`, `dyn_excess_pc_*` | tracking | §3.2 | Concentration index (§4) |
| `get_eg/dyn_zscores` | tracking | §3.2 | Share-above-uniform counts; uncalibrated |
| `dyn_report` verdicts | tracking `dyn_excess_pc` plus `dyn_eff` (participation ratio) (`python/dyn_analysis.py`) | reads low `dyn_eff` as "silent churn" | Structurally anti-aligned with `dyn_excess` (guide 07 §2.5). The drift shadow gives `dyn_eff` a neutral reference (§6.4), which settles the question "is concentration selection or silence?" |
| strip charts (`evoca_n_activity_render_col`) | Channon | p90 line | Display only |

---

## 4. The tracking shadow: why its excess is a concentration index

### 4.1 Reduction (guide 07 §2.1)

Fix a space of size `K`. Write `o_k(t)` for the observed bucket count, `T(t) = Σ_k o_k(t)`,
and `M(t) = Σ_{τ≤t} T(τ)`.

1. For `m > 0`, every bucket gets `q_k ≥ m/K > 0` (up to the tiny drift nudge). After a
   few steps every bucket has `nact_k ≥ 1`, so **`D_N = K`**.
2. `Σ_k nact_k(t) = Σ_τ Σ_k q_k(τ)·T(τ) = M(t)` (to rounding; `dyn_shadow_total/M` was
   verified at 1.000 in guide 07). The shadow's `ΣN` in `fs_excess_pc` sums over buckets
   with `q > 0` and `nact > 0`, which is all `K`. So **`ΣN/D_N = M/K`**.
3. Hence

   `excess_pc(t) = A_occ(t)/D_G(t) − M(t)/K`,

   where `A_occ` is the observed cumulative activity summed over currently occupied
   buckets. The flux magnitude, the mean-field drift term and `Np` have all dropped out.
   Nothing about any neutral process remains, except the constant `1/K`.
4. **Steady state.** With `D_G` roughly constant and little mass lost to extinct buckets,
   `A_occ ≈ M`, so `d(excess)/dt ≈ T·(1/D_G − 1/K)`. That is extensive in `T` and
   decreasing in occupancy. Guide 07 found it within a factor of about 1.4 of the actual
   dyn slope in 7 configs. At `D_G = 1` it is `T(1 − 1/K)`: the evo_focus eg blow-up
   (≈1.6e4 per step).

### 4.2 Why the frozen null passes by construction

With all mutation rates zero, `g_*_bitflips_last = 0`, so `rate = 0`, `m = 0`, the drift
nudge `m/Np = 0`, and `q_k = o_k/T` exactly. Then:

- `nact_k(t) = Σ_τ o_k(τ) = obs_act_k(t)` for every bucket, up to fractional carry;
- `D_N` is the set of buckets ever observed with `q > 0` *now*, which is `D_G`;
- so `ΣN/D_N = ΣG/D_G`, and the excess is 0;
- `nvar = 0`, so every z is 0 by the guard.

This holds **whatever selection is doing**. A frozen population with standing variation
under strong selection (for example a random-LUT start with `mu_lut = 0`) reads exactly 0.
So the null validates the wiring (flux reaches `m`), not the statistic. The reciprocal-freeze
"validation" of causal-control v2 is this identity (guide 07 §2.1 consequences).

### 4.3 Why the baselines would be wrong even with memory

- **dyn.** Mutation moves mass only between `(i,0)` and `(i,1)`. A neutral baseline must
  keep each input's mass `n_i(t)` fixed and spread only *within* the input. The tracking
  smear moves mass across inputs toward 1/500. Any v-field that uses few inputs (a
  quiescent field, `D_G` small) therefore reads as "concentrated", and so as "excess".
- **eg.** Under the actual operator, each orbit digit is a Markov chain on the hidden pair
  (mask, value). The 12 bits each flip at the same rate, so stationarity is uniform over
  the 4 (mask, value) states: digit 0 (wildcard) has probability ½, digits 1 and 2 have
  ¼ each. The neutral stationary key distribution is the **product measure**
  `π(key) = Π_j (½ if d_j=0 else ¼)`:
  - the all-wildcard key has π = 1/64;
  - each fully specified key has π = 1/4096;
  - the ratio between them is 64.

  A uniform 1/729 baseline is wrong by up to that factor. More importantly, at the
  measured rates the population is nowhere near stationarity within a run. The per-bit
  relaxation rate is `2·mu_egene·(births/Np)` ≈ 7e−5 per step, a time constant of about
  14,000 steps. So the correct neutral reference is a *transient* diffusion out of the
  founder key. Only a stateful shadow can supply that.

---

## 5. Design of the drift shadow

### 5.1 Requirements (from N's brief), and how the design meets each

| Requirement | Design element | § |
|---|---|---|
| Carries state across steps | R replicate populations of individuals, persistent | 5.2 |
| Same fixed bucket spaces | Projection of shadow individuals into eg keys and dyn (input, output) | 5.3, 5.4 |
| Same demography (Np, births, deaths) | Per step: `D_real` uniform victims, then `B_real` uniform parents among survivors | 5.2 |
| Same realized mutation flux | Replay of each real birth's `(nf, na, ne)` onto one shadow birth; positions from the shadow stream | 5.2 |
| No selection | Parent and victim choice independent of genotype | 5.2 |
| Stochastic (finite-population sampling) | Individuals, not mean-field frequencies. Drift and drift-induced concentration emerge, including linkage | 5.2, 5.6 |
| Shadow RNG stream | `shadow-rng` machinery (`xorshift32`, `poisson_sample_from`), with per-replicate sub-streams | 5.7 |
| Honest null tests | Synthetic neutral; frozen-clonal + marker; marker in evolving runs | 5.8, 7 |

### 5.2 The shadow population (shared by both channels)

**Individual.** One full genome: `lut[32]`, `egenes[8]`, `egenes_mask[8]`, `active`. That
is 49 bytes, stored as structure-of-arrays like the real grid. No position, food, `v`, age
or alive flag. One replicate holds up to `cap = N²` individuals and keeps a live count
`n_r`.

**Why individuals and not per-bucket frequencies** **[reasoned]**. Both of the obvious
"closed-form" alternatives lose something the null needs.

- **Independent two-allele processes per dyn input** (guide 07 §5 Option B's sketch).
  The real run is asexual: all 250 loci and the egenome share **one genealogy**. Under
  that genealogy, allele frequencies at different loci co-vary: when a lineage drifts up,
  every locus it carries moves together. Per-locus marginals are right, but the variance
  of any aggregate (a sum over inputs, a z-tail count, an excess) is understated. That is
  exactly the kind of quiet miscalibration that produced 113–173 false z>3 buckets
  before.
- **Independent per-key processes for eg.** The drift unit is the cell, not the egene. A
  cell's (on average 1.1) active keys share one parent, so treating `T_eg` active egenes
  as `T_eg` independent alleles sets the drift size at the wrong value. It also cannot
  represent hidden wildcard values, pseudogene slots or duplication.

Individuals represent all of this by construction, at a cost of 49 bytes each.

**Per-step update**, run after Phase 4 next to `neutral_apply_demography`, for each
replicate `r`:

1. **Deaths.** Remove `D = g_deaths_last` individuals chosen uniformly without replacement
   (swap-remove, O(1) each). Clamp to `n_r`.
2. **Births.** For `b = 1..B` with `B = g_births_last`:
   - choose parent `p` uniformly from the **survivors of step 1** (the pre-birth set,
     indices `< n_r^surv`), so that newborns cannot reproduce in their birth step, which
     matches EvoCA;
   - copy the genome into slot `n_r`;
   - apply mutation replay (below);
   - `n_r++`.
   - **Default: with replacement.** A Wright–Fisher draw, matching Channon. **Option:
     without replacement**, matching EvoCA's at most one child per parent per step.
     Section 7(a) says when the difference matters: the offspring-number variance differs
     by a factor of about (1 − B/n), negligible at 1.2% births per step.
3. **Mutation replay.** During Phase 4 the real run appends `(nf, na, ne)` for each birth
   to a per-step buffer `mut_log[B]` (three small ints per birth; reset at the top of the
   step). Shadow birth `b` takes record `π(b)`, where `π` is a uniformly random
   permutation from the shadow stream, and applies **the same operator code** with these
   counts:
   - `nf` LUT flips at uniform positions in `[0,250)`, or over the real run's current
     `active_bits` when `restricted_mu` is on, exactly as `neut_mutate_inplace`;
   - `na` presence toggles with the Negene→0 rejection, then dup-on-activate at
     `p_dup_egene`, then `ne` value/mask flips over 8 slots × 12 bits.

   **Refactor requirement.** Factor the Phase 4 operator into one function,
   `genome_mutate(genome*, nf, na, ne, next)`, called by both the real run (with counts it
   just drew from `rng_next`) and the shadow (with replayed counts and `shadow_next`). One
   operator, two callers. Then the kernels cannot drift apart. (CLAUDE.md: build the
   argument set once and pass it to both.)
4. **Incremental tallies** (so the projection costs O(events), not O(Np)). On every death
   subtract, and on every birth add, the individual's contribution to:
   - `c_r[i]`, the count of individuals with LUT bit `i` = 1 (250 counters);
   - `e_r[key]`, the count of active egenes with that ternary key (729 counters).

**Why replay realized counts instead of the nominal rates.** It is literally "the same
realized mutation flux", and it costs nothing. It stays correct if mutation rates become
per-cell or heritable (the `mu_lut_gene` scaffold). It also removes one source of
real-vs-shadow variance. What it conditions on: if selection acts on *mutation rate*
(mu-genome evolution), the realized flux is itself a product of selection, and the null
absorbs it. State that as part of the estimand when mu-genome runs use this shadow.

### 5.3 eg projection and kernel

**Projection.** Real side: `fs_eg_obs_pop[key]`, as today, from all active egenes of alive
cells. Shadow side: `e_r[key]`. Both are per-step masses with `T_eg^real` and `T_eg^shadow`,
which are **not** forced equal. The shadow's Negene drifts neutrally (presence toggles,
no egene tax). The real Negene is under selection (`tax_per_egene`).

Excess is computed on **shares** `e_r[key]/T_eg^shadow`, rescaled to the real mass:

`a^S_r,key(t) += T_eg^real(t) · e_r[key](t)/T_eg^shadow(t)`

This puts both sides in the same activity units, and it conditions on the real `T_eg`.
Selection on Negene itself is then **not** in the eg key statistic. It is reported
separately as `negene_excess` (real mean Negene minus the replicate mean) because, unlike
key identity, it is a distinct heritable trait.

**The implied kernel** (documentation and Python reference; the C uses the operator
directly). Per active slot per birth, every one of the 12 (mask, value) bits flips with
probability `ρ = P(odd Poisson(mu_egene) count) = (1 − e^(−2·mu_egene))/2 ≈ mu_egene`,
independently. Per orbit digit, on the hidden 4-state chain:

| from \ to | wild (m=0,v=0) | wild (m=0,v=1) | off (1,0) | on (1,1) |
|---|---|---|---|---|
| wild, v=0 | 1−2ρ+ρ² | ρ(1−ρ) | ρ(1−ρ) | ρ² |
| wild, v=1 | ρ(1−ρ) | 1−2ρ+ρ² | ρ² | ρ(1−ρ) |
| off | ρ(1−ρ) | ρ² | 1−2ρ+ρ² | ρ(1−ρ) |
| on | ρ² | ρ(1−ρ) | ρ(1−ρ) | 1−2ρ+ρ² |

- Projected onto the key, value flips under a wildcard are **silent**.
- `wild → off/on` happens at total rate `ρ`, split by the hidden value. In the
  hidden-state stationary approximation that is ½ each; it is exact in the individual
  model.
- `off ↔ on` and `off/on → wild` each happen at rate `ρ`.
- The key kernel is the 6-fold product.
- The slot-level events sit on top:
  - presence toggles at `mu_egenome` per slot, which add or remove a key from the cell's
    multiset;
  - duplication, which copies an active key into a newly active slot;
  - activation of a pseudogene, which brings in its independently drifted key.

  These are why the individual model is preferred to any closed-form kernel.

**Seeding.** At enable time `t0`, copy every alive cell's full egenome into every
replicate. The `run_sim` founder is monoclonal (value `0b000011`, mask `0x3F`, one active
slot). Over a 15k-step run the shadow is a transient diffusion out of that key (§4.3), and
the stochastic replicates matter most early, while few keys are occupied.

### 5.4 dyn projection: three designs and a recommendation

The question is what, in a population with no positions and no `v`, plays the role of
"cells querying input `i`".

| | Real side | Shadow side | Conditions on | Problem |
|---|---|---|---|---|
| **A.** Per-input independent two-allele WF (guide 07 §5 sketch) | `dyn_pop[2i+o]` (queried transitions) | per input `i`, a scalar `p_i` with WF drift at size `n_i` or `Np`, mutation at the realized flip rate, mass `n_i·p_i` | `n_i(t)` | Loses linkage (§5.2). "Effective number of carriers" is ill-defined: every individual carries an allele at every input, and only expression is limited to queriers. Inherits problem B below |
| **B.** Individual shadow, projected through observed inputs | `dyn_pop[2i+o]` (queried transitions) | `n_i(t)·c_r[i]/n_r` and `n_i(t)·(1 − c_r[i]/n_r)` | `n_i(t)` | **Dynamical self-selection** (below) |
| **C. (recommended)** Individual shadow vs real **population** allele frequencies, both weighted by observed query mass | `n_i(t)·P_i(t)`, `P_i` = fraction of **all alive cells** with LUT bit `i` = 1 | `n_i(t)·c_r[i]/n_r` | `n_i(t)` | Weight–frequency covariance (below); needs real-side allele counts |

**Dynamical self-selection (why B, and A, are confounded)** **[reasoned; measurable as
B − C]**. In the real run, *which genotypes query input `i`* is not independent of their
output at `i`, because the output determines the next state, and so the next input.

- Example: input 0 (v=0, all rings empty). Cells whose output there is 0 stay at input 0
  for as long as their neighborhood stays empty. Cells whose output is 1 leave it at the
  next step.
- So the queriers of input 0 are enriched for output 0 relative to the population
  frequency, purely from the CA's fixed-point structure, with no difference in
  reproduction.
- Generally, outputs that keep a cell at the same input are overrepresented among that
  input's queriers.
- B compares those enriched queriers against a shadow that only knows population
  frequencies. It reads the enrichment as excess, a non-Darwinian false positive whose
  size depends on how quiescent the v-field is. That is the same bias, rewarding simple
  v-fields, that guide 07 §2.1 found in the tracking shadow's uniform baseline, entering
  by a different route.

**C** removes it: both sides use frequencies over the whole population, and only the
*weights* come from the queries.

- C is a classical population-genetic test of whether allele frequencies at the expressed
  loci depart from neutral drift plus mutation, with each locus weighted by how much it
  is expressed.
- Selection that maintains an output at a frequently used input (purifying), or sweeps a
  new one in (directional), moves `P_i` away from the shadow's `p_i`, which relaxes toward
  ½ under mutation. That is the effect C is built to detect.
- Never-queried loci have weight 0. That retires the "silent churn" problem by
  construction: a LUT flip at an unused input has no weight on either side.

**C's residual issue.** `n_i(t)` and `P_i(t)` can covary through the dynamics: an input
can become frequent *because* of the population's output at it. The excess then contains
a term `Cov_t(n_i, P_i)` that the shadow, whose `p_i` is independent of `n_i`, does not
match. Its sign is not fixed. Measure it in validation (c) by comparing C with a variant
that uses lagged weights `n_i(t − W)`, or run-mean weights `n̄_i`. If the difference is
material, use run-mean weights. That makes the dyn statistic a fixed-weight locus test, at
the cost of non-causal (end-of-run) weights.

**Recommendation.** C as the dyn drift-shadow statistic. Also emit B's real side minus C's
real side, `Σ_i n_i·(P̂_i^queriers − P_i)` in activity units, as **`dyn_expression_bias`**,
a descriptor of dynamical self-selection. This is new information that neither shadow
gives today, and it shows how much of the old dyn excess was this effect.

**Real-side cost of C.** It needs `P_i(t)` for the real population, which means 250
counts over alive cells. Recounting every step is about Np·250 bit tests, roughly 2 ms at
N=256, more than the step itself. Maintain them **incrementally** in the real step:
subtract on tax death (Phase 2c) and on eviction, add on birth after mutation. That is
O(250) per event with no RNG, so it is a pure observer. It is still a hook inside the real
step, which is a structural-risk moment (section 10).

### 5.5 Seeding, enabling, windows and cadence

- **Enable** `evoca_fs_drift_enable(R)` copies the alive genomes into all R replicates at
  `t0`, zeroes all drift accumulators, and starts a **separate real-side accumulator**
  from `t0`. Today's `dyn_act` and `fs_eg_obs_act` are cumulative since init; the drift
  statistics must use the same window on both sides. `run_sim` enables it at t=0 when
  `fs_shadow` includes `'drift'`.
- **Re-enable mid-run** resets and reseeds. This is the "re-baseline at the switch" fix
  that memory note `dyn-shadow-frozen-rich-substrate` asks for in #3c-style phase-switch
  designs.
- **Resync period `W`** (optional, default ∞). Every `W` steps, reseed the replicates
  from the real population and close a window. `W = ∞` measures excess accumulated since
  `t0`. A finite `W` measures ongoing selection within windows. The memoryless tracking
  shadow is the mean-field `W = 1` limit, which is the cleanest way to say what it lacks.
  For long runs, report the last window as well as the cumulative value, so that early
  founder-sweep transients do not dominate a ranking.
- **Cadence.** The shadow advances **every step** (demography is per step). Activity
  accumulates per step on both sides (500 + 729 multiply-adds per replicate). Only the
  *statistics* (§6) are computed at `sample_every`, so cadence affects resolution, not
  correctness. The Channon convention (activity counted per sample) is not used here,
  because the fixed-space real sides are already per step.

### 5.6 Replicates, cost and memory

**Why replicates.** A drift-shadow statistic needs a null distribution. Analytic WF
variance does not survive linkage, the transient diffusion from the founder, or the
cumulative-activity functionals. Replicates give the Monte Carlo null for free: under H0
the real run is exchangeable with the R shadows, so the rank of the real statistic among
the R+1 values gives an exact p-value **if the shadow is the right neutral model**.
Section 5.8 covers the "if".

**Cost** (N=256, JOINT base: Np ≈ 15.4k, B ≈ 181, D ≈ 178 per step, 1.54 ms per step):

| Work per replicate per step | Ops | Estimate |
|---|---|---|
| 178 victims, swap-remove, tally subtract (250 LUT bits + ≤8 keys) | ~178 × 260 | |
| 181 births: copy 49 B, replay ≈2.5 + 0.33 flips, tally add | ~181 × 270 | |
| activity accumulation (500 + 729) | ~1.2k | |
| **total** | ≈ 1e5 simple ops | **≈ 0.05–0.1 ms** |
| real-side allele counts (C), once, not per replicate | ≈ 1e5 | ≈ 0.05–0.1 ms |

- So R = 8 adds about 0.5–0.9 ms per step (+30–60%). R = 16 roughly doubles step time.
- A bit-sliced vertical counter (carry-save adds of the 256-bit LUT into 8 bit-planes)
  cuts the 250-op tally to a few dozen word ops if this becomes the bottleneck.
- At N=128, everything scales by about ¼.
- These are estimates from op counts, not measurements. The pre-merge test notebook must
  time them.

**Memory.** 49 B × N² per replicate is 3.2 MB at N=256. R = 8 is about 26 MB; R = 16
about 51 MB. Plus 8 × (250 + 729 + 4,096 optional raw (mask, value) counts) per replicate.
That is fine for pool workers.

**Recommendation.** R = 8 by default. A 1-in-9 Monte Carlo p-value has a floor of 0.11,
so per-run binary tests at α = 0.05 need R ≥ 19. But ranking runs needs a z-type effect
size, which R = 8 estimates adequately, and calibration (section 7) uses many synthetic
runs, not many replicates per run. Offer R = 19 or 39 for confirmatory single-run claims.

### 5.7 RNG

- The brief says the drift shadow uses the `shadow-rng` stream. I recommend the
  **machinery** (`xorshift32`, `poisson_sample_from`, `shadow_seed_from`) with **its own
  state words**, not literally the Channon shadow's `g_shadow_rng`.
- **Reason.** If both shadows draw from one stream, enabling the drift shadow changes the
  Channon shadow's *realization*, so `excess_pc` and `a_new` would stop being
  bit-reproducible across `fs_shadow` settings. That is the same kind of coupling
  `shadow-rng` was written to remove between shadow and dynamics.
- **Proposed derivation.** Replicate `r` uses `s_r = shadow_seed_from(seed ^ TAG_DRIFT ^
  r·0x9E3779B9)`, with a tag distinct from "SHADWRNG". All are reseeded in
  `evoca_set_seed`.
- The marker trait (§5.8) needs a stream too, because its mutations must not touch
  `rng_next`.
- xorshift32's period (2³² − 1) is ample: at about 500 draws per replicate per step, a
  15k-step run uses about 7.5e6 draws.
- **Dependency.** The drift-shadow branch builds on `shadow-rng`, which must merge first.
  `shadow-rng` was cut from `530c334`, which is older than `main`'s Tier-2 merges, so
  expect a rebase.

### 5.8 Null tests, honestly

**What a frozen run now gives.** All mutation off means replayed counts are all zero, and
the shadow is pure neutral sampling from the `t0` composition. The real run is selection
plus the same sampling.

- **Clonal founder** (the `run_sim` default: GoL LUT, single egene): both sides stay
  monoclonal. Every frequency is 0 or 1 on both sides, and the excess is **exactly 0 for
  a real reason**: there is nothing to select. Keep this as a wiring test.
- **Founder with standing variation** (random LUTs, or an evolved population frozen
  mid-run): the excess is **expected to be positive**, since selection sorts standing
  variation. A positive reading is a finding, not a failure. The old "frozen ⇒ 0"
  assertion must be dropped for this case. The memory note's frozen-rich-substrate case
  now reads selection on a frozen evolved substrate. Whether that is "interesting" is a
  scientific question, no longer an artifact.

**The right null checks** (details in §7):

1. **Synthetic neutral, exact.** Run the shadow's own process as "real". This checks the
   implementation and the statistic.
2. **Frozen-clonal + neutral marker in EvoCA.** All genome mutation off and a clonal
   genome, so nothing is selectable. Add a **neutral marker**: a per-cell heritable byte
   pair (12 bits, projected through `egene_ternary_key` so that it lives in the same
   729-key space with the same kernel). It is mutated at a chosen rate from its own
   stream, and read by nothing in the dynamics.
   - The marker experiences the **real spatial demography**: local reproduction,
     food-driven birth sites, evictions. The drift shadow, panmictic by design, then
     reads marker excess.
   - Under H0 this is the empirical null for the misspecification the synthetic test
     cannot see: spatial structure and its effect on effective size.
   - The `mu_lut_gene` scaffold shows the pattern: per-cell, inherited, read nowhere.
     Recommend a dedicated `marker` array rather than overloading it.
3. **The marker inside evolving runs measures the draft floor.** The genome is asexual,
   so a neutral marker hitchhikes on selected lineages. Its excess in an evolving run is
   the *linked-selection baseline*: the excess any neutral trait shows because selection
   elsewhere reshapes the genealogy. This is not a false positive of the test, since
   selection is present. But it is the floor below which a channel's excess cannot be
   *attributed* to selection on that channel. Report `eg` and `dyn` excess together with
   the marker excess from the same run.

### 5.9 What the drift shadow does not see (state these as the estimand)

1. **Selection through input frequencies.** dyn conditions on `n_i(t)`. Rules that make
   useful inputs common are not credited. The test is conservative in this direction.
2. **Selection on population size or on Negene.** Conditioned out of the key statistics.
   Negene is reported separately.
3. **Channel attribution under linkage.** eg excess can be LUT selection hitchhiking, and
   the reverse. Only the marker floor and the cross-channel pattern can separate them.
   Per-bucket z-scores are hitchhiking-contaminated even when the run-level statistic is
   right.
4. **Spatial structure.** The shadow is panmictic. Section 7 measures the consequence.

**Alternative considered: a spatial shadow.** Shadow genomes live on the real lattice, die
where real cells die, and each shadow birth at site `y` copies a uniformly random alive
neighbor of `y`. It reproduces local structure. It also **leaks selection**: at an
expansion front an empty site often has one alive neighbor, the real parent, so the
shadow inherits the real, selection-driven genealogy. That is not neutral. It was
rejected in favor of a panmictic shadow plus the marker calibration.

---

## 6. Excess and z-scores under the drift shadow

### 6.1 Why Ā = ΣA/D is the wrong functional here

Bedau–Packard's mean cumulative activity per present component discriminates in **open**
spaces, where components go extinct and `D` carries the signal.

In these small fixed spaces at EvoCA's mutation supply, every component is present almost
always on both sides:

- stylized dyn-like test: `D_real` = 79.5 vs `D_shadow` = 79.5 of 80 at θ ≈ 2, and 80.0
  vs 80.0 at θ ≈ 10;
- so Ā ≈ M/D on both sides, and the excess is about 0 under selection as well **[run]**.

(The tracking shadow "worked" only because its `D_N` was pinned at `K`.) For eg it
partly survives (`D_G` 1–250 of 729), but it is still dominated by occupancy. **Do not
carry `ΣG/D_G − ΣN/D_N` over to the drift shadow.**

### 6.2 Per-component quantities (both channels)

For each component `k` (eg key; dyn `(i,o)`):

- `A^G_k(t)`, the real cumulative activity from `t0`;
- `A^S_{r,k}(t)`, the activity of replicate `r`, on the same mass scale (§5.3 rescaling
  for eg; for dyn, both sides already have mass `n_i(t)` per input);
- `μ_k = mean_r A^S_{r,k}`, `σ_k = sd_r A^S_{r,k}`;
- **excess** `e_k = A^G_k − μ_k`. Note `Σ_k e_k = 0` identically: both sides carry the
  same total mass, so excess is a *redistribution*. This is why aggregates must be
  distributional, not sums;
- **z** `z_k = e_k / σ_k`, defined when `σ_k > 0`, with a floor `σ_k ≥ max(1, 0.01·μ_k)`
  so that near-deterministic components do not divide by noise. The floor is needed
  because at θ ≫ 1 the replicate spread can be tiny (§6.5).

These replace `get_eg_zscores` / `get_dyn_zscores` under the drift shadow. They are
calibrated in the sense of section 7.

### 6.3 Run-level statistics (recommended set)

The names below assume "drift" (section 9).

| Key | Definition | Reads as |
|---|---|---|
| `dyn_drift_consistency` | `Z_dyn = Σ_i w̄_i·ζ_i / sqrt(Σ_i w̄_i²)`. Here `ζ_i` is the replicate-standardized cumulative per-input homozygosity `F_i = (A_i0² + A_i1²)/(A_i0 + A_i1)²`, and `w̄_i` is input `i`'s share of query mass over the window | Outputs at expressed inputs are held more consistently than drift plus mutation would hold them. Primary dyn statistic |
| `eg_drift_concentration` | `log Simpson(A^G) − mean_r log Simpson(A^S_r)`, with `Simpson(A) = Σ A_k²/(Σ A_k)²` over 729 keys | Real eating-key activity is more concentrated than neutral diffusion from the founder. Primary eg statistic |
| `*_drift_z` | the above standardized by the replicate sd of the same statistic (each replicate scored against the others, leave-one-out) | Effect size in null-sd units |
| `*_drift_p` | `(1 + #{r: stat_r ≥ stat_G})/(R+1)` | Monte Carlo p-value (floor 1/(R+1)) |
| `*_drift_tv` | `½ Σ_k |e_k| / M` (total-variation distance between real and mean-shadow activity shares) | Any departure, two-sided. Diagnostic for misspecification |
| `*_drift_nhi` | `#{k: z_k > 3}` minus its replicate expectation (each replicate scored against the others) | Calibrated per-bucket tail count. Replaces the uncalibrated z>3 count |
| `negene_drift` | real mean Negene minus the replicate mean | Selection on egene number |
| `marker_drift_*` | the eg statistics computed on the marker | Draft floor (§5.8) |

**Intensive by construction (P3).**

- `F_i`, Simpson index and TV distance are scale-free in `M`.
- `Z_dyn` and the z versions are in null-sd units.
- So population and run-length effects enter only through power, not through magnitude.
  That answers guide 07's third requirement.
- Rank configurations on the *intensive* effect, `eg_drift_concentration` and
  `Z_dyn`/`ζ`-based, not on z alone, which grows with power (√M).

### 6.4 Relation to existing readouts

- **`dyn_eff`** (participation ratio of `dyn_act`) is `1/Simpson` of the real *queried*
  dyn activity. Under the drift shadow its neutral reference is the same functional on
  the shadow projection. "Is low `dyn_eff` selection or silence?" (guide 07 §2.5) becomes
  a comparison with a number attached.
- **`A_new` for fixed spaces.** The per-component band is `[q_{1−α}(A^S_{·,k}), ∞)`.
  `A_new,fs(t) = (1/D_G) Σ_{k: A^G_k > band_k} A^G_k`, in Channon's normalization; a
  recent-crossing window is optional. With R = 8 the per-component quantile is crude, so
  use `μ_k + 3σ_k` instead. This is an optional consumer and is not needed for [ACT-D1](../Research/Activity-metrics/ISSUES.md).

### 6.5 What the stylized check says about calibration and power [run]

Model:

- Np = 500 (484 on a 22×22 torus for the spatial case), L = 40 loci, weights `w_i ∝ 0.7^i`
  (skewed accessibility; the top input carries 30%);
- monoclonal founder, 25 deaths and 25 births per step (5%), 1,500 steps (75
  generations);
- per-locus flip probability μ per birth;
- R = 19 panmictic neutral replicates.

"Real" is either the neutral process itself, a spatial voter-type neutral process (birth
copies a random Moore neighbor), or a selected process: parent weight (1 − s)^(number of
allele-1 at the first k loci), purifying.

| Real process | μ | θ = 2Np·μ | generations | Reject at 0.05, Simpson stat | Reject at 0.05, `Z_dyn`-type stat | runs |
|---|---|---|---|---|---|---|
| neutral panmictic | 0.002 | 2 | 75 | 0.065 | 0.055 | 200 |
| neutral panmictic | 0.01 | 10 | 75 | 0.035 | 0.045 | 200 |
| neutral panmictic | 0.01 | 10 | 150 | 0.070 (sd z 1.14) | 0.080 (sd z 1.15) | 100 |
| neutral **spatial** (voter, 22×22) | 0.002 | ≈2 | 75 | 0.030 (sd z 0.85) | — | 100 |
| neutral **spatial** (voter, 22×22) | 0.01 | ≈10 | 75 | 0.020 (sd z 0.88) | 0.040 (sd z 0.96) | 100 |
| selected, s=0.05, k=1 | 0.01 | 10 | 75 | 0.52 (mean z 1.8) | 0.43 | 60 |
| selected, s=0.05, k=1 | 0.01 | 10 | 150 | **0.83** (mean z 3.4) | 0.73 | 40 |
| selected, s=0.02, k=1 | 0.01 | 10 | 150 | 0.15 | 0.20 | 40 |
| selected, s=0.02, k=5 | 0.01 | 10 | 150 | 0.58 | **0.75** (mean z 3.4) | 40 |
| selected, s=0.01, k=1 | 0.01 | 10 | 75 | 0.13 | 0.12 | 60 |
| selected, s=0.005, k=5 | 0.01 | 10 | 75 | 0.12 | 0.13 | 60 |
| selected, s=0.02, k=5 | 0.002 | 2 | 75 | 0.08 | 0.15 | 60 |

Readings:

- **Calibration.** Panmictic neutral rejection rates fall inside the binomial 95% bands
  (200 runs: [0.02, 0.08]; 100 runs: [0.01, 0.09]). The design is calibrated when its
  model is right. The z sd of 1.09–1.15 is slightly above 1. Part of that is expected:
  real-vs-replicate-mean has variance σ²(1 + 1/R), and the replicate sd is a t-type
  estimate. Part may be skew. Use the rank p-value, not a Gaussian reading of z, for any
  binary claim.
- **Spatial misspecification is conservative in this model.** The structured neutral
  population drifts *less* in global frequency than the panmictic shadow, so the null sd
  of z is 0.85–0.88 and the rejection rate is 0.02–0.03. A one-sided selection test stays
  valid, with lost power. EvoCA's spatial process is not a voter model (birth sites are
  food-driven and evictions are directed), so this must be re-measured with the marker
  null (§7c) before it is relied on.
- **Power is modest per run and grows with generations.**
  - One selected input carrying 30% of the query mass, at `Np·s` = 25, is detected in
    about 45–50% of runs after 75 generations, and about 75–85% after 150.
  - Selection spread over 5 inputs at `s` = 0.02 is detected in about 75% of runs after
    150 generations, where the per-locus `Z_dyn` statistic beats the global Simpson.
  - Weak selection (`Np·s` ≤ 10, one locus) is not detectable in a single run of this
    size.
  - EvoCA runs are larger (Np ≈ 15k, about 180 generations in 15k steps), which should
    raise power substantially, but that is an extrapolation.
  - Two-statistic summary: Simpson is better for one strong locus, `Z_dyn` for several
    moderate ones. Report both until EvoCA data say which matters.
- **θ ≫ 1 makes the shadow nearly deterministic** **[reasoned]**. At mutation–drift
  balance the per-locus frequency variance falls roughly as 1/θ, and EvoCA's θ ≈ 300 per
  locus is 30× the stylized θ = 10. I did not measure the replicate spread in EvoCA. Expect
  very large `z` values in EvoCA from *any* systematic departure, including
  misspecification and draft. That is why the effect size against the marker null, not
  replicate z, is the ranking currency.

The stylized model is mine, chosen for speed, not fitted to EvoCA. It is an existence
check of calibration and a rough power scale, not an EvoCA calibration.

---

## 7. Validation plan and acceptance criteria

All tests enter through the public boundary (`run_sim` and the `EvoCA` wrapper). Each
regression test is shown to fail with its fix reverted (CLAUDE.md).

### (a) Synthetic neutral (exact null), Python reference implementation

- **Reference.** A pure-numpy implementation of §5.2–5.4 and §6.3, taking an event stream
  (D, B, per-birth `(nf, na, ne)`, `n_i(t)`, `T_eg(t)`). The C implementation exposes a
  test hook that consumes the same stream. **Parity:** with xorshift32 ported to Python
  and the same seeds, the shadow tallies `c_r[i]` and `e_r[key]` are **bit-identical** at
  every step for 1,000 steps (N=32-scale streams).
- **Calibration runs.** "Real" = one extra replicate of the shadow process, driven by the
  same stream. 400 runs per cell of {θ ∈ {0.5, 2, 10, 50}} × {accessibility uniform,
  geometric 0.7, geometric 0.9} × {eg-like 729-key product kernel, dyn-like 250 loci}.
- **Acceptance criteria:**
  - rejection rate at α = 0.05 in **[0.028, 0.072]** (binomial 95% band for n=400) in
    every cell;
  - KS test of the p-values against uniform on {1/(R+1), …} with p > 0.01;
  - mean `*_drift_z` within ±0.10 (2 SE); sd within [0.90, 1.20] (expected about
    sqrt(1 + 1/R) plus t-tails; the stylized check gave 1.05–1.15);
  - per-bucket `#{|z_k| > 3}` ≤ 1% of present components, **against 23–35% for the
    tracking shadow on the same skewed model** (guide 07 §2.2). This is the test the
    current construction fails.
- **Demographic-variant robustness.** Real = without-replacement parents (EvoCA-like),
  shadow = with replacement. Same criteria.

### (b) Synthetic imposed selection (power)

- Same generator, with purifying or directional selection `s` at k = 1 or 5 of the
  top-weighted inputs (dyn-like), or a fitness advantage for one key family (eg-like).
- **Acceptance criteria**, anchored to the §6.5 stylized numbers:
  - at Np = 500, θ = 10, 150 generations, R = 19, α = 0.05: power ≥ 0.70 for `Np·s` =
    25 on the top input (stylized 0.73–0.83), and ≥ 0.65 for `s` = 0.02 on 5 top inputs
    (stylized 0.58–0.75), each on the better of the two statistics;
  - power monotone in `s` and in generations;
  - a power curve recorded for R = 8 for the record;
  - a C or Python implementation that falls below these on the same generator has a bug
    or a weaker statistic. Investigate before loosening anything.
- **Constructed EvoCA fixture.** Force output `o*` at the most-queried input in every
  birth (a test-only hook) under otherwise neutral marker conditions. `dyn_drift_z` must
  exceed 5 within 2,000 steps at N=64.

### (c) EvoCA runs

| Test | Setup | Accept |
|---|---|---|
| **P4 observer** | Same seed, `fs_shadow` off / tracking / drift / both | Lattice, `alive`, `lut` and every tracking and Channon key **bit-identical** across all four (needs §5.7 sub-streams) |
| **P5 reproducible** | Same seed twice, drift | All drift keys identical |
| **Clonal wiring** | GoL founder, all μ = 0 | All drift statistics exactly 0 |
| **Frozen-clonal marker null** | GoL founder, genome μ = 0, marker μ ∈ {0.001, 0.003, 0.01}; N ∈ {64, 128}; 20 seeds each; 5,000 steps | Marker `Z` p-values: rejection ≤ 0.10 at α = 0.05 (20 seeds; looser band, since the misspecification is real). Record the marker effect-size distribution: its 95th percentile is the **EvoCA null band** |
| **Draft floor** | Causal-control JOINT and LUT_only bases with the marker on, 4 seeds, N=128, 5,000 steps | Report `marker_drift_*` next to `eg_`/`dyn_drift_*`. No pass/fail, but [ACT-D1](../Research/Activity-metrics/ISSUES.md) rankings use channel effect sizes **minus** the marker floor |
| **Standing-variation freeze** | Random-LUT founder, μ = 0, selective base | `dyn_drift_z` > 0 in ≥ 18/20 seeds. This is the test that inverts the old frozen null |
| **Old-key parity** | `main` vs branch, `fs_shadow='tracking'` | Every pre-existing `run_sim` key bit-identical (`tests/test_metrics_regression.py`) |
| **Dynamical self-selection** | GoL-like quiescent vs active configs | Report `dyn_expression_bias`. Expected larger in quiescent fields; a hypothesis, recorded rather than asserted |
| **Timing** | N=256, R ∈ {1, 8, 16} | Record ms per step. Accept ≤ +60% at R = 8, or revisit the bit-sliced tally |

---

## 8. Open design questions for N (each with a recommendation)

| # | Question | Options | Recommendation |
|---|---|---|---|
| Q1 | dyn real side | B (queried transitions, continuity with `dyn_act`) vs C (population allele frequencies, weighted by query mass) | **C**, and emit B − C as `dyn_expression_bias`. B is confounded by dynamical self-selection (§5.4) |
| Q2 | RNG | literally the Channon `g_shadow_rng` vs the same machinery with per-replicate sub-streams | **Sub-streams.** Sharing makes Channon realizations depend on `fs_shadow` (§5.7) |
| Q3 | Mode surface | single `fs_shadow` ∈ {tracking, drift} vs a set {tracking, drift, both} | **{'tracking', 'drift', 'both'}**, with `run_sim` default **'both'**. The drift shadow is a pure observer on its own stream, so reproducibility costs nothing: the tracking keys stay bit-identical. Default the C library to drift **off** (explicit enable), so direct C and notebook users see no change |
| Q4 | Replicates | 1 / 8 / 19+ | **R = 8** default for scans (effect sizes); 19 or 39 for confirmatory single-run claims |
| Q5 | Neutral marker trait | add one (small C plus its own stream) vs rely on synthetic tests only | **Add it.** It is the only null that sees EvoCA's spatial demography, and the only measure of the draft floor |
| Q6 | eg mass conditioning | condition on real `T_eg` (key shares only) vs let shadow Negene drift into the key statistic | **Condition on real `T_eg`.** Report `negene_drift` separately |
| Q7 | Parent sampling | with replacement (WF, Channon) vs without (EvoCA's at most one child per step) | **With replacement** by default, matching Channon. The variant is a validation (a) robustness cell |
| Q8 | Window | cumulative from `t0` vs resync every `W` | Cumulative by default, plus a last-window value with `W = n_steps/2`. Resync is opt-in |
| Q9 | dyn weights | contemporaneous `n_i(t)` vs run-mean `n̄_i` | Contemporaneous, with run-mean as a check (§5.4). Switch if validation (c) shows a material `Cov(n_i, P_i)` term |
| Q10 | Channon fixes in the same branch | full-genome hash, isolated-occurrence screen, banded `A_new` | **No.** Channon stays as is (N's decision). Note that the drift population, carrying full genomes, could later supply a full-genome hash projection if wanted |
| Q11 | `_PARAM_KEYS` | add `fs_shadow` / `fs_reps` there vs keep them as `run_sim` keywords | **Keep them out of `_PARAM_KEYS`.** That tuple is forwarded to `EvoCA.init()` as dynamics metaparameters (`init_kwargs`), and these are observer settings. Add them as `run_sim(..., fs_shadow='both', fs_reps=8, fs_resync=None, marker_mu=0.0)`, echo them in the output dict, add them to `save_scan_config`, and put a `'measurement'` block in the `.evoca` recipe (`export_recipe` version 4) |
| Q12 | Rename tracking keys | keep `eg_excess_pc`, `dyn_excess_pc` vs alias them as `*_concentration` | Keep the old names for continuity. Add aliases `eg_concentration_*` and `dyn_concentration_*` and document them as such (guide 07 Option E) |

---

## 9. Naming

| Candidate | For | Against |
|---|---|---|
| **drift shadow** | Population-genetic "drift" is the standard name for the neutral sampling process. It contrasts cleanly with "tracking" (follows the observed composition) | At EvoCA's θ ≫ 1 the visible content is mutation memory more than sampling drift. A purist could object, though "neutral drift" conventionally includes mutation |
| dynamic shadow | Says "has dynamics" | **Collides with the `dyn` channel:** "the dyn dynamic shadow", and `dyn_dynamic_*` keys. Reject |
| stateful shadow | Exactly the defect fixed | Describes the implementation, not the model. Stays apt only until the next stateful variant |
| population shadow / replica shadow | Emphasizes individuals, or replicates | The Channon shadow is also a population, so "population" does not distinguish. "Replica" names an accessory feature |
| Moran / WF shadow | Precise model name | Commits to a sampling scheme that is an option (Q7) |
| neutral shadow | Literature term | Every shadow here is meant to be neutral; it does not distinguish |

**Recommendation: "drift shadow"**, with the key prefix `<channel>_drift_*` and C prefix
`fs_drift_*`. Rename the old one "tracking shadow" in docs and comments.

---

## 10. Structural risk and implementation notes

These are the moments CLAUDE.md asks to flag:

- shadows become **pluggable** (`fs_shadow` modes);
- a **second call site** of the Phase 4 mutation operator (the shadow replay);
- **observer hooks inside the real step** (the mutation log, real-side allele counts);
- a **second measurement path** for the same quantity.

Suggested branch order:

1. Merge `shadow-rng`.
2. Refactor `genome_mutate`, with a bit-parity test of real trajectories.
3. Add the mutation log and real allele counts (observer only; parity test).
4. Add the drift-shadow core and the Python reference, with parity.
5. Add the `run_sim` surface.
6. Add the marker trait.
7. Pre-merge test notebook (memory: pre-merge protocol).

This is a natural point to offer `/code-review ultra` (user-triggered) before merge.

---

## 11. Glossary

| Term | Meaning here |
|---|---|
| Activity `A_k(t)` | Cumulative per-step mass of component `k` since `t0` (Bedau–Packard) |
| Component | The unit whose activity is counted: genome hash (Channon), ternary egene key (eg), (input, output) transition (dyn) |
| `D_G`, `D_N` | Number of components currently present (real, shadow) |
| Excess | Real minus shadow, for a stated functional of `A` (§6) |
| `T`, `M` | Mass per step; cumulative mass |
| Flux | Realized mutation counts this step (`g_*_bitflips_last`) |
| θ | 2·Np·u, mutation supply per locus per generation; θ ≫ 1 means drift is weak relative to mutation |
| Tracking shadow | Memoryless fixed-space construction (`fs_space_step`): composition reset to observed every step |
| Drift shadow | Stateful replicate populations with mirrored demography and replayed flux, projected into fixed spaces (§5) |
| Channon shadow | Individual LUT population mirrored on births and deaths, hash-keyed (`neutral_*`) |
| Dynamical self-selection | Enrichment of an output among an input's queriers caused by CA dynamics, not reproduction |
| Draft floor | Excess of a neutral marker in an evolving run, caused by linkage to selected sites |
| Marker | A heritable, phenotypically silent per-cell trait used to calibrate the null |
| Monte Carlo p | Rank of the real statistic among `R` shadow replicates |
| `A_new` | Channon 2001: `(1/D)·Σ` activity in a band `[a0, a1]`. In `run_sim`: a count above p90 (different) |

---

## 12. Evidence ledger and references

**Runs this pass** (scratch, outside the repo):

- `rates.py`: births, deaths, flux, Negene and step time at N=128/256, causal-control
  JOINT base (§2.2).
- `drift_mc.py` and `drift_mc2.py`: the stylized calibration and power check (§6.5).

**Code read:** `C/evoca.c` on `main` (sections cited inline); the `shadow-rng` diff
(`git diff 530c334 12f4775 -- C/evoca.c`); `python/evoca_explore.py`,
`python/evoca_py.py`, `python/dyn_analysis.py`; `tests/test_fixed_space_shadows.py`,
`tests/test_zscores.py`.

**Recorded:** `Research/Activity-metrics/plans/2026-10-18_activity-objective-guide.md` (§2.1–2.5, §5);
`Docs/EvoCA_research_directions.md` §2; memory notes `shadow-models-lut-only-mismatch`
and `dyn-shadow-frozen-rich-substrate`.

**References** (only those I am confident of):

- Bedau, M. A., Packard, N. H. (1992). Measurement of evolutionary activity, teleology,
  and life. *Artificial Life II*, 431–461.
- Bedau, M. A., Snyder, E., Packard, N. H. (1998). A classification of long-term
  evolutionary dynamics. *Artificial Life VI*, 228–237.
- Channon, A. D. (2001). Passing the ALife test: activity statistics classify evolution
  in Geb as unbounded. *ECAL 2001*, LNCS 2159. (`Docs/references/channon_ad_ecal2001.pdf`: shadow
  recipe, `Anew` band definition, isolated-occurrence screening.)
- Channon, A. D. (2006). Unbounded evolutionary dynamics in a system of agents that
  actively process and transform their environment. *Genetic Programming and Evolvable
  Machines* 7:253–281.
- Background for §5–6 (standard population genetics, not specific papers): Wright–Fisher
  and Moran neutral models; mutation–drift balance (θ); genetic draft and hitchhiking
  under linkage.
