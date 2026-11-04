# A full-genome neutral reference for evolutionary activity

> **2026-11-03, later (read first).** This plan was written without reading `Docs/shadows.md` or the `drift-shadow` branch. The full-genome individual shadow with the shared mutation operator (§3 here) and the fixed-space projections (§13) already exist there as the drift shadow, validated synthetically, and its panmictic null fails EvoCA's neutral-marker test (marker $z$ from $-3.9$ to $-18$; `plans/2026-10-18_drift-shadow-findings.md`). A panmictic full-genome Channon shadow would inherit that failure. Before any build, re-scope against that branch and the 2026-10-18 D1 (pedigree-conditioned shadow). See `Docs/shadows.md` §14.5. The Channon-specific parts (open-space full-genome keys, the seeding split, a shadow flux history) are the only pieces not covered there.

_2026-11-03. Plan, not yet built ([ACT-N9](../ISSUES.md)). N asked for both this and the
egene-cost stratification of Pareto pair B (done: notebook §4a′); this is the larger piece._

> **Headline.** EvoCA's individual-level neutral shadow (Channon style) is keyed on the LUT
> only, while the real activity table is keyed on the full genome (LUT plus egenes). The plan
> replaces it with a shadow that carries and mutates the **full genome**, under the **same
> mutation operator** as the real run, and records the population history the activity flux
> needs. Because the shadow copies the real run's births and deaths one for one, every
> statistic computed on both sides compares the same demography, so the per-step vs.
> per-generation ambiguity of the N7 comparisons cancels. If it holds up, it becomes the
> standard reference for EvoCA's activity statistics.

## Outline

| § | Section | Purpose |
|---|---|---|
| 1 | What exists | The three references in the code today and what each can and cannot compare |
| 2 | Defects in the current shadow | Key mismatch, seeding inconsistency, no flux history |
| 3 | Design | What a shadow member carries, how it mutates, how it is keyed |
| 4 | Flux on the shadow | History, band choice, the excess statistic |
| 5 | Alternatives considered | Freezing egene mutation; fixed-space shadows only; patching keys |
| 6 | Tests | What must hold bit-exact, the controls that make a test discriminate |
| 7 | Validation runs | Torque runs, with temperature logging |
| 8 | The old-vs-new notebook | What it compares, on which runs |
| 9 | Making it the standard | Outputs, naming, the transition for recorded numbers |
| 10 | Cost | Memory, time, effort |
| 11 | Critique of the first draft | What an opposing expert would object to, and the revisions |
| 12 | Decisions for N | |

## 1. What exists

| reference | level | key | mutation | gives | limits |
|---|---|---|---|---|---|
| Channon neutral shadow (`neut_*`, `nact_*` in `C/evoca.c`) | individuals, 1:1 with the alive count | **LUT hash** | LUT bit flips (restricted to the real run's active bits when `restricted_mu`) | `N_total_activity`, excess, `A_new` threshold (shadow p90) | LUT only; §2 |
| `dyn` drift shadow (`fs_shadow_*`) | closed form over 500 LUT (input, output) buckets | fixed space | realized LUT flip count | `dyn_excess_pc` | no individuals, so no clones, no flux; frozen-substrate artifact ([ACT-Q10](../ISSUES.md)) |
| `eg` drift shadow (`fs_shadow_*`) | closed form over 729 ternary egene keys | fixed space | realized egene flip count | `eg_excess_pc` | as above; void at `D_G = 1` ([ACT-Q3](../ISSUES.md)) |

The real side: `act_*`, keyed by `lut_hash_cache[i] = cell_genome_hash(i)` = FNV-1a over the LUT
bytes followed by the sorted active (value, mask) egene pairs (`cell_rehash`, about line 1120).

Only the Channon shadow is individual-level, so only it can say whether a **genome's clone**
rises faster than drift. That is what the flux needs ([CHEM-N7](../../Chem/ISSUES.md) rerun: at `mu_lut` 0.001 up to
70% of the flux is carried by clones), and what `A_new` was meant to measure.

## 2. Defects in the current shadow

1. **Key mismatch (known, memory `shadow-models-lut-only-mismatch`).** Real keys include
   egenes; shadow keys do not. In any run where egenes mutate, real diversity includes egene
   variation the shadow cannot have. Excess and `A_new` are then undefined; the causal-control
   EGENE_only arm shows it (excess about -26 with the LUT frozen). The N7 and flux runs have
   `mu_egene = mu_egenome = 0.003`, so the clonal-flux comparison would inherit it.
2. **Seeding inconsistency (new, 2026-11-03).** `evoca_neutral_enable` seeds each shadow member
   with `lut_hash_cache[i]`, the **full-genome** hash, but children are keyed with
   `lut_hash_fn(child_lut)`, the **LUT-only** hash. So an unmutated child of a seeded member
   lands in a different bucket from its parent. Every seeded lineage is split in two at its
   first birth, even with egenes fixed. The effect is a transient (the seeded buckets die
   out), but it inflates shadow diversity early and distorts `A_new` thresholds in short runs.
   Confirmed 2026-11-03 (N = 64, all mutation off, random LUTs and egenomes, 30 steps): all
   33 shadow births produced keys the real run never had (the real run had no new keys at
   all). The size of the effect on recorded statistics is to be measured in the notebook (§8).
3. **No population history.** The real table keeps `act_pop_hist` (64 ticks per bucket) for the
   flux; the shadow table does not, so the flux cannot be computed on it.

## 3. Design

**A shadow member** carries the full heritable genome and nothing else:

| field | bytes | as in the real cell |
|---|---|---|
| LUT | 32 | `lut + i*LUT_BYTES` |
| egene values | 8 | `egenes[i*NEGENOME_MAX ...]` |
| egene masks | 8 | `egenes_mask[...]` |
| active-slot byte | 1 | `active[i]` |
| cached hash | 4 | `lut_hash_cache[i]` |

No position, food, age, or tax: the shadow is neutral by construction (uniform parents,
uniform kills, counts mirrored from the real step, as now).

**Mutation.** The real child-mutation code (LUT pass, active-bit pass with dup-on-activate,
egene bit pass; `C/evoca.c` about lines 2240-2320) is factored into one function

```c
static int genome_mutate(uint8_t *lut, uint8_t *eg, uint8_t *egm, uint8_t *active,
                         uint32_t (*rng)(void), int *nf, int *na, int *ne);
```

called by the real step with `rng_next` and by the shadow with `shadow_rng_next`. One operator,
stated once: the shadow cannot drift from the real mutation rules. The refactor must not change
the real run's RNG call order (golden-hash test, §6).

**Key.** `genome_hash_static(lut, egene_buf, nb)` with the egene buffer built by the same helper
as `cell_egene_buf` (factored to take arrays, not a cell index). Seeding uses the same function,
which fixes defect 2.

**Restricted mutation.** As now: the shadow uses the real run's current active-bit set
(`active_bits`, `n_active`). This is a deliberate coupling (the shadow mutates where the real
population's rules are exercised) and is kept.

## 4. Flux on the shadow

- **History.** `nact_pop_hist[nact_cap * FLUX_W_MAX]`, written in `evoca_n_activity_update`
  (one slot per tick), carried through `nact_resize` and `nact_compact` exactly as `act_*` does.
- **One walk.** `act_crossings_walk` (added 2026-11-03 with the cap fix) is generalized to take
  the table (keys, values, history, cap, tick counter), so real and shadow flux run the same code.
- **The band.** Options:
  - (a) each side in its own p20-p30 band. Simple, but the two bands differ, so a difference in
    flux mixes band position with clone dynamics.
  - (b) **both sides in the shadow's band** (the neutral distribution defines "low activity"),
    as `A_new` uses the shadow's p90. Recommended as primary.
  - (c) the real band for both. Reported for symmetry.
- **The statistic.** `clonal_flux_excess = real clonal flux - shadow clonal flux`, in the shared
  band, per tick. Since births and deaths are identical on both sides, per tick and per birth
  give the same sign and ratio; there is no turnover normalization to choose. Also reported:
  total flux on both sides, and the excess of `A_new` and excess activity on the new key.

## 5. Alternatives considered

| alternative | verdict |
|---|---|
| Freeze egene mutation in experiments that use the shadow | Rejected as the standard (egene evolution is half the genome and pair B's second axis). **Kept as a control:** with `mu_egene = mu_egenome = 0` and uniform egenes, old (after the seeding fix) and new shadows must agree exactly in distribution; this is a validation test. |
| Rely on the fixed-space drift shadows | They have no individuals, so no clones and no flux. They stay, as the per-component nulls. |
| Patch only the keys (hash egenes from the real parent into the shadow child) | Wrong: shadow lineages are not real lineages, so there is no real parent's egene to borrow. |
| Lineage-based reference (real parent hashes, `lin_parent_hash`) | A different question (net return per lineage, Lineage-co-evolution); not a neutral null. |

## 6. Tests

All through the public wrapper (`EvoCA`, `run_sim`), each shown to fail without its fix:

1. **Real trajectory unchanged by the refactor**: golden state hashes of a mutating run
   (LUT and egene mutation on, restricted on and off), recorded before the refactor.
2. **Shadow is still a pure observer**: shadow-on and shadow-off real trajectories identical
   (`tests/test_shadow_rng.py` extended to egene-mutating runs).
3. **Seeding consistency**: with all mutation off, an enabled shadow keeps exactly the real
   run's bucket set and counts after many births (fails today: defect 2).
4. **Egene variation reaches the shadow**: `mu_lut = 0`, egene mutation on: shadow diversity
   grows above 1 (today it stays at the seeded set).
5. **Operator parity**: the shadow's realized flip counts per birth (LUT, active, egene) match
   the real run's distribution (`g_*_bitflips_last`) within sampling error at matched rates.
6. **Shadow flux**: shadow crossings and flux uncapped and equal to the walk's sum (as
   `tests/test_activity_flux.py` does for the real side); history survives resize and compaction.
7. **Egene-frozen control** (§5): old and new shadows agree in distribution.
8. **Reachability**: `run_sim(..., shadow=True)` emits the new columns and they change when the
   shadow key changes.

## 7. Validation runs (torque)

Conditions from the flux rerun: W1 and C0 × `mu_lut` {0.03, 0.003, 0.001} × 4 seeds, N = 512,
20000 steps, both shadows on, flux in the three bands, plus the egene-frozen control at
`mu_lut` 0.001. About 30 runs, an hour or so on 24 performance cores.

Temperature: `macfanpro log --rate 1 --duration <run length> --no-expire` alongside; report the
peaks of TCDX and of the hottest `Tp*` key, and fan rpm. No fan commands.

## 8. The old-vs-new notebook

A dedicated notebook, `Research/Activity-metrics/neutral-reference.ipynb` (decision 2), so the
comparison reads start to finish without the long Activity-metrics history. Sections:

1. **Why** (the three defects), each shown on a small run: the seeding split (bucket counts
   before and after the first births), the key mismatch (shadow diversity in an EGENE_only run,
   old vs. new), the flux cap (old reading 8192, new).
2. **Old vs. new on the same runs**: excess activity, `excess_pc`, `A_new`, and the flux
   (total and clonal) with old and new shadows, side by side, on
   - the causal-control arms (the harness doubles as a bug detector: EGENE_only should become
     interpretable),
   - the pure-evo joint winner and its frozen control,
   - the N7 / flux conditions (W1, C0, the `mu_lut` sweep, from §7's torque runs).
3. **What changes for recorded results**: which headline numbers in Activity-metrics and Chem
   move, and by how much; which conclusions survive.
4. **The new statistic**: clonal flux excess against the shadow, by condition and `mu_lut`;
   whether it separates W1 from C0 and selection from drift.
5. **Verdict**: recommend (or not) making it the standard; consequences for [ACT-D1](../ISSUES.md) and the
   Pareto search's pair B.

## 9. Making it the standard

- The new shadow becomes the default; the LUT-only shadow stays selectable
  (`neutral_enable(key='lut')`) for reproducing recorded results, with the seeding fix applied
  to both.
- `run_sim` output keeps the old column names for the old key during a transition, and adds
  the new ones with a prefix (`fg_` for full genome) until the notebook's verdict; then the
  defaults switch and `scan_config.json` records `shadow_key`.
- Recorded scan results are not rewritten; the notebook (§8.3) states what moved.

## 10. Cost

- Memory at N = 512: shadow members 100k × 53 B = 5 MB; shadow flux history
  `nact_cap × 64 × 4 B`, about 64 MB at 250k slots (the real table already pays the same).
  If that matters, `FLUX_W_MAX` can drop to 16 (the flux uses `window = 8`).
- Time: the shadow's per-birth cost rises from one LUT mutation to the full operator, and the
  history write adds one pass per tick over the shadow table; expected below 10% of a step.
- Effort: C refactor and shadow (about a day with tests), wrapper and `run_sim` columns (half a
  day), torque runs (an hour), notebook (a day).

## 11. Critique of the first draft (opposing expert)

1. *"A neutral shadow of the full genome is not neutral for egenes: the real egene count is held
   down by `tax_per_egene`, the shadow's drifts freely, so the shadow will always look more
   diverse."* That is the intended contrast (the tax is selection), but it means the full-genome
   excess mixes "egene count is constrained" with "particular genomes rise". **Revision:** the
   notebook reports the excess within egene-count classes as well as overall, and the Pareto
   use keeps egene-cost strata (§4a′).
2. *"Comparing each side in its own band confounds band position with dynamics."* The first
   draft used own bands. **Revision:** shared shadow band as primary (§4b), own bands reported.
3. *"Per-tick vs. per-generation normalization was the N7 problem; why is it gone?"* Because
   the shadow copies births and deaths one for one; only within one run, not across conditions.
   **Revision:** cross-condition comparisons use the excess (real minus its own shadow), never
   raw flux; the notebook says so.
4. *"Restricted mutation couples the shadow to the real run's exercised bits, so the null is
   not fully independent."* True and already the case; the alternative (unrestricted shadow
   mutation) mutates silent bits the real run never sees and inflates shadow diversity.
   **Kept**, and stated in the doc and notebook.
5. *"Defect 2 may have moved recorded numbers."* It affects early shadow statistics only.
   **Revision:** the notebook measures it on the recorded scan settings before claiming that.

## 12. Decisions for N

1. **Go ahead with the C refactor** (one mutation operator shared by real and shadow, §3)? It
   touches the real step; the golden-hash test guards it. Recommended.
2. **Notebook**: a new `neutral-reference.ipynb` (recommended) or a section in
   `Activity-metrics.ipynb`? *Settled 2026-11-03: N said "put it in old-vs-new"; the notebook
   exists with sections 1-6, and section 7 waits for this build.*
3. **Primary band**: the shadow's (recommended) or each side's own?

## 13. Addendum (2026-11-03): the fixed-space shadows need this too (superseded: see the note at the top; this re-derives `Docs/shadows.md` §4–§5)

Computing N's paired readout (total variation, $\tfrac12\sum_k |G_k - N_k| / \sum G$; C
`evoca_eg_tv`/`evoca_dyn_tv`, commit a5ac683) showed that the fixed-space drift shadows are
memoryless. Each tick they rebuild their composition from the real run's
(`q = (1-m) observed + m uniform`), so both the old excess and TV are real-run statistics
([ACT-Q11](../ISSUES.md); `neutral-reference.ipynb` section 3). The repair belongs in this build:

- **Egene keys from shadow individuals.** The full-genome shadow carries each member's
  slots, so its 729-key histogram is a neutral egene composition with memory and drift.
- **LUT buckets from shadow LUTs, real inputs.** Shadow members do not run the lattice. For
  each LUT input $i$, take the real run's input use $w_i$, and take the fraction $p_i$ of
  shadow members whose LUT outputs 1 at $i$. The shadow bucket mass is
  $w_i (1-p_i), w_i p_i$, so the input marginal is shared by construction. That is the right
  null: mutation changes outputs, not which neighborhoods occur.
- **Statistic.** The paired TV between real and projected shadow, against the TV between two
  independent shadows (replicate RNG streams) as the null spread.
- **Tests.** Two that the memoryless shadow fails: (a) the identity
  `TV = m x TV(real, uniform)` must *not* hold; (b) a neutral egene (eating independent of the
  egene) must read TV(real, shadow) within the shadow-shadow spread. (a) is in
  `tests/test_fixed_space_tv.py` as a strict xfail.
- **Cost.** About half a day on top of §10. The projections are histograms of arrays the
  shadow already holds.
