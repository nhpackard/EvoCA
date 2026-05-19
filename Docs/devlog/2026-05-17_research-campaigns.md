# 2026-05-17 — research campaigns (pure-evolution program + others)

Compute: torque.protolife.com (32 cores). Local: causal-control v2.
Conventions: this file = findings/why; `commit-log.md` = ledger;
`research_board.md` = status/decisions.

---

## Campaign #1 — pure-evo, LUT + egene  (`Scans/2026-05-17_pure_evo_joint`)

**Setup.** 300 random configs, productive box with `m_scale` and
`mu_lut` brackets pushed upward, N=256, 8000 ticks, neutral shadows
on. Ranked post-hoc by the *corrected* turnover-invariant composite:
`excess_pc_slope` (component-normalised whole-LUT-genome, §2),
`eg_excess_pc_slope` (fixed-space egene shadow, S2d),
`dyn_excess_pc_slope` (fixed-space LUT-entry shadow, S2d),
`n_distinct_genomes_temporal_std`, `unique_top_genomes`. 284/300
candidates (16 extinct, 0 saturated).

**Headline finding (significant).** Across **every** top config:
- `excess_pc_slope` ≈ **0.000–0.001** — the component-normalised
  *whole-genome-hash* excess shows essentially **no signal** anywhere.
- `dyn_excess_pc_slope` strongly **positive** (~23–60) and
  `eg_excess_pc_slope` strongly **positive** (~40–135).

Interpretation: at these parameters selection is real and measurable
in the **effective dynamics** (which LUT-entry transitions the alive
population exercises) and in **egene cognition**, but is *not*
detectable as whole-genome-hash excess. This directly vindicates the
§2 thesis: whole-genome FNV-hash churn ≠ effective selection;
`dyn_activity`/`eg_activity` are the sharper probes. Corollary: the
2026-04-27 scans, which optimised raw ΣG−ΣN whole-genome excess, were
chasing a quantity that under correct normalisation is ~flat — their
"`mu_lut` wants 0.06" optimum was an artdefact of the broken metric.

**New optimum (contradicts the old scans).** Under corrected metrics
the top configs favour **low `mu_lut` (0.01)** + **high `m_scale`
(2.5)** — the opposite of the old high-`mu_lut` regime. Low mutation
preserves coherent rule structure while strong eating drives the
effective-dynamics/cognition selection that the corrected metrics
actually see.

**Open design question (→ board R1).** `m_scale = 2.5` is again at
the **grid ceiling** — still unbracketed after pushing it from 1.5 to
2.5. A follow-up needs `m_scale ∈ {2.5, 3.5, 5.0}` to bracket it.

**Status.** Done (491 s on torque). `results.csv` versioned at
`e6525b5`. Feeds #2/#3 (LUT-only → egene-on-winners) and #4
(does spatial RD collapse under this objective?).

---

## Campaign #2 — pure-evo, LUT-only (egene frozen)  (`Scans/2026-05-17_pure_evo_lut_only`)

**Setup.** 180 configs, same grid minus egene axes; `mu_egene =
mu_egenome = tax_per_egene = 0`. 258 s on torque.

**Findings.**
- **0/180 extinct** — with the fully-specified egene frozen, the
  population is extremely robust across the whole box.
- `dyn_excess_pc_slope` is **positive but markedly smaller** than #1
  joint (top configs ~3–20 vs #1's ~23–60). Freezing egene
  *reduces* effective-rule-transition excess. Cross-campaign causal
  inference (the causal-control logic at scan scale): **egene
  co-evolution amplifies LUT effective-dynamics selection** — the two
  subsystems are not independent; cognition evolving alongside the
  rule makes the rule's effective selection stronger.
- `excess_pc_slope` again ≈ 0 (occasionally tiny +ve) — whole-genome
  excess remains a non-signal, consistent with #1.
- **Optimum (contrasts #1):** LUT-only favours **high `mu_lut`
  (0.15, ceiling)** + **high `m_scale` (2.0–2.5, ceiling)** + **low
  `food_inc` (0.010, low edge)**. The mu_lut optimum is therefore
  *conditional on egene*: with egene frozen, high LUT mutation is
  favoured (raw rule exploration must substitute for the adaptive
  cognition that's switched off); with egene co-evolving (#1), low
  mu_lut wins and cognition does the adaptive work. This
  mu_lut↔egene interaction is exactly the kind of metaparameter
  relationship the program is after.
- Spatial: top LUT-only configs keep `correlation_length` ~55–75 —
  RD-scale structure persists under LUT-only pure-evo (relevant to
  #4).
- Multiple axes unbracketed at edges (`m_scale`, `mu_lut`,
  `food_inc`-low) — feeds R1 and suggests a low-`food_inc` extension
  too.

**Status.** Done. `results.csv` versioned. Winners feed #3.

---

## Campaign #3 — egene-on-#2-winners, paired 3a/3b  (`Scans/2026-05-17_pure_evo_egene_on_winners`)

288 cfg (12 bases × 12 egene × 2 modes), 805 s on torque.

| mode | viable/total | eg_excess_pc (mean/med) | dyn_excess_pc | extinct |
|------|---:|---:|---:|---:|
| 3a seeded-joint (LUT evolves) | 127/144 | 20.3 / 16.7 | 14.0 | 17 |
| 3b frozen-LUT (mu_lut=0)      | 78/144  |  5.6 / 2.0  | 0.00 | **66 (46%)** |

3a > 3b eg_excess on **9/12 bases** (+13…+29).

**Finding.** Egene evolution is far weaker (eg_excess 20→6) and
**46% lethal** when the LUT cannot co-adapt. Strong co-evolution /
co-dependence signal: cognition needs the rule to co-evolve for both
efficacy and viability.

**Design limitation (flagged, not buried).** `run_sim` always seeds
`lut='gol'`; with `mu_lut=0`, 3b freezes the LUT at **GoL**, *not* at
an evolved good substrate. So 3b ≈ causal-control EGENE_only relocated
to the winning param regimes — it does **not** isolate "egene on a
*good evolved* substrate", so #3 **strengthens but does not cleanly
resolve** the frozen-substrate confound. The orchestrator over-stated
what 3b would test.

**Decisive resolution (supersedes the #3-as-disambiguator claim):**
1. **#8** (direct lineage joint-retention vs neutral) needs no
   freezing — promoted from corroborator to *the* decisive test.
2. **#3c (new, queued):** evolve LUT-only → `extract_patch` the
   *evolved* LUT → freeze that → evolve egene. The clean
   "good evolved substrate" test 3b was meant to be (uses S2c).

`dyn_excess_pc`=0 in 3b reconfirms the fixed-space null holds
(mu_lut=0 ⇒ no rule evolution ⇒ exactly zero).

---

## Campaign #4 — RD winners under pure-evo: does spatial structure vanish?  (`Scans/2026-05-18_rd_under_pure_evo`)

**Setup.** The 30 `Runs/2026-04-27_top_*` evo+spatial balanced-composite
winners, deduped by metaparam set → **10 distinct recipes** × 3 seeds =
30 runs, N=256, 12 000 ticks, `init='halfplane'`, shadows on. The
direct test of NP's central intuition: robust evolution needs
life-death turnover, so pure-evo should *not* re-select the static
spatial structure the spatial-composite objective rewarded; if it
collapses to static, suspect the (spatial) metric. 89 s on torque,
**0/30 extinct**.

**Headline finding (significant).** RD-scale spatial structure is
**not robust under pure-evo**. `correlation_length` final/mean across
all 30 runs: median **0.23**, mean 0.35, max 0.89 — *every* run's
spatial correlation decays to a fraction of its run-mean. The RD
structure the 2026-04-27 winners exhibit is an artifact of the
**spatial-composite objective they were selected under**, not an
attractor of pure evolutionary turnover. (This refines, not
contradicts, #1/#2's "correlation_length persists": that was a *box
winner / mean* statement; #4 shows the *2026-04-27 spatial winners'*
structure does not survive as an *endpoint* once the spatial objective
is removed — the mean was inflated by a structured early transient.)

**Two regimes (the interesting part).**
- **Dynamically-alive, spatially-washed-out (8/10).** Spatial
  correlation collapses (cL f/m 0.06–0.64) but the population keeps
  strong effective-dynamics turnover (`dyn_excess_pc_slope` ≈ 12–27).
  These are *not* collapsing to static structure — they stay
  dynamically selected, they simply don't preserve the particular
  correlation length the spatial objective rewarded.
- **Static colony-filling blob (2/10):** `cfg49` (top_1) & `cfg30`
  (top_2), both **m_scale 0.6 / mu_lut 0.01 / food_inc 0.025**.
  `largest_patch` final/mean ≈ **1.00–1.01** with mean patch ≈ 45 600
  cells (~70 % of the 65 536-cell grid) — a single frozen blob —
  *and* the **lowest `dyn_excess_pc_slope` in the set (~5.8)**. These
  are precisely the degenerate "suspect the metric" cases NP
  predicted: static structure ⟺ weakest effective-dynamics
  selection. Cf. the local artifact
  `Runs/2026-04-27_R-D_to_static_colony_filling_world.evoca`.

**Interpretation.** NP's intuition is corroborated with a sharpened
form: pure-evo does **not** preserve the spatial-objective's RD
structure; the healthy pure-evo regime keeps life-death turnover
(high `dyn_excess`) while letting spatial correlation wash out, and
the *only* runs that lock into static structure are exactly those
with the **least** dynamic selection signal. "RD spatial structure" is
therefore an objective-induced phenotype, not an intrinsic property of
the evolving dynamics; `dyn_excess`/`eg_excess` (not
`correlation_length`) remain the load-bearing probes — consistent with
the §2 / #1 / causal-control-v2 through-line.

**Status.** Done (recovered: the launching SSH session dropped on a
local-machine DNS failure/reboot *after* the scan completed
autonomously on torque; `run.log` = "Done in 89s", 30 rows, 0 errors;
results pulled to local 2026-05-18). Unblocks **R1** (m_scale bracket
{2.5,3.5,5.0} + low-food_inc, pre-approved to run after #4).

---

## Campaign #8 (QUEUED, fires after #3) — direct lineage co-evolution metric

Goal: replace the *indirect* excess-comparison evidence for LUT↔egene
co-evolution with a *direct* measurement on the S2a lineage field.

Design (ready to execute):
- Run lineage ON at a #3a winner regime (LUT-pure-optimal substrate
  with egene co-evolving) + the JOINT productive corner, multi-seed,
  post-transient window.
- Per parent→child edge, decompose change into ΔLUT (LUT-only hash
  changed?) and Δegene (active egene key changed?) — snapshot per-cell
  (birth_id, parent_id, lut-hash, egene-key) periodically; rebuild
  chains via parent_id.
- Statistic: correlation / excess joint-retention of (ΔLUT, Δegene)
  across *surviving* edges **vs the neutral product of the two
  marginal mutation frequencies**. Mutations are independent by
  construction at the mutation step, so any correlation in surviving
  edges beyond the mutation-rate product is *selection jointly
  retaining co-adapted LUT+egene* — the clean, direct co-evolution
  signal.
- Decisive cross-check with #3: if 3a (good substrate) shows positive
  joint-retention while 3b (frozen LUT) cannot, synergy is real
  co-evolution, not the GoL-frozen-substrate handicap confound.

**Enabler built (2026-05-18).** The "snapshot per-cell lut-hash"
step above needs a per-cell *LUT-only* hash; the activity tracker /
`lineage_parent_hash` only expose the *combined* LUT‖egene genome
hash (`lut_hash_cache`; the [[shadow-models-lut-only-mismatch]]
issue). Added `evoca_get_lut_hashes(uint32_t *out)` — reuses the
existing pure-LUT `lut_hash_fn` over each cell's 32 LUT bytes; purely
additive, no dynamics change (the 40 pre-existing tests pass
byte-for-byte; +2 new in `test_lut_hashes.py` pin egene-invariance —
the property the combined hash fails — and correctness vs an
independent Python FNV mirror; 42/42, clean `-Wall`). On branch
`lut-hashes-accessor` (`09ef368`), local, awaiting the M1 human merge
gate. Once merged, the #8 campaign harness can be built directly:
periodically read `get_lut_hashes()` + active-egene-key + lineage
(birth_id/parent_id) and apply the joint-retention statistic above.

Status: enabler in-review (M1); campaign build unblocked on merge.
(Original "auto-launches once #3 results are in" superseded — #3
done; #8 promoted decisive, gated only on M1.)

---

## Campaign #3c — egene on a frozen *evolved* substrate  (`Scans/2026-05-19_3c_coevolution_substrate`)

**Setup.** Spec `Docs/plan-3c-coevolution-substrate.md`. Single
continuous sim, in-situ phase switch (no patch/S2c). 80 rows
(10 seeds × 3 phase-1 `mu_lut` × {A,C} + 10 B + 10 B′), 0 errors,
650 s torque, launched detached (rode through the usual torque
flakiness). **`determinism_ok=True` on all 30/30 A/C pairs** — the
S1 shared-fork (D-4), the entire basis of the clean contrast, held
exactly. Premise confirmed via the merged `get_lut_hashes()`:
A/C `sub_frac_evolved=1.000`, ~38 000 distinct LUTs at the switch
(phase-1=3000 is amply evolved — D-2's "1000 risks null" worry was
right to be conservative; 3000 is solidly past it). B/B′ frac 0.0.

**Result is mixed and partly methodological — reported honestly.**

1. **The #3 3a≫3b gap was mostly GoL-substrate *viability
   collapse*, not clean co-evolution synergy.** Plain frozen-GoL
   (B), after a 3000-tick settle, is **10/10 extinct** (final_pop 0,
   eg_slope 0 trivially) — exactly NP's "GoL dies to a dull fixed
   point fast" intuition. So a large part of 3b's poor eg_excess was
   *the population dying*, not *egene failing to evolve on a frozen
   substrate*. Use **B′ (viable short-settled GoL) as the real
   3b-analog anchor**: 0/10 extinct, final_pop ~55 600, but
   `eg_excess_pc_slope ≈ 0` (−0.00006) — **frozen GoL ⇒ no egene
   selection even when alive.** B itself mainly proves the regime is
   lethal to settled GoL.

2. **Residual genuine co-evolution effect, but modest and weak in
   absolute terms.** Substrate held byte-identical, the only paired
   contrast that isolates "freezing per se": **C−A `eg_excess_pc`
   slope > 0 in 22/30 (whole phase-2) and 17/30 (late window), mean
   +0.6 / +1.0.** Consistent direction (co-evolving rule helps egene
   selection beyond a frozen evolved rule) — but the *magnitudes are
   tiny* (per-arm |slopes| ≲ 1 vs #1/#3/causal-control's 20–135), and
   it strengthens with substrate quality: by phase-1 `mu_lut`,
   C−A = +0.03 (0.01) / +1.02 (0.03) / +0.80 (0.06). Net mapping to
   the spec inference table: **C ≳ A ≫ B(dead)/B′(flat-zero)** — i.e.
   *substrate-death dominates the 3b deficit, with a small real
   ongoing-co-evolution component on top, on a weak-signal regime.*

3. **Long-timescale (NP point 2): egene selection is NOT robustly
   proceeding at 8000 ticks here.** Early-vs-late phase-2 windows:
   eg_slope early strongly −ve (A −5.7, C −8.1 — switch transient as
   egene mutation comes on), late ≈ 0 (A +0.03, C +1.06). It
   **plateaus**, not grows, by ticks 7–8k in every arm. A clean
   negative answer to "is evolution still proceeding robustly at that
   timescale" *at this regime* — worth knowing.

**Methodological caveat (load-bearing, must not be skipped).** A's
`dyn_excess_pc_slope` ≈ **+335 ± 317** (vs C's ≈ +1) is an
**artifact, not biology**. With `mu_lut=0` freezing a *rich evolved*
substrate, realised LUT flux → 0, so the S2d dyn fixed-space neutral
baseline degenerates toward 0 while the alive population keeps
exercising the rich frozen rule → observed/≈0 → excess explodes. The
causal-control v2 dyn-null (`dyn_excess_pc`=0 for EGENE_only) held
*only because it froze GoL* (a sparse substrate); freezing an
*evolved* substrate breaks that null. **`dyn_excess_pc` is
uninterpretable for any frozen-evolved-substrate arm** (A here; also
relevant to #8 and the ring-stability work). Use it only where
`mu_lut>0` (C: ≈0 — itself consistent with R1 that this regime is not
a strong rule-selection corner). New instance of the
[[shadow-models-lut-only-mismatch]] family — logged to memory.

**Honest bottom line.** #3c does *not* deliver a clean
"co-evolution synergy is real and strong." It shows the headline 3b
deficit was largely GoL **death**, leaves a **small but consistent**
genuine freezing penalty (C>A) that scales with substrate quality,
finds the egene signal **weak and plateaued** at this regime, and
surfaces a **new dyn-shadow artifact** for frozen-rich substrates.
The decisive clean co-evolution question now properly falls to **#8**
(direct lineage joint-retention — no freezing, no shadow, immune to
both the viability and the dyn-shadow confounds here).

**Status.** Done & recovered; notes/devlog/board written. Suggest a
follow-up board item: re-baseline the eg/dyn fixed-space shadows at a
mid-run switch (small C addition) before any future phase-switch
campaign leans on absolute excess.

---

## Causal-control v2 (local)  (`Scans/2026-05-17_causal_control_v2`)

Bug-corrected rerun of the 2026-05-16 v1 control (which had exposed
the shadow-scope bug, EGENE_only excess ≈ −26). 3 arms × 4 seeds,
N=256, 15000 ticks. Means over seeds (cog_spec, intake, n_distinct,
**excess_pc**, **eg_excess_pc**, **dyn_excess_pc**):

| arm | cog_spec | intake | n_distinct | excess_pc | eg_excess_pc | dyn_excess_pc |
|-----|---:|---:|---:|---:|---:|---:|
| JOINT      | 11.1 | 0.044 | 16182 | 0.0002 | **64.8** | **33.0** |
| LUT_only   | 25.0 | 0.055 |  9377 | 0.0001 | **0.000** | 17.3 |
| EGENE_only | 20.4 | 0.062 |   583 | −46.4  | 31.0 | **0.000** |

**Methodological validation (the headline).** The fixed-space shadows
pass their null controls *exactly*:
- LUT_only (egene frozen) → `eg_excess_pc` = **0.000** (no egene
  evolution ⇒ zero egene excess).
- EGENE_only (LUT frozen) → `dyn_excess_pc` = **0.000** (no rule
  evolution ⇒ zero effective-dynamics excess).
The neutral-null invariant S2d was built for now holds on the *actual*
causal controls at full scale. The metric is sound.

**Bug closed.** v1 EGENE_only excess ≈ −26 (artifact). v2: the
*fixed-space* `eg_excess_pc` = **+31** (egene evolution genuinely
beats neutral), while the *whole-genome* `excess_pc` = **−46** for the
same arm — i.e. whole-genome-hash excess remains actively pathological
for egene-driven runs, exactly the §2 prediction and precisely why
the fixed-space shadow was necessary. Use `eg_/dyn_excess_pc`; never
`excess_pc` for egene-driven evolution.

**Science.** Cognition is individually selectable (EGENE_only
`eg_excess_pc` = +31, LUT frozen). Effective rule dynamics are
individually selectable (LUT_only `dyn_excess_pc` = +17, egene
frozen). In JOINT *both* are larger than their isolated values
(eg 65 > 31; dyn 33 > 17) ⇒ **LUT and egene co-evolution mutually
amplify selection in both subsystems** — synergistic, not
independent. This corroborates the scan-scale #1-vs-#2 contrast with
a clean controlled experiment (two independent methods, same
conclusion).

---

## R1 — bracket the #1 unbracketed levers (m_scale, food_inc)  (`Scans/2026-05-18_R1_mscale_bracket`)

**Setup.** User-approved 2026-05-18, gated on #4. Reruns the #1 joint
design *identically* (300 cfg, 8000 ticks, N=256, same corrected
metrics/seed/init) with only two axes changed: `m_scale ∈
{2.5,3.5,5.0}` (bracket past #1's ceiling-pinned 2.5) and `food_inc ∈
{0.005,0.008,0.010,0.013,0.018}` (extend below #2's 0.010 low edge).
Rows therefore compare directly to #1. 300 rows, 0 errors, 743 s on
torque (≈2.5 s/cfg — heavier than #1's ~1.6 s, expected at higher
m_scale). Launched detached; survived a multi-minute VPN/torque
connectivity blackout with zero impact (the #4 detachment lesson
applied — `nohup`+`disown`, `.venv` python).

**Headline 1 — m_scale is NOT runaway; productive optimum is
interior ≈2.5–3.5.** Top-20 (corrected composite) m_scale frequency:
**2.5×10, 3.5×8, 5.0×2**. m_scale 5.0 is the *most viable* (see
Headline 3) yet barely appears among metric winners. So #1's
2.5-at-ceiling was the *lower edge of a productive plateau*, not a
true peak, and "more m_scale is better" is **disproved** — the lever
is now bracketed with an interior optimum near 3. eg-excess is also
markedly higher in R1 winners (top configs eg_excess ≈ 160–195 vs
#1's ≈ 115–135): the higher-m_scale regime amplifies egene-cognition
selection, consistent with #1's m_scale→cognition link, but
saturates/declines by 5.0.

**Headline 2 — food_inc optimum is HIGH; the low extension is
closed.** Top-20 food_inc: **0.018×10, 0.013×8, 0.010×2,
0.008/0.005 absent**. Sub-0.010 food_inc produces *zero* metric
winners even at m_scale 5.0 where it is survivable. #2's LUT-only
low-`food_inc` lean therefore **does not generalise to the joint
corrected objective** — same conditional-on-egene/objective pattern
already seen for the mu_lut optimum. The low edge is firmly bracketed
shut.

**Headline 3 — non-monotone U-shaped viability + a new
mu_lut×m_scale interaction.** Extinction (final_pop=0) is **107/300
(36 %)** vs #1's 16/300 (5 %), but structured:

| m_scale \ food_inc | 0.005 | 0.008 | 0.010 | 0.013 | 0.018 |
|---|---|---|---|---|---|
| **2.5** | 23/23 | 10/12 | 6/17 | 2/23 | 0/21 |
| **3.5** | 24/26 | 18/20 | 13/18 | 6/19 | 0/26 |
| **5.0** | 2/19  | 2/20  | 1/19  | 0/15  | 0/22  |

Extinction is *worst at intermediate* m_scale (2.5–3.5) under food
scarcity and *lowest at m_scale 5.0 across all food_inc* — a U-shaped
robustness landscape. Mechanism: intake `M=(m/25)·matches·F` must
clear `tax`; at m_scale 2.5–3.5 with low F intake falls below tax →
starvation collapse, while m_scale 5.0 clears tax even at scarce
food. The earlier scans never saw this (m_scale capped at 2.5,
food_inc never below 0.010). Consequently #1's **"low `mu_lut`
wins" optimum is conditional on m_scale ≤ 2.5**: R1 top-20 mu_lut
mode shifts to **0.15** (0.15×7, 0.03×6, 0.10×3, 0.01×2) vs #1's
low-skewed (0.01×5, 0.03×6). mu_lut×m_scale is a genuine
metaparameter interaction — exactly the class of relationship the
program targets.

**Through-line intact.** `excess_pc_slope` ≈ 0.000–0.001 in every top
config (both R1 and #1) — whole-genome-hash excess remains a
non-signal; `dyn_/eg_excess_pc` are the load-bearing probes. Fully
consistent with #1/#2/causal-control-v2/§2.

**Caveat (self-critique).** The composite used here is a z-normalised
sum of `dyn/eg_excess_pc_slope`, `n_distinct_genomes_temporal_std`,
`unique_top_genomes` — a faithful reconstruction of #1's ranking, not
the canonical `evoca_from_scan_top` call. Two top-20 rows (m_scale
2.5, mu_lut 0.15, `nDstd`≈14k–17k, modest fpop) look like
near-collapse population *thrashing* inflating `nDstd`, not healthy
diversity — the `nDstd` term can reward dying-population churn. The
load-bearing evidence is therefore the **param-frequency shifts and
the dyn/eg-excess-dominated rows (1,2,3,5,6,7: m_scale 2.5–3.5,
food_inc 0.013–0.018, mu_lut 0.03–0.15)**, which are robust to
ranking weights; the absolute composite ordering is not.

**Status.** Done & recovered. `results.csv`/`scan_config.json` pulled
to local; `notes.md` written. Resolves board R1.
