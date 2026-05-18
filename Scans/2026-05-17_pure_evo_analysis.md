# EvoCA Pure-Evolution Campaign Analysis: corrected excess metrics, LUT↔egene co-evolution

**Date**: 2026-05-17 (living document — extended as #3/#4–#8 land)
**Scans**: `Scans/2026-05-17_pure_evo_joint/`,
`Scans/2026-05-17_pure_evo_lut_only/`,
`Scans/2026-05-17_pure_evo_egene_on_winners/` (in flight),
`Scans/2026-05-17_causal_control_v2/`
**Reasoning trail / per-campaign narrative**:
`Docs/devlog/2026-05-17_research-campaigns.md`
**Compute**: torque.protolife.com (32 cores) for scans; local for the
causal control.

## Question

Two coupled questions from `Docs/EvoCA_research_directions.md`
§optimization and §activity:

1. If we optimise **purely for evolutionary metrics** (no spatial
   terms), what regime is selected — and does it differ from the
   2026-04-27 balanced campaigns once the **excess-activity metric is
   corrected**?
2. Are LUT-rule evolution and egene (cognition) evolution
   **independent**, or do they **co-evolve** — and can we measure it
   cleanly rather than infer it?

The motivating prior: the 2026-04-27 campaigns optimised raw
ΣG−ΣN whole-genome excess activity, which a session-earlier audit
showed is **turnover-fragile and outright pathological for
egene-driven evolution** (the shadow modelled LUT-only while activity
hashed the full genome). The corrected metrics below are the
methodological core of this campaign.

## Methodology

Each scan: random sample from a discrete grid, run headless via
`evoca_explore.run_sim` on a `multiprocessing.Pool`, N=256, 8000
ticks (vs 5000 in 2026-04-27 — longer for a cleaner OEE read),
`sample_every=100`, fixed `seed=0`, neutral shadows on, left-halfplane
init. The `evoca_set_seed` (S1) fix makes each config reproducible.

| scan | configs | runtime (32c) |
|------|---:|---:|
| `pure_evo_joint` (#1) | 300 | 491 s |
| `pure_evo_lut_only` (#2) | 180 | 258 s |
| `pure_evo_egene_on_winners` (#3) | 288 | *in flight* |
| `causal_control_v2` (local) | 12 (3 arms × 4 seeds) | — |

**Candidate filter** (as 2026-04-27): not extinct AND time-mean
alive fraction ≤ 0.95.

### The corrected metrics (the key advance over 2026-04-27)

| metric | space | what it measures |
|--------|-------|------------------|
| `excess_pc_slope` | whole-LUT-genome FNV hash, open | component-normalised ΣG/D_G − ΣN/D_N (Channon-2006 correction, §2) |
| `dyn_excess_pc_slope` | 500 fixed (LUT-input,output) buckets | excess over a realised-flux Wright–Fisher drift baseline (S2d) — selection on **effective rule transitions** |
| `eg_excess_pc_slope` | 729 fixed ternary egene keys | same construction over egene space — selection on **cognition** |

The two fixed-space shadows have an exact neutral-null property
(zero realised mutation flux ⇒ excess ≈ 0), **validated on the
causal controls below**. Pure-evo composite = unweighted mean of
min-max-normalised, positive-clipped {`excess_pc_slope`,
`eg_excess_pc_slope` (joint only), `dyn_excess_pc_slope`,
`n_distinct_genomes_temporal_std`, `unique_top_genomes`}.

### Grid

Productive box from 2026-04-27 with the two unbracketed levers pushed
upward: `food_inc∈{0.010,0.013,0.018}`, `tax∈{0.030,0.035,0.040}`,
`gdiff∈{0.05,0.06,0.08}`, `m_scale∈{1.0,1.2,1.5,2.0,2.5}`,
`mu_lut∈{0.01,0.03,0.06,0.10,0.15}`. Joint adds
`mu_egene∈{0.003,0.01}`, `mu_egenome∈{0.003,0.01,0.03}`,
`tax_per_egene∈{0,8e-4}`. LUT-only freezes those three at 0.

## Campaign #1 — pure-evo, LUT + egene (300 cfg, 284 candidates)

Top configs (corrected composite):

| rk | f_inc | tax | gdiff | m_sc | mu_lut | exPC | egPC | dynPC | ndiv_sd | corrL |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.018 | 0.035 | 0.08 | 2.5 | 0.01 | 0.001 | 134.6 | 37.3 | 5717 | 21.1 |
| 2 | 0.010 | 0.030 | 0.08 | 1.2 | 0.03 | 0.001 | 91.9 | 60.5 | 3093 | 9.4 |
| 3 | 0.013 | 0.030 | 0.08 | 1.5 | 0.06 | 0.000 | 77.7 | 56.9 | 3798 | 9.0 |

**Finding.** `excess_pc_slope` ≈ **0** across *every* candidate;
`dyn_excess_pc` (~23–60) and `eg_excess_pc` (~40–135) strongly
positive. Selection is real and measurable in effective dynamics and
cognition, **invisible as whole-genome-hash excess**. The corrected
optimum is **low `mu_lut` (0.01) + high `m_scale` (2.5)** — the
*opposite* of the 2026-04-27 high-`mu_lut` optimum, which is hereby
attributable to the broken raw metric, not biology.

## Campaign #2 — pure-evo, LUT-only / egene frozen (180 cfg, 180 candidates)

**0/180 extinct** — with the fully-specified egene frozen, the box is
uniformly viable. Top configs want **high `mu_lut` (0.15, ceiling) +
high `m_scale` (2.0–2.5, ceiling) + low `food_inc` (0.010, low
edge)**; `dyn_excess_pc` positive but **markedly smaller** than #1
joint (~3–20 vs ~23–60). RD-scale spatial structure persists
(`corrL` ~55–75).

**Finding.** Freezing egene *reduces* effective-rule excess ⇒ egene
co-evolution **amplifies** LUT effective-dynamics selection. And the
`mu_lut` optimum is **conditional on egene**: egene off ⇒ high
`mu_lut` (raw rule exploration substitutes for adaptive cognition);
egene on (#1) ⇒ low `mu_lut` (cognition does the adaptive work). A
concrete metaparameter↔metaparameter interaction.

## Causal-control v2 (controlled corroboration, 3 arms × 4 seeds)

| arm | cog_spec | intake | n_distinct | excess_pc | eg_excess_pc | dyn_excess_pc |
|-----|---:|---:|---:|---:|---:|---:|
| JOINT      | 11.1 | 0.044 | 16182 | 0.0002 | **64.8** | **33.0** |
| LUT_only   | 25.0 | 0.055 |  9377 | 0.0001 | **0.000** | 17.3 |
| EGENE_only | 20.4 | 0.062 |   583 | −46.4 | 31.0 | **0.000** |

1. **Metric validated.** Frozen subsystem ⇒ its fixed-space excess is
   *exactly* 0 (LUT_only→eg=0; EGENE_only→dyn=0). The neutral-null
   invariant holds on the actual controls at full scale.
2. **Bug closed.** v1 EGENE_only excess ≈ −26 (artifact). v2
   fixed-space `eg_excess_pc` = +31 (real), while the old whole-genome
   `excess_pc` = −46 for the same arm — whole-genome excess remains
   pathological for egene-driven runs, exactly as §2 predicted.
3. **Synergy.** Each subsystem is individually selectable; in JOINT
   *both* exceed their isolated values (eg 65>31; dyn 33>17) ⇒ LUT
   and egene co-evolution **mutually amplify** selection. Independent
   method, same conclusion as #1-vs-#2.

## Cross-campaign synthesis

1. **The 2026-04-27 "mu_lut wants 0.06" optimum was a metric
   artifact.** Under correct normalisation whole-genome excess is a
   non-signal; the real selection signal lives in `dyn`/`eg` space.
   Re-rank any prior conclusion that leaned on raw excess.
2. **`mu_lut` has no context-free optimum** — it depends on whether
   egene co-evolves (low when it does, high when frozen). Future
   scans must report `mu_lut` *conditioned on the egene regime*.
3. **LUT↔egene co-evolution is real and synergistic** by two
   independent methods — but see the confound below; #3 is the
   decisive disambiguator and #8 the direct measurement.
4. **`m_scale` is still unbracketed** (peaks at the 2.5 ceiling in
   both #1 and #2) → board **R1**: follow-up with
   `m_scale∈{2.5,3.5,5.0}` (and a low-`food_inc` extension, from #2).

## Caveats / threats to validity

- **Frozen-substrate confound.** "Mutual amplification" compares
  JOINT against isolated arms whose frozen subsystem sits at a
  *specific* value (LUT=GoL, possibly maladaptive; egene=uniform
  fully-specified). Part of the amplification could be a handicap on
  the isolated arms, not true interaction. **#3's paired 3a
  (seeded-joint) vs 3b (frozen-LUT) on a LUT-*optimal* substrate is
  the clean test**; #8 (direct lineage ΔLUT×Δegene joint-retention vs
  neutral) is the direct measurement.
- **Run length.** 8000 ticks. The `excess_pc≈0` vs `dyn/eg≫0`
  contrast is a *within-run* differential (same length for all three)
  so it is genuine, but the absolute long-t behaviour of the metrics
  warrants a dedicated long run.
- **Single seed for scans.** #1/#2 use `seed=0` (now reproducible via
  S1). The causal control used 4 seeds and agrees, which partly
  covers this, but scan-level seed variance is unquantified.
- `unique_top_genomes` likely still sample-cap-bound (a 2026-04-27
  observation); it is one of five composite terms, not decisive here.

## Open / in flight

- **#3** (DONE): 3a≫3b (eg_excess 20 vs 6; 3b 46% extinct) — strong
  co-dependence. BUT 3b's `mu_lut=0` froze the LUT at *GoL* (run_sim
  always seeds gol), not at an evolved substrate, so #3 strengthens
  but does NOT cleanly resolve the confound. #8 promoted to decisive;
  #3c (evolve→extract evolved LUT→freeze→evolve egene, via S2c) queued.
- **#4** (queued): seed pure-evo from the 2026-04-27 balanced
  winners; track spatial `correlation_length`/`largest_patch` *over
  the run* — does RD structure collapse under pure-evo? (#1/#2 show
  `corrL` is *not* trivially destroyed, so the answer is likely "no,
  not automatically" — to be confirmed with time series.)
- **#5–#7**: genelife ring-ladder A/B; bifurcation sweep;
  patch-transfer 2×2.
- **#8** (auto after #3): direct lineage co-evolution metric.
- **R1**: bracket `m_scale` upward (+ low-`food_inc` extension).

## Recommended next steps

1. Let #3 land; read the **3a-vs-3b** contrast first — it is the
   single most decisive number for the co-evolution claim.
2. Run #8 (direct lineage measurement) regardless of #3 sign — it
   converts indirect inference into a measured statistic.
3. #4 with spatial time-series, then R1's `m_scale` bracket.
4. A targeted long-t run (≥50k ticks) at the #1 optimum to
   characterise the absolute (not just differential) behaviour of
   `dyn_/eg_excess_pc`.

## Reproducing

```bash
# scans live in Scans/2026-05-17_*/ ; each has scan.py + scan_config.json + results.csv
cd Scans/2026-05-17_pure_evo_joint && python3 scan.py     # (on a many-core box)
# ranking: see the analysis snippets in
#   Docs/devlog/2026-05-17_research-campaigns.md
```
