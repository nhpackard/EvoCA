# Campaign #8 — direct lineage co-evolution metric

**Status: ON HOLD (2026-05-19) — metric design flaw caught in
smoke-testing; user rethinking the observable.** Decisive clean test
of LUT↔egene co-evolution. Enabler `get_lut_hashes()` merged (M1,
`f8d6d66`). Lineage = S2a (merged). Pure Python on `main`
(analysis-only, no C change); run detached on torque. The
*machinery* (lineage capture + parent-child edge reconstruction +
permutation null) is built and validated in a local WIP harness; the
*observables* below (D1/D3 binary indicators) are wrong, see
"Design flaw" section.

---

## ⚠ Design flaw caught in smoke-testing (2026-05-19)

The binary indicators specified below (D1/D3) — `ΔLUT = "lut-hash
changed?"`, `Δegene = "egene-key changed?"` — are **degenerate at any
productive `mu_lut`**.

**Evidence.** Tiny-scale smoke (N=64, T_WARM=400, T_MEAS=1500):
- coevo arm (`mu_lut=0.03`): `p(ΔLUT) = 0.9993` — the 2×2 collapses
  (`n11=301, n10=669, n01=0, n00=1`); essentially every surviving
  edge has ΔLUT=1. `φ ≈ 0.02, z = 0.6, p = 0.73` — *trivially* near
  zero, not because co-evolution is absent.
- frozen-LUT arm (`mu_lut=0`): `p(ΔLUT) = 0.0` by construction.

**Root cause.** With `LUT_BITS = 250`, per-reproduction
`P(≥1 LUT bit flipped) = 1 − exp(−mu_lut · 250)`. That's **0.918 at
`mu_lut = 0.01`**, **0.9995 at `0.03`**, **1.0 at `≥0.06`**. The
smoke's empirical 0.9993 matches the analytic value exactly — so
"did the lut-hash change" has essentially **no variance** at any
regime where LUT actually evolves. φ is ~0 because the indicator
saturates, not because the population lacks (ΔLUT, Δegene)
association. Between a saturated coevo arm and a null frozen arm,
the binary contrast carries **no information either way**. The
devlog framing ("any correlation in surviving edges beyond the
mutation-rate product …") remains correct in *spirit*; the binary
indicator is the wrong observable for it on a 250-bit LUT.

**What is salvageable.** The harness machinery — lineage ON via
`init(..., lineage=True)` (NOT `enable_lineage()` separately —
`init()` resets the flag), warm-up + periodic snapshot of
`(birth_id, parent_id, lut_hash, egene_key, alive)`, incremental
`birth_id → genome` map exploiting genome-fixed-over-life, edge
reconstruction via `parent_id`, surviving-edge filter via observed
lifespan, within-run label-permutation null — is sound, validated,
and reusable. Only D1/D3 need revision.

---

## Possible fixes (for user consideration — no decision pinned)

Each keeps the devlog logic (independent at mutation step ⇒
association among *surviving* edges = joint retention under
selection) but on observables with variance.

**Option A — Magnitude (Hamming) correlation.** Replace binary
indicators with magnitudes: `ΔLUT = Hamming(parent.lut_bytes,
child.lut_bytes)` over the 32 LUT bytes (0–250 bit distance);
`Δegene = edit distance(parent.egene, child.egene)` (e.g. sum of
active-mask Hamming + per-active-slot value/mask Hamming). Spearman
correlation among surviving edges; permutation null on the egene
magnitudes vs LUT magnitudes. Full variance, directly implements the
devlog logic at the right granularity. **Prerequisite:** bulk
`get_lut` snapshots (32 B/cell × N² × n_snapshots ≈ 240 MB/run at
N=256, 120 snapshots) — heavier but feasible. Egene magnitudes:
already computable from `get_egenes` / `get_egenes_mask` /
`get_active`.

**Option B — Magnitude + coarse phenotype.** A as primary, plus a
secondary categorical `ΔLUT_phen = lut_complexity level (1/2/3)
changed?` — low-variance but phenotypically meaningful, and it links
the co-evolution metric to the **ring-stability** program. Two
readouts, cross-checked. (`lut_complexity` exists in C, exposed via
`evoca_lut_complexity_counts`; per-cell access may need a tiny
additive accessor analogous to `get_lut_hashes`.)

**Option C — Information-theoretic.** Mutual information
`I(ΔLUT_bits ; Δegene_bits)` on a per-edge basis, computed from
binned magnitudes; permutation null. Distribution-free; useful if
the LUT-egene coupling is non-monotonic.

**Option D — Generation-spanning joint persistence.** Reframe "joint
retention" as *joint persistence across surviving lineages*. Walk
parent_id chains back `k` generations; measure whether specific
(LUT-class, egene-key) combinations persist longer than neutral
drift would predict. Tests retention as a *temporal* property of
lineages, not a per-edge property. Heavier; needs longer runs.

**Option E — Conditional / phenotypic bit-flip rate.** Among
surviving edges, is the *effective* per-bit LUT-mutation rate
suppressed when the egene also changed (or vice versa)? Selection
filtering jointly-incompatible (rule-change, egene-change) pairs
would show as conditional rate suppression. Closer to a fitness-cost
read-out than a correlation; needs careful baseline.

**Option F — Phenotypically-active bits only.** Define ΔLUT on bits
that were actually *queried* this step (intersect with the
`restricted_mu` `lut_active` mask). Silent flips don't count;
reduces saturation if active footprint ≪ 250. Caveat: active
footprint at N=256 is often most of 250 → may not fully de-saturate;
useful as a refinement of A.

Likely a combination (A + B as cross-checks; A is the cleanest
direct analog of the devlog statistic). Open until you decide.

---

## Original spec (frozen — under revision per "Design flaw" above)

**Why #8 is now decisive.** #1/#2/causal-control = *indirect*
(excess-comparison). #3c = confounded (GoL-substrate death + a
dyn-shadow artifact, see [[dyn-shadow-frozen-rich-substrate]]). #8 is
**direct** and **shadow-free**: it reads the S2a lineage field, so it
is immune to both the viability confound and every fixed-space-shadow
pitfall.

## Idea

At the mutation step, ΔLUT and Δegene are drawn by **independent**
Poisson processes (separate passes, `C/evoca.c:1349/1373/1410`). So
across reproduction edges they are independent *by construction*. Any
(ΔLUT, Δegene) association that appears among **surviving** edges is
therefore selection jointly retaining co-adapted LUT+egene — the clean,
direct co-evolution signal, with no shadow model and no freezing.

## Definitions (pinned)

**D1 — egene-key.** Per cell: FNV-1a over `active_byte ‖ (value,mask)
bytes of active slots in slot order` (`get_active`/`get_egenes`/
`get_egenes_mask`). Captures the full ternary egene identity (mirrors
the C genome hash's egene fold); inactive slots excluded. ΔLUT uses
`get_lut_hashes()` (LUT-only, M1). ΔLUT=1 iff child lut-hash ≠ parent
lut-hash; Δegene=1 iff child egene-key ≠ parent egene-key.

**D2 — edge reconstruction.** An organism's genome is fixed birth→
death (only the *child* mutates, at its birth). Snapshot the lattice
every `S` ticks in a post-transient window; build
`birth_id → (lut_hash, egene_key, first_t, last_t, parent_id)`. Edge =
`parent_id → birth_id` usable iff *both* genomes were captured.
**Surviving edge** = child observed alive over a span ≥ `T_surv`.

**D3 — statistic (2×2 over surviving edges).** Contingency of
(ΔLUT ∈ {0,1}) × (Δegene ∈ {0,1}). Report, per run:
- `p11` observed joint-change frequency;
- independence expectation from the *surviving* marginals
  `p(ΔLUT)·p(Δegene)`; **excess = p11 − product**;
- **φ coefficient** and **odds ratio** (headline association);
- continuity cross-check vs the devlog's analytic *mutation-rate*
  product `P_mut(ΔLUT)=1−e^{−mu_lut·250}` (LUT marginal is exact);
- **label-permutation null**: shuffle Δegene across surviving edges,
  recompute φ over `N_PERM` perms → z / empirical p. This is the
  assumption-light significance test and it cancels survivor-capture
  bias (operates *within* the captured survivors).

**D4 — controls.**
- **Frozen-LUT arm** (`mu_lut=0`): ΔLUT≡0 ⇒ metric null *by
  construction* — the "3b cannot show it" check, here clean (no
  shadow, no viability confound).
- Per-run permutation null (above); multi-seed.

**D5 — regime.** R1-informed JOINT productive corner (continuity with
#3c): `m_scale=3.0, food_inc=0.015, tax=0.035, gdiff=0.06,
mu_lut=0.03, mu_egene=0.01, mu_egenome=0.01`, N=256. Warm-up
`T_WARM=4000`, measurement window `T_MEAS=6000`, snapshot `S=50`,
`T_surv=200`, `N_PERM=200`, seeds 0–7. Arms: {co-evo, frozen-LUT} ×
8 seeds = 16 runs.

## Readout / inference

- Co-evo arm: φ > 0 and permutation-z ≫ 0, excess > 0  ⇒ **direct
  evidence of real LUT↔egene co-evolution** (selection jointly retains
  co-adapted pairs beyond independent mutation).
- φ ≈ 0 / z ≈ 0 ⇒ no detectable joint retention at this regime
  (mutations independent *and* independently retained).
- Frozen-LUT arm: ΔLUT≡0 (sanity: confirms the metric is null when it
  must be) — and isolates that any co-evo arm signal needs a mutating
  rule.

## Caveats (surfaced)

- **Survivor-capture bias:** organisms living `< S` ticks are missed.
  The permutation null is computed *within captured survivors*, so the
  bias cancels for the significance test; marginals reported both
  empirically and (LUT) analytically for transparency. Lower `S` if
  capture fraction is low (reported per run).
- φ measures association, not direction of causation; #8 establishes
  *whether* joint retention exceeds independence, complementing the
  excess-based #1/#2/causal-control.
- No shadow models used ⇒ immune to
  [[shadow-models-lut-only-mismatch]] / [[dyn-shadow-frozen-rich-substrate]].

## Build

`Scans/2026-05-19_8_lineage_coevolution/scan.py` — lineage ON,
warm-up, snapshot loop, incremental `birth_id` map, 2×2 + φ + odds +
permutation null, per arm/seed. Pure Python on `main`. Smoke-test
locally (tiny scale; assert frozen arm ΔLUT≡0), then launch detached
on torque (nohup; SSH-only health check).
