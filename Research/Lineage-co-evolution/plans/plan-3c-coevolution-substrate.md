# Campaign #3c — egene evolution on a frozen *evolved* substrate

**Status:** spec finalized 2026-05-19 (user-approved). Harness to be
built on `main` (Python-only, no C change — runtime mutation setters
already exist). Queued for torque.

Supersedes the patch/S2c-transfer version (user dropped that 2026-05-18)
and the original "queued, uses S2c" board entry.

---

## 1. The question #3c must resolve

#3 found 3a (LUT co-evolving) ≫ 3b (LUT frozen) for egene evolution,
but **3b froze GoL**, not an evolved substrate, so the gap conflates
two hypotheses. #3c separates them with a single continuous sim:
evolve LUT-only, then either freeze or keep evolving the rule while
egene evolves.

| Pattern across arms | Conclusion |
|---|---|
| A ≈ C ≫ B | #3 gap was **substrate quality** (GoL is just bad); egene evolves fine on any good substrate, frozen or not. Weakens co-evolution-synergy. |
| C ≫ A ≈ B | **Freezing per se** kills egene evolvability regardless of substrate quality → genuine ongoing co-evolution synergy. |
| C ≫ A ≫ B | Both effects, graded. |

**A** = evolved-frozen · **C** = evolved-still-co-evolving (≈3a) ·
**B** = GoL-frozen (≈3b, anchor).

The **decisive, confound-free contrast is A vs C** from a *byte-
identical shared phase-1 substrate* (only "frozen?" differs). B/B′ are
calibration anchors to the prior #3/3b scale, not the primary
inference.

## 2. Resolved design decisions (rationale)

**D-1 Regime: single R1-informed regime, both phases.** *Not* the
top-10 #2 LUT-only winners. #2's LUT-only optimum wants low
`food_inc` (~0.010); R1 (2026-05-18) showed the **egene-driven**
objective wants `food_inc` **high (0.013–0.018)** and that #2's
low-food lean **does not generalize to egene-driven dynamics**.
Phase-2 is egene-driven, so #2-winner substrates would make arm A
lose to C *because the phase-2 environment is egene-hostile*, not
because freezing hurts — exactly the confound R1 predicts. Using one
R1-productive regime for both phases puts R1's finding in the core and
keeps A-vs-C clean. (User 2026-05-19: "R1-informed regime good.
alternative makes the water a bit muddy.")

**D-2 Phase-1 = 3000 ticks (not 1000/500).** Marginally-evolved
substrate ⇒ A ≈ B trivially ⇒ null *by construction*. Real LUT
structure took thousands of ticks in #1/#2. The just-merged
`get_lut_hashes()` (#8 enabler, M1) lets us **measure** substrate
evolvedness per run rather than assume it; shorten in a follow-up only
if the diagnostic saturates earlier.

**D-3 Both B and B′.** B = plain frozen GoL (direct 3b analog). B′ =
GoL pre-run 100 steps then frozen (settled-static GoL, addresses "GoL
dies to fixed points fast"). Both cheap; both included (user
2026-05-19: "both").

**D-4 Shared phase-1 via S1 determinism.** A and C must start from the
*same* evolved substrate. Free via the S1 deterministic RNG: run
phase-1 with seed *s* for A; re-run phase-1 with the *same* seed *s*
for C — bit-identical, asserted in-harness by `get_lut_hashes()`
equality. No snapshot/restore API needed. Cost: phase-1 computed twice
per (substrate,seed) A/C pair (modest).

## 3. Final spec

| Item | Value |
|---|---|
| Regime (both phases) | `m_scale=3.0, food_inc=0.015, tax=0.035, gdiff=0.06` |
| Phase-1 (LUT-only) | **3000** ticks, `mu_lut=0.03`, egene mutation off |
| Phase-2 | **8000** ticks, `mu_egene=0.01, mu_egenome=0.01, tax_per_egene=0` |
| Arm A (evolved-frozen) | phase-1 evolve → `mu_lut→0` → phase-2 egene |
| Arm C (evolved-co-evo) | phase-1 evolve → `mu_lut` stays 0.03 → phase-2 egene |
| Arm B (GoL-frozen) | frozen GoL, identical phase-2 |
| Arm B′ (settled GoL) | GoL run 100 ticks, then frozen, identical phase-2 |
| Diversity | 10 seeds × phase-1 `mu_lut ∈ {0.01, 0.03, 0.06}` |
| Budget | ≈ 10 × 3 × {A,C,B,B′} ≈ 120 runs × ~11k ticks; est. ~30–45 min torque (R1 ref 300×8k ≈ 740 s) |

## 4. Readout

**Primary:** `eg_excess_pc_slope` (A vs C vs B/B′), `dyn_excess_pc_slope`,
viability/extinction, n_distinct egene.

**Substrate-evolvedness diagnostic (per run, via `get_lut_hashes()`):**
fraction of alive cells with lut-hash ≠ GoL-hash; # distinct evolved
lut-hashes; phase-1 accumulated `dyn_excess`. Validates D-2's premise
that the substrate is genuinely evolved.

**Long-timescale ongoing-evolution robustness (user point 2 — first-
class):** 8000 phase-2 ticks is long; a single whole-phase slope can
hide a process that stalled early. So split phase-2 into early vs late
windows (e.g. ticks 0–2000 vs 6000–8000) and report, for each arm:
`eg_excess` slope late-window vs early-window; egene
turnover/`n_distinct` rate late vs early; whether selection is *still
active* at ticks ~7000–8000 (late-window slope significantly > 0).
This directly answers "is evolution still proceeding robustly at that
time scale," and additionally sharpens A-vs-C: if C's late-window
activity stays high while A's decays, ongoing co-evolution is what
sustains long-run evolvability.

**A-vs-C extra:** C's phase-2 LUT drift (how far the rule keeps moving
once unfrozen), via `get_lut_hashes()` distinct-count trajectory.

## 5. Implementation risks to pin in the harness

1. **S1 full-state determinism**, not just metric determinism: assert
   replayed phase-1 yields identical `get_lut_hashes()` before
   diverging A/C. If it fails, the shared-fork premise (D-4) is void.
2. **Shadow scoping:** the S2d fixed-space eg/dyn neutral shadows must
   be (re)baselined at the **phase-2 switch**, scoped to phase-2 only,
   so `eg_excess` isn't diluted/biased by LUT-only phase-1.
3. Comparability: phase-2 = 8000 matches #3's egene-phase scale so
   `eg_excess` numbers sit on the same axis as 3a/3b.

## 6. Build / run plan

- Harness `Scans/2026-05-19_3c_coevolution_substrate/scan.py`
  (single continuous sim per run; phase switch via
  `update_mu_lut`/`update_mu_egene`; deterministic phase-1 replay for
  the A/C fork; corrected-metric + diagnostic instrumentation).
- Pure Python on `main` (no C change). Launch detached on torque
  (`nohup`+`disown`, `.venv` python — the recovered #4/R1 lesson;
  SSH is the only valid torque health check, ICMP is filtered).
- Findings → `Docs/devlog/2026-05-17_research-campaigns.md` (#3c
  section), `notes.md`, board.
