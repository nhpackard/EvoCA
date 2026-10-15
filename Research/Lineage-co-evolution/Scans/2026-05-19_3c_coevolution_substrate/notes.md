# 3c_coevolution_substrate (#3c)

Egene evolution on a frozen *evolved* substrate vs co-evolving vs
GoL-frozen anchors. Spec: `Docs/plan-3c-coevolution-substrate.md`.
80 rows, 0 errors, 650 s torque. **`determinism_ok=True` on all
30/30 A/C pairs** — the shared-fork crux (D-4) held exactly.
Premise confirmed: A/C substrate `frac_evolved=1.000`, ~38k distinct
LUTs at the switch (phase-1=3000 is amply evolved; D-2 vindicated);
B/B′ = GoL (frac 0.0).

**Finding (1-line, honest/mixed):** The #3 3a≫3b gap is **mostly a
GoL-substrate viability collapse, not a clean co-evolution synergy** —
plain frozen-GoL (B) goes **10/10 extinct** after a 3000-tick settle,
and viable frozen-GoL (B′) shows **flat-zero egene excess**. Holding
the substrate byte-identical, keeping the rule co-evolving (C) beats
freezing it (A) only **modestly but consistently** (paired C−A
eg_excess_pc>0 in 22/30; mean +0.6) — and the *absolute* egene signal
is weak at this regime (|slopes|≲1 vs #1's 20–135) and **plateaued by
ticks 7–8k** in every arm (answers the user's long-timescale question:
no, egene selection is not robustly still proceeding at 8000 here).

**Methodological caveat (load-bearing):** A's `dyn_excess_pc_slope`
≈ **+335** is an **artifact, not biology**: with `mu_lut=0` freezing a
*rich evolved* substrate, realised LUT flux → 0 so the S2d dyn neutral
baseline degenerates to ~0 while the population still exercises the
rich frozen rule → excess explodes. The causal-control dyn-null held
only because it froze *GoL* (sparse). **`dyn_excess_pc` is
uninterpretable for frozen-evolved-substrate arms.** New instance of
the [[shadow-models-lut-only-mismatch]] family.

Full analysis + adversarial caveats → `../../Docs/devlog/2026-05-17_research-campaigns.md` (#3c section).
