# CHEM-N1: wave and spiral survey under two-way waste exchange (2026-10-23)

**Why.** Two-way exchange (CHEM-D4, `waste_exchange=1`) produced traveling phase waves and wave sources in one seed (`Experiments/2026-10-23_twoway_exchange`). Does it give spirals, and where?

**Design** (driver `scan.py`):
- Recipe `Runs/2026-10-23_chem_fw_spiralpoint.evoca`, with the standard rule, mutation on and a halfplane start; 6000 steps.
- **Grid** (N=512, 108 runs): `waste_sink` {0.013, 0.015, 0.017, 0.019} x `wdiff` {0.01, 0.03, 0.08} x `waste_tox` {1.0, 0.6, 0.4} x seeds 0-2.
- **Large** (N=1024, seed 0): `wdiff` 0.01 and 0.08 at `waste_sink` 0.017.
- **Control** (one-way, N=512, seeds 0-2): `waste_sink` 0.015, `wdiff` 0.01.
- Seed-0 runs save a space-time plot and the final frame (`kymo/`, gitignored).

**Diagnostics** (`python/evoca_waves.py`, 8x8 blocks):
- **Phase** relative to each block's time-mean state (steps 2000-3000).
- **Kuramoto R**: near 1 for synchronized bulk pulsing, low for waves.
- **Defect survival** at lags of 60, 120, 240 and 480 steps (half a period to four), against a half-lattice-shifted chance level. A defect survives if a same-sign defect lies within 1 block. The excess is (observed - chance) / (1 - chance).

**Why these measures** (smoke test, seed 0, N=512):
- The N5 test (60-step persistence ratio with a 2-block window) is saturated at this defect density (about 9% of plaquettes).
- A large-loop winding filter does not thin noise defects either (`evoca_waves` note; measured on random fields).

| Case | R | Excess at 1 period | Excess at 4 periods |
|---|---|---|---|
| two-way, `waste_sink` 0.017, `wdiff` 0.01 | 0.11 | 0.50 | 0.17 |
| two-way, 0.015, 0.08 | 0.12 | 0.51 | 0.22 |
| one-way control | 0.39 | 0.30 | 0.095 |

Survival peaks at a lag of one period: the defect pattern repeats each cycle.

**Results.** Pending (torque).
