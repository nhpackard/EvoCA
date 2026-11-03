# Porting N's EvoCA spiral recipe to ChemCA

For N, 2026-11-02. Goal (N): see the dynamics of `Research/Chem/Runs/2026-11-02_spiral_wave_ish.evoca`
(EvoCA food/waste, N = 512) in ChemCA, with appropriately ported parameters.

## Summary

- **The EvoCA pattern is ecological.** With all mutation off it keeps its period (142–160 steps), its
  look and its levels (alive 0.18 vs 0.27). So a chemistry-only model can in principle carry it.
- **The EvoCA pattern is a bulk oscillation with turbulent fronts, not persistent spirals.** Its
  phase field has 300+ defects per frame, and they persist no better than chance at any coarse-graining
  (8, 16, 32 cells). The arcs N sees are transient.
- **A parameter-only ChemCA port gives real spirals, at EvoCA's period, but not EvoCA's pattern.**
  Calibrated ports (s2) oscillate at period 148–175 (EvoCA ≈ 146–160) and show persistent spiral cores
  (s3: 46–59% of cores survive four periods, against 4–21% by chance). The movies (s4) show textbook
  rotating spirals for chem109.
- **But the ports cannot fill the lattice, because ChemCA has no food cap.** EvoCA regrows food to 1 and
  discards only the overflow. ChemCA's open mode has a proportional sink, so:
  - a weak sink (`sink_F` 0.004) lets uninvaded food reach about 3, and a cell entering it takes up
    a large waste burst and dies of toxicity, so rich food becomes a barrier: large regions are never
    invaded (movies);
  - a sink that holds F at 1 (`sink_F` 0.012) wastes most of the supply (it removes 0.012·F per
    site per step), and every surviving configuration (86 of 150) is static, with no oscillation.
- **Recommendation:** a small ChemCA extension, an open-mode source cap per compound (`cap_c`: the
  source stops at the cap and the excess is discarded and booked as outflow), default off. With
  `cap_F = 1`, `sink_F = 0`, EvoCA's food rule ports exactly. N to decide.

## 1. The two step rules, side by side

| | EvoCA (n_chem = 2) | ChemCA (fw2, open mode) |
|---|---|---|
| Food source | +food_inc per site, F capped at 1, overflow discarded | +inc_F per site, no cap; proportional sink `F·(1−e^−sink_F)` |
| Waste sink | constant: min(W, waste_sink) per site | proportional: `W·(1−e^−sink_W)` |
| Diffusion | `core_diffuse(gdiff)`, `(wdiff)` | same rule (`core_diffuse_precise`) |
| Uptake | mouthful = min(m/25 · matches · F, F, 1 − f): pattern-matched, capped by headroom | transport u_F·(E_F − a_F), no headroom cap |
| Energy store | f (private food is the energy) | e (energy), from converting store F → W |
| Waste made | as tax is paid: w += tax per step | as food is converted: w += converted F (all at once with fast k) |
| Cost | tax (+ per-egene) per step; death at f ≤ 0 | m_e per step; death when e < m_e |
| Excretion | one-way: waste_disp · w · (1 − W) | two-way: u_W · (w − W) (= EvoCA `waste_exchange=1`) |
| Toxicity | w ≥ waste_tox kills | store W ≥ tox_W kills |
| Reproduction | f ≥ 1, child to the lowest-f Moore neighbor (dead first), split | e ≥ 1, same rule on e |

Three differences are structural: the food cap, waste made at intake, and uptake without a headroom
cap. The first is the one that matters most (above).

## 2. The port

| Tier | EvoCA | ChemCA |
|---|---|---|
| Exact rule | gdiff 0.08, wdiff 0.01 | diff_F 0.08, diff_W 0.01 |
| Exact rule (with `waste_exchange=1`) | waste_disp 0.09 | u_W 0.09 |
| Measured | effective tax 0.050 (0.0301 + 0.0008 · 25; s1 nomut) | m_e 0.0501, c_g 0 |
| Same meaning | food_inc 0.012 | inc_F 0.012, light 0, open_mode 1 |
| Same meaning | f_init 0.5, F_init 1, halfplane | e_init 0.5, E_F 1, store 0, halfplane |
| Calibrated (s2) | eating (pattern, headroom) | u_F |
| Calibrated (s2) | waste_tox 0.792 (waste at tax rate) | tox_W (waste at intake: bursts) |
| Calibrated (s2) | F cap at 1 | sink_F (no equivalent: see Summary) |
| Calibrated (s2) | waste_sink 0.019 constant | sink_W ≈ 0.019 / W̄ = 0.1 (W̄ ≈ 0.19) |
| Off | mutation (LUT, egenes) | mu_rgene = mu_ractive = 0 |

## 3. Runs

- **s1** (`s1_evoca_characterize.py`, local): EvoCA recipe at N = 256 (base, nomut, exch1) and 512 (base),
  6000 steps. Late means:

  | run | alive | tax/cell | eat fraction | F | W | w (alive) | tox share of deaths | period |
  |---|---|---|---|---|---|---|---|---|
  | base N256 | 0.27 | 0.039 | 0.34 | 0.48 | 0.20 | 0.36 | 0.13 | 150 |
  | base N512 | 0.27 | 0.039 | 0.33 | 0.48 | 0.20 | 0.36 | 0.13 | 147 |
  | exch1 N256 | 0.26 | 0.040 | 0.35 | 0.50 | 0.18 | 0.39 | 0.16 | 142 |
  | nomut N256 | 0.18 | 0.050 | 0.56 | 0.58 | 0.19 | 0.38 | 0.14 | 144 |

  Evolution lowers the effective tax (0.050 → 0.039) by shedding egene positions; the pattern does not
  depend on it. Figures: `figs/s1_frames_*.png`, `figs/s1_ts_*.png`.
- **s2** (`s2_scan.py`, torque, 450 runs, N = 256, 5000 steps): u_F × tox_W × u_W × sink_W × sink_F.
  319 extinct. The direct numeric port (u_F small, tox_W 0.792) never reproduces: a cell drains only its
  own site and peaks at e ≈ 0.95. Ranked against the EvoCA target in `s2_report.py`; best eight in
  `figs/s2_top8.png`.
- **s3** (`s3_spirals.py`, torque, N = 512, 2 seeds each): spiral-core persistence. Matched fraction of
  cores after a lag, against chance (later frame shifted by half the lattice), 8-cell blocks, seed 0:

  | run | alive | F | W | cores/frame | 150 steps | chance | 600 steps | chance |
  |---|---|---|---|---|---|---|---|---|
  | evoca | 0.27 | 0.48 | 0.20 | 320 | 0.82 | 0.66 | 0.73 | 0.67 |
  | evoca nomut | 0.18 | 0.58 | 0.19 | 367 | 0.85 | 0.72 | 0.75 | 0.72 |
  | chem109 | 0.14 | 1.14 | 0.08 | 84 | 0.71 | 0.17 | 0.49 | 0.18 |
  | chem307 | 0.05 | 2.36 | 0.03 | 104 | 0.76 | 0.11 | 0.59 | 0.12 |
  | chem418 | 0.12 | 1.02 | 0.17 | 66 | 0.55 | 0.15 | 0.21 | 0.12 |
  | chem13 | 0.06 | 2.24 | 0.03 | 57 | 0.64 | 0.06 | 0.45 | 0.04 |

  Seed 1 agrees. At 16- and 32-cell blocks EvoCA's excess stays ≤ 0.06; ChemCA's stays +0.27 to +0.59 at
  16 cells.
- **s4** (`s4_movie.py`, local): side-by-side movies, EvoCA left, ChemCA right, both rendered as red =
  alive, green = F (0–1), blue = W (0–0.5); steps 2000–4000, every 4th step, 30 fps. Files
  `figs/s4_chem{109,307,418}_s0.mp4` (30–40 MB each, not in git); frame sheets
  `figs/s4_sheet_*.png`.

## 4. Critique (what an opposing expert would say)

1. **"A matching period is not a spiral."** The first ranking scored the period, which mostly reflects the
   bulk oscillation. → Added the core-persistence test (s3) before claiming spirals. It changed the
   conclusion about the target: EvoCA's arcs are not persistent spirals.
2. **"The port is a fit, not a port."** Four of the knobs (u_F, tox_W, sink_F, sink_W) have no numeric
   counterpart, and the scan picked them on one seed at N = 256. → The mapping is split into tiers (§2),
   and the candidates were re-run at N = 512 with two seeds (they hold).
3. **"The mismatch is in the model, not the parameters."** The population levels (alive 0.05–0.14 vs
   0.18–0.27, F 1–2.4 vs 0.5) and the uninvaded regions do not move with any calibrated knob, and every
   configuration that keeps F at 1 is static. → The recommendation is a model extension (food cap), not
   more scanning.
4. **"ChemCA's spirals are too clean."** EvoCA's eating depends on the local CA pattern, which adds noise
   at every site; ChemCA's uptake does not. The clean spirals are what the chemistry does without that
   noise, so the comparison isolates what the CA layer contributes: it breaks spirals into turbulence.

## 5. Next (N to decide)

1. **Food cap in ChemCA open mode** (`cap_c`, default 0 = off, bit-identical otherwise): then
   `inc_F = 0.012, cap_F = 1, sink_F = 0` is EvoCA's food rule exactly. Re-run s2 (fewer knobs) and s3.
2. If waste timing still matters after the cap: an option where maintenance (m_e) is paid by converting
   store F → W, as EvoCA's tax does, so waste is made at the tax rate instead of at intake.
3. A noise control: add per-site multiplicative noise to ChemCA uptake, to test whether noise alone turns
   the clean spirals into EvoCA's turbulence.

## 6. With ChemCA's food cap (2026-11-02, later)

ChemCA gained `cap_<c>` (commit `678695e`): `inc_F 0.012, cap_F 1, sink_F 0` is EvoCA's food rule
exactly. **s5** (`s5_capscan.py`, torque, 150 runs, N = 256): u_F × tox_W × u_W × sink_W.
40 extinct. Best by the same score (0.19, against 0.84 before the cap):

| | alive | F | W | period |
|---|---|---|---|---|
| EvoCA target | 0.22 | 0.53 | 0.19 | 146 |
| cap137 (u_F 1, tox_W 1.6, u_W 0.4, sink_W 0.05) | 0.20 | 0.52 | 0.21 | 147 |
| cap105 (u_F 0.75, tox_W 1.6, u_W 0.2, sink_W 0.05) | 0.22 | 0.52 | 0.22 | 144 |

Toxicity never binds in these (tox share 0, identical results for tox_W 1.6, 2.4, off): the
refractory phase is food depletion alone, where EvoCA has 13% toxic deaths.

Spiral test (s3, N = 512, 2 seeds), 8-cell blocks:

| run | cores/frame | 150 steps | chance | 600 steps | chance |
|---|---|---|---|---|---|
| evoca | 320 | 0.82 | 0.66 | 0.73 | 0.67 |
| cap105 | 641 | 0.94 | 0.90 | 0.91 | 0.90 |
| cap137 | 170 | 0.73 | 0.37 | 0.61 | 0.38 |

So the cap fixes the levels, and the lattice fills (movie `figs/s4_cap137_s0.mp4`, sheet
`figs/s4_sheet_cap137.png`: the right panel now looks like EvoCA). cap105 is EvoCA-like turbulence
(dense defects, no persistence); cap137 keeps some persistent cores (+0.23 over chance at four
periods; EvoCA +0.06). One ChemCA knob (u_F, waste export) moves the port between EvoCA's
turbulence and persistent spirals.

## 7. N's question: how to widen the spacing between fronts (s6)

`s6_wavelength.py`, torque, 32 runs of the recipe (mutation on), N = 512, 6000 steps: food_inc ×
gdiff × waste_sink. Scale = first zero of the row autocorrelation of F (patch half-size); no run has
a positive second peak, so there is no regular wave train, only a typical patch size.

| food_inc | gdiff 0.08: scale | period | alive |
|---|---|---|---|
| 0.012 | 32 | 147 | 0.27 |
| 0.008 | 38 | 200 | 0.16 |
| 0.006 | 44 | (236 at waste_sink 0.01) | 0.11 |
| 0.004 | 50 | — | 0.07 |

(waste_sink 0.019; 0.010 is similar except at food_inc 0.012, where it breaks the pattern.)

- **Lower food_inc widens the fronts' spacing** (32 → 50 cells from 0.012 to 0.004, roughly
  food_inc^-0.4) and lengthens the period (147 → 236), at the cost of population (0.27 → 0.07).
- **gdiff has a narrow window.** Only 0.08 (the recipe's value) gives fronts. At 0.02 every run
  goes extinct (in the initial boom and bust); at 0.3 and 1.0 the food is mixed so fast that no
  fronts form: a fine, uniform scatter of cells (scale 3–9).
- Figures: `figs/s6_sweep_wsink0.019.png`, `figs/s6_sweep_wsink0.01.png`.

**N's verdict (2026-11-02, from watching the runs): food_inc 0.012 (the recipe value) stays best.**
Lower food_inc does lengthen the wavelength (more green between fronts), but the fronts become
structured, with fingers/worms reaching into the food, which distorts the spiral/horseshoe dynamics.
Higher food_inc supports more live cells per front, so the fronts are fatter, and the fat front
averages the fine structure out. So front spacing and front smoothness trade off through food_inc;
widening the spacing needs a different knob that leaves front thickness alone (open).
