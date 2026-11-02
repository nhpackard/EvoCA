# Closed-loop Gray-Scott regimes (2026-10-24)

For N. Pure numpy/scipy/numba; no EvoCA code touched. All scripts are in this folder, figures in `figs/`. Rerun order is by file number; whole study is about 25 min on a 16-core laptop.

## Summary

Chemistry: u + 2v -> 3v (k1 u v^2), v -> p (k2 v), p -> u (r p); u+v+p = m̄ conserved. Discrete step used throughout: reaction per channel with exponential integrator (u->v converts u(1-exp(-k1 v^2)), v->p converts v(1-exp(-k2)), p->u converts p(1-exp(-r))), then per-species blur x <- (1-d)x + d*mean3x3(x). The blur's effective diffusion constant is D = d/3 (verified in `test_gs.py`).

**Recommended parameter sets** (dt=1, lattice spacing 1, m̄ = 1; for other m̄ keep a = k1 m̄^2/k2 fixed, i.e. k1 -> k1/m̄^2 with u,v,p measured in units of m̄). Diffusion weights (du, dv, dp) = (0.6, 0.3, 0.3) in all four.

| Set | k1 | k2 | r | a = k1 m̄²/k2 | rho = r/k2 | What it does | Period / speed / wavelength | Min. lattice |
|---|---|---|---|---|---|---|---|---|
| S1 oscillatory | 14.15 | 0.1 | 0.0094 | 141.5 | 0.094 | Well-mixed relaxation oscillation (sharp v spikes). In 1D a seeded site becomes a pacemaker emitting pulse trains. In 2D: irregular wave/phase patches, no clean target pattern. | well-mixed period 117 steps (v 0.024 to 0.28); 1D emission period about 88; pulse speed 0.40 sites/step; wavelength about 35 sites | 128 |
| S2 excitable pulse | 8.68 | 0.1 | 0.0092 | 86.8 | 0.092 | One-way solitary pulses; medium returns to the resting state (u=m̄, v=p=0). Free ends retract: no spiral. | speed 0.32, FWHM 5 sites, recovery (min. pulse spacing) about 165 steps, so wake about 50 sites | 64 |
| S3 spiral | 6.0 | 0.1 | 0.005 | 60 | 0.05 | Excitable medium with long refractory time; broken front gives a persistent spiral pair. Survives at N=128, 256, 512. | spiral period 342-347 steps; wavelength about 46-48 sites; spiral wave speed lambda/T about 0.14 (solitary pulse 0.26); about 20-24% of sites active | 128 |
| S3b slower spiral | 3.0 | 0.1 | 0.003 | 60 | 0.03 | Same, larger and slower. Died at N=128 (too small); fine at 256, 512. | period about 553; wavelength about 76 | 256 |
| S4 invasion front (fallback) | 2.56 | 0.1 | 0.0346 | 25.6 | 0.346 | Front invading the resting state, leaving a uniform active state with damped oscillation. No sustained waves. | front speed 0.093 sites/step; plateau v 0.20, u 0.22, p 0.58; damped ring-down period about 66 | any |

Time-scale knob: for S3, k2 = 0.03 / 0.05 / 0.1 / 0.2 (a, rho fixed) gives spiral period 978 / 612 / 342 / 230 and wavelength about 106 / 80(coarse) / 48 / 39 (N=256). So period scales about 1/k2 and wavelength about k2^(-0.5 to -0.7). Keep k2 <= 0.1: at k2=0.1 the pulse is only 3-6 sites wide (lattice-scale); smaller k2 resolves it better.

Headline findings:

1. Mass conservation + 2 species: local dynamics are 1D, no oscillation (confirmed, Task 4). The third species (p) makes the reduced system 2D, and limit cycles do occur (Task 1). They live on a **thin ridge** in (a, rho) next to the Hopf locus, a_H ≈ 2.2 rho^-1.8 (rho <~ 0.2). This is a narrow window, not a broad region: evolving rates can leave it easily.
2. The trivial state (u=m̄, v=p=0) is **always stable** (eigenvalues -k2, -r): v = 0 is absorbing. Nothing happens without a seed of v. All waves here are excitation waves in a resting medium.
3. Excitable pulses and spirals are the **robust** regime (broad in a, rho; a ≳ 4 for any excursion, a from about 30 to 100 and rho from 0.03 to 0.1 for spirals). Spirals need a long refractory time (rho = r/k2 ≲ 0.07): the shortest sustainable wave period is about 1/rho x 1.7 steps/k2-scaled (S3: spiral period equals the minimum stimulation period, 340-370).
4. Diffusion: pulse speed grows about as sqrt(dv) and falls mildly with du; propagation **fails when du/dv ≳ 3-4 or dv ≲ 0.1 with du ≥ 0.2**. dp (waste diffusion) has no effect on 1D pulse speed or width (0 to 0.9). In 2D, dp = 0 doubles spiral period (about 600 vs 342) and wavelength. This is the opposite of open Gray-Scott, where Du > Dv is needed for spots.
5. **No spots / Turing patterns found.** Linear Turing instability of the active steady state needs Du/Dv about 25-60 (at Dp = Dv), outside the sensible range of box-blur weights; attempts at du/dv = 30 gave ragged invasion fronts and lattice-scale speckle, not spots.

## Method

- `gslib.py`: scaled ODE (time in 1/k2, concentration in m̄): du = -a u v^2 + rho p, dv = a u v^2 - v, dp = v - rho p. Two dimensionless numbers, **a = k1 m̄²/k2** and **rho = r/k2**, replace (k1, k2, r, m̄). Also continuous-time explicit-Euler and discrete-step integrators.
- Fixed points: v from (1 + 1/rho) v² - v + 1/a = 0 (m̄=1), u = 1/(a v), p = v/rho; exist iff a ≥ 4(1 + 1/rho) (saddle-node, i.e. k1 m̄² ≥ 4 k2 (1 + k2/r)). Trivial point always stable.
- `gs2d.py`: numba version of the discrete step (checked equal to the numpy version).
- `test_gs.py`: mass conservation to 1e-9, non-negativity, numba == numpy, D_eff = d/3, reachability (changing dv changes the result through the public step), discrete well-mixed vs ODE peak timing at small rates. All pass.
- `s0_open_gs_check.py`: open Gray-Scott, explicit Euler, dt = 1, Du = 0.16, Dv = 0.08. Reproduces spots (F,k = 0.035, 0.065), labyrinth/worms (0.037, 0.060 and 0.030, 0.057) and isolated solitons (0.014, 0.054) at t = 20000 (`f0_open_gs_check.png`). Integrator sanity check passed.
- `s1_wellmixed.py`, `s1b_cycles.py`: Jacobian regime map, kick experiments, limit-cycle amplitude/period (`f1*.png`).
- `s2_1d_scan.py`, `s2b_gallery.py`, `s3_sources.py`: 1D discrete scans and space-time gallery (`f2*.png`, `f3`).
- `s4_2d.py`, `s4b_spiral_scan.py`, `s4c_spiral_map.py`: 2D, broken-front seed (half-line kick with refractory p behind it), spiral survival map (`f4*.png`).
- `s5_diffusion.py`: Turing dispersion of the 3-species active fixed point (`f5a`). `s10_pulse_props.py`: clean 1D pulse measurement (first-arrival times at +150, +250 from a two-sided kick on a ring of 1000) and the (du, dv) map (`f10`).
- `s6_two_species.py`: Task 4. `s7_spots.py`: spot search. `s8_recommended.py`: min. stimulation period, spiral vs N, k2, diffusion. `s9_osc_discrete.py`: discrete well-mixed oscillation map and 2D oscillatory medium. `s11_figures.py`: final figures for the table.

## Results

### Task 1: well-mixed (f1, f1b, f1c, f9a)

- Regimes in (a, rho), upper (high-v) fixed point: none below a = 4(1+1/rho); unstable focus just above the saddle-node; stable focus (bistable with rest) above the Hopf curve; stable node for rho >~ 1.
- Hopf locus (trace = 0 on the upper branch): rho = 0.0236 -> a_H = 1925; 0.0677 -> 265; 0.194 -> 44.0; 0.557 -> 11.3 (merges with the saddle-node at rho about 0.56). Fit a_H ≈ 2.2 rho^-1.8 for rho <~ 0.2.
- Limit cycles (ODE, scaled): found only on the ridge between the saddle-node and a_H near a_H; amplitude of v up to about 0.35, period 7-13 /k2 (so 70-130 steps at k2 = 0.1). Discrete-time well-mixed version (k2 = 0.1) has the same ridge, a ≈ 56 (rho 0.17) to 1600 (rho 0.026), period 69-108 steps, amplitude from 0.001 to 0.17, very thin in a (a few %). Whether the Hopf is sub- or supercritical was not determined.
- Excitability: from the resting state, a v kick above a threshold gives a large excursion (peak v 0.3-0.95) and returns to rest for all a ≳ 4-10 where the active state is not a stable attractor (yellow area in `f1c`); threshold kick v0 falls from about 0.2 at a ~ 10 to 0.005 at a ~ 200. Where the active focus is stable (blue), the same kick ends in the active steady state (bistable, not excitable). Recovery time scale about 1/r.

### Task 2: spatial (f2, f3, f4, f9b)

- 1D gallery (`f2b`): single one-way pulses (a = 57.7, 87; rho = 0.02-0.09), pulse-train sources (a = 58-300, rho = 0.06-0.15), invasion fronts into a uniform active state (rho >= 0.3, a = 17-40).
- Pulse-train sources (`f3`, 12x12 scan): only on the Hopf ridge: (a, rho) = (449, 0.051) T = 79, speed 0.59, wl 46; (252, 0.069) T = 82, 0.49, 40; (141.5, 0.094) T = 88, 0.40, 35; (79, 0.128) T = 107, 0.31, 33. Source periods agree with well-mixed periods (70-120 steps).
- 2D (`f4`): at S1-like parameters the wake oscillates and the pattern becomes irregular wave turbulence; from a noisy uniform oscillating state (`f9b`) it becomes large patches in near synchrony with phase waves, never targets or spirals. At S2 a pulse passes and nothing remains.
- Spirals (`f4b`, `f4c`, `f11b`): survive at a = 30-100 with rho = 0.03-0.07 (survival fractions 0.08-0.3 of active sites); a = 40 gives ragged, fragmenting spirals; a >= 90 pulses leave no tip. At rho >= 0.15 the whole domain fills (uniform active). Spiral period and wavelength: see table. In `f11b` u and p show the spiral too, with u+v+p = 1 on average and 0.73-1.19 locally.
- Diffusion: see finding 4 and `f10`. Spiral (a = 60, rho = 0.05) vs (du, dv, dp) at N = 256: (0.6, 0.3, 0.3) T = 342; (0.45, 0.3, 0.3) T = 323; (0.3, 0.3, 0.3) T = 292; (0.9, 0.3, 0.3) T = 426; (0.6, 0.3, 0.0) T = 599; (0.6, 0.3, 0.9) T = 364. All still spiral. (0.9, 0.15, 0.3) gave a fine-scale irregular pattern (spectral wavelength about 10).

### Task 4: closed two-species analog (f6)

u + w = m̄, u + w -> 2w (k u w), w -> u (r w): dw/dt = w (k (m̄ - w) - r), a 1D flow with a single stable fixed point w* = m̄ - r/k (when k m̄ > r): monotone relaxation, no oscillation possible (the linearization has eigenvalues 0 for mass and -k w* for the other). Diffusion: the max real part of the eigenvalues of J - D q² over Du, Dw in [0.01, 100]² and q is <= 1e-6 (stable) and eigenvalues are never complex. Diffusion alone gives Fisher-type invasion fronts and then relaxation to the homogeneous state; with Du = 100 Dw a w-spike spreads out slowly (transient, relaxing over about 1e5 time units). General statement (mass-conserving two-component systems can show only stationary mass-redistribution patterns when the reactive-equilibrium curve slopes suitably): this rate law does not satisfy it.

## Caveats and what could not be found

- Only the **discrete-time** model was explored in space. Continuous-time spatial simulations were run only for the open Gray-Scott check and the two-species analog; no explicit-Euler stability analysis for the 3-species closed system was needed (stability limit for 5-point Euler is D dt <= 0.25 in 1D/2D with dx = 1; the box blur is unconditionally stable for d <= 1).
- At k2 = 0.1 the pulses are lattice-scale (FWHM 3-6 sites); results (especially speeds and the exact survival boundaries) depend on k2 and the splitting order of the reaction channels (1: autocatalysis, 2: decay, 3: regeneration, applied in sequence). I did not test other orders.
- Wavelengths come from radial spectra (±15%; values for S3 ranged 45-57 depending on bin choice); periods from probe peaks. Speed in 2D spirals is lambda/T.
- Only one broken-front seed type was tested; spirals do not form spontaneously from the resting state or from noise (v=0 absorbing, and an active uniform state is only reached for large rho). Spiral survival maps used T = 4000 only; slow decay beyond that is possible (S3 checked to t = 6000-8000 at N = 256, 512).
- The oscillation ridge is thin. In an agent model where rates mutate, S1 will be left quickly; S2/S3 (excitable) are far more tolerant.
- No Turing spots, no stationary patterns, no target patterns that persist in 2D, no chaotic regime was characterized (S1 in 2D looks irregular; not analyzed quantitatively).
- Sub- vs supercritical Hopf, existence of a limit cycle outside the thin ridge, and the exact refractory/dispersion relation were not mapped; the 'min. stimulation period' values come from three parameter sets only (S2: about 165; S3: about 340-370; S1 source: about 130 or less).
- The first-pass `s5_diffusion.py` pulse test and `s2_1d_scan` fate classifier were flawed (wrap-around and argmax tracking); the pulse part was removed from s5 and replaced by `s10_pulse_props.py`. `f2a_1d_scan.png` is a coarse classifier map; treat `f2b` and `f3` as the reliable 1D figures.
