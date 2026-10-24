# CHEM-N6: layered fronts and excitation waves (2026-10-24)

**Question (N).** N saw "noisy fronts, almost spiral waves" at one-way exchange, `waste_tox` 0.4,
`waste_sink` 0.019 (N=256, interactive). The CHEM-N1 statistics had never measured that point. N
suggested detecting the layered fronts (the spiral arms) with a crude spatial power spectrum,
rather than tracking moving cores.

**Answer.** Yes, these are **excitation waves**, in the one-way medium at `waste_tox` 0.4 as well as
in two-way exchange.
- The space-time plots (`kymo_view.png`) show the textbook signatures: fronts leave sources (V
  shapes), travel at about 1-2 cells/step (by eye) and annihilate on collision (inverted V).
- Sites fire with one period, 115-120 steps, consistent across N = 256, 512 and 1024 and across
  seeds.
- They are not periodic in space. The fronts come from many sources, so the spatial spectrum has
  **no ring**, only a knee at the patch scale (about 100-200 cells). That is why the crude spatial
  spectrum, and the CHEM-N1 core statistics, did not see them.
- They are not spirals in any resolved sense. Spirals and target sources both give layered fronts,
  and nothing here separates them.

**This corrects CHEM-N5's "no propagation under one-way".** N5 ran at `waste_tox` 1.0, and at 1.0
the one-way control in N1 shows no period either. It is the lower toxicity that makes the one-way
medium excitable. The notebook's toxicity sweep ranked this point only by alive std, which
measures bulk fluctuation and is blind to waves.

## Runs (`scan.py`, torque, 24 runs, 0 errors)

All runs use recipe `2026-10-23_chem_fw_spiralpoint`, 6000 steps, 3 seeds:

- one-way (tox 0.4, sink 0.019): wdiff {0.01, 0.08} x N {512, 1024}, plus N=256 at wdiff 0.01;
- two-way at N=1024: (sink 0.019, wdiff 0.03, tox 0.6), (0.015, 0.08, 1.0) and (0.019, 0.01, 0.4).

Saved per run (`frames/`, not committed):

- 4 x 4 block means of alive and W every 20 steps over steps 5000-6000;
- the middle row every step over steps 5400-6000;
- the final colormode-6 frame.

## Statistics (`python/evoca_waves.py`; offline in `analyze.py`, so they can change without a re-run)

| Statistic | What it measures | Verdict on these runs |
|---|---|---|
| `ring_spectrum` contrast | how far a ring rises above the dip on its low-k side | no ring anywhere: contrast <= 1.6 (a real ring is > 10 on constructed fields) |
| `ring_spectrum` moving | ring power in the fluctuation about the time-mean field | wave runs 0.99 (moving); wdiff 0.08 runs 0.87-0.95 alive, 0.3-0.5 for W (partly frozen) |
| `site_period` | site-averaged temporal autocorrelation, 4-cell coarse-graining | **the discriminator**: wave runs 115-120 steps, peak 0.29-0.37; speckle runs none |
| `kw_spectrum` | power at w != 0, strongest moving mode | moving share 0.95-0.99 everywhere, noise included; the peak share is tiny (broadband): uninformative |

| configuration | period (alive) | ACF peak | ring contrast | look |
|---|---|---|---|---|
| one-way N=256, wdiff 0.01 | 116 | 0.33 | 1.07 | waves |
| one-way N=512, wdiff 0.01 | 117 | 0.36 | 1.46 | waves |
| one-way N=1024, wdiff 0.01 | 117 | 0.33 | 1.42 | waves |
| one-way N=512/1024, wdiff 0.08 | none | 0.00 | 1.1 | speckle, no structure |
| two-way N=1024 (0.015, 0.08, 1.0) | 120 | 0.37 | 1.62 | waves |
| two-way N=1024 (0.019, 0.03, 0.6) | 115 | 0.34 | 1.56 | waves |
| two-way N=1024 (0.019, 0.01, 0.4) | 87 | 0.29 | 1.37 | waves, shorter period |

**On the 39 saved CHEM-N1 runs** (`n1_fronts.csv`, `analyze_n1.py`, from the saved middle rows):

- the one-way control (tox 1.0) has no period;
- two-way runs have one in 50-100% of cases per cell, except wdiff 0.08 at tox 1.0 (20%);
- the ACF peak correlates 0.60 (Spearman) with N1's 4-period defect-survival excess, and 0.04 with
  Kuramoto R.

## Measures tried and dropped (kept here so they are not retried blindly)

- **Ring prominence alone.** It ranks large-scale variance: a reddened spectrum peaks at k = 1. Hence the contrast check, and a red-noise test.
- **A fixed-speed ridge** `|w| = c|k|` in the space-time spectrum. It found speed 2.0 on constructed irregular fronts, but on real runs it locked onto slow drift near w = 0 and scored the speckle runs highest. Whitening each frequency row drowned the signal in per-cell speckle.
- **A structure-tensor slope.** Frozen and moving were separated well. But the speed was biased by the noise structure (about 9c on constructed fronts with band-scale smoothing), and on real runs it was 1-5 cells/step from seed to seed.
- Consequence: **the front speed is not measured**. "1-2 cells/step" is by eye from `kymo_view.png`.

## Figures

- `kymo_view.png`: space-time plots, the evidence.
- `spectra.png`: radial spectra. The knee is there; the ring is not.
- `frames.png`: final frames, seed 0.

## Infrastructure note (torque)

The run took about 1 h 25 min of wall time for about 16 min of compute. 15 minutes in, torque went
into **"Thermal Emergency Sleep"** (pmset log, 00:37:38) and slept until 01:20. From outside that
looked like a network outage. The pmset log shows 15 such events since 2026-10-20, all preceded
by "Ignored DarkWake thermal emergency signal". Between jobs torque cycles between sleep and dark
wake; a `caffeinate -u` assertion made it fully reachable at once. Recommendations for N:

- launch scans under `caffeinate -is`;
- consider fewer workers;
- check vents and ambient temperature.

## Next

- Evolution in vs out of the wave regime (Track A): egene and LUT turnover, activity, lineage return.
- The front speed and spiral-versus-target remain open. A 2D phase from `site_period`-filtered block time series (Hilbert) could revisit cores at the right scale.
