# Chem model v2: scan statistics, Pareto search, and a metabolism-based agent

For N, 2026-10-23. Deep analysis of N's message of 2026-10-23:

> Those parameters have noisy fronts. Almost spiral waves. Let's review the stats you are using
> to evaluate the scans. Whatever happened to the idea of replacing grid scans with a few
> generations of Pareto search? [...] Not sure the Vanilla EvoCA dynamics that we have developed
> are readily adaptable to chem. [seven points: drop the tax for a food→waste reaction; rate
> genes; separate code?; tension with toxic waste; mass conservation with an environmental
> waste→food reaction; catalysis (A,B,C,k); how to organize new code]

Nothing here is implemented. Section 13 lists the decisions for N.

---

## Status and corrections (2026-10-24, after N's reading)

- **Status: design only, Track B, later.** N: "let's not gallop ahead". Current work stays in vanilla EvoCA's food/waste mode (Track A): understand its dynamics and their evolutionary properties, including wgenes ([CHEM-N3](../ISSUES.md), which is *not* paused; §8e's recommendation is withdrawn). The one thing vanilla cannot do cleanly is a genetic metabolic rate (§4), which is why chemca exists at all.
- **Gray-Scott is agreed as the reference starting point for chemca**, with a correction to §0 point 4 and §6.
  - Gray-Scott has three species: U + 2V -> 3V, V -> P, with U fed and everything washed out. P is inert, so the dynamics are usually written with two variables.
  - The autocatalysis is **cubic**. My F + W -> 2W is quadratic, and that mostly gives simple invasion fronts. The threshold behavior (spots, pulses, spiral-like waves) needs the cubic term.
  - A closed-mass version **needs the third species**. With two species and conserved mass, each site's well-mixed chemistry is one-dimensional and cannot oscillate. The closed loop is: F + 2A -> 3A inside agents (yielding energy), A -> P inside agents (P exported), P -> F in the environment by light.
    - A is the activator; it is inherited at division.
    - P is the waste in the toxic sense.
    - Toxicity on P gives the refractory zone; on A it gives the "need some, die of too much" tension.
  - The first chemca is therefore three compounds, not two.
- **§2 (statistics) is superseded by [CHEM-N6](../ISSUES.md)** (`Research/Chem/Scans/2026-10-24_N6_layered_fronts/notes.md`).
  - Following N's suggestion, the layered fronts were measured instead of the cores.
  - The result: excitation waves, also one-way at tox 0.4. They are irregular in space, so they show no spectral ring, but each site fires with one period (115-120 steps). The site period is the working discriminator.
  - Core tracking (§2c) was not built.

---

## 0. Summary

**Scan statistics (section 2).**
- The point N watched (one-way exchange, `waste_tox` 0.4, `waste_sink` 0.019) was **never measured with the wave statistics**. [CHEM-N1](../ISSUES.md) ran one-way only as a control, at tox 1.0 / sink 0.015. The toxicity sweep in the notebook recorded alive std only.
- The [CHEM-N1](../ISSUES.md) spiral test is also biased against what N saw. It matches phase defects within 1 block (8 cells) at lags of 60-480 steps. A spiral tip that meanders, which real tips do, moves farther than that between frames. A noisy front with a meandering core therefore scores as "no spiral".
- R, alive std and the global period are bulk measures; they can rule out synchronized pulsing, but they cannot detect spirals.
- Proposal:
  - track defects frame to frame instead of a fixed-lag match;
  - take the local phase from each block's own time series (Hilbert transform);
  - add a phase-coherence length and the space-time spectrum;
  - calibrate every statistic against N's visual labels on a set of films before it ranks anything.

**Pareto generations (section 3).**
- What exists is post-hoc selection only. The cone/box selection was ported from pdt-server as `~/Projects/conepareto` (`select`, layers, regions, wedge-angle derivation). The generational sampler around it was explicitly not ported (report `Research/Activity-metrics/plans/2026-10-15_cone-pareto-package.md`, "Not ported").
- EvoCA uses `select` only to pick galleries from finished grids (Reaction-diffusion, Activity-metrics). [ACT-D2](../../Activity-metrics/ISSUES.md) and [RXD-N1](../../Reaction-diffusion/ISSUES.md) still wait on it.
- No EvoCA board item ever proposed replacing grids with generations. The idea lives only in the pdt-server comparison (`~/Projects/pdt-server/pdt-param-optimization.md`: "the thing PDT has that EvoCA lacks: sequential [...] active learning").
- Recommendation: build a small model-agnostic generational driver, but only after the statistics are calibrated (otherwise it optimizes a miscalibrated objective faster), and use it first on the v2 model.

**Chem v2 (sections 4-11).** N's instinct is right: the metabolism N describes does not fit EvoCA's energy bookkeeping. The core facts:
1. **The tax already is a food→waste reaction**, at a fixed rate (`g_metab_product`: the tax converts compound 0 into waste). What N proposes changing is the *rate law* (genetic, proportional to stock) and, implicitly, *what the reaction gives the agent*.
2. **That second point conflicts with pinned decision B** (`Research/Chem/plans/plan-chem.md` §2). Under it, the reproduction currency is stored free energy, `f = Σ energy · amount`. Converting food into energy-free waste then *lowers* `f`, so any rate gene evolves to zero. A reaction is only worth running if it *pays*, so the redesign needs a separate energy currency `e` that is produced by the reaction and spent on reproduction. Mass (food + waste) is conserved; energy is not, because it enters with the environmental waste→food regeneration (a "light" field) and leaves as heat.
3. **The tension N worries about is partly already there.** Stoichiometry fixes the waste made per unit of energy, so a rate gene cannot make "less waste per offspring"; it can only set the *tempo*. Fast cells win the eviction race and slow cells lose it, while fast cells poison themselves unless they export. The result is an interior optimum `k*`. But in a uniform world a single optimum means convergence and stasis, the vanilla result again. Sustained evolutionary activity needs `k*` to depend on the neighbors.
4. **The most natural way to make cells need waste** is product activation: the waste catalyzes the agent's own reaction, `F + W → 2W`. A cell then needs some waste to run its metabolism and dies from too much. Its optimum depends on the waste its neighbors export, which is the density dependence point 3 asks for. Two analogies support it:
   - glycolysis is activated by its own products (Sel'kov's 1968 model of glycolytic oscillations);
   - with feed and removal, `F + 2W → 3W` is the Gray-Scott system, a standard source of spots, waves and spiral-like patterns.

   Recommendation: make the catalyst a gene field (`C = none` or `C = W`), so that evolution chooses between plain and autocatalytic metabolism rather than having it imposed.
5. **Catalysis rate law: mass action**, `v = k · a · c`, integrated exponentially so that a step never overdraws a store.
   - A presence gate ("only if C present") is meaningless with continuous stores once any C diffuses in, unless a threshold is added.
   - Michaelis-Menten is derived for enzymes at trace concentration, which is not this case.
   - Rates are discrete (log-spaced index), so genomes hash and activity statistics stay meaningful.
   - Reaction genes `(A, B, C, k_idx)` reuse the egene slot machinery: 8 slots, an active mask, duplication on activation.
   - Per-compound free energies `G_c` make every reaction's yield `G_A − G_B`. That rules out free-energy cycles in larger networks.
6. **Code: a separate model in the same repository and library, not a branch and not flags in `evoca.c`.**
   - A new `C/chemca.c` with its own step function and genome.
   - A small shared core (PRNG, Poisson, diffusion, neighbor choice) extracted from `evoca.c`, guarded by the existing bit-exact goldens. This retires [CHEM-P1](../ISSUES.md) along the way.
   - A `ChemCA` Python class behind a common sim interface, with `"model"` in recipes (missing means `"evoca"`, a load-time default), and `build_sim` dispatching on it.
   - The Research tree stays on `main`.
   - The main cost is `controls.py`, which makes 56 direct `_lib.` calls and must move to sim methods. That is a layer-boundary move and a new pluggable point, so I flag it and suggest `/code-review ultra` there.
7. **Conflict to resolve:** [CHEM-D3](../ISSUES.md)/N3 (wgenes) is decided, but wgenes are CA-pattern genes for the current model. If v2 goes ahead without a CA, wgenes are superseded by transport genes. I suggest pausing [CHEM-N3](../ISSUES.md) until N decides.

---

## Outline

| § | Section | Purpose |
|---|---|---|
| 1 | Starting facts | What the current chem actually does, and the prior results the design must answer to |
| 2 | Scan statistics review | What each statistic measures, what it misses, why "almost spirals" went undetected, replacements |
| 3 | Pareto generations | History, status, why not yet, a driver sketch |
| 4 | What a reaction must give an agent | The energy-currency conflict with pinned decision B |
| 5 | Mass and energy accounting | Closed mass, open energy, the light field; what replaces `food_inc` and `waste_sink` |
| 6 | Why agents would need waste | Options for the tension, the tempo argument, the recommendation |
| 7 | Rate laws and catalysis | Gate vs mass action vs Michaelis-Menten; integration; energy levels |
| 8 | Genome and mutation | Reaction-gene format, discrete rates, gene cost, CA coupling, wgenes |
| 9 | Code organization | Options, recommendation, C and Python layout, recipes, controls, scans |
| 10 | Tests and guards | What must be bit-exact, conservation identities, reachability |
| 11 | Staged plan | Stages 0-4 with deliverables |
| 12 | Adversarial critique | Flaws found in the draft and the revisions they produced |
| 13 | Decisions for N | The calls only N can make |

---

## 1. Starting facts

**Current energy and mass flow** (`C/evoca.c` Phase 2c and 2d; recipe `2026-10-23_chem_fw_spiralpoint`):

- **Intake.** Eating moves `F → f` through the egene match. The mouthful is `(m_scale/25) · matches · F`, capped by `F` and by the headroom `1 − f`.
- **Tax.** The tax is `tax + tax_per_egene · positions` per alive cell per step (0.035 + 0.0008 · positions at the spiral point). It removes energy from `f`. When `n_chem = 2`, the amount removed becomes private waste `w` (`g_metab_product = 1`). **The tax is already a fixed-rate food→waste reaction.**
- **Reproduction.** A cell reproduces at `f ≥ 1`, and the child evicts the Moore neighbor with the lowest `f`.
- **Excretion.** One-way: `w → W` at `disp · w · max(0, 1−W)`. Two-way: `disp · (w − W)`.
- **Toxic death** at `w ≥ waste_tox`.
- **Environment.** `F` regenerates at `food_inc`, capped at 1. `W` has no source; it loses a uniform `waste_sink` per site and diffuses with `wdiff`.
- **Energy currency.** `f = Σ chem_energy · amount`, with food energy 1 and waste energy 0 (pinned decision B).

**Prior results the redesign must answer to:**

| Result | Source | What it implies for v2 |
|---|---|---|
| The balance regime is only 0.001-0.002 wide in `waste_sink` | [CHEM-N5](../ISSUES.md) | A per-site sink against density-proportional production is a fine-tuning problem. A closed mass budget changes its form (§5). |
| W piles up to 5-12 when the sink is slower than production, and acts as uniform aging | [CHEM-D2](../ISSUES.md), [CHEM-Q3](../ISSUES.md) | Unbounded W is an artifact of an open system without a regeneration path. Mass conservation bounds it (§5). |
| One-way exchange: local cycling, no propagation (soft refractory zone) | [CHEM-N5](../ISSUES.md) | The exchange law matters more than toxicity. v2 starts with two-way (passive) transport. |
| Two-way exchange: traveling waves with stationary sources over a broad region; no spirals resolved; phase noisy at 32-64 cells because colonies are speckle | [CHEM-N1](../ISSUES.md), [CHEM-Q4](../ISSUES.md) | The speckle plausibly comes from CA-pattern eating (patchy intake). A model without CA-gated intake is a direct test of that hypothesis (§8). |
| N's eye: one-way, tox 0.4, sink 0.019, N=256 shows noisy fronts and "almost spirals" | this message | The statistics must be checked against the eye (§2). |
| Vanilla EvoCA evolves to stasis (accepted); activity moves into chem | memory, 2026-10-23 | The v2 design must provide density- or frequency-dependent selection, or it reproduces stasis (§6). |

---

## 2. Scan statistics review

### 2a. What the scans compute now

[CHEM-N1](../ISSUES.md) (`Research/Chem/Scans/2026-10-23_N1_twoway_waves/scan.py`, module `python/evoca_waves.py`) uses 8×8 blocks and a burn-in to step 3000. [CHEM-N5](../ISSUES.md) and the toxicity sweep used only the first three statistics in the table below.

| Statistic | Definition | What it detects | What it misses |
|---|---|---|---|
| `alive_mean`, `alive_std` | global alive fraction over time | bulk density and bulk oscillation | All spatial structure. Traveling waves and spirals average out in a global mean, so they score *low* std, like a stationary state. The toxicity sweep's "coarsest mosaic" ranking used this, so it measured bulk fluctuation, not mosaic size. |
| `period`, `acf_peak` | first ACF peak of the global alive series | a bulk period | Waves: the global signal is weakly periodic, so the period is often nan where waves are clean. |
| `R_mean` (Kuramoto) | `|mean e^{iθ}|` over blocks with amplitude ≥ 0.1 | synchronization (R → 1) | Low R is shared by spirals, plane waves and noise, so R excludes bulk pulsing and nothing else. |
| `amp_mean` | mean distance of a block from its time-mean (F, W) | cycle strength | Nothing about organization. |
| `n_defects` | plaquettes with ±2π winding of θ | phase singularities | Dominated by noise pairs at speckle density. |
| `surv_lag`, `surv_lag_excess` | fraction of defects at t+lag with a same-sign defect within **1 block** at t, against a half-lattice-shifted chance | persistent cores | **Meandering cores** (below). |
| `toxic_share` | toxic deaths / all deaths | the refractory mechanism | — |

### 2b. Why "almost spirals" would score as "no spirals"

1. **Coverage gap, not detector failure, for the point N watched.**
   - [CHEM-N1](../ISSUES.md) sampled one-way exchange only at its control point (tox 1.0, sink 0.015, wdiff 0.01).
   - The tox 0.4 / sink 0.019 one-way point exists only in the notebook's toxicity sweep, which recorded alive std.
   - That point has never been through the wave statistics. Separately, "tox 0.4 destroys structure" was a two-way result and should not be read as covering one-way.
2. **Matching tolerance against tip motion.**
   - The wave speed is about 1.5 cells/step, and a spiral tip in an excitable medium meanders at a comparable speed (Winfree's meander).
   - At a lag of 60 steps a meandering tip moves tens of cells, several blocks, so a live core fails the 1-block test.
   - The long-lag survival test can only detect *stationary* cores, which is exactly the "stationary sources" N1 found (target sources, not spirals). It is blind to the noisy-front, wandering-core spiral N describes.
3. **Phase from a fixed center.** θ = atan2(W − W0, F − F0) about each block's time-mean over steps 2000-3000. If the pattern drifts after step 3000 (nonstationary sources, coarsening), the center is stale and θ is distorted. Noise in F and W at the 8-cell scale enters θ directly; there is no temporal smoothing.
4. **One block size.** b = 8 is below the 32-64-cell noise scale reported in [CHEM-Q4](../ISSUES.md). Defects at b = 8 are mostly speckle. b was never varied.
5. **One-dimensional film.** The kymograph is the middle row only. A spiral appears in a row cut as a source with alternating-direction waves, which looks much like a target source.

### 2c. Replacements

Ordered by value per cost; all are cheap on stored frames.

| New statistic | How | Distinguishes |
|---|---|---|
| **Local time-series phase** | Per block, the Hilbert transform of the band-passed W(t) (or alive(t)) around the dominant period. Standard in cardiac reentry phase mapping. | It does not depend on a global center and smooths over time, so θ is less noisy than atan2 on one frame. |
| **Tracked defects** | Frames every 10 steps; link each defect to the nearest same-sign defect within `v_max · Δt` (about 2 blocks); output trajectories and **lifetimes**; chance from the same linking on time-shuffled frames. | Spiral cores are long-lived trajectories, even when they meander. Noise defects live about one frame. This retires the fixed-lag survival test. |
| **Phase-coherence length** ξ | The spatial autocorrelation of `e^{iθ}` and the distance where it falls to 1/e. | The quantitative form of [CHEM-Q4](../ISSUES.md). A spiral needs ξ much larger than the core size; noise has ξ of about one block. |
| **Radial structure factor S(k)** of W or alive | 2D FFT, radially averaged; peak position and width | The pattern wavelength, and mosaic coherence (a narrow peak) against speckle (a broad one). It replaces alive std as the "mosaic" measure. |
| **Space-time spectrum** | FFT over (x, y, t) of stored frames | Bulk oscillation (a peak at k = 0), traveling waves (a ridge ω = c·k), and stationary patterns (ω = 0) are separated cleanly. |
| **Front roughness** | Interface width of the alive/dead or high/low-W front against front length | Measures "noisy fronts" directly, which is the visual feature N named. |
| **Scale scan** | Repeat the defect and ξ computation at b = 8, 16, 32 | Topological cores persist across b; noise defects vanish as b grows. |

### 2d. Calibrate against the eye before ranking

Global rule: a plausible result is the moment to check against something external. Here the external instrument is N's eye. The procedure:

1. Pick about 30 runs spanning the [CHEM-N1](../ISSUES.md) grid, [CHEM-N5](../ISSUES.md) and the tox 0.4 / sink 0.019 one-way point. Save 300-step films (colormode 6) at N=512.
2. N labels each film: bulk / stationary / target / spiral-like / noise.
3. Compute every statistic in §2c on the same runs.
4. Keep the statistics that separate the labels, and set thresholds there.
5. Only then use them in a scan or as Pareto objectives.

Cost: about 30 runs (minutes on torque) plus about 30 minutes of N's viewing. The film viewer can be a gallery `.show()` per run, which already exists.

**Immediate cheap check (offered, not run):** run the N1 diagnostics plus tracked defects on the one-way tox 0.4 / sink 0.019 point, 3 seeds at N=512, and at wdiff 0.01 and 0.08. This closes the coverage gap in §2b.1.

---

## 3. Pareto generations

### 3a. History

| Date | Event |
|---|---|
| 2026-05-21 | `pdt-server/pdt-param-optimization.md` compares PDT's generation loop (explore/exploit sampling, refit per generation) with EvoCA's one-shot grids: "the thing PDT has that EvoCA lacks: sequential, surrogate-driven active learning." |
| Oct 2026 | `prompt1-Oct2026.md`: N asks to abstract the restricted (cone) Pareto approach from pdt-server for any dual-parameter optimization. |
| Oct 2026 | Review `Docs/Archive/Oct2026/00-review-of-prompt1.md` §6-§8. EvoCA's "optimization" is scalarization over a static grid, with weaknesses: compensation, concave regions unreachable, pool-relative normalization, silent weighting. It proposes a package. |
| Oct 2026 | `~/Projects/conepareto` built, with R parity (`Research/Activity-metrics/plans/2026-10-15_cone-pareto-package.md`): `select`, layers, regions, hypervolume, and later wedge-angle derivation and `on_shortfall`. **"Not ported: [...] the lambda-sweep sampler. They sit around `select` rather than inside it; say if you want them."** No answer was recorded. |
| Oct 2026 | Used post-hoc: `rd_pareto_picks` / `rd_pareto_layers` galleries (Reaction-diffusion), and Activity-metrics notebook 4a. |
| open | [ACT-D2](../../Activity-metrics/ISSUES.md) ("adopt Pareto selection for (dyn, eg)?", recommendation yes, box region) and [RXD-N1](../../Reaction-diffusion/ISSUES.md) (re-select the 04-27 winners with box/wedge) are still open on the board. |

So the generational search was never built. The selection half exists; the sampling half does not. I found no board item for it in EvoCA. It fell between the pdt comparison and the port's "not ported" list.

### 3b. Why not run it on the current chem scans right away

- **The objective is not calibrated** (§2). A generational search climbs whatever statistic it is given. With the current spiral statistic it would climb toward stationary target sources, the thing that statistic actually rewards.
- **Single-seed noise.** pdt-server's own critique (§2b of the comparison) applies: an arg-max over noisy evaluations selects the lucky ones. The [CHEM-N1](../ISSUES.md) seeds disagree enough (defect excess 0.15-0.21) to matter.
- **The parameter space is about to change.** If v2 is adopted, the v1 parameters (`waste_sink`, `food_inc`) disappear (§5). A driver written generically survives that; a v1-specific campaign does not.

### 3c. Driver sketch (model-agnostic)

`python/evoca_search.py`:

- `search(run_fn, space, objectives, region, gens, per_gen, seeds)`.
- `space`: per-parameter bounds and scale (log or linear); sampling happens in normalized coordinates.
- **Gen 0:** a Latin hypercube of 60-100 points.
- **Each generation:**
  - evaluate every point with 2 seeds;
  - combine the seed values per objective as the mean minus half the seed spread, a cheap robustness penalty;
  - `conepareto.select(s1, s2, region='box' | 'qbox' | 'wedge', k)`;
  - spawn `per_gen` children around the selected points (Gaussian in normalized space, width halved each generation, plus 20% fresh LHS points against premature collapse);
  - **re-evaluate the survivors with fresh seeds** each generation, so a lucky survivor cannot persist.
- Outputs: one `results.csv` per generation plus a lineage column, so the existing `Scan` and `Gallery` tooling reads each generation as a scan.
- Budget: 5 generations × 50 points × 2 seeds = 500 runs, about the cost of four [CHEM-N1](../ISSUES.md) grids, over a 5-8-parameter space a grid cannot cover.
- Two objectives must be named first. For chem, candidates are (pattern coherence ξ or tracked-core lifetime, persistence or viability), and in v2 (pattern, evolutionary activity on reaction genes).

The surrogate model PDT uses is left out on purpose: EvoCA evaluates the true objective cheaply enough, and §1 of the pdt comparison argues that the surrogate is where PDT's weakness lies.

---

## 4. What a reaction must give an agent

N's first point is to eliminate the metabolic tax and replace it with a food→waste reaction. As §1 shows, under `n_chem = 2` the tax *is* that reaction at a fixed rate. The change N is proposing therefore has two parts:

1. **The rate law.** The rate is proportional to the stock and genetically set (`v = k · f`), rather than a constant draw.
2. **The payoff.** Today the reaction is pure loss. It lowers `f`, which is the reproduction currency, and buys nothing.

Under pinned decision B (`f = Σ energy_c · amount_c`, food energy 1, waste 0), any reaction F→W lowers `f`. Selection on a heritable `k` then drives `k → 0`: cells that metabolize less keep more `f` and reproduce sooner. **N's own design rules this out**, so decision B has to change. The options:

| Option | Currency for reproduction | What the reaction does | Verdict |
|---|---|---|---|
| (i) Keep B | stored free energy `f` | destroys currency | k → 0. Rejected. |
| (ii) **Energy currency `e`** | `e`, produced by reactions | `v` units of A→B yield `v · (G_A − G_B)` of `e` | **Recommended.** It is ATP-like: food is stored potential, metabolism realizes it, growth spends it. |
| (iii) Biomass compound | a body compound `B` built from F | F → B + W | More realistic (growth needs mass). It adds a third compound and a stoichiometry choice now, so defer it to a later stage. |

With (ii):
- A cell's state is its compound stores `a_c`, its energy `e`, its genome, and `alive`.
- **Reproduction** at `e ≥ 1`: split `e` and every store 50/50.
- **Eviction** targets the neighbor with the lowest `e`, not the lowest food.
- **Death:** toxicity (`w ≥ tox`), eviction, and, only if a maintenance cost `m_e > 0` is enabled, `e < 0`.
- On any death the stores return to the site. The cell's `e` is lost; it is heat, not mass.

**On "eliminate the tax".** A *mass* tax (removing food as a fixed cost) goes. An *energy* maintenance cost is a separate thing, and physically ordinary (basal metabolism). It has one job a reaction cannot do: removing cells that have stopped metabolizing. Without it, a cell with an empty food store and `e < 1` never dies unless evicted. In an empty region nothing evicts it, so inert "zombies" occupy space indefinitely. Recommendation: `m_e` as a metaparameter, default 0 in the first experiments, and scan it once to see whether zombies matter. The per-gene cost in §8 is also an energy maintenance cost.

---

## 5. Mass and energy accounting

N's point 5: conservation of mass for food + waste, with an environmental waste→food reaction complementing the agents' food→waste.

**Design.**
- **Mass:** `M = Σ_x (F + W) + Σ_alive (f + w)` is exact and conserved every step. It is set by the initial state, and the control parameter is the mean mass per site `m̄ = M / N²`.
- **Environment reaction** W → F at rate `r(x)` per step, the "light", integrated as `ΔF = W · (1 − e^{−r(x)})`.
  - It is the only energy input to the world.
  - `r(x)` is a field: uniform, or patterned like `env_mask` (gradients for the checkpoint §B structured-gradient experiment).
- **Energy** enters via `r`, is carried as F's potential `G_F`, is realized by agents into `e`, and leaves as dissipation (reproduction overhead, death, maintenance). Energy is booked as a flux, not conserved.
- **Gone:** `food_inc` (replaced by `r` acting on the existing W), `waste_sink` (W is not destroyed, only recycled), and the F ≤ 1 regeneration cap. `F` can still be capped at 1 per site by limiting the regeneration to `min(ΔF, 1 − F)`, so that a site does not hold more food than a cell can eat.

**What this does to the [CHEM-N5](../ISSUES.md) and [CHEM-D2](../ISSUES.md) problems** (a hypothesis to test, not a result):

- In v1, production scales with density, `ρ · tax`, while the sink is a fixed per-site removal `waste_sink`. Balance therefore needs `waste_sink ≈ ρ* · tax`, and nearby values either pile W up (aging regime, W of 5-12) or starve the field. That is the 0.001-wide band.
- In v2, total W is bounded by `M`, so W cannot exceed `m̄` on average. The pile-up of [CHEM-D2](../ISSUES.md) cannot happen in the mean, though local accumulation can.
- Consumption and regeneration feed back on each other: agents die when W is high, stop producing, and `r` restores F. That loop is a relaxation oscillator by construction, the predator-prey-like cycle that the v1 balance regime only reached by tuning.
- **Fine-tuning moves, it does not vanish:** the regime is set by `(m̄, r)` against the agents' rates. The rates evolve, though, so the population can tune itself to the environment, which v1's fixed tax could not.
- **Test:** a 2D scan over `(m̄, r)` at fixed genes, looking for the width of the oscillating regime, against the v1 band.

**Keep an "open" mode** (`food_inc`, `waste_sink`, `r = 0`) as an option so that v2 can be compared with the v1 results directly.

---

## 6. Why agents would need waste

N's point 4: if waste is toxic, agents will evolve rates that make little of it; there must be tension that makes agents need some waste.

### 6a. The tempo argument (what is already there)

- With one reaction F → W, every unit of food converted yields `G_F − G_W` of energy and exactly one unit of waste. **Waste per unit of energy, and so per offspring, is fixed by stoichiometry.** A rate gene cannot evolve "less waste per offspring".
- What `k` sets is **tempo**: energy per unit time, and therefore reproduction rate, and with it waste production per unit time.
  - High `k`: faster reproduction, wins the eviction race, but `w` builds up faster than export can remove it. The cell reaches toxicity, or empties its food store and has to wait on uptake.
  - Low `k`: safe, but out-reproduced and evicted.
- So there *is* tension, between tempo and toxicity, with an interior optimum `k*(export, uptake, local W)`.
- **The problem:** in a uniform world `k*` is one number. The population converges to it and stops, which is stasis again. For sustained activity, `k*` has to depend on what neighbors are doing (frequency- or density-dependent selection).

### 6b. Mechanisms that make waste *needed*

| Mechanism | Form | Per-agent tension? | Density dependence | Naturalness | Cost |
|---|---|---|---|---|---|
| **A. Product activation / autocatalysis** | F + W → 2W at rate `k · f · w` | Yes: some `w` is needed to metabolize, too much kills, so an interior optimum `w*` | Strong: neighbors' exported W primes metabolism, so a colony's success depends on local W | High: glycolysis has product activation (Sel'kov 1968 glycolytic oscillator); Gray-Scott is this reaction with feed and removal | One gene field (`C`); no new compound |
| B. Cross-feeding (producer niche) | some agents run W → F with light | Only for producers | Strong: two niches, ecology | High, but needs light inside agents (an "autotroph gene") | A second reaction and a light coupling for agents |
| C. Waste as building material | growth needs W | Yes | Weak | Low (ad hoc) | A biomass compound |
| D. Collective dependence only | closed mass: today's waste is tomorrow's food | **No** (a common good, no individual payoff) | Global | High | None (comes with §5) |

**Recommendation: A, offered as a gene choice rather than imposed.**
- The reaction gene's catalyst field `C` can be `none` (plain F → W, rate `k · f`) or `W` (autocatalytic, rate `k · f · w`). A mutation can switch it.
- The first experiment then asks whether evolution picks autocatalytic metabolism, and where. This is a cleaner question than building it in.
- The plain variant is the control.
- B arrives naturally in stage 3, once there are more compounds and reactions. D comes for free with §5, but it does not create individual-level selection; spatial structure could turn it into a group effect, which is worth watching but not something to rely on.

### 6c. What A predicts (for checking stage 2 against)

- **Ignition.** A newborn inherits half its parent's `w`, a metabolic inheritance of catalyst. A cell with `w = 0` cannot start unless it takes up W from its site. Colonies therefore spread into waste-primed ground and stall in clean ground, which is a front-propagation mechanism with a built-in excitability threshold.
- **Refractoriness.** Behind a front, W is high (toxic) and F depleted. Regeneration `r` converts W back to F, so the refractory time is about `1/r`. The v1 refractory zone was "soft" ([CHEM-N5](../ISSUES.md)); here the refractory variable is the same W that ignites, which is the Gray-Scott structure in which waves and spiral-like patterns are known to occur.
- **Selection on k** is density dependent: the best `k` is lower in high-W regions (toxicity close) and higher in fresh ground. A traveling front could carry a different `k` distribution than its wake, which is the coupling between pattern and evolution that the activity direction wants.

---

## 7. Rate laws and catalysis

N's point 6: beyond two compounds, catalysis. Michaelis-Menten, or simply A→B with rate k only if C is present? The gene would be (A, B, C, k).

| Law | Rate per step | Pros | Cons |
|---|---|---|---|
| Presence gate | `k · a` if `c > 0` | Simplest to state | With continuous stores, "present" is true almost everywhere once C diffuses or is inherited, so the gate stops discriminating. A threshold `θ` fixes that but adds a free parameter and a discontinuity. |
| Threshold gate | `k · a` if `c ≥ θ` | Switch-like, CA-flavored | Extra parameter; the rate is insensitive to how much C there is above `θ`; discontinuous fitness landscape. |
| **Mass action** | `k · a · c` | One parameter; smooth; autocatalysis is the special case C = B; it is the FKP-style network kinetics | Unbounded in `c`, but stores are bounded (toxicity, capacity) |
| Michaelis-Menten | `k · c · a / (K + a)` | Saturation in substrate | MM assumes an enzyme at trace concentration with fast complex formation. Here the catalyst is a compound in amounts comparable to the substrate, so MM is not the derived law. The substrate is already capped by store capacity, which saturates anyway. It adds `K`. |

**Recommendation: mass action, `v = k · a · c`, with an exponential integrator:**

```
converted = a · (1 − exp(−k · c))        # c = 1 for an uncatalyzed gene
```

- A step can never convert more than the store `a`, so no clamps are needed. The current exchange code needs explicit `d > a` clamps for exactly this reason.
- For small `k · c` it reduces to `k · a · c`.
- With several genes using the same substrate, the per-step rates sum, and the converted amount is split in proportion to each gene's rate, so that the order of the genes does not matter.
- Michaelis-Menten can be added later as an optional rate law if store sizes become large.

**Energy levels.**
- Each compound has a free energy `G_c`, set by the world. A reaction A → B yields `G_A − G_B` of energy to the cell.
- A reaction with `G_B > G_A` costs energy: it runs only while `e` lasts, and it is the hook for anabolic reactions (making catalysts) later.
- Because the yield is a difference of potentials, any closed reaction cycle yields exactly zero. Free-energy loops are impossible however large the network grows, which an arbitrary per-reaction yield could not guarantee. For two compounds, `G_F = 1` and `G_W = 0`, matching today's energy weights, so the units carry over.

**Reaction direction.** In stage 1, F → W only. A W → F gene inside an agent is possible in the representation but costs energy (`G_F > G_W`) and buys nothing, unless agents can harvest light (§6b, mechanism B, stage 3).

---

## 8. Genome and mutation

### 8a. Reaction genes ("rgenes")

- **Slot format:** `(A, B, C, k_idx)` in 4 bytes. `A, B ∈ compounds`; `C ∈ {none} ∪ compounds`; `k_idx ∈ 0..31` with `k = k_min · 2^(k_idx/4)`, about 2.4 decades in quarter-octave steps.
- **Storage:** 8 slots plus an active mask, **exactly the egene layout** (`NEGENOME_MAX = 8`, `active` byte, dup-on-activate). The machinery is reused, not reinvented.
- **Mutation:**
  - a Poisson number of point mutations per child (`mu_rgene`);
  - each picks a slot and a field: `k_idx ± 1` (most of the weight), or `C`, `A`, `B` redrawn;
  - an active-bit pass with duplication on activation, as for egenes;
  - mutations producing A = B are rejected.
- **Why discrete `k`:**
  - Activity statistics need genome *identity*, and continuous Gaussian mutation makes every birth a new genome, so cumulative activity saturates trivially.
  - A log grid with ±1 steps keeps neutral and near-neutral variants identical, keeps the genome hashable (lineage tables, activity, shadows), and makes "how far k moved" an integer.
- **Two compounds:** the useful genes are F → W with C = none or C = W. The format already holds stage 3.

### 8b. Gene cost

The current `tax_per_egene` is a mass tax per egene position, so it goes. Without some cost, dup-on-activate fills all 8 slots, and genome bloat would dominate activity. Recommendation: an energy cost `c_g` per active rgene per step (enzyme maintenance). It is the only energy tax, and it is physically motivated. Default small, scanned once.

### 8c. Transport (uptake and export)

- v1 has egene eating for F and `disp` for W.
- Recommendation for stage 1: passive two-way transport for both compounds, `d_c = u_c · (E_c − a_c)`, with fixed `u_F` and `u_W` as metaparameters. Two-way is what propagated in [CHEM-N1](../ISSUES.md).
- Transport genes `(compound, u_idx)` can be stage 2 or 3. They supersede wgenes (§8e).

### 8d. Does a chem agent keep the CA (LUT and egene eating)?

This is N's call. My recommendation is **no in stages 1-2, then optional CA gating in stage 4**:

- **Without the CA**, the model isolates metabolic evolution and RD pattern formation. If [CHEM-Q4](../ISSUES.md)'s speckle comes from CA-pattern intake (a hypothesis, not established), v2 will show coherent phase where v1 could not, and that test is worth having on its own.
- **The case for keeping it:** without the CA, v2 is an individual-based reaction-diffusion model, and EvoCA's distinctive feature, heritable local dynamics coupled to resources, is gone.
- **The middle path** (stage 4): the CA state `v` gates expression. An rgene runs only when `v = 1`, so LUT evolution sets the temporal program of metabolism (dormancy under toxicity, alternation). For that to be adaptive, the LUT inputs need some chemical sensing (for example, a bit for `w > tox/2`). That is a separate design step, best taken after stage 2 shows what the chemistry does alone.

### 8e. Conflict: [CHEM-D3](../ISSUES.md) / [CHEM-N3](../ISSUES.md) (wgenes)

- [CHEM-D3](../ISSUES.md) decided on wgenes: egene-like D4 patterns that scale excretion, in the *current* model.
- In v2 without a CA, a pattern gene has no input. The heritable-exchange role passes to transport genes (`compound, u_idx`).
- [CHEM-D3](../ISSUES.md)'s motivation was that eating and excreting would demand different CA patterns and so select for temporal alternation. That maps onto the stage-4 CA gating, not onto stage 1.
- Resolution options:
  - (a) pause [CHEM-N3](../ISSUES.md) until N decides on v2 (my recommendation: building wgenes into v1 now would be dev in the model being superseded);
  - (b) build wgenes in v1 anyway, as the CA-coupled baseline for comparing against v2 at stage 4.

---

## 9. Code organization

N's point 7: branches, or a separate model (`sim_build_chem`), possibly with separate C code?

### 9a. Options

| Option | Description | For | Against | Verdict |
|---|---|---|---|---|
| A. Long-lived branch | v2 lives on `chem2` | isolation | the Research tree diverges, shared-infrastructure fixes need two merges, notebooks cannot hold v1 and v2 side by side (N's own objection) | Rejected |
| B. Mode flag in `evoca.c` | `if (g_model == CHEM)` in the step | one file | `evoca.c` is 4865 lines; the hot loops already carry `n_chem` branches; `f`, eviction and death semantics differ, so nearly every phase forks; the bit-exact goldens would guard a file being rewritten throughout | Rejected |
| **C. Separate model, same repo and library** | `C/chemca.c` with its own state, step and genome; a shared core; Python `ChemCA` | the Research tree stays on `main`; v1 is untouched (its goldens stay meaningful); the shared infrastructure stays single | needs a shared sim interface and a `controls.py` refactor | **Recommended** |
| D. Separate repo | a new project | total isolation | duplicates display, controls, scans, galleries and recipes; the Research tree splits | Rejected |

### 9b. C layout (option C)

- **`C/core.c`, `C/core.h`.** Extract only what both models need and what is pure:
  - the PRNG and `poisson_sample`;
  - one `diffuse(field, tmp, d)`, the [CHEM-P1](../ISSUES.md) refactor, which unifies the three copies of the box blur;
  - Moore-neighbor helpers.

  These must be byte-for-byte equivalent in behavior; the existing goldens (`chem_parity`, food goldens) are the guard. The alternative of copying the helpers into `chemca.c` leaves two copies of the PRNG to drift, against the "state it once" rule.
- **`C/chemca.c`, `C/chemca.h`.** Prefix `chemca_` (avoiding the existing `evoca_set_chem_*` names, which belong to v1's food/waste). File statics hold its state: `N`, `alive`, `E[c][N²]`, the per-cell stores, `e[N²]`, `rgenes[N²·8]`, `ractive[N²]`, `r_field[N²]`, `G[c]`, and the metaparameters.
- **Step order (stage 1):**
  1. environment: W → F by `r(x)`;
  2. diffusion of each `E_c`;
  3. transport (two-way, every compound);
  4. reactions (rgenes; mass action, exponential integrator) into `e`;
  5. costs (`c_g` per active gene, `m_e`);
  6. toxicity and energy death (stores return to the site);
  7. reproduction (`e ≥ 1`; evict the lowest `e`; split; mutate);
  8. bookkeeping: the mass identity and the energy flux.
- **Store layout.** Dense `[N²·n_chem]` for `n_chem ≤ 8`: 4 MB at N=512 and 4 compounds, simpler and faster than slots. Pinned decision A (sparse slots) was made for universes of about 1000 compounds; it applies in stage 3 at the earliest, and the store is internal to `chemca.c`, so the layout can change without touching the API. This is a deliberate departure from decision A for the small case; I am flagging it rather than doing it silently.
- **Build.** One library: `gcc -O2 -Wall -fPIC -dynamiclib -o C/libevoca.dylib C/evoca.c C/core.c C/chemca.c`. The CLAUDE.md build line, the torque build and the Makefile all change at once (and the Makefile's zsh failure gets fixed while there).
- **Display.** `chemca_colorize(buf, mode)` fills the same ARGB buffer, so `sdl_worker.py` needs no change: it only sees pixels.

### 9c. Python

- **`python/sim_base.py`.** The interface both models implement: `init(N, **params)`, `reinit()`, `state(...)`, `step(n)`, `set_seed`, `get_alive`, `colorize(buf, mode)`, `colormode_count`, `colormode_names`, `params()`, `export_recipe`, `free`, and `SLIDERS`, a declarative slider spec of (name, min, max, step, log, setter).
- **`python/chemca_py.py`.** The `ChemCA` class with ctypes bindings and getters (`get_E(c)`, `get_store(c)`, `get_e`, `get_rgenes`, the mass and energy fluxes).
- **Recipes** carry `"model": "evoca" | "chemca"`. A missing key is read as `"evoca"` at load time (data corrections at load time), so every existing recipe and scan row loads unchanged.
- **One entry point.** `build_sim(params, ...)` and `make_sim(N, **params)` dispatch on `params.get('model', 'evoca')`. A `build_chem_sim` alias is fine if N likes the name, but scans, galleries and recipes should route by the recipe's `model` key so that there is one source of truth, not a second argument list.
- **Run harness.** `run_sim` stays EvoCA's. Its roughly 30 metrics are mostly LUT- or egene-specific (41 references to activity, shadow, dyn and eg). A `run_chemca` gets its own metric set: the §2c wave statistics, the mass and energy budget, the rgene diversity and the activity on rgene hashes. `build_sim` dispatches; the metric functions do not pretend to be shared.
- **`controls.py`** is the real cost.
  - It makes 56 direct `sim._lib.` calls and has 12 hand-built sliders. It must call sim methods and build its sliders from `sim.SLIDERS`.
  - This **moves a layer boundary** (C calls migrate from the UI into the model classes) and **makes the panel pluggable** (a slider registry). Both are moments where expensive defects cluster.
  - Guards: the existing headless harness (`tests/test_controls_pause_label.py`) extended to both models, plus a reachability test that opens the panel for each model and moves each declared slider.
  - At that point I suggest `/code-review ultra`, which only N can launch.

### 9d. Research tree

- v2 work stays in `Research/Chem/`. Scans record `model` in `scan_config.json`.
- The v1 results ([CHEM-N1](../ISSUES.md), [CHEM-N5](../ISSUES.md)) remain reproducible because v1 is untouched.
- Notebook sections label their model.
- Development happens on a short-lived feature branch per stage with the usual merge gate (test notebook, review), not a long-lived model branch.

---

## 10. Tests and guards

| Guard | What | Fails without |
|---|---|---|
| Core extraction is behavior-neutral | all existing goldens bit-exact after `core.c` (chem parity, food goldens, waste identities) | any drift in the PRNG or diffusion |
| Mass identity, every step | `ΣE_F + ΣE_W + Σ stores` constant to float tolerance, through reactions, transport, deaths, eviction and reproduction | any leak (each death path) |
| Energy flux identity | `Δ(Σe) = produced − reproduction overhead − costs − lost at death` | a missing booking |
| Constructed one-step fixtures | one cell: uncatalyzed and autocatalytic conversion against `a(1 − e^{−kc})`; transport both directions; toxic death returns the stores; an uphill reaction stops when `e` runs out | each rate law |
| Cycle yield zero | a constructed A→B→A gene pair yields `e` change 0 | per-reaction yields |
| Mutation reachability | through `init(..., mu_rgene=...)` and `run_chemca`: `k_idx` moves, `C` toggles, A = B rejected | params not forwarded |
| Recipe round trip | export → load → identical state hash, with `model: chemca`; an old recipe without `model` loads as EvoCA | the dispatch default |
| Controls, both models | the headless panel opens; every `SLIDERS` entry reaches its setter | the hardcoded `_lib` calls |

Each regression test is shown failing with its fix reverted (standing rule).

---

## 11. Staged plan

| Stage | Content | Deliverable | Gate |
|---|---|---|---|
| **0a** | Statistics: §2c in `evoca_waves.py` (Hilbert phase, tracked defects, ξ, S(k), the space-time spectrum); constructed-field tests (meandering spiral pair, plane wave, target source, noise) | module, tests | the tests separate the four constructed cases |
| **0b** | Eye calibration (§2d) plus the one-way tox 0.4 / sink 0.019 coverage check | labeled set, a table of which statistics separate labels | N's labels |
| **0c** | `C/core.c` extraction ([CHEM-P1](../ISSUES.md)) | refactor | all goldens bit-exact |
| **1** | `chemca` skeleton: two compounds, mass-conserving, `r(x)` regeneration, two-way transport, rgenes F → W with C ∈ {none, W}, energy currency, `c_g`, `m_e`, toxicity; `ChemCA`; recipes `model`; colormodes; `controls.py` refactor | branch, test notebook, `/code-review ultra` offered | §10 table green; mass identity exact |
| **2** | Experiments: (`m̄`, `r`) regime map; does evolution choose autocatalysis; does `k` track pattern phase; the generational Pareto driver (§3c) on (pattern, activity) | scans, Chem notebook section | the §2d-calibrated statistics |
| **3** | `n_chem > 2`: free energies `G_c`, the full (A, B, C, k) genes, transport genes, the light-harvesting gene (the producer niche); sparse store if the universe grows | dev plan | cycle-yield and mass tests |
| **4** | Optional CA coupling: `v` gates expression, LUT chemical sensing | design doc first | N's call |

Stages 0a-0c are independent and can run in parallel. Stage 1 needs 0c. Stage 2 needs 0b.

---

## 12. Adversarial critique of the first draft

I wrote a first draft and then critiqued it from the position of a rival expert (an artificial-chemistry and RD modeler) holding an opposing draft. These are the substantive flaws found, each with the revision it produced.

1. **The currency flaw.**
   - *Draft:* "replace the tax with `v = k · f` food→waste; the rest unchanged".
   - *Objection:* under decision B, conversion lowers the reproduction currency, so `k → 0` and the rate gene is dead on arrival.
   - *Revision:* §4. A separate energy currency `e`, eviction by `e`, and the conflict with pinned decision B stated explicitly rather than patched.
2. **Accepting the premise of point 4 uncritically.**
   - *Draft:* "toxicity drives rates toward low waste; add a need for waste to counter it".
   - *Objection:* stoichiometry fixes waste per unit of energy, so the gene sets tempo, not waste per offspring. The tension (tempo against toxicity) already exists. What is missing is density dependence, without which the interior optimum gives stasis again.
   - *Revision:* §6a, and the recommendation reframed as density dependence. Autocatalysis is offered as a *gene choice* (C = none or W) with the plain variant as a control, rather than built into the rate law.
3. **Continuous rates.**
   - *Draft:* `k` a float with Gaussian mutation.
   - *Objection:* activity statistics key on genome identity, and continuous mutation makes every birth novel, so activity saturates trivially, the opposite of what moving activity into chem needs.
   - *Revision:* discrete log-spaced `k_idx` (§8a).
4. **Overclaiming the mass-conservation fix.**
   - *Draft:* "a closed mass budget removes the narrow balance band of [CHEM-N5](../ISSUES.md)".
   - *Objection:* unjustified generalization. The band comes from matching a fixed sink to density-proportional production, and in a closed system the matching moves to `(m̄, r)`.
   - *Revision:* §5 states it as a hypothesis with a test. The claim kept is the structural one: W is bounded by `M`, which does follow. An open mode is kept for comparison.
5. **Pareto generations now.**
   - *Draft:* "run a generational Pareto search on the chem parameters next".
   - *Objection:* the objectives are miscalibrated (§2b), and single-seed arg-max selects the lucky (pdt-server's own critique).
   - *Revision:* calibrate first, use two seeds with a spread penalty, re-evaluate survivors with fresh seeds, and build the driver model-agnostic for v2 (§3b, §3c).
6. **Copying helpers into `chemca.c`.**
   - *Draft:* a separate library with copied PRNG and diffusion.
   - *Objection:* two copies of one fact drift.
   - *Revision:* extract `core.c` under the goldens (§9b), which also retires [CHEM-P1](../ISSUES.md).
7. **Silent departure from decision A.**
   - *Draft:* dense stores, no comment.
   - *Revision:* flagged in §9b, with the reason (decision A targets large universes) and the condition for switching back.

A further objection I considered and did not adopt: "drop the CA permanently; it only adds noise". The speckle-from-CA claim is a hypothesis, and EvoCA's identity is the CA-resource coupling. So §8d keeps CA gating as stage 4 and makes it N's decision.

---

## 13. Decisions for N

1. **Currency.** Adopt an energy currency `e` produced by reactions, replacing pinned decision B for v2 (§4)? Keep `m_e` (energy maintenance) as an optional metaparameter, default 0?
2. **Mass.** Closed mass with W → F regeneration by a light field `r(x)` as the default, with an open mode kept for comparison (§5)?
3. **The need for waste.** Autocatalysis offered as a gene choice (C = none or W), with plain metabolism as the control (§6b)?
4. **Rate law.** Mass action with the exponential integrator; discrete log-spaced rates; free energies `G_c` per compound (§7, §8a)?
5. **The CA.** Leave the CA and egene eating out of stages 1-2, with CA gating as a stage-4 option (§8d)?
6. **Wgenes.** Pause [CHEM-N3](../ISSUES.md) (wgenes in v1) pending v2, or build them as the CA-coupled baseline (§8e)?
7. **Code.** A separate model `chemca` in the same repo and library, with a shared `core.c`, a sim interface, recipe `model` keys and the `controls.py` refactor (§9)? The name `chemca` is a placeholder.
8. **Statistics.** Do stages 0a-0b first (new statistics plus your labels on about 30 films), and the one-way tox 0.4 coverage check now (§2c, §2d)?
9. **Pareto generations.** Build the model-agnostic driver (§3c) after calibration, and answer [ACT-D2](../../Activity-metrics/ISSUES.md) at the same time (the same `select` call)?
