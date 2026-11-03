# Pareto search: generations of cone/box selection instead of grids

_2026-11-03. Moved from `Research/Chem/plans/2026-10-23_chem-model-v2.md` §3 (written 2026-10-24),
because the search is model-agnostic. Status: design only, not built ([ACT-N8](../ISSUES.md))._

N's idea: replace one-shot parameter grids with a few generations of Pareto search on two
objectives. The selection half exists (`~/Projects/conepareto`,
[`2026-10-15_cone-pareto-package.md`](2026-10-15_cone-pareto-package.md)); the sampling half
does not.

## 0. Objective pairs (N, 2026-11-03)

The search needs two named objectives. N's candidates:

| pair | objective 1 | objective 2 | notes |
|---|---|---|---|
| **A (the original)** | evolutionary activity | spatial structure (e.g. RD patterns, fronts) | the pair behind the 04-27 RD scans and [RXD-N1](../../Reaction-diffusion/ISSUES.md); structure needs a calibrated statistic (chem v2 plan §2; the s6 front-spacing and layered-front measures are candidates) |
| **B** | LUT activity | egene activity | two activity channels; the (`dyn`, `eg`) drift-shadow excesses already proposed for selection in [ACT-D2](../ISSUES.md) (box region, recommended yes) |

For either, the activity measure must reward adaptation, not persistence: the N7 evolution
scan found that the excess statistics mostly reward long-lived genotypes, which is why the
activity flux is being re-tested (`Research/Chem/Scans/2026-10-24_N7_evolution_waves`, rerun
2026-11-03). Which activity readout fills these slots is [ACT-D1](../ISSUES.md).

## 1. History

| Date | Event |
|---|---|
| 2026-05-21 | `pdt-server/pdt-param-optimization.md` compares PDT's generation loop (explore/exploit sampling, refit per generation) with EvoCA's one-shot grids: "the thing PDT has that EvoCA lacks: sequential, surrogate-driven active learning." |
| Oct 2026 | `prompt1-Oct2026.md`: N asks to abstract the restricted (cone) Pareto approach from pdt-server for any dual-parameter optimization. |
| Oct 2026 | Review `Docs/Archive/Oct2026/00-review-of-prompt1.md` §6-§8. EvoCA's "optimization" is scalarization over a static grid, with weaknesses: compensation, concave regions unreachable, pool-relative normalization, silent weighting. It proposes a package. |
| Oct 2026 | `~/Projects/conepareto` built, with R parity (`Research/Activity-metrics/plans/2026-10-15_cone-pareto-package.md`): `select`, layers, regions, hypervolume, and later wedge-angle derivation and `on_shortfall`. **"Not ported: [...] the lambda-sweep sampler. They sit around `select` rather than inside it; say if you want them."** No answer was recorded. |
| Oct 2026 | Used post-hoc: `rd_pareto_picks` / `rd_pareto_layers` galleries (Reaction-diffusion), and Activity-metrics notebook 4a. |
| open | [ACT-D2](../ISSUES.md) ("adopt Pareto selection for (dyn, eg)?", recommendation yes, box region) and [RXD-N1](../../Reaction-diffusion/ISSUES.md) (re-select the 04-27 winners with box/wedge) are still open on the board. |

So the generational search was never built. The selection half exists; the sampling half does not. I found no board item for it in EvoCA. It fell between the pdt comparison and the port's "not ported" list.

## 2. Why it was not run on the chem scans (2026-10-24)

- **The objective is not calibrated** (chem v2 plan §2). A generational search climbs whatever statistic it is given. With the current spiral statistic it would climb toward stationary target sources, the thing that statistic actually rewards.
- **Single-seed noise.** pdt-server's own critique (§2b of the comparison) applies: an arg-max over noisy evaluations selects the lucky ones. The [CHEM-N1](../../Chem/ISSUES.md) seeds disagree enough (defect excess 0.15-0.21) to matter.
- **The parameter space is about to change.** If v2 is adopted, the v1 parameters (`waste_sink`, `food_inc`) disappear (chem v2 plan §5). A driver written generically survives that; a v1-specific campaign does not.

## 3. Driver sketch (model-agnostic)

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
- Budget: 5 generations × 50 points × 2 seeds = 500 runs, about the cost of four [CHEM-N1](../../Chem/ISSUES.md) grids, over a 5-8-parameter space a grid cannot cover.
- Two objectives must be named first. For chem, candidates are (pattern coherence ξ or tracked-core lifetime, persistence or viability), and in v2 (pattern, evolutionary activity on reaction genes).

The surrogate model PDT uses is left out on purpose: EvoCA evaluates the true objective cheaply enough, and §1 of the pdt comparison (`~/Projects/pdt-server/pdt-param-optimization.md`) argues that the surrogate is where PDT's weakness lies.
