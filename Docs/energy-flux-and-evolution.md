# Free-energy flux and evolution

**Anchor document.** N, 2026-10-24:

> We seek to find relationships between (free) energy flux through an evolving system and its
> evolutionary properties. Ideally we want to articulate laws relating energy flux to evolution
> and evolvability.

Every workstream, campaign and model extension, chemca included, should be able to say in one
line what it contributes to this goal. When it cannot, it is a rabbit hole. That can be fine,
but only as a deliberate choice. This document makes the goal operational. It builds on the
thermodynamic framing already in `Docs/checkpoint 3 July 2026.md` §L (Prigogine, Schneider & Kay,
Lotka, Chaisson) and `Docs/EvoCA_research_directions.md` §6, and on measured data, rather than
restating them.

## 0. Summary

- **A law needs three things we do not yet fully have:**
  1. an energy-flux variable that can be varied independently, and in the end evolves;
  2. evolutionary measures that are valid in the regime being measured;
  3. a collapse test: different parameter routes to the same flux give the same evolutionary outcome.
- **In the open models (vanilla, food/waste), total throughput is imposed.** At steady state it equals the source (`food_inc` x sites); it is identical, about 4158 per step, across all five [CHEM-N7](../Research/Chem/ISSUES.md) conditions. Evolution can change only:
  - **per-capita** throughput Φ;
  - the **partition** of the flux into metabolism (tax), death and unharvested overflow;
  - its spatial and temporal **structure** (waves).

  The law-bearing variables in the open models are partition and structure, not total flux.
- **Prior data already show partition and structure tracking evolution.**
  - Wave regime ([CHEM-N7](../Research/Chem/ISSUES.md)): 12% of the flux goes unharvested and deaths rise 35-fold. Per-capita Φ is up 18%, turnover is 5 times faster, genotype persistence above drift collapses, and a trait (the egene count) still rises.
  - Resource-driving pilot: A_new is an inverted U against `food_inc/tax`, while Φ is nearly degenerate (Φ = max(tax, food_inc)). Per-capita flux is pinned by the imposed tax.
- **chemca is the model in which flux becomes endogenous**, which is the reason to build it.
  - Mass is closed, and free energy enters only by light acting on waste (P -> F).
  - Reaction rates are genes, and the free-energy yield is G_A - G_B per conversion.
  - So throughput depends on what the population does, and Lotka's maximum-power hypothesis (selection increases the throughput a population captures) becomes testable. It is untestable in vanilla, where throughput is fixed.
- **The current evolutionary measures are the weakest link.**
  - The activity statistics reward persistence.
  - The activity flux clips at its buffer and is about D/10 when genomes are singletons.
  - A_new is not well posed when egenes evolve.

  Laws need a valid measure of evolutionary change and of evolvability. That is a first-class task, not housekeeping.

## 1. Outline

| § | Section | Purpose |
|---|---|---|
| 2 | Energy flux, operationally | what is imposed and what is emergent in each model; the accounting identities |
| 3 | Evolution and evolvability, operationally | which measures are valid where; known defects |
| 4 | Candidate laws | testable forms, each with its prediction and falsifier |
| 5 | Prior data against the candidates | what we already know |
| 6 | Why chemca serves the goal | endogenous flux; free-energy accounting |
| 7 | Standards for a law-grade result | collapse tests, decorrelation, scaling, seeds |
| 8 | Map of the work onto the goal | every workstream's one line; the rabbit-hole test |
| 9 | Adversarial critique | flaws found in the draft and the revisions they produced |
| 10 | Decisions for N | |

## 2. Energy flux, operationally

| Model | Free-energy input | Sinks (dissipation channels) | Imposed | Emergent (evolution can change it) |
|---|---|---|---|---|
| vanilla | `food_inc` per site (capped at F <= 1) | tax (metabolism), death (stores lost), overflow (input past the cap), eviction | total input; tax per cell | density; per-capita Φ (pinned at about tax when food-limited, [RES-Q2](../Research/Resource-driving/ISSUES.md)); harvested fraction; spatial structure |
| food/waste (`n_chem=2`) | same | same; the tax is also the source of waste mass | total input; tax; waste sink, toxicity, exchange | partition (tax, death, overflow); waves; refractory structure |
| chemca closed | light on waste: r(x) x P (or W) | agent metabolism (reaction yield into `e`, then costs and reproduction), death (e lost), env relaxation | light rate r; total mass m̄; G_c | **total throughput** (depends on how much waste the population produces and where); per-capita yield; reaction rates themselves (genes) |

Accounting identities (built and tested):
- vanilla/food-waste: the L1 dissipation accountant (`get_dissipation`; source = tax + death + overflow + Δstorage);
- chemca: the mass identity, plus an energy budget (light input, agent yield, costs, death loss; engine in progress).

Identities are what turn "energy" from a metaphor into a variable.

**Free energy versus energy.** In vanilla and food/waste, "food" is free energy only in the
source-sink sense (`checkpoint` §L: held in disequilibrium by the tax sink). There is no
entropy term. In chemca, G_c makes reaction free energy explicit: a closed cycle yields zero, and
uphill reactions cost. Only chemca can support a statement about *free* energy in the strict
sense.

## 3. Evolution and evolvability, operationally

| Measure | What it measures | Valid where | Known defect |
|---|---|---|---|
| `dyn_excess`, `eg_excess` (activity above a drift shadow) | persistence of components above neutral | moderate turnover, components with populations | rewards persistence: stasis scores high ([CHEM-N7](../Research/Chem/ISSUES.md)); dyn artifact when `mu_lut = 0` freezes a rich substrate |
| A_new (Bedau-Packard) | count of components above the neutral p90 | LUT-only evolution | not well posed when egenes evolve (shadow is LUT-only) |
| activity flux | population of components crossing a low activity band | genomes with populations (D << pop) | clips at the 65536-sample buffer; about D/10 when genomes are singletons ([CHEM-N7](../Research/Chem/ISSUES.md) notes) |
| trait trajectories (egene count, specificity, ring level, k in chemca) | directional phenotypic change | everywhere | needs a neutral expectation (the drift of a trait under the same mutation) |
| invasion / competition assays (`patch_transfer.py`) | fitness of evolved against ancestral types in a given medium | everywhere | per-assay cost; the choice of medium matters |
| lineage tracing (lineage-return reducer) | which lineages persist and where (wave sources, Savill et al. 1997) | everywhere | lineage identity under mutation |
| **evolvability** | rate at which heritable variation produces fitness-increasing change | undefined in the code so far | needs a definition: candidates are the variance of offspring fitness per mutation (a mutational-neighborhood scan) and the response to a change of environment (re-adaptation time) |

**The law-relevant evolutionary variables** should be rates per unit time *and* per unit energy:
adaptations per generation and per unit of dissipated energy. The second normalization is what
makes a statement about energy and evolution, rather than about time.

## 4. Candidate laws (hypotheses, each with a falsifier)

- **L1, the window (inverted U; Schneider & Kay, Chaisson, `research_directions` §6).** Selectable complexity and adaptive change peak at an intermediate source/sink gradient (`food_inc/tax`). There is extinction below the window and drift (gradient collapse, fitness variance near 0) above it.
  - Falsifier: curves at different `tax` do not superimpose against the ratio, or complexity keeps rising on the rich flank ([RES-N1](../Research/Resource-driving/ISSUES.md)).
- **L2, maximum power (Lotka 1922).** Where throughput is endogenous, selection increases the free-energy flux the population captures, up to a limit set by a trade-off (in chemca: tempo against toxicity).
  - Prediction in chemca: mean k and captured light flux rise over evolutionary time until waste toxicity or ignition limits bind.
  - Falsifier: throughput falls or wanders while fitness rises, i.e. selection favors efficiency (yield per unit flux) over power.
  - The alternative, minimum dissipation or maximum efficiency, is the competing law. Which one wins, and when, is itself the law-level result.
- **L3, partition.** The rate of evolutionary change scales with the fraction of the flux routed through **turnover** (death plus birth), not with total flux.
  - [CHEM-N7](../Research/Chem/ISSUES.md) is suggestive: waves raise the death and overflow shares, and turnover rises 5 times.
  - Falsifier: equal turnover fractions reached by different routes (tax, toxicity, imposed mortality) give different evolutionary rates.
- **L4, structured flux.** Spatial and temporal structure in the flux (waves, lumpiness, refractory zones) changes the *level* at which selection acts (pattern-level selection, Boerlijst & Hogeweg; Savill et al. 1997) and the evolvability, at fixed mean flux.
  - Tests: [RES-N5](../Research/Resource-driving/ISSUES.md) (lumpy food), the [CHEM-N7](../Research/Chem/ISSUES.md) wave and no-wave contrast with a fitness measure, S1/S2 of the Hogeweg survey.
  - Falsifier: no difference at matched mean flux and turnover.
- **L5, evolvability and flux variance.** Evolvability (mutational neighborhoods with fitness variance) is selected up when the environment, here the flux, fluctuates on the time scale of a few generations. This is the energetic version of the "evolution of evolvability under fluctuating environments" results (Crombach & Hogeweg 2008, in `Docs/references/hogeweg/`).
  - Falsifier: no rise of evolvability measures under periodic flux (`food_per` > 1) against a constant flux.

Dimensionless groups to look for: source/sink ratio (`food_inc/tax`), turnover per unit flux
(births / dissipated energy), structure period / generation time (the wave period of 115 steps
against generation times of 42-217 steps in [CHEM-N7](../Research/Chem/ISSUES.md)), and in chemca a = k1·m̄²/k2 and ρ = r/k2
from the regime study.

## 5. Prior data against the candidates

| Source | Finding | Bears on |
|---|---|---|
| Resource-driving pilot ([RES-Q1](../Research/Resource-driving/ISSUES.md)/Q2) | A_new is an inverted U in `food_inc/tax` (peak at a ratio of about 0.6; collapse at 2.9); Φ is nearly degenerate (= max(tax, food_inc)); alive density = ratio below 1 | L1 supported (pilot only); in vanilla, per-capita Φ is not an independent axis |
| R1 rescan ([RES-Q4](../Research/Resource-driving/ISSUES.md)) | extinction AUC is slightly better for the ratio than for `food_inc` alone, but the `tax` range is too narrow | L1 collapse test pending ([RES-N1](../Research/Resource-driving/ISSUES.md)) |
| [CHEM-N7](../Research/Chem/ISSUES.md) (2026-10-24) | equal total flux in every condition; waves shift the partition (overflow 12%, death x35, per-capita Φ +18%), turnover x5, persistence above drift collapses, egene count rises | L3, L4 (suggestive); the measure validity problem in §3 |
| [CHEM-N6](../Research/Chem/ISSUES.md) | excitation waves; period about 115 steps | gives the structure variable for L4 |
| GS regime study | closed-loop chemistry: excitable pulses and spirals at a about 30-100, ρ about 0.03-0.07; an oscillatory ridge only | sets chemca's starting point for L2 |
| vanilla stasis (memory, 2026-10-23) | evolution settles; activity is persistence | a law has to cover the stasis endpoint (L1 high flank?) |

## 6. Why chemca serves the goal (and what to protect in its design)

In vanilla, total throughput is imposed and per-capita throughput is pinned by the tax, so
evolution cannot change how much energy flows. Lotka's question cannot be asked. chemca changes
three things:

1. **Throughput is endogenous.** Light acts on waste, so the energy entering the system depends on the waste the population makes and where. Total flux becomes a dependent, evolvable variable.
2. **Free energy is exact.** G_c, closed cycles yield zero, and uphill reactions cost. The energy budget is an identity, not a proxy.
3. **The rate of exploitation is a gene** (k), which is the very quantity Lotka's principle is about.

Design choices to protect for the goal's sake:

- keep the energy budget identity exact and tested;
- record light input, captured yield and dissipation per channel as standard run outputs;
- keep an open mode, so the same chemistry can run with imposed or endogenous flux. That comparison is a direct test of whether endogeneity matters.

The mass-gate option (Kaneko style) removes energy from the reproduction gate. It is useful as
a *contrast*: whether evolution tracks energy flux when energy does not gate reproduction.

## 7. Standards for a law-grade result

1. **An independent axis.** Vary flux by at least two routes (for example `food_inc` and `tax`; in chemca, r and m̄) and test **collapse**: outcomes superimpose against the flux variable, not against the route.
2. **Decorrelate the confounds.** Density, turnover and food per cell co-vary with flux. Report them, and use matched-turnover or matched-density contrasts (an imposed-mortality control).
3. **A valid evolutionary measure in the regime tested** (§3). For now, trait trajectories and invasion fitness are the safest; activity statistics only where their conditions hold.
4. **Scale and replication.** Several N (finite size), at least 3 seeds, runs long enough relative to the generation time.
5. **Per unit energy.** Report evolutionary rates per generation *and* per unit dissipated energy.
6. **Prefer a dimensionless form.** State laws in dimensionless groups (§4) so they can transfer between models.

## 8. Map of the work onto the goal

| Workstream / item | One line on the goal |
|---|---|
| Resource-driving ([RES-N1](../Research/Resource-driving/ISSUES.md) inverted U, [RES-N5](../Research/Resource-driving/ISSUES.md) lumpy food) | L1 directly (imposed flux, collapse test); L5 via lumpiness. **Central, overdue.** |
| Chem Track A (food/waste waves, N6/N7, wgenes, S1/S2) | L3 and L4: how flux structure and partition change evolution at fixed total flux |
| Chem Track B (chemca: fw2, gs3, mass gate) | L2: endogenous flux, maximum power or efficiency; exact free energy |
| Activity-metrics | the measures that §3 says are the weakest link; fix the flux clipping; define evolvability |
| Lineage-co-evolution (net return) | per-lineage energy capture: the individual-level version of L2 |
| Criticality-bifurcation, Ring-complexity | supporting: transitions in the evolutionary dynamics as functions of flux variables |
| Reaction-diffusion | structure of the flux field (L4) |

**The rabbit-hole test**, to apply before starting anything new: "What will this tell us about how
free-energy flux relates to evolution or evolvability, and which candidate law does it bear on?"
If the answer is "none, but it unblocks X", name X.

## 9. Adversarial critique of the draft

An imagined opponent: a statistical physicist who works on nonequilibrium thermodynamics and is
skeptical of "energy laws of evolution". The substantive flaws found, and the revisions they
produced:

1. **"Vary the energy flux" is ill-posed in the open models.**
   - Draft: treat `food_inc` as the flux axis.
   - Objection: at steady state, total throughput equals the source regardless of evolution ([CHEM-N7](../Research/Chem/ISSUES.md): about 4158 in every condition). The axis is imposed, and evolution cannot respond along it.
   - Revision: §0 and §2 separate total, per-capita, partition and structure. The open-model laws are restated in partition and structure terms (L3, L4). L2 is assigned to chemca, where flux is endogenous.
2. **The evolutionary measures were taken as given.**
   - Draft: plot activity against flux.
   - Objection: the activity statistics reward persistence, the flux clips and degenerates, and A_new is ill-posed with egenes. A law built on a defective measure is an artifact law, which is our documented failure mode: plausible output, no loud failure.
   - Revision: §3 is a validity table. §7 requires a valid measure per regime. Activity-metrics is promoted to "weakest link".
3. **Confounding of flux with turnover and density.**
   - Draft: correlations of flux and evolution across conditions.
   - Objection: every route to a different flux also changes density and turnover, so a correlation does not separate "energy" from "demography".
   - Revision: §7 items 1-2 (collapse across routes, matched-turnover controls). L3 is stated as a turnover law, so it competes with an energy law instead of hiding inside one.
4. **"Free energy" was used loosely.**
   - Objection: vanilla has no entropy, so "free energy" is only a source-sink metaphor there.
   - Revision: §2 says plainly that only chemca supports strict free-energy statements, and that the open models support energy-throughput statements.
5. **Maximum power was presented as the law.**
   - Objection: competing principles exist (minimum dissipation, maximum efficiency), and which applies depends on constraints.
   - Revision: L2 is stated with its competitor, and the deciding result is the outcome, not an assumption.

## 10. Decisions for N

1. **Adopt this as the anchor document**, linked from `Research/README.md` and checked at each plan ("rabbit-hole test", §8)?
2. **Promote [RES-N1](../Research/Resource-driving/ISSUES.md)** (the inverted-U campaign with the `tax`-collapse test) to the top of the queue? It is the most direct existing test of L1, and it has waited since October.
3. **Define evolvability operationally** (§3, last row) as an Activity-metrics task, and fix the flux clipping as part of it?
4. **Make energy-budget outputs standard** in chemca's `run_chemca`: light input, captured yield, dissipation by channel.
