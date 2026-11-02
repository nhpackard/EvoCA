# The energy gradient in EvoCA

_Response to N's [`energy-analog.md`](energy-analog.md), 2026-10-25. It bears on candidate laws
L2 (maximum power against efficiency) and L4 (structured flux) in
[`energy-flux-and-evolution.md`](energy-flux-and-evolution.md)._

## 1. Summary

- **N's worry is right, but "missing" is slightly too strong.** Every energy flow in EvoCA has a
  source and a sink. In vanilla EvoCA, though, the sink is **perfect and stateless**: the tax
  disappears into nothing, and a unit of food is worth the same everywhere, always. The gradient
  exists, but it is infinite and constant, so nothing that evolves can act on its low side.
- **The food/waste model already gives the sink a state.** Waste piles up where the sink is slow,
  excretion stops as local waste approaches 1, and private waste kills.
  - The waste field is what made the medium refractory and produced the waves ([CHEM-N6](../Research/Chem/ISSUES.md), [CHEM-N7](../Research/Chem/ISSUES.md)).
    The most interesting dynamics found so far came from making the sink finite.
  - But the sink acts only as a **hazard**, not as a **yield**. A unit of food eaten next to a waste
    dump is worth exactly as much as one eaten in clean surroundings.
- **`chemca` does not have it either.** The yield of A → B is the fixed difference `G_A − G_B`
  (`Research/Chem/plans/2026-10-23_chem-model-v2.md` §7). There is no dependence on the
  concentrations of A and B, which is the part of thermodynamics that makes a gradient a gradient.
- **The main consequence: there is no rate-yield trade-off anywhere in EvoCA.** Taking food
  faster never lowers the energy captured per unit. Without that trade-off, L2 has nothing to
  bite on: maximum power and maximum efficiency point the same way, so "evolution maximizes
  power" and "evolution maximizes efficiency" cannot be told apart.
  - Our own data show the symptom. In `chemca`, `gs3` evolved the F → P shortcut, which captures
    twice the energy per conversion (1 against 0.5), and nothing opposed it.
- **Recommendation (§6):** a single heritable **coupling** trait `ε`, the energy captured per
  unit of food, with an eating rate that falls as `ε` approaches the local driving force
  `Δμ(x) = Δμ0 + ln(F/W)`.
  - This puts the gradient (both reservoirs, locally) and the rate-yield trade-off into the
    food/waste model with one formula.
  - It stays in Track A, is off by default (bit-identical), and makes L2 testable without `chemca`.

## 2. Where the gradient is now

| Model | Source (high side) | Sink (low side) | Does the sink's state matter? | Energy captured per unit food |
|---|---|---|---|---|
| Vanilla (`n_chem=1`) | `food_inc` regenerates `F` (capped at 1) | the tax: food simply vanishes | no: the sink is infinite and stateless | fixed, 1 |
| Food/waste (`n_chem=2`) | as vanilla | the tax turns private food into private waste; excretion `∝ (1 − W)`; environmental `waste_sink` | **yes, as a hazard**: excretion stalls, toxic death | fixed, 1 |
| `chemca` (branch) | light regenerates F from W | waste W, closed mass | yes: light needs W, mass is closed | fixed, `G_A − G_B` per conversion |

N's two examples are both cases where the **low side limits the power**:
- **The runner** is limited by rejecting heat. A cooler sink means more sustainable power.
- **The vent** is a gradient between two reservoirs. Life sits *between* them, not inside the
  hot one.

In EvoCA's food/waste model the low side limits *survival* (toxicity) but never *power*.

## 3. What a gradient does that a stock of energy does not

1. **It sets how much work a conversion can yield, from both ends.**
   - For a reaction F → W the free energy released is `Δμ = Δμ0 + kT ln([F]/[W])`.
   - As waste builds up, each conversion is worth less, and at equilibrium it is worth nothing.
   - A heat engine is the same: Carnot efficiency `1 − T_c/T_h` depends on the sink temperature.
2. **It makes speed cost efficiency.**
   - A conversion that captures `ε` units of the available `Δμ` runs at a net rate that falls as
     `ε → Δμ`. In flux-force form, `J ∝ 1 − exp(−(Δμ − ε)/kT)`.
   - Capturing everything means running at equilibrium, infinitely slowly. Running fast means
     capturing little.
   - So power `ε · J(ε)` peaks at an intermediate `ε`: **maximum power is not maximum efficiency**.
     Classic forms of this:
     - Odum and Pinkerton's "time's speed regulator" (maximum power at about 50% efficiency);
     - the Curzon–Ahlborn efficiency at maximum power, `1 − sqrt(T_c/T_h)`, for heat engines;
     - the rate-yield trade-off of ATP-producing pathways (Pfeiffer, Schuster and Bonhoeffer 2001).
3. **It makes the sink a resource, and space matters.**
   - Waste dumped locally lowers the neighbors' yield. That is a public bad, so fast, wasteful
     strategies are a tragedy of the commons.
   - Pfeiffer et al. showed that **spatial structure** favors the efficient, high-yield strategy:
     the wasteful cell poisons its own kin.
   - EvoCA is spatial, and its diffusion constants (`gdiff`, `wdiff`) set the mixing, so this
     prediction is directly testable here.
4. **Life organizes along gradients and dissipates them faster** (Schneider and Kay 1994,
   "Life as a manifestation of the second law"). This is the vent picture: a separated source and
   sink, with structure at the interface.

## 4. What this means for the program

- **L2 is untestable in the open models as they stand,** not only because total flux is imposed
  (the anchor doc's point) but because no trait trades rate against yield.
  - The anchor doc says L2 needs `chemca`. **That needs a correction:** `chemca` as designed also
    lacks the trade-off, since its yields are concentration-independent. With the coupling trait
    of §6, L2 becomes testable in Track A.
- **The prior data fit this reading:**
  - **`gs3` evolved the higher-yield F → P shortcut** (twice the energy per conversion). With a
    fixed yield per conversion, more yield is free, so selection takes it.
  - **`fw2` collapsed when F → W genes duplicated** until the waste was toxic (fixed by a gene
    cost). That is a rate escalation with no rate penalty except death.
    - A graded penalty (yield falling as local waste rises) would replace the cliff with a slope.
    - It would make that run a measurement rather than a failure.
  - **The waves ([CHEM-N6](../Research/Chem/ISSUES.md)/N7) came from the sink side.** The refractory zone is local waste. Making
    yield depend on `W` couples the waves to energy capture directly, not only through deaths.
- **Evolvability.** N asks whether a gradient would make the model "more evolvable". A trade-off
  creates a real fitness landscape where there was a ramp:
  - an interior optimum `ε*`;
  - an optimum that depends on the environment (local `F`, `W`, mixing);
  - frequency dependence (the commons).

  Those are the conditions for diversification and for ongoing adaptation when the environment
  moves. That makes greater evolvability plausible, but it is an **unproven assumption**, and it
  should be tested rather than designed in (§7).

## 5. Options

| | Change | Track | Cost | Bears on | Notes |
|---|---|---|---|---|---|
| G1 | Yield per unit food falls with local waste: `y = clip(Δμ(x)/Δμ0, 0, 1)` | A (food/waste) | small (C, one factor in eating) | L2 weakly, L3 | Gradient without a trade-off trait: yield is environmental, not heritable |
| **G2** | **Heritable coupling `ε`**: yield `ε`, rate factor `1 − exp(−(Δμ(x) − ε))` | A (food/waste) | moderate (one per-cell float trait, mutation, recipe and shadow plumbing) | **L2 directly**, L4 | Recommended; contains G1 as its environmental part |
| G3 | Separated source and sink (vent geometry): `F` regenerates in one region, `W` is removed in another (masks exist: `env_mask`; source/sink masks for waste are listed as a next step in `plan-chem.md`) | A | small (masks) | L1, L4 | Complements G2; Schneider–Kay test: does evolution increase gradient dissipation? |
| G4 | Concentration-dependent free energies in `chemca`: yield `G_A − G_B + ln([A]/[B])` | B | one-line change in the yield, but Track B | L2 | The `chemca` version of G1/G2 once Track B resumes |
| — | Tax depending on the local environment (the runner's ambient temperature) | A | small | L3 | Already partly present: toxicity plays this role |

## 6. Recommended design (G2, for discussion, not yet a plan)

All quantities are per site, in units where the vanilla model has yield 1.

- **Driving force:** `Δμ(x) = Δμ0 + ln((F(x) + δ) / (W(x) + δ))`.
  - `Δμ0` is a metaparameter, the standard free energy of eating.
  - `δ` regularizes the logarithm.
- **Trait:** each cell carries `ε ∈ [0, ε_max]`, the energy captured per unit food. It is
  inherited and mutates like the egene rates.
- **Eating:**
  - mouthful `M' = M · max(0, 1 − exp(−(Δμ(x) − ε)))`, where `M` is today's mouthful (egene match,
    caps);
  - private energy gain `ε · M' / ε_ref`;
  - waste produced per unit of energy, `1/ε` relative to today, so efficient cells pollute less per
    offspring.
- **Parity:** with `Δμ0 → ∞` (a metaparameter switch `gradient=0`) and `ε = ε_ref`, the factor is 1
  and the yield is 1. That is today's model, bit for bit, and the first test.
- **Where it acts:**
  - near clean, full sites, the driving force is large and the factor is close to 1 for any
    modest `ε`;
  - in waste-laden or depleted sites, only low-`ε` (fast, wasteful) cells can eat at all.

  So `ε*` should vary in space and track the waves.

**Measurements** (all exist or are small additions):
- the distribution of `ε` over time and space;
- power per cell `ε · M'`;
- efficiency, as energy captured per unit food;
- the partition (tax / death / overflow);
- energy per unit of waste.

All are reported per unit energy, as the anchor doc requires.

## 7. Predictions and falsifiers

1. **Interior optimum.** Evolved `ε` settles strictly inside `(0, Δμ)`.
   - Falsified if it runs to either bound in every environment.
2. **Mixing decides the strategy** (Pfeiffer et al.). Lower `wdiff` (more structured) gives
   higher evolved `ε` (efficient); higher `wdiff` (well mixed) gives lower `ε` (fast).
   - Falsified if `ε*` is independent of `wdiff` at matched total flux.
3. **Maximum power or maximum efficiency (L2).** Compare evolved `ε*` with the `ε` that
   maximizes power at the evolved local `Δμ`, and with the `ε` that maximizes efficiency.
   - Lotka/Odum predicts the former.
   - The collapse test of the anchor doc applies: reach the same total flux by different routes
     (`food_inc`, `Δμ0`, `waste_sink`) and check that `ε*` depends only on the flux variables.
4. **Gradient dissipation (Schneider–Kay, with G3).** Under evolution, the source-to-sink
   difference across the colony shrinks over time against a mutation-off control.

## 8. Critique of the first draft, and the revisions it produced

1. **"The gradient is missing" was the wrong framing.** A defender of the current model would
   say that food is free energy relative to a perfect sink, so the gradient is there.
   - Revision: the claim is now that the sink is *degenerate* (infinite, stateless), and that in
     food/waste it is a hazard but not a yield. That locates exactly what to change.
2. **The first draft proposed going straight to thermodynamics in `chemca`.** That ignores N's
   "let's not gallop ahead" (Track A first).
   - Revision: G2 is designed for the food/waste model on `main`, and the `chemca` version (G4)
     is deferred to Track B.
   - This also exposed a real error in the anchor doc: `chemca` as designed would not have made
     L2 testable.
3. **"More realistic, therefore more evolvable" was an unstated assumption.** The prior data cut
   both ways: `fw2` shows rate escalation can destroy the population.
   - Revision: evolvability is a hypothesis with falsifiers (§7), not a design goal. G2 is
     justified by L2 alone.
4. **Confounding.** The coupling trait changes per-capita flux and population density at once.
   - Revision: every comparison is at matched total flux, with the density/turnover decorrelation
     the anchor doc requires (prediction 3).
5. **A numerical hazard.** `ln(F/W)` diverges as `W → 0` and is undefined at `F = 0`.
   - Revision: the regularizer `δ`, and a clip at zero mouthful when `Δμ ≤ ε` (no reverse
     reaction in stage 1).

## 9. References

Pfeiffer et al. (2001) is in `Docs/references/`:
[`pfeiffer2001_atp_pathways_rate_yield.pdf`](references/pfeiffer2001_atp_pathways_rate_yield.pdf)
(N found a copy, 2026-10-26). A correction to it was published in Science 293:1436 (2001); that
correction is not in the repo. The others are not downloaded.

- A. J. Lotka (1922), "Contribution to the energetics of evolution", PNAS 8:147.
- H. T. Odum and R. C. Pinkerton (1955), "Time's speed regulator: the optimum efficiency for
  maximum power output in physical and biological systems", American Scientist 43:331.
- F. L. Curzon and B. Ahlborn (1975), "Efficiency of a Carnot engine at maximum power output",
  American Journal of Physics 43:22.
- T. Pfeiffer, S. Schuster and S. Bonhoeffer (2001), "Cooperation and competition in the
  evolution of ATP-producing pathways", Science 292:504. DOI 10.1126/science.1058079.
- E. D. Schneider and J. J. Kay (1994), "Life as a manifestation of the second law of
  thermodynamics", Mathematical and Computer Modelling 19:25.
- H. J. Morowitz (1968), *Energy Flow in Biology*, Academic Press.
