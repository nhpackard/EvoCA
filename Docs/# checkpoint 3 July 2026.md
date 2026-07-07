# checkpoint 3 July 2026

High level comments for Claude's consideration.

The purpose of this document is to take a step back from the details of the current research agenda, and consider some big picture issues that might suggest structural changes we might consider.  If we implement changes here, these will have effects on the details of our campaigns, to be sorted out afterward.

This should be considered deep analysis.  Claude should take his time, consult literature when appropriate (and include citations), and certainly check out the local repos indicated.

Claude should push back on expressed intuitions, whenever appropriate, especially with reasoned responses containing alternative perspectives.  Claude's response should be appended to this document.

Each section should be considered a potential dev path, at least with a series of commits, possibly with a branch, commits, PR.

# metaparam structure

Metaparams should be divided into two groups, and this division should be reflected in two indpendent metaparam structures:  one for environment, one for agents. E.g. a partial example:
```
env_metaparams = {
    'N':512,
    'food_inc':20,
    'gdiff: 0.8,
}
agent_metaparams = {
    'mu_lut':...,
    'mu_egene':...,
    ...
}
```
The point:  I want to think about them differently, and explore them (somewhat) independently.  I am thinking of the different environmental metaparams as possible different niches in Nature, ranging from arctic cold to equatorial jungle.  Also:  agent_metaparams are all candidates for being themselves mutated as part of the evolutionary process.  Agent metaparams should be contained in each agent, even if some, for the moment, are constant across all agents.   Besides variations caused  by mutation, we might also explore simply constructing populations with a distribution of values for one or more of these params, to see if any values seem to be dominant.

# Evolutionary patterns to seek

## Evolution as a function of driving by resource

Roughly, I would like to have a resource parameter that has the properties
- low resource => everything dies, or a small population barely survives => not much evolution is going on.
- high resource => everything lives, not much evolution because everything is similarly fit since getting resource is so easy.
  - an exception:  even with lots of resources, there could still be nontrivial evolution of egenomes.
  - also, even without evolving egenomes, there might still be evolution of LUT for faster eating => reproductive edge.
  - but in the high resource regime, I'm thinking the fitness landscape is rather flat, with mild hills, rather than having lots of high peaks and low valleys.

So the intuition is that 
- low resource driving =>  low evolution
- high resource driving => low evolution
- intermediate resource driving => high evolution

Is this intuition valid??


## Resource driving

It may be that our current resource driving with its two metaparams, food_inc and gdiff, are not rich enough  to provide resource parameterization intuitively described above.

I have the feeling that global increment of all food at every time step is not creating enough opportunity for agents to be stressed, die, and be replaced, enabling evolutionary opportunity.  

One simple additional parameter could be a food update period, e.g. `food_per`.  Now, food_per=1, updated every time step.  If we made it > 1, then patterns in the food field would emerge due to the differential eating capabilities of the agents.  The would add stochastic contingency, but might also provide structure for evolution of dynamics to work for/against.  Worth noting, in this regard that in ~/Projects/Bugs, the sensory-motor agents are consuming food (as here), and they definitely create spatial structure that affects their evolution.

## Innovation + optimization

There is, in the literature (Rasmussen et. al.) an intuition that evolution might have large innovative steps, followed by optimization of parameters governing the innovation.  This is essentially a form of punctuated evolution (one of the many possible explanations). 

Seems like it should be possible to observe this pattern in EvoCA, especially in the evolution of the egenome.  If all the egenes are turned off initially, turning one of them on makes a big difference.  Subsequently, adjustment of the exact value of the egene is more like the optimization phase.  

An essential question for such studies is what serves as fitness, or a proxy thereof?  Could be population level itself, the usual definition of fitness in population biology.  Could also be direct measures of food intake, etc.

This (egenome example) might be best studied, i.e. focused on, in the extreme case mentioned above: the constanly resource rich environment.  Possibly with all eating patterns present in the egenome, so it is really a task of discovering the best egenes without stochasticity introduced by exploring the LUT, since all LUT results yield good food because of the resource rich environment.

Another special case for egenome study:  Fix the LUT, with mu_LUT=0, and evolove only the egenomes.  It would be tempting to fix LUT with GoL, but this is not good because GoL dynamics converge to uninteresting fixed point + local oscillations.  Maybe we need to find a slightly more active proxy than GoL.  A quick check indicates that 100 steps of GoL initial LUT with mutation, then mu_LUT=0, seems to result in ongoing dynamics (to be confirmed).

# Resource complexity

We currently have a single  variable (scalar) resource field.  One way to open up the space of evolutionary paths is to make the resource field richer, with more possible values, i.e. to have a vector valued resource.  The multiple values could represent compounds interacting with each other via chemical reactions.

This additional complication entails significant input on the design of the chemistry.  Possibilities are 
- simplified autocatalytic network of reactions (e.g. Kauffman style)
- random reactions from Swarm Chemistry (cf ~/Projects/swarmchem)
- either of the above with addition of energy activation with a special compound (analogous to phosphorus in real chem), to boost reaction rates.


## implementation

I would be inclined to make a semi-permanent branch of the code containing the vectorial resource (food).  Maybe call it simply chem.  Eventually if this branch becomes dominant in the research, it could be forked to it's own repo.  I have some collaborators that might be interested in this branch more than the current main branch.  

I am thinking of the chemistry as being specified *a priori* in the model, and being fixed for all time.  I see 2 possibilities for chemistry definition
- stocahstic growth of chemistry with new compounds and new reactions defined on the fly, analogous to 'mutations' of the chemistry, with a corresponding mutation param.  This is the simulation style of Farmer Kauffman Packard paper on evolution of autocatalytic networks.
- pre-defined chemistry, choosing a universe of $N_\text{chem}$ compounds with a network of reactions and corresponding rates defined at the outset.  'Growth' of the chemistry within this framework would be obtained by initializing a sim with a small number of compounds with nonzero concentration, and allowing the pre-defined reactions to gradually populate the space of compounds and reactions.

Efficiency might lead us to consider agent-local food vectors, that would grow as reactions happen, but I would be inclined to allocate memory for a resource plane for each compound.  For a world with N=512, one resource lattice would be 512x512 = 2^18 = 262,144 numbers, so 100 compounds => 26,214,400, 1000 compounds => 252,144,000 values, etc.  Seems like we don't run into memory problems, even for 10,000 compounds.

# Activity

## Relevant activity

I am worried that we are seeing a bunch of activity that does not reflect actual phenotypic evolution, i.e. change in how agents react (functionally) to their environment.

I wonder if there might be a form of activity that has the same spirit as our "restricted_mu"; only count those LUT entries that actually get used, not simply all those that are present.

## Use of neutral model

Currently, we are, in the style of Alistair Channon, computing average activity of a neutral model and counting activity only above that threshold (which changes over time).

But isn't it the case that each gene present has a neutral shadow?  So we could threshold each gene separately.  There would, of course, be more fluctuations in the estimate of the neutral threshold, but the increase in statistical fluctuation must be weighed against sharper estimate per gene.

We have, perhaps discussed this in the past, during the implementation of our current neutral model version.  Perhaps I just need to be reminded of that conversation.

---

# Claude analysis (2026-07-07)

Framing note: this is a deep-analysis response to the checkpoint above.
It is grounded in the *current* code (I re-read `C/evoca.c`,
`python/evoca_py.py`, `python/evoca_explore.py`, `Docs/probes.md` —
the 50-day-old project memories were stale in several places, flagged
below), in the two local repos you pointed at (`~/Projects/Bugs`,
`~/Projects/swarmchem`), and in the literature (citations inline, new
entries appended to `Docs/references/references.md`). Where I think an
expressed intuition is wrong or under-specified I say so with a reason,
not a hedge. Each lettered section is a candidate dev path.

## Shape of this pass

| § | Checkpoint topic | What this section does |
|---|---|---|
| A | metaparam structure | Endorse the split; show it is *two* changes (a zero-risk refactor + the per-cell mutable genes); fix where the boundary actually falls |
| B | resource-driving inverted-U | Validate the *shape*, correct the *mechanism* — high resource does not stop turnover, it makes evolution drift-dominated |
| C | resource driving richness / `food_per` | Endorse; identify the real mechanism (depletion-structure lifetime) and the deeper lever (regen toward a spatial template — the Bugs mechanism) |
| D | innovation + optimization (punctuation) | Endorse; Coreworld precedent; separate the two egene experiments you conflated; fix the fitness proxy and the frozen-substrate metric trap |
| E | resource complexity / chemistry | Endorse the branch; swarmchem is a *negative* reference; adopt your own "fixed universe + sparse seeding" option and go sparse-vector not dense-plane |
| F | relevant activity (used entries) | You already built this — it is `dyn_activity`; promote it to the primary phenotypic signal |
| G | per-gene neutral model | Remind you of the prior decision; the canonical per-component statistic is still missing; full per-gene shadows are ill-posed in the open LUT space, well-posed (do it) in the fixed spaces |
| H | cross-cutting | The lineage field (now opt-in in code) is the shared backbone |
| I | dev logistics | Branch-vs-main per path |
| J | structural & scaling review | Ownership, boundaries, scaling risks named |
| K | roadmap + self-critique | Prioritised list; the flaws I caught revising this |
| L | energy gradient (added later) | Verdict: unifies §B/§C/§E and does **not** perturb them; the gradient already exists as a scalar through-flux, and is missing only as *structure* and as *measurement* |

## Load-bearing claims (read these even if you skip the body)

1. **The env/agent split is really two changes, and only one is hard.**
   Splitting the *namespace* (two dicts) is a zero-risk refactor; do it
   now. Storing agent-params *per cell so they mutate* is the valuable,
   invasive change — and it is exactly the already-planned `mu-genome`
   branch (research_directions §10). The SoA memory layout makes it
   cheap (one `calloc` + one copy line per gene). The niche metaphor
   assigns the boundary cleanly *except* the tax family, which
   genuinely splits: baseline `tax` is an environmental cost-of-living;
   `tax_per_egene`/`tax_lut`/`tax_ring` are agent-side self-imposed
   metabolic costs.

2. **Your inverted-U intuition is right in shape but the mechanism at
   the high-resource flank is not "no evolution" — it is "drift-
   dominated evolution."** In EvoCA the reproduction rule (replace the
   lowest-food Moore neighbour whenever `f≥1`) keeps turnover high even
   when food is superabundant; what collapses is the *fitness variance*
   (Fisher's fundamental theorem: response to selection ∝ additive
   variance in fitness), so genomes turn over by drift, not selection.
   That matters because raw activity *stays high* on the flat flank —
   which is precisely why you cannot read adaptive progress off raw
   activity, and why §F/§G are not a side-quest but the instrument this
   whole intuition needs. Your scans already located the productive
   peak (`food_inc ≈ 0.013`).

3. **`food_per` is worth adding; its real job is to lengthen the
   lifetime of self-generated food-depletion structure.** Bugs
   (`~/Projects/Bugs`) generates its evolutionarily-relevant spatial
   structure by grazing craters that regrow slowly toward a fixed
   source field. EvoCA's every-step spatially-uniform increment erases
   depletion structure as fast as agents create it. `food_per>1` is one
   fix; the deeper lever is regrowth toward a spatial *template*
   (Bugs' `F_source`), for which EvoCA already has the scaffolding
   (`env_mask`, `F` clamped to [0,1]).

4. **The chemistry branch is a good idea, but swarmchem will not help
   and your own second definition-option is the right one.**
   `~/Projects/swarmchem` is a kinematic flocking model with *no
   reaction network* — a negative reference. The anchor is your own
   Farmer–Kauffman–Packard (1986) autocatalytic-sets model. Your
   "fixed universe of N_chem compounds, seed sparse, let reactions
   populate" option unifies pre-defined vs grown chemistry *and*
   justifies sparse agent-local food vectors over the dense
   per-compound planes you defaulted to — which is the difference
   between a tractable and an intractable first cut.

5. **On activity: half of what you asked for already exists, and the
   other half you should not build the way you propose.** "Count only
   used LUT entries" is `dyn_activity` (500 used-transition buckets,
   updated only for the entry that actually fired) plus the existing
   per-step `lut_active` mask — promote it to the primary phenotypic-
   evolution signal. "Give each gene its own neutral shadow" is
   ill-posed in the open LUT-hash space (a just-born genome has no
   drift history to estimate an expectation from) but well-posed in the
   fixed spaces (eg=729, dyn=500) — do per-bucket z-scoring there. The
   canonical middle option — Bedau–Packard per-component *new activity*
   against one shared neutral threshold — is still not implemented and
   is the cheapest real upgrade.

6. **The energy-gradient worry (§L, added later) does not overturn any
   of the above — it unifies them.** EvoCA already has the essential
   thermodynamic ingredient: a source→sink free-energy through-flux
   (`food_inc` creates food, tax + death + cap-overflow destroy it,
   agents extract reproduction from the flow). What is missing is that
   the gradient is nowhere *measured* and is *unstructured*. Reading the
   model this way sharpens §B (the real control parameter is the
   source/sink *ratio* ≈ `food_inc/tax`, and the high-resource flat
   flank is literally gradient-collapse toward equilibrium), and
   retroactively justifies §C's source-field template (a spatial
   gradient) and §E's activation compound (an internal ATP-like
   gradient). Caveat: EvoCA agents are sessile, so the "navigate to the
   vent" intuition maps only weakly — the internal/chemical gradient is
   the better fit.

---

## A. Metaparam structure: env vs agent — endorse, with a sharper boundary

I endorse the split; the niche metaphor (arctic ↔ jungle as points in
env-metaparam space) is a clean organizing idea and matches how these
knobs actually function. Two things to get right.

**It is two independent changes, and conflating them will cost you.**
Ground truth (current code): all 14 metaparams are global `static` C
variables (`C/evoca.c:150-179`: `gfood_inc, gm_scale, ggdiff, gmu_lut,
gmu_egene, gmu_egenome, gp_dup_egene, gtax, gtax_per_egene, gtax_lut,
gtax_ring, grestricted_mu, …`), and the Python side is a single flat
namespace — `init(...)` takes them as flat kwargs and `params()` emits
a flat dict (`python/evoca_py.py`). So:

- **Change 1 — namespace split (do now, zero dynamical risk).** Group
  the *Python* representation into `env_metaparams` / `agent_metaparams`
  dicts, with `params()`/`params_str()`/`.evoca` export deriving from
  them. No C change, no behaviour change. This is pure ergonomics and
  makes the niche-vs-organism exploration you describe natural to
  script. Low cost, immediate payoff.

- **Change 2 — per-cell agent genes (the real one).** "Agent
  metaparams should be contained in each agent, even if some are
  constant for now" is a genome-structure change: today no agent-param
  is stored per cell. This *is* the mutation-of-mutation direction
  already scoped as the `mu-genome` branch (research_directions §10),
  and the good news from re-reading the code is that it is cheap: the
  per-cell layout is Struct-of-Arrays (separate `calloc`'d arrays,
  `C/evoca.c:801-825`), so each new per-cell gene is one array + one
  copy line in the Phase-4 reproduction block — not a struct refactor.
  Your "even if constant for now" instinct is exactly right: put the
  storage in first with mutation rate 0, so turning on evolvability
  later is a one-line change, not a data-model migration.

**Where the boundary actually falls — and the one place it doesn't.**
The niche metaphor makes most assignments obvious: environment =
`{N, food_inc, gdiff, food_per(new), env_mask/template(new)}`; agent =
`{mu_lut, mu_egene, mu_egenome, p_dup_egene, m_scale, restricted_mu,
LUT, egenes}`. The exception worth designing around is **the tax
family, which splits across the boundary**:

- baseline `tax` is environmental — the cost of living imposed by the
  niche (a harsh/arctic niche = high tax); it belongs in `env`.
- `tax_per_egene`, `tax_lut`, `tax_ring` are *agent-caused* costs — the
  organism's own genome complexity is what it pays for. Conceptually
  these are agent-side (metabolic overhead of the phenotype the agent
  carries), even though the *rate* is currently a global policy.

That distinction is not pedantic: if you later make tax-rates evolvable
(a plausible extension), only the agent-caused ones make sense as
per-cell genes; the baseline cost-of-living should stay an
environmental control. So the split should be by *causal owner*
(who imposes the cost), not by whether the word "tax" appears.

**The distribution-seeding experiment you mention needs Change 2
first.** "Construct populations with a distribution of values for one
param, see if any value dominates" requires per-cell storage of that
param — it is not doable with the current globals. Once the per-cell
array exists (even with mutation off), seeding a spatial or random
distribution of, say, `m_scale` and watching which value's carriers
win is a zero-mutation selection assay — and a clean, cheap precursor
to the full evolvable-rate experiment. This is the natural *first*
use of the `mu-genome` machinery and I'd run it before enabling
rate mutation.

## B. The resource-driving inverted-U: right shape, wrong high-flank mechanism

The intuition — low resource ⇒ little evolution, high resource ⇒
little evolution, intermediate ⇒ much — is a real and well-supported
*shape*. It is the ecological hump-shaped productivity–diversity
relationship (Rosenzweig 1995), and it is the dissipative-structure
"window" (Schneider & Kay 1994; Chaisson 2001) already in the
references: complexity is selectable only where flux is high enough to
pay for it and low enough that it still confers a differential. Your
scans have already *found* this window empirically — `food_inc`
bracketed tightly at ≈0.013. So the qualitative answer to "is this
intuition valid?" is: yes, and you have already measured the peak.

**But the two flanks are not symmetric, and the high flank is where
your own phrasing is misleading.** The low flank fails by *extinction /
no standing variance* (nobody survives to vary; the population is a
handful of cells barely clearing tax). The high flank does **not** fail
by "everything lives so no evolution." Look at the reproduction rule:
a cell reproduces whenever `f ≥ 1` and *always* kills the
lowest-food Moore neighbour. That means **turnover and spatial
competition persist at every resource level** — superabundant food does
not switch off death, it just makes everyone reach `f=1` quickly. What
collapses at high resource is the *variance in reproductive success*:
when even a mediocre egene/LUT gets enough food, the fitness
differential between genomes → 0, and by Fisher's fundamental theorem
the adaptive response → 0. Evolution does not stop; it becomes
**drift-dominated**. Genomes keep turning over, hashes keep churning,
raw activity stays high — but none of it is adaptive.

This is not a nitpick; it reframes the whole measurement problem and
connects directly to §F/§G:

- Your two "exceptions" are correct and now have a mechanism. "Even
  with lots of resource there could still be egene evolution" and "even
  without egene evolution there could be LUT evolution for faster
  eating" are both statements that *some* fitness variance survives at
  high resource — via cognitive efficiency or eating speed even when
  raw survival is trivial. Whether they actually do is the
  experiment, and the right readout is **per-lineage net metabolic
  return** (cumulative mouthful − cumulative tax; research_directions
  §6), not population count.
- "In the high-resource regime the landscape is flat with mild hills"
  is exactly the drift regime. The danger is that raw activity will
  *look* alive there. So the high-flank test is the sharpest possible
  motivation for the neutral-model work: **the flat landscape is
  detectable only as raw-activity-high but excess-activity≈0.** If you
  run a `food_inc` sweep and plot excess (component-normalised, and the
  fixed-space `dyn_excess`) against `food_inc`, the dissipative-window
  prediction is an inverted-U in *excess* that the raw activity will not
  show. That single plot both tests your intuition and validates the
  metric — do it early.

Predicted result, stated so it can be falsified: raw activity roughly
monotone-increasing or flat in `food_inc` above the extinction
threshold; *excess/dyn_excess* inverted-U peaking near 0.013; per-egene
cognitive load and per-lineage net return also inverted-U, collapsing
(not merely plateauing) on the high side. If instead cognitive load
keeps *rising* at high `food_inc`, that falsifies the "flat landscape"
claim and would be the more interesting outcome (cognition as a
resource sink that pays for itself even under slack survival selection).

## C. Resource driving richness / `food_per`: endorse; the mechanism is depletion-structure lifetime

Strong endorse — I think this is the highest value-to-cost structural
change in the document, and the Bugs comparison you gestured at is
exactly the right evidence.

**What Bugs actually does (from the code), and why it evolves richly.**
In `~/Projects/Bugs`, a bug eats *only from its own cell* (`bugs.c`
eat step), and food regrows by slow additive relaxation toward a
*fixed source field* `F_source` (`F ← min(F+food_inc, F_source)`),
optionally diffused. Because regrowth (`food_inc≈0.01`, ~100 ticks to
refill) is far slower than a bite, grazing carves persistent depletion
craters and trails, and the food field develops low-wavenumber spatial
structure. The bugs' *perception* (per-direction food thresholds) then
co-adapts to that self-generated structure over ~10k ticks. The
spatial structure is the substrate evolution works on.

**Why EvoCA under-generates this.** EvoCA regenerates food every step
by a spatially *uniform* increment. Uniform + every-step means
depletion structure is topped back up almost as fast as agents create
it — the field is kept near-flat, so there is little spatial substrate
for dynamics to exploit. (Notably: Bugs *also* regenerates every step —
so the difference is not the period, it is that Bugs regrows toward a
capped *template* while EvoCA adds uniformly. Worth being precise about
this, because it changes the recommendation.)

**So the real lever is the eating-rate ÷ regrowth-rate ratio and
whether regrowth is spatially structured — not the period per se.**
`food_per>1` helps, and its precise mechanism is that it **lengthens
the lifetime of self-generated depletion patterns**: between
replenishments, agents deplete the field and the resulting spatial
structure *persists* for `food_per` steps instead of being erased next
tick. That is genuinely useful and cheap (one modulo in Phase 2 of the
step loop). But I'd pair it with the deeper knob:

- **`food_per` (temporal, trivial):** regen only every K steps. Adds
  temporal contingency and extends depletion-structure lifetime. One
  line in `evoca_step`.
- **Regrowth-toward-template (spatial, the Bugs mechanism):** replace
  the uniform additive increment with `F ← min(F + food_inc·rate,
  F_source(x))` for a per-cell source field. EvoCA already has the two
  ingredients — `env_mask` (a per-cell 0/1 regen gate) generalises to a
  float `F_source`, and `F` is already clamped to [0,1] (a ceiling).
  This is what actually makes depletion craters *stable* rather than
  transient, because a grazed cell in a low-`F_source` region stays
  poor. It is also how you would later paint niches (arctic patch vs
  jungle patch) into a single world — which ties back to §A's
  environment metaparams.

Recommendation: add `food_per` first (it is nearly free and directly
tests "does temporal lumpiness create evolutionary opportunity?"), then
the source-field generalisation of `env_mask` as the mechanism that
makes spatial food structure persistent. Both belong on the same small
branch. Test the claim with the `ts` `F_env` variance and the spatial
correlation-length probe: the prediction is that mean `F_env` can be
held fixed while its *spatial variance* rises with `food_per` and with
source-field contrast — and that turnover / excess-activity rise with
that variance. If evolution does **not** respond to food-field
structure, that is an important negative result about how coupled the
CA dynamics really are to the resource field.

## D. Innovation + optimization (punctuation): endorse, with three corrections

This is a well-posed and literature-anchored direction. The
"large innovative step, then optimization of the parameters governing
it" pattern is punctuated equilibrium (Eldredge & Gould 1972), and it
has direct artificial-life precedent: the Coreworld (Rasmussen,
Knudsen, Feldberg & Hindsholm 1990, *Physica D* 42:111-134) evolved
through *seven successive epochs* of cooperative structure, each a
punctuation; Lindgren (1991) saw the same epochal structure in evolving
iterated-game strategies. Your egene model is arguably a cleaner
punctuation generator than either, because activating a previously-off
egene slot is a discrete, identifiable innovation event. Three
corrections.

**(1) You described two different experiments as one; keep them
separate.** "Turn one egene on ⇒ big difference, then tune its value"
is the *punctuation* experiment: it *requires starting with genes off*
so activation is the innovative step. "Resource-rich, all eating
patterns present, discover the best egenes without LUT stochasticity"
is the opposite setup — with all slots already on there is no
activation event to punctuate, only optimization. Both are worth
running, but they answer different questions:
- *Punctuation study:* seed `set_egenome_pair_all` with one active
  centre-only slot (minimal cognition), watch the `egene` `spec`/`load`
  strips and the `eg_activity` 729-space for step-changes when a new
  slot activates. The signature is a staircase, not a ramp.
- *Optimization-only study:* seed with many/all slots active in a
  resource-rich, LUT-frozen world and measure how fast max-match climbs
  to optimum. This isolates the "tuning" phase.

**(2) The fitness proxy: use two, and prefer the lineage one.**
Population level is the classical proxy and cheap, but it is confounded
(population also tracks resource, not just adaptedness). The direct
proxy is **per-lineage net metabolic return** (cumulative intake −
cumulative tax), which the now-merged opt-in **lineage field**
(`lin_parent_hash/birth_id/parent_id`, `C/evoca.c:230-232`, enabled via
`evoca_enable_lineage()`) finally makes measurable — a staircase in
lineage net-return is the crisp punctuation signal. Use population as
the coarse readout and lineage net-return as the causal one.

**(3) The GoL-proxy problem — your instinct is right, but mind the
frozen-substrate metric trap.** You correctly note GoL freezes to
fixed-point + local oscillators, an uninteresting substrate for
egene-only evolution, and that "100 GoL steps with mutation, then
`mu_lut=0`" gives ongoing dynamics. That "evolve-then-freeze" recipe is
exactly campaign #3c, which found it *does* yield a viable active
frozen substrate (unlike freezing raw GoL, which went 10/10 extinct).
But there is a trap I have to flag because it will bite this experiment:
on a *frozen evolved* substrate (`mu_lut=0`, rich rule) the
`dyn_excess_pc` metric becomes an **artifact** — with LUT flux driven
to zero the neutral baseline collapses while the alive population keeps
exercising the rich frozen rule, so `dyn_excess_pc` explodes (#3c arm A
≈ +335, not biology). Rule: for any LUT-frozen arm, measure egene
evolution with the **eg (729) shadow and lineage net-return**, and do
**not** interpret `dyn_excess_pc` (or the LUT-hash excess, which is
LUT-only and undefined for egene-driven runs anyway). This is a
known result — I'm reminding rather than re-deriving.

## E. Resource complexity / vector chemistry: endorse the branch, redirect the sources

I endorse a semi-permanent `chem` branch and the eventual fork — a
vector-valued resource genuinely opens the evolutionary path space, and
a reaction network is the natural way to make food *earned through
transformation* rather than merely *harvested*. Three redirections.

**swarmchem is a negative reference — do not mine it for machinery.**
I read it: `~/Projects/swarmchem` is a port of Sayama's Swarm Chemistry,
a *kinematic* multi-agent flocking model. "Chemistry" there is a
metaphor for mixing recipes of behavioural parameter sets; there are no
compounds, no reactions, no rates, no network, and the evolutionary
operators aren't even ported. It has nothing to contribute to a
reaction-network resource field. The real anchor is your own
**Farmer, Kauffman & Packard (1986), "Autocatalytic replication of
polymers," *Physica D* 22:50-67** — catalyzed cleavage/condensation
with a critical-diversity threshold for autocatalytic-set emergence —
and the Coreworld (§D) for how such a chemistry behaves under spatial
selection.

**Adopt your own second definition-option; it dominates the first for a
first cut.** You offered (i) stochastic on-the-fly growth of chemistry
(FKP-style) and (ii) a pre-defined universe of `N_chem` compounds +
reactions, "grown" by seeding a few compounds and letting the fixed
reactions populate the space. Option (ii) is strictly better to build
first: it is reproducible (critical for the exponent/criticality work),
has no meta-mutation-rate to tune, and — importantly — its "growth by
population of a fixed universe" gives you the emergence phenomenology of
(i) without the bookkeeping of dynamically allocating new
species/reactions mid-run. You can always add on-the-fly chemistry
mutation later as a second knob. (This also mirrors the EvoCA design
philosophy already in the repo: fixed finite spaces — the 729 egene
keys, the 500 dyn buckets — bought you closed-form neutral baselines;
a fixed chemistry universe buys the same rigor.)

**Go sparse-vector, not dense-plane — the scaling argument.** You
defaulted to "a resource lattice per compound" and reasoned from
capacity (10k compounds × 512² ≈ 2.6·10⁹ floats ≈ 10 GB — fits).
Capacity is not the binding constraint; **per-step compute and memory
bandwidth are.** A dense scheme costs O(N_chem·N²) memory *touched
every step* and O(R·N²) for R reactions — at N_chem=1000, N=512 that is
a quarter-billion cells streamed per step, which will dominate the CA
core (currently ~150 fps at N=256) and make the fps target unreachable.
But under option (ii)'s "seed sparse, let it populate" dynamics, **most
compounds are absent at most sites almost all the time** — the
occupancy is sparse. So store food per cell as a *sparse vector*
(active-compound list, like the egene slots already are), and iterate
reactions only over co-present reactant pairs. This is O(active
compounds) not O(N_chem), and it is the same sparse-slot pattern the
egene model already uses. Memory drops from "plane per compound" to
"a few slots per cell." This is the single design decision that decides
whether `chem` is tractable at N=512.

**One design gap to resolve before coding:** you have not said what an
agent *does* with a vector resource — how multi-compound concentration
maps to the scalar private food `f` that gates reproduction, and how
the egene match generalises (does an egene now specify a target
*compound profile* in the neighbourhood, or still a spatial pattern?).
That mapping is the actual scientific core of the branch (it defines
what "eating well" means chemically), and it is under-specified in the
checkpoint. I'd pin it down first — everything else is plumbing.

## F. Relevant activity (used entries): you already built it — promote it

Your worry — "we may be seeing activity that does not reflect
phenotypic change" — is correct, well-founded, and *already
instrumented*. The full-genome `activity` probe keys on an FNV-1a hash
of the whole LUT (`C/evoca.c:747-788`), so it counts silent churn:
flipping a LUT bit for an input the cell never encounters changes the
hash (a "new species") without changing behaviour. That is the churn
you rightly distrust.

What you're asking for already exists as **`dyn_activity`**: 500 buckets
over `(LUT-input, output)` transitions, incremented in Phase 1 **only
for the specific entry that actually fired, and only for alive cells**
(`C/evoca.c:1155,1162-1167`). It is a per-*used-entry* activity by
construction — exactly your "count only entries that actually get used"
proposal, at the transition granularity. The per-step `lut_active`
mask (`C/evoca.c:1154-1177`) is the same idea at the whole-mask level
and is what `restricted_mu` already uses to keep mutations phenotypic.

So the recommendation is not "build it" — it's **promote `dyn_activity`
(and its fixed-space `dyn_excess`) from a side-probe to the primary
phenotypic-evolution readout**, and demote full-genome `activity` to a
churn diagnostic (its value is precisely the *comparison*: `activity`
high while `dyn_activity` concentrated ⇒ silent churn; both moving ⇒
real functional change).

The one genuinely new artifact worth building is the **genome-level**
analog of the transition-level idea: an activity that keys on the cell's
LUT *masked to the entries it actually exercises* (the "effective rule"
over the inputs that occur in that cell's realised neighbourhood
distribution), so two genomes differing only in never-used entries
collapse to one species. That is a strict upgrade over full-genome
`activity` and would directly answer your worry at the species level
rather than the transition level. It reuses the existing per-step
`lut_active` machinery; the only new part is per-cell (not global)
used-masking, which is a modest change.

## G. Per-gene neutral model: the reminder you asked for, and where your instinct holds

You asked to be reminded of the prior conversation. Here it is, and it
bears directly on your proposal.

**What EvoCA computes today (three things, none of them per-gene).**
(1) The Channon-style **LUT-hash shadow** (`C/evoca.c:1778-2036`) is a
*single* shadow population that mirrors the real run's birth/death
counts under random (non-selective) inheritance; excess is computed in
Python as a *pooled scalar* — raw `ΣG − ΣN` or component-normalised
`ΣG/D_G − ΣN/D_N` (`evoca_explore.py:284-308`). (2) & (3) The two
**fixed-space drift shadows** (eg=729 ternary keys, dyn=500 transitions)
maintain a per-bucket Wright–Fisher drift occupancy driven by the
*realised* mutation flux, and expose a per-component-normalised scalar
each. None of the three applies a *per-gene threshold*.

**Your proposal — "each gene has its own neutral shadow, threshold each
gene separately" — is the right instinct in fixed spaces and an
ill-posed one in the open LUT space.** The reason is statistical, and
it is the reason the current design is pooled:

- In the **open LUT-hash space**, a genome that just arose has *no
  history* — you cannot estimate its neutral expected activity from its
  own drift, because there is exactly one realisation of it and it is
  brand new. The pooled/shared-threshold estimator exists precisely to
  *borrow statistical strength across genes*: it estimates one neutral
  activity distribution from the whole shadow ensemble and thresholds
  every gene against it. Per-gene shadows in an open space trade a
  little bias for catastrophic variance (each gene's "shadow" is a
  single sample). This is the crux of the earlier conversation, and
  it's why the LUT shadow was built pooled.

- In the **fixed finite spaces (eg 729, dyn 500)**, your instinct is
  exactly right and *should* be implemented. Each bucket is present for
  the whole run, so it has a well-defined drift trajectory, and the
  Wright–Fisher baseline gives a **closed-form per-bucket variance** —
  so you can compute a per-bucket z-score (observed activity − drift
  mean)/drift-sd and threshold each key/transition individually. Today
  the fixed-space shadows only expose an *aggregate* per-component
  scalar; adding per-bucket z-scoring is the concrete, well-posed
  version of your proposal, with the extra fluctuation you anticipated
  ("more fluctuations in the estimate... weighed against sharper
  estimate per gene") bounded analytically rather than guessed. This is
  a real upgrade and I recommend it.

**And the canonical option that sits between pooled and per-gene is
still missing.** The Bedau–Packard *new-activity* statistic A_new —
count the components whose cumulative activity exceeds a *single*
neutral-derived threshold — is a per-component test with a shared
threshold, and it is the actual literature statistic (research_directions
§2.1 already recommended it; Bedau, Snyder & Packard 1998). EvoCA has
the pooled difference and the fixed-space scalars but not A_new. It is
the cheapest rigorous improvement and it is turnover-invariant, which
directly de-fangs the high-resource drift problem from §B. Order of
work: A_new first (cheap, canonical), per-bucket z-scoring in the fixed
spaces second (your instinct, well-posed), per-gene open-space shadows
**not at all** (ill-posed — the pooled shadow is the correct estimator
there).

## H. The cross-cutting backbone: the lineage field

Three of the sections above (§B per-lineage net return, §D punctuation
fitness proxy, §C/§F "who actually out-reproduced") reduce to
parent→child bookkeeping. That field now exists in code as an opt-in
SoA addition (`lin_parent_hash/birth_id/parent_id`, zero-cost when off).
It remains the single highest-leverage instrument; most of the
"is this adaptive or drift?" questions this checkpoint raises are only
answerable through it. Nothing to build — just *use* it (enable it in
the §B/§D campaigns and add a lineage-net-return reducer alongside the
existing metrics).

## I. Dev logistics (branch vs main)

| Path | Where | Rationale |
|---|---|---|
| §A.1 namespace split | **main** | Python-only, no dynamics change, immediately useful |
| §A.2 per-cell agent genes | **branch** `mu-genome` (already scoped) | genome-structure change; a half-tuned rate-gene contaminates every run |
| §B inverted-U sweep + excess plot | **main** | pure config + analysis on existing code |
| §C `food_per` + source-field | **branch** `resource-driving` | changes Phase-2 of the step loop → alters every run's dynamics |
| §D punctuation / optimization campaigns | **main** | config + lineage reducer; no C change |
| §E vector chemistry | **branch** `chem` (semi-permanent, fork-candidate) | large, invasive, and you want it shareable with collaborators separately from main |
| §F promote `dyn_activity`; effective-rule species | **main** (analysis) / small **main** C add for effective-rule masking | reuses existing machinery |
| §G A_new + per-bucket z-scoring | **main** | additive metric code; old runs stay reproducible |

Same discipline as before: dynamics-altering paths (§A.2, §C, §E) get
branches with a human merge gate and a pre-merge test notebook;
analysis/metric additions go to main because they re-interpret existing
data rather than change it. §A.1, §B, §D, §F-analysis, §G are all
low-risk and could proceed in parallel immediately; §C and §E are the
two that need branch isolation and careful validation.

## J. Structural & scaling review (zoom-out)

- **Artifacts / single owner.** `params()`/the `.evoca` recipe is
  currently the single serialization owner of metaparams — *keep it
  that way through the §A split*. The risk in splitting into
  `env_metaparams`/`agent_metaparams` is creating two sources of truth;
  derive both dicts from (or serialize both into) the one canonical
  recipe, don't fork the storage. This is the one place §A could
  introduce a duplication bug.
- **Boundaries.** `food_per` and `chem` correctly belong in the C
  library step loop, not the orchestration layer; the distribution-
  seeding and niche-painting belong in Python init. That's the right
  split — flagging only so the `chem` branch doesn't drift
  reaction-network logic up into `evoca_explore.py`.
- **Scale.** The dominant scaling risk in the whole document is §E's
  dense-per-compound-plane default: O(N_chem·N²) streamed/step. The
  sparse-vector redesign is not an optimisation to defer — it is the
  precondition for the branch to hit the fps target at N=512. Every
  other proposed change is O(N²) or cheaper. `food_per` actually
  *reduces* mean per-step work (regen fires 1/K of the time).
- **Reuse.** §F should reuse the existing `lut_active`/`dyn_activity`
  machinery, not add parallel used-entry tracking; §C's source field is
  a float generalisation of the existing `env_mask` (don't add a third
  env-gating array); §D/§B should reuse the merged lineage field, not
  re-derive lineage from activity snapshots.

## K. Prioritised roadmap

**Do now (low/zero code, high information, all on main or trivial):**
1. §A.1 namespace split (env/agent dicts) — ergonomics, unblocks the
   niche framing.
2. §B `food_inc` sweep with the *excess/dyn_excess* inverted-U plot —
   tests your central intuition and validates the metric in one figure.
3. §G A_new (canonical per-component new-activity) — cheapest rigorous
   metric upgrade; turnover-invariant, defuses §B's drift problem.
4. §F promote `dyn_activity`/`dyn_excess` to primary phenotypic signal
   in the analysis (pure reinterpretation).

**Do next (scoped, mostly main):**
5. §D punctuation + optimization campaigns using the lineage field and
   per-lineage net return (mind the §D.3 frozen-substrate metric trap).
6. §G per-bucket z-scoring in the fixed eg/dyn shadows — your per-gene
   instinct, made well-posed.
7. §F effective-rule (used-mask) genome species — the genome-level
   answer to your churn worry.

**Do on branches (invasive, human merge gate, pre-merge notebook):**
8. §C `resource-driving`: `food_per` then `env_mask`→`F_source`.
9. §A.2 `mu-genome`: per-cell agent genes (storage first, mutation
   off), then the distribution-seeding selection assay, then evolvable
   rates.
10. §E `chem`: fixed compound universe + sparse agent-local vectors;
    resolve the eating/egene-mapping design *before* coding.

**Single highest-leverage cheap item:** the §B excess-vs-`food_inc`
plot — it simultaneously tests the checkpoint's central intuition,
exposes whether raw activity is misleading you in the drift regime, and
tells you whether the productive window you already found is a genuine
selection peak or a turnover artifact.

**Single most consequential structural decision:** §E's sparse-vector
vs dense-plane choice — it determines whether the chemistry branch is
tractable at research scale, and it should be made before any `chem`
code is written.

### Note on how I pressure-tested this (self-critique)

Four claims changed under an adversarial re-read before I posted this:

1. *High-resource flank.* My first draft said "high resource keeps
   selection on because reproduction always kills a neighbour." An
   adversary correctly objects: constant *death* is not constant
   *selection*. Revised to the drift-domination mechanism (turnover
   persists, fitness *variance* collapses) — which is both more correct
   and ties §B to the activity sections. (§B)
2. *`food_per` mechanism.* First draft said period-batching "only adds
   temporal lumpiness." Wrong: batching lengthens the *lifetime* of
   self-generated depletion structure, so it does interact with spatial
   structure. And the sharper point — Bugs regrows every step *too*, so
   the real difference is template-vs-uniform, not period — reshaped the
   recommendation. (§C)
3. *Chemistry definition.* First draft ranked pre-defined over
   stochastic-growth as if they were rivals. They aren't: your own
   option (ii) *is* growth-within-a-fixed-universe, and it's also what
   makes sparse vectors valid. Merged them. (§E)
4. *Tax placement in the env/agent split.* First draft put "tax" wholly
   on one side. The family actually splits by causal owner (baseline =
   environment, complexity-taxes = agent), which is the non-obvious part
   worth designing around. (§A)

## L. Energy gradient: not missing in essence — missing in structure and in measurement

**Verdict first (you asked whether this perturbs the other findings):
it does not. It unifies §B, §C, and §E under one thermodynamic reading
and sharpens §B's control parameter.** But the framing needs one
correction, because "where is the energy gradient — is it missing?"
has a two-part answer: the *scalar* gradient is already there; the
*structured* gradient, and any *measurement* of either, is what's
missing.

### The gradient already exists as a source→sink through-flux

The energy-gradient perspective is not a new lens on top of the model —
it is inherent in Prigogine's notion of a **dissipative structure**
(Nicolis & Prigogine 1977): a system held far from equilibrium by a
through-flux self-organizes, "order through fluctuations." A dissipative
structure needs a low-entropy source, a sink to dump degraded output,
and it lives in the flux between them; Schneider & Kay (1994) carry this
into biology ("nature abhors a gradient" — living systems as especially
effective gradient-degraders). EvoCA has exactly this structure, and I
confirmed it in the step loop:

- **Source:** `F_food[i] += gfood_inc` (Phase 2, `C/evoca.c:1183`) —
  food (free energy) created from nothing each step.
- **Sinks:** tax destroys private food every step (Phase 2c); a dying
  cell's remaining `f` is zeroed (death dissipates); and `food_inc`
  added past the `F=1` cap is discarded. These are the dissipation
  channels.
- **Agents are the dissipative structures in between:** they route food
  from source to sink and extract work — reproduction — from the flow.

So food is not "just energy." *With the tax sink present* it is free
energy held in a maintained disequilibrium. The cleanest proof that the
tax-sink is what makes it a gradient: delete the tax and food
accumulates to the `F=1`/`f=1` caps, the system reaches equilibrium,
the flux stops, and selection collapses — which is *exactly* the §B
high-resource flat flank. Your worry is therefore half-right: the
through-flux gradient is present; what is absent is (i) any
*measurement* of it, and (ii) any *structure* to it (it is scalar and
spatially uniform, with no chemical disequilibrium for agents to
specialize to).

### Where EvoCA sits in the lineage — Prigogine plus heredity

Being precise about *which* layer of the thermodynamic-evolution stack
EvoCA occupies matters, because it is what makes the model more than a
Bénard cell. Three layers:

1. **Prigogine (dissipative structures):** far-from-equilibrium systems
   held by a flux self-organize. This is the physics of *why* structure
   appears in a flux at all — and EvoCA's source→sink flux and its
   self-organizing reproduction fronts are already Prigoginian in
   exactly this sense.
2. **Schneider & Kay (life as gradient-degrader):** Prigogine carried
   into biology — living matter as an especially effective gradient
   dissipator.
3. **Lotka (maximum power / energetics of evolution, 1922):** the
   *Darwinian* layer — selection tunes the *rate* of gradient
   exploitation. This is the one Prigogine's framework leaves out.

The distinction that matters for us: a classic dissipative structure
**does not inherit or select** — a Bénard cell has no genome, so its
gradient-use cannot ratchet across generations. EvoCA's distinctive
content is dissipative-structure physics **plus heredity plus
selection**; it sits at that intersection. That is why the §B inverted-U
is predicted by the *Lotka* layer (a selected optimum in gradient
steepness) and not by Prigogine's alone (which would give only a
physical steady state). The energy-gradient perspective is inherent in
Prigogine; the *evolvability* of gradient exploitation — the whole point
of the model — is the layer he leaves out.

### Why this unifies rather than perturbs

- **§B (inverted-U) — sharpened.** The high-resource flat flank is
  **gradient collapse**: as `food_inc ≫ tax` the system approaches
  equilibrium, per-capita dissipation → 0, fitness variance → 0,
  evolution → drift. So the right control parameter for the §B sweep is
  not `food_inc` alone but the **source/sink ratio ≈ `food_inc / tax`
  = the gradient steepness**. Lotka's (1922) maximum-power principle —
  selection proceeds so as to maximize the sustainable energy *flux*
  through the system — predicts the productive regime is the steepest
  gradient the population can survive, i.e. an interior optimum. That
  *is* the inverted-U, now with a thermodynamic mechanism. Concrete
  revision to §B: re-plot the sweep against `food_inc/tax`, and expect
  the productive peak to track the ratio, not the absolute food input.
- **§C (spatial structure) — reinforced.** A source region plus a sink
  region is literally a spatial energy gradient; the `F_source`
  template I recommended in §C is precisely how you impose one. The
  thermodynamic reading gives that template a *purpose* beyond "make the
  food field lumpy": it is a standing disequilibrium.
- **§E (chemistry) — reinforced.** A high-free-energy compound relaxing
  to a low-free-energy product, with the agent catalyzing the reaction
  and capturing the released work — plus the activation compound you
  already proposed as a phosphorus/ATP analog — is a genuine *internal*
  free-energy gradient. Chemistry is the natural home for real
  free-energy accounting; the activation compound is retroactively
  justified here as the ATP gradient.

### The honest limitation: sessile agents cannot navigate a gradient

Two of your three examples — the racer in cool air, swimming to the
vent — are about an organism *positioning itself* in a gradient. EvoCA
agents do not move; they can exploit a spatial gradient only by
*reproducing across it* (colonization fronts flowing from source toward
sink). So the spatial-navigation intuition maps weakly. The intuitions
that map *cleanly* to a sessile CA are (a) the internal ATP/redox
gradient (→ chemistry, §E) and (b) the source/sink through-flux
gradient (→ §B/§C). Worth stating plainly so we don't build a spatial
gradient expecting foraging behavior the model structurally cannot
express.

### What to actually build (cheapest → deepest)

1. **Measure it (nearly free — do first).** Add a dissipation-rate
   readout: energy in (Σ `food_inc` over regen cells), energy dissipated
   (Σ tax collected + death losses + cap-overflow), and net throughput
   per capita. EvoCA *already* sums clamped mouthfuls
   (`g_sum_mouthful`, `C/evoca.c:375`) — the intake side is done; the
   tax/death side is a one-line accumulator. That per-capita throughput
   is Chaisson's (2001) energy-rate-density Φ; plot it against
   complexity (cognitive load, ring level) to test the
   dissipative-structure claim directly (this is research_directions §6
   made concrete). It also hands §B its principled x-axis.
2. **Impose a spatial gradient (env-side, on the §C branch).** Source
   field high on one region, tax (or a second sink) higher elsewhere —
   or a monotone `food_inc` gradient across the lattice. This creates a
   persistent standing disequilibrium and a reproduction flux from
   source to sink; it is evolvable because a more efficient lineage
   extends the viable distance from the source. This is the
   thermodynamically-motivated version of the §C template.
3. **Chemical free-energy currency (the `chem` branch, §E).** The
   deepest and most realistic: free energy carried by compounds,
   released by agent-catalyzed reactions, rate-gated by an
   activation/ATP compound. This is where "energy gradient" stops being
   an analogy and becomes an accounting identity in the model.

### One caution against over-claiming

A gradient is a *sustaining* condition — it holds the system away from
the equilibrium that kills selection — but it is not *sufficient* for
open-endedness. A steep but *uniform* gradient just drives faster
turnover (more drift, per §B). The innovation ladder still has to come
from *structure*: chemical niches, spatial niches, the egene/LUT
complexity axes. So "energy gradient → more evolvable" holds in the
weak sense (prevents gradient-collapse stagnation) but not the strong
sense (guarantees sustained innovation) — the strong sense needs the
gradient to be *structured*, which routes right back to §C and §E.
Adding energy-gradient realism is worth doing; its payoff is that it
unlocks and motivates those two, not that it is a standalone win. And
because Morowitz (1968) showed that energy flow through a system is what
organizes it into cycles, the highest-value first step is the cheap one
— *measure the flux you already have* (item 1) — before adding any new
mechanism.

## References added (see Docs/references/references.md)

- **Farmer, J. D., Kauffman, S. A., Packard, N. H. (1986).**
  *Autocatalytic replication of polymers.* Physica D 22:50-67.
- **Rasmussen, S., Knudsen, C., Feldberg, R., Hindsholm, M. (1990).**
  *The coreworld: emergence and evolution of cooperative structures in
  a computational chemistry.* Physica D 42:111-134.
- **Eldredge, N., Gould, S. J. (1972).** *Punctuated equilibria: an
  alternative to phyletic gradualism.* In Models in Paleobiology.
- **Lindgren, K. (1991).** *Evolutionary phenomena in simple dynamics.*
  In Artificial Life II, Addison-Wesley.
- **Rosenzweig, M. L. (1995).** *Species Diversity in Space and Time.*
  Cambridge University Press. (hump-shaped productivity–diversity.)
- **Fisher, R. A. (1930).** *The Genetical Theory of Natural Selection.*
  Oxford. (response to selection ∝ additive variance in fitness.)
- **Lotka, A. J. (1922).** *Contribution to the energetics of
  evolution.* PNAS 8:147-151. (Maximum-power / maximum-energy-flux
  principle — the §L thermodynamic reading of the inverted-U.)
- **Morowitz, H. J. (1968).** *Energy Flow in Biology.* Academic Press.
  (Energy flow through a system organizes it into cycles — the §L
  "measure the flux first" argument.)
- **Nicolis, G., Prigogine, I. (1977).** *Self-Organization in
  Nonequilibrium Systems.* Wiley. (Dissipative structures — §L.)
- (already in references: Bedau–Snyder–Packard 1998; Channon 2006;
  Schneider & Kay 1994; Chaisson 2001; Kauffman & Levin 1987.)

