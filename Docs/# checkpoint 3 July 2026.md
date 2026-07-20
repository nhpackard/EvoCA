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

# N Responses 19 July 2026

## Load-bearing claims:
I give responses to these because they are sometimes the most direct path to specifying an action.

### 1. 
- Action:  So let's do the split into two dicts.  
- Re. mu-genome branch, let's continue the plan to launch the per agent distinction in the development there.
- thx for tax distinction.  
    - Action: model.md needs to be updated; Phase 2c: Tax and Death mentions only tax, tax_per_egene, and tax_lut
### 2
- agree that upper resource regime tends toward neutral drift.  That is what I meant by "in the high resource regime, I'm thinking the fitness landscape is rather flat, with mild hills"  
- Action:  Scans must have more summary info.  I'm thinking a notebook with figures that would, e.g. clearly display peaks or transitions in the data, e.g. food_inc ~ 0.013.  E.g. scatter plot of run result as function of food_inc with other scan params encoded in marker (color, size, etc).
### 3
- Action:  add food_per.
- Use of env_mask can also use food_per; instead of 'F clamped to [0,1]', can have F restored to 1 at env_mask values every food_per (defaulting to current when food_per=1).
### 4
- Re. swarmchem:   True that swarmchem spec is not based on a reaction network directly.  But I'm not clear why it is a completely "negative reference".  If I create a swarm chem spec with many many components and random interactions between components, shouldn't I expect a reaction network to emerge?  Emergence of a reaction network from chemical bond priors is actually interesting, and to my knowledge unexplored in swarmchem literature (or correct me!).
- Action: place chemistry branch on the research map, with a plan for its own branch.
- Cf. 6 below, with comment on simple energy chemistry.  Should this be a branch off of chemistry, or a branch of its own?
### 5
- Action: plan scan using dyn_activity
- When you say '"Give each gene its own neutral shadow" is ill-posed in the open LUT-hash space (a just-born genome has no drift history to estimate an expectation from)', this seems a pretty weak objection.  What's wrong with simply eliminating the ill-posed case by defining it to be zero (in case of no drift history)?
### 6 and §L
- Yes: "the gradient is nowhere *measured* and is *unstructured*".  
- Action:  implement §C's source-field spatial gradient.  I like this, and have used it many times in the past.  Especially in Bugs context
- Good point about EvoCA agents being sessile, and this meaningfully constrains our evolutionary narrative.  Should we consider making a motile version???  (would be a new major branch, like chemistry).  What could possibly be the move algorithm?  1) swap places with neighbor that has least food.  2) displace neighbor with least food, leave dead cell behind.  3) implement an LUT brain as in Bugs.  4) other?
- Here is a model that would capture an explicit gradient dependence:
  - Expand food chemistry to *two chemical species*:  food and waste.
  - Food is gathered as currently
  - Waste is manufactured by the metabolization of food:  any current tax does not simply decrement food, but changes it to waste.
  - Waste must be discarded into the environment.  Roughly, I'm thinking about an inverse of the feeding algorithm; environmental food low=> easy to discard waste, environmental food high => hard to discard waste.
  - open question: have an analog to mouthful algorithm using a waste analog to egenes (i.e. wgenes)?  Maybe we could start with a simple no-wgene version.
  - Penalty for waste reaching max (1.0): agent dies.
  - the point:  gradient = food-waste comprises an explicit gradient 'measurement' that conditions agent evolution.
  - what do you think?  would this change evolution in a material way, or wind up giving pretty much same as current observations?

## Sections

Any comments and/or questions and/or actions left over after Load Bearing sections

### §A

- interesting to think of baseline tax as environmental.  I usually regard it as agent-based, with harsh environments (lower resource) evolving agents that have lower tax.  This framing requires, of course, evolvable per-agent tax.
- Re. distribution seeding:  agree that "Construct populations with a distribution of values for one param, see if any value dominates" requires per-cell storage of that param.  Plan to explore in the mu-genome branch.

### §B

- when you say "**the flat landscape is
  detectable only as raw-activity-high but excess-activity≈0.**", doesn't this suggest a new activity probe, raw-activity - excess-activity ?

### §C
 
- Anything wrong with simply having food_per param work on restoration of env_mask food source?  In fact, we could have env_mask active always, with default being all sites 1, shaped food being env_mask zero except specifically shaped food source.  env_mask taking values between 0, 1 (with restoration of food up to env_mask values every food_per steps) could be used to implement the gradient environments.

- Definitely want to test "does temporal lumpiness create evolutionary opportunity?".  
  - Action: Must go on the research somewhere explicitly.

### §D

- (1) separate punctuation / optimization or not:  
  - Not sure I agree that punctuation and optimization can be cleanly separated.  Your suggestion for punctuation has a 'punctuation event' hapening when a new slot activates.  I suspect that such a punctuation event would yield negligible effect.  Maybe a punctuation event could be defined as adding a slot filled with random egene values.
  - Anyway, in nature, there is no explicit separation.  Yet punctuation + optimization is observed.  So we should be able to observe it in a good model for evolution, without explicit separation.
- (2)  How do I see "per-lineage net metabolic return"?  Do we have a probe for that?
- (3)  Good point about frozen-substrate metric trap.  

### §E
- agree with predefined universe over stochastic on-the-fly.
- You say "meta-mutation-rate to tune", but actually there will be metaparams to define the "predefined universe".  
- agree with sparse-vector implementation.  But we may want to scan the resource vectors to populate additional resource planes for visualization.
- design gap:  yes, we must decide what to do with all these resources.  My first thought is to treat them like the current food resource, with the food plane becoming a vector valued plane.  We must have a simplified version of the eating algorithm.  Having egenomes for every resource is untenable. We need a mechanism for resource communication between agents.  This could be via the invention of a new agent interaction algorithm which would be added model complication, and a pain.  Maybe have the resource gathering (eating) algorithm look at Moore neighborhood.  Resource communication between agents could than be via this enlarged eating neighborhood.  Open to other suggestions.

### §F

- Action:  implement "genome-level analog of the transition-level idea: an activity that keys on the cell's LUT masked to the entries it actually exercises".  What should we call it? dyn_genome_activity ??

### §H

You say:  "Three of the sections above (§B per-lineage net return, §D punctuation fitness proxy, §C/§F "who actually out-reproduced") reduce to parent→child bookkeeping. That field now exists in code as an opt-in SoA addition (lin_parent_hash/birth_id/parent_id, zero-cost when off)."

What is "That field" that now exists?

You say "It remains the single highest-leverage instrument"

What is "it"?

Yes: let's use the lineage field. Please expand your comment in §H a bit, with explicit research program to distinguish between adaptive and drift.

### §I, §K

- Expand table to include all actions identified in responses above.
- Revise prioritized roadmap to include all actions identified in responses above.


### §L

- appreciate "Where EvoCA sits..." comments.  I mostly agree.  Must re-examine Lotka.
- Re. "what to actually build":  chem branch has a special priority because a potential collaborator might be interested.  
  - Action:  create the chem branch with spawned agents, perform the natural tests to retrieve current functionality with one resource field.

## Big picture
I have a couple of comments regarding the big picture of this research program:
- We need to be more careful to digest research as it is produced.  This should entail going beyond forming md files to report the results, and should include forming notebooks to summarize results, including both md cells and graphics to clearly illustrate the major discoveries.  (this is an elaboration of my comment on enhancing result reporting of scans, above)

- I think we need to consider what the aim of this endeavor should be, regarding communication with the world.  Are we aiming for a giant article?  a series of articles?  a book?  Or perhaps publication through structuring of the github repo, with notebooks reporting results?  In this last case, we would need to be a bit more rigorous about structure of the github repo.  Or consider a publication fork.

---

# Claude analysis (2026-07-19): response to N's replies

Framing: this is a deep-analysis pass responding to your `# N Responses 19 July 2026`
block. I re-verified every load-bearing claim against the *current* code and the
two local repos before answering (evidence inline; the file:line citations are
fresh this pass, not carried over). Where you pushed back I either concede with a
reason, refine, or hold the line with a stronger argument than the first pass gave
— never a polite split-the-difference. Where you raised a new dev path (motility,
waste chemistry, publication structure) I treat it as a first-class candidate with
its own verdict. The consolidated action table (§I′) and revised roadmap (§K′)
capture **every** action you named plus the ones this pass adds; nothing is dropped.

## Digest of this pass (read this if you read nothing else)

1. **The three pointed factual questions, answered.** "That field"/"it" in §H = the
   **lineage field** (`lin_parent_hash/lin_birth_id/lin_parent_id`, opt-in via
   `evoca_enable_lineage`, `C/evoca.c:215-232,1628-1685`). The per-lineage
   net-return probe you ask about in §D(2) **does not exist yet** — we have the
   *intake* side (`g_sum_mouthful`, `C/evoca.c:375`) and a per-egene mean-intake
   probe, but no lineage-keyed cumulative (intake−tax) reducer. It is a small,
   well-scoped build on top of the lineage field, and it is the single instrument
   most of this checkpoint depends on.

2. **Your "define the ill-posed case to zero" (LB5) is not a weak-objection patch —
   it silently changes the estimand.** In the open LUT-hash space, zero-filling the
   no-history genomes doesn't rescue a per-gene neutral threshold; it turns *excess
   activity* back into *raw activity* for exactly the churning majority of genomes
   the neutral model exists to discount. Full argument below — this is the one place
   I most firmly hold the line, with a sharper reason than the first pass.

3. **Your punctuation pushback (§D-1) is right and I over-engineered the first
   answer.** A good evolutionary model should show punctuation *emergently in one
   run*, not via two contrived setups. I concede the two-experiment framing and
   replace it with a single unforced run + a punctuation *detector*. And your
   mechanistic instinct is sharper than you gave it credit for: activating an egene
   slot is initially near-neutral-to-slightly-deleterious (it adds tax before it
   wins any max-score), so the *signature* is a dip-then-climb in lineage net-return
   — which is a cleaner, more literature-faithful punctuation marker than a positive
   jump.

4. **swarmchem (LB4): you're right that I overstated, but the correction is
   narrower than your rebuttal.** Random *pairwise* interactions among many
   components do produce emergent higher-order structure — that's real (particle
   life / primordial-particle systems). But stock swarmchem cannot host it: its
   coupling is *one-sided and type-blind* (an agent responds to the neighbor
   *aggregate* via its own genome only; there is no type×type matrix — verified in
   the code). And emergent self-**assembly** is not an emergent **reaction network**:
   the latter needs species-identity transitions (A+B→C), which no force-based swarm
   gives you without an explicit reaction rule. The literature for your actual
   question ("reaction network from bond priors") is *artificial chemistry*
   (Dittrich–Ziegler–Banzhaf; Fontana–Buss), not swarm chemistry. Net: swarmchem
   stays a negative reference *for the chem field*, but your instinct points at a
   genuinely under-explored question that belongs to a *motile particle* model —
   which connects straight to your motility question.

5. **The waste model (LB6) would change evolution materially — in the exact
   direction your own §B intuition wants.** Food→waste conversion with
   disposal-difficulty rising in crowded/rich regions installs a *density-dependent
   cost that stays selective precisely where the food axis goes slack*. That directly
   attacks the high-resource drift-collapse: abundance no longer means trivial
   survival, because waste management becomes the live selective axis. It is best
   built as the **first milestone of the chem branch** (N_chem = 2, food+waste) —
   which also answers your "branch off chem or its own?" question. One refinement:
   couple disposal to a *waste* field, not to the food field.

6. **Motility is not "another branch like chemistry" — it is potentially a change to
   what EvoCA *is*.** Chemistry enriches the resource while leaving the CA-on-a-fixed-
   lattice core intact. Motility breaks the core invariant (cells sit still and play
   a CA defined on their fixed neighbourhood). Your four move-algorithms split
   cleanly: options 1–2 are undirected/anti-foraging as literally stated and won't
   climb a gradient; only option 3 (the Bugs LUT-brain) gives *evolvable gradient
   navigation*, which is the whole point of §L's motile motivation. But adopting it
   turns EvoCA into "Bugs with CA-genomes." That may be exactly what you want — but
   it should be a conscious identity decision, not a branch we drift into.

---

## Load-bearing responses

### LB1 — env/agent split, mu-genome, tax placement

- **Two-dict split: confirmed, do now, on `main`.** Zero dynamical risk; derive both
  dicts from the single canonical recipe (do *not* fork storage — see §J′).
- **mu-genome per-agent distinction: continue on the branch, storage-first.** Agreed.
- **tax framing — you and I are not actually in conflict; the reconciliation is the
  interesting part.** I called baseline tax "environmental"; you call it agent-based
  (harsh/low-resource niches evolve *lower-tax* agents), which needs per-agent
  evolvable tax. Both are right about *different layers*, and the resolution has
  teeth: the environment sets a **cost-of-living floor** `tax_env` (a niche
  property — arctic = high floor), and the agent carries an **evolvable
  tax-modifier gene** above it, `tax_eff = tax_env + tax_gene`. Your harsh-niche
  intuition is then a *prediction of the model*: in a high-`tax_env` niche,
  selection on `tax_gene`↓ is strong. Crucially, §L says tax is the sink that makes
  food a *gradient*; if tax could evolve all the way to 0 the gradient collapses and
  selection dies — so the environmental **floor is thermodynamically load-bearing**,
  not just bookkeeping. Design consequence: when tax becomes evolvable in the
  mu-genome branch, keep a small env floor and evolve only the modifier.
- **Action (yours, confirmed):** `model.md` Phase 2c currently lists tax,
  tax_per_egene, tax_lut but the code's tax has a **fourth term** `tax_ring`
  (`C/evoca.c:1197-1200`, `tax = gtax + gtax_per_egene·Negene + gtax_lut·popcount +
  gtax_ring·ring-level`). Queued as a doc fix (A3).

### LB2 — scans need real summary figures

Endorse without reservation; this is the same instrument §B needs. The figure you
describe (run-outcome vs `food_inc`, other scan params encoded in marker
color/size/shape) is right, and I'd add three things so it does double duty as the
§B test: (a) put **`food_inc/tax`** (the §L gradient-steepness ratio) on the x-axis
as a second panel, since the peak should track the *ratio*, not absolute food; (b)
overplot *raw activity* and *excess/dyn_excess* on the same axis so the drift regime
shows up as the raw–excess divergence (LB/§B below); (c) mark the ≈0.013 window you
already found. This is a `main` analysis notebook, and it is the natural first
artifact of the "digest research as notebooks" big-picture point.

### LB3 — food_per + env_mask unification

Strong endorse, and your unification is cleaner than my first-pass "add a separate
`F_source` array." Ground truth: `env_mask` today is **`uint8` binary**, default all
1, and regen is *additive to a hard 1.0 cap every step* (`C/evoca.c:1180-1184`:
`if(!env_mask[i]) continue; F+=food_inc; if(F>1) F=1`). Your proposal — make the mask
**float-valued in [0,1]** and restore *up to the mask value* every `food_per` steps —
collapses three features (temporal lumpiness, spatial template, gradient niche) into
one knob, and reduces *exactly* to current behaviour at `food_per=1`, `mask≡1`. One
design fork worth deciding consciously:

- **(a) additive relaxation toward the float ceiling, fired every `food_per` steps**
  (`F ← min(F + food_inc·food_per, mask[i])`). Keeps `food_inc` as the *rate* and the
  mask as the *ceiling/template*; preserves the slow-relaxation mechanism that lets
  depletion craters persist (the Bugs behaviour — confirmed: Bugs does `F ←
  min(F+food_inc, F_src)` every step, `bugs.c:1417`). **Recommended.**
- **(b) hard set `F ← mask[i]` every `food_per` steps.** Sawtooth pulses; also
  interesting, but erases sub-`food_per` depletion structure at each pulse and drops
  `food_inc` as a knob.

Go with (a): it is the strictly more general one (b falls out as `food_inc→∞`), and
it's the variant whose *mechanism* (depletion-structure lifetime, §C) we actually
theorised. Concrete change: `env_mask` `uint8`→`float`, one modulo gate in Phase 2.
Answering your §C sub-question directly: **nothing is wrong with routing everything
through `env_mask`** — a single always-active float mask (default 1.0 everywhere,
shaped values for niches/gradients, restored every `food_per`) is the right unified
design, and it retires the idea of a third env-gating array.

### LB4 — swarmchem: I overstated "negative reference"; here is the precise correction

You pushed on "completely negative reference," and you're partly right, so let me be
exact about *what* is and isn't there, because the distinction is scientifically
load-bearing.

**What stock swarmchem cannot give you (verified in the code this pass).** Its
inter-agent coupling is **one-sided and type-blind**: an agent's acceleration depends
only on *its own* genome (`c1..c5`, `nbhd_radius`) applied to the *aggregate* mean of
its neighbours' positions/velocities — there is no type×type interaction, no
per-pair parameter, no way for "A near B" to differ from "A near C". So there is no
substrate in stock swarmchem for "random interactions *between* components": the
model has no place to *put* a pairwise interaction. That is why it's a negative
reference for a *reaction network*, and that claim stands.

**Where your intuition is right, and it's a real phenomenon.** If you *build* a
different model — many particle types with a **random type×type interaction matrix**
(attraction/repulsion depending on the specific pair) — you *do* get emergent
higher-order structure: stable clusters, membranes, self-maintaining and even
self-replicating spatial "organisms." This is the *particle life* / *primordial
particle systems* result (Schmickl, Stefanec & Crailsheim 2016, *Sci. Rep.*;
Ventrella's *Clusters*), and Sayama's own *heterogeneous* Swarm Chemistry (2009)
moves toward it. So "many components + random pairwise priors ⇒ emergent structure"
is correct and well-precedented.

**But — the distinction that decides the design.** Emergent self-*assembly* is not an
emergent *reaction network*. A reaction network requires **species-identity
transitions** (a particle/quantity of A *becomes* B; A+B→C with conservation and
catalysis). Force-based particle systems produce persistent *morphologies* you could
relabel as "molecules," but no particle ever changes *type* — there is no
stoichiometry, no catalysis, unless you add an explicit reaction rule on
contact. So random *bond-force* priors give emergent **shape**, not emergent
**chemistry**. Your specific question — "does a reaction network *emerge* from bond
priors?" — is genuinely under-explored, but its home is the **artificial chemistry**
literature (Dittrich, Ziegler & Banzhaf 2001, *Artificial Life* 7:225-275; Fontana &
Buss's AlChemy, 1994 — organisations/reaction networks emerging from a λ-calculus
"bond" prior), *not* swarm chemistry, and not the flocking substrate.

**Net for EvoCA.** For the **chem** branch (a resource *field* consumed by *sessile*
agents), the swarm substrate doesn't transfer at all and FKP-1986 remains the anchor.
Your reaction-network-from-priors question is real but lives in a **motile
particle** model — i.e. it's actually a question for the *motility* branch, not the
chem branch. If motility happens, revisiting "type×type priors ⇒ emergent
reaction/assembly" there would be novel and worth it. I've added the particle-life
and artificial-chemistry citations so the distinction is on the record.

### LB5 — "give each gene its own neutral shadow, define no-history = 0"

This is the one place I most firmly disagree, and I owe you a stronger reason than
"ill-posed," because as stated that *did* sound like a hand-wave.

**The zero-fill doesn't patch a boundary case — it changes what the statistic
measures, for the majority of genes.** Recall `excess = ΣG − ΣN`
(`evoca_explore.py:284`): the neutral term ΣN is what a *non-selective* process
would produce, and excess is the part above that. The whole point is that a neutral
genome *still accumulates activity* just by persisting and drifting. Now, in the
**open LUT-hash space**, almost every genome that matters is *young and rare* — the
churn you distrust in §F is precisely a stream of short-lived, freshly-minted hashes.
A per-gene neutral shadow would need *that specific hash* to also exist in the shadow
population with its own drift history; it generally does not (the shadow drifts to
*different* hashes). So "no drift history" is not the exceptional case — **it is the
generic case for the exact genes the metric is supposed to discipline.**

Set those to zero and `excess = ΣG − 0 = ΣG`: you have quietly turned *excess
activity* back into *raw activity* for the churning majority. The neutral correction
then does nothing for the population that most needs it, while still "working" for
the handful of long-lived genomes — which will *look* like a working neutral model
while being raw activity in disguise. That's worse than no per-gene threshold,
because it's a silent bias toward calling churn "adaptive."

**Why pooling is the *correct* estimator here, not a workaround.** The neutral
expectation for "a gene like this, arising now" is a property of the *ensemble*
birth/death/drift dynamics, not of the individual gene's (nonexistent) history. The
pooled/shared-threshold shadow estimates one neutral distribution from the whole
shadow ensemble and borrows statistical strength across genes — which is exactly the
information the zero-fill throws away. This is a bias–variance point: per-gene
open-space shadows trade a little bias for catastrophic variance (n=1 per gene), and
zero-filling the n=0 genes converts that into *systematic* bias.

**Where your instinct is exactly right and *should* be built: the fixed spaces.** In
the 729-key egene space and the 500-bucket dyn space, every bucket is present for the
whole run, has a well-defined Wright–Fisher drift trajectory, and a **closed-form
per-bucket variance** — so a per-bucket z-score `(obs − drift_mean)/drift_sd` is
well-posed and gives you exactly the "sharper per-gene threshold, at the cost of more
fluctuation" you asked for, with the fluctuation *bounded analytically* rather than
guessed. That's G2 in the roadmap and I recommend it. So the honest split is: your
per-gene instinct is *correct in the fixed spaces* (build it) and *ill-posed in the
open space* (the pooled shadow is right; zero-fill makes it worse, not fixed). And
the canonical middle rung — Bedau–Packard **A_new** (per-component new-activity vs a
*single shared* threshold, G1) — is the cheapest rigorous version of "threshold each
gene" and is turnover-invariant, which is why it also defuses the §B drift problem.

- **Action (yours):** dyn_activity scan — planned as D1 (see roadmap); pairs with
  promoting `dyn_excess` to the primary phenotypic readout.

### LB6 / §L — measure the flux, structured gradient, motility, and the waste model

Four distinct things here; I'll take the two you flagged as actions first, then the
two open questions.

- **Measure the flux (L1): confirmed, cheapest item, do first.** Intake side already
  exists (`g_sum_mouthful`); the sink side (Σ tax collected + Σ death losses + Σ
  cap-overflow) is a one-line accumulator per step. Per-capita net throughput is
  Chaisson's Φ and hands §B its principled x-axis. `main`.
- **Spatial source-field gradient (§C, L2): confirmed** — and it *is* the float-mask
  generalisation from LB3, so L2 and C2 are the **same change** viewed twice (env
  measurement + env structure). Build once.

**Open question — motility (your 4 options). This is the most consequential
architectural decision in your replies, and I want to slow you down on it.**
Chemistry enriches the *resource* and leaves the CA core untouched. Motility touches
the **core invariant**: EvoCA is "every cell carries a CA rule and the CA plays out
on a *fixed* lattice." If agents move, the CA-on-a-lattice picture is no longer
coherent — what is `v(x)`'s neighbourhood update when the occupants shuffle between
steps? Motility doesn't extend EvoCA; it converts it toward **"Bugs with
CA-genomes."** That might be the right destination, but name it as an identity choice.
On the four algorithms specifically:

| Option | Verdict |
|---|---|
| 1. swap with least-food neighbour | As stated this moves the agent *toward* the poorest cell — **anti-foraging**; it cannot climb a food gradient. Only makes sense as a mixing/diffusion rule, not adaptation. If the intent is foraging, it must be "move toward *highest*-food neighbour." |
| 2. displace least-food neighbour, leave dead cell behind | Non-conservative; leaves vacancy "trails" (interesting spatial structure) but again undirected unless the target is chosen by gradient. A churny diffusion coupled to reproduction filling vacancies. |
| 3. LUT brain as in Bugs | **The only evolvable, gradient-sensing option**, and the only one that makes §L's "navigate the gradient" meaningful. Verified Bugs mechanism: a **512-entry LUT** keyed on the 9-bit thresholded-food Moore neighbourhood → one of **120 moves** (8 dirs × 15 magnitudes), thresholds themselves evolvable (the egenome). This is a *second* genome alongside the CA-rule LUT — a large but principled addition. |
| 4. other | Two I'd add: **(4a)** a *non-evolvable greedy climb* (move to highest-food Moore neighbour) as a **baseline to test whether motility helps at all** before committing to the evolvable brain; **(4b)** a single evolvable "taxis gene" (probability of stepping up-gradient) — a graded rung between passive and full-brain. |

Recommendation: if you pursue motility, do it as a *separate major branch* (`motile`),
start with **4a** (does mobility even change outcomes on a gradient world?), and only
build the **option-3** Bugs-brain if 4a shows motility matters — because option 3 is
the change that redefines the model, and it should be earned by evidence, not
assumed. And note the payoff couples to LB4: a motile particle world is where your
"reaction-network-from-priors" question actually lives.

**Open question — the food+waste model. My verdict: yes, materially different, and
in the direction your own §B intuition wants — build it as chem's first milestone.**
Reasoning, not applause:

- The current tax is a scalar sink: `f -= tax`. Your waste model converts food→waste
  and makes *disposal* the hard part, with difficulty rising in rich/crowded regions.
  That installs a **density-dependent cost that stays selective exactly where the
  food axis goes slack.** In the current model, high resource ⇒ everyone survives ⇒
  drift (the §B high-flank collapse). With waste, high resource ⇒ waste piles up ⇒
  waste-management becomes the live selective axis ⇒ **selection does not
  collapse at abundance**, it *changes character*. This is precisely the "even in a
  resource-rich world there is still nontrivial evolution" exception you posited in
  the original checkpoint — the waste model is a mechanism that *guarantees* it.
- It is a genuine **structured internal gradient** (§L) without full vector
  chemistry: a 2-vector (food, waste) is the minimal chem, so it doubles as the
  validation milestone for the chem plumbing (N_chem=2 before N_chem=many). That
  answers "branch off chem or its own?": **make food+waste the first concrete
  deliverable *on* the chem branch.**
- **One refinement (push-back on the coupling).** You proposed disposal-difficulty
  keyed to the *local food* field ("inverse of feeding"). Cleaner: give waste its
  own environmental field `W_env`, dump waste into it, disposal rate ∝ `(1 − W_env)`
  (hard where the neighbourhood is *already waste-saturated*). Then (i) waste has its
  own spatial gradient and plumes, (ii) the crowding penalty *emerges* (crowds make
  local waste) instead of being imposed by conflating it with food, and (iii) it's a
  true second compound, which is what you want the chem machinery to exercise. Start
  without wgenes (uniform disposal) exactly as you suggest; add a wgene max-score
  disposal later by symmetry with eating.
- Honest uncertainty: whether this produces *sustained* new evolution or just a
  one-time shift to a waste-limited steady state is exactly the experiment — but even
  the null result ("waste just moves the drift-collapse to a waste-poisoning
  collapse") would be informative about whether a second conserved currency changes
  the selection structure.

## Section-level responses (leftovers)

**§A.** tax-as-agent reconciled above (LB1). Distribution-seeding on mu-genome:
agreed, it's the natural *first* use of per-cell storage (a zero-mutation selection
assay) and I'd run it before enabling rate-mutation.

**§B — your "raw − excess" probe.** Nice instinct, and let me sharpen it with a fact
that makes it *free*: since `excess = ΣG − ΣN`, we have **`raw − excess = ΣN` — the
neutral shadow's own activity, which we already log** (`N_total_activity`). So the
quantity you want isn't a new measurement; it's the shadow baseline, already in
hand. The *diagnostic* you're really reaching for is the **adaptive fraction**
`excess/raw` (equivalently `1 − ΣN/ΣG`): in the drift regime it →0 while raw stays
high; at the productive peak it's maximal. That single normalised scalar is the crisp
"is the activity adaptive or drift?" readout — a graded version of your idea. Trivial
to add to the §B notebook; no C change.

**§C.** Unified float-`env_mask` answered in LB3 (yes, route everything through it).
"Does temporal lumpiness create evolutionary opportunity?" is logged as an explicit
campaign (C3) on the research map, with the falsifiable prediction from the first
pass: mean `F_env` held fixed while spatial variance rises with `food_per`, and
turnover/excess rise with that variance (or the informative null: they don't).

**§D — punctuation vs optimization. You're right; I concede the framing and
strengthen the mechanism.** Two contrived setups was the wrong instinct — a good
model should show punctuation *emergently in a single run*, which is exactly how
Coreworld's seven epochs and Lindgren's strategy-epochs appear, and what
Bedau–Packard evolutionary-activity waves *are*. So: **one unforced run + a
punctuation detector**, not two experiments. And your mechanistic doubt is sharper
than you claimed: with the current multi-slot egene (verified — `Negene` active
ternary slots, mouthful = **max** score across active egenes,
`C/evoca.c:1052,1198`), activating a slot adds a pattern to the max-pool *and* adds
`tax_per_egene` — so a fresh random slot is **near-neutral-to-slightly-deleterious at
birth** (it costs tax before it wins any max), and only becomes beneficial once
optimized to out-score the incumbents. That is *textbook* punctuated equilibrium
(innovations near-neutral at origin, adaptive after refinement), and it gives a
concrete detectable **signature: a slot-activation event → a dip/plateau in
lineage net-return → a climb as the slot is tuned.** Your "slot filled with random
values" *is* what activation already does, and it confirms the dip. So the
experiment is: single run, lineage net-return time series, detect the
dip-then-staircase around `Negene` increments. This is better science than my
first-pass separation — thank you for the push.

**§D(2) — "how do I see per-lineage net metabolic return? do we have a probe?"** No.
Verified: we have `g_sum_mouthful` (global intake) and the per-egene mean-intake
probe, but **no lineage-keyed net-return reducer**. It needs building (D2), and it's
small given the lineage field already exists: at each death/step, attribute
`(intake − tax)` to the cell's `lin_birth_id`, reduce by lineage root (walk
`parent_id` chains). This is the instrument §D, §B, and §H all depend on — highest-
leverage small build in the document.

**§E.** Conceded: (i) "no meta-mutation-rate" was imprecise — a predefined universe
still has universe-defining *metaparams* (compound count, reaction density, rate
distribution); the point that survives is *no per-run stochastic species allocation*.
(ii) Sparse vectors as storage, **plus** optional dense per-compound planes rendered
*on demand for visualization* — agreed, that's a viz projection, not the hot-path
representation. (iii) Your eating design — vector food plane, simplified eating,
resource communication via an **enlarged (Moore) eating neighbourhood** rather than a
new agent-interaction protocol — is the right instinct (reuse the neighbourhood you
already scan; don't invent a messaging layer). The one thing still to pin *before
coding* remains the mapping from the resource vector to the scalar `f` that gates
reproduction, and whether an egene targets a *compound profile* or still a *spatial
pattern*. The food+waste milestone (LB6) is the ideal place to settle it concretely
with N_chem=2.

**§F — naming the used-mask genome activity.** `dyn_genome_activity` is
serviceable but overloads "dyn" (which currently means the 500 *transition* buckets);
a reader will expect a fixed 500-space. Since the new probe is a *genome-level
species* count where the species key is the LUT **masked to entries the cell actually
exercises**, I'd name it for what it keys on: **`eff_activity`** ("effective-rule
activity") or **`used_lut_activity`**. My pick: **`eff_activity`**, with the docstring
"full-genome activity over the *effective rule* (LUT restricted to exercised
entries)," and keep `dyn_activity` for the transition-level probe. Reuses the existing
per-step `lut_active` mask; the only new part is *per-cell* used-masking.

**§H — expanded, with an explicit adaptive-vs-drift research program.** "That field"
= the lineage field (`lin_parent_hash/lin_birth_id/lin_parent_id`); "it" = that same
field, which I called the highest-leverage instrument because parent→child
bookkeeping is what turns correlational activity into causal selection claims. The
program to separate adaptation from drift, concretely:

1. **Instrument:** enable lineage; build the D2 net-return reducer; log per-lineage
   `(intake, tax, births, extinction step)`.
2. **Neutral control by construction:** the drift shadow already gives the null. The
   test is per-lineage: does a lineage's realised growth/net-return exceed what its
   *shadow* (same birth/death counts, non-selective inheritance) produces? Lineages
   that beat their shadow are under selection; those tracking it are drifting.
3. **Two orthogonal readouts, cross-checked:** (a) **A_new** (per-component
   new-activity vs shared threshold, G1) — turnover-invariant, so it stays honest in
   the high-resource drift regime where raw activity lies; (b) **lineage net-return
   staircases** — a lineage whose net-return steps *up* and *persists* is adaptive; a
   lineage whose hash churns while net-return is flat is drift. Agreement between the
   population-level A_new and the lineage-level staircase is the strong evidence.
4. **The decisive sweep:** run the §B `food_inc/tax` sweep with all of the above.
   Prediction: adaptive-fraction (`excess/raw`) and count of shadow-beating lineages
   are **inverted-U** in the ratio; raw activity is monotone/flat. If raw is
   inverted-U too, my drift-collapse story is wrong — that's the falsifier.

This is the backbone that makes §B, §D, and §L *measurable* rather than narrated.

**§L.** Lotka re-examination: the piece that actually delivers your interior optimum
is **Odum & Pinkerton (1955)** — *maximum power at intermediate efficiency* — which
is the sharper form of Lotka's principle and the one that predicts the inverted-U as
a selected optimum in gradient-exploitation rate rather than a monotone "steeper is
better." I've added it; it's the citation to lean on for the `food_inc/tax` peak.
chem-branch priority (collaborator): acknowledged — the branch-creation +
single-field-parity task (L3) is in the roadmap as a near-term branch item, and the
food+waste milestone gives the collaborator something concrete and novel early.

## Big picture

- **Digest-as-notebooks: endorse as standing process, not a one-off.** Concretely:
  every campaign closes with a `Research/<campaign>/` subdir co-locating its summary
  notebook (md narrative + the figures that show the discovery) with any `.md`
  write-up — a result and its narrative never live in separate directories. The §B/LB2
  scan-summary notebook is the first instance and the template. I'd make "campaign
  isn't done until its digest notebook exists" a workflow rule (it pairs with the
  existing pre-merge test-notebook protocol).
- **Publication aim — my recommendation, with a reason.** Don't pick "one big
  article" now. The structure of the work (a core model + semi-independent
  branches: resource-driving, chem, mu-genome, possibly motile) maps naturally onto a
  **repo-as-publication with a small series of focused papers/notebooks**, not a
  monograph. Concretely: (i) make the **GitHub repo the primary artifact** —
  `Research/` as the results layer (per-campaign subdirs, notebook + narrative
  co-located), `Docs/` as the spec/theory layer, a top-level
  README that reads as a paper index; (ii) each mature branch yields *one* focused
  methods+results notebook that can become a short paper (the neutral-model/activity
  methodology is itself a publishable methods contribution — it's the sharpest,
  most-defensible piece); (iii) reserve the "book/large synthesis" for after two or
  three of those exist and the throughline (dissipative structure + heredity +
  open-ended activity metrics) is demonstrated, not asserted. A **publication fork**
  is premature; a *disciplined main repo* gets you the same shareability without
  maintaining two histories. This also raises the bar on repo structure now, which
  the digest-notebook rule already pushes toward. (Flag: this is the one area where I
  have the least project-specific ground truth — it's a judgement call about your
  audience and collaborators, so treat it as a recommendation to react to, not a
  finding.)

## §I′ — consolidated action table (every action, old + new)

| # | Action | Source | Where | Risk |
|---|---|---|---|---|
| A1 | env/agent two-dict split (derive from single recipe) | LB1 | **main** | none |
| A2 | mu-genome: per-cell agent-gene storage (mutation off), then distribution-seeding assay, then evolvable rates (env-floor + agent-modifier for tax) | LB1,§A | **branch** mu-genome | high |
| A3 | `model.md` Phase 2c: add `tax_ring` term | LB1 | **main** (doc) | none |
| B1 | scan-summary notebook: outcome vs `food_inc` *and* `food_inc/tax`, markers encode scan params, overlay raw vs excess, mark 0.013 | LB2,§B | **main** | none |
| B2 | adaptive-fraction readout `excess/raw` (raw−excess = ΣN already logged) | §B | **main** | none |
| C1 | add `food_per` (temporal regen period) | LB3 | **branch** resource-driving | med |
| C2/L2 | `env_mask` uint8→float ceiling; additive relaxation toward mask every `food_per`; spatial gradient niches | LB3,§C,§L | **branch** resource-driving | med |
| C3 | log "temporal lumpiness ⇒ evolutionary opportunity?" as explicit campaign w/ falsifiable prediction | §C | campaign | — |
| D1 | dyn_activity scan; promote `dyn_excess` to primary phenotypic signal | LB5,§F | **main** | none |
| D2 | **build** per-lineage net-return reducer on the lineage field | §D(2) | **main** (small C+py) | low |
| D3 | single-run punctuation detector: lineage net-return dip-then-staircase at `Negene` increments | §D | **main** | none |
| E1 | chem branch on research map + plan | LB4,§E | doc + **branch** chem | — |
| E2 | food+waste as chem **milestone 2** (first scientific instance, N_chem=2; milestone 1 is single-food parity); needs minimal reaction-execution engine; waste has own `W_env` field, disposal ∝ (1−W_env) | LB6,§E | **branch** chem | high |
| E3 | chem design: sparse vectors (hot path) + optional dense planes (viz); enlarged Moore eating neighbourhood for resource communication; pin vector→`f` mapping first | §E | **branch** chem | high |
| F1 | `eff_activity` (effective-rule / used-mask genome species) | §F | **main** (C add) | low |
| G1 | Bedau–Packard **A_new** (per-component new-activity, shared threshold) | LB5 | **main** | low |
| G2 | per-bucket z-scoring in fixed eg(729)/dyn(500) spaces (your per-gene instinct, well-posed) | LB5 | **main** | low |
| G3 | *decision:* do **not** build per-gene open-space shadow w/ zero-fill (changes estimand) | LB5 | — | — |
| H1 | adaptive-vs-drift program (§H above): lineage + A_new + shadow-beating test on §B sweep | §H | **main** | low |
| L1 | dissipation accumulator (Σ tax + death + overflow); per-capita throughput Φ; §B x-axis | §L | **main** | low |
| L3 | create `chem` branch w/ spawned agents; parity tests for single-field functionality | §L | **branch** chem | med |
| M1 | **DEFERRED 2026-07-19.** Motility on hold: characterize sessile EvoCA (baseline) first, so motility's contribution is attributable. Open sub-question: start motility from *EvoCA* or from *Bugs*? (see note below) | LB6 | **branch** motile (not now) | very high |
| P1 | digest-as-notebooks standing rule; `Research/<campaign>/` (notebook + md co-located) per campaign | BigPic | process | — |
| P2 | publication: repo-as-primary-artifact + focused notebook-papers; defer book; no fork yet | BigPic | strategic | — |

## §K′ — revised prioritized roadmap

**Do now (main, low/zero code, high information):**

1. A1 env/agent split · A3 model.md tax_ring fix — trivial, unblock the niche framing.
2. L1 flux accumulator + B1/B2 scan-summary notebook with the `food_inc/tax` axis and
   raw-vs-excess overlay — this one figure tests the central intuition *and* validates
   the metric (still the single highest-leverage cheap item).
3. G1 A_new — cheapest rigorous, turnover-invariant metric; defuses drift regime.
4. D1 promote dyn_activity/dyn_excess to primary phenotypic readout (reinterpretation).

**Do next (mostly main, scoped):**

5. D2 lineage net-return reducer — the instrument §B/§D/§H all need.
6. H1 adaptive-vs-drift program on the §B sweep (needs D2 + G1).
7. D3 punctuation detector (single unforced run) · G2 fixed-space per-bucket z-scoring
   (your per-gene instinct, made well-posed) · F1 `eff_activity`.

**Do on branches (invasive, human merge gate, pre-merge notebook):**

8. resource-driving: C1 `food_per` → C2/L2 float-`env_mask` gradient. Campaign C3.
9. mu-genome: A2 storage-first → distribution-seeding assay → evolvable tax
   (env-floor + agent-modifier).
10. chem (L3 branch + collaborator priority) — **general N_chem machinery first**
    (E3 design: sparse storage + vector→`f` mapping pinned up front), then milestones
    on it: **(1)** single-food parity → **(2)** food/waste binary (E2; needs the
    minimal reaction-execution engine + §L structured-gradient experiment) →
    **(3)** larger resource spaces + autocatalytic reaction-network design (FKP-style)
    and the egene→compound-profile redesign. Ordering per N 2026-07-19 (§E′ note).

**Deliberate decisions before any code:**

- **M1 motility — DEFERRED 2026-07-19** (N's call): explore sessile EvoCA first to
  establish what motility would *add*; and the prior question — EvoCA-as-trunk vs
  Bugs-as-trunk — is itself unsettled (see note). Revisit after the baseline
  campaigns (§B/§D/§H) are in hand.

  *Note — EvoCA-start vs Bugs-start (framing for when we revisit).* The two models are
  near-duals: **EvoCA** = sessile, rich genome (CA-rule LUT + eating egenome), resource
  *field*; **Bugs** = motile, thin genome (a food-sensing LUT-brain), no CA. Motility
  machinery (LUT-brain, `place_or_bump` occupancy, sensing) *already exists in Bugs* —
  so "Bugs + EvoCA's CA-rule/egenome richness" may be *less* work than "EvoCA + motility
  + resolving CA coherence." But the CA dynamics are EvoCA's distinctive scientific
  content, which Bugs lacks entirely. The real question is which distinctive feature is
  the trunk — the evolvable CA rule (EvoCA) or motile sensing (Bugs) — and the two are in
  genuine tension (a CA is defined on a *fixed* neighbourhood; if occupants move, "what is
  `v(x)`'s update?" must be answered either way). Deciding this needs the sessile baseline
  first, which is exactly why the deferral is correct.
- **E2 vs E3 ordering — RESOLVED 2026-07-19** (N's call): build the **general N_chem
  universe machinery first** (the trunk), then instantiate milestones small→large on
  it. Agreed milestone sequence:
  1. **Parity** — retrieve current single-chem (food) results at N_chem=1 (= L3
     regression gate).
  2. **food/waste binary** (N_chem=2) — validates sparse storage + vector→`f` mapping;
     runs the §B high-flank / structured-gradient experiment.
  3. **Larger resource spaces** — where the **autocatalytic reaction-network** design
     (FKP-style) lands; deferred to here by N.
  Distinction to hold in the design: food→waste *is a reaction*, so a minimal
  **reaction-execution engine** (transform A→B at a rate) is needed by milestone 2,
  while the **reaction-network specification** (autocatalytic universe) waits for
  milestone 3. The egene→compound-profile redesign also defers to 3 (milestones 1–2
  have ≤2 compounds and reuse current eating on the food component). Pin **sparse
  storage** and the **vector→`f` mapping** in the general design before milestone 2.
- **P2 publication shape — RESOLVED 2026-07-19, structure refined 2026-07-19**
  (N's call): **repo-as-artifact.** Layout, as settled:
  - **`Research/` is the results layer**, organized as **per-campaign subdirs**
    (`Research/<campaign>/`). Each subdir **co-locates** its summary notebook with
    its `.md` write-up (and a short per-campaign README), so a result and its
    narrative never live in separate directories. This deliberately avoids the
    disliked pattern of a results `.md` in one place linking to a notebook in
    another (`../Notebooks`). A top-level **`Research/README.md`** is the index /
    reading-order that summarizes the research endeavor.
  - **Raw artifacts stay where they land now** (`Runs/`, `Scans/`, `ProbeLogs/`).
    `Research/` holds *curated* results and may **copy** selected raw artifacts in;
    campaign notebooks read raw data from those dirs in place. The distinction that
    matters: **md → notebook-in-another-dir is what we avoid** (solved by
    co-location); **notebook → raw data** is just code reading its inputs and is
    fine.
  - **Notebooks keep their outputs** so GitHub renders figures inline (no run
    needed by a visitor). Note (verified 2026-07-19): the repo does **not** strip
    notebook outputs today — no `nbstripout` filter, `.gitattributes`, pre-commit,
    or hook — and existing notebooks already commit *with* outputs, so `Research/`
    needs **no special git treatment now**. If repo-wide output-stripping is added
    later (to tame noisy *working* notebooks), exempt `Research/` with a nested
    `Research/.gitattributes` (`*.ipynb filter=` / `*.ipynb diff=`). For a
    publication repo the robust posture is to *not* rely on filters for `Research/`
    and just commit outputs directly (the current default).
  - **`Docs/` stays the hodge-podge for now** (spec/theory + planning + logs mixed);
    the code-doc vs planning-doc split is a **publish-time** cleanup, not now.
  - No publication fork, no book for now; top-level README reads as a paper-index.
  This *raises the bar on repo structure now* — the P1 rule ("a campaign isn't done
  until its `Research/<campaign>/` digest notebook exists") is the enforcement
  mechanism, and the §B/LB2 scan-summary notebook is the first artifact built to
  that standard.

## §K″ — parallel execution plan (what runs concurrently)

The §K′ lists are **priority-ordered, not sequential.** Almost everything
in "Do now" and "Do next" is mutually independent and should run as
**concurrent agents**, not one-after-another. The only hard serialization
points are the human merge gate and a handful of explicit data
dependencies (called out below). This section makes the parallelism
manifest so agents can be spawned in waves.

**Isolation rule (why worktrees).** Concurrent agents that edit
`C/evoca.c` / rebuild `libevoca.dylib` on a single working tree collide.
So **every code change runs in its own `git worktree`+branch** (the proven
S2a–e batch pattern from the 2026-05 integration), and only analysis /
notebook / doc work runs directly on `main`. Merge stays a human gate; an
orchestrator digest to `Docs/research_board.md` collects each agent's
result + any decision it needs, so review is one digest, not N
transcripts. Standing preference (N): **spawn parallel agents whenever the
work is independent** — the default is concurrency, serial is the
exception that must be justified by a dependency or the merge gate.

### Wave 1 — launch now, all concurrent (no cross-dependencies)

| Agent / branch | §I′ actions | Isolation | Kind |
|---|---|---|---|
| `env-agent-split` | A1 + A3 | worktree | py + doc, trivial |
| `metric-flux` | L1 | worktree | small additive C+py |
| `metric-anew` | G1 | worktree | small additive C+py |
| `metric-zscore` | G2 | worktree | small additive C+py |
| `eff-activity` | F1 | worktree | small additive C+py |
| `lineage-return` | D2 | worktree | small additive C+py |
| `resource-driving` | C1 (→ C2/L2) | worktree (semi-permanent) | feature branch |
| `mu-genome` | A2 storage-first | worktree (semi-permanent) | feature branch |
| `chem` | E1/E3 machinery + milestone-1 parity gate | worktree (semi-permanent) | feature branch |
| `dyn-promote` | D1 | **main** | analysis / reinterpretation |
| `research-migrate` | P1 first digest notebooks (see candidates below) | **main** | notebook / doc |

That is **11 tracks in parallel.** The five `metric-*` worktrees are
exactly the shape of the S2a–e batch that already rebased and merged with
zero cross-branch conflicts, so this is a known-good fan-out, not a
gamble. Each Wave-1 code branch ships with its pre-merge test notebook per
protocol.

### Wave 2 — gated on Wave-1 merges (dependencies explicit)

| Agent / action | Waits on | Why |
|---|---|---|
| `H1` adaptive-vs-drift program on the §B sweep | `lineage-return` (D2) + `metric-anew` (G1) | the shadow-beating test needs both instruments merged |
| `B1/B2` final scan-summary figures | `metric-flux` (L1) + `metric-anew` (G1) | the `food_inc/tax` axis and raw-vs-excess overlay need those metrics — but scaffold the notebook in Wave 1 and drop figures in when they land |
| `D3` punctuation detector | `lineage-return` (D2) | needs per-lineage net-return |
| `C3` resource-driving campaign | `resource-driving` branch validated | runs on the new dynamics |
| chem **milestone 2** (food/waste, E2) | `chem` machinery (E3) + minimal reaction engine | food→waste is a reaction |

### Deliberately serial (not parallelizable)

- The **human merge gate** between waves (branches never auto-merge;
  dynamics-altering `resource-driving`/`mu-genome`/`chem` especially).
- The **chem milestone ordering** E-1→2→3 (parity → food/waste → larger
  spaces + autocatalytic network) — each milestone gates the next.
- Anything on the **`motile`** branch — deferred (M1).

### First Research/ entries — candidates (2026-07-20)

`research-migrate` (Wave 1) seeds `Research/` from completed campaign
work. Candidates below, strongest first; the WIP index table lives in
`Research/README.md` (Campaign + Question filled now; Headline + Notebook
fill as each digest notebook is built). Ordering favors the most
publication-defensible pieces.

1. **`neutral-model-methodology`** — *Can excess-activity metrics separate
   adaptive evolution from drift/turnover?* The flagship methods
   contribution: fixed-space eg(729)/dyn(500) shadows, reciprocal-freeze
   causal controls, and the two shadow-scope bugs those controls caught.
   (The checkpoint flags this as the sharpest, most-defensible piece.)
2. **`resource-driving-inverted-U`** — *Does evolvability peak at
   intermediate resource, and is the true axis `food_inc/tax` (gradient
   steepness)?* The §B/LB2 flagship; **science in progress** (the template
   digest notebook), so Headline stays open.
3. **`pure-evo-regime`** (campaigns #1/#2) — *What metaparameter regime
   maximizes open-ended activity under joint vs LUT-only evolution?*
   Headline: whole-genome excess ≈ 0 everywhere; dyn/eg excess strongly
   positive; low `mu_lut` + high `m_scale` wins; the `mu_lut` optimum is
   conditional on egene co-evolution.
4. **`coevolution-substrate`** (#3/#3c) — *Does egene co-evolution amplify
   rule selection; is a freezing penalty real?* Headline: egene-freeze
   reduces dyn-excess; #3c mostly GoL-substrate death confound + the
   `dyn_excess_pc` frozen-rich artifact (methodological result).
5. **`RD-robustness`** (#4) — *Is reaction-diffusion spatial structure
   robust under pure-evolutionary optimization?* Headline: no — corrL
   washes out; the static-colony cases are the lowest-`dyn_excess` ones
   (corroborates the suspect-the-metric intuition).
6. **`viability-brackets`** (R1) — *Where are the productive optima in
   (`m_scale`, `food_inc`)?* Headline: `m_scale` interior ≈2.5–3.5,
   `food_inc` high 0.013–0.018, U-shaped viability.

## Self-critique (adversarial pass — flaws caught and fixed before posting)

1. *Motility, first draft understated the stakes.* I originally filed motility as
   "another major branch like chemistry." An adversary rightly objects: chemistry
   leaves the CA-on-fixed-lattice core intact; motility *breaks* it (what is the CA
   neighbourhood update when occupants move?). Revised M1 to a *model-identity*
   decision with a cheap 4a gate before the identity-changing option-3, rather than a
   co-equal branch. (Changed a recommendation.)
2. *Waste model, first draft said "materially different" without saying which
   direction.* Vague endorsement is applause, not analysis. Revised to the specific
   claim that waste installs a **density-dependent cost that stays selective at
   abundance** — i.e. it targets the *exact* §B high-flank drift-collapse — and
   flipped the coupling from food-field to a dedicated `W_env` field so the crowding
   penalty *emerges*. (Changed an assumption and a design axis.)
3. *"raw − excess" probe, first draft treated it as a new probe to build.* On
   checking the formula (`excess = ΣG − ΣN`), raw − excess = ΣN, which we *already
   log*. So the honest answer is "it's the shadow baseline you already have," and the
   *useful* new scalar is the adaptive fraction `excess/raw`. Corrected from "build
   it" to "you have it; here's the sharper graded version." (Removed a spurious
   action, added a better one — the binary→graded move.)
4. *Punctuation, first draft defended my two-experiment separation.* Under
   adversarial re-read your single-run emergent framing is simply more correct
   (Coreworld/Lindgren/Bedau all emergent), so I conceded and *replaced* the
   recommendation with a detector, then found the *stronger* mechanistic point
   (slot-activation is near-neutral-at-birth ⇒ dip-then-staircase) that my separated
   version would have hidden. (Reversed a recommendation.)

## §J′ — structural & scaling review (zoom-out)

- **Single source of truth for metaparams (the §A split's one real risk).** Derive
  both `env_metaparams`/`agent_metaparams` dicts *from* the one canonical `.evoca`
  recipe; do not fork storage. This is the only place A1 can introduce a duplication
  bug.
- **`env_mask` is now doing three jobs (regen gate, spatial template, gradient).** That
  consolidation (LB3/C2/L2) is *good* — it retires two would-be arrays — but it means
  the float-mask semantics must be documented in one place (model.md) as the single
  env-structure primitive, or three campaigns will each reinvent its meaning.
- **Scaling: the dominant risk remains chem's per-step cost, not memory.** Sparse
  vectors (E3) are the precondition for N=512, unchanged from the first pass. New
  entry: **motility (M1) has a different scaling risk** — random-order sequential
  movement with occupancy conflicts (Bugs' `place_or_bump` does up to 64 probes per
  bug) is inherently serial and cache-hostile; at N=512 this could dominate the step
  loop far more than the CA core. Prototype cost at N=256 before committing.
- **Reuse, not parallel machinery:** F1 reuses `lut_active`; C2 generalises
  `env_mask` (no third array); D2/D3/H1 reuse the lineage field; the waste field
  reuses the food-field regen path with a sign flip. The one place to *resist* reuse
  is naming: `eff_activity` must not be folded into `dyn_activity` (different key
  space — see §F).

## References added (to Docs/references/references.md)

- **Odum, H. T., Pinkerton, R. C. (1955).** *Time's speed regulator: the optimum
  efficiency for maximum power output in physical and biological systems.* American
  Scientist 43:331-343. (Maximum power at *intermediate* efficiency — the interior
  optimum behind §B/§L's inverted-U; the sharper form of Lotka you asked to
  re-examine.)
- **Sayama, H. (2009).** *Swarm Chemistry.* Artificial Life 15(1):105-114.
  (Kinematic recipe model; heterogeneous variant gestures at type-specific coupling —
  LB4.)
- **Schmickl, T., Stefanec, M., Crailsheim, K. (2016).** *How a life-like system
  emerges from a simple particle motion law.* Scientific Reports 6:37969.
  (Primordial particle systems — emergent self-replicating structure from random
  pairwise motion priors; the "particle life" precedent — LB4.)
- **Dittrich, P., Ziegler, J., Banzhaf, W. (2001).** *Artificial chemistries — a
  review.* Artificial Life 7(3):225-275. (The correct literature home for
  "reaction network from bond priors" — LB4.)
- **Fontana, W., Buss, L. W. (1994).** *The arrival of the fittest: toward a theory of
  biological organization.* Bulletin of Mathematical Biology 56:1-64. (AlChemy —
  reaction networks/organizations emerging from a λ-calculus interaction prior —
  LB4.)
- **Packard, N. H. (2019).** *Intrinsic adaptation in a simple model for evolution.*
  Artificial Life 25(1). (The Bugs model — LUT-brain motility, the option-3
  reference for M1.)

