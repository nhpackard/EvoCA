# EvoCA research: stability analysis

This is an alternative approach aiming to replace the "patch transfer" thread (listed as research campaign 7 in research_board.md).

# the basic idea

Instead of extracting and inserting patches, let's think of dividing the world (lattice) into two vertical bands, with different populations in each band.  The key requirement for the setup of the stability measurement is that the two sub-populations in the two bands have some identifiable characteristic that distinguishes a cell being in one or the other.  Call left sub-population $A$ and right sub-population $B$.

Then we can define a dominance measure, per row:
$$
D_i = N_A - N_B
$$
and over all rows, we can then compute 
$\mu_D = $ mean over rows and $\sigma_C = $ the standandard deviation over rows.

Further, we can observe how $\mu_D$ and $\sigma_D$ vary over time, denoted $\mu_D^t$ and $\sigma_D^t$, where $t=0$ is the preparation of the initial population in bands covering half the lattice size, so $\mu_D^0 = 0$  and $\sigma_D^0 = 0$

If $\mu_D^t$ goes positive, $A$ dominates.  If $\mu_D^t$ goes negative, $B$ dominates.

# different preparations

We now discuss other possible definitions of $A$ and $B$.

## different rings

One natural wat to define $A$ and $B$ is by the complexity of the local rules implemented, i.e. whether the local rules depend on ring 1, 2, or 3.   So we could test whether more complex rules tend to dominate less complex rules, e.g. 

If $A$ has ring $3$ rules and $B$ has ring $2$ rules, and $\mu_D^t$ is seen to increase with time, we would conclude that ring $3$ rules dominate ring $2$ rules.

Naturally such a result is contingent on the settings of all other metaparams.  For the particular case that all mutation is set to zero, we obtain a basic result for competition between rules of different complexity, but these are still contingent on other metaparams, e.g. having to do with food and reproduction.  We will want to look at superposition of $\mu_D^t$ graphs as other metaparams (e.g. food) are varied.

We will want to compare stability between different initial preparations

- ring 1 vs ring 2
- ring 2 vs ring 3
- ring 1 vs ring 3

## inclusion of evolution

We currently do not have the ability of restricting mutation to act within a ring (or to be ring bounded).  We have discussed the possibility of making this enhancement for possible future work flows.

For the present, however, we may still include the possibility of evolutionary dynamics in $A$ / $B$ preparation.

### evolved vs. wild type

Consider a generic evolutionary experiment, e.g. for metaparams, use top winner for one of the recent scans, and initialize to random LUT.  Evolve for $N$ ticks.

Take $A$ to be half of the evolved population, $B$ to be half of the initial population.

One subtlety:  following $A$ and $B$ after the definition above will require a bit of book-keeping:  when reproduction happens reproduction from $A$ will need to be labeled with the $A$ tag, and likewise for $B$.  

Then we can address the question:  *does the evolved population dominate the unevolved initial population?*

### gradual vs rapid evolution

Another easy evolutionary experiment:

$A$:  half a lattice of an evolved population from initial population of 'simple' ring 1 LUT's

$B$:  half a lattice of an evolved population from initial population of random LUTs (mix of ring 1, ring 2, ring 3)

To address the question:  *is gradual evolutionary learning from simple initial condition more effective than evolutionary learning from a complex initial condition?"*

---

# Claude — agent critique, refinements, and additional A/B preparations (2026-05-19)

*Deep-analysis pass appended at user request. Integrates R1 and #3c
results + the #8 metric-degeneracy lesson into the framework. Style:
the framework above is sound and worth pursuing — what follows is a
set of substantive refinements (each changes a metric, assumption, or
experimental step, not just wording) and additional preparations to
broaden the design surface. Order: where the framework is strong →
flaws and revisions → additional A/B preparations → adversarial pass →
operational starting plan → open questions.*

## (from chat)

Headline shape of what I added

Section: A. Strengths
What it does: Brief — what the framework already gets right (tag separation, generality, spatial
dominance, replaces #7 cleanly)
────────────────────────────────────────

Section: B. Geometry
What it does: Surfaces what the doc leaves implicit: vertical bands → two interfaces (periodicity);

tabulates 6 alternative geometries (band, stripes, mosaic, single-cell invasion, quadrants) —
many
    questions are unanswerable in the default geometry
────────────────────────────────────────

Section: C. Ten refinements (each changes a metric / assumption / step)
What it does: C.1 add frequency-based F_i alongside D_i (raw count conflates competition with
total-pop change — exactly R1's regime-comparison problem); C.2 σ_D is weak → add front position,

roughness (KPZ-testable), contact fraction; C.3 A/B tag is free for  mu_lut=0 via get_lut_hashes,

needs a small additive a_b_tag C addition for evolutionary preps (queue alongside the
ring-bounded
    operator); C.4 reframe "evolved vs wild-type" as seeded invasion of descendants; C.5 viability
screen is a  prerequisite — without it ring-class dominance reads are dominated by which seed
survives at all (#3c's lesson); C.6 uniform-rule vs class-mosaic preps answer different questions

(Simpson's-paradox risk); C.7 pin symmetric ICs outside the rule axis; C.8 5-way outcome
classifier (extinct-A, extinct-B, A-dom, B-dom, co-existence);  C.9 early-vs-late window check
(#3c's plateau lesson); C.10 null control is mandatory  (#8's lesson)
────────────────────────────────────────

Section: D. Eight additional A/B preparations
What it does: Feasible now (no C): density, initial-state-structure, egenome-class,
co-evolved-pair,
    diversity, invasion variant of all the above; Needs small additive C: hostile-environment (may
already be feasible via env_mask — flagged to verify), per-region mutation rate
────────────────────────────────────────

Section: E. Adversarial peer pass
What it does: Three rival objections and the revisions they force — promote front metrics to
primary; mandate Screen-0 viability; mandate null spread baseline
────────────────────────────────────────

Section: F. Operational starting plan
What it does: Concrete first-five-steps if you want to begin: dominance harness → null control →
viability screen → ring-2-vs-ring-3 mu_lut=0 competition → queue a_b_tag C addition
────────────────────────────────────────

Section: G. Open methodological questions
What it does: θ/T_obs thresholds; ring-class representative selection; mosaic-front redefinition;
A/B egenome matching default; >2-tag generalisation
────────────────────────────────────────

Section: H. Program integration
What it does: Where this sits vs #7, Track-1/Track-2, #8, R1 & #3c lessons

Three load-bearing claims worth flagging up-front for your review

1. The A/B tag is already free whenever mu_lut = 0 — get_lut_hashes() is a stable per-cell tag
because hash never drifts. The "bit of book-keeping" in the doc only needs a real implementation in
the LUT-mutating regime, and that one is best done as a small additive C a_b_tag (parallel to S2a
lineage) — the same lightest-weight-merge-gated category as the ring-bounded operator. Both
Track-2-ish, both parallelizable.
2. The viability screen is not optional, it's a prerequisite. #3c just showed us that at the
R1-productive regime, a frozen GoL substrate goes 10/10 extinct in 3000 ticks. Most random LUTs
share that fate. Without screening, "A dominates B" readings collapse into "A is viable and B
isn't" — which is information, but only if you measured it on purpose. Step 0 of any prep using
fixed LUTs.
3. The null control is the part of the #8 lesson that transfers directly. #8 was killed by an
observable with no variance; this framework's analog is "what does μ_F^t look like when A and B are
identical preps with different tags?" That stochastic spread is the baseline against which all
real preps must be quoted. Without it, a μ_F ≈ +0.15 reading is uninterpretable.

I did not commit the file — it's your working document, so I'd rather you review/edit/integrate the
appendix and commit yourself (or tell me to). If you want me to commit it as-is, say the word. Two
factual hedges to flag: I noted env_mask write-API needs a one-line check before promoting
hostile-environment from "queue C add" to "feasible now"; and the KPZ β ≈ 1/3 claim is a standard
1+1D result but predicts a testable roughening exponent, not a guaranteed one for EvoCA dynamics.


## A. What the framework gets right

1. **Tag separated from genome.** Defining A and B as a *tag* (rather
   than as "the genome of A") cleanly accommodates the case where A
   and B start as the same genome (a null control, §C.10) and the
   case where they diverge under mutation (lineage erosion, §C.4).
2. **General — no ring-bounded mutation needed.** The whole framework
   is orthogonal to whether we ever build the ring-bounded mutation
   operator (Track-2). It is the natural measurement engine for
   Track-1 (`mu_lut=0` prepared competition) and remains the engine
   for Track-2 (evolution within a ring class) once that exists.
3. **Spatial dominance, not just population frequency.** Per-row
   `D_i` carries information that a global ratio loses (e.g.,
   front position, roughness). Worth keeping even when refined.
4. **Replaces the patch-transfer thread cleanly.** The patch program
   (#7) was always asking how much of "patch viability" was patch-
   stamp-artifact vs real competition; the half-lattice / band prep
   sidesteps stamping entirely.

## B. Geometry of the prep — what the doc leaves implicit (and why it matters)

Vertical bands with periodic boundaries have **two A/B interfaces**
(mid-lattice and the wraparound at the left/right edge), not one.
Many of the "interesting" dynamics — front roughening, finger
formation, KPZ-like growth — will be dominated by interface physics
at those two surfaces, not by bulk competition. This affects what
`σ_D^t` actually measures.

→ **Recommend explicit geometry options** rather than fixing
"vertical bands":

| Prep geometry | Interfaces | Primary use |
|---|---|---|
| **Half-lattice band** (current) | 2 vertical | Standard competition; bulk vs front mixed |
| Horizontal band | 2 horizontal | Same as above; choose to align with any anisotropy in `init='halfplane'` already used by `run_sim` |
| **Stripes** (alternating columns, width w) | ~N/w | High interface density — fast-equilibrating; tests local stability |
| **Mosaic** (random per-cell A/B at given fraction) | many, dispersed | No coherent fronts; pure frequency competition; the natural null geometry for a mean-field-like ratio |
| **Single-cell or small-patch invasion** | 1 closed boundary | Invasion-stability question (is A a successful rare invader?) — *complementary* to the 50/50 competition question |
| **Quadrants** | 4 perpendicular interfaces | Three-or-four-way comparisons in one run |

The doc's prep is the half-lattice band — fine as a default — but
the framework should *explicitly* accept multiple geometries; some
questions (invasion stability, frequency-only dominance) are
unanswerable in the half-band prep.

## C. Substantive refinements (each changes a metric, assumption, or step)

### C.1 `D_i` as raw count conflates competition with total-population change → add a frequency-based companion metric.

`D_i = N_A − N_B` shifts both because the *balance* shifts and
because the *total* row population shifts. Two regimes can give
identical `D_i` for different reasons. **Add the frequency-based
companion:**
$$F_i^t = \frac{N_A(i,t) - N_B(i,t)}{N_A(i,t)+N_B(i,t)} \in [-1, +1]$$
(defined only when row population ≥ a threshold; else NaN.) Report
$\mu_F^t$, $\sigma_F^t$ alongside $\mu_D^t$, $\sigma_D^t$.

`D` captures *absolute* dominance (a dying B is more dominated than
a halved B); `F` captures *competitive* dominance (60% vs 80% A is
unambiguously stronger dominance and is invisible to `D` if total
population also rises). Reporting both keeps the absolute-vs-relative
distinction explicit. This matters especially when comparing across
metaparameter regimes where total population differs (R1 showed
36 % extinction at the wrong end of the m_scale/food_inc box; #3c
showed B's 10/10 extinction at the productive regime). `F` is
comparable across regimes; `D` is not.

### C.2 `σ_D` is a weak second moment for the half-band geometry — add front-position and contact-fraction measures.

For a vertical-band prep with mostly-coherent fronts, `σ_D` over
rows is dominated by finite-sample noise of an essentially 1-D
process; it is not informative about *what the front is doing*.
**Add three more informative second moments:**

- **Front position** $\bar x_F(t)$: the centroid (along the column
  axis) of the A/B interface (median of the position where
  $F_i^t$ crosses 0 along each row); how does it drift?
- **Front roughness** $w_F(t) = \mathrm{std}_i\{x_F(i,t)\}$: how
  ragged is the interface? KPZ-like roughening scales as
  $w \sim t^\beta$ with $\beta \approx 1/3$ in 1+1 D — testable.
- **Contact fraction**: fraction of alive cells with at least one
  opposite-tag Moore-neighbour. Tracks how much A-vs-B interaction
  is happening (a *static* mosaic-equilibrium has high contact;
  domain-coarsening has falling contact).

For `mosaic` geometry the front is meaningless and *contact
fraction* / coarsening exponents become the natural second moment.

### C.3 The A/B tag through reproduction is *free* in some cases, needs a small additive C change in others.

The doc says "a bit of book-keeping" for evolutionary preps. There
are three regimes:

1. **`mu_lut = 0`, mu_egene/egenome = 0 (Track-1 pure-rule
   competition).** The LUT-hash itself is a stable per-cell tag —
   `get_lut_hashes()` distinguishes A's seed-LUT from B's
   automatically and forever. **No tag bookkeeping needed.**
2. **`mu_lut = 0`, evolutionary egene.** Same as (1) for the rule;
   egene drift doesn't change the LUT-hash. Still free.
3. **`mu_lut > 0` (LUT mutates).** LUT-hash drifts; lineage
   `parent_id` chains preserve descent but founders all have
   `birth_id = 0`, so we can't tell A-founder descendants from
   B-founder descendants from the lineage field alone. **This is
   the case that needs a real tag.**

The cleanest minimal solution for case 3 is **a small additive C
addition** that parallels the S2a lineage field: a per-cell `uint8`
`a_b_tag`, set at `set_alive_*` time (or by a new prep helper), and
copied parent→child at reproduction (alongside the existing
`lin_birth_id`/`lin_parent_id`/`lin_parent_hash` propagation in
`C/evoca.c:1455`). It never mutates. ~15 lines C + accessor + ctypes
binding, gated by `enable_ab_tag(True)` analogous to lineage. This
is **dynamics-altering only in the sense of "needs new C code" — its
runtime effect is purely a tag copy**, so it should be classified
as the lightest-weight kind of code change but still goes through
worktree + merge gate + pre-merge notebook (the protocol).

Until that tag exists: case 3 preps can still be **approximated** by
walking lineage chains back to founders + matching to the initial
A/B band by founder position. The reconstruction is exact for cells
whose entire ancestry was captured in the periodic lineage snapshots
(the #8 machinery already does this) and lossy for cells whose
chain was broken (a founder ancestor reproduced and then died
between snapshots). For initial work this approximation is enough.

→ **Practical recommendation:** start Track-1 cases (1)–(2) now
(no new code); queue the additive `a_b_tag` C change for the
evolutionary preps and treat it as the same kind of merge-gated work
as the ring-bounded mutation operator (own worktree, pre-merge
notebook). They are independent and could parallelize.

### C.4 Reframe "evolved vs wild-type" — under mutation, A and B are *descendants of* seeds, not the seeds.

Once mu > 0 acts on both sides, A is *no longer* "the evolved
population" by tick 1000 — it is "the descendants of the evolved
seed, now further mutated under the same dynamics as B's descendants
of the random seed." The dominance measure reports lineage propagation,
which can decouple from current trait quality: A could dominate
because the *evolved seed survived longer* (head-start effect) even
if A's current rule complexity has fallen back into B's range.

→ **Recommend two complementary readouts** for evolutionary preps:
- *Lineage dominance*: $\mu_F^t$ on the A/B tag (the doc's measure).
- *Trait dominance*: at each sampling time, compute the mean
  rule-complexity (`lut_complexity` level), mean cog_specificity,
  mean intake, etc. **separately on A-tagged and B-tagged cells**.
  This distinguishes "A's seed-derived rules are winning by being
  better" from "A-tagged lineages happen to dominate but their
  rules have converged toward B's".

The pair makes the experimental conclusion robust to the seed-
identity erosion that pure lineage dominance otherwise hides.

### C.5 R1 + #3c bite here too — viability screen is a prerequisite, and the regime choice is not free.

#3c showed that at the R1-productive regime
(`m_scale=3.0, food_inc=0.015`), a frozen-GoL substrate goes
10/10 extinct in 3000 ticks. That is directly analogous to what
will happen in Track-1 ring-stability with `mu_lut=0`: **most
random LUTs do not sustain population at any given regime** (Conway's
LUT-rule landscape is mostly non-viable). If A's seeded LUTs die
quickly and B's don't (or vice versa), $\mu_F^t \to \pm 1$ trivially
— not because of ring-class dominance, but because of viability
mismatch. The "dominance" reading is meaningless without controlling
for solo viability.

→ **Add a viability screen as a required Step 0** for every
prep that uses fixed LUTs:

> *Screen-0 (viability):* for each candidate LUT (and each ring
> class' representative set), seed a *uniform-LUT* full-lattice run
> at the chosen regime, run T_screen ticks (e.g. 3000), record final
> population, min population over the run, and the get_lut_hashes
> diagnostic. **Reject** any LUT whose solo viability fails a
> threshold (e.g., min_pop > some floor, final_pop > another floor).
> Only LUTs that pass Screen-0 go into A/B competitions. Report
> Screen-0 results alongside the competition results.

This converts "viable" into a measurable property, and an extinct A
or B arm becomes an *informative* outcome ("A is non-viable at this
regime") rather than noise.

Related: at `mu_lut = 0` (Track-1), **the regime choice cannot just
be inherited from R1 (which optimized for *egene-driven* dynamics
with `mu_lut > 0`).** Track-1 might need its own regime sweep, or at
least its own viability screen, before any A/B comparison is
interpretable. This is the same lesson #3c surfaced: *the regime
that's productive for one objective is not necessarily productive
for another* (cf. R1 disproving #2's low-food lean for egene-driven
runs). Pin a Track-1-appropriate regime before scaling.

### C.6 Ensemble vs uniform within a ring class — the doc is ambiguous, and the answer matters.

"A has ring-3 rules" is two distinct experiments:

- **Uniform-rule prep:** A is *one specific* ring-3 LUT applied to
  every A cell (B likewise). Tests the class-representative
  competition. Replicate over many randomly-drawn representatives.
- **Class-mosaic prep:** A is a *random mix* of ring-3 LUTs (a
  different ring-3 LUT per cell). Tests class-vs-class dominance
  *averaged over within-class diversity* — a different scientific
  question (the latter asks "averaging over the class, who wins";
  the former asks "what's the typical pair outcome").

These give different answers in general (Simpson's-paradox-like
effects are possible — class-mosaic A can lose even if every
individual representative would beat every individual B
representative).

→ **Recommend running both as paired experiments** for at least one
class pair (e.g. ring-2 vs ring-3) to characterize the
within-class-vs-between-class variance ratio.

### C.7 Symmetric initial conditions outside the rule axis.

The doc doesn't state which initial conditions (state v(x), private
food f(x), environmental food F(x), egenome) are matched between A
and B. **Pin them explicitly** — any asymmetry biases the test:

- `v(x)` symmetric: e.g. independent draws from the same
  distribution in both bands, or mirror-symmetric.
- `f(x)`, `F(x)`: identical across bands (or independent draws
  from the same distribution).
- egenome: identical (or with `egene` evolution turned off and a
  fixed seed value).
- `alive` mask: identical fraction in both bands.

This is implicit in `set_alive_halfplane`/`set_egenome_all` etc., but
should be a stated invariant of the prep so subtle asymmetries don't
sneak in (analogous to how the #3c shared-fork required asserting
byte-identical phase-1 substrates; the analog here is "asserting
asymmetry is only in the rule axis").

### C.8 Outcome categorization — replace "μ_D goes positive ⇒ A dominates" with a 5-way classifier.

A single sign-of-μ_D reading collapses real qualitative differences:

| Outcome | Definition (suggested) |
|---|---|
| A-extinct | A population reaches 0 before T_obs |
| B-extinct | B population reaches 0 before T_obs |
| A-dominance | `μ_F` at T_obs > +θ, both alive |
| B-dominance | `μ_F` at T_obs < −θ, both alive |
| Co-existence | `|μ_F|` < θ at T_obs and stable in last window |

Pin θ (e.g. 0.3) and T_obs (e.g. 8000 ticks for evolutionary preps,
shorter for `mu_lut=0` because no drift extending the time horizon).
This converts a continuous measure into a categorical readout that
is straightforward to tabulate across many runs.

### C.9 Time-scale convention pinned + long-timescale check (the user's point-2 lesson from #3c).

#3c taught us that egene/dyn signals can *plateau* before T_obs
without it being obvious. **Bake the early-vs-late window check
in from day one** of stability experiments: report `μ_F^t` over
windows [0, T/4], [3T/4, T] to detect plateau / reversal /
non-monotonicity. A monotone `μ_F^t` is qualitatively different
from a non-monotone one (transient invasion followed by reversal).

### C.10 The null control — *required* for significance.

The doc proposes informative experiments but no null experiment.
**Run a tag-only null first**: A and B are *identical* preps
(same LUT, same v(x), same everything) but tagged differently.
$\mu_F^t$ must hover near 0 with a finite stochastic spread; that
spread (over seeds) is the *null distribution* against which all
real A-vs-B differences must be benchmarked. Without this, a small
positive $\mu_F$ at T_obs is interpretively ambiguous.

This is the analogue of the #8 permutation null (which was the right
instinct there too): *measure how much "dominance" appears from
finite-population stochasticity alone*, then test real preps against
it.

## D. Additional A/B preparations

Building on the doc's three (rings, evolved-vs-wild-type, gradual-vs-
rapid), here are eight more, sorted by feasibility *under current
capability* vs *needing a small addition*. Each is one named prep —
the user can fold the useful ones into the doc.

**Feasible now (no new C code):**

1. **Density A vs B** — same LUT/egene; A seeded at fraction 1.0
   (full density in its band), B at 0.5. Tests how much of "dominance"
   is initial-population-mass alone. Important *null variant*.
2. **Initial-state structure A vs B** — same LUT; A seeded with
   structured v (e.g. glider population), B with random v. Tests
   whether prepared spatial structure carries forward into population
   advantage independent of rule.
3. **Egenome class A vs B** — same LUT; A has high-specificity
   egenome (low active count, narrow match), B has broad (high
   active count). Tests selectivity-vs-generalism dominance in the
   eating axis. Connects to #1/causal-control's cognition story.
4. **Co-evolved-pair A vs B** — A = (LUT, egenome) jointly evolved
   from a long joint-productive run (#1-winners); B = (LUT-only-
   evolved, wild-type egenome). Tests whether co-adaptation gives a
   real *competitive* advantage, complementing the joint-retention
   reading that #8 (when reborn with a real metric) will give.
5. **Diversity A vs B** — A = single uniform LUT across its band
   (low intra-A diversity); B = mosaic of LUTs within the same ring
   class (high intra-B diversity). Tests whether within-band
   diversity helps or hinders a band against a homogeneous
   competitor — directly relevant to evolvability literature.
6. **Invasion variant of all of the above** — same preps but A in a
   small patch (e.g. 20×20) immersed in a sea of B. Asks the
   *invasion-stability* question (does rare A take over?), which is
   not the same question as 50/50 competition: a stable co-existence
   at 50/50 can mask total invasion-failure at rare. (Cf. the
   classic Maynard-Smith ESS distinction: invadability ≠ resident-
   stability.)

**Needs a small additive C change (queue alongside `a_b_tag`):**

7. **Hostile-environment A vs B** — symmetric `env_mask` modulation:
   A region has nominal `food_inc`, B region has reduced (or zero)
   `food_inc`. Tests rule robustness to environmental scarcity.
   `env_mask` already exists per `evoca_get_env_mask`; we may already
   have the API to do this without new C — worth checking before
   queueing. If yes, promote to "feasible now."
8. **Per-region mutation rate A vs B** — A operates under high
   `mu_lut`, B under low (same regime otherwise). Tests
   exploration-vs-exploitation dominance. Needs per-region mutation,
   which we don't have (`mu_lut` is global) — small additive C
   change ('per-cell or per-region mutation rate') to make it
   feasible.

## E. Adversarial peer pass — three rival objections, each producing a revision

Per the project critique discipline ("what would a strong domain
expert with an opposing draft object to?"):

1. **"Your dominance metric is just a 2D extension of population
   genetics — what does the CA layer add beyond noise? Show me one
   *spatial* observable that wouldn't appear in a well-mixed
   ODE model."** → Revision: §C.2's *front-position and contact-
   fraction* aren't reproducible from a frequency-only model; KPZ
   roughening (if found) is a 2D-specific signature. Promote front
   metrics to *primary*, not secondary, to justify the spatial CA
   framing.

2. **"You're conflating fitness with viability. Most random rules
   are unviable; your A-vs-B differences are dominated by which
   class has more viable seeds, not by competitive interaction."** →
   Revision: §C.5's Screen-0 viability filter is mandatory before any
   ring-class result is published; report Screen-0 viability rates
   alongside dominance.

3. **"Without a null, your finite-population statistics will lie:
   any 256×256 stochastic run with two equal tagged populations
   shows nonzero μ_F by t=T just from drift, and you have no way to
   tell that from real selection."** → Revision: §C.10's identical-
   prep null run is mandatory. Establish the stochastic null spread
   first, then quote *real preps as standard deviations above the
   null*, not as absolute μ_F values.

## F. Operational starting plan (what I would do first)

If asked to operationalize this tomorrow, the smallest set of steps
that unblock the rest:

1. **Implement the dominance harness** as a Python helper on `main`
   (no C change for Track-1 (1)–(2) regimes). Inputs: an A/B tag
   computed each snapshot from `get_lut_hashes()` (case 1/2) or
   from an explicit tag (case 3 once the C addition lands). Outputs:
   `D_i^t`, `F_i^t`, `μ`/`σ` of each, front position, front roughness,
   contact fraction, outcome category (§C.8). Pin sampling rate.
2. **Run the null control (§C.10).** Identical A and B preps,
   varying only the tag. 5–10 seeds. Characterize the stochastic
   spread of `μ_F^T`.
3. **Run the viability screen (§C.5).** Each candidate uniform LUT
   (a handful per ring class) → solo full-lattice run at the chosen
   regime → keep the viable ones for competition.
4. **Run the simplest informative A-vs-B experiment**: ring-2 vs
   ring-3 uniform-LUT (one viable representative each) at `mu_lut=0`,
   multi-seed, multi-food_inc. Confirm the framework works on the
   case that requires no new code at all.
5. **Queue** the `a_b_tag` C addition (small additive, own worktree,
   pre-merge notebook) to unlock the evolutionary preps cleanly.

## G. Open methodological questions to pin before scaling

- **θ, T_obs values** for the 5-way outcome classifier (§C.8).
- **Ring-class representative selection**: from which population /
  random distribution are ring-N LUTs drawn for §C.6? Random over
  the reduced basis? Random + viability-screened?
- **What counts as a "front" in mosaic geometry** (§B): probably
  redefine the spatial observable for that prep (contact fraction
  replaces front position).
- **Match between A's and B's egenomes**: identical (controls for
  cognition axis) vs symmetrically randomized (matched distribution
  but distinct realizations) — both defensible; pin the default.
- **Whether to extend to >2 tags** (e.g., ABC three-way) — the
  framework generalizes, but the dominance measure needs a different
  form (pairwise vs simplex).

## H. How this relates to the rest of the program

- **Replaces #7** (patch transfer), as the doc states — and does so
  cleanly because the band/mosaic preps need no stamping.
- **Track-1 (`mu_lut=0`)** is fully reachable now under current
  workflows; needs only the dominance harness + viability screen +
  null control (all Python, on `main`). Connect to the
  [[ring-stability-two-track-plan]] memory.
- **Track-2 (ring-bounded mutation)** remains its own merge-gated
  workstream; once it lands, evolutionary preps within a fixed-class
  region become possible.
- **#8** (when reborn with a non-degenerate observable) and this
  framework are *complementary*: #8 measures co-evolution as
  joint-retention along lineages within a single population; the
  A-vs-B dominance framework measures *competition between
  populations* prepared with different rules/genomes. Both are
  needed; neither subsumes the other.
- **R1 / #3c lessons** baked in throughout (§C.1, §C.5, §C.9):
  population dynamics ≠ competitive dynamics; regime choice is
  conditional on objective; long-timescale plateaus must be checked.
- **#8 metric-flaw lesson** baked into §C.1 and §C.10: choose
  observables with variance, pin a null before drawing conclusions.

End appended critique.
