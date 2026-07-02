
# EvoCA research directions

This document serves to outline research objectives of EvoCA, with the goals of

- sharpening EvoCA as a research tool by adding features that enable making progress toward the objectives, or even better, attaining them.

- creating a map of research paths that are clear enough to be pursued and managed by Claude, through running sequences of simulations in parallel, to make progress on the paths.

Each section will include some suggestions for enhancements of the platform, with research motivations.

This document will serve as a prompt to Claude, to be regarded as deep analysis.

I expect recommendations, including feedback on my suggestions for enhancements and suggestions for changes and/or other enhancements.

## dev logistics

We will have to discuss dev strategy.  Some of the research directions might benefit from creating their own branch, then maybe or maybe not merging in the future.

# Overall big picture issues

## Open ended evolution

One of the fundamental motivations of EvoCA is to construct a model that is rich enough to display some form of open-ended evolution, and yet to be simple enough to understand, and to measure relationships between model metaparameters and open-endedness.

## Evolutionary landscapes

The nature of open-endedness, when present, may be strongly dependent on the evolutionary landscape.  A couple of different possibilities (not meant to be exhaustive)

- The population is on the slope of a Mt Fuji fitness peak, its evolution is open-ended in the sense that constantly varying genomes can claw their way up before they are replaced.  But:  is the population moving up the mountain?  If so how is its movement measured?

- The fitness landscape is rugged, with many independent peaks, and sub-populations surviving because they exploit one of many peaks.  But:  is the population moving toward higher and higher peaks?  Or does it typically get stuck?

Questions:

- How can different landscapes be distinguished?

- Do we need to look at how properties of evolution (like the landscape) scale with world size, i.e. lattice size?

# Different issues generating research paths

## Activity

Many of the analytics of EvoCA (in particular the different versions of activity) are aiming to get a bead on evolutionary innovation and progress, and to detect and measure ongoing evolutionary innovation and progress.

We have various versions of activity
- activity (genome activity)
- egenome
- dyn_activity
- pat_activity
and for activity, we have
- n_activity to measure excess over neutral model

One possible path re. activity:
**make neutral model variants of all activity versions, enabling the measurement of excess activity for all versions.**

## Resource flow

The model is structured to have resources flow from the environment into live cells, with various kinds of taxes imposed, so the evolutionary learning task, writ large, is to learn how to get food from the environment to overcome the taxes to such a degree that reproduction can take place.

Resources represent a sort of energy.  Can we observe flow of energy into genomic innovation that may be taxed, but generates more food intake than the additional tax?

Terrestrial life is generally regarded as a dissipative structure living in the flux of energy coming from the sun.  Can we see a concrete parallel in EvoCA?  Can we measure a relationship between organismic energy usage and organismic complexity?

### egenomes

current dev on egenomes is aimed at exploring a version of resource flow.  If egenomes are taxed by number of non * egenes, live cells with more complex egenomes must perforce be harvesting food more efficiently than less complex egenomes.

Are all the pieces in place to measure this relationship directly, and ascertain when successful egenome evolution is actually happening? 

Experiments should include 

- starting with simple egenomes, and seeing if they grow to more complex ones
- fixing egenomes with particular spatial patterns (e.g. checkerboard), set mu_egene = 0, and see if LUT evolution leads to local patterns that lead to enhanced feeding.

### LUT genenomes

We need to impose a tax structure on egenomes related to the complexity of the egenome, as measured by what neighborhood ring it depends on.  Similar to taxing egenomes for their complexity.

Then experiments should include

- starting with simple genomes (inner ring), does evolution proceed to more complex genomes.
(genelife experiments seem to indicate that in that context, without resource flow, the answer is not necessarily).
- what is the relationship between mu_lut and mu_egene?

## optimization experiments

We have benifitted greatly from optimization experiments that helped find beautiful reaction diffusion spatial dynamics, with evolution taking place on the fronts.

Remaining questions:
- what happens if we really push pure evolutionary metrics without spatial pattern metrics?  Do we always optimize to uninteresting spatial dynamics (static structures)?

My intuition is that robust evolution requires life-death turnover, so I don't expect to see static structures.  If I see static structures, I will suspect the metrics driving the optimization.

- what happens when we take some of the nice RD generating metaparams and use them as seeds for pure evolution oriented optimization?  We have done one round of that, but I would like to see what happens when we really push

- What happens when we optimize
    - egenome evolution keeping lut genome constant
    - lut genome keeping egenome constant
    - combining evolutionary metrics for simultaneous egenome and lut evolution

- should we consider larger worlds for optimization?  I would like to scope optimization campaigns that could take days on a 32-core machine with 128GB memory.  But obviously I want to burn cpu time for results, not just to take more time.

- there will be phase transitions as metaparams are varied.  Observables often diverge, or have a discontinuity at such transitions.  I would like to identify such transitions and measure critical exponents when appropriate.

# Mutation of mutation

I would like to include mutation parameters in the genome, so that they are themselves subject to mutation.

more or less obvious questions:
- does the distribution of mutation params in the population converge, and if so to what shape?
- how do mutation rates in the population respond to environmental stress (e.g. sudden perturbation in global food for a time).
- how do the distribution or mutation rates depend on resource environment in general?

# Evolution success 

if a population evolves for a while, what happens if

- a patch of it is placed in an environment of random live cells
- a patch of it is placed in an environment of independently evolved cells

---

# Claude analysis (2026-05-16)

This section responds to the document above as deep analysis. It is
grounded in the current code (`C/evoca.c`, `python/evoca_explore.py`,
the probe set in `Docs/probes.md`) and in the four completed scan
campaigns (`Scans/2026-04-27_*`, `Scans/2026-04-27_scan_analysis.md`),
not on generic OEE advice. Literature anchors are in
`Docs/references/references.md` (open-access papers downloaded there;
paywalled ones linked). Where my reading of your framing differs from
yours I say so explicitly rather than smoothing it over.

Section map (purpose in one line each):

1. **Reframing** — three places the prompt's framing should change before any code is written.
2. **The neutral-model program — and the bug in the one you have** — `excess_activity_slope` is broken in exactly the regime your scans optimise in; root cause and fix.
3. **One primary OEE order parameter** — what the scan harness should optimise instead of the current composite.
4. **Distinguishing landscapes as a continuous map, not a binary** — ruggedness as a scannable scalar; Mt-Fuji-vs-rugged operationalised.
5. **Scaling and criticality are one program** — finite-size scaling; the correlation-length plateau is probably a criticality signature, not a grid ceiling.
6. **Resource flow → complexity** — the dissipative-structure measurable, and why your scans already found the resource-flux window.
7. **Egenome evolution: instrumentation audit** — what is and isn't in place to answer "is cognition evolving?"
8. **The LUT-complexity-tax design gap** — `tax_lut` does not tax what you asked it to tax.
9. **Optimization program** — wiring cognitive metrics, separable vs joint optimisation, the determinism bug, large-world budget.
10. **Mutation of mutation** — C-side design, the error-threshold prior, stress protocol.
11. **Patch-transfer success assay** — concrete API, scoring, the controls that make it mean something.
12. **Dev logistics** — branch strategy per direction.
13. **Prioritised roadmap** — table + a "do these first" short list.

### Notes from claude code window:

⏺ Appended as a new section "Claude analysis (2026-05-16)" to Docs/EvoCA_research_directions.md (13 sections + section
   map), and created Docs/references/references.md with two open-access PDFs downloaded into Docs/references/:

- taylor2016_OEE_perspectives_york.pdf (CC-BY)
- banzhaf2016_defining_simulating_OEE_novelty.pdf

Channon 2001 was already in Docs/. The Bedau–Packard series and the paywalled landscape/criticality/mutation
references are cited with publisher links (SFI page was access-denied at fetch time — noted, not guessed). No
fabricated citations.

The analysis pushes back on the framing in three places and surfaces a concrete defect, rather than just elaborating
your list:

1. A real bug, not a quirk. Your excess_activity_slope inverts at mu_lut≈0.06 because the neutral shadow in
C/evoca.c:1715 reproduces with the same mutation operator as the real run — at high μ that's a no-heredity model,
not a neutral one. This makes the evo_focus "mu_lut=0.06 optimum" rest on two of four non-informative metrics. Fix
is ~10 lines (per-component normalisation, the Channon-2006 correction). Flagged as the most urgent item because
it's currently producing a wrong optimisation conclusion.
2. The lineage field is the backbone. EvoCA tracks genome presence but not parent→child lineage; the landscape map,
resource→complexity trajectory, lineage-resolved cognition, and patch-transfer success criterion all reduce to
lineage bookkeeping. Single highest-leverage addition.
3. Three reframings: Mt-Fuji-vs-rugged → a continuous ruggedness coordinate (scannable); scaling and
critical-exponents → one finite-size-scaling program (your ξ≈N/3 plateau is plausibly a criticality signature, not a
grid ceiling — and the two are distinguishable via susceptibility/Binder cumulant); the resource window the
dissipative-structure framing predicts is the food_inc≈0.013 bracket your scans already found.

Also surfaced: tax_lut taxes bit-count, not the ring-depth complexity you asked to tax (would confound the genelife
A/B); COG_METRICS is defined but wired into nothing; the per-process RNG non-determinism is a correctness blocker
for exponent work, not a caveat.

The roadmap ends with a "do this week, zero-code" short list — the reciprocal egene/LUT isolation controls need no
new code and bound every later cognition claim.

Nothing committed — the working tree has the doc edit, Docs/references/, and the earlier untracked .evoca/scan
files. Want me to commit the analysis doc + references (and only those)?


---

## 1. Reframing (read this before the rest)

Three framing changes that propagate through everything below:

- **"Is the population moving up the mountain?" needs a fitness proxy
  EvoCA does not yet expose.** There is no explicit fitness in the
  model; fitness is realised, not assigned — it is the lineage's
  reproductive success and persistence. Every landscape and
  open-endedness question therefore reduces to *lineage bookkeeping*
  (who descends from whom, how many offspring, how long a genome's
  descendants persist). EvoCA currently tracks genome **presence**
  (activity) and **population**, but not **parent→child lineage**. The
  single highest-leverage instrumentation addition in this whole
  document is a per-cell parent-hash / lineage-ID field written in
  Phase 4. Almost every "are we making progress" question becomes
  measurable once that exists, and stays unanswerable (or only
  proxied) without it. I treat this as the backbone dependency in the
  roadmap.

- **"Mt Fuji vs rugged" is a false binary.** Ruggedness is a
  continuous, *tunable* property (Kauffman & Levin 1987; in NK terms,
  K is a knob), and in EvoCA the knob is the metaparameter vector
  itself plus the spatial structure of the food field. The right
  deliverable is not a classifier that stamps a run "Fuji" or
  "rugged"; it is a **ruggedness coordinate** computed per run and
  mapped across metaparameter space — same machinery as the existing
  scans, one more column.

- **"Scaling with lattice size" and "critical exponents" are the same
  experiment.** Critical exponents are *defined* through finite-size
  scaling. Your scan analysis already reports a correlation-length
  plateau at ≈ N/3 and calls it a probable "grid-size ceiling." The
  competing — and more interesting — hypothesis is that the system
  sits near a critical point where ξ is cut off by L, which is *also*
  exactly ξ ≈ N/const. These two readings are experimentally
  distinguishable and the test is one protocol (§5), not two
  sections.

## 2. The neutral-model program — and the bug in the one you have

Your proposal: *"make neutral model variants of all activity versions,
enabling the measurement of excess activity for all versions."* Agreed
in direction, but there is a prior problem that has to be fixed first
because it is silently distorting the optimisation you have already
run.

**The bug.** `Scans/2026-04-27_scan_analysis.md` §"One Channon-metric
anomaly" and the `evo_focus` notes both record that
`excess_activity_slope` *fell* (top-5 mean 448 → 228, and went
negative for individual top configs) precisely as `mu_lut` was pushed
to 0.06 — the regime your last scan declared optimal. The notes
attribute this to "a known quirk … the excess metric is being
saturated by raw turnover on both sides." It is worth being precise
about the root cause, because it is not a quirk of the Bedau–Packard
metric, it is a defect in *this particular shadow construction*:

`neutral_apply_demography()` in `C/evoca.c` reproduces shadow members
with **the same `gmu_lut` LUT-bit-flip procedure as the real run**
(`evoca.c:1715` "Mutate one shadow LUT in place: same procedure as
evoca_step's"). At `mu_lut = 0.06`, Poisson(0.06 × 250) ≈ 15 bit
flips per birth: every birth in *both* populations is essentially a
fresh genome. The shadow is then not a *neutral* model (heredity +
drift, selection removed) — it is a *no-heredity* model. ΣN inflates
because random parent choice in the shadow never lets a lineage
consolidate, so it manufactures new hash buckets at the maximum
possible rate; ΣG − ΣN flattens or inverts not because real evolution
stopped but because the baseline became wrong. The Bedau–Snyder–Packard
(1998) construction you co-authored assumes the shadow preserves the
**same diversity/turnover statistics under non-selective
inheritance**; this implementation preserves the same *mutation
operator*, which at high μ is a different thing.

**Fixes, in order of preference:**

1. **Switch the optimisation target to a turnover-invariant
   statistic.** Raw ΣG − ΣN is the wrong summary at high μ. New
   activity per component A_new and mean cumulative activity Ā (your
   1998 paper's own preferred statistics) are designed to be invariant
   to the bucket-creation rate. Concretely: report
   `excess = (ΣG/D_G) − (ΣN/D_N)` (per-component-normalised, where D
   is diversity = count of buckets with pop>0), not `ΣG − ΣN`. This is
   a ~10-line change in `run_sim` and is exactly the
   "component-normalised" correction Channon (2006) applied to fix the
   identical failure mode in Geb. The downloaded Banzhaf (2016) and
   Taylor (2016) papers both flag raw cumulative activity as
   turnover-fragile.
2. **Make the shadow actually neutral.** Decouple the shadow's
   reproduction from `gmu_lut`: in the shadow, copy the parent LUT
   *without mutation* and instead inject new buckets at the rate that
   matches the real run's *realised* novel-genome rate (measured, not
   assumed). This is the literal Bedau–Packard shadow. More code, but
   it makes excess activity meaningful across the whole μ range
   instead of only below μ ≈ 0.03.
3. **Re-examine the `evo_focus` conclusion.** The headline finding
   "`mu_lut = 0.06` wins" rests on a composite where one of four
   metrics (`excess_activity_slope`) was inverted by this bug and the
   other (`unique_top_genomes`) was saturated at the sample cap (your
   own analysis, observation 8). Effectively two of four EVO metrics
   were non-informative in that scan. The mu_lut=0.06 optimum may be
   real, but it should be re-derived with the corrected metric before
   it drives a multi-day campaign.

**Which activity versions need a shadow, and how:**

| Version | Shadow needed? | Construction |
|---|---|---|
| `activity` (LUT genome) | yes — exists, fix per above | per-component-normalised excess |
| `eg_activity` (729 ternary egene keys) | **yes — missing** | fixed 729-bucket space → shadow is trivial: a multinomial-drift baseline (Wright–Fisher on 729 alleles at the realised mutation flux). No hash table needed; closed-form variance even available. |
| `dyn_activity` (500 LUT-entry buckets) | yes — missing | same idea: 500 fixed buckets, drift baseline from realised per-entry flip rate. This is the *sharpest* one: it directly measures whether *effective* rule transitions are under selection, independent of silent hash churn. |
| `pat_activity` (local spatial patterns) | yes — missing | shadow = same patterns under a phase-randomised (Fourier-surrogate) alive mask preserving spatial spectrum but destroying lineage structure. |

Because `eg_activity` and `dyn_activity` live in **fixed finite bucket
spaces**, their neutral baselines are far cheaper and more rigorous
than the open-ended LUT-hash shadow — no realloc, no random-kill
bookkeeping, and an analytic drift variance you can z-score against.
**Recommendation: do `dyn_activity` and `eg_activity` shadows first
(closed-form, cheap, rigorous), and treat the LUT-hash shadow fix as a
parallel track.** This inverts the prompt's implied order ("all
versions" as if uniform); the fixed-space ones are both easier and
more diagnostic.

## 3. One primary OEE order parameter

You want metaparameter↔open-endedness *relationships* and a loop Claude
can manage. That requires a single scalar the headless harness computes
per run, robust to the turnover artefact of §2. The current scan
composite (four metrics, two of them shown non-informative at the
chosen optimum) is not that.

Proposed primary objective, computable in `run_sim` today plus the
lineage field of §1:

> **Unbounded-activity diagnostic U\*** = sign and magnitude of the
> long-run slope of *component-normalised* cumulative excess activity
> over a long window, AND the long-run slope of diversity D(t).
> Class-3 (unbounded, per Bedau–Snyder–Packard 1998): both stay
> bounded away from zero as t → ∞. Class-2 (bounded): excess > 0 but
> diversity-normalised new activity decays. Class-1: excess ≈ 0.

Operationally: run long (≥ 50 k ticks, not 5 k — see §9), fit the
slope on the second half only (discard transient), use the §2.1
normalisation. This gives the scan a single principled axis with a
literature-grounded interpretation, replacing the ad-hoc unweighted
4-metric mean. Keep the spatial/cognitive metrics as *reported
context*, not as score components — your own cross-scan finding
(observation 6: "no Pareto trade-off in the regime we're searching")
says you do not need them in the objective.

## 4. Distinguishing landscapes as a continuous map

Reframed per §1: compute a ruggedness coordinate per run, map it over
metaparameter space. Three measurable coordinates, in increasing order
of what they need:

- **R1 — realised-fitness autocorrelation along lineages** (needs the
  lineage field). Define a genome's realised fitness as its
  per-capita offspring rate over its lifetime. Walk parent→child
  chains; the lag-1 autocorrelation of realised fitness along chains
  is the classic landscape-correlation statistic. High autocorrelation
  + a persistent advancing frontier ⇒ Fuji-like (smooth, directional).
  Low autocorrelation + many coexisting mid-fitness lineages ⇒ rugged.
  This is a single scalar per run.
- **R2 — beneficial-mutation fraction vs adaptedness** (needs lineage
  field). The Kauffman–Levin signature: on a smooth landscape the
  fraction of mutations that are beneficial declines monotonically and
  predictably as the lineage climbs; on a rugged landscape it stays
  noisy and roughly constant because every genome sits near some local
  optimum. Bin children by parent realised fitness, plot Δfitness
  distribution per bin; the slope is the coordinate.
- **R3 — the activity-wave profile you already render.** The
  `q_activity` probe (deciles of per-genome cumulative activity,
  modelled on genelife-2 Fig. 3) already discriminates these two
  cases without lineage tracking, today: a single advancing front
  shows as a coherent high-decile band pulling away from the low
  deciles (Fuji-like directional sweep); a rugged landscape shows many
  comparable mid-deciles with no breakaway band. **Use `q_activity`
  for immediate qualitative classification before the lineage field
  lands; R1/R2 for the quantitative scan column after.**

The deliverable is a heatmap: ruggedness coordinate over the
(food_inc, m_scale, mu_lut, tax) box your scans already bracket.

> **[Superseded by NP-reply R2 below: the "ruggedness rises with
> mu_lut, falls with food_inc" prediction is withdrawn. It smuggled
> back the a priori metaparam→ruggedness link NP correctly rejects.
> Ruggedness in EvoCA is emergent with no a priori functional tie to
> the metaparams; the map is purely empirical/exploratory, with no
> expected monotone form, and whether it even varies smoothly (vs
> jumping at a bifurcation — R3) is itself the result.]**

Your existing scan data may already let you compute R3 retrospectively
if the q_activity decile traces were logged; if not, this is a re-run,
not new code.

## 5. Scaling and criticality are one program

> **[Corrected by NP-reply R3 below: these are NOT one experiment.
> Finite-size scaling is only one of three distinct exponent programs;
> the parameter-driven bifurcation one is the more relevant. Read §5
> together with NP-reply R3.]**

The correlation-length plateau at ξ ≈ N/3 across scans 2–4 is
currently interpreted (scan_analysis observation 2) as a grid-size
ceiling. That is one of two hypotheses:

- **H_ceiling**: ξ is intrinsically larger than N/3 but clipped by the
  half-period of a periodic box. Prediction: ξ_measured grows ∝ N
  with the *same* prefactor at every N, and the *physical* correlation
  length (some run-intrinsic scale) is N-independent.
- **H_critical**: the productive regime sits near a continuous phase
  transition; the true correlation length diverges, so ξ_measured is
  pinned to ≈ const·N at *every* N (finite-size-limited at all sizes).
  This is the textbook finite-size signature of criticality.

**Both give ξ ≈ N/const, so the plateau alone cannot distinguish
them** — which is why the current "grid-size ceiling" conclusion is
under-determined. The discriminator is the *susceptibility* and the
*Binder cumulant*, not ξ itself:

Finite-size-scaling protocol (one campaign):
1. Pick the transition-suspect axis. Best candidates from your data:
   `food_inc` (tightly bracketed at 0.013 — a sharp optimum bordering
   extinction is a classic transition) and `mu_lut` (the diversity
   metrics jump sharply across it). Sweep the axis finely through the
   suspect value.
2. For each of N ∈ {128, 256, 512, 1024}: order parameter = alive
   fraction (or diversity D); measure its variance × N² (the
   susceptibility χ) and the Binder cumulant U₄ = 1 − ⟨m⁴⟩/(3⟨m²⟩²).
3. **Signatures:** χ peaks that grow with N and a Binder-cumulant
   crossing at a common point ⇒ genuine transition (H_critical) ⇒
   extract β/ν, γ/ν by data collapse. χ peaks that *don't* grow with
   N ⇒ no transition (H_ceiling); the plateau is then genuinely a box
   artefact and large-N scanning has diminishing returns.

This directly answers the prompt's "do we need to look at how
properties of evolution scale with lattice size" (yes, but as a
criticality test, not a generic robustness check) and "identify
transitions and measure critical exponents when appropriate" (the
"when appropriate" is decided by the Binder crossing, step 3). It also
tells you whether the multi-day large-world campaigns you want to
scope (prompt §optimization) are worth the CPU: only if H_critical, in
which case N=1024 reveals new structure; if H_ceiling, N=512 is the
point of diminishing returns and the CPU is better spent on longer
runs at N=256.

## 6. Resource flow → complexity

This is the strongest place to integrate prior data, because your
scans have *already* measured the key quantity without naming it.

The dissipative-structure framing (Schneider & Kay 1994; Chaisson's
energy-rate-density Φ_m) predicts a **window**: complexity-bearing
structure can only be selected for when the resource flux is high
enough to pay the metabolic cost of complexity but low enough that
complexity still confers a reproductive differential. Below the
window: extinction. Above: selection is slack, any genome reproduces,
no pressure toward cognition. EvoCA makes this concrete and falsifiable
because the flux is literally `food_inc`.

Your scans found `food_inc` **tightly bracketed at 0.013** (scan
analysis: "tight peak in both scans"; the "interesting middle" box is
food_inc ∈ [0.008, 0.018]). That bracket *is* the predicted resource
window — you have already located it empirically. The research move is
not to discover it but to *test the mechanism*: across the food_inc
sweep, plot mean cognitive load (the `egene` probe's `load` strip) and
the §4 ruggedness coordinate against food_inc. The dissipative-structure
prediction is a non-monotone (inverted-U) relationship between
food_inc and selectable complexity, peaking inside [0.008, 0.018], with
complexity *collapsing* (not just survival failing) on both sides.

Concrete per-genome measurable for "energy into genomic innovation
that is taxed but pays for itself" (your prompt's exact question):

> **Net metabolic return of a lineage** = (cumulative mouthful eaten
> by all members of the lineage) − (cumulative tax paid). Innovation
> "pays for itself" iff a more complex descendant has higher net
> return than its simpler ancestor.

`eg_food` already gives cumulative food per egene key; `tax_per_egene`
and `tax_lut` give the cost side analytically from the genome. With
the lineage field this becomes a directly plottable
return-vs-complexity trajectory — the concrete EvoCA parallel to
"energy usage vs organismic complexity" you asked for. Without the
lineage field it is only a population-level scatter (still useful, less
sharp).

## 7. Egenome evolution: instrumentation audit

Your prompt asks: *"Are all the pieces in place to measure this
relationship directly, and ascertain when successful egenome evolution
is actually happening?"* Honest audit:

**In place:** `egene` probe (spec/load/intake), `eg_activity` (729
keys), `eg_food` (cumulative food per egene), `egenome_stats`,
per-slot introspection, `tax_per_egene`. The validation table in
`egene_discussion.md` (V1–V3) shows the knobs behave sensibly.

**Missing for "ascertain when it is *actually* happening":**

1. **A neutral baseline for `eg_activity`** (§2) — without it, rising
   egene activity cannot be distinguished from drift through the
   729-key space. This is the single missing piece that most directly
   blocks your stated question, and it is the *cheap* shadow (fixed
   space, closed-form drift variance).
2. **Lineage-resolved intake** (§1) — "successful egenome evolution"
   means a more-specified egene lineage out-reproduces its less-
   specified ancestor *because* of intake, not coincidentally. The
   `egene` probe shows population means; it cannot separate "cognition
   improved" from "cognitively-poor cells died for unrelated reasons."
3. **A causal control**: your proposed experiment "fix egenomes
   (checkerboard), set mu_egene=0, let LUT evolve" is exactly the
   right causal isolation and it needs *no new code* — `set_egenome_pair_all`
   + `mu_egene=mu_egenome=0` already does it. The complementary
   control (fix LUT, evolve only egene: `set_lut_all` +
   `mu_lut=0`) likewise exists. **These two reciprocal controls should
   be run before any joint-evolution interpretation**; they bound how
   much of any observed cognitive gain is attributable to each genome.
   This is the highest-value zero-code experiment in the document — do
   it this week.

Your "start simple, see if it grows complex" experiment: seed with
`set_egenome_pair_all(0, 0b000001)` (centre-only, minimal cognition,
already supported) and watch the `egene` `spec` strip. But note the
initial-viability concern in `egene_discussion.md` is real and
under-tested: a centre-only egene scores at most +1 per match, so
mouthful is tiny; at the scan-productive corner with `tax_per_egene>0`
the seed may not clear tax before going extinct. **Recommend a
viability pre-screen**: a short headless sweep of
(food_inc, tax_per_egene) at fixed minimal egene to find where the
minimal-cognition seed survives at all, *then* run the
complexity-growth experiment only inside that survivable subregion.
Otherwise "simple did not grow complex" is confounded with "simple
died immediately," which is not a statement about cognitive evolution.

## 8. The LUT-complexity-tax design gap

Your prompt: *"impose a tax structure on egenomes related to the
complexity of the egenome, as measured by what neighborhood ring it
depends on."* (I read "egenomes" here as the LUT rule, per the section
heading "LUT genenomes" — flagging the ambiguity rather than guessing
silently.)

**The gap:** `tax_lut` already exists, but it taxes the **popcount of
'1' bits in the LUT**, not **ring-dependence depth**. These are
different. A rule can have many '1' bits yet condition only on n1
(complexity level 1 in the `lut_complexity` probe's own
classification); another can have few '1' bits but genuinely depend on
n3 (level 3). Taxing bit-count penalises a property that is largely
orthogonal to the "what ring does it depend on" complexity you
actually asked to tax. Using `tax_lut` as a proxy for ring-complexity
would silently confound the genelife-comparison experiment you want to
run.

**Principled alternatives (pick one):**

- **Ring-dependence tax** (closest to your stated intent): per cell,
  compute the minimum ring set the LUT actually conditions on — this
  is *exactly the computation `lut_complexity` already performs* for
  the probe (level 1/2/3). Tax = `tax_ring × level`. Reuses existing
  code; the probe and the tax then measure the same thing, which is
  what makes the "does evolution climb the ring ladder under resource
  flow?" experiment interpretable.
- **Effective-dependence tax** (sharper, more code): mutual information
  between each ring's count and the LUT output, summed. Penalises
  *used* dependence, not declared structure — the LUT analogue of
  what the ternary egene mask already does for cognition. This is the
  more elegant design and aligns the LUT and egene cost models
  conceptually, but it needs an MI estimator over the realised input
  distribution per cell per step (expensive; restricted-mutation's
  `lut_active` mask gives a cheap approximation of the support).

The genelife null you cite ("without resource flow, simple→complex
does not necessarily happen") becomes a clean A/B: same ring-tax,
`food_inc = 0` (genelife-like, no resource flow) vs `food_inc` in the
window. The hypothesis "resource flow is what enables the climb up the
ring ladder" is then a single controlled comparison with a
literature-anchored null. This is the experiment that most directly
addresses the EvoCA-vs-genelife thesis in your project's premise — but
only if the tax measures ring-depth, not bit-count.

## 9. Optimization program

Point-by-point on the prompt's optimisation questions, grounded in
`evoca_explore.py`:

- **"Push pure evolutionary metrics, do we always get static
  structures?"** Your intuition (robust evolution needs life–death
  turnover, so static structure ⇒ suspect the metrics) is correct and
  the §2 bug is evidence *for* that intuition: the metric most likely
  to be gamed by a static-but-hash-churning attractor is exactly the
  raw-excess one that broke. With the §2 fix and the §3 U\* objective
  (which requires sustained diversity *and* turnover), a static
  attractor scores low by construction — the objective can no longer
  be satisfied by a frozen lattice that merely permutes silent LUT
  bits. So: adopt §3, then push pure-evo freely; the failure mode you
  fear is designed out.
- **"Optimise egenome-only / LUT-only / joint."** The harness is
  *almost* ready: `COG_METRICS` is defined in `evoca_explore.py:445`
  but **wired into nothing** — `evoca_from_scan_top` and the scan
  scripts never use it. Three line-level changes give you all three
  campaigns: (a) egene-only: scan with `mu_lut=0`, score by
  `COG_METRICS`; (b) LUT-only: `mu_egene=mu_egenome=0`, score by
  EVO_METRICS; (c) joint: score by a *stated* combination (do not
  unweighted-mean across incommensurable axes — see next point).
- **Multi-objective, not scalarised mean.** The current composite is
  an unweighted mean of min-max-normalised metrics — this implicitly
  weights by each metric's empirical range and hides trade-offs. For
  the joint egene+LUT campaign specifically, report the **Pareto
  front** over (U\*, cognitive-return) rather than a scalar; your own
  data says there is no trade-off in the *current* box, but the joint
  campaign deliberately leaves that box, so the trade-off may reappear
  and a scalar would mask it.
- **The determinism bug is a correctness issue, not a caveat.**
  `scan_analysis.md` observation + `evo_focus` "Caveats": the C
  xorshift RNG is per-process and never reset, so `multiprocessing.Pool`
  worker reuse makes per-config results non-reproducible; the 50
  baseline reruns gave different numbers than their originals. This
  means *no individual scan row is trustworthy*, only aggregates —
  which is acceptable for ranking but **not** for the
  finite-size-scaling and critical-exponent work (§5), which needs
  per-config reproducibility to do data collapse. Add `evoca_set_seed()`
  and call it at the top of `run_sim` *before* §5 or §9-large-world
  campaigns. This is a prerequisite, not an optional cleanup.
- **Run length.** All scans are 5 000 ticks. The §3 unbounded-activity
  diagnostic is a statement about t → ∞; 5 k ticks cannot distinguish
  bounded from unbounded (Bedau–Snyder–Packard explicitly require long
  windows). Before scaling N, scale t: a few configs at 50–100 k ticks
  will tell you more about open-endedness than a wide N-scan at 5 k.
- **Large-world budget.** Defer until §5 resolves H_ceiling vs
  H_critical. If H_ceiling, the 32-core/128 GB budget is best spent on
  long-t × multi-seed at N=256 (open-endedness is a temporal question
  at fixed size). If H_critical, spend it on the N ∈ {512, 1024} ×
  fine-axis sweep for exponent extraction. Either way, the budget
  decision is *downstream* of a cheap N=256 test — do not commit
  multi-day CPU before that gate.

## 10. Mutation of mutation

The cleanest, most literature-grounded direction in the document, and
mechanically well-scoped.

**C-side design:** add two per-cell float genes `mu_lut_gene`,
`mu_egene_gene` (8 bytes/cell; ~13% memory bump at the current 60
B/cell — acceptable). In Phase 4 the *child's* mutation rates are read
from the *parent's* genes; the rate genes themselves mutate
multiplicatively (log-normal step, so the rate stays positive and
scale-free) before being used. Global `mu_lut`/`mu_egene` become the
*initial* population values and an optional hard ceiling. This is
~40 lines in `evoca_step`'s reproduction block plus accessors; it does
not touch the CA core.

**Priors that make this sharp (not just "let's see"):**

- **Eigen's error threshold (1971):** there is a maximum mutation rate
  above which heredity collapses for a given genome length. EvoCA's
  LUT is 250 bits → a concrete predicted ceiling the population should
  evolve to stay below. *Prediction:* the population-mean `mu_lut_gene`
  converges to just under the error threshold for whatever effective
  genome length is under selection (which, via the restricted-mutation
  and dyn_activity machinery, is the *active* bit count, typically
  20–50, not 250 — so the predicted ceiling is higher than the naive
  1/250). This links directly to your dyn_activity probe and is a
  non-obvious, falsifiable quantitative prediction.
- **Sniegowski et al. (2000):** mutation rates evolve *down* under
  stable conditions (selection for fidelity) and *up* transiently
  under stress/maladaptation (hitchhiking of mutators with the
  beneficial variants they generate). EvoCA can test this cleanly
  because you control the environment: the prompt's "sudden food
  perturbation" is the canonical mutator-hitchhiking experiment.

**Protocol for your three questions:**
- *Does the distribution converge, and to what shape?* Long run,
  stable food; expect a right-skewed distribution pinned below the
  active-length error threshold (log-normal-ish, not delta — spatial
  structure maintains rate variance). The shape itself is the result.
- *Response to stress?* Equilibrate, then step `food_inc` down for a
  bounded interval and back. Predicted signature: a transient spike in
  population-mean rate during the maladapted interval, decaying back
  after — the mutator-hitchhiking curve. Measure spike amplitude vs
  perturbation depth.
- *Dependence on resource environment?* Cross the rate-gene experiment
  with the §6 food_inc sweep: predicted monotone decrease of evolved
  mean rate with resource abundance (richer environment ⇒ stronger
  fidelity selection ⇒ lower evolved μ).

## 11. Patch-transfer success assay

Status: **not implemented** (`egene_discussion.md` Item 2: "start when
ready"). It is a small, well-contained feature with a subtle
experimental-design core — the code is easy, the *controls* are where
it succeeds or fails as science.

**API (genome-level, not pixel-level — matches the recipe philosophy):**

```
sim.extract_patch(r0, c0, side) -> PatchGenomes
    # snapshot LUTs + egene slots + active masks for the square
    # region; positions normalised to patch-local coords.
sim.stamp_patch(patch, r0, c0, mode='overwrite')
    # write those genomes into the target sim's region; set alive,
    # default v/f. mode='overwrite' | 'invade' (only into dead cells).
```

**Scoring:** patch alive-fraction trajectory relative to t=0, *plus*
the boundary flux (births across the patch perimeter minus deaths) —
the perimeter flux is the leading indicator; total fraction is the lag
indicator. "Success" = sustained positive boundary flux, not merely
"patch didn't shrink in the first 100 ticks" (a patch can transiently
hold then collapse).

**Controls that make it mean something (this is the part the prompt
under-specifies):**

1. **Self-transfer control.** Transplant a patch of population P into a
   fresh field initialised from P's *own* recipe. This is the null:
   any growth/shrinkage here is boundary/initialisation artefact, not
   competitive superiority. *Every* invasion result must be reported
   *relative to* its self-transfer control or it is uninterpretable.
2. **Reciprocal transfer (mutual invasibility).** Adaptive-dynamics
   (Geritz et al. 1998) makes the right question precise: it is not
   "does A grow in B" but the 2×2 sign structure of {A-in-B, B-in-A}.
   Both grow → coexistence; one-way → competitive exclusion (genuine
   "evolutionary success"); neither → mutual resistance. Reporting
   only the one-way assay your prompt describes can mislabel
   coexistence as success. **Recommend the reciprocal 2×2 as the unit
   of measurement.**
3. **Size-series.** Run at patch side ∈ {N/6, N/4, N/3} — a critical
   patch size below which even a superior population fails to
   establish (a propagule-pressure threshold) is itself a key
   measurable and connects back to the §5 correlation length (a patch
   smaller than ξ cannot self-sustain its own spatial structure).

## 12. Dev logistics

| Direction | Branch? | Rationale |
|---|---|---|
| §2 neutral-shadow fix + fixed-space shadows | **main** | additive; everything downstream depends on it; no behaviour change to existing runs |
| §3 U\* objective + §9 metric wiring | **main** | analysis-side only (`evoca_explore.py`); old scans still reproducible |
| §1 lineage field | **main**, behind a compile flag / opt-in | core data addition many directions need; keep it zero-cost when off so it never slows the CA benchmark |
| §9 `evoca_set_seed()` | **main**, immediately | correctness fix; blocks §5 |
| §5 finite-size-scaling campaign | **main** | pure config + analysis |
| §8 ring-complexity tax | **branch** `ring-tax` | changes the step's tax computation → alters dynamics of every run; merge only if the genelife A/B validates it |
| §10 mutation-of-mutation | **branch** `mu-genome` | genome-structure change (per-cell rate genes); invasive; precedent is the existing `bigLUT` branch; merge after the error-threshold prediction is confirmed |
| §11 patch transfer | **branch** `patch-transfer` then merge | self-contained API; merge once the reciprocal-control protocol is validated, since the assay is generally useful |

General: keep the §1 lineage field opt-in and benchmarked — the
project's performance budget (N=256 ≈ 150 fps) is a stated constraint
and a parent-hash write per birth must not silently erode it. Branch
discipline: the two dynamics-altering changes (§8, §10) get branches
because a half-tuned tax or rate-gene scheme contaminates every run's
results; the analysis-side changes go straight to main because they
re-interpret existing data rather than changing it.

## 13. Prioritised roadmap

**Do first (this week, mostly zero or low code, high information):**

1. **Reciprocal egene/LUT isolation controls** (§7.3) — *zero code*,
   uses existing setters; bounds every later cognition claim.
2. **§2.1 per-component normalisation of excess activity** (~10 lines
   in `run_sim`) — unblocks every open-endedness measurement and
   un-breaks the metric your last scan optimised on.
3. **`evoca_set_seed()` + reset in `run_sim`** (§9) — correctness;
   prerequisite for §5.
4. **N=256 finite-size pre-test on food_inc** (§5 steps 1–3 at two
   sizes, N ∈ {128, 256}) — cheap; decides whether the multi-day
   large-world budget is worth spending at all.

**Do next (1–2 weeks, scoped code):**

5. Lineage field, opt-in + benchmarked (§1) — backbone for §4/§6/§7.
6. Fixed-space `dyn_activity` and `eg_activity` neutral baselines
   (§2) — cheap, rigorous, directly answers "is egene evolution real."
7. U\* objective + `COG_METRICS` wiring + long-t runs (§3, §9).

**Do after the above gate (branches):**

8. Ring-complexity tax + genelife A/B (§8, branch).
9. Mutation-of-mutation + error-threshold/stress protocol (§10,
   branch).
10. Patch-transfer assay with reciprocal 2×2 controls (§11, branch).

**Single highest-leverage item:** the lineage field (§1). It is the
common dependency of the landscape map (§4), the resource→complexity
trajectory (§6), the lineage-resolved cognition test (§7), and the
patch-transfer success criterion (§11). Everything else is sharper
once it exists; several questions are simply unanswerable without it.

**Single most urgent fix:** the §2 neutral-shadow defect — not because
it is the biggest feature but because it is currently producing a
*wrong optimisation conclusion* (the mu_lut=0.06 "optimum") that, left
unaddressed, would misdirect the multi-day campaigns the prompt wants
to scope.

*Sources for this section's literature anchors:
`Docs/references/references.md` (Taylor 2016 and Banzhaf 2016
downloaded full-text; Channon 2001 already in Docs; Bedau–Packard
series and paywalled landscape/criticality/mutation references cited
with publisher links there).*

# NP response

On reframing issues:

> "Is the population moving up the mountain?" needs a fitness proxy EvoCA does not yet expose. There is no explicit fitness in the model; fitness is realised, not assigned — it is the lineage's reproductive success and persistence. Every landscape and open-endedness question therefore reduces to lineage bookkeeping (who descends from whom, how many offspring, how long a genome's descendants persist). EvoCA currently tracks genome presence (activity) and population, but not parent→child lineage. The single highest-leverage instrumentation addition in this whole document is a per-cell parent-hash / lineage-ID field written in Phase 4. Almost every "are we making progress" question becomes measurable once that exists, and stays unanswerable (or only proxied) without it. I treat this as the backbone dependency in the roadmap.

I like the idea of focus on lineage.

>"Mt Fuji vs rugged" is a false binary. Ruggedness is a continuous, tunable property (Kauffman & Levin 1987; in NK terms, K is a knob), and in EvoCA the knob is the metaparameter vector itself plus the spatial structure of the food field. The right deliverable is not a classifier that stamps a run "Fuji" or "rugged"; it is a ruggedness coordinate computed per run and mapped across metaparameter space — same machinery as the existing scans, one more column.

NK model is built to have K be a tunable ruggedness parameter.  It is not very useful to say that the analog is "the metaparameter vector itself plus the spatial structure of the food field".  Yes: ruggedness is determined by these aspects of EvoCA, but not as directly as varying K in the NK model.  For EvoCA, ruggedness is more of an emergent property with no obvious a priori links to the metaparams.

>"Scaling with lattice size" and "critical exponents" are the same experiment. Critical exponents are defined through finite-size scaling. Your scan analysis already reports a correlation-length plateau at ≈ N/3 and calls it a probable "grid-size ceiling." The competing — and more interesting — hypothesis is that the system sits near a critical point where ξ is cut off by L, which is also exactly ξ ≈ N/const. These two readings are experimentally distinguishable and the test is one protocol (§5), not two sections.

No.  You are right that critical exponents may be associated with finite size scaling.  You are wrong that that is the only place they occur.  In physics, exponents are also measured for things like applied magnetic field for models of magnets (e.g. Ising model).  More relevant for us is critical behavior in bifurcations, changes in global dynamics as parameters are varied.  The most famous example is probably the period doubling bifurcation, with its famous Feigenbaum scaling.

# Claude reply to NP response (2026-05-16)

**R1 — lineage.** Noted; no disagreement. It stays the backbone
dependency in the roadmap (§13). Nothing to revise.

**R2 — ruggedness is emergent, not an NK-style knob. Conceded.** The
NK analogy was loose in exactly the way you identify: in NK, *K is
ruggedness by construction* — set a priori, before any dynamics. EvoCA
has no such dial. Ruggedness is an emergent property of the coupled
(CA rule ⊗ food field ⊗ selection ⊗ spatial structure) system with no
a priori functional tie to the metaparameters. My "the knob is the
metaparameter vector plus the food field" sentence dressed an emergent
property up as a control parameter, and — worse — I then smuggled the
same error back as a "falsifiable prediction" (ruggedness ↑ with
mu_lut, ↓ with food_inc). That prediction is **withdrawn** (marked in
§4). Consequences for the research path, which actually sharpen it:

- The metaparam → ruggedness map (§4 heatmap) remains the right
  deliverable, but it is **purely empirical and exploratory**: no
  expected monotone form, and no claim that ruggedness is even a
  smooth function of the metaparams. Whether it varies smoothly or
  jumps discontinuously is itself a primary result — and that is the
  bridge to R3 (a discontinuous jump in an emergent property as a
  metaparam is tuned *is* a bifurcation).
- Because the "landscape" in EvoCA is co-constructed by the population
  (neighbours and the food field both change as the population
  evolves), ruggedness is plausibly **non-stationary within a single
  run**. So R1 (lineage realised-fitness autocorrelation) must be
  measured in sliding windows and its drift tracked, not collapsed to
  one per-run scalar. The non-stationarity of ruggedness is a more
  faithful EvoCA observable than a static coordinate — this is a
  genuine refinement your objection forces, not a concession of
  scope.

**R3 — critical exponents are not only finite-size scaling.
Conceded; I over-collapsed §5.** Folding "scaling with lattice size"
and "critical exponents as metaparams vary" into one experiment was
wrong. Correct taxonomy — three *distinct* exponent programs, of which
your original prompt was asking mostly for the second and third:

- **(a) Finite-size scaling.** Exponents from L-dependence at a
  continuous transition. This is the *only* thing my original §5
  covered. Its legitimate, narrow use here is the specific test of
  whether the ξ ≈ N/3 plateau is criticality vs a periodic-box
  artefact (the H_critical/H_ceiling discriminator). Keep §5's
  Binder-cumulant/susceptibility protocol — but only for *that*
  question.
- **(b) Control-parameter (field-like) exponents at fixed L.** Order
  parameter ~ |p − p_c|^β as a metaparameter p is tuned through a
  transition — the Ising magnetisation-vs-temperature/field analogue.
  This is the literal reading of your prompt's "phase transitions as
  metaparams are varied … measure critical exponents." It needs no
  large worlds at all: fix N, sweep p finely (food_inc and mu_lut are
  the prime suspects, from the scan brackets), fit β, γ, δ from the
  order-parameter and susceptibility singularities at p_c.
- **(c) Bifurcations in the global dynamics — the one you flagged.**
  Changes in the *temporal* dynamics of global observables as a
  metaparam is varied: fixed point → limit cycle → period-doubling
  cascade → chaos, with Feigenbaum δ ≈ 4.6692 as the canonical
  scaling (and intermittency / quasiperiodic routes as alternatives).
  This is dynamical, not spatial-equilibrium; the relevant
  observables are EvoCA's *global time series* — `ts` probe's `pop`,
  `F_env`, `lut_div`, total/excess activity — not the spatial
  correlation length. The right instrument is a **bifurcation
  diagram**: slowly ramp one metaparam (or run a fine grid of fixed
  values), and for each value plot the long-time recurrent set of an
  observable (e.g. local maxima of the population time series after
  transient). Period-doubling shows as the classic pitchfork cascade;
  Feigenbaum δ is then extracted from the geometric convergence of
  successive doubling thresholds.

Practical EvoCA consequence: program (c) is mostly **analysis on swept
runs you can already produce** — `ProbeLogs/` already persists the
`ts` scalar series; a bifurcation-diagram harness is a driver that
sweeps one metaparam, runs long (post-transient), and reduces each
run's `pop`/activity series to its recurrent extrema. No C changes,
no large worlds. This makes (c) cheap and high-value, and it is the
program that most directly answers your phase-transition question —
so it should be promoted in the roadmap above the finite-size-scaling
work, which I had over-weighted. Revised priority: **(c) bifurcation
diagrams on the `ts` series ≻ (b) control-parameter exponents ≻ (a)
finite-size scaling** (the last only as the ξ-plateau discriminator).

Net effect on the roadmap (§13): R3 promotes a new low-code item —
"bifurcation-diagram harness over `ProbeLogs` `ts` series, sweeping
food_inc / mu_lut / tax, looking for period-doubling and Feigenbaum
scaling" — into the *Do next* tier, ahead of the large-world
finite-size campaign, which drops to conditional (gated on the
ξ-plateau test still being of interest after (c) characterises the
dynamical transitions).


# On parallelization and dev logistics

## NP

❯ Thanks for clarification and agreement on framing.  Your added comments are useful.

  Re proceeding with dev:

  I'm happy that you like my suggestion of 'causal control' (your 7.3).  Let's definitely pursue this.  Do you have
  enough specificity to launch runs on this right away?

  So your roadmap 2 fixes the bug you identified?  Seems like we want to do that right away.

  Sidebar:  I would like your introspection on how you found the bug noted in 2nd paragraph of your section 2.  I am
  constantly worried that the model is complex enough, and the dynamical behavior is complex enough so that it is
  difficult to tell if the code is bug free.  I would like to know:
  - does your discovery of this bug suggest any further checks we might implement to search for other bugs?
  - could we begin a test suite (using pytest or similar)?
     - check that GoL works as it should.
     - maybe identify some test cases where tax>XX => death tax<XX => life, and check these as a test.
     - other?

  Proceeding: my normal mode of interacting with you seems to have my personal interaction with you as a bottleneck.
   I would like your suggestions on parlaying your manegerial powers to have multiple threads of the research program
   going in parallel (e.g. everywhere you suggest a branch in the repo could be handled by independent agents).  I
  read that ex-Anthropic Boris Cherny says he has 'thousands of agents running overnight' with a fairly vanilla
  Claude code setup.  I am not suggesting that we push to thousands of agents, but I am suggesting we might move in
  the direction of running more agents in parallel.

  If I were to do things manually, I would perhaps start multiple claude code sessions.  But that would put me always
   in the critical path of progress, and I am trying to find ways to offload these critical paths to you.

  I await your suggestion here.

## Claude

⏺ The causal control already paid off — on its first run it surfaced a second instrumentation bug. Saving that
finding, then answering your questions.

Wrote 2 memories (ctrl+o to expand)

⏺ Everything below is done and verified; the campaign and tests have already run.

1. §7.3 causal control — launched, finished, and it already found a second bug

Specificity was sufficient: the reciprocal controls need zero new code — run_sim already seeds lut='gol' +
egenome='uniform' every run, so "freeze LUT" is just mu_lut=0 and "freeze egene" is just mu_egene=mu_egenome=0.
Driver: Scans/2026-05-16_causal_control/driver.py, 3 arms × 3 seeds at the egene_discussion V3 corner (N=128, 3000
ticks), each run in its own process (maxtasksperchild=1) so the per-process RNG is fresh — these results are
reproducible (variance across seeds is tiny, confirming it).

Results (means over 3 seeds, all non-extinct):


| arm | spec | load | intake | n_distinct | excess_pc_slope |
| ---| ---| ---| ---| --- | --- |
| JOINT |  16.2|  19.3|   0.050|        2842|           +0.003 | 
| LUT_only (egene frozen)|  25.0|  25.0|   0.054|        2056|           +0.002 | 
|  EGENE_only (LUT frozen)| | 20.6|  29.2|   0.061|         147 |            −25.9| 



Two findings:

- Substantive science: with the substrate frozen (EGENE_only), cognition evolves the highest food intake (0.061 vs
JOINT 0.050) and the most cognitive load — egene evolution genuinely improves feeding when the LUT isn't shifting
underneath it. This is the cleanest "cognition pays" signal yet, and it's only visible because the control isolates
it.
- A second bug, found by the control itself: EGENE_only's excess_pc_slope ≈ −26 (consistent across seeds) is not
"egene evolution fails." The neutral shadow in C/evoca.c carries LUT bytes only, but get_activity hashes the full
genome (LUT ‖ egene). So excess activity compares full-genome diversity against an LUT-only baseline — it is
undefined for any egene-driven run. This directly confirms §2's recommendation to build the fixed-space eg_activity
shadow, and it means the EGENE_only/JOINT excess numbers must not be interpreted until that exists. Saved to project
memory.

2. Roadmap-2 fix — done

python/evoca_explore.py now emits excess_pc_slope / excess_pc_final (component-normalised ΣG/D_G − ΣN/D_N, the
Channon-2006 correction) alongside the legacy raw metric, which is retained and flagged in-code rather than silently
dropped. One decision I deliberately did not make unilaterally: EVO_METRICS still lists excess_activity_slope.
Switching it to excess_pc_slope changes nearest_evo's behaviour and the meaning of any future "evo" scan — that's
your call. Recommendation: switch it, and re-rank the evo_focus scan with the corrected metric before any campaign
trusts the mu_lut=0.06 optimum. Say the word and I'll flip it.

Note this fix is necessary but not sufficient — finding #1 above shows the deeper shadow-scope bug is independent of
the normalization.

Sidebar: how I found the §2 bug, and what it implies

The reasoning path, honestly: I didn't find it by reading C line-by-line. I found it by a consistency cross-check
between two documents you'd written. scan_analysis.md reported excess_activity_slope falling while the other three
evo metrics rose, and labelled it a "known quirk." That phrase was the trigger: a metric moving opposite to every
correlated metric, explained away rather than explained, is a smell. I then traced the metric back to its definition
(run_sim), then to the shadow construction in C, and the bug was visible at evoca.c:1715 — the shadow uses the same
mutation operator as the real run, which stops being "neutral" at high μ. The generalizable move: a metric whose
anomaly has been rationalized in prose is a prime suspect; chase it to its lowest-level definition. The causal
control found bug #2 by a different mechanism — a designed discrepancy (freeze one subsystem, see if a metric does
something impossible).

That gives the two highest-yield bug-hunting strategies for a model this complex:

1. Differential/invariant tests (cheap, deterministic): freeze a subsystem and assert a metric cannot move; assert
conservation laws (food is transferred, not created — sum of f+F changes only by food_inc minus what tax destroys;
reproduction splits food 50/50; dead cells are fully zeroed). These catch the "metric measures something other than
you think" class — exactly bug #2.
2. Analytic-oracle tests (cheap, exact): GoL is a closed-form oracle; the tax deadline ceil(f0/tax) is closed-form.
Anywhere the model reduces to arithmetic you can compute by hand, assert the exact value.

Further checks worth adding (not yet done): a food-conservation invariant test (instrument total f+F per step,
assert the only sinks/sources are food_inc, tax, and the eating transfer — this would catch a whole class of
accounting bugs in Phases 2–4); a reproduction invariant (child genome == parent genome when all mutation rates are
0; food exactly halved); and a shadow-scope assertion (the shadow's diversity space must match the activity hash's
space — the bug #2 guard).

The test suite is started (tests/, runs with python3 -m pytest -q, 8 passing):
- test_gol.py — blinker period-2, block still-life, glider (1,1)-translation, and LUT-stability-under-zero-mutation.
The GoL oracles passing is strong evidence the CA core is correct.
- test_tax_survival.py — high-tax kills by the exact ceil(f0/tax) deadline and not before; zero-tax persists;
dead-cell zeroing invariant.
- test_metrics_regression.py — locks in the roadmap-2 fix so a refactor can't silently drop it (the exact failure
mode the original bug exploited).

One of my own tests failed first run — it called get_activity without evoca_activity_update(). That's instructive:
it's the same class of error as bug #2 (assuming a probe reflects state without the update that populates it), and
the test suite caught it immediately. The harness works.

Parallelization — proposal (awaiting your go-ahead before spawning anything)

The honest constraints first, because they shape the design:

- Scientific runs are CPU-bound and already parallel. Your scans saturate 12+ cores via multiprocessing.Pool. Ten
agents each launching a scan would contend for the same cores and run slower, not faster. Agent parallelism helps
independent code/analysis workstreams, not throughput of one CPU-bound campaign.
- The determinism bug (evoca_set_seed) is a hard prerequisite for any parallel scientific campaign — without it,
parallel runs aren't reproducible and the results don't compose. This gates the parallel-research vision; it should
be the very first parallelized task's output.
- Shared mutable repo state (the single C/libevoca.dylib, Runs/, Scans/) means concurrent agents editing/building
collide unless isolated.

Given that, the design that actually offloads you:

Tier A — parallel now, independent, isolated (safe to run many): the branch-tagged directions from §12 are
embarrassingly parallel because each is a separate branch with no shared files. Each runs as a background agent in
its own git worktree (isolation: "worktree"), so builds and edits never collide:
- ring-tax (§8) — C change + genelife A/B
- mu-genome (§10) — C change + error-threshold protocol
- patch-transfer (§11) — new API + reciprocal-control harness
- lineage-field (§1) — opt-in C field + benchmark

Plus analysis-only agents on main (read-mostly, low collision): the bifurcation-diagram harness (R3), the
eg_activity/dyn_activity fixed-space shadows (§2), the evoca_set_seed fix (do this one first and alone — everything
else depends on it).

Tier B — coordination so you're not the bus. The bottleneck isn't compute, it's integration decisions. The mechanism
that removes you from the critical path:
- A single tracked Docs/research_board.md — one row per workstream: branch, owner-agent, status, last result,
blocking-question, decision-needed. Each agent appends its result and any decision it needs; it does not block
waiting for you.
- A nightly orchestrator (via /schedule or a cron routine) that: pulls each worktree's status, runs the test suite
against each branch, writes a digest to the board, and only flags the specific decisions that need a human. You read
one digest in the morning instead of N agent transcripts.
- Branches never auto-merge. Agents propose; merge stays a human gate (the §8/§10 dynamics-altering changes
especially — a half-tuned tax contaminates everything).

Concretely, what I'd spawn if you approve, in dependency order:
1. Alone, first: one agent — evoca_set_seed + determinism fix on main. Everything waits on this.
2. Then, parallel worktree agents: lineage-field, ring-tax, patch-transfer, plus a main analysis agent for the
bifurcation harness and the eg/dyn shadows.
3. Nightly: orchestrator digest to research_board.md.

This keeps you out of the per-step loop while keeping you the merge authority. I have not spawned anything — the
system constraint is that I should only launch sub-agents when you explicitly direct it. Tell me to proceed (and
with which of the steps 1–3), and whether you want the nightly orchestrator on /schedule.

Nothing committed. The working tree now has: the research-doc edits + Docs/references/, the evoca_explore.py fix,
tests/, Scans/2026-05-16_causal_control/. Want me to commit these (logically grouped — doc, fix+tests, campaign)
now?

