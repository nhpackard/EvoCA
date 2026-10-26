# On Activity

I write this as a reflection on the endeavor of measuring evolutionary activity.

Claude:  please read this as deep analysis

# the basic idea

Basic idea is to attach counters to elements of individuals within a population, elements that change over evolutionary time scales, typically during reproductive events.  Activity measures cumulative counts of particular elements (typically particular values of an element) within the entire population, as a measure of how successful that particular evolutionary innovation is.  If an element contributes to the propagation of individuals, its increased presence contributes to the increase of observed activity.

The idea of activity is that innovations in the population (new elements) that confer reproductive advantage are flagged by new activity waves.

# problems

# need for a neutral shadow model

Elements have zero contribution to propagation of individuals, but still accumulate activity, because they are not actively expelled from the population.  They can just ride along.  This is relatively well understood, and activity may be normalized by subtracting off activity counts of parallel neutral shadow models.  Subtleties include whether an activity threshold is computed globally, or element by element.

This issue may be put aside for the moment, unless it can be argued that it impacts the further issues below.

# problem of many elements leading to equivalent functionality.

This is the essential problem I would like to address.

Observationally, I commonly see a population generating tons of activity waves, indicating that many new elements are being created and are contributing to propagation of individuals, conferring the necessary advantage to create an observed new activity wave.  But if I look at the dynamics of the population, I see no discernable change.

This indicates that the 'innovations' causing the new waves generate the same (statistically, as filtered by my visual apparatus) basic dynamic functionality.  The evolution seems to reach a plateau of many many innovations that are 'neutral' in the sense that they lead to the same spatial dynamics.

This is te problem: we want activity to indicate significant evolutionary progress, but the population seems functionally and visually static (on average).

We have tried various tweaks on activity, attaching activity counters to different aspects of individuals (varying the definition of 'element').  but so far, no real joy.

What I want from activity that I'm not getting at the moment:  
- If I see no new discernable dynamics, I should see (new) activity go to zero.
- If I see ongoing new activity I should see new dynamics, as reflected in new domain structure, new textures, new fronts of advantageous subpopulations, etc.

What I want from you (Claude):
- an audit of all the variations of activity (definitions of elements and counts) we have in our current EvoCA implmentation.
- suggestions for new approaches to activity that could possibly yield the desires I articulate
- clear statements of limitations to the current activity concept that might prevent the desires I articulate.

---

# Claude's reflections (2026-10-23, deep analysis)

## Summary

- **Activity measures something real, but not what you want.**
  - Activity counts the *persistence* of heritable variants. Its neutral shadows test whether a variant persists *better than drift*.
  - Your desiderata need something else: whether a variant is *different in what it does* from what came before, judged at the level where you look, which is the lattice's spatial dynamics.
  - No choice of "element" fixes this as long as the count stays keyed on the variant's *identity*. The key has to become a distance in a *phenotype* or *dynamics* space, with an archive of what has been seen.
- **Why you see many waves and no change.** In EvoCA the many-to-one map from genotype to dynamics is extreme. I measured it in today's standard-rule pure-evo winner (N=96, t=3000):
  - A living cell has used a median of 105 of its 250 LUT entries over its lifetime.
  - The 25 most frequent inputs carry 70% of all transitions, and the top 50 carry 91%.
  - So most mutations change the rule only where the population almost never looks.
  - Lineages carrying such mutations still sweep, through drift or through small genuine fitness differences, and each sweep is an activity wave.
  - The waves are honest about *succession*, but blind to *function*.
- **The current variants each fix part of this, and none fixes all of it.**
  - `eff_activity` was meant to collapse silent churn. Its population-union mask saturates (248 of 250 entries here), so it collapses nothing.
  - `dyn_activity` lives in the right space, the transitions actually used. But it is a fixed, finite space (500 buckets) that saturates early, so it can register *concentration*, not *novelty*.
  - Local-pattern activity lives in phenotype space too, but it is not heritable and counts dead sites.
- **What I suggest, in priority order:**
  1. **Measure your observable directly.** Build a *dynamics-novelty* signal on population-level spatial descriptors (textures, domain sizes, fronts, periods), against an archive, with a frozen-genome control run as the null. That alone tests both your requirements.
  2. **Link it to heredity** with frequency-weighted *phenotypic* mutation classification at birth (cheap, in C) and a phenotype-distance key for genomes.
  3. **Count a wave as "new activity" only when it is both selected** (beats the shadow) **and novel** (outside the archive by more than a calibrated distance).
  4. **Validate** on a known evolutionary change of dynamics: under the standard rule, evolution damps the boom-bust waves to a fixed point around t≈4000–6000. The metric must flag that event and stay quiet in your plateau runs.
- **The limits are stated plainly in the "Limitations" section below.** The deepest one is that in EvoCA "dynamics" is an ecological, population-level property, while activity is an individual-level count. Your second requirement can only be met by measuring at the population level and *attributing* change back to lineages, not by any per-individual element.

## Outline

1. Your two requirements, stated as testable conditions.
2. Audit: every activity variant in the current implementation.
3. What the prior data already tell us.
4. Diagnosis: why "many waves, no visible change".
5. Limitations of the activity concept (clear statements).
6. New approaches: what could meet the requirements.
7. A validation design that uses events we already know.
8. Critique of my first draft, and the revisions it produced.
9. Recommended plan and costs.
10. Questions for you.

## 1. Your two requirements, stated as testable conditions

- **R1 (no false positives):** if the population's dynamics show no discernible novelty, new activity goes to zero.
- **R2 (no hollow positives):** sustained new activity implies new dynamics, visible as new domain structure, textures, or fronts of advantageous subpopulations.

Two remarks shape everything below.

- **These are conditions on a map from genotype events to phenotype change**, and the phenotype you judge is *population-level and spatial*, not a property of one cell. A cell's behavior is its rule applied to its neighbors, who carry other rules. Textures and fronts are ecological products.
- **"Discernible" needs a null and a scale.** Frozen genomes (`mu=0`) still produce fluctuating dynamics. So "no new dynamics" has to mean "no departure larger than the frozen-genome system's own variability", and the scale should be calibrated to your eye: you labeling a small set of runs is a cheap and legitimate calibration.

## 2. Audit: every activity variant in the current EvoCA implementation

| Variant | Element (key) | Key space | Counting | Neutral baseline | Status and known problems |
|---|---|---|---|---|---|
| `activity` (genome activity) | FNV-1a hash of LUT bytes plus sorted active egenes | open | cumulative presence per alive genome per update | Channon shadow population (random parent/death, same demography) | Silent churn: any bit flip mints a species. The shadow is LUT-only while the key is the full genome, so it is invalid for egene evolution ([ACT-Q1](../ISSUES.md)). Demoted to a churn diagnostic (`dyn_report`) |
| `excess_activity` (raw) | ΣG − ΣN over the above | — | pooled | Channon shadow | Inverts at high `mu_lut` (the shadow becomes a no-heredity model). Legacy; do not rank on it |
| `excess_pc` | ΣG/D_G − ΣN/D_N | — | per-component | Channon shadow | Turnover-invariant. About 0 in pure-evo runs (this is the "whole-genome excess ≈ 0" result) |
| `A_new` (Bedau–Packard) | count of genomes above the shadow's p90 | — | count | Channon shadow p90 | Turnover-invariant. Blind in egene runs (0 in EGENE_only while `eg` excess is positive). A frozen control gives 1, not 0 |
| `n_activity` | the shadow's own genome activity | open | as `activity` | — | Baseline only |
| `eff_activity` (F1) | LUT masked to the population-union used mask, plus egenes | open (coarsened) | as `activity` | none | The union saturates (248/250 here, 250/250 in the notebook's runs), so it equals full-genome activity. The per-cell mask splits identical genomes (3–10x more species). Collapses nothing ([ACT-Q2](../ISSUES.md)) |
| `eg_activity` | egene ternary key (value, mask) | fixed, 729 | cumulative per key | closed-form Wright–Fisher drift driven by realized egene flux | Valid for egenes. `eg_excess_pc` explodes when one key is present ([ACT-Q3](../ISSUES.md)) |
| `eg_food` | food intake per egene key | fixed, 729 | cumulative food | none | A *function*-weighted readout (what the egene earns), unexploited as activity |
| `dyn_activity` | (LUT input, output) transition actually fired by an alive cell | fixed, 500 | cumulative per bucket, per-step population | closed-form drift on realized LUT flux | The right space (used behavior). Finite, so support saturates early. `dyn_excess_pc` is artifact-prone at `mu_lut=0` (the #3c explosion turned out to be legacy-rule, [ACT-Q10](../ISSUES.md)) |
| z-scores (G2) | per bucket of `eg`/`dyn` | fixed | (obs − drift mean)/sd | Wright–Fisher | Well posed, exactly 0 under a freeze. Bidirectional (many z < −60), needs aggregation and an FP calibration ([ACT-Q9](../ISSUES.md)) |
| `dyn_report` (D1) | participation ratio `dyn_eff` over the 500 buckets, plus genome churn | — | level gate (60) | none | Measures *concentration*, not novelty. Its verdict on the winner flipped with the eating rule (silent-churn 38 legacy → functional 70 standard) and sits near its gate ([ACT-Q7](../ISSUES.md)) |
| Pattern activity and entropy (`pat_*`) | local v-configuration (VN/Moore/n3 window) at *every* site | fixed, ≤ 8192 | cumulative per pattern, Shannon entropy | none | The only *phenotype-side* (dynamics) activity. Not heritable, counts dead sites, no baseline, not in `run_sim` |
| Activity crossings (flux probe) | activity band crossings per genome | open | slopes in a band | — | A shape probe of waves; rarely used |
| Lineage net return (D2) | per-genome-lineage intake − tax | open | cumulative | none | A *fitness* readout per lineage, not activity. It can separate "adaptive sweep" from "drift sweep" |
| `unique_top_genomes`, `n_distinct_genomes` | identity of the dominant genome; diversity | open | counts | none | Turnover, not innovation. Saturate at the sample count |

Pattern: every *genotype-keyed* variant either churns (open key) or saturates (fixed key). The only *dynamics-keyed* variant (patterns) is not connected to heredity. That gap is exactly your problem.

## 3. What the prior data already tell us

- **Selection is real, but mostly invisible.** In the pure-evo scans (standard-model rescans of 2026-10-23), `dyn` and `eg` excess are strongly positive (median 22.6 and 16.1) while whole-genome excess is about 0 (mean 0.78). The *components* (transitions, egene keys) are under selection; the *genomes* turn over like drift. That is the many-equivalent-genomes picture, quantified.
- **The genotype-to-behavior map is extremely many-to-one** (the measurement above). About 105 of 250 entries are ever used by a cell, and 50 entries carry 91% of use. A uniformly random LUT bit flip lands on an entry carrying under 1% of the population's input mass 86% of the time, and on one with under 0.1% about two times in three (62–66%; same run).
- **The coarsening we built does not coarsen.** `eff_activity`'s union mask is 248–250 of 250 in every evolving run, because over a long run every input eventually occurs somewhere. Rare is not never, and a binary used/unused test cannot express "rarely".
- **Evolution in EvoCA tends to remove dynamical novelty, not add it.**
  - Campaign #4: reaction-diffusion spatial structure washes out under pure evolution (correlation-length final/mean median 0.23).
  - Today, under the standard rule, the boom-bust oscillation is damped to a fixed point by evolution (window std 0.22 → 0.005), while with mutation off it persists.
  - Your visual plateau may often be the *attractor* of selection, not a failure of measurement.
  - A metric meeting R1 would then rightly read near zero most of the time. Its value lies in catching the rare transitions, including transitions *toward* stasis like the damping.
- **Readouts flip with model details.** The `dyn_report` verdict and the #3c `dyn_excess` explosion both turned out to depend on the eating rule. Any new metric needs validation on the standard model, with frozen and shadow controls built in from the start.

## 4. Diagnosis: why "many waves, no visible change"

Four mechanisms, roughly in order of importance.

1. **Silent and near-silent mutations ride lineage succession.** In a saturated, high-turnover population every lineage eventually dies out, and its replacement is a variant. If the variant differs only on rarely used entries, its sweep is a new wave with an unchanged phenotype. The neutral shadow removes the *drift* part of this, but not:
2. **Adaptation without innovation.** Some sweeping variants are genuinely fitter at *the same dynamics*: slightly better eating, a slightly better response on a common input, slightly lower tax. That is real selection, so it beats the shadow and properly counts as adaptive activity. Visually nothing changes, because the improvement is quantitative and spread through the existing texture. Your R1 asks to *exclude* this, which is a legitimate choice, but it means the target quantity is not "adaptation" but "*phenotypic innovation that is selected*".
3. **Phenotype is ecological.** Textures and fronts are produced by *mixtures* of genotypes interacting through `v`. A new texture can arise from a shift in the *composition* of existing types (no new element), and a new element can sweep without changing the texture (its behavior matches the type it replaces in the contexts that occur). Individual-keyed counts cannot see the first case, and they over-count the second.
4. **Finite phenotype spaces saturate.** The phenotype-side keys we have (`dyn` 500 buckets, egene 729 keys, patterns ≤ 8192) fill early. After that, activity there measures *re-weighting*, not novelty. Novelty is then a property of the *distribution* over the space, not of any key's identity.

## 5. Limitations of the activity concept (clear statements)

- **L1. Identity-keyed activity cannot satisfy R1.** Whatever the element, a count of how long distinct keys persist will register succession among keys that are behaviorally equivalent, unless the key *is* the behavior class. Defining behavior classes requires a distance and a tolerance, which activity as defined does not have.
- **L2. The neutral shadow answers a different question.** It tests "fitter than drift", not "different from what was there". Making shadows better (fixed-space, per-bucket, z-scored) sharpens the adaptation signal, and cannot by itself produce an innovation signal.
- **L3. Open-ended novelty needs an open phenotype space and a metric on it.** Fixed phenotype spaces saturate (point 4 above). Open ones (genomes) are many-to-one. The way through is a *metric* phenotype space plus an *archive*: novelty means "far from everything seen", not "a new key".
- **L4. R2 is a population-level claim and cannot be met by per-individual elements alone.** EvoCA's dynamics are collective. The best an individual-level measure can do is be *attributed* to population-level change, and the change has to be measured where it happens.
- **L5. "Discernible" is observer-relative and noisy.** Dynamics fluctuate without evolution. R1 must be stated against a frozen-genome null, at a scale calibrated to human discrimination. Otherwise the metric either fires on noise or misses what you see.
- **L6. Statistical power is finite.** At N=96–256 many dynamical descriptors have large seed-to-seed spread (e.g. `eg` slopes 9 vs 21 across two seeds under the legacy rule). Novelty detection needs replicate controls, or it will mistake one realization's wander for innovation.
- **L7. Evolution here may mostly *remove* dynamical structure.** A metric that meets R1 and R2 will often read zero. That is a correct result about EvoCA, not a failure of the metric. If you want sustained innovation, the *model* may need a driver of open-endedness, such as a non-stationary environment, the chem food/waste coupling, or a richer resource space.

## 6. New approaches: what could meet the requirements

Ordered by how directly each serves R1/R2, with cost.

**A. Dynamics-novelty activity (primary; tests R1 and R2 directly).**
- **The descriptor.** Every K steps, compute a descriptor vector Φ(t) of the lattice: the local-pattern distribution (`pat_update`, n3 window), the spatial power spectrum (dominant wavelength, anisotropy), domain/patch size distribution, correlation length, alive-fraction period and amplitude, front speed (displacement of genotype-class boundaries), food/waste field statistics.
- **The novelty signal.** Keep an archive of past descriptors. Novelty(t) is the distance from Φ(t) to its k-nearest archived neighbors, normalized per feature.
- **The null.** Fork the run at time t into an evolving branch and a frozen-genome branch (`mu=0`; same state, same seed). The frozen branch's novelty is the dynamical noise floor.
- **The signal.** New dynamics activity = (evolving novelty − frozen-branch novelty), thresholded at δ.
- **How it meets the requirements.** R1 holds by construction (no departure beyond the frozen null gives zero). R2 is the definition.
- **Cost.** Python plus existing probes; forks double the compute, but periodic short forks suffice.

**B. Frequency-weighted phenotypic mutation classification at birth (cheap, at the source).**
- **What it computes.** In Phase 4, for each child mutation, the *input-frequency mass* of the changed LUT entries, weighted by the population's recent input distribution (`dyn_pop`, windowed). That is the fraction of encounters on which the child behaves differently from its parent.
- **What it gives.** Births carry a graded "phenotypic effect" e ∈ [0, 1]. Lineages inherit the cumulative effect since a reference ancestor.
- **The activity it enables.** Count activity only for lineages whose cumulative e exceeds ε, i.e. a "phenotypic species" keyed by the *behavior-relevant* part of the genome.
- **What it fixes.** This repairs `eff_activity` the right way (graded and windowed, not a binary union over all time) and directly removes near-silent churn (mechanism 1).
- **Cost.** Small C change in the reproduction loop, plus the lineage field already on main.

**C. Phenotype-distance key and archive for genomes (connects heredity to function).**
- **The distance.** For genomes g, h: d(g, h) = the input-frequency-weighted fraction of the population's current inputs on which their rules disagree, with the egene term weighted by intake differences, using `eg_food`.
- **The archive.** Cluster with tolerance ε into phenotype classes, and keep an archive of classes that ever reached significant activity.
- **The count.** **Novel adaptive activity** = waves of a class that (i) beats the fixed-space shadow (z-tail) and (ii) sits farther than δ from every archived class.
- **Caveat.** The input distribution moves as the population evolves, so d is context-dependent. Freeze the reference distribution per epoch and report sensitivity.
- **Cost.** Python, moderate; uses LUT hashes, per-cell masks and `dyn_pop`.

**D. Front-based sweep detection (your "fronts of advantageous subpopulations").**
- **The test.** For each rising phenotype class (from B or C), track its spatial geometry. A *front* is a compact domain whose radius grows roughly linearly; *diffuse spread* is a rise scattered across the lattice.
- **Why it matters.** Visible sweeps are fronts. Hitchhiking and near-neutral succession are diffuse. The ratio of front-like to diffuse sweeps is a direct, cheap observable of "visible advantage".
- **Cost.** Python, using the lineage field and connected components (the patch code exists).

**E. Knock-out validation of individual waves (causal, expensive, decisive).**
- **The test.** For a sampled wave, use the patch-transfer API (`extract_patch` / `stamp_patch`) to replace the wave's genome with its parent's in a forked run, and measure the change in Φ.
- **What it gives.** The causal phenotypic effect of that innovation. Use it to *validate* A–D, not as a running metric.

**F. What to keep from the current set.**
- **Keep** the fixed-space `eg`/`dyn` z-scores as the *selection* test (component (i) in C). Keep `A_new` and full-genome activity as churn diagnostics only.
- **Retire** `eff_activity` in favor of B, and `dyn_report`'s concentration gate as a novelty proxy.

## 7. A validation design that uses events we already know

- **Negative controls, which must read zero:**
  - frozen genomes;
  - the neutral shadow (shadow vs a second independent shadow, which also gives [ACT-Q9](../ISSUES.md)'s false-positive rate);
  - your plateau runs, the pure-evo winner after burn-in.
- **Positive controls, which must fire:**
  - **Boom-bust damping** under the standard rule. Mutation-on vs mutation-off is a known evolutionary change of dynamics around t≈4000–6000 (the Reaction-diffusion notebook, `rd_boom_variants`).
  - **The onset of the waste-balance patch regime** in the Chem workstream, when a parameter or evolutionary shift crosses the balance.
  - **A planted invasion:** stamp a known different-behavior genome patch into a plateau population, using E.
- **Calibrate δ two ways:**
  - **Frozen-branch spread:** set δ above it.
  - **Your eye:** you label about 20 run segments as "changed / not changed", and δ is set where the metric agrees with you best.
- **Report** per run: novelty rate, front fraction, selected-and-novel wave count, and the classic metrics alongside, so the new and old readings can be compared.

## 8. Critique of my first draft, and the revisions it produced

I drafted sections 4–6 first and then critiqued the draft as an opposing domain expert would. Five substantive flaws:

1. **I first treated "many waves, no change" as purely neutral equivalence.**
   - *The objection:* the prior data show real selection on components (`dyn`/`eg` excess strongly positive).
   - *Revision:* I added mechanism 2 (adaptation without innovation) and stated that your R1 deliberately excludes it. The target became "*selected AND novel*", and I added lineage net return (D2) as the way to tell adaptive sweeps from drift sweeps.
2. **I first proposed a genome "fingerprint in a standard environment" as the core novelty key.**
   - *The objection:* in EvoCA phenotype is context-dependent (rules act on neighbors carrying other rules), so a standard environment does not exist, and a moving reference makes the fingerprint drift and invent novelty.
   - *Revision:* population-level dynamics novelty (A) is now the primary judge of R1/R2, and genome-level keys (B, C) are demoted to *attribution*, with an epoch-frozen reference distribution and a sensitivity report.
3. **I first proposed a binary silent/phenotypic classification of mutations.**
   - *The objection:* that is exactly the binary that broke `eff_activity` (the union saturates because rare is not never), and the measured input-frequency distribution is strongly graded.
   - *Revision:* B is now a graded, frequency-weighted, windowed effect size with a threshold, not a used/unused flag.
4. **I first stated R1 as an absolute.**
   - *The objection:* dynamics fluctuate with frozen genomes, and "discernible" is observer-relative.
   - *Revision:* frozen-branch nulls built into A, and δ calibrated against your own labels (L5, section 7).
5. **One axis from your prompt was missing:** "fronts of advantageous subpopulations".
   - *Revision:* added D (front-based sweep detection), which is cheap and uses machinery we already have.

I also cut a re-explanation of neutral shadows, which you already know (you set that issue aside). I kept only the point that shadows answer a different question (L2).

## 9. Recommended plan and costs

| Step | What | Kind | Cost | Unblocks |
|---|---|---|---|---|
| 1 | Descriptor Φ and archive novelty in Python; frozen-branch fork helper | analysis code | ~1 day | A, the validation |
| 2 | Validation on the known events (boom-bust damping, planted invasion) and the negative controls; you label ~20 segments to calibrate δ | analysis | ~1 day + your hour | any claim |
| 3 | Graded phenotypic-effect classification at birth (C, Phase 4) with a lineage-inherited cumulative effect | dev (C + ctypes + tests) | ~1 day | B, the repaired eff species |
| 4 | Phenotype-distance classes, the selected-and-novel count, and front detection | analysis code | ~1–2 days | C, D |
| 5 | Re-rank the pure-evo and Chem runs on the new readouts, 3+ seeds | campaign (torque) | hours | [ACT-D1](../ISSUES.md) |

Steps 1–2 alone answer whether R1/R2 are achievable at all in EvoCA. I would do them before any C change.

## 10. Questions for you

1. **Do you want "adaptation without innovation" counted, excluded, or reported separately?** I recommend reporting it separately: it is real evolution, and seeing how much of the activity it accounts for is informative.
2. **Which dynamical features matter most to your eye?** Textures, domain sizes, fronts, oscillation and period, or something else. These set Φ's components and weights.
3. **Will you label ~20 run segments to calibrate δ?**
4. **If the honest answer is that EvoCA mostly evolves *toward* stasis (L7), is the research goal to measure that well, or to change the model to sustain innovation** (non-stationary environment, chem coupling)? This decides whether the activity work or the chem work is the main road.

# Further thoughts

## pattern activity

Just want to point out that local spatial pattern activity is actually a form of "dynamics-novelty signal on population-level spatial descriptors".  Doesn't refer to an archive as in your suggestion  (i)

## use flexibility in visualization

I see the point that adaptation can occur without changing visual perception of dynamics; you say
"
Adaptation without innovation. Some sweeping variants are genuinely fitter at the same dynamics: slightly better eating, a slightly better response on a common input, slightly lower tax. That is real selection, so it beats the shadow and properly counts as adaptive activity. Visually nothing changes, because the improvement is quantitative and spread through the existing texture. Your R1 asks to exclude this, which is a legitimate choice, but it means the target quantity is not "adaptation" but "phenotypic innovation that is selected".
"

This brings another perspective I forgot to mention:  visualization.

I ask for a link between adaptation as signalled by activity and change in visual dynamics.  You point out that real adaptation can indeed happen with no significant change in visual dynamics.  my discussion and cosequently yours was focused on trying to develop a different form of activity that would be linked to change of *current visuals*.  But another approach is to adapt the visuals to see current (valid) adaptations more sharply.  E.g. we could display eating efficacy, etc.

## Re. your questions:

1. **Do you want "adaptation without innovation" counted, excluded, or reported separately?** I recommend reporting it separately: it is real evolution, and seeing how much of the activity it accounts for is informative.

cf. comments on visualization.  adaptation without (current) visual effect is *not* to be excluded.  Reported separately (or with visualization choice) is fine

2. **Which dynamical features matter most to your eye?** Textures, domain sizes, fronts, oscillation and period, or something else. These set Φ's components and weights.

No preference.  anything visible would be good.

3. **Will you label ~20 run segments to calibrate δ?**

4. **If the honest answer is that EvoCA mostly evolves *toward* stasis (L7), is the research goal to measure that well, or to change the model to sustain innovation** (non-stationary environment, chem coupling)? This decides whether the activity work or the chem work is the main road.

I think this is the real point.  Vanilla EvoCA may simply evolve to dynamic/functional stasis. That's what it does.  I was trying to use mechanisms like patchy and boom regimes to give the population space to explore functional alternatives, but nothing drastically different seems to appear.

"whether the activity work or the chem work is the main road" is not how I would phrase the question.  For significant functional adaptation,  I think chem work is the direction to go.  But going in this direction is not exclusive of activity!  I want to see activity in that workflow as well, and hone the activity concept there to link activity signalled adaptation with visually perceived change.

## novelty detection

I resist going down the rabbit hole of exotic novelty detection using an archive.  I think if it is really the case that "EvoCA mostly evolves *toward* stasis" this would end up being a waste of time and effort.  We should just accept the result and innovate in the (meta) space of models and measurements, as in the chem direction.


# Claude's response to "Further thoughts" (2026-10-23)

## Summary

- **Pattern activity:** you are right, and my audit undersold it. Cumulative local-pattern activity *is* a dynamics-novelty signal on population-level spatial descriptors. Its pattern table is an archive already, a fixed and exhaustive one, so no learned archive is needed. I am revising the audit row (below). The one thing it still lacks is a null: a frozen-genome fork run gives one cheaply, with no new machinery.
- **Visualization:** adopted as a first-class route. Instead of forcing activity to track what the current display shows, add displays that show what activity is tracking. Concrete color modes are listed below. They fit your Q1 answer: adaptation without visible change is kept, and a color mode makes it visible.
- **Novelty archive:** dropped, as you ask. Proposals B and C (per-mutation phenotypic effect, phenotype-distance archive) and the δ-labeling exercise are dropped with it. Q3 needs no answer now.
- **Direction:** vanilla EvoCA evolves to dynamic/functional stasis. That is accepted as a result, not treated as a measurement failure. Functional adaptation is pursued in the chem model, and activity moves there with it. The aim there is the one you set: activity-signaled adaptation that you can see.
- **A design point this raises:** in the chem model as it stands, *no chem trait is heritable*. Excretion, toxicity threshold and metabolic product are all global. So activity in chem runs would still count LUT/egene variants, and adaptation to waste could only show up indirectly (where cells live, the evolvable tax terms). For activity to signal *chem* adaptation, something chem-related must be inherited. [CHEM-D3](../../Chem/ISSUES.md) (wgenes, evolvable excretion) is the existing item for this. I now think it should move ahead of "first characterize the no-wgene model". Your call.

## Revised audit row: pattern activity

| Variant | Element | Key space | What it counts | Baseline | Status |
|---|---|---|---|---|---|
| Pattern activity and entropy (`pat_*`) | local v-configuration (VN / Moore / n3 window) at every site | fixed, ≤ 8192: an exhaustive archive | cumulative presence per pattern, Shannon entropy | none yet; a frozen-genome fork is the natural null | **The archive-free form of proposal A.** Counting dead sites is correct for a descriptor of the lattice's dynamics. Not tied to heredity, and that is fine for R1. Not wired into `run_sim` |

What changes from my first reading: I listed "not heritable" and "counts dead sites" as defects. For a measure of what the *lattice* does, they are not defects. Heredity enters by comparison with the frozen fork: novelty in the evolving run minus novelty in the frozen fork is the part that evolution contributed.

## Color modes that make current adaptations visible

Existing modes: 0 LUT hash, 1 state + food, 2, 3 births, 4 age, 5 Moore neighborhood, 6 food/waste (chem). None shows *how well* a cell is doing. Candidates, roughly cheapest first:

| Mode | What it shows | What it makes visible | Data source and cost |
|---|---|---|---|
| Intake | food eaten this step, `M(x)`, or the match score that drives it | eating efficacy: a better egene/LUT combination lights up before it changes the texture | needs a per-cell intake buffer written in the eating loop (the lineage accumulators `lin_intake_cell` already do this when lineage is on) |
| Energy | private food `f(x)` for alive cells | who is near reproduction and who is near tax death; boom-bust phases | `f_priv`, no new state |
| Reproduction tempo | steps since last birth or reproduction | the fast-cycling subpopulation; fronts of advantageous lineages | `last_event_step`, no new state |
| Lineage growth | each cell colored by its genome's recent population growth rate | your "fronts of advantageous subpopulations" directly; this is proposal D as a display | per-hash population over a short window, from the activity table |
| Waste exposure (chem) | local `W` against the toxicity threshold, per alive cell | which cells live at the toxic margin; any later waste-tolerance adaptation | `chem_env`, no new state |

Cost per mode: a `case` in `evoca_colorize` and an entry in `colormodes.COLOR_MODES`. The intake and lineage-growth modes each need one buffer. Adding modes changes the colormode count again, so this should land after `chem` merges (mode 6 is on that branch).

## Revised plan

| Item | Status | Cost |
|---|---|---|
| A. Novelty on population descriptors | **Kept in archive-free form:** pattern activity plus a frozen-fork null, wired into `run_sim` | small |
| B. Graded per-mutation phenotypic effect | dropped | — |
| C. Phenotype-distance archive | dropped | — |
| D. Front detection | kept, as the lineage-growth color mode first; a quantitative front count only if the display shows fronts worth counting | small |
| E. Knock-out validation | deferred; only relevant if a new measure is proposed | — |
| F. Keep/retire guidance | stands as written | — |
| δ calibration labeling | dropped (Q3 withdrawn) | — |
| Adaptation color modes | **new, first:** intake, energy, reproduction tempo, lineage growth, waste exposure | small–medium |
| Activity in chem runs | **new:** the standard activity probes (`A_new`, `dyn_excess`, pattern activity) run in chem; [CHEM-Q2](../../Chem/ISSUES.md) already asks the main question | small |
| Heritable chem trait ([CHEM-D3](../../Chem/ISSUES.md), wgenes) | **raised in priority**, pending your decision | medium (C mutation and hash changes) |
