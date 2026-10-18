# 10 — Drift shadow: what was found, why the EvoCA null fails, and what to do (2026-10-18)

For N. This is the long form of the drift-shadow status: what the statistics are (with
equations), what the validation showed, why the frozen-clonal marker null fails on EvoCA, and the
decisions, with recommendations.

**Sources.**
- Branch `drift-shadow` (`16f0807`), checked out for review at **`/tmp/evoca_drift`**.
- Its pre-merge notebook is `tests/notebooks/drift-shadow-test.ipynb`. It exists only on the
  branch, which is why it is not on `main`. Open it from `/tmp/evoca_drift/tests/notebooks/`.
  - It does not run the simulations itself. It reads `tests/notebooks/drift-shadow-validation.json`,
    written by `tools/drift_validation.py`.
- The design is in `Docs/shadows.md` §5–§8. The implementation report is §13 (on the branch).
- Code: `C/evoca.c` (`fs_drift_*`, `genome_mutate`, `marker_mutate`) and `python/drift_stats.py`.

Tags: **[measured]** = in the validation JSON or notebook outputs; **[code]** = read from the
implementation; **[reasoned]** = my argument, not yet tested.

---

## 0. Summary

**What the drift shadow is.** R replicate neutral populations, R = 8 by default.
- Each replicate individual carries a full EvoCA genome.
- Every step it receives the real run's death and birth counts. Victims and parents are chosen
  uniformly at random.
- Each newborn replays one real newborn's realized mutation counts, through the same mutation
  code the real run uses.
- It is a neutral population with the real demography's *totals* and the real mutation flux, but
  no selection and **no population structure**: it is panmictic, so anyone can be anyone's parent.

**What works.**
- On synthetic panmictic populations the statistics are calibrated (§3a) and have the designed
  power against imposed selection (§3b).
- Turning the shadow on does not change the real run, bit for bit.
- It costs +33% per step at R = 8.

**What does not.** On EvoCA the frozen-clonal neutral-marker null fails.
- Setup: a 12-bit "marker" trait that nothing in the dynamics reads, riding the real genealogy.
  With the marker the *only* thing that varies, the marker's concentration statistic should be
  ≈ 0 in units of the null sd.
- Observed: it comes out **−3.9 to −18 null-sd** (§4.1).
- **Cause, measured in part:** real surviving lineages accumulate **1.1–2.2×** the mutations per
  unit time that panmictic lineages do (§5).
- So the shadow's null is not EvoCA's null. Every evolving run reads negative, and no selection
  can yet be claimed from the `eg` or `dyn` drift statistics.

**My diagnosis goes one step past the branch report** (§5–§6).
- **The bias comes from EvoCA's genealogy, not from genotypes.** Children are placed on the
  neighbor with the least private food, typically a dead site with uneaten food. So the *child*
  edge of a birth tends to have higher reproductive value than the parent edge, and ancestral
  lines run through child edges (and hence mutations) more often than 1-in-2. That is spatial and
  ecological structure, and any genotype-blind panmictic null will misread it.
- **The branch's proposed fix,** an age-structured mirror, matches only part of this. It could
  flip the residual bias to the *anti-conservative* side (§6.2).

**Recommendation for the gating decision (D1).**
- Make the null **pedigree-conditioned**: neutral shadow genomes that ride the *real* parent→child
  events. Do not try to rebuild the demography synthetically.
- Score selection with an **expression-weighted substitution contrast**. This is the dN/dS idea:
  the pedigree shadow plays the role of synonymous sites (§7, Option P).
- It is exact under the null by construction, but not tautological in the way the tracking shadow
  was. The reason is that it copies genealogy, not allele content.
- It is also cheap to test first. EvoCA already has a free internal reference: LUT entries that
  are never queried are neutral sites on the same genealogy (§7, Option S).

**Decisions** (details in §8):

| # | Decision | Recommendation |
|---|---|---|
| D1 | How to repair the null | First run the cheap diagnostic (§9). Then build the **pedigree shadow + expression-weighted contrast** (Option P). Use the never-queried-input contrast (Option S) as a zero-cost preview now. Do not build the age-structured mirror as the main fix |
| D2 | Primary selection test | **No to the two-sided TV as primary** on any panmictic null; it detects demography. Under a pedigree null, use the expression-weighted contrast as primary and TV as a secondary any-departure check |
| D3 | `marker_mu` default 0.003 | **Yes.** Free, non-perturbing, and the marker null is the instrument that found this |
| D4 | Merge `drift-shadow` now? | **Yes, as infrastructure, with `run_sim` defaulting to `fs_shadow='tracking'`** (drift off; this reverses my earlier "'both'" recommendation), so scans don't pay +33% for an uncalibrated statistic. Offer `/code-review ultra` first |

---

## 1. Setup and notation

**Real run.**
- $n_t$ alive cells at step $t$.
- $B_t$ births and $D_t$ deaths in step $t$.
- Each real birth $b$ has realized mutation counts $(n_f, n_a, n_e, n_m)_b$: LUT flips, egene
  presence toggles, egene value/mask flips, marker flips.

**Shadow replicate $r = 1..R$.**
- A list of $n^{(r)}_t$ genomes.
- Per step: remove $D_t$ uniformly chosen individuals. Then add $B_t$ children, each a copy of a
  uniformly chosen survivor.
- Child $b$ is mutated with the real counts $(n_f, n_a, n_e, n_m)_b$ at freshly drawn positions
  (`genome_mutate` in replay mode) **[code]**.

**Channels and components $k$.**
- **eg:** the 729 ternary egene keys (value bits with wildcard mask). An individual contributes
  one count per active egene.
- **dyn:** the 500 components $(i, o)$, with LUT input $i \in \{0..249\}$ and output $o \in \{0,1\}$.
  - The real side uses the **whole-population** allele frequency
    $P_i = \#\{\text{alive cells with LUT bit } i = 1\}/n$.
  - It is weighted by the real query count $q_i(t)$ (how many cells used input $i$ this step).
  - Real per-step mass: $(q_i(1-P_i),\; q_i P_i)$. The shadow uses its own $P^{(r)}_i$ with the
    same $q_i(t)$ (`Docs/shadows.md` §5.4, design C).
- **marker:** a 12-bit trait, 6 value and 6 mask bits, mapped into the same 729-key space as eg.
  - It is inherited at real births, flipped $n_m \sim \text{Poisson}(12\,\mu_m)$ times from its own
    RNG stream, and **read by nothing** in the dynamics.
  - In the shadow, the marker is inherited and mutated like every other replicate gene (replayed
    $n_m$).

**Cumulative activity** of component $k$ from the enable step $t_0$:

$$A^{G}_k(t) = \sum_{s=t_0}^{t} m^{G}_k(s), \qquad A^{(r)}_k(t) = \sum_{s=t_0}^{t} m^{(r)}_k(s),$$

where $m_k(s)$ is the per-step mass of component $k$ (count of carriers, or the query-weighted
allele mass for dyn). $G$ is the real run. Both sides carry the same total mass per step, so
$\sum_k A^G_k = \sum_k A^{(r)}_k = M$.

---

## 2. The statistics

All statistics are computed for $R+1$ **units**: the real run ($u = 0$) and each replicate
($u = 1..R$). Each unit is scored against the other $R$ units (leave-one-out). Under the null the
$R+1$ units are exchangeable, so the rank p-value is exact **[code: `drift_stats._rank_stats`]**.

**2.1 Concentration** (eg, marker; secondary for dyn). With the Simpson index
$S(A) = \sum_k A_k^2 \big/ \big(\sum_k A_k\big)^2$:

$$c_u \;=\; \log S(A^{(u)}) \;-\; \frac{1}{R}\sum_{v \neq u} \log S(A^{(v)}).$$

$c_0 > 0$ means real activity is more concentrated on few keys than neutral diffusion from the
founder would leave it. The reported `*_drift_concentration` is $c_0$.

**2.2 Rank p and the "z".** These are the quantities behind the phrase "the marker's z":

$$p \;=\; \frac{1 + \#\{r \ge 1 : c_r \ge c_0\}}{R+1}, \qquad
z \;=\; \frac{c_0 - \bar c_{1..R}}{\operatorname{sd}(c_1,\dots,c_R)} .$$

- $z$ is the real unit's statistic in units of the replicate spread. Under a correct null it is
  approximately standard (sd ≈ 1; ≈ 1.25 for R = 8 because of $t_7$ tails).
- **"The marker's z"** is this $z$ computed from the *marker* channel's concentration
  (`marker_drift_z`). The marker is neutral by construction, so under a correct null
  $z_{\text{marker}} \sim (0, 1)$. A mean of −12 says the null is wrong. The marker is not the
  problem.

**2.3 dyn consistency** (primary dyn statistic). Per input $i$, the cumulative output
homozygosity and its standardized deviation are

$$F^{(u)}_i = \frac{(A^{(u)}_{i0})^2 + (A^{(u)}_{i1})^2}{(A^{(u)}_{i0} + A^{(u)}_{i1})^2}, \qquad
\zeta^{(u)}_i = \frac{F^{(u)}_i - \bar F^{(-u)}_i}{\max\!\big(\operatorname{sd}^{(-u)}_i,\,10^{-6},\,0.01(1-\bar F^{(-u)}_i)\big)} .$$

They are aggregated with the real query shares $w_i = q^{\text{cum}}_i / \sum_j q^{\text{cum}}_j$:

$$Z^{(u)}_{\text{dyn}} = \frac{\sum_i w_i\, \zeta^{(u)}_i}{\sqrt{\sum_i w_i^2}} .$$

$Z^{(0)}_{\text{dyn}} > 0$ means outputs at expressed inputs are held more consistently than drift
plus mutation would hold them (purifying selection).

**2.4 Two-sided TV** (added on the branch):

$$\mathrm{TV}_u = \tfrac12 \sum_k \Big| A^{(u)}_k - \tfrac1R \sum_{v\ne u} A^{(v)}_k \Big| \Big/ M ,$$

ranked like $c_u$. It is large for *any* difference between the real and neutral activity
distributions.

**2.5 Other keys.**
- $\Delta\text{Negene}$ (`negene_drift`) = time-mean real Negene − replicate mean.
- `*_nhi` = count of components with $z_k > 3$ minus its replicate expectation.

---

## 3. Synthetic validation (panmictic "real" population) [measured]

**(a) Calibration.** Np = 500, 25 births and deaths per step, 75 generations, 400 runs per cell,
mutation supply θ ∈ {0.5, 2, 10, 50} × input accessibility {uniform, geometric 0.7, 0.9}, R = 19.

| Criterion | Result |
|---|---|
| Rejection at α = 0.05 in [0.028, 0.072], primary statistics | 22/24 cells inside (0.085 and 0.075 outside; about 1.3 expected outside by chance) |
| Mean $z$ within ±0.10, sd in [0.90, 1.20] | max \|mean\| 0.08; sd 1.01–1.17 |
| Per-component $\lvert z_k\rvert > 3$ ≤ 1% | **fails**: up to 4.6% dyn, 5.2% eg (activities are heavy-tailed). The tracking shadow gave 23–35%. Use the calibrated `nhi` instead |

**(b) Power** (Np = 500, θ = 10, R = 19, α = 0.05):
- 0.85 on purifying selection at one locus with s = 0.05 over 150 generations;
- 0.79 with s = 0.02 at five loci;
- power rises with s and with generations.

**The sweep blind spot.** A directional sweep that completes inside the window reads *negative* on
the one-sided statistics (power 0.00–0.01). The reason, at one locus ($K = 2$ components):
- Before the sweep all activity sits on the old allele; after it, all of it on the new one.
- So $A^G = (\tau, T-\tau)$ for a sweep at time $\tau$, giving
  $S(A^G) = \big(\tau^2 + (T-\tau)^2\big)/T^2$. That is as low as ½ for a mid-window sweep.
- A neutral replicate that never moved has $S \approx 1$.

Cumulative activity forgets *when* mass sat where, so a completed sweep looks like
*diversity*. TV catches it (power 0.67–0.85). **[measured; the formula is mine]**

---

## 4. EvoCA results [measured]

N = 128 unless noted, 5,000 steps, R = 19.

### 4.1 The frozen-clonal marker null fails

**Base:** the pure-evo winner with all genome mutation off and a single clonal genome. Only the
marker varies.

| N, $\mu_m$ | mean $z_{\text{marker}}$ ± sd (20 seeds) | one-sided rejections | TV rejections | $\lvert c\rvert$ q95 |
|---|---|---|---|---|
| 64, 0.001 | −3.9 ± 2.0 | 0/20 | 75% | 0.30 |
| 64, 0.003 | −5.4 ± 1.8 | 0/20 | 95% | 0.67 |
| 64, 0.01 | −7.1 ± 2.6 | 0/20 | 95% | 1.30 |
| 128, 0.001 | −8.5 ± 3.8 | 0/20 | 100% | 0.32 |
| 128, 0.003 | −12.4 ± 3.1 | 0/20 | 100% | 0.67 |
| 128, 0.01 | −18.1 ± 4.0 | 0/20 | 100% | 1.38 |

Reading:
- The intensive bias $\lvert c\rvert$ is roughly independent of N.
- $z$ grows with N because the replicate sd shrinks, so $z \propto$ bias / sd.
- So this is a fixed misspecification, made more visible by power. It is not noise.

### 4.2 Lineage mutation load

The marker's Hamming weight approximates the number of marker mutations along the ancestral
line. Frozen-clonal winner base, $\mu_m = 5\cdot10^{-4}$, R = 8:

| N, seed | real | shadow (mean ± sd) | panmictic expectation | real / shadow |
|---|---|---|---|---|
| 64, 1 | 0.269 | 0.157 ± 0.055 | 0.12 | 1.7 |
| 64, 2 | 0.136 | 0.119 ± 0.035 | 0.11 | 1.1 |
| 128, 1 | 0.280 | 0.129 ± 0.013 | 0.12 | 2.2 |
| 128, 2 | 0.212 | 0.119 ± 0.025 | 0.12 | 1.8 |

**Mechanism of the sign.** More mutations along real lineages means each marker lineage has
wandered further through key space, so real marker activity is *less* concentrated than the
shadow's. Hence $c_0 < 0$ and a negative $z$.

### 4.3 Evolving runs (4 seeds each)

| Run | eg concentration $c_0$ | marker $c_0$ | eg − marker | $Z_{\text{dyn}}$ | dyn TV $z$ | $\Delta$Negene |
|---|---|---|---|---|---|---|
| JOINT | −0.45 to −1.31 | −1.40 to −2.01 | **+0.38 to +1.56** (4/4) | −12 to −14 | +25 to +39 | −0.77 to −0.86 |
| LUT_only | 0 (egenes frozen) | −2.3 to −3.0 | — | −11 to −13 | — | — |
| pure-evo winner | −1.29 to −1.91 | −1.38 to −2.16 | −0.39 to +0.69 (2/4) | −39 to −43 | — | −0.46 to −0.49 |
| frozen control | 0 | −0.52 to −0.65 | — | 0 | — | — |
| standing-variation freeze (N = 64, random LUTs, μ = 0) | — | — | — | > 0 in **20/20** | — | — |

Reading:
- **Everything evolving reads negative, which is the bias direction.** So nothing can be claimed
  from the one-sided statistics.
- **The standing-variation freeze (20/20 positive) is *not* a clean positive control.** With no
  mutation the load bias (§5, H1/H2) is switched off. But the structural bias H3 is not:
  - clumped, autocorrelated fecundity makes real allele frequencies drift *faster* than the
    panmictic $1/n$;
  - faster drift means higher homozygosity, which reads positive on $Z_{\text{dyn}}$ even with no
    selection.

  So 20/20 is consistent with selection sorting standing variation *and* with H3 alone. It shows
  the statistic can go positive, not that it measures selection. **[reasoned]**
- **$\Delta$Negene < 0 is a second positive control.** `tax_per_egene` imposes a known cost per
  active egene, and the shadow sees real Negene held *below* neutral.
  - The load bias should push in the other direction, because neutral toggling with
    dup-on-activate drifts Negene upward and real lineages accumulate more of it.
  - So the reading is conservative. This is the one EvoCA result I would quote. **[reasoned for
    the bias direction]**
- **"eg − marker > 0" is not evidence yet.** The marker's 12-bit kernel and its rate $\mu_m$ are
  not matched to the egene kernel (8 slots × 12 bits, presence toggles, duplication) or to the
  egene rates. So the two channels feel the load bias differently, and the difference has no
  known null.
- **There is no dyn floor** at all, because the marker lives in eg key space.

---

## 5. Why the panmictic null fails: lineage structure in EvoCA

### 5.1 The rate of mutations along an ancestral line

Follow the ancestral line of a randomly sampled present-day cell backward in time. Mutations enter
the line only at **child edges**: births in which the line passes from the newborn copy (the
mutated one) to its parent. The parent copy is unmutated. Let:
- $\lambda$ be the rate at which the line passes through birth events, either as the parent or as
  the child;
- $f_c$ be the fraction of those passages that use the child edge.

The expected number of mutation events on the line per step is

$$\rho \;=\; \lambda\, f_c\, \bar\nu ,\qquad \bar\nu = \text{mean mutations per child}.$$

**In the shadow** parents are uniform and deaths are uniform, so parent and child copies are
exchangeable after the birth. Then $f_c = \tfrac12$ and $\lambda = 2B/n$, giving
$\rho_{\text{shadow}} = \bar\nu\, B/n$. That is the "panmictic expectation" column of §4.2
(0.11–0.12), and the measured shadow values match it **[measured]**.

**In EvoCA** both factors can exceed their panmictic values. The measured ratio
$\rho_{\text{real}}/\rho_{\text{shadow}} = 1.1$–$2.2$ bounds their product.

### 5.2 Three EvoCA mechanisms (hypotheses, with predicted direction)

| | Mechanism | Factor | Effect on real concentration vs panmictic | Status |
|---|---|---|---|---|
| H1 | **Child-site advantage.** The child goes to the neighbor with the *lowest* private food: a dead site, or a starving cell, typically sitting on uneaten environmental food. The parent stays on the patch it has depleted. The child's line therefore has higher reproductive value | $f_c > \tfrac12$ | more load → **lower** concentration (negative bias) | [reasoned]; directly testable (§9) |
| H2 | **Short generation intervals along successful lines.** Reproduction needs $f$ to climb from ≈ 0.5 back to 1, which is fast in rich patches. Lines that persist run through rich patches, whose rich-patch generation interval is shorter than the population mean $n/B$ | $\lambda > 2B/n$ | more load → lower concentration (negative) | [reasoned]; testable |
| H3 | **Clumped, autocorrelated fecundity** (patches persist; neighbors share fate). This raises offspring-number variance and its correlation along lines, so $N_e < n$ and coalescence is faster | — | **higher** concentration (positive bias, anti-conservative) | [reasoned]; currently masked by H1/H2 |

The measured net bias is negative, so H1 and H2 dominate H3 *in the runs tested*. H3 does not go
away when H1 and H2 are corrected. That point drives the recommendation in §6.2.

---

## 6. What the null should condition on

### 6.1 Two possible estimands

- **$E_{\text{pan}}$:** "Does real evolution depart from a panmictic neutral population with the
  same totals and mutation flux?" This is what the drift shadow tests now.
  - Any population structure counts as a departure, including spatial, ecological and life-history
    structure that is independent of genotype.
  - In EvoCA that structure is intrinsic (food patches, local placement of the child). So
    $E_{\text{pan}}$ cannot isolate selection on genome content.
- **$E_{\text{ped}}$:** "Given the realized genealogy, does genome content depart from neutral
  mutation along that genealogy?"
  - This is the dN/dS question.
  - Genealogy-shaping forces that are not about genome content (spatial structure, life history,
    and even hitchhiking on selected lineages) are absorbed into the null.
  - What remains is selection's footprint *on the content itself*: which sites change, and how
    often relative to neutral sites on the same tree.

For EvoCA, **$E_{\text{ped}}$ is the defensible estimand.**

### 6.2 Why the age-structured mirror is not enough

The branch proposes mirroring the life history: shadow individuals carry a birth step, and each
real death or real parent is mirrored by a uniform draw within the same age class.
- It targets H2, and partly H1, since a young cell often sits on a fresh site.
- It does **not** reproduce:
  - the child-versus-parent edge asymmetry *within* an age class (H1 proper);
  - the spatial and temporal correlation of fecundity (H3).
- If it removes the negative H1/H2 bias but leaves the positive H3 bias, the corrected null
  becomes **anti-conservative**: neutral content would read as "selected".
- The marker null would catch that, but only after the build. A demography that can only be
  approximated synthetically leaves this risk permanently. Using the real pedigree removes it by
  construction.

---

## 7. Options for repairing the null

| Option | What it is | Null exact? | Detects | Cost |
|---|---|---|---|---|
| **A** (branch's (a)) | Age-structured mirror in the panmictic shadow | No: approximates H2, misses H1 and H3 | $E_{\text{pan}}$ minus age effects | moderate C change; open-ended validation |
| **B** (branch's (b)) | Marker "floor": compare channel statistics with a neutral marker on the real genealogy | Only if marker kernel and rates match the channel; currently they don't | eg − marker difference | a dyn-space marker is needed for dyn |
| **C** (branch's (c)) | Accept the bias; report only floor-referenced effect sizes | No | — | none; leaves ACT-D1 blocked |
| **P** (recommended) | **Pedigree shadow.** For each real cell, R neutral shadow genomes (full 250-bit LUT + egenome) that ride the real parent→child events. At each real birth the child's shadow genomes copy the parent's and get the child's real mutation counts at fresh random positions (`genome_mutate` replay). At each real death or eviction they die with the cell | **Yes**: under content-neutrality the real genome and each shadow genome are exchangeable given the pedigree | $E_{\text{ped}}$: expression-weighted substitution contrast (below); TV becomes meaningful | 49 B × R per cell (25 MB at N=256, R=8). It generalizes the existing per-cell marker. Reuses `genome_mutate` |
| **S** (preview, no C change) | **Internal synonymous-site contrast.** LUT inputs that are never (or rarely) queried are neutral sites on the *same* genealogy with the *same* per-bit mutation rate (without `restricted_mu`). Compare turnover at expressed inputs with turnover at unexpressed ones | Approximately: needs a clean "never queried" set; the query distribution evolves | purifying selection at expressed inputs; positive selection if expressed sites diverge *faster* | analysis only, using existing data (per-input query counts, allele frequencies) |

### 7.1 The statistic for Option P (the dN/dS analog)

For input $i$, let $\Delta^{(u)}_i$ be the **divergence** at site $i$ in unit $u$: the population
mean of $\mathbb{1}[\text{bit } i \neq \text{founder bit } i]$ at the end of the window (or the
number of fixed differences). Under content neutrality, $\Delta^{(0)}_i$ and the $\Delta^{(r)}_i$ are
exchangeable for every $i$.

Why divergence and not a turnover sum such as $\sum_t |P_i(t) - P_i(t-1)|$:
- Purifying selection removes deleterious mutants, but only after they have appeared at low
  frequency.
- A turnover sum counts those transient polymorphisms, so it dilutes the signal.
- dN/dS uses fixed differences for the same reason.

Define the expression-weighted contrast

$$K^{(u)} \;=\; \frac{\sum_i w_i\,\Delta^{(u)}_i}{\sum_i \Delta^{(u)}_i}\;\Big/\;\frac{\sum_i w_i}{250}\,,$$

the share of turnover that lands on expressed sites, relative to a uniform spread.
- **Purifying selection** at expressed inputs gives $K^{(0)} < K^{(r)}$: real changes avoid
  expressed sites.
- **Adaptive change** gives $K^{(0)} > K^{(r)}$: real changes concentrate there.
- It is ranked over the $R+1$ units exactly as in §2.2.

**What it can and cannot see.**
- It sees selection that differentiates *sites*. That is what "adaptive activity on the rule" should
  mean.
- It does not see a single sweep at a random site as different from neutral. Neither can any
  pedigree-conditioned test, and that is correct: one sweep is not distinguishable from one lucky
  neutral fixation.
- Sweeps that land preferentially on expressed sites, or repeat there, are seen.
- For eg the analog uses the egene structure: value bits that the current fiducial pattern
  actually matches versus bits masked out (wildcards).

---

## 8. Decisions, with recommendations

**D1. How to repair the null.**
- **Recommendation: Option P,** preceded by the §9 diagnostic and an Option S preview.
- Do not build Option A as the main fix. It is a moving target (H1 and H3 remain) and can become
  anti-conservative.
- Option P also unifies two workstreams:
  - the pedigree is the same object as LIN-D2's barcoded clades;
  - the H1 lineage program can use the same per-cell shadow genomes.
- **What it unblocks.** ACT-D1 becomes "two channels, each scored by a pedigree-conditioned contrast,
  Pareto across channels within egene-cost strata".

**D2. Primary selection test.**
- **Under the panmictic shadow, none.** TV rejects for demography alone: 75–100% on a neutral
  marker.
- **Under Option P:** the expression-weighted contrast $K$ is primary. TV and $Z_{\text{dyn}}$ are
  secondary descriptors.
- Keep the sweep blind spot of §3 in mind for any cumulative-activity statistic.

**D3. `marker_mu` default.**
- **Yes, 0.003.** It draws from its own stream and cannot change the real run (tested).
- It is the instrument that exposed the misspecification, so it should be on wherever drift
  statistics are reported.
- Under Option P it becomes the low-cost per-cell check that the pedigree mirror is wired
  correctly: real marker and shadow markers must be exchangeable.

**D4. Merge `drift-shadow` now?**
- **Yes, as infrastructure.** Parity is clean, the refactor (`genome_mutate`) is what Option P
  needs anyway, and the synthetic machinery (`drift_synth`, `drift_ref`) is reusable.
- **But change one default.** This reverses my earlier recommendation (Q3 in `Docs/shadows.md`
  §8, which you accepted) that `run_sim` default to `'both'`.
  - `'both'` was chosen when the drift statistics were expected to be calibrated.
  - They are not calibrated on EvoCA, so `'both'` makes every scan pay +33% per step for
    statistics nobody should read yet.
  - Default to `'tracking'` until a calibrated null exists, and opt in with `fs_shadow='both'`.
- This merge touches pluggable modes, adds a second caller of the mutation operator and puts an
  observer in the real step. That makes it a good moment for `/code-review ultra`, which only you
  can start.

---

## 9. Cheap diagnostic before building anything

This is about a day of analysis with existing code (lineage fields `lin_parent_id` and
`lin_birth_id`, plus the marker):

1. **Measure $f_c$ and $\lambda$ directly (H1 vs H2).**
   - Run frozen-clonal with `lineage=True`.
   - Reconstruct the ancestral lines of the final population.
   - Count births passed through, and the fraction via child edges.
   - Prediction if H1 dominates: $f_c$ well above ½ and $\lambda \approx 2B/n$.
   - Check $\lambda f_c / (B/n)$ against the measured 1.1–2.2.
2. **Measure H3:** the variance of offspring number per parent, and its correlation with the
   parent's own parent (the lineage autocorrelation of fecundity). This sizes the positive bias
   that any synthetic-demography fix would leave.
3. **Option S preview on the existing `pure_evo_joint` / causal-control configs:**
   - per-input divergence from the founder in the *real* genome;
   - split by cumulative query share (top quartile vs never-queried);
   - against the same split under the frozen-control seeds.

   If expressed inputs diverge measurably less (purifying) or more (adaptive) than unexpressed
   ones in evolving runs, there is a selection signal that needs no shadow at all.
   - **Circularity caveat.** Which inputs are expressed depends on the genome. A mutation that
     changes an output can change which inputs get queried. Use cumulative query share over the
     window, not the final share.

---

## 10. Caveats and what is not established

- **H1–H3 are reasoned mechanisms.** Only their net effect (the 1.1–2.2× load) is measured.
- **Option S's "never queried" set may be small and drifting.** With 250 inputs and evolving
  CA dynamics, many inputs are queried at least occasionally. Use query share as a continuous
  weight rather than a binary split.
- **Option P conditions on the realized genealogy, so it cannot detect selection that acts *only*
  through genealogy shape.** An example is selection on mutation rate (mu-genome), whose realized
  flux the shadow replays. This is the same caveat as in `Docs/shadows.md` §5.2, and it should be
  stated as part of the estimand.
- **The TV and per-component tail results** are calibrated only on panmictic synthetic data.
- **Food cap.** The validation ran on base `3598e08`, which includes the food-cap merge, so its
  EvoCA runs use `food_cap=1`. The recipes carry no `food_cap` key, so they load with the
  default. The frozen-clonal base (the pure-evo winner) barely changes under the cap (guide 09).
