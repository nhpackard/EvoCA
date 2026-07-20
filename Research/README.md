# EvoCA Research

Curated results layer of the EvoCA repo-as-artifact. This directory holds
the **findings** of the research program; the model spec, theory, and
campaign planning live in [`../Docs/`](../Docs), and raw run/scan data
lives in [`../Runs/`](../Runs), [`../Scans/`](../Scans), and
[`../ProbeLogs/`](../ProbeLogs).

## Structure

One **subdirectory per campaign**, `Research/<campaign>/`, each
**co-locating** its summary notebook with its narrative so a result and
its write-up never live apart:

```
Research/
  README.md                     <- this index / reading order
  <campaign>/
    <campaign>.ipynb            <- summary notebook: narrative + figures
    <campaign>.md               <- (optional) prose write-up
    README.md                   <- (optional) one-paragraph orientation
```

Conventions:

- **Notebooks are committed with their outputs** so GitHub renders the
  figures inline — a visitor sees the results without running anything.
  (The repo does not strip notebook outputs; if repo-wide stripping is
  ever added for noisy working notebooks, exempt this directory with a
  `Research/.gitattributes` containing `*.ipynb filter=` and
  `*.ipynb diff=`.)
- **Raw artifacts stay in `../Runs/`, `../Scans/`, `../ProbeLogs/`.**
  Campaign notebooks read them in place; a campaign may also copy the
  specific curated artifacts it relies on into its own subdir. Narrative
  (`.md`) never links to a notebook in another directory — that is the
  pattern this layout exists to avoid; a notebook reading raw data is
  fine.

## Standing rule

A campaign **is not done until its `Research/<campaign>/` digest notebook
exists** (pairs with the pre-merge test-notebook protocol). See the P1/P2
decisions in [`../Docs/checkpoint 3 July 2026.md`](../Docs/checkpoint%203%20July%202026.md).

## Drafting workflow

Digest notebooks are **drafted** by the `research-migrate` agent into
`Research/_drafts/<campaign>.ipynb` — a staging area, so the top level
stays clean and no per-campaign subdir is created until a result is
approved. Review a draft there; on **your approval** it graduates to its
per-campaign home `Research/<campaign>/<campaign>.ipynb`, and only then
does its link replace the `⏳` in the index below. **Entry of a notebook
link in the table is gated by your approval** — an existing draft never
auto-populates the table.

## Index (reading order)

Rows appear when a campaign **starts** (Campaign + Question filled); the
Headline and Notebook columns fill in when its digest notebook is built.
`⏳` = digest notebook not yet migrated (science may already be done).
Candidates below are drawn from campaign work to date, strongest /
most-publication-defensible first.

| Campaign | Question | Headline result | Notebook |
|----------|----------|-----------------|----------|
| neutral-model-methodology | Can excess-activity metrics separate adaptive evolution from drift/turnover? | Fixed-space eg(729)/dyn(500) shadows + reciprocal-freeze controls; two shadow-scope bugs caught | ⏳ |
| resource-driving-inverted-U | Does evolvability peak at intermediate resource, and is the true axis `food_inc/tax` (gradient steepness)? | _in progress_ | ⏳ |
| pure-evo-regime (#1/#2) | What regime maximizes open-ended activity under joint vs LUT-only evolution? | Whole-genome excess ≈ 0; dyn/eg excess strongly +ve; low `mu_lut` + high `m_scale` wins (mu_lut optimum conditional on egene co-evo) | ⏳ |
| coevolution-substrate (#3/#3c) | Does egene co-evolution amplify rule selection; is a freezing penalty real? | Egene-freeze reduces dyn-excess; #3c mostly GoL-substrate death confound + `dyn_excess_pc` frozen-rich artifact | ⏳ |
| RD-robustness (#4) | Is reaction-diffusion spatial structure robust under pure-evo optimization? | No — corrL washes out; static-colony cases are lowest-`dyn_excess` | ⏳ |
| viability-brackets (R1) | Where are the productive optima in (`m_scale`, `food_inc`)? | `m_scale` interior ≈2.5–3.5; `food_inc` high 0.013–0.018; U-shaped viability | ⏳ |
