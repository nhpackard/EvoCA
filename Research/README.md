# EvoCA Research

Curated results layer of the EvoCA repo-as-artifact. This directory holds
the **findings** of the research program *and* the raw run/scan data that
backs them, organized by workstream; the model spec, theory, and campaign
planning live in [`../Docs/`](../Docs), and probe logs in
[`../ProbeLogs/`](../ProbeLogs).

**Two layers, two files.** This README is the *results* index — the
curated reading order of finished findings. The companion
[`research_board.md`](research_board.md) is the *coordination* surface —
the living process layer: decisions-needed, in-flight workstream/branch
status, the reconciliation with the checkpoint roadmap, and the running
log. The board tracks what is *active*; it carries campaign **status**
(headline + link back into this index) rather than the authoritative
result, so a finding never lives in two places that can drift. Read the
board for "what's happening now," this index for "what we've found."

## Structure

One **subdirectory per workstream**, `Research/<Workstream>/`, holding
everything for that line of work: its lab notebook, its open questions,
its finished digests, and the raw recipes and scans. Workstreams:
`Activity-metrics`, `Resource-driving`, `Lineage-co-evolution`,
`Reaction-diffusion`, `Criticality-bifurcation`, `Ring-complexity`
(the list is defined once, in `WORKSTREAMS` in
[`../python/evoca_paths.py`](../python/evoca_paths.py)).

```
Research/
  README.md                     <- this index / reading order
  SCANS.md                      <- how scans work (drivers, harness, promotion)
  research_board.md             <- coordination surface
  artifacts.csv                 <- provenance map: old path -> new location
  <Workstream>/
    <Workstream>.ipynb          <- lab notebook, chronological, newest section at bottom
    ISSUES.md                   <- single source of open questions / next directions
    digests/                    <- finished campaign digest notebooks (+ unapproved drafts)
    Runs/<file>.evoca           <- recipes
    Scans/<date>_<name>/        <- scan.py, scan_config.json, results.csv, notes.md
    Scans/<date>_<name>.md      <- (optional) cross-scan analysis write-up
  _unsorted/Runs/               <- default home for new recipes with no workstream
```

`ISSUES.md` is the one place open questions live: the lab notebook
displays it and the board pulls from it, so it is never copied by hand.

Resolving locations in code: use [`python/evoca_paths.py`](../python/evoca_paths.py)
(`REPO_ROOT`, `WORKSTREAMS`, `runs_dir()`, `scans_dir()`, `all_runs()`,
`all_scan_dirs()`) rather than hardcoded paths. `import_run('<file>.evoca')`
finds a recipe by bare filename, and `evoca_from_scan('<date>_<name>')` and
the other scan helpers accept a bare scan name. A **new** recipe goes to
the workstream named in the `EVOCA_WORKSTREAM` environment variable, else
to `Research/_unsorted/Runs/` (sort it later with `git mv`).

Conventions:

- **Notebooks are committed with their outputs** so GitHub renders the
  figures inline — a visitor sees the results without running anything.
  (The repo does not strip notebook outputs; if repo-wide stripping is
  ever added for noisy working notebooks, exempt this directory with a
  `Research/.gitattributes` containing `*.ipynb filter=` and
  `*.ipynb diff=`.)
- **Raw artifacts live beside the workstream that owns them**
  (`Research/<Workstream>/Runs/`, `.../Scans/`); probe logs stay in
  `../ProbeLogs/`. Digest notebooks read them in place, locating the
  repo root by walking up to the directory containing `C/evoca.h`.
  Narrative (`.md`) never links to a notebook in another workstream; a
  notebook reading another workstream's raw data is fine (e.g.
  neutral-model-methodology reads scans owned by Lineage-co-evolution).

## Standing rule

A campaign **is not done until its `Research/<campaign>/` digest notebook
exists** (pairs with the pre-merge test-notebook protocol). See the P1/P2
decisions in [`../Docs/checkpoint 3 July 2026.md`](../Docs/checkpoint%203%20July%202026.md).

## Drafting workflow

Digest notebooks are **drafted** by the `research-migrate` agent directly
into `Research/<Workstream>/digests/<campaign>.ipynb`, flagged as drafts
(unapproved). Review a draft there; on **your approval** it is treated as
a finished digest and only then does its link replace the `⏳` in the
index below. **Entry of a notebook link in the table is gated by your
approval** — an existing draft never auto-populates the table, even
though it already sits in `digests/`.

## Index (reading order)

Rows appear when a campaign **starts** (Campaign + Question filled); the
Headline and Notebook columns fill in when its digest notebook is built.
`⏳` = digest notebook not yet migrated (science may already be done).
Candidates below are drawn from campaign work to date, strongest /
most-publication-defensible first.

| Campaign | Question | Headline result | Notebook |
|----------|----------|-----------------|----------|
| neutral-model-methodology | Can excess-activity metrics separate adaptive evolution from drift/turnover? | Fixed-space eg(729)/dyn(500) shadows + reciprocal-freeze controls; two shadow-scope bugs caught | ⏳ (draft at `Activity-metrics/digests/neutral-model-methodology.ipynb`) |
| resource-driving-inverted-U | Does evolvability peak at intermediate resource, and is the true axis `food_inc/tax` (gradient steepness)? | _in progress_ | ⏳ (draft at `Resource-driving/digests/resource-driving-inverted-U.ipynb`) |
| pure-evo-regime (#1/#2) | What regime maximizes open-ended activity under joint vs LUT-only evolution? | Whole-genome excess ≈ 0; dyn/eg excess strongly +ve; low `mu_lut` + high `m_scale` wins (mu_lut optimum conditional on egene co-evo) | ⏳ (draft at `Activity-metrics/digests/pure-evo-regime.ipynb`) |
| coevolution-substrate (#3/#3c) | Does egene co-evolution amplify rule selection; is a freezing penalty real? | Egene-freeze reduces dyn-excess; #3c mostly GoL-substrate death confound + `dyn_excess_pc` frozen-rich artifact | ⏳ (draft at `Lineage-co-evolution/digests/coevolution-substrate.ipynb`) |
| RD-robustness (#4) | Is reaction-diffusion spatial structure robust under pure-evo optimization? | No — corrL washes out; static-colony cases are lowest-`dyn_excess` | ⏳ (draft at `Reaction-diffusion/digests/RD-robustness.ipynb`) |
| viability-brackets (R1) | Where are the productive optima in (`m_scale`, `food_inc`)? | `m_scale` interior ≈2.5–3.5; `food_inc` high 0.013–0.018; U-shaped viability | ⏳ (draft at `Resource-driving/digests/viability-brackets.ipynb`) |
| reaction-diffusion-gallery | What do EvoCA's photogenic reaction-diffusion regimes look like, and which recipes produce them? | _gallery topic — curating example `.evoca` recipes; the most photogenic sims to date (distinct from #4's robustness result)_ | ⏳ |
