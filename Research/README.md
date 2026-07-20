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
decisions in [`../Docs/# checkpoint 3 July 2026.md`](../Docs/#%20checkpoint%203%20July%202026.md).

## Index (reading order)

_No campaigns migrated yet. Add a row per campaign as its digest lands._

| Campaign | Question | Headline result | Notebook |
|----------|----------|-----------------|----------|
| _(none yet)_ | | | |
