# Artifact census: Runs, Scans, ProbeLogs, drafts, root notebooks (2026-10-15)

Read-only census for N. Companion data: `Research/artifacts.csv` (120 rows; columns
old_path, owner_workstream, new_path, used_by, kind, confidence, note). Nothing was moved,
added, or committed.

## 1. Method

- **Runs (87):** parsed each `.evoca` JSON for `descriptor`, `source_scan`, `source_config_idx`.
  83 have a `source_scan`; 4 hand-made recipes have none (2026-03-16, 04-25, 04-26,
  04-27_R-D...). File-name class (`top_N` balanced composite, `top_evo_N` EVO-only, `auto_N`)
  was cross-checked against the ranking rules in `Scans/2026-04-27_scan_analysis.md` and the
  filter in `Scans/2026-05-18_rd_under_pure_evo/scan.py:38-42` (which globs `top_*`, skipping
  `top_evo`/`auto`).
- **Scans:** read each `notes.md`, each `scan.py`/`driver.py` docstring, the two top-level
  analysis .md files, `Scans/README.md`, and `Docs/devlog/2026-05-17_research-campaigns.md`
  section headings; assigned by the question each scan answers (review §3a table).
- **Drafts:** `git ls-tree -r research-drafts Research/_drafts` (6 notebooks, commit 84fc29b);
  assigned by the scan each one reads.
- **ProbeLogs (193 files, untracked):** grouped by descriptor suffix. Only 6 are attributable
  to a recipe.
- **Policy applied:** exactly one owner per artifact; `used_by` lists other workstreams that
  read or compare against it. Where a scan and its Runs disagree, the Runs follow the metric
  used to select them, not the scan that produced them (see 3.1).
- **Not covered by the data (no owner):** Criticality-bifurcation and Ring-complexity own
  **zero** existing artifacts. `python/bifurcation.py` and `tests/test_bifurcation.py` belong
  to the former, but they are code, not Runs/Scans. Ring-stability results (memory note) are
  not in the repo as Scans. Their dirs will start with `ISSUES.md` and a notebook only.

## 2. Per-workstream summary

| Owner | Runs | Scans (dirs) | Scan .md | Probe-log groups | Draft digests | Root nb | Total rows |
|---|---|---|---|---|---|---|---|
| Activity-metrics | 47 | 5 | 1 | 1 (4 files) | 2 | 0 | 56 |
| Reaction-diffusion | 40 | 4 | 1 | 1 (2 files) | 1 | 0 | 47 |
| Lineage-co-evolution | 0 | 3 | 0 | 0 | 1 | 0 | 4 |
| Resource-driving | 0 | 1 | 0 | 0 | 2 | 0 | 3 |
| Criticality-bifurcation | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Ring-complexity | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| shared-infra | 0 | 0 | 1 (Scans/README.md) | 0 | 0 | 4 | 5 |
| archive | 0 | 0 | 0 | 1 (~187 files) | 0 | 4 | 5 |

Scans by owner: AM = causal_control (v1), causal_control_v2, pure_evo_joint,
pure_evo_lut_only, evo_focus. RD = initial, refined, push_edges, rd_under_pure_evo.
LC = egene_on_winners (#3), 3c_coevolution_substrate, 8_lineage_coevolution.
RS = R1_mscale_bracket.

Confidence over all 120 rows: high 73, medium 42, low 5.

## 3. Assignments needing N's decision (medium and low)

### 3.1 `top_evo` Runs: the largest ambiguity (25 Runs, medium)
`Runs/2026-04-27_top_evo_{1..20}` (25 files; ranks 1-5 appear twice, once per source scan) have
`source_scan` = `initial` or `push_edges` (both RD scans) but were ranked by the EVO-only
composite. I assigned them to **Activity-metrics** (they exist to showcase activity, and the
04-29 sets from `evo_focus` go there with high confidence), with RD in `used_by`. Alternative:
keep each Run with its source scan's workstream (RD). Costs: the Activity-metrics choice splits
a scan from the Runs it produced, which weakens the "co-location by provenance" argument of §5
of the review. Note the legacy `excess_activity_slope` ranking behind these is flagged invalid
for egene-driven runs. Suggest: Activity-metrics.

### 3.2 Hand-made / duplicate Runs (7)
| File | Proposed | Why uncertain |
|---|---|---|
| `2026-03-16_2_branch_activity` (low) | Activity-metrics | No scan; N=512, gdiff=0, food_inc=0.12, activity probe on. Pre-dates every scan. Could be `archive`. |
| `2026-04-27_auto_test` (medium) | RD | Same cfg34 as `top_4_cfg34`; looks like a helper test export. Archive candidate. |
| `2026-04-27_auto_{1,2,3}` (medium) | RD | Duplicate content of `top_N` of the same cfg (49, 30, 70) in `initial`; harmless dupes. |
| `2026-04-27_high_excess_activity_slope` (medium) | Activity-metrics | Is cfg98 of push_edges, identical to RD `top_9_cfg98`; named for the activity metric. Possibly same content as top_9; diff before dropping one. |
| `2026-05-17_best_corr_length`, `2026-05-18_best_corr_length` (medium) | RD | All three `best_corr_length` files point to initial cfg49. Re-exports? Probably duplicates; ProbeLogs refer to the 05-17/05-18 file names, so do not delete before updating logs. |

### 3.3 Scans (4 medium)
- **`2026-04-27_evo_focus` (AM, medium):** review §3a assigns it to Activity-metrics, but its
  base configs come from `push_edges` (RD) via a relative path (scan.py:32-33) and its
  ranking used the legacy metric. If it moves away from `push_edges`, that path must be
  rewritten (see section 4).
- **`2026-05-17_pure_evo_egene_on_winners` (LC, medium):** #3 asks the LUT-vs-egene interaction
  question (Lineage), but is a pure-evo composite scan reading the Activity-owned
  `pure_evo_lut_only` results. Has no `notes.md` or `scan_config.json`.
- **`2026-05-18_R1_mscale_bracket` (RS, medium):** it brackets `m_scale` and `food_inc` (the
  resource axis), as the review says, but it is a mirror of #1 (Activity). Consistent with
  the review table; flagged only because its notes.md is mainly about `mu_lut x m_scale`.
- **`2026-05-16_causal_control` (AM, high but note):** no `notes.md`; only
  `results.csv` + `driver.py`. Superseded by v2; archive-within-Activity candidate.
- `2026-05-19_8_lineage_coevolution`: only an untracked `scan.py` (+ `__pycache__`); on hold.

### 3.4 Scan .md files (both medium)
- `2026-04-27_scan_analysis.md`: RD (3 of 4 scans are RD, composite-score method) but also
  covers evo_focus and the metric definitions; content is cross-cutting.
- `2026-05-17_pure_evo_analysis.md`: Activity-metrics (#1, #2, causal v2) with pointers to #3
  (LC), #4 (RD), R1 (RS).
- `Scans/README.md`: shared-infra, stays; needs path update.

### 3.5 Draft digests (2 medium)
`pure-evo-regime` (AM; the same #1/#2 result also feeds Lineage's #3 question) and
`viability-brackets` (RS; follows the R1 scan). `resource-driving-inverted-U` is a skeleton
with no data: RS, high.

### 3.6 ProbeLogs (1 low group)
- 4 files `ts_2026-04-27_top_evo_20_cfg131...` (AM, high): produced by `Untitled.ipynb`.
- 2 files `ts_2026-05-1{7,8}_best_corr_length` (RD, high).
- **~187 `*_{ts,egene,egenome}_manual.csv`** (05-10 to 05-18, descriptor `manual`): ad hoc
  interactive sessions, no recipe, not attributable per file. Proposed: archive wholesale
  (`Research/_archive/ProbeLogs/manual_2026-05/`), low confidence. They are untracked in git, so a
  move is plain `mv`. Alternative: delete or leave in a top-level `ProbeLogs/`, which
  `python/controls.py:716` keeps writing to anyway (it should keep its location: it is an
  output sink, not a research artifact).

### 3.7 Root notebooks
| Notebook | Class | Confidence | Note |
|---|---|---|---|
| `evoca_explore.ipynb`, `evoca_test.ipynb` | shared-infra, stay | high | harness tutorial; tests/benchmarks |
| `evoca_scans.ipynb`, `evoca.ipynb` | shared-infra, stay | medium | scan-ranking scratch across RD/causal/egene scans (25 path hits); early general notebook |
| `lut-hashes-accessor-test.ipynb` | archive | medium | test nb for merged branch |
| `test.ipynb`, `Untitled.ipynb`, `Untitled1.ipynb` | archive | low | untracked/scratch |

## 4. Path coupling that a move breaks

No shell scripts and no torque launch/pull scripts exist in the repo (searched `*.sh`,
`rsync|scp|ssh` in .py/.md/Makefile). The torque side lives outside git, so **N's torque
launch/pull commands (which assume `Scans/<name>/` at the top level) need a separate check.**
`tests/` has zero references to Runs/Scans/ProbeLogs. `python/bifurcation.py`,
`python/patch_transfer.py` have none.

### 4.1 Harness (`python/`)
| file:line | What | Break mode |
|---|---|---|
| `python/evoca_explore.py:379-389` `_default_scan_dir()` | newest-mtime subdir of `<repo>/Scans` | will find nothing / raise once Scans/ is split; mtime rule is also fragile |
| `python/evoca_explore.py:392-393` `_resolve_scan_dir()` | `None` -> default; otherwise path passed through unchanged | no legacy-name lookup: this is where the `artifacts.csv` resolver belongs |
| `python/evoca_explore.py:396-397` `_load_scan()` | reads `<scan_dir>/scan_config.json`, `results.csv` | follows the resolved dir |
| `python/evoca_explore.py:423-436` `evoca_from_scan()` | resolves scan dir then loads | same |
| `python/evoca_explore.py:456-458` | default `runs_dir` = `<repo>/Runs` | exports land in old location |
| `python/evoca_explore.py:472` | writes `source_scan` = basename of scan dir (bare name) | only a writer; **nothing in the repo reads `source_scan` back**, so legacy names in the 83 recipes are provenance only. A resolver is needed just if a future tool reads them |
| `python/evoca_explore.py:575`, `600-602`, `644-645`, `735-736`, `769` | `nearest_*`, `evoca_from_scan_top` call `_resolve_scan_dir` / `evoca_from_scan` | route through the resolver |
| `python/evoca_explore.py:8,431,600,732,491-494` | docstrings/comments naming `Scans/`, `Runs/` | doc update |
| `python/evoca_py.py:1438` `export_recipe()` | writes to `<repo>/Runs` | needs a `workstream=`/`runs_dir` arg |
| `python/evoca_py.py:1672,1684,1692-1698` `import_run()` | no-arg listing globs `<repo>/Runs/*.evoca`; explicit paths pass through | listing breaks (would need to glob `Research/*/Runs/`) |
| `python/controls.py:716` | `ProbeLogs` sink = `<repo>/ProbeLogs` | works if the dir stays; filename embeds `_origin_recipe` basename (renaming recipes decouples log names) |

### 4.2 Scan drivers (`Scans/*/`)
Every `scan.py`/`driver.py` computes `PYTHON_DIR = HERE/../../python` (depth 2). **After the
move (depth 4) all 13 break**, even ones with no other coupling:
`2026-04-27_{initial,refined,push_edges,evo_focus}/scan.py:17/17/17/24`,
`2026-05-17_{pure_evo_joint,pure_evo_lut_only,pure_evo_egene_on_winners}/scan.py:32/26/31`,
`2026-05-18_{R1_mscale_bracket,rd_under_pure_evo}/scan.py:31/25` (ROOT),
`2026-05-19_{3c_coevolution_substrate,8_lineage_coevolution}/scan.py:39-40/31-32` (ROOT),
`2026-05-16_causal_control/driver.py` and `2026-05-17_causal_control_v2/driver.py` (`_HERE/../../python` near line 20).
Fix once by a shared `find_repo_root()` helper (as the draft notebooks already do) instead
of a depth count.

Cross-workstream sibling paths (break differently, depend on the final owner split):
- `Scans/2026-04-27_evo_focus/scan.py:32-33` `SOURCE_SCAN = ../2026-04-27_push_edges`
  (AM -> RD).
- `Scans/2026-05-17_pure_evo_egene_on_winners/scan.py:37` reads
  `../2026-05-17_pure_evo_lut_only/results.csv` (LC -> AM).
- `Scans/2026-05-18_rd_under_pure_evo/scan.py:38-39` globs `ROOT/Runs/2026-04-27_top_*.evoca`
  (will become `Research/Reaction-diffusion/Runs/...`; the `top_evo`/`auto` exclusion at
  line 40 must be re-expressed against the new names `*_top[0-9]*`, `*_evo*`, `*_auto*`).

### 4.3 Notebooks
- Root: `evoca_explore.ipynb` cells 31, 48 (`Runs/2026-03-16_2_branch_activity.evoca`,
  `Runs/2026-04-25_waves_in_blank_half.evoca`); `evoca.ipynb` cell 27 (same 03-16 run);
  `evoca_scans.ipynb` cell 7 L1-7 (`Scans/2026-04-27_initial/results.csv`,
  `Scans/2026-05-16_causal_control/results.csv`,
  `Scans/2026-05-17_pure_evo_egene_on_winners/`,
  `Runs/2026-04-27_top_evo_20_cfg131_score0.51.evoca`), cell 12 (`Scans/2026-04-27_initial`),
  and cell 11/17 via `evoca_from_scan(scandir...)`/`nearest_params(...)` (default scan dir);
  `Untitled.ipynb` cell 2 L2 (`Runs/...top_evo_20_cfg131...`); `test.ipynb` cell 8 comment only.
- Drafts (branch `research-drafts`), each finds the repo root by walking up for a `Scans/` dir
  (cell 1, L11-13) and then reads `ROOT/"Scans/<name>/results.csv"` (cell 2 unless noted):
  `RD-robustness` (L1 `Scans/2026-05-18_rd_under_pure_evo`), `coevolution-substrate`
  (`Scans/2026-05-19_3c_...`), `neutral-model-methodology` (cells 3, 5, 7: causal_control v2, v1, 3c),
  `pure-evo-regime` (cell 2 L1-2: #1, #2), `viability-brackets` (R1),
  `resource-driving-inverted-U` (placeholder `../../Scans/<date>_...`, cell 1 L3, L11).
  **The `Scans/`-dir root heuristic stops working once `Scans/` is empty or gone** (use
  `evoca_py.py` or `Makefile` as the root marker).
  `neutral-model-methodology` reads scans owned by AM and LC: cross-workstream path.
- Docs: `Research/README.md:6-7,42`, `Research/research_board.md:105` (glob `../Runs/*R-D*`),
  `Scans/README.md` lines 15, 23, 29, 37, 73, 86, 97, 154, 178, 189, plus the
  `Docs/devlog/2026-05-17_research-campaigns.md` scan headers and
  `Scans/2026-04-27_scan_analysis.md:456-484` "Reproducing" commands and
  `Scans/2026-05-17_pure_evo_analysis.md:193-199`.

## 5. Other findings
- Only 14 Runs were tracked at start (review §1); `git ls-files Runs` currently lists 87 so the
  untracked set has since been added or the status is out of date; re-check before `git mv`.
  `ProbeLogs/` is entirely untracked and `Scans/2026-05-19_8_*` and the `__pycache__` dirs are
  untracked.
- `2026-04-27_top_*` Runs span two source scans (initial: cfg49/30/70... ; push_edges), so
  the `rd_under_pure_evo` bases are a mixture; fine under single-owner (RD).
- Proposed shorter names: `<date>_<scan-abbrev>_<class><NN>_cfg<k>_s<score>.evoca`
  with abbrev `ini`, `pe`, `ef` and class `top`, `evo`, `auto`
  (e.g. `Runs/2026-04-27_top_1_cfg236_score0.90.evoca` ->
  `Research/Reaction-diffusion/Runs/2026-04-27_pe_top01_cfg236_s0.90.evoca`); all 87 new names
  verified unique. Scan dirs drop redundant words (`pure_evo_joint` -> `joint`,
  `3c_coevolution_substrate` -> `3c_substrate`, `rd_under_pure_evo` -> `under_pure_evo`).
  The 04-29 `evo_focus` Runs keep their 04-29 date prefix. The `descriptor` inside each JSON
  is unchanged (does not need to match the file name).
- Structural risk (global-instructions trigger): this restructure **turns a configured path
  into a derived one and adds a resolver layer**; I recommend offering `/code-review ultra`
  on the harness path refactor once it exists.

## 6. Recommended order of moves
1. Commit the pending Oct-2026 docs and any untracked Runs/Scans (`git add` first so `git mv`
   preserves history). Ignore `__pycache__`.
2. Merge Tier-1 branches first (several edit `python/evoca_explore.py`).
3. Land the harness change **before** any file moves: a `find_repo_root()` helper, a
   `workstream=`/`root=` argument on `_resolve_scan_dir`, `evoca_from_scan*`, `export_recipe`, `import_run`
   (including the no-arg listing), and a legacy-name resolver backed by `artifacts.csv`. Keep the
   old defaults working (backward compatible) so the move can be done in small steps.
4. Fix all 13 scan drivers to use the root helper (identical mechanical edit) and rewrite the 3
   cross-workstream paths in 4.2. Verify by importing each scan module (not running).
5. `git mv` per `artifacts.csv`, one workstream at a time, in order **Reaction-diffusion**
   (Runs + 4 scans + analysis md; self-contained except `evo_focus`), **Activity-metrics**,
   **Lineage-co-evolution**, **Resource-driving**; each followed by `pytest tests/` and a
   headless `import_run` / `evoca_from_scan_top` smoke check on a moved Run, and a notebook
   check (`evoca_scans.ipynb` and the drafts).
6. Move the drafts (after `research-drafts` is merged) to `digests/` and replace their
   root-finder with the shared helper.
7. Update root notebooks, `Scans/README.md`, `Research/README.md`, and the board's `../Runs/*R-D*`
   reference. Leave `ProbeLogs/` in place (sink) and archive the manual logs last.
8. Last: the 13 `source_scan` fields in recipes need no rewrite (provenance only); keep old
   names in `artifacts.csv` as the resolver.
