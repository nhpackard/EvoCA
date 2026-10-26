# 2026-10-25: Docs threads absorbed into Research/, issue-ID links, aggregate ISSUES

N's request:
1. Remember that issue-ID references are links to their ISSUES.md, and apply that to every `.md`.
2. Add a script that aggregates the workstream ISSUES files into `Research/ISSUES.md`.
3. Give thoughts on `Docs/Energy analog.md`.
4. Move the live `Docs/` threads (the plan-xxx files and others) into the `Research/<Workstream>/`
   structure, with references from the workstream ISSUES files, and archive the loose ends in
   `Docs/Archive/`.

N asked that the moves come first, so that each link is written once against its final path.

## Plan

1. Classify every file in `Docs/` as one of:
   - a workstream thread (moved to `Research/<WS>/plans/`);
   - project-wide (stays in `Docs/`);
   - a finished plan or loose end (moved to `Docs/Archive/`).
2. Move with `git mv`, so history follows. Rewrite every reference in live files: repo-rooted
   paths, relative links that resolved to a moved file, and relative links inside moved files.
   Historical records (devlog, Archive, `research_board_archival.md`) and N's open working files
   (`Research/Chem/Chem.ipynb`, `prompt1-Oct2026.md`) are left as written.
3. Add a "Plans and design docs" section at the top of each affected ISSUES.md.
4. `tools/build_board.py` renders the workstream sections once and writes two copies:
   - the board section;
   - the new `Research/ISSUES.md`.

   Relative links are rebased for files in `Research/`. This was also a latent bug in the board
   copy, where any relative link in an ISSUES file had broken.
5. `tools/link_issue_refs.py` links bare IDs to their defining ISSUES.md. It is idempotent and
   skips:
   - definition sites;
   - same-file references;
   - code and existing links;
   - undefined IDs.

   `build_board.py` applies the same linking to the generated copies.
6. Tests through the public functions, with constructed fixtures; each rule shown to fail when
   broken.
7. Respond to the energy-analog question: `Docs/energy-gradient.md`.

## What was done

**Moves.**

| From `Docs/` | To |
|---|---|
| `plan-chem.md`, `2026-10-23_chem-model-v2.md` | `Research/Chem/plans/` |
| `plan-8-lineage-coevolution.md`, `plan-3c-coevolution-substrate.md`, `Oct2026/08-lineage-decisions-guide.md` | `Research/Lineage-co-evolution/plans/` |
| `activity-thoughts.md`, `Oct2026/07-activity-objective-guide.md`, `Oct2026/10-drift-shadow-findings.md`, `Oct2026/05-cone-pareto-package.md` | `Research/Activity-metrics/plans/` |
| `2026-10-21_spirals.md`, `2026-10-22_food-thresh.md` | `Research/Reaction-diffusion/plans/` |
| `EvoCA_research_stability.md` | `Research/Ring-complexity/plans/` |
| `2026-10-21_galleries.md` | `Docs/devlog/` |
| `2026-10-24_kaneko-literature.md`, `channon_ad_ecal2001.pdf` | `Docs/references/` |
| `Energy analog.md` | `Docs/energy-analog.md` (no space in the name) |
| finished plans, `model2.md`, old prompts, the October process reviews, two non-EvoCA files | `Docs/Archive/` (see its README) |

`Docs/tmp.md` (empty, untracked) was deleted.

**References.**
- 37 files were rewritten for the moves. A link check over every live `.md` finds 728 relative
  links and none broken.
- Not rewritten:
  - `Research/Chem/Chem.ipynb`, which has N's unsaved edits; its references to moved docs are
    stale until N saves and the paths are updated;
  - the historical records.

**Tools.**
- `tools/build_board.py` now also writes `Research/ISSUES.md`.
- New: `tools/link_issue_refs.py`. It found 116 defined IDs and linked 413 references in 20 files.

**Tests.** `tests/test_issue_tools.py`, 11 tests. Six rules are each shown to fail when broken:
link rebasing, protected spans, the same-file skip, the trailing ID boundary, the duplicate
definition check, and the linking of the generated copies. The full suite has 261 tests, all
passing.

**Docs.**
- `Docs/README.md` (what stays in Docs and why) and `Docs/Archive/README.md` (what was archived
  and why).
- The `Research/README.md` structure now includes `plans/` and the aggregate.
- `Docs/energy-gradient.md`.

## Open after this

- When `chemca` merges, its `Docs/plan-chemca.md` should go to `Research/Chem/plans/`.
- `plan-chem.md` milestone 3 (the reaction network, CHEM-N4) overlaps the `chemca` design. N to
  decide whether `chemca` supersedes it.
- `Docs/energy-flux-and-evolution.md` says L2 needs `chemca`. `energy-gradient.md` §4 argues that
  is wrong as designed. The anchor doc is left unedited while N is reviewing it.
