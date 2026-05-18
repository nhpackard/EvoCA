# Devlog

Durable record of **how and why the model is evolving** — the reasoning
layer that commit messages and the research board don't capture.

- `research_board.md` = status / decisions / timeline (what state things are in).
- commit messages = what changed in a diff.
- **`commit-log.md`** = running ledger of every commit/merge that
  lands on `main` (newest first) + push ranges — the at-a-glance
  "what has changed in the repo." Maintained by the orchestrator on
  every commit/merge to `main`.
- **devlog** (this) = *why*: motivation, key findings, design
  decisions, the alternatives that were rejected and why, flagged
  modelling choices, obstacles and how they were resolved.

One file per dated work batch (or per major workstream). Each
workstream section records: **Motivation · Key findings · Design
decisions & rejected alternatives · Flagged modelling choices ·
Obstacles · Status/commit**.

Convention going forward: every spawned workstream agent writes its
devlog section as a committed deliverable on its branch (not just a
chat report), and the orchestrator consolidates per batch here.
