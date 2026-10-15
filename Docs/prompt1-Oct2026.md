
This is a prompt for claude, regarding 
- firing up the workstreams; review of research_board.md and its embedded suggestions.
- structure of the research dir with associated subdirs for each workflow.

There will be some questions to resolve, but I would ultimately like agents fired up to do as much of the work described here in parallel as is possible.

Some of the  requests are asking for text response.  Rather than simply spew stuff to the CLI, i would like text response captured in md files. (but summary report to CLI is good.)  Could be added to the bottom of this file, or could be myriad new little files, either tagged with Oct2026, or maybe we make a Oct2026 subdir of Docs and put there the reasonably named files as responses to my queries.  Open to your suggestion re structuring your response.  Main thing is I don't want it to be ephemeral data living temporarily in this session.

# cleanup of dev branches
- review of existing branches, PRs, and potential PRs
- lets merge everything that seems mergable (not buggy, and with all major issues/decision points resolved)
- in a research_board section, let's list all existing dev branches associated with each workstream, with a sketch of the plan to conclude work in the branch and merge

# the Research dir

## research_board.md

- separate 'archival non-authoritative' sections into a different file, e.g. research_board_archival.md
- List of workflows is basically good, but a little too high context.  E.g. "H1 adaptive-vs-drift (needs D2+G1); §B food_inc/tax figure (needs L1+G1)"  I find the codes like H1, D2, G1 are actions in the table above, but that too is cryptic; are these actions needing more dev? or actions that will be implemented in a merge of the indicated branches?  Also the labels presumably have some historical sense, but they aren't particularly clear in the immediate context.  This can be ok; labeling actions to be referred to can reasonably have more or less arbitrary labels.  Not clear to me where the sections referred to live (e.g. §B).  let's have a little explanation.

Re. "Campaign results → Research/" section:  Let's go further than described.

- in the Research dir, let's have a subdir for each work stream.  E.g.
  - Wave
  - Criticality-bifurcation
  - Ring-complexity
  - Lineage-co-evolution
  - Reaction-diffusion
- In each of these dirs there should be
  - A growing notebook with most current dev in the workstream in the bottom section, with comments on directions and outstanding decisions.  Each section should include example runs to illustrate explored phenomenology, with the runs tailored to have the relevant probes for illuminating the work stream.
    - would be nice to coordinate this with the research board.  Open to suggestions.  Maybe have a little md file of issues and next directions that gets imported into the notebook and also into the research board.  Not sure how easy this is.  Md doesn't have a lot of import facility.  what do you think?
  - Runs dir and Scans dir for that particular workstream
    - I would like to eliminate the top level Runs and Scans dirs, ultimately.  They become a morass.  having runs and scans per work stream will simplify naming the files.
    - I am expecting the work stream notebooks to access particular results in the respective runs and scans dirs as examples

So... the action item here is not to use this structure going forward, but to initialize the structure with existing work/examples from each of the work streams.  *and* to use it going forward.

# Scans / Runs

Some of the scans/runs were using activity as a param to optimize.  After current work in the Wave workstream is concluded, we may have a new version of activity to optimize.  Pls make a list of the candidates

# A separate dev issue
 
Remind me what the difference is between Runs and Scans.

I know that in some of the Runs/Scans we were essentially exploring multiparam optimization.  E.g. for the reaction diffusion work stream, maximizing both structure and activity.

Remind me:  how are you handling dual param optimization?  If this is explained in an existing md file, point me to the explanation. If not, let's create a little explanatory md file.

Remind me: what is the complete list of where we are using dual param optimization?

I would like to consider using the restricted (cone) pareto approach that we developed in the pdt-server repo (which review, in ~/Projects/pdt-server/backend).  This would be an alternative to whatever you have been doing.  I would like to abstract the restricted pareto frontier approach out of pdt-server and package it up to be useable for any dual-param optimization.

    