"""Where EvoCA research artifacts live.

Runs (.evoca recipes) and Scans (campaign dirs) belong to a workstream:

    Research/<Workstream>/Runs/<file>.evoca
    Research/<Workstream>/Scans/<date>_<name>/{scan.py, scan_config.json, results.csv, notes.md}

The workstream list below is the single source; directories are derived
from it. A new recipe with no workstream given goes to the workstream in
the EVOCA_WORKSTREAM environment variable, else to Research/_unsorted/Runs
(sort it later with `git mv`).
"""
import os

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
RESEARCH = os.path.join(REPO_ROOT, 'Research')

WORKSTREAMS = (
    'Activity-metrics',
    'Resource-driving',
    'Lineage-co-evolution',
    'Reaction-diffusion',
    'Criticality-bifurcation',
    'Ring-complexity',
    'Chem',
)
UNSORTED = '_unsorted'


def _check(workstream):
    if workstream not in WORKSTREAMS and workstream != UNSORTED:
        raise ValueError(f"unknown workstream {workstream!r}; one of {WORKSTREAMS}")
    return workstream


def default_workstream():
    """EVOCA_WORKSTREAM if set (validated), else the _unsorted bucket."""
    ws = os.environ.get('EVOCA_WORKSTREAM')
    return _check(ws) if ws else UNSORTED


def runs_dir(workstream=None):
    """Research/<workstream>/Runs (workstream=None -> default_workstream())."""
    ws = _check(workstream) if workstream else default_workstream()
    return os.path.join(RESEARCH, ws, 'Runs')


def scans_dir(workstream):
    return os.path.join(RESEARCH, _check(workstream), 'Scans')


def runs_dir_for_scan(scan_dir):
    """The Runs dir beside a scan's workstream (Research/<ws>/Scans/<scan>
    -> Research/<ws>/Runs), so recipes emitted from a scan land with it.
    A scan outside that layout falls back to runs_dir()."""
    scans = os.path.dirname(os.path.abspath(os.path.normpath(scan_dir)))
    if os.path.basename(scans) == 'Scans':
        return os.path.join(os.path.dirname(scans), 'Runs')
    return runs_dir()


def all_scan_dirs(workstream=None):
    """Every scan dir (one containing results.csv) in one workstream, or all."""
    wss = (_check(workstream),) if workstream else WORKSTREAMS
    out = []
    for ws in wss:
        base = os.path.join(RESEARCH, ws, 'Scans')
        if os.path.isdir(base):
            out += [os.path.join(base, d) for d in sorted(os.listdir(base))
                    if os.path.isfile(os.path.join(base, d, 'results.csv'))]
    return out


def all_runs(workstream=None):
    """Every .evoca recipe in one workstream, or all (plus _unsorted)."""
    wss = (_check(workstream),) if workstream else WORKSTREAMS + (UNSORTED,)
    out = []
    for ws in wss:
        d = os.path.join(RESEARCH, ws, 'Runs')
        if os.path.isdir(d):
            out += [os.path.join(d, f) for f in sorted(os.listdir(d)) if f.endswith('.evoca')]
    return out
