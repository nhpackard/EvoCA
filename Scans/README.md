# EvoCA Parameter Scans

Each subdirectory `<YYYY-MM-DD>_<name>/` is one scan campaign:

```
<date>_<name>/
  notes.md       — context, hypothesis, parameter ranges, post-hoc observations
  scan.py        — driver: builds param grid, fans out to multiprocessing.Pool
  results.csv    — output: one row per config, columns are summary metrics
```

The harness module is `python/evoca_explore.py` — `run_sim(params, ...)`
returns a dict of summary metrics. Drivers should import that and feed
configurations through it.

Run a scan with:
```bash
cd Scans/<date>_<name>
python3 scan.py
```
