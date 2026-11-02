# Galleries: one interface to lists of recipes and scan rows (2026-10-21)

Branch `rd-galleries`. Prompted by N's request: in `Reaction-diffusion.ipynb` the "example runs of the picks" needed variables defined in cells far above; recipe lists had a different shape in every section (`GALLERY` tuples, `WAVES` tuples, the `picks` DataFrame indexed by pool position), and a figure title such as "CB.cfg102" did not lead to a recipe without a search. N asked for a uniform way to access the libraries of scans, per-section setup cells, and noted the galleries could use parallelism.

## Plan

1. **One type for every list of configurations.** `Gallery` (ordered, named list of `Entry`) in `python/evoca_gallery.py`. An `Entry` has exactly one source, a recipe file or a scan row (`scan`, `cfg`), plus `overrides` (through `update_<key>()` setters), `title`, and free-form `meta` (scores, Pareto layer, methods).
2. **Indexing that matches the figures.** `G[i]`, `G['push_edges.cfg246']`, `G['cfg246']` (unique '.'-suffix), `G[[...]]`, slices, `where(**meta)`, `split(key)`. Evaluating a gallery shows its index table, including the metaparameters a load *actually* sets (corrections and overrides applied, read back from a sim). Figure titles carry `[i] label`.
3. **Actions on an entry.** `load` (headless, reproducibly seeded), `show` (SDL window, starts running), `scan_row()`, `effective_params()`.
4. **Libraries on disk.** Each notebook section saves the galleries it builds to `Research/<WS>/Galleries/<name>.json` (written only when the content changes; `owner` names the section). `Gallery.open(name)` works from any kernel.
5. **Parallel, cached runs.** `Gallery.run(N, n_steps, every, snaps, cms, seed, seed_v, metrics, overrides, workers)` runs entries in spawn worker processes (the C library holds one lattice per process) and caches each result in `.cache/gallery_runs/` (gitignored). The cache key is the entry's source content (recipe bytes, or scan row plus `scan_config.json`), the merged overrides and run arguments, and a hash of the C library and the Python loaders. Labels and titles are not part of the key. `GalleryRun` gives `summary(t_start)`, `tiles()` and `strip()`.
6. **Notebook restructure.** A global Setup section (imports, shared helpers `load_0427`, `osc`, a listing of saved galleries), then one `# ── section setup` cell per section that defines every name the section uses. Each figure cell calls `G.run(**SPEC)` itself, a cache hit after the first time, so any cell runs after the global setup plus its section setup.

## Decisions taken (flag for N)

- **Pareto picks load from their scan rows**, not from the `top_N` recipes, so they get the scan-row load-time corrections (`food_cap=0`, pre-rename `mu_egenome` read as `mu_egene`) without `sim.update_food_cap(0)` by hand. The recipe file for each pick is kept in `meta['recipe']` for reference.
- **RXD-D4 is still N's.** The regime gallery, decay and exploration galleries load their recipes *as written* (legacy recipes therefore run with `food_cap=1`). The table's `food_cap` column makes this visible. To decide D4 for the notebook only, add `overrides={'food_cap': 0}` to those entries.
- **Seeded initial v everywhere.** The old `run_track` used `seed_v=False` (an unseeded v draw), which a cache would freeze arbitrarily; `Gallery.run` seeds v from `seed`. This changed the regime gallery, exploration and decay figures slightly (below).

## Summary of what was done

- `python/evoca_gallery.py` (new): `Entry`, `Gallery`, `GalleryRun`, `Track`, `all_galleries`, `scan_label`, `clear_cache`.
- `tests/test_evoca_gallery.py` (new, 15 tests): lookup and ambiguity, sub-galleries, scan-row corrections through `load`, overrides through `load`, `show` forwarding to `run_with_controls` (monkeypatched), table `food_cap`, parallel equals serial, cache hit does not re-run, cache key follows overrides/spec and not label/title, run-time overrides take effect, save/open round trip without churn. Full suite: 162 passed.
- **Defect found by the tests and fixed:** the first cache key hashed entry overrides and run-time overrides separately, so the same effective run (`food_cap=0` set on the entry vs passed to `run`) was cached twice. Now keyed on the merged overrides, as `_run_job` applies them; `test_run_overrides_take_effect` failed before the fix (3 cache files for 2 distinct runs).
- `Reaction-diffusion.ipynb` rebuilt (59 cells) from N's working copy: seven galleries (`rd_explore`, `rd_waves`, `rd_boom_variants`, `rd_regimes`, `rd_decay`, `rd_pareto_picks`, `rd_pareto_layers` with `PARETO_LAYERS = PARETO.split('layer')`); interactive cells are one commented line each (`# PICKS.show(3, N=512)`).
- **Checks against external evidence.** The food_cap 0/1 wave table, the 10k-step persistence table, the 29-pick table and the per-method medians reproduce the previous execution exactly (same seeds), which also shows the scan-row door and the recipe-plus-override door give identical runs for the picks. A scripted check ran every code cell in a fresh kernel with only the global setup and its section setup: 30/30 ran.
- **Wall time.** Full notebook 211 s before; 70 s now from a cold cache (15 workers); seconds when cached.
- **Text corrected after re-execution.** (a) The decay observation said cfg39 was "not clearly worse" than cfg35 at N=128; the window means (3.8 vs 6.7 now, 3.8 vs 7.0 in the earlier unseeded execution) say it is worse after t~5000, matching the recorded N=256 ordering. (b) `fill_dead_region` is described as what the figure shows: evolving colored domains fill the dead half and invade the frozen grey GoL half. (c) cfg236's gallery description no longer claims the lowest alive density at t=600.
- `Research/README.md`: `Galleries/` in the tree and a paragraph on `evoca_gallery`. `.gitignore`: `.cache/`.
- `rd-galleries-test.ipynb` (repo root, un-run): PASS/FAIL review cells.

## Structural and scaling review

- **Artifacts.** Each gallery file has one owner (the section that builds it, recorded in `owner`). Run results have one owner (the cache). Nothing is produced in two places.
- **Opportunities (not implemented):**
  - *Other workstream notebooks* (Activity-metrics, Criticality-bifurcation, ...) keep their own tuple lists and hand-rolled run loops; adopting `Gallery` there is mechanical (a few cells each).
  - *The section-independence check* (fresh kernel per cell, global + section setup only) is a scratch script; as `tools/check_nb_sections.py` it could gate every lab notebook. Cost: about 1 s per cell with a warm cache.
  - *`Gallery.run` metrics are fixed* to `_spatial_metrics` or alive fraction; a `metrics=callable` hook would let lineage/activity notebooks use the same runner. Small change, but it makes the runner pluggable, so it would deserve its own reachability test.
  - *Cache growth* is unbounded (`clear_cache()` exists; picks at N=128 are ~100 kB per entry). Fine at this scale.
