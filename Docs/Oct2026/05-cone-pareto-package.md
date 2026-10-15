# 05 - conepareto package (port of pdt-server cone/box Pareto selection)

For N. Repo: `~/Projects/conepareto/` (git init, one commit 28c2b36, no remote). pdt-server untouched.

## Files
- `src/conepareto/`: `layers.py` (pareto_layers, pareto_front), `normalize.py` (normalize_values, resolve_ranges, normalize_with_ranges), `regions.py` (7 regions + registry, quantile floors/ceilings/band bounds, validity_mask), `relax.py` (relax_acceptance), `select.py` (select + Report), `hypervolume.py` (hypervolume_2d), `__init__.py`.
- `tests/test_parity.py` (R parity), `tests/test_api.py` (contract/validation), `tests/fixtures/make_fixtures.R` + committed CSVs (2.7 MB).
- `examples/two_metric_selection.py`, `README.md` (quick start, ported design notes nearly verbatim, R-is-reference statement), `pyproject.toml` (numpy only; pytest optional).

## R functions ported
ParetoLayer2D, ParetoFront2D (helper used by hypervolume), ResolveParetoObjectiveRanges, ParetoHypervolume2D, SymmetricWedge/Box/BoxQuantile/InverseBox/InverseBoxQuantile/Band/BandQuantile Acceptance, BoxQuantileFloors, InverseBoxQuantileCeilings, BandQuantileBounds, ValidityMask, RelaxAcceptance, SelectParetoSweepCandidates (as `select`), NormalizeValues (from InitFunctions.R; needed because the sampler normalizes with it).

## Deviations (all naming/API, no intended semantic change)
- Generic names: response -> `first`, secondObjective -> `second`, snake_case (mapping in regions.py docstring and README).
- Registry names: wedge, box, qbox, ibox, iqbox, band, qband (R function names also accepted). No asymmetric wedge, by design (README explains).
- R `quantile(type=1)` is reimplemented rather than using numpy's `inverted_cdf` (R's 4*eps fuzz). Note: that fuzz turns out to be behaviorally inert for type 1 (the `h` flag compensates), so removing it is an equivalent mutant; the `h = nppm > j` logic is what matters and is covered.
- `select` returns a `Selection` with a `Report` (layers_used, truncated_layer, floors_achieved = min score per axis of the selection, worst_imbalance = max |s1-s2|, relaxed, scale, admitted, cone_width, validity_admitted/relaxed). floors_achieved and worst_imbalance are computed by me from the selection (the R sampler computes them in its diagnostics, not in SelectParetoSweepCandidates); fixtures check them against the same formulas applied in R.
- R `stopifnot` errors become `ValueError`. Error paths are tested in Python only (not parity-tested).

## Parity results
Fixtures from real R 4.5.2 runs sourcing pdt-server read-only (SelectParetoSweepCandidates and NormalizeValues are extracted by parsing, since their files do other things at top level). Six pools: ties, dups, boundary (values exactly on 0.5 / diagonal / axes), 6x6 integer grid, random-300 rounded to 2 decimals (incidental ties), and the real Gen1 pool (10 000 candidates from `DualObjectiveParetoSweepDiagnosticsGen1.RData`; the R script also asserts R's ParetoLayer2D reproduces the saved layers).
- 23 region/parameter cases x 6 pools = 138 mask comparisons; 736 relaxation cases (incl. validity guard on, and guard-dropped); 1610 selection cases (3 within-layer orders, validity configs) entered through public `select`; layers, fronts, hypervolume (pooled and fixed ranges), normalization, quantile bounds per pool.
- Result: 39 tests pass (33 parity + 6 API). 581 of the selection cases are relaxed, 198+70 of the 736 relaxation cases relax the region, and 207 are cases where the validity guard had to be dropped, so the fixtures are not trivial.
- Discrimination check: (1) `>` to `>=` in the group-vs-larger-s1 comparison -> test_layers[ties] fails; (2) breaking tied-s1 grouping -> fails; (3) removing the same-s2 member filter (duplicate/tie handling) -> fails; (4) `h = nppm > j` replaced by False in the type-1 quantile -> mask tests fail. All restored; suite green again.

## Not ported
SecondObjectiveFunctions.R (domain-specific, per instructions). Also not ported: the sampler's pareto_max_imbalance / inclusion-target search (derive a wedge angle from a target imbalance or frontier inclusion) and the lambda-sweep sampler. They sit around `select` rather than inside it; say if you want them.

## Open questions for N
1. Should the imbalance-to-angle search (R ~line 458-600 of the sampler) be ported as a helper? It is the only part of the sampler that derives a region parameter from a stated goal.
2. The example shows a box asking for floors 0.4/0.4 on a 2000-pool with k=40 is relaxed (achieved 0.312). That is correct R-style behavior, but for EvoCA use you may prefer an error option instead of silent-ish relaxation (it is reported). Want a `strict=True` flag? Not in R, so not added.
3. EvoCA structural note: this is a new standalone module with its own R-parity fixtures; nothing in EvoCA imports it yet. Offering `/code-review ultra` on the package is optional; nothing pluggable beyond the region registry was introduced beyond what R has.
