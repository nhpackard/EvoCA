# Test suite

Run from the repo root:

```bash
python3 -m pytest -q tests/
```

`conftest.py` rebuilds `C/libevoca.dylib` once per session if missing or
stale, and provides a `make_sim` factory fixture (fresh `EvoCA`,
auto-freed). 40 tests / 10 files. Each line below is what the test
actually asserts; the workstream tag links it to `Docs/devlog/`.

## CA core — `test_gol.py` (GoL oracles)

- `test_blinker_period_2` — 3-cell blinker oscillates with period 2.
- `test_block_still_life` — 2×2 block invariant for many steps.
- `test_glider_translates` — glider returns to shape shifted (1,1) after 4 steps.
- `test_gol_population_stable_without_mutation` — with mu=0, exactly 1 distinct LUT genome over a long run (no silent LUT corruption).

## Tax / survival — `test_tax_survival.py`

- `test_high_tax_kills_all_by_threshold` — all cells die exactly by `ceil(f0/tax)`, and not before.
- `test_zero_tax_population_persists` — tax=0 with no death pressure keeps the whole population alive.
- `test_dead_cell_state_is_zeroed` — dead cells have v=0 and a zeroed LUT (structural invariant).

## Determinism — `test_determinism.py` (S1)

- `test_same_seed_is_reproducible` — same `seed` ⇒ identical key metrics in one process (the `evoca_set_seed` fix).
- `test_different_seeds_differ` — different seeds produce different results (seed actually does something).

## Metric regression — `test_metrics_regression.py` (§2 fix)

- `test_run_sim_emits_component_normalised_excess` — `run_sim` still emits both raw and component-normalised excess, finite (locks in the §2 fix).

## Lineage — `test_lineage.py` (S2a, §1)

- `test_default_off_is_zero_cost_and_zero_data` — lineage OFF: accessors all-zero, no added work.
- `test_parent_hash_matches_single_live_genome` — ON + mutation OFF: every child's parent-hash == the single live genome hash.
- `test_known_single_reproduction_link` — one founder in a dead field: exactly one child, correct parent-hash / birth_id / parent_id.
- `test_lineage_on_does_not_perturb_ca` — additivity: GoL blinker still period-2 with lineage ON.

## Ring-dependence tax — `test_ring_tax.py` (S2b, §8)

- `test_level1_taxed_less_than_level3` — same `tax_ring>0`: a level-1 LUT loses strictly less per step than level-3, and the level-3 world dies first.
- `test_tax_ring_zero_is_bitwise_noop` — `tax_ring=0.0` reproduces prior dynamics byte-for-byte (mutation-tested guard).

## Recipe round-trip — `test_recipe_roundtrip.py` (D2 / D4)

- `test_tax_lut_ring_roundtrip[init|final]` — `tax_lut`/`tax_ring` survive export→import for both recipe kinds (D2).
- `test_params_includes_tax_lut_ring` — `params()` reports both knobs.
- `test_full_param_set_roundtrip[init|final]` — `mu_egenome`/`p_dup_egene`/`tax_per_egene` also round-trip (D4).

## Fixed-space neutral shadows — `test_fixed_space_shadows.py` (S2d, §2)

- `test_run_sim_emits_fixed_space_excess` — `run_sim` emits finite `eg_/dyn_excess_pc` final+slope.
- `test_excess_pc_independent_of_channon_shadow_flag` — fixed-space shadows run inside `sim.step()`, independent of the Channon-shadow flag.
- `test_neutral_drift_excess_is_zero` — zero mutation (no realised flux) ⇒ excess ≈ 0 (neutral-null invariant).
- `test_selection_excess_clearly_positive_and_above_neutral` — strong selection ⇒ excess clearly > 0 and ≫ neutral.

## Patch transfer — `test_patch_transfer.py` (S2c, §11)

- `test_extract_stamp_roundtrip_genomes_exact` — extract→stamp reproduces alive-cell genomes byte-exactly.
- `test_roundtrip_preserves_species_hash` — self-stamp preserves the distinct-species count.
- `test_self_transfer_approximately_neutral` — the §11 null: a population into its own field stays ~constant.
- `test_invade_only_writes_dead_cells` — `mode='invade'` never displaces a resident.
- `test_run_assay_returns_well_formed_dict` — the assay harness returns a well-formed result dict.

## Bifurcation harness — `test_bifurcation.py` (S2e, R3)

- `test_reducer_fixed_point` / `_exact_constant` — a constant/decaying series ⇒ one recurrent value.
- `test_reducer_period_two` / `_with_transient_and_noise` — alternating signal ⇒ exactly two recurrent values (robust to transient+noise).
- `test_reducer_period_four` — period-4 signal ⇒ four recurrent values.
- `test_feigenbaum_on_synthetic_cascade` — geometric-gap thresholds ⇒ δ ≈ 4.669.
- `test_feigenbaum_too_few_thresholds` — <3 thresholds ⇒ no estimate.
- `test_doubling_thresholds_and_delta_roundtrip` — fabricated doubling sweep → thresholds → δ end-to-end.
- `test_to_scatter_flattens` — sweep results flatten to (xs, ys) scatter.
- `test_sweep_driver_smoke` — one tiny real sweep: driver/sim wiring intact.
