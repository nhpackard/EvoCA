# 2026-10-23: standard-model rescan of `2026-05-19_3c_coevolution_substrate`

**Why.** The source scan ran the legacy uncapped mouthful (`food_cap=0`, defect A9) at `m_scale` 3.0 in every row. The capped rule is the standard model (N, 2026-10-23; `Docs/2026-10-23_food-cap-standard.md`).

**What.** `scan.py` is a copy of the source driver, unchanged except for a shim after its imports that forces `food_cap=1` in `EvoCA.init`, and `N_WORKERS` from the environment. Arms, seeds, phases and protocol are as in the source. Today's code was used (the neutral shadow has its own RNG stream since 2026-10-15). Run on torque with `tools/rescan_standard_2026-10-23.sh`: 80 rows, 303 s, `determinism_ok=True` on every row that records it.

**Used by:** the Lineage-co-evolution notebook and the `coevolution-substrate` digest (main thread). The source scan is kept as the recorded legacy-rule result.
