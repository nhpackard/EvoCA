# 2026-10-23: standard-model rescan of `2026-05-16_causal_control`

**Why.** The source runs used the legacy uncapped mouthful (`food_cap=0`, defect A9) at `m_scale` 1.2. The capped rule is the standard model (N, 2026-10-23; `Docs/2026-10-23_food-cap-standard.md`).

**What.** `driver.py` is a copy of the source driver, unchanged except for a shim after its imports that forces `food_cap=1`. Arms, seeds, N and steps are as in the source. Today's code was used (the neutral shadow has its own RNG stream since 2026-10-15). Run on torque with `tools/rescan_standard_2026-10-23.sh`; see `run.log`.

**Used by:** the Activity-metrics notebook and the `neutral-model-methodology` digest (main thread). The source run is kept as the recorded legacy-rule result.
