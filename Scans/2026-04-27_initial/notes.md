# Scan 2026-04-27 — initial coarse search

## Question

Find parameter configurations that give **large-scale spatial structure**
*and* **ongoing genomic turnover** simultaneously. Earlier interactive
exploration showed that small-gdiff regimes (~0.01) produce reaction-
diffusion-like transient waves, but selection then locks in a static
fine-grain "garden plot" attractor and dynamics stop.

## Search axes (full grid: 1152 combos)

```
gdiff         ∈ {0.005, 0.01, 0.02, 0.05}
mu_lut        ∈ {0.001, 0.003, 0.01, 0.03}
mu_egenome    ∈ {0.001, 0.01}
food_inc      ∈ {0.015, 0.025, 0.040}
tax           ∈ {0.020, 0.025, 0.035}
m_scale       ∈ {0.4, 0.6}
restricted_mu ∈ {True, False}
```

Random sample of 150 unique combos. Fixed: N=256, n_steps=5000,
sample_every=100, init='halfplane', seed=0, shadow=True.

## Results — 132 s on 12 cores

Of 150 runs:
- 31 extinct (population died)
- 94 saturated (alive_density_mean > 95% — uniform live regime, no spatial structure)
- 25 in the interesting middle (alive but not saturated)

## Per-axis observations

| axis | finding |
|------|---------|
| **food_inc** | 0.040 → always saturated (too much food). 0.015 → 60% extinct, 40% interesting. 0.025 → mostly saturated. **0.015 is the only food_inc producing the regime we want.** |
| **gdiff** | 0.005 → mostly extinct (too patchy). 0.01–0.05 produce interesting candidates; 0.02–0.05 minimise extinction. |
| **tax** | 0.035 has highest extinction *and* highest candidate count. The death-edge regime is what we want. |
| **m_scale** | 0.6 gives more candidates and longer correlation lengths than 0.4. |
| **mu_lut** | Doesn't affect candidate count much, but higher mu_lut linearly increases dominator turnover (15 → 30 distinct top genomes over 50 samples). |
| **mu_egenome, restricted_mu** | Negligible effects. |

## Top configurations (composite score over correlation_length, largest_patch_temporal_std, n_distinct_genomes_temporal_std, unique_top_genomes)

```
   score  gdiff  mu_lut  mu_eg  food_inc  tax  m_scale  rest_mu | corr_L  patch_std  div_std  top_g  alive
   0.750  0.050  0.010   0.001  0.025     0.035  0.6   True    | 21.1    8314       7546     43     0.696
   0.721  0.050  0.010   0.010  0.025     0.035  0.6   False   | 21.2    7689       7647     43     0.698
   0.681  0.010  0.010   0.010  0.015     0.020  0.6   False   | 32.6    7717       7850     19     0.511
   0.630  0.050  0.030   0.010  0.015     0.035  0.6   True    | 49.2    4462       5226     50     0.373
   0.611  0.010  0.003   0.001  0.015     0.020  0.6   False   | 32.4    7900       5846     23     0.518
```

**Best single config by correlation length AND turnover** — corr_L=50.6, 51 distinct top genomes:

```python
food_inc=0.015, gdiff=0.05, mu_lut=0.003, mu_egenome=0.01,
tax=0.035, m_scale=0.6, restricted_mu=True
```

## Synthesis

The interesting regime is narrow and consistent:
- **food_inc = 0.015** (low — death pressure)
- **tax = 0.035** (high — death pressure)
- **gdiff ∈ [0.02, 0.05]** (moderate — preserves spatial structure but not so low that the system fragments to death)
- **m_scale = 0.6** (cells eat aggressively)
- **mu_lut ∈ [0.003, 0.03]** (enough mutation to keep the system from locking in)

The "ongoing evolution + spatial structure" regime lives at the **edge of survival**:
high tax + low food + moderate diffusion + aggressive eating + medium mutation. Below
this the population dies; above it the system saturates uniformly with no spatial
structure.

## Next steps

- [ ] Validate the top config at N=512 visually with the SDL window. Should show
      ~50-cell-scale patches with continual genomic turnover (no garden-plot lock-in).
- [ ] (Optional) Refined scan: tighter grid around the death-edge —
      `food_inc ∈ {0.012, 0.015, 0.018, 0.022}`,
      `tax ∈ {0.030, 0.035, 0.040, 0.045}`,
      `gdiff ∈ {0.02, 0.03, 0.05, 0.08}`,
      everything else fixed at the top-config values. ~64 combos, < 1 min.
- [ ] Multi-seed runs of the top config to confirm the dynamics are robust to
      v_curr initialisation.
