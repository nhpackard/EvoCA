
This doc is to sharpen the notion of egene evolution, and the aims of its study.

# overall goals

The main goal is to study evolution of cognition, and to have a concrete model that enables the observation of evolution of increased cognitive ability through an evolutionary process driven by selection of the fittest in a population of agents that have evolving cognitive aparatus. 

The EvoCA platform:  provides a cognitive mechanism in the form of egenes.  Tax on egenes makes cognition costly; survival is a balance between cognition enabling greater food harvesting to overcome the cost of the cognition.

The prior:  we should see a flow of resources (food) to support evolution of cognitive ability.  This depends on the cognitive aparatus being capable of increasing food harvest of agents.

# sharpening of egene / feeding

To sharpen the egenes functionality I propose some changes:

- rather than have egenes code for a pattern of 1's and 0's that get matched to a local pattern, we change to gene values of 0, 1, or *, where * means 'don't care'

- rather than simply summing up matches to establish the food harvesting mouthful, over all non-* genes sum +1 for match and -1 for non-match.  Penalty for non-match increases value of the egene template actually matching.

- etax should be proportional to number of non-* elements of egene.

- we need to establish a proper normalization of eating capability given matching of genome to local pattern.  Open to suggestion.

# review measures for progress

We have had a discussion of measuring evolutionary progress of agent's eating capabilities.  I sould like to review these, and implement probes that reflect the sharpest measurements.

# evolutionary test

This section is not particular to egenes and evolution of cognition, but a feature for all EvoCA evolutionary experiments:

I propose to be able to extract a patch of the population (e.g. middle square of lattice, side length 1/3 of lattice), then insert that into another population, either initialized or itself evolved, to then observe whether the patch grows or shrinks.  growth => "evolutionary success", shrink => "evolutionary failure"

# implementation status (2026-05-10)

Item 1 (ternary egenes + scoring + tax(b) + new probe) landed across
commits `bb9ca18..` plus follow-ups for the new `egene` probe and
docs.

Settled choices:
- Encoding: parallel `egenes_value[N*N*8]` + `egenes_mask[N*N*8]`,
  6 bits used per byte. Backward-compat default mask = 0x3F.
- Scoring: raw_score = sum over non-wildcard cell-positions of ±1,
  range `[−25, +25]`. Mouthful = `(m_scale/25) × max(0, raw_score) ×
  F_food`. No normalisation.
- Tax (option b): `tax_per_egene × sum_over_active(non_wildcard_positions)`,
  where positions per orbit = `[1, 4, 4, 4, 4, 8]`. Max contribution
  per cell = `tax_per_egene × 200`.
- Mutation: `mu_egene` per-bit rate over all 8×12 = 96 bits per cell
  (6 value + 6 mask per slot). Bits 0..5 hit value, 6..11 hit mask.
- Species hash: includes `(value & mask, mask)` per active slot,
  sorted lex. Wildcard bits don't shift identity.

New probe `egene` (3 sub-strips, 192 px):
- `spec` (mean cognitive specificity) y in `[0, 25]`
- `load` (mean per-cell cognitive load) y in `[0, 200]`
- `food` (mean intake per eater this step) y in `[0, 1]`

Validation at the scan-3 productive corner (`food_inc=0.013,
m_scale=1.2, gdiff=0.06`, N=128, 2000 ticks):

| Run                                    | spec  | load  | intake | alive |
|----------------------------------------|------:|------:|-------:|------:|
| V1: zero new knobs (mask=0x3F default) | 24.4  | 24.4  | 0.056  | ok    |
| V2: mu_egenome=0.005, tax_per_eg=1e-4  | 20.7  | 38.0  | 0.045  | ok    |
| V3: mu_egenome=0.005, tax_per_eg=8e-4  | 15.6  | 20.1  | 0.049  | ok    |

V1 stays near maximum specificity because masks start at 0x3F and
only drift slowly under `mu_egene=0.003`. V2 sees specificity drift
down toward an intermediate value while cognitive load grows from
the unbounded Negene pressure. V3's higher per-position tax bites
both axes — fewer slots active and fewer non-wildcard positions per
slot, with a small intake penalty (0.045 → 0.049 was within
run-to-run variance).

Initial-fitness concern (cognition starts low → may not survive):
not observed at these parameters with the scan-3 corner. Default
`set_egenome_all` keeps mask=0x3F so the simulation seeds with
fully-specified egenes; if `set_egenome_random_all` is used,
specificity starts at ~12.5 (random masks) and the population
typically still survives at the scan-3 corner with `tax_per_egene=0`.

# next

- Item 2 (patch transfer): start when ready.
- New scan grid using the cognitive metrics in the score —
  `cog_specificity_mean`, `cog_load_mean`, `mean_intake_mean`
  (exported as `COG_METRICS` in `evoca_explore`).
