
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

# Claude:

Good direction — a few reactions and questions before any code lands.

  Ternary egenes (0/1/*). Sharpening agreed. Storage cleanest as two uint8 per egene: value[6 bits] and mask[6 bits],
  where mask bit = 1 means non-*. Doubles per-cell egene storage from 8 → 16 bytes (still trivial). With mask=0x3F
  (all non-*), the new model behaves identically to the current one except for the −1 penalty term, so it's not
  strictly backward-compat — old scan winners will eat differently.

  +1/−1 scoring with normalization. I'd push back gently on normalising — the natural pressure of the system is that a
   richer egene (more non-*) has both higher max reward and higher tax, and forcing them to map to the same [0, 1]
  range erases the "complexity costs more, pays more when right" gradient that makes it interesting. My
  recommendation:

  - raw_score = sum over the 25 cell-positions covered by non-* orbits of ±1, range [−25, +25]
  - mouthful = (m_scale / 25) × max(0, max_over_active_egenes(raw_score)) × F_food

  Clamp at 0 because negative eating is biologically odd, but keep the un-normalised range so a 1-orbit egene tops out
   at +1 and a 6-orbit egene tops out at +25. Open to discussion if you have a specific normalisation in mind.

  etax. Two readings of "non-* elements":

  - (a) literal: count of non-* orbits (0–6) — tax cap at tax_per_egene × 6.
  - (b) count of non-* cell positions (0–25), so the corner orbit (4 cells) costs 4× the centre orbit (1 cell) —
  closer to actual cognitive cost.

  I'd default to (b); centre-orbit-only egenes with perfect match are the cheapest possible "cognition" then, which
  feels right.

  Measures review. With the ternary change, I'd add to the egenome probe (or a paired probe):

  - mean non-* count per egene (cognitive specificity)
  - mean non-* count per cell (sum over active egenes — total cognitive load)
  - mean eaten food per alive cell per tick — this is the measure that actually answers "are eaters getting better?"
  and we don't have a clean version of it yet (the closest is f_priv, but that's the stash, not the intake rate).

  Channon excess and lut_div/eg_ent stay useful but become coarser indicators next to those.

  Patch transfer. Clean concept — implementation details to nail down:

  - Patch contents: alive, LUT, egenes, active, v_curr. Should f_priv go too (the patch carries its energy reserves)?
  I'd say yes — otherwise the transplant starts hungry and dies before it can demonstrate fitness.
  - F_food stays from the destination (it's the environment, not the organism).
  - Same N required for the patch source and destination (different N would need stretching/cropping — separate
  problem).
  - API sketch: sim.save_patch(cx, cy, side) → dict; sim.load_patch(patch_dict, cx, cy). Both can be done on a paused
  sim.
  - Persistence: save/load to .evopatch JSON files under a new Patches/ dir.

  Order I'd take it. Items are independent, but #1 (ternary) changes the meaning of the existing "match" metric. So:

  1. Ternary egenes + scoring + tax(b) + new ternary-aware probe metrics.
  2. Eaten-food-per-tick metric (small, drops cleanly into egenome probe or a new one).
  3. Patch transfer (independent, do anytime).

  Worth noting: the current scan winners will not behave the same under the ternary change. We should plan a fresh
  scan after #1 lands rather than expecting old recipes to keep working.

  Want me to start on #1 with the (mask, value) encoding + clamp-at-0 + tax(b), or push back further on any of this?

## response

1. re. +-1 and normalization:

Not sure I completely understand 

   'the natural pressure of the system is that a
   richer egene (more non-*) has both higher max reward and higher tax, and forcing them to map to the same [0, 1] range erases the "complexity costs more, pays more when right" gradient that makes it interesting.'

I guess there is a point:  using egene complexity defined as non-* is clean.  making +-1 distinction muddies the situation.

Anyway I'm happy to try your suggestion as a starting point.  

I'm a little worried that normalizing to max 25 means that initial 'dumb' egenes will be getting very little food.  maybe too little to survive.  May have to adjust taxes accordingly.

2. etax:

definitely (b).

3. probes:  

Let's make an egene probe with all 3 traces.

4. order:

sounds good.  let's do it.
make sure we have commited all work to date, before beginning new commits.

