

**N**

❯ I basically agree with your contention:  we have successfully found a regime of ongoing interesting evolution.

But...

I still think we should be able to see waves somehow.  You called them predator-prey waves.  I have called them
(in other model contexts, eg genelife) reaction diffusion waves.  Anyway, such waves are intuitively to be
expected; I'm surprised we can't see them easily.

**Claude**

⏺ I think you're right that the intuition is correct — consumer-resource + diffusion + finite recovery time should
produce R-D waves. The fact that we can't elicit them easily is probably telling us about specific structural
choices in EvoCA that damp the oscillation. Three candidates:

1. No refractory delay. Classical R-D oscillations (BZ, FitzHugh-Nagumo, predator-prey waves in genelife if I
remember the lattice variant correctly) require a delay between consumption and re-availability — either a slow
reaction intermediate, or a recovery period during which the resource genuinely cannot be re-consumed. In EvoCA,
food regen and eating happen in the same step on the same field. No delay means perturbations decay rather than
propagate as oscillations.

2. Reproduction-as-migration, not local growth. When a cell hits f_priv ≥ 1, the child goes to the minimum-f_priv
neighbor — by construction, offspring move from rich → poor regions. That's a Fisher-KPP-style flux operator:
dense regions monotonically seed sparse ones, but there's no mechanism for a dense region to crash in place and
refill. Cells just emigrate as conditions deteriorate. Compare to most classical R-D systems where the "predator"
stays put and overshoots locally before crashing.

3. Selection actively prefers stable solutions. This is the deepest factor. The garden-plot pattern you saw isn't
accidental — once a non-trivial mutation rate is on, evolution finds attractors that are stable (= survive).
Oscillating regimes are by definition less stable than fixed-point ones, so over evolutionary time the system
slides off them. R-D waves in pre-evolutionary models (where genome doesn't adapt) don't have this issue. In
genelife, I'd bet either the genome doesn't directly modulate reproduction efficiency, or the timescales of
evolution and oscillation are well separated.

What I'd try, in roughly increasing-effort order:

- Sparse random init (set_alive_fraction(0.02) instead of halfplane). Spatially heterogeneous starting density may
seed local fronts that interact rather than smoothing into a uniform sea.
- Disable evolution for a probe run: mu_lut=0, mu_egenome=0. Lock the LUT to GoL and see whether non-evolutionary
R-D dynamics work in this model. If you see waves there but not with evolution on, point #3 is dominant. If you
see no waves there either, points #1–2 are.
- Add a refractory mechanism — small C-level change: a per-cell "regen cooldown" counter such that a cell that
just had its food eaten or a cell that just died can't gain food for K steps. This directly inserts the missing
delay. Easily reversible if it doesn't help.

If you've seen R-D waves cleanly in genelife with comparable primitives, it'd be useful to look back at what those
primitives actually were — odds are decent that one of #1 or #2 is the structural difference. Want me to sketch
the refractory cooldown patch as a switch you can toggle on/off, so we can find out empirically?