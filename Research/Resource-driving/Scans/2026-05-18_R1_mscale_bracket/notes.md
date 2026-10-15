# R1_mscale_bracket (R1)

Bracket the two levers #1 left pinned at edges: `m_scale ∈
{2.5,3.5,5.0}` (was ceiling-pinned at 2.5) and `food_inc ∈
{0.005…0.018}` (extend below #2's 0.010 low edge). Mirrors #1 joint
design exactly (300 cfg, 8000 ticks, N=256, same metrics) so rows
compare directly. 300 rows, 0 errors, 743 s on torque. Launched
detached (`nohup`+`disown`) — survived a multi-minute VPN/torque
connectivity blackout untouched; recovered & analysed.

**Finding (1-line):** both levers are now **bracketed, not
runaway** — the `m_scale` *productive* optimum is **interior ≈2.5–3.5**
(5.0 is maximally *viable* but under-represented among metric winners),
the `food_inc` optimum is at the **HIGH** end 0.013–0.018 (the
sub-0.010 extension is closed: lethal at m_scale 2.5/3.5, no metric
winners anywhere), and #1's "low `mu_lut` wins" is revealed to be
**conditional on m_scale ≤ 2.5** — with m_scale bracketed up the
mu_lut optimum shifts high (top-20 mode 0.15). New mu_lut×m_scale
interaction; non-monotone U-shaped viability (worst extinction at
*intermediate* m_scale under food scarcity, lowest at m_scale 5.0).

Full analysis → `../../Docs/devlog/2026-05-17_research-campaigns.md`
Reasoning trail / caveats → same devlog section
