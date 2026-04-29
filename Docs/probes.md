# EvoCA Probes Reference

Probes are optional strip-chart windows enabled via the `probes` dict
parameter to `run_with_controls()`.  Each probe renders in its own SDL
window, stacked to the left of the main lattice window.

Call `available_probes()` (from `python/controls.py`) to list all probe
names and their one-line descriptions at runtime:

```python
from python.controls import available_probes
available_probes()
```

## Probe Summary

| Probe name       | Window size   | What it shows                                          |
|------------------|---------------|--------------------------------------------------------|
| `env_food`       | 512 x 128     | Mean +/- std of environmental food F(x) over time      |
| `priv_food`      | 512 x 128     | Mean +/- std of private food f(x) over time            |
| `births`         | 512 x 128     | Mean +/- std of births array over time                 |
| `activity`       | 512 x 256     | LUT genome activity (scrolling hash-colored strip)     |
| `eg_activity`    | 512 x 256     | Egenome activity (scrolling hash-colored strip)        |
| `eg_food`        | 512 x 256     | Cumulative food per egene byte (mouthful split across max-match-tied winners) |
| `egenome`        | 512 x 256     | 4-strip Negene/eating diagnostics: mean ± std band, distinct egene values, mean max-match, frac at Negene_max |
| `lut_complexity` | 512 x 128     | Stacked area: green=n1 only, yellow=n1+n2, red=all     |
| `eg_pop`         | 512 x 128     | Stacked area: population fraction per egenome value    |
| `q_activity`     | 512 x 128     | Activity quantile deciles (log-scaled strip chart)     |
| `entropy`        | 512 x 128     | Local-pattern Shannon entropy over time                |
| `pat_activity`   | 512 x 256     | Local-pattern activity (scrolling hash-colored strip)  |
| `ts`             | 512 x 256     | Grouped scalar time-series (two 3-trace strips)        |

Example:

```python
run_with_controls(sim, probes={
    'activity': True,
    'eg_activity': True,
    'lut_complexity': True,
    'eg_pop': True,
    'env_food': True,
})
```

---

## Activity Tracking

### LUT Activity

Tracks cumulative presence of each distinct LUT genome over time.
Each genome is identified by its FNV-1a hash (cached in `lut_hash_cache[N*N]`).
An open-addressing hash table maps hash -> `{activity, pop_count, color}`.

- `evoca_activity_update()`: clears pop_counts, scans alive cells, increments activity
- `evoca_activity_render_col(col, height)`: renders one column of the scrolling
  strip chart. Alive genomes in full color; extinct genomes dimmed (RGB x 0.15).

**Y-axis saturation formula** (from genelife):

    y = (H-1) - (H-1) * act / (act + ymax)

This is a hyperbolic saturation curve: `act=0` maps to the bottom, `act=ymax`
maps to mid-chart, and `act->inf` approaches the top.  Tunable via `act_ymax`
halve/double buttons (default 2000).

### Egenome Activity

Mirrors LUT activity for the 6-bit egenome (fiducial eating pattern).  Since
there are only 2^6 = 64 possible egenomes, uses fixed-size arrays instead of
a hash table.  Wild-type egenome is colored white; mutants get FNV-1a hash
colors.  Separate `eg_act_ymax` halve/double buttons.

---

## LUT Complexity Probe

Classifies each alive cell's LUT by the minimum ring set it depends on:

- **Level 1** (green): rule depends only on (v_x, n1) -- constant across n2 and n3
- **Level 2** (yellow): rule depends on (v_x, n1, n2) -- constant across n3
- **Level 3** (red): rule depends on all three rings (v_x, n1, n2, n3)

Dead cells are excluded from counts.  The stacked area chart shows population
fractions at each complexity level over time.  Useful for studying whether
evolution discovers higher-ring dependencies starting from simple n1-only rules.

---

## Egenome Population Probe (eg_pop)

Stacked proportional bar chart showing the population fraction of each of the
64 possible egenome values over time.  Each time step renders one column.
Band height is proportional to population share; every populated egenome gets
at least 1 pixel; the last band absorbs any rounding remainder.

**Fixed band order**: egenomes are stacked in a fixed order (odd indices
descending, then 0, then even indices ascending) so that the wild-type
egenome (index 0) sits at the vertical centre.  Band positions are stable
across frames -- bands do not re-sort by population each step.

Colors match those used by the `eg_activity` probe (white = wild-type,
FNV-1a hash colors for mutants).  Internally reads `eg_pop[]` which is
populated by `evoca_eg_activity_update()` -- the probe ensures this update
runs even when `eg_activity` is not enabled.

---

## Egene Food Intake Probe (eg_food)

Same scrolling-strip layout as `eg_activity`, but the per-egene quantity
is **cumulative food obtained**, not presence.  At every eating step, a
cell's mouthful (after the headroom clamp) is split equally across the
egenes that tied for the best match — option-c attribution.  If a cell's
slot 3 and slot 5 both achieved match=18, each of their egene-byte
buckets gets `mouthful / 2`.

Stored internally as `uint64_t eg_food[64]` scaled by `1e6`, so 6 decimal
digits of resolution per food unit (mouthful is bounded by 1.0 by
construction).  Rendered with the same hyperbolic Y-axis saturation as
`eg_activity`:

    y = (H-1) - (H-1) * f / (f + ymax)

Tunable via `evoca_set_eg_food_ymax()` (default `1e6`, ≈ 1 unit of
food).  Use `sim.get_eg_food()` from Python to retrieve the raw table.

This probe answers "how much food has each egene actually fed into the
population?" — complementary to `eg_activity` ("how often is each egene
present?").  An egene that's common but in poor matching neighbourhoods
will have high activity but low food; an egene that fits the local food
field well shows the opposite.

---

## Egenome Stats Probe (egenome)

Four stacked sub-strips (each 64 px tall, 256 px total) charting Negene
distribution and eating efficiency over time.  All sub-strips have
**fixed y-ranges** so values can be read directly without auto-scale
context.

| Sub-strip | Trace                              | Y range  | Colour |
|-----------|------------------------------------|----------|--------|
| 0 (top)   | Mean Negene with ± std band        | [0, 8]   | green line, dim-green band |
| 1         | Distinct egene values (out of 64)  | [0, 64]  | orange |
| 2         | Mean max-match (across alive cells eating this step, 0..25) | [0, 25] | blue |
| 3 (bot)   | Fraction of alive cells at Negene = 8 | [0, 1]  | white |

Source: `evoca_egenome_stats(out)` fills a 5-element float buffer with
`[mean_negene, std_negene, distinct_egene_values, mean_max_match,
frac_at_max]`; `sim.egenome_stats()` returns the same data as a dict.

Reading the strips together:

- Mean Negene rising while distinct egenes plateau ⇒ population is
  consolidating around a fixed set of egene patterns and adding
  duplicates / shifts within that set, not exploring new ones.
- Mean max-match climbing without Negene rising ⇒ the LUT is evolving
  toward neighbourhoods that better match the *existing* egenes
  (eating gets better at fixed breadth).
- `frac at max` > 0 with `tax_per_egene > 0` is a sign the per-egene
  tax isn't strong enough to bound breadth; with `tax_per_egene = 0`
  it's expected to climb toward 1.

---

## Activity Quantile Probe (q_activity)

Shows the distribution of LUT genome activity over time as 9 decile curves
(p10 through p90) on a log-scaled Y axis.  Inspired by Figure 3 ("Activity
Profile") of `Docs/genelife-2.pdf`.

At each time step, `evoca_q_activity_deciles()` collects the cumulative
activity of every currently-alive genome, normalises each by the diversity D
(number of alive species), sorts the values, and extracts the 10th through
90th percentiles.  The 9 decile values are written to shared memory as
float32 tracks.

The SDL window renders each decile as a colored line using a blue-to-red
gradient: p10 (blue), p20, p30 (cyan), p40 (teal), p50 (green/median),
p60 (yellow), p70 (orange), p80, p90 (red).  The Y axis is log10-scaled
and auto-ranged from the data, so both very young (low activity) and ancient
(high activity) genomes are visible simultaneously.

Requires `activity` probe data (the function calls `evoca_activity_update()`
if `activity` is not separately enabled).

---

## Click-to-Identify (eg_activity, eg_pop)

Clicking on either the `eg_activity` or `eg_pop` probe window identifies
the egenome under the cursor and displays a 5x5 pattern overlay in the
top-right corner of the window.

**How it works**: the click handler reads the pixel color at the clicked
position in the shared-memory pixel buffer and matches it against the
known egenome color table.  For `eg_activity`, both full (alive) and
dimmed (extinct, RGB x 0.15) colors are checked.  If the clicked pixel
is background, the handler walks upward row by row until it finds a
colored pixel.

**Overlay**: a 25x25 pixel (5x5 cells at 5px each) rendering of the
matched egenome's D4-symmetric fiducial pattern, drawn with a 1px white
border using the egenome's assigned color.

**Console output**: each click prints the matched egenome index, its
binary representation, current population fraction, and alive/extinct
status to the subprocess stdout (relayed to the notebook cell output).

**Implementation note**: click coordinates are captured inside the
`SDL_PollEvent` loop but processing is deferred to the main render loop,
outside the event poll.  This avoids a macOS-specific segfault triggered
by accessing numpy shared-memory arrays during SDL event dispatch.

---

## Grouped Time-Series (ts)

Six scalar traces in two stacked strips of three traces each:

| Strip | Trace           | Source                                              |
|-------|-----------------|-----------------------------------------------------|
| top   | `pop`           | alive-cell count                                    |
| top   | `F_env`         | mean of environmental food F over the whole lattice |
| top   | `f_priv`        | mean of private food f over alive cells only        |
| bot   | `lut_div`       | distinct live LUT genomes (FNV-1a hash buckets)     |
| bot   | `eg_ent`        | Shannon entropy (bits) of the 64-bucket egenome distribution |
| bot   | `activity_flux` | flux-probe slope sum in the `[p20, p30]` activity band |

Each trace auto-scales independently within its strip: the observed
[min, max] over the visible window maps to the middle 80% of the
strip's height (10% padding top and bottom).  Zero samples are treated
as uninitialised and excluded from both autoscale and plotting.
Adjacent samples are connected by vertical fill so traces render as
continuous lines.  A small color-swatch legend sits at the top-left of
each strip; the right edge of each strip is the cursor (newest data).

**Note on `eg_ent`**: the egenome is a 6-bit D4-symmetric fiducial with
only 64 possible values, so a raw distinct-count saturates almost
immediately under any mutation pressure (all 64 buckets fill and the
trace flatlines).  Shannon entropy over the 64-bucket distribution
ranges 0..6 bits and stays responsive to mutation rate even after all
buckets are occupied — higher entropy means a more even distribution
across egenome space, lower entropy means one or a few egenomes
dominate.

---

## Y-axis Scale Controls

The `activity`, `eg_activity`, and `pat_activity` probes each have
halve/double buttons (`<| name |>`) to adjust their Y-axis saturation
scale (`act_ymax`, `eg_act_ymax`, `pat_act_ymax`).  These replace the
previous slider widgets.  Halving makes waves rise faster; doubling
spreads out low-activity entries.

---

## Reproduction Age Histogram

Tracks the distribution of time between successive reproduction events
(or birth-to-first-reproduction).  A per-cell timestamp `last_event_step[N*N]`
records the step of each cell's most recent birth or reproduction.  At each
reproduction, `age = step - last_event_step[parent]` is binned into
`repro_age_hist[1024]`.

A configurable `repro_age_t0` (default 0) skips transient: only events
where both the current step and the parent's last event are >= t0 are
counted.  `reset_repro_age_hist()` clears the histogram.
