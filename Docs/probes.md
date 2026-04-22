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
| `lut_complexity` | 512 x 128     | Stacked area: green=n1 only, yellow=n1+n2, red=all     |
| `eg_pop`         | 512 x 128     | Stacked area: population fraction per egenome value    |
| `entropy`        | 512 x 128     | Local-pattern Shannon entropy over time                |
| `pat_activity`   | 512 x 256     | Local-pattern activity (scrolling hash-colored strip)  |

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
