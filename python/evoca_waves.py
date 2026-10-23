"""Wave and spiral diagnostics for an oscillating medium (Chem CHEM-N1;
shared with Reaction-diffusion RXD-N6).

The cycle phase of a coarse block is theta = atan2(W - W0, F - F0), from its
mean food F and waste W minus a center (F0, W0), in absolute units (not
standardized). In the food -> colony -> waste -> food cycle each block goes
around the (F, W) plane once per period, so theta is its position in the
cycle. The center should be each block's TIME-averaged state (pass it as
`center`). The default, the frame's spatial mean, makes a synchronized bulk
oscillation invisible (every block sits at the mean), so it sees only
spatial heterogeneity.

- `phase_field`: theta and the cycle amplitude on the b x b coarse grid.
- `defects`: plaquettes where theta winds by +-2 pi (phase singularities).
  Noise makes many, in close +- pairs.
- `loop_winding`: the winding of theta around a square loop of half-size L.
  Measured on smooth random fields, it does NOT separate noise defects from
  cores once defects are dense: a large loop encloses many, and the net is
  +-1 about half the time. Kept for inspection, not used as a filter.
- `cores`: defects in surroundings with cycle amplitude >= amp_min.
- `matched`: defects of a later frame with a same-sign defect of an earlier
  frame nearby. A spiral core persists for many periods while the phase
  turns around it; noise defects come and go. The discriminator is the
  survival of matches at lags of several periods against a chance level
  (the earlier frame shifted by half the lattice).
- `order_parameter`: Kuramoto R = |mean exp(i theta)| over amplitude-masked
  blocks. R near 1 is bulk (synchronized) pulsing; waves and spirals give
  low R with large coherent phase domains.
- `period_of`: dominant period of a time series from its autocorrelation.
"""
import numpy as np


def _wrap(d):
    return (d + np.pi) % (2 * np.pi) - np.pi


def block_means(A, b):
    """Mean of A over each b x b block: the (N // b)^2 coarse grid."""
    n = A.shape[0] // b
    return A[:n * b, :n * b].reshape(n, b, n, b).mean(axis=(1, 3))


def phase_field(F, W, b, center=None):
    """(theta, amplitude) on the (N // b)^2 coarse grid. center = (F0, W0),
    scalars or coarse-grid arrays (e.g. time-averaged block_means); None
    uses the frame's spatial means (see the module note)."""
    Fc, Wc = block_means(F, b), block_means(W, b)
    F0, W0 = (Fc.mean(), Wc.mean()) if center is None else center
    x, y = Fc - F0, Wc - W0
    return np.arctan2(y, x), np.hypot(x, y)


def defects(theta):
    """Integer winding (+1, -1, 0) of every plaquette (i, j)-(i+1, j+1), periodic."""
    t01 = np.roll(theta, -1, axis=1)
    t11 = np.roll(t01, -1, axis=0)
    t10 = np.roll(theta, -1, axis=0)
    wind = _wrap(t01 - theta) + _wrap(t11 - t01) + _wrap(t10 - t11) + _wrap(theta - t10)
    return np.rint(wind / (2 * np.pi)).astype(int)


def loop_winding(theta, i, j, L):
    """Winding of theta around the square loop of half-size L enclosing
    plaquette (i, j), i.e. the corners (i - L + 1 .. i + L) x (j - L + 1 .. j + L)."""
    n = theta.shape[0]
    lo_i, hi_i, lo_j, hi_j = i - L + 1, i + L, j - L + 1, j + L
    path = ([(lo_i, c) for c in range(lo_j, hi_j)] + [(r, hi_j) for r in range(lo_i, hi_i)] +
            [(hi_i, c) for c in range(hi_j, lo_j, -1)] + [(r, lo_j) for r in range(hi_i, lo_i, -1)])
    th = np.array([theta[r % n, c % n] for r, c in path])
    return int(np.rint(_wrap(np.roll(th, -1) - th).sum() / (2 * np.pi)))


def cores(theta, amp, amp_min=0.1):
    """Positions (k x 2) and signs of the defects whose 2 x 2 corner blocks
    have mean cycle amplitude >= amp_min."""
    q = defects(theta)
    a = (amp + np.roll(amp, -1, 0) + np.roll(amp, -1, 1) + np.roll(np.roll(amp, -1, 0), -1, 1)) / 4
    q[a < amp_min] = 0
    return np.argwhere(q != 0), q[q != 0]


def matched(pos_a, sgn_a, pos_b, sgn_b, n, tol=2):
    """Number of entries of a with a same-sign entry of b within tol coarse
    cells (periodic Chebyshev distance)."""
    if len(pos_a) == 0 or len(pos_b) == 0:
        return 0
    d = np.abs(pos_a[:, None, :] - pos_b[None, :, :])
    d = np.minimum(d, n - d).max(axis=2)
    return int(((d <= tol) & (sgn_a[:, None] == sgn_b[None, :])).any(axis=1).sum())


def order_parameter(theta, amp, amp_min=0.1):
    """Kuramoto R over blocks with amplitude >= amp_min (nan if none)."""
    m = amp >= amp_min
    return float(np.abs(np.exp(1j * theta[m]).mean())) if m.any() else float('nan')


def period_of(series, dt):
    """(period in steps, autocorrelation peak): the first local maximum of
    the autocorrelation after its first zero crossing; period nan if that
    peak is below 0.2."""
    x = np.asarray(series, float) - np.mean(series)
    if len(x) < 4 or x.std() == 0:
        return float('nan'), 0.0
    acf = np.correlate(x, x, 'full')[len(x) - 1:]
    acf /= acf[0]
    half = len(acf) // 2
    neg = np.where(acf[:half] < 0)[0]
    if len(neg) == 0:
        return float('nan'), 0.0
    for k in range(neg[0] + 1, half - 1):
        if acf[k] >= acf[k - 1] and acf[k] >= acf[k + 1]:
            peak = float(acf[k])
            return (float(k * dt) if peak >= 0.2 else float('nan')), peak
    return float('nan'), 0.0
